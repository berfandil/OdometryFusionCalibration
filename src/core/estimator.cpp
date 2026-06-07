// ofc/core/estimator.cpp — the caller-pumped fusion facade (Slice 2 + Slice 5 time-sync).
//
// STRICT CORE: all working memory lives in Impl and is sized from Config at init();
// step() allocates nothing; no exceptions; bounded loops.
//
// TIME-SYNC (Slice 5, D16). When cfg.timesync_enabled, each step() samples every
// source's ‖ω‖ over a small sub-interval ending at the frontier, feeds it (keyed to
// base time) into TimeSync, and runs a bounded update() (‖ω‖ xcorr -> per-source offset
// histogram). A source's committed offset (histogram concentration >= commit_concentration)
// then SHIFTS its fusion query interval earlier by `off` so its internal +off re-shift
// lands back on true base time (sign per D21 — REMOVES a planted offset); below the gate,
// the configured prior is used. The offset + confidence are surfaced in CalibSnapshot.
// With time-sync OFF the per-source offset is exactly the prior (default 0) — Slice-2
// behavior unchanged.
//
// Slice 2 pipeline (no extrinsic/scale calibration, no absolute corrections — DESIGN §§4-5):
//   per step(now): frontier t1 = now - fusion_delay
//     integration interval [q0, t1]:
//       q0 = last_t1        (steady state — picks up where the last fuse ended)
//       q0 = t1 - window_s  (bootstrap of the FIRST fuse — no prior frontier yet)
//     dt = (t1 - q0) seconds                   (ACTUAL elapsed motion, not a fixed step)
//     for each registered source:
//       B = query(q0 - off, t1 - off)         (source-frame delta; off = committed/prior
//                                               clock offset, removed per Slice 5)
//       B_corr = { B.R, B.t / prior_scale }   (de-scale the reported translation, D20)
//       A = X o B_corr o X^-1                  (frame-align to base; X = sensor->base
//                                               prior extrinsic — see convention below)
//       w = clamp(weight_prior * sigma_confidence)
//     median of {A_i} weighted by {w_i}        (split-metric Weiszfeld)
//     ESKF predict on the median delta (dt = t1 - q0), adaptive Q from the spread
//     on success: last_t1 <- t1               (gap/overlap-free across ticks)
//     populate Result.frontier and Result.tip
//   Returns NotReady until >= 1 source covers the interval (full lifecycle is Slice 3).
//   NOTE: `window_s` is the bootstrap/lookback interval for the first fuse only; the
//   integrator no longer assumes one window per tick, so tick cadence may differ from
//   window_s without opening integration gaps/overlaps (review fix — DESIGN §7).
//
// Frame-align convention (documented): each SensorConfig::prior_extrinsic X is the
// sensor->base transform (a point/twist expressed in the sensor frame maps to base
// via X). A motion B measured in the sensor frame becomes, in the base frame, the
// conjugation A = X o B o X^-1. With X = identity (aligned mount) A == B.
#include "ofc/core/estimator.hpp"

#include "ofc/core/buffer.hpp"
#include "ofc/core/calibration.hpp"
#include "ofc/core/eskf.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/median.hpp"
#include "ofc/core/smoother.hpp"
#include "ofc/core/timesync.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
constexpr Scalar kNanosPerSec = Scalar(1e9);
constexpr int    kMaxSourcesCap = 32;   // matches Result::calib/health array sizes
// Absolute-reference corrections cap (Slice 11, D22). Fixed-capacity registry (strict
// core: no post-init heap, bounded loops). 8 is generous for the GPS/INS/map-match
// plugin count a single rig carries; past it add_correction returns CapacityExceeded.
constexpr int    kMaxCorrectionsCap = 8;

// Variance-EMA reliability (Slice 9, D17). A source's residual-to-consensus EMA mean +
// variance accrue per fused step; once it has at least kRelWarmup samples the reliability
// multiplier is recomputed from the ratio of a robust baseline variance (the MEDIAN of the
// warmed-up participants' variances) to its own variance. Before warmup — and whenever
// fewer than 2 participants have warmed up — the reliability holds at 1.0 (Slice-2 weights
// unchanged), so a short run never perturbs fusion.
constexpr int    kRelWarmup = 20;       // min residual samples before reliability applies
constexpr Scalar kRelVarEps = Scalar(1e-12);   // guards the ref_var / resid_var ratio

Timestamp secs_to_ns(Scalar s) {
    return static_cast<Timestamp>(std::llround(s * kNanosPerSec));
}

// Sigma-confidence of a windowed delta: an inverse-covariance scalar. We use
// 1 / (mean diagonal + eps) so a tighter (smaller-Sigma) source weighs more. This
// is the prior x Sigma-confidence weight of DESIGN §4.
// PLACEHOLDER (unit-mixing): the mean over all 6 diagonal entries blends
// translation (m^2) and rotation (rad^2) variances into one scalar, so a source
// with tiny rad var but huge m var gets a misleading confidence. This is an
// accepted Slice-2 placeholder until the Slice 9 reliability EMA replaces it with a
// properly unit-separated weight.
Scalar sigma_confidence(const Mat6& cov) {
    const Scalar mean_diag = cov.diagonal().mean();
    constexpr Scalar kEps = Scalar(1e-9);
    return Scalar(1) / (std::max(Scalar(0), mean_diag) + kEps);
}

// ‖ω‖ extraction for time-sync (DESIGN §6 time-offset row, D16): the body angular-rate
// MAGNITUDE over a small sub-interval [t-h, t]. From the source's reported SE(3) delta
// over that window, ‖ω‖ ≈ ‖log(ΔR)‖ / h. This is extrinsic-invariant (a rotation
// magnitude survives the sensor->base conjugation), so it can be sampled per source
// BEFORE any spatial calibration and cross-correlated against the reference. dt is the
// sub-interval length in seconds (> 0).
Scalar omega_norm_from_delta(const SE3& motion, Scalar dt) {
    if (!(dt > Scalar(0))) return Scalar(0);
    return so3::log(motion.R).norm() / dt;
}

// Commit gate for an EXTRINSIC/SCALE DOF that gets RE-ANCHORED on commit (Slice 8). Same
// τ_commit ∧ N_min first-commit + τ_drop hysteresis as commit_gate(), but the re-open
// (committed -> uncommitted) additionally requires the histogram to be RE-POPULATED
// (votes >= N_min). Without this guard, the re-anchor's histogram RESET drops vote mass to
// ~0, which momentarily reads confidence 0 < τ_drop and would falsely re-open — then re-fill
// and re-commit, RE-ANCHOR-ing again every few steps (commit thrash). Requiring a full
// histogram before honoring a confidence drop makes the re-open mean a GENUINE drift (a full,
// low-concentration histogram), not the post-reset transient. The reference / offset DOF keep
// using the plain commit_gate().
bool commit_gate_reanchor(bool prev_committed, Scalar confidence, Scalar votes,
                          Scalar commit_concentration, int commit_min_votes,
                          Scalar commit_drop) {
    if (prev_committed) {
        // Stay committed through a sparse (just-reset / still-refilling) histogram; only
        // re-open on a FULL histogram that has genuinely lost concentration.
        if (votes < static_cast<Scalar>(commit_min_votes)) return true;
        return confidence >= commit_drop;
    }
    return confidence >= commit_concentration &&
           votes >= static_cast<Scalar>(commit_min_votes);
}
} // namespace

// All working memory lives here and is sized from Config at init().
struct Estimator::Impl {
    Config cfg{};
    Result result{};
    bool   inited = false;

    // Registered sources (pointers owned by caller). Fixed array sized by the cap;
    // the live count is bounded by cfg.max_sources.
    const ISource* sources[kMaxSourcesCap] = {};
    int            source_count = 0;

    // Registered absolute-reference corrections (Slice 11, D22). Pointers owned by the
    // caller; fixed-capacity registry (strict core). Evaluated in step() AFTER predict and
    // BEFORE publishing, each Mahalanobis-gated by cfg.mahalanobis_chi2. Reset in init().
    const ICorrection* corrections[kMaxCorrectionsCap] = {};
    int                correction_count = 0;

    // Per-source prior extrinsic (sensor->base), prior translation scale, and prior
    // weight, looked up from the SensorConfig matching the source id (identity / 1.0 /
    // 1.0 if no config provided). prior_scale de-scales the source's reported
    // translation before the frame-align (D20 — inverts the sim's per-source scale).
    SE3    prior_extrinsic[kMaxSourcesCap];
    Scalar prior_scale[kMaxSourcesCap];
    Scalar weight_prior[kMaxSourcesCap];
    // Per-source prior clock offset (seconds). Sign per D21 / CONFIG §9: positive =>
    // source clock ahead of base; its reported [t0,t1] reads true [t0+off, t1+off].
    Scalar prior_time_offset[kMaxSourcesCap];

    // Time-sync (Slice 5, D16). When cfg.timesync_enabled, each step() feeds every
    // source's current ‖ω‖ sample into `timesync` and periodically runs update(); the
    // committed per-source offset shifts that source's fusion query interval so the
    // planted clock skew is removed before fusion. Heap-resident inside Impl (the
    // strict-core no-post-init-heap contract: Impl is allocated once in init()).
    TimeSync timesync;
    bool     timesync_active = false;
    int      ticks_since_sync = 0;

    // Phase-1 calibration (Slice 6, D5/D11/D20). Driven from the straight-regime calib
    // stage each fused step: per source the DE-SCALED reported delta B_corr + the fused
    // consensus twist/translation -> 3-channel so(3) histogram @ prior basepoint (yaw,
    // pitch) + magnitude-ratio scale histogram. Slice 6 only ESTIMATES + exposes into
    // CalibSnapshot.{extrinsic, scale, confidence}; feeding it back into fusion is Slice 8
    // (NOT wired). Heap-resident inside Impl (allocated once in init()).
    Phase1Calibrator calib1;
    // Phase-2 calibration (Slice 7, D5/D10/D21). Driven from the TURN-regime calib stage:
    // per source, hand-eye A·X = X·B recovers the roll (about forward, 1-D rotation
    // residual -> circular S¹) + the xyz lever arm (linear normal-equation LS), on top of
    // Phase 1's yaw/pitch. Slice 7 only ESTIMATES + exposes into CalibSnapshot.{extrinsic
    // (full rotation + translation), extrinsic_confidence, translation_confidence}; no
    // feedback into fusion (Slice 8). Heap-resident inside Impl (allocated once in init()).
    Phase2Calibrator calib2;
    // Per-step scratch for the calibration observe() (no heap in step()).
    SourceId calib_ids[kMaxSourcesCap];
    // Registered slot of each covered (calibration) source, aligned with calib_ids — lets
    // the Slice-10 smoother push/look-up the per-source slot without an id re-scan.
    int      calib_slot[kMaxSourcesCap];
    SE3      calib_reported[kMaxSourcesCap];
    // Per-source Σ-confidence aligned with calib_ids — the same inverse-covariance scalar
    // used as the fusion weight (PRE clamp/prior, the raw confidence). Fed into both
    // calibrators' observe() so the D5 vote_weight Confidence/Combo modes can weight each
    // vote by how trustworthy the source's window is.
    Scalar   calib_conf[kMaxSourcesCap];

    // ---- Per-sensor fixed-lag RTS smoother (Slice 10, D18) ----------------------------
    // The DEEPER calibration frontier. When ANY source has SensorConfig::per_sensor_kf, the
    // calibration feed runs at now − delay_fast − L (L = calib_lag_s · tick_rate_hz steps):
    // each step's calibration inputs are buffered in a depth-(L+1) ring, and the lag-L-OLD
    // entry drives the calibrators — with each smoothing-enabled source's twist replaced by
    // the TWO-SIDED RTS-smoothed twist (sharper, zero-phase). FUSION is untouched (causal).
    //
    // smoother_active is false (and the whole deeper path is skipped — calibration runs at
    // the causal frontier exactly as Slice 6/7/8) UNLESS at least one source enables
    // per_sensor_kf. So the DEFAULT (every per_sensor_kf == false) is byte-identical to
    // pre-Slice-10. Heap-resident inside Impl (allocated once in init()).
    TwistSmoother smoother;
    bool          smoother_active = false;            // any source has per_sensor_kf
    bool          smooth_source[kMaxSourcesCap] = {}; // per registered slot: per_sensor_kf
    int           calib_lag_steps = 0;                // L (clamped), == smoother.lag_steps()

    // The deeper-frontier calibration-input ring (depth L+1). Each filled slot is one past
    // step's full calibration input set: every covered source's de-scaled B_corr + its
    // Σ-confidence + measured twist + the source id, the covered count, the consensus
    // (fused_omega/fused_trans/med SE3), and the step dt. Sized at the compile-time caps so
    // init() is the sole allocator (strict core). Logical order: 0 = oldest .. fill-1.
    static constexpr int kCalibRing = TwistSmoother::kMaxLag + 1;
    struct CalibFrame {
        int      nc = 0;
        SourceId ids[kMaxSourcesCap];
        SE3      reported[kMaxSourcesCap];     // de-scaled B_corr (raw, delayed)
        Scalar   conf[kMaxSourcesCap];
        int      slot[kMaxSourcesCap];         // registered slot of each covered source
        // Time-alignment stamp (Slice-10 review MAJOR). For each covered entry i that is a
        // smoothing source, the smoother's MONOTONIC push-count for that slot AFTER this
        // frame's push (push_seq[slot] post-increment). The smoother's emitted sample after
        // N pushes is the (N-L)-th push, so at consumption the smoothed twist time-aligns
        // with THIS frame iff (current push_seq[slot] - push_at[i]) == L. Under per-source
        // dropout that equality fails and we fall back to the raw delayed delta (a clean
        // variance loss, NOT a calibration bias). -1 for non-smoothing / un-pushed entries.
        long     push_at[kMaxSourcesCap];
        Vec3     fused_omega = Vec3::Zero();
        Vec3     fused_trans = Vec3::Zero();
        SE3      fused_motion;
        Scalar   dt = Scalar(0);
    };
    CalibFrame calib_ring[kCalibRing];
    int        calib_ring_count = 0;           // filled entries (saturates at L+1)
    // Per-registered-slot MONOTONIC smoother push count (Slice-10 review MAJOR): incremented
    // once each time this slot is pushed into the smoother (covered AND smoothing-enabled).
    // Unlike TwistSmoother::Slot::count (which saturates at the ring depth) this never
    // saturates, so push_seq - push_at is the true number of pushes between a frame and its
    // consumption — the exact L-step alignment test. Reset to 0 in init().
    long       push_seq[kMaxSourcesCap] = {};
    // Per-step scratch the deeper path writes the lag-L entry's (possibly smoothed) deltas
    // into before calling observe() (no heap in step()).
    SE3        calib_smoothed[kMaxSourcesCap];

    // Per-source time-offset commit state (DESIGN §2/§6). A source's offset is COMMITTED
    // (driving fusion) once its histogram concentration clears commit_concentration AND
    // its vote count clears commit_min_votes (N_min); once committed it STAYS committed
    // (hysteresis) until concentration falls below commit_drop. Fixed-capacity bool array
    // (strict core); advanced once per step() in update_commit_state(). Indexed by the
    // registered-source slot. Always false for the reference / when time-sync is off.
    bool     offset_committed[kMaxSourcesCap] = {};

    // Per-source, per-DOF EXTRINSIC/SCALE commit state (Slice 8 feedback loop, DESIGN §6).
    // Same N_min + hysteresis gate (reuse commit_gate) on each DOF's own confidence + vote
    // count. On the RISING EDGE of a commit the freshly-committed value is swapped into the
    // fusion prior (prior_extrinsic / prior_scale) AND the calibrator is RE-ANCHORED to that
    // value with the DOF's histogram RESET (so the stale votes — cast @ the old basepoint —
    // do not pin the now-residual mode; subsequent votes refine the residual around the new
    // basepoint — the two-rate iterate, DESIGN §6). Once committed, every step keeps the
    // fusion prior tracking the calibrator's latest committed estimate (basepoint ∘ small
    // residual). Indexed by the registered-source slot. Reference never commits its extrinsic
    // (it is the pinned gauge). Advanced once per step() in apply_calib_feedback() AFTER the
    // result is published — so a step() never sees a half-updated prior (atomic swap).
    bool     ext_committed[kMaxSourcesCap]   = {};   // yaw/pitch (so(3) direction)
    bool     roll_committed[kMaxSourcesCap]  = {};   // roll about forward
    bool     lever_committed[kMaxSourcesCap] = {};   // xyz lever arm
    bool     scale_committed[kMaxSourcesCap] = {};   // per-source scale

    // Scratch for the per-step fuse (no heap in step()).
    SE3    aligned[kMaxSourcesCap];
    Scalar weights[kMaxSourcesCap];
    // Maps median entry k -> registered source slot (cold-start may exclude sources from
    // the median, so the median index no longer matches the covered order). Used by the
    // per-source residual diagnostics to attribute the consensus distance correctly.
    int    med_slot[kMaxSourcesCap];

    // Variance-EMA reliability state (Slice 9, D17). Per registered-source slot, accrued
    // each fused step from that source's residual-to-consensus d (the split_distance used
    // by the residual diagnostics). The D17 BIAS/VARIANCE SPLIT:
    //   * resid_mean = EMA mean of d  -> the SYSTEMATIC component (a biased source's
    //                  residual is large but CONSISTENT, so it lands here, NOT in the
    //                  variance). Surfaced as SourceHealth::bias; routed to the calibrator,
    //                  never folded into the weight.
    //   * resid_var  = EMA variance of d AROUND its running mean -> the zero-mean SCATTER.
    //                  A genuinely NOISY source has large scatter; a biased-but-low-noise
    //                  source has SMALL scatter (its bias is in resid_mean) -> NOT
    //                  penalized. reliability = clamp(ref_var / max(resid_var, eps), floor,
    //                  cap) downweights only true scatter.
    //   * resid_n    = sample count (for the kRelWarmup gate).
    //   * reliability= the applied multiplier (init 1.0; held at 1.0 pre-warmup / for
    //                  non-participants). Reset in init().
    Scalar resid_mean[kMaxSourcesCap];
    Scalar resid_var[kMaxSourcesCap];
    int    resid_n[kMaxSourcesCap];
    Scalar reliability[kMaxSourcesCap];
    // Fixed stack scratch for the robust baseline (median of warmed-up variances). Sized at
    // the cap so the per-step median selection is heap-free + bounded (strict core).
    Scalar rel_scratch[kMaxSourcesCap];

    // Lifecycle state (Slice 3, D23). Persisted across steps; the phase is
    // RECOMPUTED every step from the current fuse signals, so an upgrade
    // (DEGRADED->NOMINAL) and a graceful downgrade (NOMINAL->DEGRADED on source
    // loss) both fall out automatically. `ever_fused_` latches on the first
    // successful fuse and distinguishes the pre-first-fuse WARMUP from the
    // post-fuse DEGRADED of a total source loss. Reset in init().
    Phase  phase_      = Phase::Init;
    bool   ever_fused_ = false;

    // Set both result.phase and result.readiness from one Phase (the coarse [0,1]
    // readiness map of DESIGN §7/§8). Centralizes the mapping so the two fields
    // never drift apart.
    void publish_phase(Phase p) {
        phase_ = p;
        result.phase = p;
        switch (p) {
            case Phase::Init:     result.readiness = Scalar(0.0);  break;
            case Phase::Warmup:   result.readiness = Scalar(0.25); break;
            case Phase::Degraded: result.readiness = Scalar(0.6);  break;
            case Phase::Nominal:  result.readiness = Scalar(1.0);  break;
            // Keep the phase->readiness map TOTAL: a future Phase value falls back to a
            // safe "not ready" readiness rather than silently keeping the last value.
            default:              result.readiness = Scalar(0.0);  break;
        }
    }

    Eskf   eskf;
    bool   eskf_started = false;

    // End of the last successfully-fused window. The next predict integrates the
    // actual motion over [last_t1, t1] (dt = t1 - last_t1), so integration is
    // gap/overlap-free regardless of how the caller spaces step() ticks. Only
    // advanced on a successful fuse; a step with no covering source leaves it put.
    Timestamp last_t1     = 0;
    bool      has_last_t1 = false;

    // Resolve the SensorConfig for a given source id (linear scan, bounded by
    // sensor_count). Returns nullptr if none.
    const SensorConfig* sensor_for(SourceId id) const {
        if (cfg.sensors == nullptr) return nullptr;
        for (int i = 0; i < cfg.sensor_count; ++i) {
            if (cfg.sensors[i].id == id) return &cfg.sensors[i];
        }
        return nullptr;
    }

    // Advance every source's time-offset commit state for this step (DESIGN §2/§6).
    // Runs the N_min + hysteresis gate on the current histogram reading; the resulting
    // per-source bool then selects committed-estimate vs prior in effective_offset() and
    // is surfaced as CalibSnapshot::committed. No-op (all false) when time-sync is off.
    // Bounded by source_count; no heap.
    void update_commit_state() {
        for (int i = 0; i < source_count; ++i) {
            if (sources[i] == nullptr) { offset_committed[i] = false; continue; }
            const SourceId id = sources[i]->id();
            if (!timesync_active || id == cfg.reference_sensor_id) {
                offset_committed[i] = false;          // reference / time-sync off never commits
                continue;
            }
            offset_committed[i] = commit_gate(
                offset_committed[i], timesync.confidence(id), timesync.vote_count(id),
                cfg.commit_concentration, cfg.commit_min_votes, cfg.commit_drop);
        }
    }

    // Effective per-source clock offset (seconds) applied to its fusion query. When the
    // source's offset estimate is COMMITTED (the N_min + hysteresis gate in
    // update_commit_state has latched), use the estimated offset; otherwise fall back to
    // the configured prior (default 0). Slot `i` indexes the registered-source arrays.
    Scalar effective_offset(int i) const {
        if (timesync_active && offset_committed[i]) {
            return timesync.offset(sources[i]->id());
        }
        return prior_time_offset[i];
    }

    // Whether a source participates in the fusion median this step (COLD-START switch, D6).
    //   * MedianFromStart : every covered source participates from the first tick (the
    //                       Slice-2 behaviour — used as the bootstrap-loop default since the
    //                       median needs ≥3 sources to reject an outlier and to give the
    //                       calibrators a consensus to vote against).
    //   * ReferenceOnly   : before a source's extrinsic is confident, only the reference
    //                       (the pinned gauge — extrinsic = prior by construction) drives the
    //                       output (reference-only dead-reckon); a non-reference source joins
    //                       the median once its extrinsic has COMMITTED (its prior is now
    //                       trustworthy). The reference always participates.
    // Slot `i` indexes the registered-source arrays. (Full lifecycle is Slice 3; this is the
    // minimal cold-start interplay the brief asks for — it must not regress Slice-2 when
    // cold_start == MedianFromStart.)
    bool participates(int i) const {
        if (cfg.cold_start == ColdStart::MedianFromStart) return true;
        const SourceId id = sources[i]->id();
        if (id == cfg.reference_sensor_id) return true;
        return ext_committed[i];   // ReferenceOnly: join once the extrinsic is committed
    }

    // Drive BOTH calibrators (Phase 1 + Phase 2) for one set of calibration inputs — the
    // shared observe core used by the causal path (Slice 6/7) AND the deeper-frontier path
    // (Slice 10). `reported` are the (possibly RTS-smoothed) de-scaled B_corr deltas, `omega`
    // / `trans` / `motion` the consensus paired with them (causal step's consensus, or the
    // lag-L-old consensus on the deeper path — kept TIME-ALIGNED with `reported`). This is
    // exactly the sequence the Slice-6/7 call site ran inline; factoring it lets the smoother
    // swap in the lagged + smoothed inputs without duplicating the observe wiring. Heap-free.
    void run_calibration_observe(int nc_in, const SourceId* ids, const SE3* reported,
                                 const Vec3& omega, const Vec3& trans, const SE3& motion,
                                 const Scalar* conf) {
        calib1.observe(nc_in, ids, reported, omega, trans, conf);
        for (int i = 0; i < nc_in; ++i) {
            calib2.set_yaw_pitch(ids[i], calib1.extrinsic(ids[i]).R);
        }
        calib2.observe(nc_in, ids, reported, motion, omega, conf);
    }

    // Apply the calibration->fusion feedback (Slice 8) for one step, AFTER the result is
    // published — the atomic between-step swap (a step() never sees a half-updated prior).
    // Per source + DOF: advance the commit gate (commit_gate_reanchor — N_min + hysteresis,
    // with a re-fill guard so a post-reset sparse histogram does not falsely re-open). Once
    // committed, swap the latest committed value into the fusion prior every step (extrinsic
    // rotation, lever arm) so fusion tracks the refined estimate; the SCALE DOF additionally
    // re-anchors on its rising edge (folds the residual into prior_scale + resets its
    // histogram — the one DOF whose votes read the fusion prior). The reference's
    // extrinsic/scale are the pinned gauge — never committed. Bounded by source_count; no heap.
    void apply_calib_feedback() {
        for (int i = 0; i < source_count; ++i) {
            if (sources[i] == nullptr) continue;
            const SourceId id = sources[i]->id();
            // The reference is the gauge: its extrinsic + scale stay at the prior, never
            // commit, never re-anchor.
            if (id == cfg.reference_sensor_id) {
                ext_committed[i] = roll_committed[i] = false;
                lever_committed[i] = scale_committed[i] = false;
                continue;
            }

            // NOTE on commit_min_votes (Slice-8 review MINOR): the gate compares it against
            // vote_count()/roll_vote_count()/etc., which are histogram TOTALS = vote MASS, not
            // a count, under any non-`One` vote_weight (Confidence/Rotation/Combo scale each
            // vote by the Σ-confidence and/or ‖ω‖). So commit_min_votes is a MASS threshold
            // that must be RE-TUNED per vote_weight; it is a literal vote COUNT only when
            // vote_weight == One (which the feedback tests force — set_hists()). See config.hpp.

            // --- yaw/pitch (so(3) direction) -----------------------------------------
            const bool was_ext = ext_committed[i];
            ext_committed[i] = commit_gate_reanchor(
                was_ext, calib1.extrinsic_confidence(id), calib1.vote_count(id),
                cfg.commit_concentration, cfg.commit_min_votes, cfg.commit_drop);

            // --- roll ----------------------------------------------------------------
            const bool was_roll = roll_committed[i];
            roll_committed[i] = commit_gate_reanchor(
                was_roll, calib2.extrinsic_confidence(id), calib2.roll_vote_count(id),
                cfg.commit_concentration, cfg.commit_min_votes, cfg.commit_drop);

            // --- xyz lever arm -------------------------------------------------------
            const bool was_lever = lever_committed[i];
            lever_committed[i] = commit_gate_reanchor(
                was_lever, calib2.translation_confidence(id), calib2.xyz_vote_count(id),
                cfg.commit_concentration, cfg.commit_min_votes, cfg.commit_drop);

            // --- scale ---------------------------------------------------------------
            const bool was_scale = scale_committed[i];
            scale_committed[i] = commit_gate_reanchor(
                was_scale, calib1.scale_confidence(id), calib1.vote_count(id),
                cfg.commit_concentration, cfg.commit_min_votes, cfg.commit_drop);

            // === ROTATION + LEVER feedback (yaw/pitch ∘ roll, xyz) ===================
            // CONTRACTIVE EXTRINSIC RE-ANCHOR (Slice 8). calib1 recovers the forward AXIS as
            // the INVERSE minimal rotation δRᵀ @ its own basepoint, so the recovered extrinsic
            // satisfies  X.R·dir_B = e_x  exactly (calibration.cpp extrinsic()). Folding that
            // recovered yaw/pitch back into calib1's basepoint and resetting the so(3) votes is
            // therefore a CONTRACTIVE map: from a LARGE prior error the next round's g_obs lands
            // on e_x, δR→I, and the estimate converges to the planted forward axis (the genuine
            // two-rate iterate — same shape as the scale re-anchor below). calib2 recovers roll
            // @ the fed R_yp and the lever arm as an ABSOLUTE linear-LS solve; neither reads the
            // FUSION prior, so the lever swap is a value publish (no stale votes).
            //
            // Two parts, mirroring the scale DOF:
            //  (1) EVERY committed step — PUBLISH the latest committed value into the fusion
            //      prior so fusion frame-aligns with the refined estimate.
            //  (2) RISING EDGE of the yaw/pitch commit — fold the recovered yaw/pitch into
            //      calib1's basepoint + reset its so(3) histogram, so the residual re-votes
            //      ≈ 0 around the new (corrected) basepoint. Rising-edge only: re-anchoring
            //      every step would never let the histogram re-concentrate (collapsing
            //      confidence -> hysteresis re-open thrash).
            if (ext_committed[i] || roll_committed[i]) {
                // calib1.extrinsic().R is the recovered yaw/pitch (the contractive minimal
                // rotation — roll about forward leaves the forward axis fixed, so Phase 1 never
                // carries roll). Once yaw/pitch has committed AND been re-anchored,
                // calib1.extrinsic().R tracks the CORRECTED yaw/pitch, so a roll commit on top
                // composes on the corrected basepoint, not a stale one (the MINOR is resolved by
                // the contractive re-anchor). Before yaw/pitch commits, hold the fusion prior R.
                const Mat3 R_yp_committed = ext_committed[i] ? calib1.extrinsic(id).R
                                                             : prior_extrinsic[i].R;
                // calib2.extrinsic().R = R_yp · Rx(roll) folds the recovered roll on top of the
                // R_yp fed to Phase 2 every step (the same calib1 rotation), so it carries both.
                prior_extrinsic[i].R = roll_committed[i] ? calib2.extrinsic(id).R
                                                         : R_yp_committed;
            }
            if (lever_committed[i]) {
                prior_extrinsic[i].t = calib2.lever_arm(id);
            }
            // (2) Yaw/pitch RE-ANCHOR (rising edge of the yaw/pitch commit): move calib1's
            // so(3) basepoint to the recovered MINIMAL rotation (extrinsic().R — the gauge that
            // satisfies X.R·dir_B = e_x, so this composition is CONTRACTIVE) and drop the
            // now-stale residual votes (cast @ the OLD basepoint). The post-reset so3() falls
            // back to the basepoint until new votes land, so the residual re-concentrates ≈ 0
            // around the new (corrected) basepoint and calib1 tracks the planted forward axis
            // from a large initial error. Rising-edge only: re-anchoring every step would never
            // let the histogram re-concentrate (collapsing confidence -> hysteresis re-open
            // thrash) — the same shape as the scale re-anchor below. One contraction lands the
            // forward axis on truth (to the run's mode/segment-dynamics floor); the every-step
            // publish then keeps refining the residual around the fixed basepoint.
            if (ext_committed[i] && !was_ext) {
                SE3 bp = prior_extrinsic[i];
                bp.R = calib1.extrinsic(id).R;   // MINIMAL rotation (contractive basepoint)
                calib1.set_basepoint(id, bp);
                calib1.reset_so3(id);
            }
            // Keep calib2's prior (lever-arm fallback for an unobservable axis + the
            // PairwisePinnedRef gauge) tracking the committed extrinsic — value sync only,
            // NOT a histogram reset.
            if (ext_committed[i] || roll_committed[i] || lever_committed[i]) {
                calib2.set_basepoint(id, prior_extrinsic[i]);
            }

            // === SCALE feedback (the one DOF that DOES re-anchor) =====================
            // calib1 votes the magnitude ratio off the DE-SCALED B (B.t / prior_scale), so
            // its scale() is the RESIDUAL scale vs the CURRENT fusion prior_scale — the only
            // calibrated quantity that reads the fusion prior. The committed absolute scale =
            // prior_scale * residual. On the RISING EDGE, fold the residual into prior_scale
            // and RESET the scale histogram so it re-votes the (now ≈ 1) residual around the
            // new prior (the two-rate iterate — and the old ratios, taken vs the old
            // prior_scale, are genuinely stale). Rising-edge only: resetting every step would
            // never let the histogram re-concentrate (collapsing confidence -> hysteresis
            // re-open thrash).
            if (scale_committed[i] && !was_scale) {
                const Scalar residual = calib1.scale(id);
                if (residual > Scalar(0)) {
                    prior_scale[i] *= residual;
                    if (!(prior_scale[i] > Scalar(0))) prior_scale[i] = Scalar(1);
                    calib1.reset_scale(id);
                }
            }
        }
    }
};

Estimator::Estimator() : impl_(nullptr) {}

Estimator::~Estimator() { delete impl_; }   // sole heap free; nothing after init

Status Estimator::init(const Config& cfg) {
    const Status s = validate(cfg);
    if (!ok(s)) return s;
    if (impl_ == nullptr) impl_ = new Impl();   // allocate-once at init
    impl_->cfg          = cfg;
    impl_->result       = Result{};
    // Lifecycle reset (Slice 3): INIT is observable before the first step. publish_phase
    // sets both result.phase (== Init, as Result{} already does) AND result.readiness
    // (== 0.0), and resets the persisted phase_. ever_fused_ clears so the first
    // un-fusable step reads WARMUP, not the post-fuse DEGRADED.
    impl_->ever_fused_  = false;
    impl_->publish_phase(Phase::Init);
    impl_->source_count = 0;
    // Absolute-ref corrections reset (Slice 11): no plugins until re-registered.
    impl_->correction_count = 0;
    for (int i = 0; i < kMaxCorrectionsCap; ++i) impl_->corrections[i] = nullptr;
    impl_->eskf_started = false;
    impl_->last_t1      = 0;
    impl_->has_last_t1  = false;
    impl_->ticks_since_sync = 0;
    for (int i = 0; i < kMaxSourcesCap; ++i) {
        impl_->sources[i]           = nullptr;
        impl_->prior_extrinsic[i]   = SE3{};
        impl_->prior_scale[i]       = Scalar(1);
        impl_->weight_prior[i]      = Scalar(1);
        impl_->prior_time_offset[i] = Scalar(0);
        impl_->offset_committed[i]  = false;
        impl_->ext_committed[i]     = false;
        impl_->roll_committed[i]    = false;
        impl_->lever_committed[i]   = false;
        impl_->scale_committed[i]   = false;
        impl_->med_slot[i]          = -1;
        // Variance-EMA reliability reset (Slice 9, D17): clean EMA state + neutral
        // multiplier so the first kRelWarmup steps weight exactly as Slice 2.
        impl_->resid_mean[i]        = Scalar(0);
        impl_->resid_var[i]         = Scalar(0);
        impl_->resid_n[i]           = 0;
        impl_->reliability[i]       = Scalar(1);
        // Per-sensor smoother gate (Slice 10): bound in add_source(); default off.
        impl_->smooth_source[i]     = false;
        impl_->calib_slot[i]        = -1;
        // Monotonic smoother push count per slot (Slice-10 review MAJOR): empty.
        impl_->push_seq[i]          = 0;
    }
    // Deeper-frontier calibration ring + smoother (Slice 10): empty until the first step.
    impl_->calib_ring_count = 0;
    impl_->smoother.reset();

    // Time-sync: configure once at init (preallocates its buffers/histograms). Disabled
    // configs skip it entirely (behave exactly as before — offsets stay at the priors).
    impl_->timesync_active = false;
    if (cfg.timesync_enabled) {
        const Status ts = impl_->timesync.configure(cfg, cfg.reference_sensor_id);
        if (!ok(ts)) return ts;
        impl_->timesync_active = true;
    }

    // Phase-1 calibration: configure once at init (preallocates its per-source histograms).
    // Per-source priors are bound in add_source(). Always active (it self-gates on the
    // straight regime); it never feeds back into fusion in Slice 6.
    {
        const Status cs = impl_->calib1.configure(cfg, cfg.reference_sensor_id);
        if (!ok(cs)) return cs;
    }
    // Phase-2 calibration: configure once at init (preallocates its roll + xyz histograms
    // and per-source LS accumulators). Self-gates on the turn regime; never feeds back into
    // fusion in Slice 7 (Slice 8).
    {
        const Status cs = impl_->calib2.configure(cfg, cfg.reference_sensor_id);
        if (!ok(cs)) return cs;
    }

    // Per-sensor fixed-lag RTS smoother (Slice 10, D18). Scan the SensorConfigs: if ANY has
    // per_sensor_kf, the calibration feed switches to the DEEPER (lag-L) frontier and the
    // smoother is configured. With NONE enabled (the default for every existing test) the
    // smoother stays INACTIVE and the calibration call site is byte-identical to Slice 6/7/8
    // (the hard non-regression guard). L = round(calib_lag_s · tick_rate_hz), clamped to
    // [1, TwistSmoother::kMaxLag] by configure(). Process noise = the MAX per-source
    // kf_process_noise among enabled sources (one shared CV model strength — a per-source q
    // is a possible later refinement; reported as a simplification).
    impl_->smoother_active = false;
    impl_->calib_lag_steps = 0;
    {
        bool   any_kf = false;
        Scalar q_kf   = Scalar(0);
        if (cfg.sensors != nullptr) {
            for (int i = 0; i < cfg.sensor_count; ++i) {
                if (cfg.sensors[i].per_sensor_kf) {
                    any_kf = true;
                    q_kf = std::max(q_kf, cfg.sensors[i].kf_process_noise);
                }
            }
        }
        if (any_kf) {
            // L bakes the NOMINAL tick rate into a STEP count. Both the smoother ring and the
            // calib_ring advance in steps, so they stay mutually aligned at ANY cadence; but the
            // time-lag SEMANTIC (the deeper frontier trails by calib_lag_s seconds) holds only
            // at the nominal rate — off-cadence, L steps != calib_lag_s s (see smoother.hpp).
            int L = static_cast<int>(std::llround(cfg.calib_lag_s * cfg.tick_rate_hz));
            if (q_kf <= Scalar(0)) q_kf = Scalar(1);   // guard a non-positive knob
            const Status ss = impl_->smoother.configure(cfg.max_sources, L, q_kf);
            if (!ok(ss)) return ss;
            impl_->smoother_active = true;
            impl_->calib_lag_steps = impl_->smoother.lag_steps();
        }
    }

    impl_->inited = true;
    return Status::Ok;
}

Status Estimator::add_source(const ISource* src) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    if (src == nullptr)                     return Status::InvalidConfig;
    if (impl_->source_count >= impl_->cfg.max_sources ||
        impl_->source_count >= kMaxSourcesCap) {
        return Status::CapacityExceeded;
    }
    const int slot = impl_->source_count;
    impl_->sources[slot] = src;
    // Bind the prior extrinsic + translation scale + weight for this source from its
    // SensorConfig (defaults: identity / 1.0 / 1.0 when no config matches).
    const SensorConfig* sc = impl_->sensor_for(src->id());
    impl_->prior_extrinsic[slot] = (sc != nullptr) ? sc->prior_extrinsic : SE3{};
    impl_->weight_prior[slot]    = (sc != nullptr) ? sc->weight_prior    : Scalar(1);
    // prior_scale must be positive (validate is expected to enforce this, but be safe:
    // a non-positive scale would blow up the de-scale division in step()).
    const Scalar ps = (sc != nullptr) ? sc->prior_scale : Scalar(1);
    impl_->prior_scale[slot] = (ps > Scalar(0)) ? ps : Scalar(1);
    impl_->prior_time_offset[slot] =
        (sc != nullptr) ? sc->prior_time_offset_s : Scalar(0);
    // Per-sensor smoother gate (Slice 10): this slot feeds SMOOTHED twists to the
    // calibrators iff its SensorConfig::per_sensor_kf is set (else its raw delayed B_corr
    // is used on the deeper path; off entirely if no source enables it).
    impl_->smooth_source[slot] = (sc != nullptr) ? sc->per_sensor_kf : false;
    // Bind the Phase-1 calibrator's per-source prior (histogram basepoint) + scale_calib.
    {
        const SE3  Xp = (sc != nullptr) ? sc->prior_extrinsic : SE3{};
        const bool sccal = (sc != nullptr) ? sc->scale_calib : true;
        impl_->calib1.set_prior(src->id(), Xp, sccal);
        // Phase-2 prior (pinned-ref gauge + translation fallback). yaw/pitch is fed each
        // step from the live Phase-1 estimate; roll/t start from the prior.
        impl_->calib2.set_prior(src->id(), Xp);
    }
    ++impl_->source_count;
    return Status::Ok;
}

Status Estimator::add_correction(const ICorrection* corr) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    if (corr == nullptr)                    return Status::InvalidConfig;
    if (impl_->correction_count >= kMaxCorrectionsCap) return Status::CapacityExceeded;
    impl_->corrections[impl_->correction_count] = corr;
    ++impl_->correction_count;
    return Status::Ok;
}

Status Estimator::step(Timestamp now) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    Impl& s = *impl_;
    const Config& cfg = s.cfg;

    // Causal frontier at t1 = now - delay. We integrate the ACTUAL motion since the
    // last fused frontier: the query/predict interval is [q0, t1] with
    //   q0 = last_t1          (steady state — gap/overlap-free between ticks)
    //   q0 = t1 - window_s    (bootstrap of the first fused step — no prior frontier)
    // and dt = (t1 - q0) in seconds. So `window_s` is now the BOOTSTRAP / LOOKBACK
    // interval for the very first fuse, not a fixed per-tick step: the integrator is
    // correct regardless of tick cadence (DESIGN §7). last_t1 only advances on a
    // successful fuse, so a skipped (uncovered) step never opens an integration gap.
    const Timestamp t1 = now - secs_to_ns(cfg.fusion_delay_s);
    const Timestamp q0 = s.has_last_t1 ? s.last_t1 : (t1 - secs_to_ns(cfg.window_s));

    // Degenerate / non-causal interval guard: nothing to integrate. Lifecycle
    // (Slice 3): before any fuse this is WARMUP; after a fuse it is a graceful
    // DEGRADED (the last good frontier in `result` is retained) — never a crash.
    if (t1 <= q0) {
        s.publish_phase(s.ever_fused_ ? Phase::Degraded : Phase::Warmup);
        s.result.tip_valid = false;
        return Status::NotReady;
    }
    const Scalar dt = static_cast<Scalar>(t1 - q0) / kNanosPerSec;

    // ---- Time-sync sampling (Slice 5, D16) --------------------------------------
    // Feed each source's current ‖ω‖ sample (over a small sub-interval ending at the
    // frontier t1) into the cross-correlator, keyed to BASE/query time, then run a
    // bounded update() once per tick. The ‖ω‖ sub-interval is one common-grid period
    // (sample_dt) so samples land on the time-sync grid. No-op when time-sync is off.
    if (s.timesync_active) {
        const Scalar h_s  = s.timesync.sample_dt();
        const Timestamp h = secs_to_ns(h_s);
        if (h > 0 && t1 - h > 0) {
            for (int i = 0; i < s.source_count; ++i) {
                const ISource* src = s.sources[i];
                if (src == nullptr) continue;
                // ‖ω‖ at the frontier is a BACKWARD DIFFERENCE: the rotation magnitude of
                // the delta over the single sub-interval [t1 - h, t1] (h = one common-grid
                // period), divided by h. It is therefore the mean ‖ω‖ over the window
                // ending at t1, not the instantaneous value at t1 — on a piecewise-constant
                // trajectory a window straddling a turn-rate step reads a ramp, not the
                // step. Harmless for xcorr: BOTH channels see the same backward-difference
                // smoothing, so the time-OFFSET between them is preserved.
                const Expected<Delta> wq = src->query(t1 - h, t1);
                if (!wq.ok()) continue;
                // NOTE: ‖ω‖ is sampled at the source's RAW timeline (no offset applied)
                // — the planted skew between two sources' ‖ω‖ streams is exactly what
                // the cross-correlation measures. Applying the offset here would cancel
                // the very signal we are trying to recover.
                const Scalar wn = omega_norm_from_delta(wq.value().motion, h_s);
                s.timesync.push(src->id(), t1, wn);
            }
        }
        s.timesync.update();
        ++s.ticks_since_sync;
    }

    // Advance the per-source time-offset commit state (N_min + hysteresis, DESIGN §2/§6)
    // once per step, after the vote. effective_offset() and the CalibSnapshot below read
    // the latched flags. No-op (all false) when time-sync is off / for the reference.
    s.update_commit_state();

    // Collect each covered source's frame-aligned delta + weight. The CALIBRATION arrays
    // (calib_*) take EVERY covered source — calibration must see all sources to vote the
    // consensus disagreement (Slice 6/7). The FUSION median arrays (aligned/weights) take
    // only the sources that PARTICIPATE this step per the cold-start switch (D6):
    // MedianFromStart = all (Slice-2 behaviour); ReferenceOnly = the reference + any source
    // whose extrinsic has COMMITTED (a trustworthy prior). `n` indexes the median set;
    // `nc` indexes the calibration set. `med_slot[k]` maps a median entry back to its
    // registered slot for the residual diagnostics.
    int n  = 0;     // participating sources (fusion median)
    int nc = 0;     // covered sources (calibration)
    for (int i = 0; i < s.source_count; ++i) {
        const ISource* src = s.sources[i];
        if (src == nullptr) continue;
        // Apply the source's committed clock offset to its fusion query: shift the
        // interval EARLIER by `off` so the source's internal +off re-shift lands back
        // on true base time [q0, t1] (sign per D21 — removes a planted offset). When
        // time-sync is off / not yet confident this is the configured prior (default 0).
        const Timestamp off = secs_to_ns(s.effective_offset(i));
        const Expected<Delta> q = src->query(q0 - off, t1 - off);

        // Diagnostics (best-effort; coverage flag set regardless of fuse outcome). Reset
        // the per-step fields to neutral so a covered-but-non-participating source (cold
        // start) reads sensible defaults rather than a stale prior-step value; the EMA
        // reliability/bias are re-published below for the participating slots only.
        SourceHealth& h = s.result.health[i];
        h.id          = src->id();
        h.in_window   = q.ok();
        h.weight      = Scalar(0);
        h.residual    = Scalar(0);
        h.reliability = Scalar(1);
        h.bias        = Scalar(0);

        if (!q.ok()) continue;

        const Delta& d = q.value();
        // De-scale the reported translation by the per-source prior_scale BEFORE the
        // frame-align (D20). The source reports B with B.t already multiplied by `scale`
        // (the sim applies scale to B.t AFTER the X-conjugation), so dividing B.t by the
        // matching prior_scale and then conjugating inverts it exactly:
        //   B_corr = { B.R, B.t / prior_scale } ;  A = X o B_corr o X^-1.
        // prior_scale is guaranteed > 0 by add_source (guarded), so this never /0.
        SE3 B_corr = d.motion;
        B_corr.t   = d.motion.t / s.prior_scale[i];
        // Frame-align: A = X o B_corr o X^-1  (X = sensor->base prior extrinsic).
        const SE3& X = s.prior_extrinsic[i];
        const SE3 A  = se3::compose(se3::compose(X, B_corr), se3::inverse(X));

        // Weight = prior x reliability x Sigma-confidence, clamped to [floor, cap]
        // (DESIGN §4, D17). CAUSAL ORDERING: reliability[i] is the value accumulated from
        // PRIOR steps (init 1.0) — this step's residual updates it only AFTER the median is
        // solved (below), so the weight never uses this step's own residual to weight this
        // step's median (that would be circular). The Σ-confidence keeps the Slice-2
        // unit-mixing placeholder; reliability is the NEW quality factor (see sigma_confidence).
        const Scalar sigma_conf = sigma_confidence(d.cov);
        Scalar w = s.weight_prior[i] * s.reliability[i] * sigma_conf;
        w = std::max(cfg.weight_floor, std::min(cfg.weight_cap, w));

        // Capture the de-scaled reported sensor-frame delta for the calibrators
        // (direction -> yaw/pitch + magnitude ratio -> scale; D11/D20). EVERY covered
        // source feeds calibration. Also carry the raw Σ-confidence (D5 vote_weight).
        s.calib_ids[nc]      = src->id();
        s.calib_slot[nc]     = i;            // registered slot (Slice-10 smoother lookup)
        s.calib_reported[nc] = B_corr;
        s.calib_conf[nc]     = sigma_conf;
        ++nc;

        // Fusion median: only participating sources (cold-start switch, D6).
        if (s.participates(i)) {
            s.aligned[n]  = A;
            s.weights[n]  = w;
            s.med_slot[n] = i;
            h.weight      = w;
            ++n;
        }
    }

    s.result.source_count = s.source_count;

    // Lifecycle (Slice 3, D23): need at least one PARTICIPATING source covering the
    // window to fuse. Total source loss this step is a graceful downgrade, not a
    // crash: before any fuse it is WARMUP; after a fuse it is DEGRADED + NotReady
    // (the last good frontier stays in `result`; no new frontier is published).
    // Under ReferenceOnly cold-start before any commit, n is the reference alone
    // (reference-only dead-reckon) — that is n == 1, handled by the fuse path below.
    if (n == 0) {
        s.publish_phase(s.ever_fused_ ? Phase::Degraded : Phase::Warmup);
        s.result.tip_valid = false;
        return Status::NotReady;
    }

    // Robust consensus delta over the window.
    median::Params mp;
    mp.max_iters = cfg.weiszfeld_max_iters;
    mp.tol       = cfg.weiszfeld_tol;
    mp.eps       = cfg.weiszfeld_eps;
    mp.lambda    = cfg.metric_lambda;
    const median::Result med = median::solve(s.aligned, s.weights, n, mp);

    // Per-source residual to the consensus (diagnostics). Indexed by the median set via
    // med_slot (cold-start may exclude covered-but-non-participating sources from the
    // median); a covered source not in the median keeps its residual at 0.
    //
    // VARIANCE-EMA RELIABILITY (Slice 9, D17). Each fused step, attribute the consensus
    // distance d to the participating slot and accrue an EMA mean (bias) + EMA variance
    // (zero-mean scatter) of d. Then recompute the reliability multiplier from a robust
    // baseline (the median of the warmed-up participants' variances). The update is
    // strictly AFTER the median + weights are formed (causal — this step's residual never
    // weights this step's median). Bounded by n (<= max_sources); no heap.
    const Scalar alpha = cfg.reliability_ema_alpha;
    for (int k = 0; k < n; ++k) {
        const int i = s.med_slot[k];
        SourceHealth& h = s.result.health[i];
        const Scalar d = se3::split_distance(med.value, s.aligned[k], cfg.metric_lambda);
        h.residual = d;

        // West's incremental EMA mean/variance of the per-source residual. resid_mean is the
        // BIAS diagnostic: it is the EMA mean of d (= split_distance), an UNSIGNED residual
        // MAGNITUDE (>= 0) — "mean residual distance to consensus", NOT a signed per-DOF
        // offset (a source biased high or low both read positive; even a zero-mean noisy
        // source accrues positive bias). It is left for the calibrator to absorb via the
        // existing Slice-6/7/8 observe path and is exposed for DIAGNOSTICS ONLY, NOT folded
        // into the weight (D17 "bias → calibrator, not weight"). resid_var is the zero-mean
        // SCATTER around the running mean — a biased-but-consistent source has a large mean
        // but a SMALL variance, so it is NOT penalized; resid_var (not the mean) is what
        // distinguishes noise from a systematic offset (the D17 split).
        if (s.resid_n[i] == 0) {
            s.resid_mean[i] = d;
            s.resid_var[i]  = Scalar(0);
        } else {
            const Scalar delta = d - s.resid_mean[i];
            s.resid_mean[i] += alpha * delta;
            s.resid_var[i]  += alpha * (delta * delta - s.resid_var[i]);
        }
        ++s.resid_n[i];
    }

    // Robust baseline ref_var = the MEDIAN of resid_var[j] over THIS step's median
    // participants that have warmed up (resid_n >= kRelWarmup). Copy the qualifying values
    // into the fixed stack scratch and median-select in place (bounded, heap-free). If
    // fewer than 2 qualify, SKIP the recompute (hold prior reliability) — warmup / low-data
    // never thrashes the weights.
    int q = 0;
    for (int k = 0; k < n; ++k) {
        const int i = s.med_slot[k];
        if (s.resid_n[i] >= kRelWarmup) s.rel_scratch[q++] = s.resid_var[i];
    }
    if (q >= 2) {
        // Median of the q values (q <= max_sources <= kMaxSourcesCap). std::nth_element on
        // a fixed stack array is heap-free; for even q take the lower-middle (deterministic,
        // a reference variance only needs to be representative, not the exact midpoint).
        const int mid = q / 2;
        std::nth_element(s.rel_scratch, s.rel_scratch + mid, s.rel_scratch + q);
        // Floor BOTH the baseline and each source's variance by the same kRelVarEps so the
        // ratio is well-conditioned. Critically this makes the DEGENERATE near-zero case
        // (a noiseless oracle with correct priors -> all resid_var ~ machine-zero) read
        // ratio ~ eps/eps = 1 (reliability ~ 1) instead of 0/0 noise that would collapse
        // every clean source to the floor. Equal-noise sources -> equal resid_var ->
        // ratio ~ 1; a NOISIER source -> larger resid_var -> ratio < 1 (reliability floored);
        // a BIASED-but-low-noise source -> small resid_var (its bias is in resid_mean, NOT
        // here) -> ratio ~ 1 (NOT collapsed) — the D17 split.
        const Scalar ref_var = std::max(s.rel_scratch[mid], kRelVarEps);
        for (int k = 0; k < n; ++k) {
            const int i = s.med_slot[k];
            if (s.resid_n[i] < kRelWarmup) continue;   // pre-warmup keeps reliability 1.0
            const Scalar ratio = ref_var / std::max(s.resid_var[i], kRelVarEps);
            s.reliability[i] = std::max(cfg.reliability_floor,
                                        std::min(cfg.reliability_cap, ratio));
        }
    }
    // Surface the reliability + bias diagnostics for every participating slot (non-
    // participants keep the neutral 1.0 / 0.0 the result struct defaults to).
    for (int k = 0; k < n; ++k) {
        const int i = s.med_slot[k];
        SourceHealth& h = s.result.health[i];
        h.reliability = s.reliability[i];
        h.bias        = s.resid_mean[i];
    }

    // ---- Calibration stage (Phase 1 Slice 6 + Phase 2 Slice 7; D5/D10/D11/D20/D21) ----
    // The fused body angular SPEED is log(med.R)/dt (rad/s); the fused translation is the
    // per-step DISPLACEMENT med.t (m, NOT ÷dt). NOTE the gate operands carry different
    // time-normalization: straight_omega_max gates a speed, straight_trans_min a per-step
    // displacement — so straight_trans_min is CADENCE-DEPENDENT (tune to the tick rate; see
    // calibration.hpp). Phase 1 self-gates on the straight regime, Phase 2 on the turn
    // regime. NO calibration feedback into fusion here (Slice 8 runs after publish).
    const Vec3 fused_omega = so3::log(med.value.R) / dt;
    const Vec3 fused_trans = med.value.t;

    if (!s.smoother_active) {
        // --- CAUSAL frontier (Slice 6/7/8 — the DEFAULT, byte-identical path) ----------
        // No source enables per_sensor_kf: drive the calibrators IMMEDIATELY with THIS
        // step's consensus + per-source de-scaled deltas, exactly as before Slice 10. This
        // is the hard non-regression guard — the observe() sequence is unchanged.
        s.run_calibration_observe(nc, s.calib_ids, s.calib_reported, fused_omega,
                                  fused_trans, med.value, s.calib_conf);
    } else {
        // --- DEEPER (lag-L) frontier (Slice 10, D18) -----------------------------------
        // Calibration tolerates latency, so it runs at now − delay − L. Push each
        // smoothing-enabled source's MEASURED body twist log(B_corr)/dt into the smoother
        // (it maintains its own internal two-sided lag-L window), and buffer THIS step's
        // full calibration-input set into the depth-(L+1) ring. Once the ring is full, the
        // OLDEST (lag-L) entry drives the calibrators — with each smoothing-enabled source's
        // delta REPLACED by exp(smoothed_twist · dt_old) (the de-scale convention preserved:
        // the smoothed twist is itself the twist of the de-scaled B_corr) and the consensus
        // taken from that SAME lag-L entry (so the source-vs-consensus vote stays
        // time-aligned). Non-smoothed sources use their (delayed) raw B_corr.
        const int depth = s.calib_lag_steps + 1;

        // (a) Push the measured twist of every smoothing-enabled covered source. Advance the
        // per-slot MONOTONIC push count (Slice-10 review MAJOR) so the frame below can stamp
        // the exact push sequence number this source reached at this wall-clock step.
        for (int i = 0; i < nc; ++i) {
            const int slot = s.calib_slot[i];
            if (slot < 0 || !s.smooth_source[slot]) continue;
            const Vec6 zt = se3::log(s.calib_reported[i]) / dt;   // body twist of de-scaled B
            // Slot is pre-validated (slot in [0, max_sources) and smooth_source set at
            // add_source) so push() cannot fail here; ignore the Status (strict-core ethos).
            (void)s.smoother.push(slot, zt, dt);
            ++s.push_seq[slot];
        }

        // (b) Append THIS step's calibration inputs to the deeper ring (drop oldest if
        // full). Write DIRECTLY into the destination ring slot (no large stack temporary —
        // strict core). This is a LOGICAL SHIFT, not an O(1) circular buffer: the full-ring
        // case copies every entry down by one (O(L) per step, depth <= kCalibRing — bounded,
        // heap-free), matching the smoother's own ring (smoother.cpp push()).
        if (s.calib_ring_count >= depth) {
            for (int k = 0; k < depth - 1; ++k) s.calib_ring[k] = s.calib_ring[k + 1];
        }
        const int dst = (s.calib_ring_count < depth) ? s.calib_ring_count : (depth - 1);
        Impl::CalibFrame& frame = s.calib_ring[dst];
        frame.nc           = nc;
        frame.fused_omega  = fused_omega;
        frame.fused_trans  = fused_trans;
        frame.fused_motion = med.value;
        frame.dt           = dt;
        for (int i = 0; i < nc; ++i) {
            frame.ids[i]      = s.calib_ids[i];
            frame.reported[i] = s.calib_reported[i];
            frame.conf[i]     = s.calib_conf[i];
            frame.slot[i]     = s.calib_slot[i];
            // Stamp the post-push sequence number for smoothing sources (Slice-10 review
            // MAJOR); -1 marks a non-smoothing / un-pushed entry (never substituted).
            const int slot    = s.calib_slot[i];
            frame.push_at[i]  = (slot >= 0 && s.smooth_source[slot]) ? s.push_seq[slot]
                                                                     : static_cast<long>(-1);
        }
        if (s.calib_ring_count < depth) ++s.calib_ring_count;

        // (c) Once the ring is full, drive the calibrators with the lag-L-OLD entry. The
        // smoother's emitted twist is the two-sided estimate of THAT same lag-L sample
        // (both trail the live frontier by exactly L steps), so they are time-aligned.
        if (s.calib_ring_count >= depth) {
            const Impl::CalibFrame& old = s.calib_ring[0];
            for (int i = 0; i < old.nc; ++i) {
                const int slot = old.slot[i];
                // Substitute the smoothed twist ONLY when it is TIME-ALIGNED with this lag-L
                // frame (Slice-10 review MAJOR). The smoother's emitted sample after N pushes
                // is the (N-L)-th push; it lines up with THIS frame iff the source has been
                // pushed exactly L times since the frame was created, i.e.
                //   current push_seq[slot] - old.push_at[i] == L.
                // Under per-source dropout the source missed pushes inside (t-L, t], so the
                // equality fails and the emitted twist belongs to an EARLIER wall-clock step
                // than this frame's consensus — feeding it would inject a phase-error bias on
                // a turning trajectory. We then fall back to the raw delayed delta (a clean
                // variance loss instead of a calibration bias). No-dropout: equality always
                // holds, so this path is byte-identical to the pre-fix substitution.
                const bool aligned =
                    slot >= 0 && s.smooth_source[slot] && s.smoother.ready(slot) &&
                    old.push_at[i] >= 0 &&
                    (s.push_seq[slot] - old.push_at[i]) == static_cast<long>(s.calib_lag_steps);
                if (aligned) {
                    // Reconstruct the de-scaled SE3 from the smoothed body twist over the
                    // lag-L entry's dt: B_corr_smoothed = exp(w_smoothed · dt_old). This
                    // preserves the de-scale convention (the twist IS that of the de-scaled
                    // B_corr), so the calibrators' frame-align/magnitude logic is unchanged.
                    const Vec6 w = s.smoother.smoothed(slot);
                    s.calib_smoothed[i] = se3::exp(w * old.dt);
                } else {
                    s.calib_smoothed[i] = old.reported[i];   // raw delayed delta
                }
            }
            s.run_calibration_observe(old.nc, old.ids, s.calib_smoothed, old.fused_omega,
                                      old.fused_trans, old.fused_motion, old.conf);
        }
        // Before the ring fills, calibration simply has not started yet (the deeper frontier
        // is not reachable until L steps have elapsed) — a no-op, documented.
    }

    // Calibration snapshot (Slice 5 fills the time-offset DOF; extrinsic/scale stay at
    // their priors until Slices 6-8). Per source: the offset estimate currently APPLIED
    // to its query (committed time-sync value, else the prior) plus the time-sync
    // histogram confidence. `committed` marks that the estimate cleared the N_min +
    // hysteresis commit gate (so it is actually driving fusion rather than the prior) —
    // it is the latched per-source state advanced in update_commit_state() above.
    // Phase-1 calibration (Slice 6) now fills the extrinsic yaw/pitch + scale DOF. The
    // reported extrinsic is the prior with the recovered yaw/pitch correction (roll +
    // translation stay at the prior — Slice 7). `scale` is the residual scale vs the prior
    // (calib1.scale ~ 1 until a straight regime is observed) multiplied back by the prior
    // (the calibrator votes the ratio off the DE-SCALED B, so the absolute scale =
    // prior_scale * residual). `confidence` is the so(3) direction concentration (yaw/pitch
    // reliability). The time-offset DOF + its commit flag are unchanged (Slice 5). Slice 6
    // does NOT feed any of this back into fusion (Slice 8).
    for (int i = 0; i < s.source_count; ++i) {
        const SourceId id = s.sources[i]->id();
        CalibSnapshot& cs = s.result.calib[i];
        cs.id            = id;
        // Full extrinsic: Phase-1 yaw/pitch ∘ Phase-2 roll (rotation) + Phase-2 lever arm
        // (translation). calib2.extrinsic() composes the yaw/pitch rotation fed in above
        // with the recovered roll, and uses the LS-mode lever arm (the prior until a turn
        // regime is observed — so before any turning this equals the Phase-1 extrinsic).
        cs.extrinsic     = s.calib2.extrinsic(id);
        cs.scale         = s.prior_scale[i] * s.calib1.scale(id);
        cs.time_offset_s = s.effective_offset(i);
        // confidence stays the TIME-OFFSET confidence (Slice 5); per-DOF confidences are
        // reported in their own fields (each DOF converges in its own regime).
        cs.confidence    = s.timesync_active ? s.timesync.confidence(id) : Scalar(0);
        // Rotation confidence: Phase-1 yaw/pitch (so(3)) AND, once Phase 2 has roll votes,
        // the roll (S¹) concentration — combined as the MIN (the weakest rotational DOF
        // bounds the joint reliability). Before any roll votes, fall back to Phase 1 alone.
        const Scalar yp_conf   = s.calib1.extrinsic_confidence(id);
        const Scalar roll_conf = s.calib2.extrinsic_confidence(id);
        cs.extrinsic_confidence = (s.calib2.roll_vote_count(id) > Scalar(0))
                                      ? std::min(yp_conf, roll_conf)
                                      : yp_conf;
        cs.scale_confidence       = s.calib1.scale_confidence(id);
        cs.translation_confidence = s.calib2.translation_confidence(id);
        // Per-DOF commit flags (Slice 8): the latched state from the PREVIOUS step's
        // apply_calib_feedback() — i.e. which DOF is currently DRIVING this step's fusion
        // prior. Each underlying per-DOF latch is monotone (commit then hysteresis-hold), so
        // each flag is non-thrashing. `committed` is the time-offset commit (Slice 5);
        // `extrinsic_committed` is the yaw/pitch (so(3) rotation) commit; the roll commit is
        // a separate refinement on top (it does not gate this flag — combining them would
        // make the flag oscillate during the window where roll is observed but not yet
        // committed). `translation_committed` is the xyz lever-arm commit; `scale_committed`
        // the per-source scale commit.
        cs.committed             = s.offset_committed[i];
        cs.extrinsic_committed   = s.ext_committed[i];
        cs.scale_committed       = s.scale_committed[i];
        cs.translation_committed = s.lever_committed[i];
    }

    // Adaptive process noise from the inter-source spread (DESIGN §4, D4).
    // q_scale / q_floor come from Config (CONFIG §3). Adaptive: floor + q_scale*spread^2.
    // Non-adaptive: just the per-axis floor (no spread term).
    Mat6 q_pose;
    if (cfg.adaptive_q) {
        q_pose = Eskf::adaptive_q(med.spread, cfg.q_scale, cfg.q_floor);
    } else {
        q_pose = Eskf::adaptive_q(Scalar(0), cfg.q_scale, cfg.q_floor);
    }

    // ESKF predict on the median delta over the integrated interval [q0, t1] (dt is
    // the actual elapsed time, NOT a fixed window). Anchor the odom frame at the
    // first fused tick (pose starts at identity — gauge anchored at init, DESIGN §7).
    if (!s.eskf_started) {
        s.eskf.init(SE3{}, Mat12::Identity());
        s.eskf_started = true;
    }
    s.eskf.predict(med.value, dt, q_pose);

    // This fuse succeeded: advance the integration frontier so the next step picks
    // up exactly where this one ended (no gap, no overlap).
    s.last_t1     = t1;
    s.has_last_t1 = true;

    // ---- Absolute-reference corrections (Slice 11, D22) -------------------------
    // AFTER predict, BEFORE publishing the frontier/tip: for each registered ICorrection,
    // evaluate its measurement at the CURRENT (post-predict) frontier state, then run the
    // Mahalanobis-gated ESKF update. Accepted updates SHRINK P toward the measurement and
    // correct the fused pose (e.g. a GPS-like position fix removing odometry drift); the
    // published frontier + the predicted tip are re-read from the CORRECTED state below.
    // No-op (byte-identical to the predict-only path) when no correction is registered.
    //
    // Diagnostics: count availability (evaluate true), acceptance vs gate-rejection, and the
    // last NIS (computable now that there is an innovation — closes the Slice-14 NIS
    // deferral). The evaluate() State carries the real frontier stamp t1 (the absolute ref
    // needs it to sample the reference at the frontier), so we stamp a working copy first.
    CorrectionDiag cdiag;
    if (s.correction_count > 0) {
        for (int c = 0; c < s.correction_count; ++c) {
            const ICorrection* corr = s.corrections[c];
            if (corr == nullptr) continue;
            State x = s.eskf.state();      // post-predict (+ any earlier accepted update)
            x.stamp = t1;                  // give the plugin the real frontier time
            Measurement m;
            if (!corr->evaluate(x, m)) continue;   // no fix available this step
            ++cdiag.corr_evaluated;
            const bool applied = s.eskf.update(m, cfg.mahalanobis_chi2);
            if (applied) ++cdiag.corr_applied; else ++cdiag.corr_rejected;
            cdiag.last_nis = s.eskf.last_nis();     // last NIS, accepted or rejected
        }
    }
    s.result.correction = cdiag;

    // Frontier state. Drive the published stamp from the real frontier t1. Re-read AFTER
    // the corrections so the published pose/covariance reflect any accepted update.
    s.result.frontier        = s.eskf.state();
    s.result.frontier.stamp  = t1;

    // Lifecycle (Slice 3, D23): a successful fuse latches ever_fused_; the phase is
    // NOMINAL when >= min_sources_warn sources PARTICIPATED, else DEGRADED (reduced
    // outlier-rejection margin — e.g. ReferenceOnly dead-reckon at n == 1, or a
    // partial source loss below the warn threshold). Recomputed every step from `n`,
    // so a drop below min_sources_warn downgrades NOMINAL->DEGRADED automatically.
    s.ever_fused_ = true;
    s.publish_phase((n >= cfg.min_sources_warn) ? Phase::Nominal : Phase::Degraded);

    // Predicted tip: const-velocity extrapolation from the frontier to `now`.
    if (cfg.emit_predicted_tip) {
        const Scalar dt_ahead = cfg.fusion_delay_s;   // frontier (t1) -> now
        s.eskf.predict_tip(dt_ahead, cfg.tip_cov_inflation, s.result.tip);
        s.result.tip.stamp = now;
        s.result.tip_valid = true;
    } else {
        s.result.tip_valid = false;
    }

    // ---- Calibration -> fusion feedback (Slice 8, DESIGN §6) --------------------
    // ATOMIC between-step swap: run AFTER the result is published so this step saw a
    // consistent prior set; the per-DOF commit + value swap (+ the scale re-anchor) take
    // effect on the NEXT step's fuse. This closes the calibration->fusion loop — a committed
    // scale / time-offset / extrinsic is swapped into the value fusion frame-aligns/de-scales
    // with, improving the fused trajectory. No-op while no DOF has cleared the commit gate (so
    // calibration-off / unconfident == the prior-driven Slice-6/7 behaviour).
    s.apply_calib_feedback();

    return Status::Ok;
}

const Result& Estimator::latest() const { return impl_->result; }

Expected<int> Estimator::serialize(unsigned char* buf, int cap) const {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)buf; (void)cap;
    return Status::NoData;   // TODO Slice 12
}

Status Estimator::deserialize(const unsigned char* buf, int len) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)buf; (void)len;
    return Status::NoData;   // TODO Slice 12 (incl. config-hash guard)
}

// --- validate() lives with the config it checks --------------------------
Status validate(const Config& cfg) {
    if (cfg.max_sources < 1 || cfg.max_sources > 32) return Status::OutOfRange;
    if (cfg.buffer_capacity < 16)                    return Status::OutOfRange;
    if (cfg.sensor_count > cfg.max_sources)          return Status::CapacityExceeded;
    if (cfg.tick_rate_hz <= 0.0)                     return Status::OutOfRange;
    if (cfg.window_s <= 0.0)                         return Status::OutOfRange;
    if (cfg.weiszfeld_max_iters < 1)                 return Status::OutOfRange;
    if (cfg.metric_lambda <= 0.0)                    return Status::OutOfRange;
    if (cfg.commit_drop >= cfg.commit_concentration) return Status::InvalidConfig;
    if (cfg.reference_sensor_id >= cfg.max_sources)  return Status::OutOfRange;
    // Lifecycle (Slice 3): lower-bound only. NOT tied to max_sources — a
    // min_sources_warn > max_sources is a legitimate "never NOMINAL" config.
    if (cfg.min_sources_warn < 1)                    return Status::OutOfRange;

    // Median / fusion knob ranges (CONFIG §§1-4).
    if (cfg.confidence_blend < 0.0 || cfg.confidence_blend > 1.0)
        return Status::OutOfRange;                       // blend factor in [0, 1]
    if (cfg.weiszfeld_tol <= 0.0 || cfg.weiszfeld_tol >= 1.0)
        return Status::OutOfRange;                       // tol in (0, 1)
    if (cfg.weiszfeld_eps <= 0.0 || cfg.weiszfeld_eps >= 1.0)
        return Status::OutOfRange;                       // eps in (0, 1)
    if (cfg.fusion_delay_s < 0.0 || cfg.fusion_delay_s > 2.0)
        return Status::OutOfRange;                       // delay in [0, 2] s
    if (cfg.weight_floor <= 0.0 || cfg.weight_floor >= 1.0)
        return Status::OutOfRange;                       // floor in (0, 1)
    if (cfg.weight_cap < 1.0)
        return Status::OutOfRange;                       // cap >= 1
    // Variance-EMA reliability clamp (Slice 9, D17): floor in (0, 1]; cap >= 1.
    if (cfg.reliability_floor <= 0.0 || cfg.reliability_floor > 1.0)
        return Status::OutOfRange;                       // reliability floor in (0, 1]
    if (cfg.reliability_cap < 1.0)
        return Status::OutOfRange;                       // reliability cap >= 1
    if (cfg.tip_cov_inflation < 1.0)
        return Status::OutOfRange;                       // inflation >= 1
    if (cfg.q_scale <= 0.0)
        return Status::OutOfRange;                       // q_scale > 0
    for (int i = 0; i < 6; ++i) {
        if (cfg.q_floor[i] < 0.0) return Status::OutOfRange;   // each q_floor[i] >= 0
    }

    // Time-sync knob ranges (CONFIG §5). Enforced here (config-standalone) as well as
    // inside TimeSync::configure() — so a caller validating config with timesync OFF
    // still catches them. match_metric is an enum (no range check). offset_hist is a
    // HistogramConfig validated by Histogram1D::configure() at TimeSync setup.
    if (cfg.excitation_min_var < 0.0)
        return Status::OutOfRange;                        // ‖ω‖ variance gate >= 0
    if (cfg.max_lag_s <= 0.0 || cfg.max_lag_s > 2.0)
        return Status::OutOfRange;                        // max_lag_s in (0, 2]

    // Absolute-ref Mahalanobis gate (Slice 11, D22): the chi-square threshold must be
    // positive — a non-positive threshold would reject EVERY measurement (NIS >= 0 always),
    // silently disabling the correction path.
    if (cfg.mahalanobis_chi2 <= 0.0)                     return Status::OutOfRange;

    // TODO: per-sensor + histogram range checks.
    return Status::Ok;
}

} // namespace ofc

