// ofc/core/calibration.hpp — online extrinsic/scale self-calibration (Slice 6: Phase 1).
//
// PHASE 1 — STRAIGHT-REGIME yaw/pitch + per-source scale (DESIGN §3/§6, DECISIONS
// D5/D11/D20). In straight motion the base travels along its forward axis (±e_x);
// a source reports that direction rotated into its own frame, which reveals 2
// rotational DOF (yaw, pitch) of its extrinsic. Roll about the forward axis is
// UNOBSERVABLE here (Slice 7, turning). Pure translation also means a source's
// translation-magnitude ratio to the pinned reference is its scale with no
// lever-arm confound. Both are observable ONLY in straight motion — hence the gate.
//
// THE MAPPING (direction -> so(3) @ prior -> mode -> forward axis -> yaw,pitch).
// The sim/oracle pins it (sim/synthetic_source.cpp): a source reports
//     B = X_true^{-1} o A o X_true ,  then  B.t *= scale ,
// where A is the TRUE base delta and X = prior_extrinsic is sensor->base. In
// straight forward motion A.R = I and A.t = d * e_x_base (d>0 fwd, d<0 reverse), so
//     B.t = scale * R_Xtrue^{-1} * A.t   =>   dir_B = B.t / ||B.t|| = R_Xtrue^{-1} * (±e_x).
// We map the source's observed forward direction into the BASE frame through its
// PRIOR extrinsic:  g_obs = R_Xprior * dir_B  (≈ e_x when the mount matches the prior).
// REVERSE-FOLD: when the CONSENSUS (fused) motion is reverse (sign of fused_trans on the
// base-forward axis is negative) negate g_obs, so forward and reverse straight segments both
// land in the SAME hemisphere (unimodal). Folding by the consensus sign — not the source's
// own g_obs.x — keeps this correct for a sideways / far-off-prior mount (whose g_obs.x is
// noise). A folded g_obs that is still ≥90° off +e_x is outside Phase-1's small-deviation
// regime and is SKIPPED (it would otherwise need a near-π so(3) log, which is singular).
// With the fold OFF the reverse samples are not negated; staying backward they trip that
// same ≥90° guard and are DROPPED — so fold-OFF simply discards the reverse half (it does
// NOT deposit a 180°-off second peak; in the non-circular so(3) histogram such a near-π log
// would only have edge-clamped into the boundary bins anyway, never a true bin at π).
// The candidate rotation is the MINIMAL rotation in the base frame taking
// e_x -> g_obs; its so(3) log φ = log(δR) is voted into THREE Histogram1D channels
// (basepoint = the prior, i.e. φ ≈ 0 when the mount matches the prior — data sits near
// the basepoint, no pole/antipode, isotropic resolution; per D11 store all 3 channels,
// φ_x is ~0 only up to noise). Estimate = per-channel mode -> φ_mode -> δR = exp(φ_mode)
// -> recovered forward axis f = δR * e_x -> yaw = atan2(f_y, f_x),
// pitch = atan2(-f_z, hypot(f_x, f_y)).
// CONTRACTIVE EXTRINSIC (Slice-8 re-anchor). δR records the forward minimal rotation
// (δR·e_x = g_obs), so the forward axis / yaw / pitch read δR·e_x directly. The reported
// extrinsic rotation, however, is the INVERSE composition  R = δRᵀ · R_Xprior  — chosen so
// the recovered map satisfies  R·dir_B = δRᵀ·g_obs = e_x  EXACTLY (not just at the fixed
// point). That makes re-anchoring the basepoint to a freshly committed extrinsic a
// CONTRACTIVE iterate: from a LARGE prior error the next round's g_obs lands on e_x, δR→I,
// and the estimate converges to the planted forward axis. At the fixed point (prior ==
// truth) g_obs = e_x ⇒ δR = I ⇒ R = R_Xprior. Roll and the translation are kept at the
// prior (roll is Slice 7, xyz lever-arm is Slice 7).
//
// SCALE (D20). With both the source and the reference seeing the same base
// translation magnitude over the window, ||B_i.t|| / ||B_ref.t|| = scale_i / scale_ref.
// The pinned reference is the gauge (its scale is held at its prior, default 1), so the
// magnitude ratio vs the reference's reported magnitude is voted into a per-source scale
// Histogram1D; estimate = mode. A source with scale_calib == false is pinned at its prior.
//
// CONFIDENCE = histogram concentration (the Slice-4 peak-concentration), exposed per
// source. Slice 6 ONLY estimates + exposes (into Result::CalibSnapshot). Feeding these
// back to correct fusion is Slice 8 — NOT wired here.
//
// STRICT CORE: no heap after configure(); fixed-capacity per-source histograms sized at
// configure(); bounded loops; status-code returns; no exceptions; double math. Reuses
// Histogram1D storage (3 so(3) channels + 1 scale per source).
#ifndef OFC_CORE_CALIBRATION_HPP
#define OFC_CORE_CALIBRATION_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/histogram.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

namespace ofc {

// Phase-1 (straight-regime) calibrator: per-source yaw/pitch + scale from odometry
// disagreement. Configure once from a Config + the pinned reference id, then feed each
// fused STRAIGHT-regime step to observe(); read per-source estimates + confidence.
class Phase1Calibrator {
public:
    // Compile-time maxima (strict core, no heap). Matches the Result arrays / Config cap.
    static constexpr int kMaxSources = 32;

    Phase1Calibrator() = default;

    // Bind the active configuration + the pinned reference source id. Sizes the per-source
    // histograms (3 so(3) channels from cfg.so3_hist + 1 scale from cfg.scale_hist),
    // captures the straight gate thresholds / reverse-fold / ref-cross-check flags, and
    // clears all state. Validates:
    //   reference_id < kMaxSources                 -> else OutOfRange
    //   so3_hist / scale_hist pass Histogram1D::configure -> propagated
    // No heap occurs here nor in any later call.
    Status configure(const Config& cfg, SourceId reference_id);

    // Drop all votes (keeps the configuration). Estimates fall back to "no data".
    void reset();

    // Register the per-source PRIOR extrinsic (sensor->base) — the histogram basepoint
    // and the rotation the recovered correction is applied on top of. Also captures
    // whether the source calibrates scale (scale_calib). Call once per source after
    // configure(), before observe(). A source never registered uses identity / scale_calib
    // = true defaults. Returns OutOfRange if id is out of range, CapacityExceeded at cap.
    Status set_prior(SourceId id, const SE3& prior_extrinsic, bool scale_calib = true);

    // --- Re-anchor API (Slice 8 feedback loop) -------------------------------
    // RE-ANCHOR the so(3) basepoint to a NEW committed extrinsic (sensor->base). The
    // direction votes are cast as the so(3) log of the minimal rotation e_x -> g_obs,
    // where g_obs = R_basepoint * dir_B; so when fusion swaps a freshly-committed
    // extrinsic into the source's prior, the OLD votes (cast @ the old basepoint) are
    // STALE. Calling set_basepoint(id, X) moves the basepoint to X.R so subsequent votes
    // refine the RESIDUAL direction around the new basepoint (φ ≈ 0). The caller should
    // also reset_so3(id) so the stale votes do not pin the now-residual mode. Heap-free;
    // OutOfRange / CapacityExceeded as set_prior. extrinsic() = δR_modeᵀ * X.R thereafter
    // (the inverse minimal rotation — the contractive composition; see THE MAPPING above).
    Status set_basepoint(SourceId id, const SE3& basepoint);
    // Drop only the so(3) (yaw/pitch) histogram votes for `id` (keeps scale + basepoint).
    // The post-reset estimate falls back to the basepoint until new votes land. Used on a
    // yaw/pitch commit+re-anchor so the histogram re-concentrates on the residual.
    void reset_so3(SourceId id);
    // Drop only the scale histogram votes for `id` (keeps so(3) + basepoint). Used on a
    // scale commit+re-anchor (after fusion's prior_scale absorbs the committed scale, the
    // residual ratio re-votes ≈ 1). The post-reset scale() falls back to 1 until new votes.
    void reset_scale(SourceId id);

    // Consume one FUSED step. `fused_omega` is the consensus body angular RATE over the
    // step (rad/s, ÷dt) used by the straight gate's rotation test. `fused_trans` is the
    // consensus translation DISPLACEMENT over the step (m, NOT a rate) — it is both the
    // gate's translation test operand and the source of the consensus forward/reverse sign
    // (its projection on the base-forward axis ±e_x) used by the reverse-fold. Per
    // registered source i in [0, n): `ids[i]` is its id and `reported[i]` is its DE-SCALED
    // reported sensor-frame delta B_corr (B with the translation divided by prior_scale,
    // matching the estimator's pre-median de-scale — so the magnitude ratio recovers the
    // RESIDUAL scale error vs the prior). The reference source must be present among `ids`
    // for a scale vote to occur.
    //
    // GATE: votes are deposited ONLY when the fused motion is straight —
    //   ||fused_omega|| < straight_omega_max  AND  ||fused_trans|| > straight_trans_min
    // NOTE the two gate operands have DIFFERENT time-normalization: straight_omega_max is a
    // speed (rad/s) vs ||fused_omega||, while straight_trans_min is a per-step DISPLACEMENT
    // (m) vs ||fused_trans|| — so straight_trans_min is CADENCE-DEPENDENT (halving dt halves
    // the per-step displacement at fixed speed). Tune it for the deployment tick rate.
    // — and, when ref_cross_check is on, the reference source's reported delta also reads
    // "straight forward/back" (small rotation + sizable translation). Off-gate steps are a
    // no-op (the observability spine: yaw/pitch/scale are unobservable unless straight).
    //
    // REVERSE-FOLD: a far-off-prior / sideways mount's per-source g_obs.x is noise, so the
    // fold uses the CONSENSUS forward/reverse sign (from fused_trans), not g_obs.x — forward
    // and reverse straight segments collapse to ONE peak for any mount. A folded direction
    // ≥90° off +e_x (outside Phase-1's small-deviation regime) is SKIPPED, never voted as a
    // near-π so(3) log (which would be singular / NaN). With reverse_fold == false the
    // reverse samples are NOT folded; staying backward they hit that same ≥90° guard and are
    // dropped (fold-OFF discards the reverse half — no edge-clamped boundary-bin second peak).
    //
    // VOTE WEIGHT (cfg.vote_weight, D5 "votes weighted by ‖ω‖"). Each histogram vote's
    // mass is scaled by a per-source factor (see vote_weight_factor()): One -> 1;
    // Confidence -> the source's Σ-confidence (the same inverse-covariance scalar fusion
    // weights with, supplied per source in `confidences[i]`); Rotation -> the turn
    // magnitude ‖fused_omega‖ (floored to a small positive so straight-regime Phase-1
    // votes — where ‖ω‖≈0 by gate — still land); Combo -> the product of the Rotation and
    // Confidence factors. `confidences` is OPTIONAL (nullptr -> unit confidence, so
    // Confidence/Combo collapse to the Rotation factor); when supplied it must have n
    // entries aligned with `ids`. So(3) (3 channels) and scale share the same per-source
    // factor.
    //
    // Returns Ok when the step was voted, NotReady when gated out (not an error),
    // NotInitialized if unconfigured.
    Status observe(int n, const SourceId* ids, const SE3* reported,
                   const Vec3& fused_omega, const Vec3& fused_trans,
                   const Scalar* confidences = nullptr);

    // --- Per-source readouts (committed = histogram mode) --------------------
    // Recovered forward axis in the BASE frame (unit vector). Returns the prior forward
    // (R_Xprior * e_x) for an unvoted / unknown source.
    Vec3   forward_axis(SourceId id) const;
    Scalar yaw(SourceId id)   const;        // atan2(f_y, f_x)   [rad]
    Scalar pitch(SourceId id) const;        // atan2(-f_z, hypot(f_x,f_y))  [rad]
    // Recovered extrinsic: rotation = δR_modeᵀ * R_Xprior (the INVERSE minimal rotation, so
    // R·dir_B = e_x exactly and re-anchoring the so(3) BASEPOINT to it is CONTRACTIVE — the
    // estimator folds this back into calib1's basepoint on an extrinsic commit, converging the
    // recovered forward axis to truth from a large prior error. See THE MAPPING above; ROLL
    // and the TRANSLATION kept at the prior — Slice 7 fills roll + xyz). Returns the prior for
    // an unvoted / unknown source.
    SE3    extrinsic(SourceId id) const;
    // Per-source scale estimate (magnitude ratio vs the reference, mode). Returns 1.0 for
    // the reference / a scale_calib==false source / an unvoted source.
    Scalar scale(SourceId id) const;

    // Direction (so(3)) histogram CONCENTRATION in [0,1] — the yaw/pitch confidence.
    // Combined across the 3 channels as the MINIMUM concentration (the weakest channel
    // bounds the joint reliability). 0 for an unvoted / unknown source.
    Scalar extrinsic_confidence(SourceId id) const;
    // Scale histogram concentration in [0,1]. 0 for the reference / unvoted / unknown.
    Scalar scale_confidence(SourceId id) const;

    bool     configured() const { return configured_; }
    SourceId reference()  const { return reference_id_; }
    // Total so(3)-channel vote mass for source `id` (the N_min driver; 0 if none).
    Scalar   vote_count(SourceId id) const;
    bool     gated_straight() const { return last_gate_straight_; }   // last observe()

private:
    // Resolve / get-or-create a per-source slot (linear scan, bounded by source_count_).
    int slot_for(SourceId id) const;
    int ensure_slot(SourceId id);

    // The minimal rotation taking unit vector `a` to unit vector `b` (Rodrigues from the
    // cross/dot). Handles the antiparallel degenerate case with an arbitrary perpendicular
    // axis. Both inputs assumed unit-norm.
    static Mat3 rotation_between(const Vec3& a, const Vec3& b);

    // Per-vote weight factor honoring vote_weight_ (D5). `omega_norm` is ‖fused_omega‖ for
    // the step; `confidence` is the source's Σ-confidence (1 if not supplied). One -> 1;
    // Confidence -> confidence; Rotation -> max(omega_norm, kRotFloor) (floored so the
    // straight-regime ‖ω‖≈0 still produces a positive vote); Combo -> Rotation × Confidence.
    // Always strictly positive (add() ignores non-positive weights).
    Scalar vote_weight_factor(Scalar omega_norm, Scalar confidence) const;

    // Active configuration (configure() is the sole initializer; readers short-circuit on
    // !configured_).
    bool     configured_       = false;
    SourceId reference_id_      = 0;
    Scalar   straight_omega_max_ = Scalar(0.05);
    Scalar   straight_trans_min_ = Scalar(0.05);
    bool     reverse_fold_      = true;
    bool     ref_cross_check_   = false;
    VoteWeight vote_weight_     = VoteWeight::Combo;
    bool     last_gate_straight_ = false;
    HistogramConfig so3_cfg_{};
    HistogramConfig scale_cfg_{};

    // Per-source state. so3_[3*slot + c] is channel c (0=x,1=y,2=z) of the direction
    // histogram; scale_[slot] is the scale histogram. prior_[slot] is the full SE3 prior
    // extrinsic (its rotation R_Xprior is the basepoint the correction is applied on top of).
    Histogram1D so3_[3 * kMaxSources];
    Histogram1D scale_[kMaxSources];
    SE3         prior_[kMaxSources];
    bool        scale_calib_[kMaxSources] = {};
    SourceId    ids_[kMaxSources] = {};
    int         source_count_ = 0;
};

// ===========================================================================
// PHASE 2 — TURN-REGIME roll (about forward) + xyz lever-arm via HAND-EYE
// (DESIGN §3/§6, DECISIONS D5/D10/D21). The two extrinsic DOF Phase 1 cannot see:
// roll about the forward axis and the xyz lever arm are observable ONLY under
// rotation. A sensor offset by a lever arm `r` sees `ω×r`, which is zero in pure
// translation — so the translation hand-eye is SINGULAR at R_A=I, which re-derives
// the turn gate. Phase 1 has already fixed yaw,pitch; Phase 2 adds the remaining roll
// + translation on top.
//
// THE HAND-EYE EQUATION (the oracle pins it, sim/synthetic_source.cpp). A rigid sensor
// and base obey  A·X = X·B,  A = base motion, B = sensor motion, X = sensor->base. The
// sim reports  B = X^{-1} A X  (then B.t *= scale), so  A X = X B  exactly. Splitting:
//   rotation:     R_A R_X = R_X R_B
//   translation:  (R_A − I) t_X = R_X t_B − t_A
//
// ROLL (1 nonlinear DOF). With yaw,pitch fixed by Phase 1 as R_yp (the roll-free
// rotation R_X with roll=0 — exactly Phase1Calibrator::extrinsic's rotation), the full
// extrinsic rotation is  R_X(roll) = R_yp · Rx(roll),  Rx = rotation about the LOCAL
// forward axis (sensor x; Phase 1 has aligned it to base-forward). Recover roll per
// accepted window by a BOUNDED 1-D search minimizing the rotation hand-eye residual
//   res(roll) = ‖ log( R_X(roll) R_B R_X(roll)^T R_A^T ) ‖
// over roll ∈ (−π, π] (coarse grid, then a parabolic refine of the grid minimizer; the
// search is bounded — fixed grid count). At the true roll, R_A R_X = R_X R_B ⇒ res = 0;
// under turning the residual is sensitive to roll (it is flat only when R_A = I, i.e.
// no rotation — the same singularity that gates Phase 2). The minimizer is voted into a
// CIRCULAR S¹ Histogram1D (cfg.roll_hist); estimate = its mode.
//
// ROT3D (Slice 17, opt-in via Config::rot3d_enabled) — the FULL rotation extrinsic from
// the ROTATION half of the hand-eye, R_A R_X = R_X R_B, equivalently the rotation-AXIS
// correspondence  a = R_X · b  with  a = so3::log(R_A), b = so3::log(R_B)  per turn-gated
// window. Accumulate the per-slot Wahba problem  Mw = Σ w·a bᵀ  and the axis Gram
// BBw = Σ w·b bᵀ (both decay-aged by so3_hist.decay_gamma per ACCEPTED window — a fixed
// per-window decay independent of the histogram aging mode, since a SlidingK ring cannot
// evict from a rank accumulator). The LS rotation is the running Kabsch solve
//   R̂ = U·diag(1, 1, det(U Vᵀ))·Vᵀ,  U S Vᵀ = SVD(Mw).
// OBSERVABILITY lives in BBw: rank ≥ 2 (two non-parallel rotation axes) fully determines
// R_X; rank 1 (ground, yaw-only) leaves the rotation ABOUT the common axis blind — exactly
// the DOF Phase-1's straight-regime direction method covers (complementary by construction).
// The votable gate is λ_mid ≥ kAxisPairMin·λ_max on BBw (NOT λ_min — Kabsch needs only two
// axes). When votable, the basepoint-relative residual  δφ = so3::log(R̂ · R_bpᵀ)  is voted
// into 3 so(3) channels (cfg.so3_hist shape, basepoint re-anchorable — the contractive
// parametrization mirroring Phase-1); readout  rot3d() = exp(φ_mode)·R_bp. On clean
// conjugated data (B = X⁻¹∘A∘X) the recovered rotation satisfies a = rot3d()·b exactly.
// Rank-1 motion NEVER votes (channels stay empty, confidence 0, never commits) — planar
// behavior unchanged; no partial-subspace voting (covered by Phase 1 + roll).
//
// JOINT LEVER + SCALE (Slice 17b, opt-in via Config::joint_lever_scale) — the 4-unknown
// hand-eye translation LS. With the OBSERVED (prior-de-scaled) t_B = s_res·t_true the
// translation identity becomes LINEAR in u = [t_X; κ], κ = 1/s_res:
//     (R_A − I)·t_X − κ·(R_X t_B) = −t_A          (row block J = [(R_A − I) | −(R_X t_B)])
// accumulated into fixed 4×4/4×1 normal equations exactly where the 3-unknown rows
// accumulate today (same turn gate, same kRotRowMin guard, same R_X source — the running
// Kabsch R̂ when the rot3d gate is open, else R_yp ∘ Rx(roll)). The production 3-unknown
// solve is the κ≡1 special case; an unmodeled residual scale lands as a lever error along
// the dominant motion axis (Slice-17 finding: scale 1.08 → lever x off 3.5 cm). The
// running vote solve is ridge-regularized CENTERED ON [t_prior; 1] (κ prior = 1: the
// residual convention — matching Phase-1's scale()); per-axis info gating extends to the
// κ diagonal, with the LEVER axes judged against the lever sub-block's largest diagonal
// and κ against the full matrix's (mixed units: the κ column carries metres, the lever
// columns are dimensionless — a dominant κ diagonal must not mask an observable lever,
// and a starved κ — near-zero translation — must pin at 1 with scale2 confidence 0).
// Observable lever axes vote into the EXISTING xyz histograms (unchanged); the κ axis
// votes s_res = 1/κ̂ into a NEW per-source scale histogram (cfg.scale_hist shape — a
// residual ratio clustering at 1), read out via scale2()/scale2_confidence()/
// scale2_vote_count() (1/0/0 for the reference / a scale_calib==false source / unvoted —
// mirroring Phase-1's scale readouts). Unlike the 3-unknown path, the joint rows carry
// the per-source vote weight w (the κ estimate benefits from the same excitation/
// confidence weighting as its histogram votes). Default OFF keeps the 3-unknown path
// BYTE-IDENTICAL (separate accumulators; the off path is untouched code).
//
// xyz (3 linear DOF given rotation). Stack the translation row-block over accepted
// (turning) windows and solve LEAST-SQUARES for t_X:
//      Stack  (R_A − I) t_X = (R_X t_B − t_A)  =>  normal equations  AᵀA t_X = Aᵀb,
// The row's R_X is the R_yp ∘ Rx(roll) composition at the per-window roll — EXCEPT when
// rot3d is enabled AND its two-axis gate is open for the slot, in which case the running
// Kabsch solve R̂ drives the row instead (the Slice-17 lever-coupling win: on turn-only 3D
// motion Phase-1 starves, freezing R_yp at the prior, which would bias every row by ~|t|·θ;
// R̂ is the best available full rotation there). Below the gate / disabled, byte-identical.
// accumulated INCREMENTALLY into a fixed 3×3 (AᵀA) + 3×1 (Aᵀb) — no growing matrix
// (strict core). Solve via Eigen LDLT (symmetric-PSD). NEAR-SINGULAR GUARD: pure yaw
// turning leaves t_X.z unobservable (the null space of R_A−I for a z-axis rotation IS
// the z-axis), so AᵀA is rank-deficient; readout is gated on a conditioning floor
// (smallest eigenvalue of AᵀA ≥ cond_min · largest). Per-window a low-rotation row
// (‖log R_A‖ below a floor) is REJECTED before accumulation (the ω×r → 0 regime). Each
// solved component is voted into the 3×1-D cfg.xyz_hist; estimate = per-channel mode.
//
// STRATEGY (cfg.phase2_strategy, both implemented + compared — D10). The per-source
// hand-eye is identical; only the base motion A differs:
//   * VsFusedBase       : A = the fused CONSENSUS base motion (med over sources). R_yp
//                         for source i comes from Phase 1's extrinsic(i).
//   * PairwisePinnedRef : A = X_ref ∘ B_ref ∘ X_ref^{-1}, the base motion IMPLIED by the
//                         pinned reference (its extrinsic held at the prior — the gauge).
//                         Solving the same hand-eye then yields the source's X relative
//                         to the reference. Reference-paired (sufficient; not all-pairs).
// Both recover the planted roll + lever arm in the oracle.
//
// CONFIDENCE. extrinsic_confidence reuses the roll histogram concentration (the roll
// reliability); translation_confidence (a non-breaking new readout) is the MIN over the
// 3 xyz channels' concentrations (the weakest axis bounds the joint t_X reliability — an
// unobservable axis, e.g. z under pure yaw, never concentrates, so its confidence stays
// low and flags the under-determined component).
//
// Slice 7 ESTIMATES + EXPOSES only — NO feedback into fusion (that is Slice 8). It is
// driven from the estimator's turn-gated calibration slice; fusion behavior is unchanged.
//
// STRICT CORE: no heap after configure(); fixed per-source roll (1) + xyz (3)
// histograms + a fixed-size 3×3 normal-equation accumulator per source; bounded loops
// (the roll grid, the per-source scan); Status returns; no exceptions; double math.
class Phase2Calibrator {
public:
    static constexpr int kMaxSources = 32;
    // Bounded 1-D roll search resolution (strict core: fixed iteration count). FIXED
    // (not config) observability threshold the slice is judged on. The 180-sample grid is
    // ~2° over (−π, π] — a PRE-REFINE grid that ASSUMES a single-well rot_residual (smooth,
    // one minimum); a residual minimum narrower than ~2° between samples (e.g. extreme
    // yaw+roll coupling) could let the grid pick the wrong neighbourhood before the parabola.
    static constexpr int kRollGrid   = 180;     // grid samples over (−π, π]

    Phase2Calibrator() = default;

    // Bind the active configuration + the pinned reference source id. Sizes the
    // per-source roll (cfg.roll_hist, circular) + 3 xyz (cfg.xyz_hist) histograms,
    // captures the turn gate / strategy, and clears all state. Validates:
    //   reference_id < kMaxSources                -> else OutOfRange
    //   roll_hist / xyz_hist pass Histogram1D::configure -> propagated
    // No heap occurs here nor in any later call.
    Status configure(const Config& cfg, SourceId reference_id);

    // Drop all votes + reset the LS accumulators (keeps the configuration).
    void reset();

    // Register the per-source PRIOR extrinsic (sensor->base) — used by
    // PairwisePinnedRef as the pinned-reference gauge, and as the translation fallback
    // for an unvoted source. `scale_calib` (Slice 17b) mirrors Phase-1's: a source with
    // scale_calib == false never votes/reads the joint-scale (scale2) channel (its scale
    // is pinned at the prior); it defaults true, so every pre-17b caller is unchanged.
    // Call once per source after configure(). Returns OutOfRange if id is out of range,
    // CapacityExceeded at cap.
    Status set_prior(SourceId id, const SE3& prior_extrinsic, bool scale_calib = true);

    // --- Re-anchor API (Slice 8 feedback loop) -------------------------------
    // RE-ANCHOR the lever-arm prior + the PairwisePinnedRef gauge to a NEW committed
    // extrinsic (sensor->base). The xyz LS is centred on prior_.t (the ridge pins an
    // unobservable axis there) and the histogram falls back to prior_.t per channel, so a
    // freshly-committed extrinsic must replace the prior. Roll is voted ABSOLUTE about the
    // forward axis (basepoint R_yp fed via set_yaw_pitch each step), so the roll histogram
    // is re-anchored implicitly when R_yp moves — but the caller should reset_lever(id) so
    // the LS re-accumulates around the new basepoint (the row b = R_X t_B − t_A and the
    // ridge centre both change). Heap-free; OutOfRange / CapacityExceeded as set_prior.
    Status set_basepoint(SourceId id, const SE3& basepoint);
    // Drop only the roll histogram votes for `id` (keeps the lever-arm LS + prior). Used
    // when the roll basepoint R_yp moves (a yaw/pitch re-anchor) so the residual roll
    // re-concentrates. The post-reset roll() falls back to 0 until new votes land.
    void reset_roll(SourceId id);
    // Reset the lever-arm LS (AᵀA, Aᵀb, row count) AND the 3 xyz histograms for `id`
    // (keeps roll + prior). Used on a lever-arm / extrinsic commit+re-anchor: the LS rows
    // were accumulated against the OLD base/extrinsic frame, so they are stale — drop them
    // and re-accumulate around the new basepoint. lever_arm() falls back to prior_.t.
    // With joint_lever_scale (Slice 17b) this clears the FULL joint state: the 4×4/4×1
    // accumulators AND the scale2 histogram (the κ votes were running solutions of the
    // same stale rows).
    void reset_lever(SourceId id);
    // Drop the joint-scale (scale2) histogram votes AND the joint 4×4/4×1 accumulator +
    // row count for `id`, KEEPING the xyz histograms (Slice 17b). Used on a SCALE
    // commit+re-anchor (either scale path's fold into prior_scale): the accumulated rows'
    // t_B were de-scaled by the OLD prior_scale, so the single-κ joint model cannot mix
    // pre-/post-fold rows (the running κ̂ would blend the two de-scale epochs and re-vote
    // a residual ≠ 1 forever); the xyz votes are KEPT because the joint solve they were
    // read from is scale-immune — they remain valid lever estimates across the fold.
    // The post-reset scale2() falls back to 1 until new votes land. No-op when the joint
    // path is disabled (nothing accumulated).
    void reset_scale2(SourceId id);
    // RE-ANCHOR the rot3d so(3) basepoint to a NEW committed rotation (Slice 17). The rot3d
    // votes are δφ = log(R̂·R_bpᵀ), so moving the basepoint onto a freshly-committed rot3d
    // makes the residual re-vote ≈ 0 around it (the contractive iterate, mirroring Phase-1's
    // set_basepoint). NOTE: deliberately SEPARATE from set_basepoint/set_prior's every-step
    // value sync — the rot3d basepoint must only move together with a reset_rot3d (moving it
    // under live votes would shift what the accumulated δφ mode means). The Mw/BBw
    // accumulators are basepoint-INDEPENDENT and are intentionally NOT touched here.
    Status set_rot3d_basepoint(SourceId id, const Mat3& R_bp);
    // Drop only the 3 rot3d so(3) channels for `id` (keeps Mw/BBw — they are basepoint-
    // independent rank/solve state and keep accumulating across re-anchors). The post-reset
    // rot3d() falls back to the basepoint until new votes land.
    void reset_rot3d(SourceId id);

    // Provide the Phase-1-recovered yaw/pitch rotation R_yp for each source (the roll
    // basepoint: R_X(roll) = R_yp · Rx(roll)). Call whenever Phase 1 updates (cheap; the
    // calibrator just stores the latest). A source with no R_yp set uses its prior
    // rotation. Returns OutOfRange / CapacityExceeded as set_prior.
    Status set_yaw_pitch(SourceId id, const Mat3& R_yp);

    // Consume one FUSED step (the turn-regime calibration slice). `fused_motion` is the
    // consensus base motion A over the step (SE3; used by VsFusedBase and for the gate).
    // `fused_omega` is the consensus body angular RATE (rad/s, ÷dt) — the turn gate's
    // operand. Per registered source i in [0, n): ids[i] is its id, reported[i] is its
    // DE-SCALED reported sensor-frame delta B_corr (same convention Phase 1 uses). The
    // reference source must be present among `ids` for PairwisePinnedRef to vote.
    //
    // GATE (D5): votes are deposited ONLY when the fused motion is TURNING —
    //   ‖fused_omega‖ > turn_omega_min.
    // Off-gate steps are a no-op (NotReady, not an error). Per-source, a window whose
    // base-motion rotation magnitude ‖log R_A‖ is below a small floor is additionally
    // skipped (the ω×r→0 near-singular row — it would add a near-zero, ill-conditioned
    // row to the LS). Returns Ok when at least one source voted, NotReady when gated out,
    // NotInitialized if unconfigured, NoData on bad args.
    //
    // VOTE WEIGHT (cfg.vote_weight, D5 "votes weighted by ‖ω‖"). Each roll + xyz vote's
    // mass is scaled by vote_weight_factor(): One -> 1; Rotation -> ‖fused_omega‖ (the
    // lever-arm benefits from more rotation); Confidence -> the source's Σ-confidence
    // (`confidences[i]`, the same scalar fusion weights with); Combo -> Rotation ×
    // Confidence. `confidences` is OPTIONAL (nullptr -> unit confidence). The xyz LS rows
    // are themselves unweighted (the normal equations stay in physical units); the weight
    // applies to the HISTOGRAM vote of the running solution, matching Phase-1.
    Status observe(int n, const SourceId* ids, const SE3* reported,
                   const SE3& fused_motion, const Vec3& fused_omega,
                   const Scalar* confidences = nullptr);

    // --- Per-source readouts (committed = histogram mode / LS solve) ---------
    // Recovered roll about the forward axis [rad], the roll-histogram mode. 0 (prior)
    // for an unvoted / unknown source.
    Scalar roll(SourceId id) const;
    // Recovered lever-arm translation t_X (sensor->base), per-channel xyz-histogram mode.
    // The prior translation for an unvoted source / an unobservable (ill-conditioned)
    // axis. (The histogram mode is the robust estimate; the raw LS solve is available via
    // solve_lever_arm() for diagnostics.)
    Vec3 lever_arm(SourceId id) const;
    // Full recovered extrinsic: rotation = R_yp · Rx(roll_mode) (yaw/pitch from Phase 1 ∘
    // recovered roll), translation = lever_arm(). The prior for an unvoted source.
    SE3  extrinsic(SourceId id) const;

    // --- rot3d readouts (Slice 17; meaningful only when cfg.rot3d_enabled) ----
    // Full recovered extrinsic ROTATION from the axis-correspondence Wahba solve:
    //   rot3d() = exp(φ_mode) · R_bp   (φ_mode = the 3 rot3d channel modes, R_bp = the
    // rot3d basepoint). Falls back to the basepoint (= the prior rotation until re-anchored)
    // for an unvoted / unknown source. CONVENTION (pinned by tests): on clean conjugated
    // data B = X⁻¹∘A∘X the returned R satisfies  log(R_A) = R · log(R_B)  exactly, i.e. it
    // IS the sensor->base R_X (a planted +yaw reads +yaw, matching config_loader's
    // "yaw pitch roll" parse and tools/inject_calib.py).
    Mat3   rot3d(SourceId id) const;
    // MIN concentration over the 3 rot3d channels in [0,1]; 0 unvoted / unknown.
    Scalar rot3d_confidence(SourceId id) const;
    // Total rot3d channel vote mass (the N_min driver; 0 if none).
    Scalar rot3d_vote_count(SourceId id) const;
    // Whether the slot's BBw currently clears the two-axis gate (λ_mid ≥ kAxisPairMin·λ_max).
    // Diagnostics / observability self-tests; rank-1 (planar) motion never clears it.
    bool   rot3d_observable(SourceId id) const;

    // --- joint-scale readouts (Slice 17b; meaningful only when cfg.joint_lever_scale) --
    // Turn-regime per-source residual scale s_res = 1/κ̂ (the scale2 histogram mode) —
    // the SAME residual-vs-prior convention as Phase-1's scale() (t_B arrives prior-de-
    // scaled, so a converged calibration re-votes ≈ 1). Returns 1.0 for the reference /
    // a scale_calib==false source / an unvoted source / when the joint path is disabled.
    Scalar scale2(SourceId id) const;
    // scale2 histogram concentration in [0,1]. 0 for the reference / scale_calib==false /
    // unvoted / unknown / joint path disabled — in particular a κ starved of information
    // (negligible translation) never votes, so its confidence stays 0 (the per-DOF spine).
    Scalar scale2_confidence(SourceId id) const;
    // Total scale2 vote mass (the N_min driver; 0 if none).
    Scalar scale2_vote_count(SourceId id) const;

    // Roll-histogram concentration in [0,1] — the roll (extrinsic) confidence. 0 unvoted.
    Scalar extrinsic_confidence(SourceId id) const;
    // Translation confidence: MIN of the 3 xyz-channel concentrations in [0,1]. The
    // weakest (e.g. an unobservable z under pure yaw) bounds the joint t_X reliability.
    // 0 for an unvoted / unknown source.
    Scalar translation_confidence(SourceId id) const;

    bool     configured() const { return configured_; }
    SourceId reference()  const { return reference_id_; }
    Scalar   roll_vote_count(SourceId id) const;     // roll-histogram total mass
    Scalar   xyz_vote_count(SourceId id) const;      // accumulated LS row count
    bool     gated_turning() const { return last_gate_turning_; }   // last observe()

    // Direct LS solve of t_X from the accumulated normal equations (AᵀA t_X = Aᵀb), with
    // the conditioning guard. Returns the prior translation if ill-conditioned / no rows.
    // `observable` (out) flags whether the solve was well-conditioned. Diagnostics path
    // (the histogram mode is the committed estimate); also used by lever_arm()'s seeding.
    // JOINT MODE (Slice 17b): the solve reads the 4×4 joint system, but LEVER observability
    // is judged on the LEVER SUB-PROBLEM — the kCondMin eigenvalue-ratio guard runs on the
    // κ-marginalized 3×3 Schur complement (with κ ridge-pinned toward its prior 1 when
    // under-informed), so a dominant/starved κ diagonal can neither mask nor fake lever
    // conditioning (planar yaw-only still reports the prior: the z null direction survives
    // the marginalization). When the guard passes, the κ-ridged 4×4 is solved centered on
    // [t_prior; 1] and the lever components returned.
    Vec3 solve_lever_arm(SourceId id, bool* observable = nullptr) const;

private:
    // Fixed-size joint-system types (Slice 17b; strict core — Eigen fixed 4×4/4×1).
    using Mat4 = Eigen::Matrix<Scalar, 4, 4>;
    using Vec4 = Eigen::Matrix<Scalar, 4, 1>;
    int slot_for(SourceId id) const;
    int ensure_slot(SourceId id);

    // Rotation about the local forward axis (sensor x) by `roll` (Rodrigues, closed form).
    static Mat3 roll_about_forward(Scalar roll);
    // Ridge-regularized solve of the lever-arm normal equations CENTERED on the prior:
    //   (AᵀA + ridge·I)(t − t_prior) = Aᵀb − AᵀA·t_prior  =>  t.
    // An unobservable axis (≈zero information in AᵀA) is pinned near t_prior by the ridge;
    // observable axes resolve. `obs[c]` (out) flags whether axis c carries enough
    // information (AᵀA(c,c) ≥ kAxisInfoMin × max diagonal). Heap-free.
    static Vec3 solve_ridge(const Mat3& AtA, const Vec3& Atb, const Vec3& t_prior,
                            Scalar ridge, bool obs[3]);
    // Joint-system analogue (Slice 17b): ridge-regularized solve of the 4×4 normal
    // equations CENTERED on u_prior = [t_prior; 1]. Per-axis info flags `obs[4]`: the 3
    // LEVER axes are judged against the LEVER sub-block's largest diagonal (so a dominant
    // κ diagonal — metres² vs the dimensionless lever columns — cannot mask an observable
    // lever axis), κ against the FULL matrix's largest diagonal (so a starved κ — near-
    // zero translation — pins at its prior 1, never voted). The ridge is likewise
    // per-block-scaled. Heap-free (fixed 4×4 LDLT).
    static Vec4 solve_ridge4(const Mat4& AtA, const Vec4& Atb, const Vec4& u_prior,
                             Scalar ridge, bool obs[4]);
    // Per-vote weight factor honoring vote_weight_ (D5). `omega_norm` is ‖fused_omega‖;
    // `confidence` is the source's Σ-confidence (1 if not supplied). One -> 1; Rotation ->
    // omega_norm (already > turn_omega_min by the gate, so always positive); Confidence ->
    // confidence; Combo -> Rotation × Confidence. Always strictly positive.
    Scalar vote_weight_factor(Scalar omega_norm, Scalar confidence) const;

    // Two-axis observability gate on the axis Gram (Slice 17): true iff the MIDDLE
    // eigenvalue clears kAxisPairMin × the largest (two non-parallel rotation axes have
    // been excited — Kabsch needs only two correspondences, so the gate reads λ_mid, NOT
    // λ_min). Rank-1 (planar yaw-only) stays below the floor forever. Heap-free (Eigen 3×3
    // SelfAdjointEigenSolver, precedented by the lever conditioning guard).
    static bool axis_pair_observable(const Mat3& BBw);
    // Running Wahba/Kabsch solve of Mw = Σ w·a bᵀ (Slice 17): R̂ = U·diag(1,1,det(U Vᵀ))·Vᵀ
    // from SVD(Mw). Returns false (R_out untouched) on a non-finite result. Heap-free
    // (Eigen 3×3 JacobiSVD, bounded).
    static bool wahba_solve(const Mat3& Mw, Mat3& R_out);

    // The rotation hand-eye residual norm for a candidate roll (see header math).
    static Scalar rot_residual(const Mat3& R_yp, Scalar roll,
                               const Mat3& R_A, const Mat3& R_B);
    // Bounded 1-D minimizer of rot_residual over (−π, π] (grid + parabolic refine).
    static Scalar best_roll(const Mat3& R_yp, const Mat3& R_A, const Mat3& R_B);

    bool     configured_        = false;
    SourceId reference_id_       = 0;
    Scalar   turn_omega_min_     = Scalar(0.20);
    Phase2Strat strategy_        = Phase2Strat::VsFusedBase;
    VoteWeight vote_weight_      = VoteWeight::Combo;
    bool     last_gate_turning_  = false;
    bool     rot3d_enabled_      = false;            // Config::rot3d_enabled (Slice 17)
    bool     joint_              = false;            // Config::joint_lever_scale (Slice 17b)
    Scalar   rot3d_gamma_        = Scalar(0.999);    // Mw/BBw per-window decay (so3 gamma)
    HistogramConfig roll_cfg_{};
    HistogramConfig xyz_cfg_{};
    HistogramConfig rot3d_cfg_{};                    // = cfg.so3_hist (Slice 17)
    HistogramConfig scale2_cfg_{};                   // = cfg.scale_hist (Slice 17b)

    // Per-source state. roll_[slot] is the circular roll histogram; xyz_[3*slot + c] is
    // channel c (0=x,1=y,2=z) of the lever-arm histogram. prior_[slot] is the full SE3
    // prior extrinsic; ryp_[slot] is the Phase-1 yaw/pitch rotation (roll basepoint).
    // ata_[slot] / atb_[slot] are the incremental 3×3 / 3×1 normal-equation accumulators
    // for the lever-arm LS; rows_[slot] counts accumulated rows.
    Histogram1D roll_[kMaxSources];
    Histogram1D xyz_[3 * kMaxSources];
    SE3         prior_[kMaxSources];
    Mat3        ryp_[kMaxSources];
    Mat3        ata_[kMaxSources];
    Vec3        atb_[kMaxSources];
    Scalar      rows_[kMaxSources] = {};
    bool        ryp_set_[kMaxSources] = {};
    SourceId    ids_[kMaxSources] = {};
    int         source_count_ = 0;

    // rot3d state (Slice 17). rot3d_[3*slot + c] is channel c of the basepoint-relative
    // so(3)-residual histogram; mw_/bbw_ the decay-aged Wahba / axis-Gram accumulators;
    // rot3d_bp_ the re-anchorable rotation basepoint (init = the slot prior's rotation).
    Histogram1D rot3d_[3 * kMaxSources];
    Mat3        mw_[kMaxSources];
    Mat3        bbw_[kMaxSources];
    Mat3        rot3d_bp_[kMaxSources];

    // Joint lever+scale state (Slice 17b, only touched when joint_). ata4_/atb4_ are the
    // 4×4/4×1 joint normal-equation accumulators (unknowns [t_X; κ] — the legacy ata_/atb_
    // stay untouched for the byte-identical off path); scale2_[slot] is the per-source
    // residual-scale histogram (s_res = 1/κ̂ votes, cfg.scale_hist shape);
    // scale2_calib_[slot] mirrors Phase-1's scale_calib gate.
    Histogram1D scale2_[kMaxSources];
    Mat4        ata4_[kMaxSources];
    Vec4        atb4_[kMaxSources];
    bool        scale2_calib_[kMaxSources] = {};
};

} // namespace ofc
#endif // OFC_CORE_CALIBRATION_HPP
