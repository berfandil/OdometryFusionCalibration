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
// pitch = atan2(-f_z, hypot(f_x, f_y)). The reported extrinsic rotation is
// R = δR * R_Xprior (corrects the prior so the forward axis aligns); roll and the
// translation are kept at the prior (roll is Slice 7, xyz lever-arm is Slice 7).
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
    // Returns Ok when the step was voted, NotReady when gated out (not an error),
    // NotInitialized if unconfigured.
    Status observe(int n, const SourceId* ids, const SE3* reported,
                   const Vec3& fused_omega, const Vec3& fused_trans);

    // --- Per-source readouts (committed = histogram mode) --------------------
    // Recovered forward axis in the BASE frame (unit vector). Returns the prior forward
    // (R_Xprior * e_x) for an unvoted / unknown source.
    Vec3   forward_axis(SourceId id) const;
    Scalar yaw(SourceId id)   const;        // atan2(f_y, f_x)   [rad]
    Scalar pitch(SourceId id) const;        // atan2(-f_z, hypot(f_x,f_y))  [rad]
    // Recovered extrinsic: rotation = δR_mode * R_Xprior (yaw/pitch corrected; ROLL and
    // the TRANSLATION kept at the prior — Slice 7 fills roll + xyz). Returns the prior
    // for an unvoted / unknown source.
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

    // Active configuration (configure() is the sole initializer; readers short-circuit on
    // !configured_).
    bool     configured_       = false;
    SourceId reference_id_      = 0;
    Scalar   straight_omega_max_ = Scalar(0.05);
    Scalar   straight_trans_min_ = Scalar(0.05);
    bool     reverse_fold_      = true;
    bool     ref_cross_check_   = false;
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
// xyz (3 linear DOF given rotation). Stack the translation row-block over accepted
// (turning) windows and solve LEAST-SQUARES for t_X:
//      Stack  (R_A − I) t_X = (R_X t_B − t_A)  =>  normal equations  AᵀA t_X = Aᵀb,
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
    // Bounded 1-D roll search resolution (strict core: fixed iteration count).
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
    // for an unvoted source. Call once per source after configure(). Returns OutOfRange
    // if id is out of range, CapacityExceeded at cap.
    Status set_prior(SourceId id, const SE3& prior_extrinsic);

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
    Status observe(int n, const SourceId* ids, const SE3* reported,
                   const SE3& fused_motion, const Vec3& fused_omega);

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
    Vec3 solve_lever_arm(SourceId id, bool* observable = nullptr) const;

private:
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
    // The rotation hand-eye residual norm for a candidate roll (see header math).
    static Scalar rot_residual(const Mat3& R_yp, Scalar roll,
                               const Mat3& R_A, const Mat3& R_B);
    // Bounded 1-D minimizer of rot_residual over (−π, π] (grid + parabolic refine).
    static Scalar best_roll(const Mat3& R_yp, const Mat3& R_A, const Mat3& R_B);

    bool     configured_        = false;
    SourceId reference_id_       = 0;
    Scalar   turn_omega_min_     = Scalar(0.20);
    Phase2Strat strategy_        = Phase2Strat::VsFusedBase;
    bool     last_gate_turning_  = false;
    HistogramConfig roll_cfg_{};
    HistogramConfig xyz_cfg_{};

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
};

} // namespace ofc
#endif // OFC_CORE_CALIBRATION_HPP
