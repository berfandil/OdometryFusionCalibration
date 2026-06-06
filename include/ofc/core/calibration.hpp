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
// REVERSE-FOLD: if g_obs points backward (dot(g_obs, fwd_sign*e_x) < 0) negate it, so
// forward and reverse straight segments both land in the SAME forward hemisphere
// (unimodal). The candidate rotation is the MINIMAL rotation in the base frame taking
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

    // Consume one FUSED step. `fused_omega` / `fused_v` are the consensus body twist
    // angular / linear parts over the step (rad/s, m/s) used by the straight gate;
    // `fused_trans` is the consensus translation over the step (m) used to fix the base
    // forward sign. Per registered source i in [0, n): `ids[i]` is its id and
    // `reported[i]` is its DE-SCALED reported sensor-frame delta B_corr (B with the
    // translation divided by prior_scale, matching the estimator's pre-median de-scale —
    // so the magnitude ratio recovers the RESIDUAL scale error vs the prior). The
    // reference source must be present among `ids` for a scale vote to occur.
    //
    // GATE: votes are deposited ONLY when the fused motion is straight —
    //   ||fused_omega|| < straight_omega_max  AND  ||fused_trans|| > straight_trans_min
    // — and, when ref_cross_check is on, the reference source's reported delta also reads
    // "straight forward/back" (small rotation + sizable translation). Off-gate steps are a
    // no-op (the observability spine: yaw/pitch/scale are unobservable unless straight).
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
    // histogram; scale_[slot] is the scale histogram. prior_R_[slot] is R_Xprior.
    Histogram1D so3_[3 * kMaxSources];
    Histogram1D scale_[kMaxSources];
    SE3         prior_[kMaxSources];
    bool        scale_calib_[kMaxSources] = {};
    SourceId    ids_[kMaxSources] = {};
    int         source_count_ = 0;
};

} // namespace ofc
#endif // OFC_CORE_CALIBRATION_HPP
