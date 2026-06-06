// ofc/core/calibration.cpp — Phase-1 (straight-regime) calibration. See calibration.hpp
// for the full direction->yaw/pitch mapping, the scale derivation, the gate, and the
// reverse-fold; this file is the bounded, no-heap implementation.
//
// STRICT CORE: all storage is fixed-capacity member arrays bound in configure();
// set_prior()/observe()/readouts allocate nothing and run loops bounded by the live
// source count or the channel count. No exceptions; status-code returns; double math.
#include "ofc/core/calibration.hpp"

#include "ofc/core/lie.hpp"

#include <Eigen/Cholesky>     // LDLT (lever-arm normal-equation solve)
#include <Eigen/Eigenvalues>  // SelfAdjointEigenSolver (conditioning guard)

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
constexpr Scalar kUnitEps = Scalar(1e-12);   // guard ||t|| ~ 0 before normalizing
constexpr Scalar kPi      = Scalar(3.14159265358979323846);

// Floor for the Rotation vote-weight factor in the STRAIGHT regime. Phase-1 gates on
// ‖ω‖ < straight_omega_max (near zero), so a raw ‖ω‖ Rotation factor would collapse all
// votes to ~0 mass; floor it to a small positive so every gated-in vote still lands (the
// straight regime has no rotation excitation to differentiate — Rotation weighting is a
// no-op signal here by design; the knob is meaningful in Phase-2's turn regime).
constexpr Scalar kRotFloorP1 = Scalar(1e-3);

// The base forward axis (the gauge convention: in straight motion the base travels
// along ±e_x of its own frame — DESIGN §3).
const Vec3 kFwd = (Vec3() << Scalar(1), Scalar(0), Scalar(0)).finished();
} // namespace

Mat3 Phase1Calibrator::rotation_between(const Vec3& a, const Vec3& b) {
    // Minimal rotation taking unit vector a -> unit vector b (Rodrigues from cross/dot).
    // CALLER CONTRACT: observe() applies the π-singularity guard (c = cos(a,b) ≥ 0) before
    // calling this, so on the vote path θ < 90° and the result's so3::log is well away from
    // the θ=π singularity. The antiparallel (c<0) branch below is a defensive fallback for
    // any other caller; it builds a true 180° rotation whose so3::log is itself singular, so
    // it must NOT be fed to so3::log without first being guarded by the caller.
    const Vec3 v = a.cross(b);
    const Scalar s = v.norm();                 // sin(theta)
    const Scalar c = a.dot(b);                 // cos(theta)
    if (s < kUnitEps) {
        // Parallel (c>0 -> identity) or antiparallel (c<0 -> 180° about any axis ⟂ a).
        if (c > Scalar(0)) return Mat3::Identity();
        // Pick an axis perpendicular to a (the smaller-component cross is well-conditioned).
        Vec3 perp = a.cross(Vec3(Scalar(1), Scalar(0), Scalar(0)));
        if (perp.norm() < kUnitEps) perp = a.cross(Vec3(Scalar(0), Scalar(1), Scalar(0)));
        perp.normalize();
        return so3::exp(perp * kPi);
    }
    const Mat3 vx = so3::hat(v);
    // R = I + [v]x + [v]x^2 * (1 - c) / s^2
    return Mat3::Identity() + vx + vx * vx * ((Scalar(1) - c) / (s * s));
}

Scalar Phase1Calibrator::vote_weight_factor(Scalar omega_norm, Scalar confidence) const {
    // D5 vote weighting. One -> 1; Confidence -> the source Σ-confidence; Rotation ->
    // ‖ω‖ floored (straight regime ω≈0); Combo -> Rotation × Confidence. Always > 0 so
    // add() never drops the vote.
    const Scalar rot  = std::max(kRotFloorP1, omega_norm);
    const Scalar conf = std::max(Scalar(0), confidence);
    switch (vote_weight_) {
        case VoteWeight::One:        return Scalar(1);
        case VoteWeight::Confidence: return std::max(kUnitEps, conf);
        case VoteWeight::Rotation:   return rot;
        case VoteWeight::Combo:      return std::max(kUnitEps, rot * conf);
    }
    return Scalar(1);
}

Status Phase1Calibrator::configure(const Config& cfg, SourceId reference_id) {
    if (static_cast<int>(reference_id) >= kMaxSources) return Status::OutOfRange;

    so3_cfg_   = cfg.so3_hist;
    scale_cfg_ = cfg.scale_hist;

    // Configure all per-source histograms up front (one set per slot). Propagate the
    // first failure so a bad so3_hist / scale_hist is surfaced at configure().
    for (int i = 0; i < 3 * kMaxSources; ++i) {
        const Status hs = so3_[i].configure(so3_cfg_);
        if (!ok(hs)) return hs;
    }
    for (int i = 0; i < kMaxSources; ++i) {
        const Status hs = scale_[i].configure(scale_cfg_);
        if (!ok(hs)) return hs;
    }

    reference_id_       = reference_id;
    straight_omega_max_ = std::max(Scalar(0), cfg.straight_omega_max);
    straight_trans_min_ = std::max(Scalar(0), cfg.straight_trans_min);
    reverse_fold_       = cfg.reverse_fold;
    ref_cross_check_    = cfg.ref_cross_check;
    vote_weight_        = cfg.vote_weight;
    last_gate_straight_ = false;

    configured_ = true;
    reset();
    return Status::Ok;
}

void Phase1Calibrator::reset() {
    for (int i = 0; i < 3 * kMaxSources; ++i) so3_[i].reset();
    for (int i = 0; i < kMaxSources; ++i) {
        scale_[i].reset();
        prior_[i]       = SE3{};
        scale_calib_[i] = true;
        ids_[i]         = 0;
    }
    source_count_       = 0;
    last_gate_straight_ = false;
}

int Phase1Calibrator::slot_for(SourceId id) const {
    for (int i = 0; i < source_count_; ++i) {
        if (ids_[i] == id) return i;
    }
    return -1;
}

int Phase1Calibrator::ensure_slot(SourceId id) {
    const int s = slot_for(id);
    if (s >= 0) return s;
    if (source_count_ >= kMaxSources) return -1;
    const int slot = source_count_++;
    ids_[slot]         = id;
    prior_[slot]       = SE3{};
    scale_calib_[slot] = true;
    return slot;
}

Status Phase1Calibrator::set_prior(SourceId id, const SE3& prior_extrinsic,
                                   bool scale_calib) {
    if (!configured_)                          return Status::NotInitialized;
    if (static_cast<int>(id) >= kMaxSources)   return Status::OutOfRange;
    const int slot = ensure_slot(id);
    if (slot < 0)                              return Status::CapacityExceeded;
    prior_[slot]       = prior_extrinsic;
    scale_calib_[slot] = scale_calib;
    return Status::Ok;
}

Status Phase1Calibrator::set_basepoint(SourceId id, const SE3& basepoint) {
    // Re-anchor the so(3) basepoint to a committed extrinsic (Slice 8). Same validation /
    // slot semantics as set_prior; only the prior SE3 is replaced (scale_calib kept).
    if (!configured_)                        return Status::NotInitialized;
    if (static_cast<int>(id) >= kMaxSources) return Status::OutOfRange;
    const int slot = ensure_slot(id);
    if (slot < 0)                            return Status::CapacityExceeded;
    prior_[slot] = basepoint;
    return Status::Ok;
}

void Phase1Calibrator::reset_so3(SourceId id) {
    if (!configured_) return;
    const int s = slot_for(id);
    if (s < 0) return;
    so3_[3 * s + 0].reset();
    so3_[3 * s + 1].reset();
    so3_[3 * s + 2].reset();
}

void Phase1Calibrator::reset_scale(SourceId id) {
    if (!configured_) return;
    const int s = slot_for(id);
    if (s < 0) return;
    scale_[s].reset();
}

Status Phase1Calibrator::observe(int n, const SourceId* ids, const SE3* reported,
                                 const Vec3& fused_omega, const Vec3& fused_trans,
                                 const Scalar* confidences) {
    if (!configured_)                       return Status::NotInitialized;
    if (n <= 0 || ids == nullptr || reported == nullptr) return Status::NoData;

    // Rotation factor for vote weighting is the consensus turn magnitude (D5). In the
    // straight regime this is near-zero (gated < straight_omega_max), so the factor is
    // floored inside vote_weight_factor(); see that helper + kRotFloorP1.
    const Scalar omega_norm = fused_omega.norm();

    // --- Straight gate (the observability spine, D5) -----------------------------
    // Accept only near-zero rotation AND sizable translation in the FUSED consensus
    // motion. ||omega|| is the angular speed (rad/s); ||fused_trans|| the translation
    // magnitude over the step (m). Gated-out steps are a no-op (not an error).
    const bool straight = (fused_omega.norm() < straight_omega_max_) &&
                          (fused_trans.norm() > straight_trans_min_);
    last_gate_straight_ = straight;
    if (!straight) return Status::NotReady;

    // CONSENSUS forward/reverse sign (DESIGN §6: "fold into the consensus hemisphere").
    // In the straight regime the base travels along ±e_x; the sign of the FUSED translation
    // projected on the base-forward axis says whether this step is forward (+) or reverse
    // (−). We fold each source's observed direction by THIS sign (not by the source's own
    // g_obs.x, which is pure noise for a sideways mount), so forward and reverse segments
    // land in the same hemisphere for ANY mount orientation. The straight gate's
    // ‖fused_trans‖ > straight_trans_min_ already excludes the near-zero (ambiguous-sign)
    // case; an exact zero x-projection ties deterministically to forward (no negation).
    const Scalar fwd_proj = fused_trans.dot(kFwd);
    const Scalar consensus_sign = (fwd_proj < Scalar(0)) ? Scalar(-1) : Scalar(1);

    // Optional reference cross-check (ice-slide / drone niche): require the reference
    // source's OWN reported delta to also read straight forward/back — a small rotation
    // and a translation that, mapped to the base frame through its prior, lies near the
    // base forward axis (±e_x). If the reference is absent from this step, the cross-check
    // cannot be satisfied. The fold is against the canonical +e_x (the prior basepoint).
    if (ref_cross_check_) {
        bool ref_ok = false;
        for (int i = 0; i < n; ++i) {
            if (ids[i] != reference_id_) continue;
            const int rs = slot_for(reference_id_);
            const SE3 Xr = (rs >= 0) ? prior_[rs] : SE3{};
            const Vec3 bt = reported[i].t;
            const Scalar bn = bt.norm();
            const Scalar rot = so3::log(reported[i].R).norm();
            if (bn > straight_trans_min_ && rot < straight_omega_max_) {
                // Reference forward in the base frame through its prior, folded to +e_x.
                Vec3 g = Xr.R * (bt / bn);
                if (g.dot(kFwd) < Scalar(0)) g = -g;
                // "Mostly forward": within ~30° of +e_x (a generous straight-read tol for
                // the niche cross-check). An off-axis (e.g. sideways) reference fails this.
                if (g.dot(kFwd) > Scalar(0.866)) ref_ok = true;
            }
            break;
        }
        if (!ref_ok) return Status::NotReady;
    }

    // Reference reported magnitude (gauge denominator for the scale ratio). Found in this
    // step's sources; if absent, scale voting is skipped this step (extrinsic still votes).
    Scalar ref_mag = Scalar(0);
    bool   have_ref = false;
    for (int i = 0; i < n; ++i) {
        if (ids[i] == reference_id_) {
            ref_mag = reported[i].t.norm();
            have_ref = (ref_mag > kUnitEps);
            break;
        }
    }

    // --- Per-source votes --------------------------------------------------------
    for (int i = 0; i < n; ++i) {
        const SourceId id = ids[i];
        const int slot = ensure_slot(id);
        if (slot < 0) continue;                       // at capacity

        const Vec3 bt = reported[i].t;
        const Scalar bn = bt.norm();
        if (bn < kUnitEps) continue;                  // no translation -> no direction

        // Observed forward direction in the BASE frame through the prior extrinsic.
        // When the mount matches the prior this sits near +e_x for FORWARD base motion
        // and near -e_x for REVERSE base motion.
        const Vec3 dir_B = bt / bn;
        Vec3 g_obs = prior_[slot].R * dir_B;
        // REVERSE-FOLD (D5, DESIGN §6): fold a reverse-segment sample into the SAME
        // hemisphere as the forward segments using the CONSENSUS sign (sign of the fused
        // translation on the base-forward axis), not the source's own g_obs.x. For a
        // sideways/far-off mount g_obs.x is noise, so a fixed g_obs.x<0 test would split the
        // forward and reverse samples across antipodal peaks; the consensus sign folds them
        // correctly for any mount orientation. With the fold OFF (reverse_fold_ == false) we
        // do NOT negate — a reverse sample stays backward (≥90° off +e_x) and is then dropped
        // by the π-guard below (never voted; not edge-clamped boundary mass as before).
        if (reverse_fold_ && consensus_sign < Scalar(0)) g_obs = -g_obs;
        const Scalar gn = g_obs.norm();
        if (gn < kUnitEps) continue;
        const Vec3 g_unit = g_obs / gn;

        // π-SINGULARITY GUARD. Phase-1 is the SMALL-deviation-from-prior regime: a valid
        // mount's folded forward direction sits near +e_x. If the folded direction is ≥90°
        // off +e_x (cos(e_x, g_unit) < 0) the candidate rotation approaches 180°, whose
        // so3::log is singular (s = θ/(2·sinθ) → ÷0 at θ=π → NaN/huge, lie.cpp:49). Such a
        // sample is outside Phase-1's regime, so SKIP it rather than vote a π-rotation log —
        // this guarantees no NaN/clamped-π mass ever enters the histogram.
        if (g_unit.dot(kFwd) < Scalar(0)) continue;

        // Per-source vote weight (D5): scale every channel by the configured factor (the
        // turn magnitude and/or this source's Σ-confidence). `confidences` is optional —
        // a null array gives unit confidence (Confidence/Combo collapse to Rotation).
        const Scalar conf = (confidences != nullptr) ? confidences[i] : Scalar(1);
        const Scalar w    = vote_weight_factor(omega_norm, conf);

        // Candidate rotation δR taking e_x -> g_unit; vote its so(3) log (3 channels).
        // The guard above keeps θ < 90°, well clear of the so3::log singularity at π.
        const Mat3 dR  = rotation_between(kFwd, g_unit);
        const Vec3 phi = so3::log(dR);
        so3_[3 * slot + 0].add(phi.x(), w);
        so3_[3 * slot + 1].add(phi.y(), w);
        so3_[3 * slot + 2].add(phi.z(), w);

        // Scale vote = magnitude ratio vs the reference's reported magnitude (same weight).
        if (have_ref && scale_calib_[slot] && id != reference_id_) {
            scale_[slot].add(bn / ref_mag, w);
        }
    }

    return Status::Ok;
}

// --- Readouts --------------------------------------------------------------

Vec3 Phase1Calibrator::forward_axis(SourceId id) const {
    if (!configured_) return kFwd;
    const int s = slot_for(id);
    const Mat3 Rp = (s >= 0) ? prior_[s].R : Mat3::Identity();
    if (s < 0 || so3_[3 * s + 0].empty()) {
        // No data -> the prior forward axis (R_Xprior * e_x).
        return Rp * kFwd;
    }
    Vec3 phi;
    phi << so3_[3 * s + 0].mode(), so3_[3 * s + 1].mode(), so3_[3 * s + 2].mode();
    const Mat3 dR = so3::exp(phi);
    return dR * kFwd;
}

Scalar Phase1Calibrator::yaw(SourceId id) const {
    const Vec3 f = forward_axis(id);
    return std::atan2(f.y(), f.x());
}

Scalar Phase1Calibrator::pitch(SourceId id) const {
    const Vec3 f = forward_axis(id);
    return std::atan2(-f.z(), std::hypot(f.x(), f.y()));
}

SE3 Phase1Calibrator::extrinsic(SourceId id) const {
    if (!configured_) return SE3{};
    const int s = slot_for(id);
    SE3 X = (s >= 0) ? prior_[s] : SE3{};
    if (s < 0 || so3_[3 * s + 0].empty()) return X;   // prior
    Vec3 phi;
    phi << so3_[3 * s + 0].mode(), so3_[3 * s + 1].mode(), so3_[3 * s + 2].mode();
    const Mat3 dR = so3::exp(phi);
    // Corrected rotation = δR * R_Xprior (yaw/pitch corrected; roll + translation = prior).
    X.R = dR * X.R;
    return X;
}

Scalar Phase1Calibrator::scale(SourceId id) const {
    if (!configured_) return Scalar(1);
    if (id == reference_id_) return Scalar(1);
    const int s = slot_for(id);
    if (s < 0 || !scale_calib_[s] || scale_[s].empty()) return Scalar(1);
    return scale_[s].mode();
}

Scalar Phase1Calibrator::extrinsic_confidence(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0) return Scalar(0);
    // The weakest channel bounds the joint reliability.
    const Scalar cx = so3_[3 * s + 0].confidence();
    const Scalar cy = so3_[3 * s + 1].confidence();
    const Scalar cz = so3_[3 * s + 2].confidence();
    return std::min(cx, std::min(cy, cz));
}

Scalar Phase1Calibrator::scale_confidence(SourceId id) const {
    if (!configured_) return Scalar(0);
    if (id == reference_id_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0 || !scale_calib_[s]) return Scalar(0);
    return scale_[s].confidence();
}

Scalar Phase1Calibrator::vote_count(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0) return Scalar(0);
    return so3_[3 * s + 0].total();
}

// ===========================================================================
// Phase2Calibrator — turn-regime roll + xyz lever-arm via hand-eye. See
// calibration.hpp for the full derivation; this is the bounded, no-heap impl.
// ===========================================================================

namespace {
// LOAD-BEARING observability thresholds — all FIXED namespace constants (NOT config knobs;
// CONFIG lists none). They govern the rank-deficiency behaviour the slice is judged on and
// are deliberately not surfaced for tuning.
// Per-window guard: skip a row whose base-motion rotation magnitude ‖log R_A‖ is below
// this floor (the ω×r → 0 near-singular regime — adds an ill-conditioned ~zero row).
constexpr Scalar kRotRowMin = Scalar(1e-3);   // rad over the step
// Readout conditioning floor for the lever-arm LS: accept the FULL solve only when the
// smallest eigenvalue of AᵀA is ≥ kCondMin × the largest (else an axis is unobservable —
// e.g. z under pure yaw, whose AᵀA null space IS the z-axis).
constexpr Scalar kCondMin = Scalar(1e-3);
// Per-axis information floor for VOTING a channel: axis c is votable when its normal-
// matrix diagonal AᵀA(c,c) is ≥ kAxisInfoMin × the largest diagonal. An unobservable
// axis (near-zero diagonal) is below the floor and is NOT voted (its histogram stays
// empty -> low translation_confidence, the under-determined flag).
constexpr Scalar kAxisInfoMin = Scalar(1e-2);
// Ridge for the per-window voting solve (relative to the largest AᵀA diagonal). Small
// enough not to bias observable axes; large enough to pin an unobservable axis at prior.
constexpr Scalar kVoteRidgeRel = Scalar(1e-3);

// Strip the roll DOF from a rotation, returning its yaw/pitch part (the minimal rotation
// taking e_x -> the rotation's forward axis R·e_x). Roll is rotation about the forward
// axis, which leaves R·e_x fixed, so aligning e_x to that axis drops exactly the roll.
// Used as the Phase-2 roll basepoint fallback when Phase 1 has not yet supplied R_yp, so
// a nonzero PRIOR roll is not double-counted into the recovered roll on the bootstrap path.
Mat3 yaw_pitch_of(const Mat3& R) {
    const Vec3 f = R * kFwd;
    const Scalar fn = f.norm();
    if (fn < kUnitEps) return Mat3::Identity();
    const Vec3 fu = f / fn;
    const Vec3 v  = kFwd.cross(fu);
    const Scalar s = v.norm();             // sin(theta)
    const Scalar c = kFwd.dot(fu);         // cos(theta)
    if (s < kUnitEps) {
        // Forward axis already ±e_x: identity (parallel) — antiparallel is outside the
        // small-deviation regime (a near-π yaw); return identity defensively.
        return Mat3::Identity();
    }
    const Mat3 vx = so3::hat(v);
    return Mat3::Identity() + vx + vx * vx * ((Scalar(1) - c) / (s * s));
}
}  // namespace

Scalar Phase2Calibrator::vote_weight_factor(Scalar omega_norm, Scalar confidence) const {
    // D5 vote weighting. One -> 1; Rotation -> ‖ω‖ (> turn_omega_min by the gate, so
    // positive — the lever-arm benefits from more rotation); Confidence -> the source
    // Σ-confidence; Combo -> Rotation × Confidence. Always > 0 so add() keeps the vote.
    const Scalar rot  = std::max(kUnitEps, omega_norm);
    const Scalar conf = std::max(Scalar(0), confidence);
    switch (vote_weight_) {
        case VoteWeight::One:        return Scalar(1);
        case VoteWeight::Confidence: return std::max(kUnitEps, conf);
        case VoteWeight::Rotation:   return rot;
        case VoteWeight::Combo:      return std::max(kUnitEps, rot * conf);
    }
    return Scalar(1);
}

Mat3 Phase2Calibrator::roll_about_forward(Scalar roll) {
    // Rotation about the LOCAL forward axis (sensor x). Closed-form Rx(roll).
    const Scalar c = std::cos(roll), s = std::sin(roll);
    Mat3 R;
    R << Scalar(1), Scalar(0), Scalar(0),
         Scalar(0), c,        -s,
         Scalar(0), s,         c;
    return R;
}

Scalar Phase2Calibrator::rot_residual(const Mat3& R_yp, Scalar roll,
                                      const Mat3& R_A, const Mat3& R_B) {
    // R_X(roll) = R_yp · Rx(roll); residual = ‖log( R_X R_B R_X^T R_A^T )‖. At the true
    // roll this is the identity (R_A R_X = R_X R_B) ⇒ 0.
    const Mat3 R_X = R_yp * roll_about_forward(roll);
    const Mat3 E   = R_X * R_B * R_X.transpose() * R_A.transpose();
    return so3::log(E).norm();
}

Scalar Phase2Calibrator::best_roll(const Mat3& R_yp, const Mat3& R_A, const Mat3& R_B) {
    // Bounded 1-D minimizer over (−π, π]: a fixed-count grid scan, then a parabolic
    // refine of the grid minimizer using its two neighbours (wrap-aware). No heap, no
    // unbounded iteration.
    const Scalar step = (Scalar(2) * kPi) / static_cast<Scalar>(kRollGrid);
    int    best_k = 0;
    Scalar best_r = -kPi + step;       // grid covers (−π, π]
    Scalar best_v = rot_residual(R_yp, best_r, R_A, R_B);
    for (int k = 1; k < kRollGrid; ++k) {
        const Scalar r = -kPi + step * static_cast<Scalar>(k + 1);
        const Scalar v = rot_residual(R_yp, r, R_A, R_B);
        if (v < best_v) { best_v = v; best_r = r; best_k = k; }
    }
    // Parabolic refine using the residual at the minimizer's neighbours (period 2π).
    const Scalar rm = -kPi + step * static_cast<Scalar>(best_k);          // k-1 sample
    const Scalar rp = -kPi + step * static_cast<Scalar>(best_k + 2);      // k+1 sample
    const Scalar vm = rot_residual(R_yp, rm, R_A, R_B);
    const Scalar v0 = best_v;
    const Scalar vp = rot_residual(R_yp, rp, R_A, R_B);
    const Scalar denom = vm - Scalar(2) * v0 + vp;
    Scalar r = best_r;
    if (std::abs(denom) > Scalar(1e-12)) {
        Scalar off = Scalar(0.5) * (vm - vp) / denom;     // in (−0.5, 0.5) for a clean min
        if (off >  Scalar(0.5)) off =  Scalar(0.5);
        if (off < -Scalar(0.5)) off = -Scalar(0.5);
        r = best_r + off * step;
    }
    // Wrap into (−π, π].
    while (r <= -kPi) r += Scalar(2) * kPi;
    while (r >   kPi) r -= Scalar(2) * kPi;
    return r;
}

Vec3 Phase2Calibrator::solve_ridge(const Mat3& AtA, const Vec3& Atb,
                                   const Vec3& t_prior, Scalar ridge, bool obs[3]) {
    // Largest diagonal (the information scale). Pure PSD normal matrix -> diagonals >= 0.
    Scalar max_diag = Scalar(0);
    for (int c = 0; c < 3; ++c) max_diag = std::max(max_diag, AtA(c, c));
    for (int c = 0; c < 3; ++c) {
        obs[c] = (max_diag > Scalar(0)) && (AtA(c, c) >= kAxisInfoMin * max_diag);
    }
    if (!(max_diag > Scalar(0))) return t_prior;     // no information at all

    // Ridge centered on the prior: (AtA + r I)(t - t_prior) = Atb - AtA t_prior.
    const Mat3 M = AtA + (ridge * max_diag) * Mat3::Identity();
    const Vec3 rhs = Atb - AtA * t_prior;
    const Vec3 dt  = M.ldlt().solve(rhs);
    Vec3 t = t_prior + dt;
    if (!t.allFinite()) return t_prior;
    return t;
}

Status Phase2Calibrator::configure(const Config& cfg, SourceId reference_id) {
    if (static_cast<int>(reference_id) >= kMaxSources) return Status::OutOfRange;

    roll_cfg_ = cfg.roll_hist;
    xyz_cfg_  = cfg.xyz_hist;

    for (int i = 0; i < kMaxSources; ++i) {
        const Status hs = roll_[i].configure(roll_cfg_);
        if (!ok(hs)) return hs;
    }
    for (int i = 0; i < 3 * kMaxSources; ++i) {
        const Status hs = xyz_[i].configure(xyz_cfg_);
        if (!ok(hs)) return hs;
    }

    reference_id_   = reference_id;
    turn_omega_min_ = std::max(Scalar(0), cfg.turn_omega_min);
    strategy_       = cfg.phase2_strategy;
    vote_weight_    = cfg.vote_weight;
    last_gate_turning_ = false;

    configured_ = true;
    reset();
    return Status::Ok;
}

void Phase2Calibrator::reset() {
    for (int i = 0; i < kMaxSources; ++i) {
        roll_[i].reset();
        prior_[i]   = SE3{};
        ryp_[i]     = Mat3::Identity();
        ata_[i]     = Mat3::Zero();
        atb_[i]     = Vec3::Zero();
        rows_[i]    = Scalar(0);
        ryp_set_[i] = false;
        ids_[i]     = 0;
    }
    for (int i = 0; i < 3 * kMaxSources; ++i) xyz_[i].reset();
    source_count_      = 0;
    last_gate_turning_ = false;
}

int Phase2Calibrator::slot_for(SourceId id) const {
    for (int i = 0; i < source_count_; ++i) {
        if (ids_[i] == id) return i;
    }
    return -1;
}

int Phase2Calibrator::ensure_slot(SourceId id) {
    const int s = slot_for(id);
    if (s >= 0) return s;
    if (source_count_ >= kMaxSources) return -1;
    const int slot = source_count_++;
    ids_[slot]     = id;
    prior_[slot]   = SE3{};
    ryp_[slot]     = Mat3::Identity();
    ata_[slot]     = Mat3::Zero();
    atb_[slot]     = Vec3::Zero();
    rows_[slot]    = Scalar(0);
    ryp_set_[slot] = false;
    return slot;
}

Status Phase2Calibrator::set_prior(SourceId id, const SE3& prior_extrinsic) {
    if (!configured_)                        return Status::NotInitialized;
    if (static_cast<int>(id) >= kMaxSources) return Status::OutOfRange;
    const int slot = ensure_slot(id);
    if (slot < 0)                            return Status::CapacityExceeded;
    prior_[slot] = prior_extrinsic;
    return Status::Ok;
}

Status Phase2Calibrator::set_basepoint(SourceId id, const SE3& basepoint) {
    // Re-anchor the lever-arm prior + the pinned-ref gauge to a committed extrinsic
    // (Slice 8). Same validation / slot semantics as set_prior.
    if (!configured_)                        return Status::NotInitialized;
    if (static_cast<int>(id) >= kMaxSources) return Status::OutOfRange;
    const int slot = ensure_slot(id);
    if (slot < 0)                            return Status::CapacityExceeded;
    prior_[slot] = basepoint;
    return Status::Ok;
}

void Phase2Calibrator::reset_roll(SourceId id) {
    if (!configured_) return;
    const int s = slot_for(id);
    if (s < 0) return;
    roll_[s].reset();
}

void Phase2Calibrator::reset_lever(SourceId id) {
    if (!configured_) return;
    const int s = slot_for(id);
    if (s < 0) return;
    ata_[s]  = Mat3::Zero();
    atb_[s]  = Vec3::Zero();
    rows_[s] = Scalar(0);
    xyz_[3 * s + 0].reset();
    xyz_[3 * s + 1].reset();
    xyz_[3 * s + 2].reset();
}

Status Phase2Calibrator::set_yaw_pitch(SourceId id, const Mat3& R_yp) {
    if (!configured_)                        return Status::NotInitialized;
    if (static_cast<int>(id) >= kMaxSources) return Status::OutOfRange;
    const int slot = ensure_slot(id);
    if (slot < 0)                            return Status::CapacityExceeded;
    ryp_[slot]     = R_yp;
    ryp_set_[slot] = true;
    return Status::Ok;
}

Status Phase2Calibrator::observe(int n, const SourceId* ids, const SE3* reported,
                                 const SE3& fused_motion, const Vec3& fused_omega,
                                 const Scalar* confidences) {
    if (!configured_)                                    return Status::NotInitialized;
    if (n <= 0 || ids == nullptr || reported == nullptr) return Status::NoData;

    // --- Turn gate (the observability spine, D5) ---------------------------------
    // Roll + lever arm are observable ONLY under rotation. Accept only sufficiently
    // turning fused motion. Gated-out steps are a no-op (not an error).
    const Scalar omega_norm = fused_omega.norm();
    const bool turning = (omega_norm > turn_omega_min_);
    last_gate_turning_ = turning;
    if (!turning) return Status::NotReady;

    // Base motion A per strategy.
    //   VsFusedBase       : A = the fused consensus motion (passed in).
    //   PairwisePinnedRef : A = X_ref ∘ B_ref ∘ X_ref^{-1}, the motion implied by the
    //                       pinned reference (extrinsic held at the prior — the gauge).
    //                       Solving the same hand-eye then yields each source's X relative
    //                       to the reference. If the reference is absent this step, we
    //                       cannot form A_impl -> no votes (NotReady).
    SE3 A_base = fused_motion;
    if (strategy_ == Phase2Strat::PairwisePinnedRef) {
        int ref_idx = -1;
        for (int i = 0; i < n; ++i) {
            if (ids[i] == reference_id_) { ref_idx = i; break; }
        }
        if (ref_idx < 0) return Status::NotReady;
        const int rs = slot_for(reference_id_);
        const SE3 Xref = (rs >= 0) ? prior_[rs] : SE3{};
        // A_impl = Xref o B_ref o Xref^{-1}.
        A_base = se3::compose(se3::compose(Xref, reported[ref_idx]), se3::inverse(Xref));
    }

    const Mat3& R_A = A_base.R;
    const Vec3& t_A = A_base.t;
    const Scalar rotA = so3::log(R_A).norm();
    // Per-window near-singular guard: a low-rotation A is the ω×r→0 regime — its
    // translation row is near-zero and ill-conditioned. Skip the whole window's votes.
    if (rotA < kRotRowMin) return Status::NotReady;
    const Mat3 RA_minus_I = R_A - Mat3::Identity();

    bool any = false;
    for (int i = 0; i < n; ++i) {
        const SourceId id = ids[i];
        // In PairwisePinnedRef the reference solves to identity-relative (Y = I); its
        // own hand-eye against A_impl is degenerate (B_ref maps to A_impl exactly), so
        // skip voting the reference. In VsFusedBase the reference is voted too (its X
        // should stay near its prior).
        if (strategy_ == Phase2Strat::PairwisePinnedRef && id == reference_id_) continue;

        const int slot = ensure_slot(id);
        if (slot < 0) continue;                       // at capacity

        const Mat3& R_B = reported[i].R;
        const Vec3& t_B = reported[i].t;

        // Roll basepoint R_yp: the Phase-1 yaw/pitch rotation if provided. FALLBACK
        // (bootstrap, before Phase 1 has fed an R_yp — see set_yaw_pitch): the prior's
        // yaw/pitch ONLY, with its prior roll STRIPPED, so the recovered roll is measured
        // about the same forward axis and does NOT double-count any prior roll. We strip
        // roll by taking the prior forward axis f = R_prior·e_x and building the minimal
        // yaw/pitch rotation that aligns e_x -> f (roll about forward leaves f fixed, so
        // this drops exactly the roll DOF). In the WIRED estimator R_yp is set every step
        // (the wired path always sets R_yp), so this fallback only matters on the
        // bootstrap / Slice-8 path — but it is now roll-safe there too.
        const Mat3 R_yp = ryp_set_[slot] ? ryp_[slot] : yaw_pitch_of(prior_[slot].R);

        // Per-source vote weight (D5): the turn magnitude and/or this source's
        // Σ-confidence. `confidences` optional (null -> unit confidence).
        const Scalar conf = (confidences != nullptr) ? confidences[i] : Scalar(1);
        const Scalar w    = vote_weight_factor(omega_norm, conf);

        // --- Roll: 1-D rotation hand-eye residual minimizer, voted into S¹ -------
        const Scalar roll = best_roll(R_yp, R_A, R_B);
        roll_[slot].add(roll, w);

        // --- xyz: accumulate the LS row (R_A − I) t_X = R_X t_B − t_A ------------
        // Use the full R_X at the just-recovered roll for the translation row (the row
        // depends on the rotation; using the per-window roll keeps it consistent).
        const Mat3 R_X = R_yp * roll_about_forward(roll);
        const Vec3 b   = R_X * t_B - t_A;
        // Incremental normal equations (fixed 3×3 + 3×1; no growing matrix).
        ata_[slot].noalias() += RA_minus_I.transpose() * RA_minus_I;
        atb_[slot]           += RA_minus_I.transpose() * b;
        rows_[slot]          += Scalar(1);

        // Vote the RUNNING ridge-regularized solution into the per-channel histograms.
        // The running solve converges as rows accumulate; the histogram mode then gives a
        // robust (outlier-window-resistant) committed estimate on top of the LS. Only an
        // OBSERVABLE axis (enough information in AᵀA) is voted — an unobservable axis (e.g.
        // z under pure yaw) stays empty -> low translation_confidence (the spine: the
        // lever arm needs rotation, per-axis).
        bool obs[3];
        const Vec3 tx_run = solve_ridge(ata_[slot], atb_[slot], prior_[slot].t,
                                        kVoteRidgeRel, obs);
        if (obs[0]) xyz_[3 * slot + 0].add(tx_run.x(), w);
        if (obs[1]) xyz_[3 * slot + 1].add(tx_run.y(), w);
        if (obs[2]) xyz_[3 * slot + 2].add(tx_run.z(), w);

        any = true;
    }

    return any ? Status::Ok : Status::NotReady;
}

// --- Readouts --------------------------------------------------------------

Vec3 Phase2Calibrator::solve_lever_arm(SourceId id, bool* observable) const {
    if (observable) *observable = false;
    if (!configured_) return Vec3::Zero();
    const int s = slot_for(id);
    if (s < 0) return Vec3::Zero();
    const SE3 prior = prior_[s];
    if (rows_[s] < Scalar(1)) return prior.t;        // no rows -> prior

    const Mat3& AtA = ata_[s];
    // Conditioning guard via the symmetric eigenvalues (AtA is symmetric PSD). If the
    // smallest is below kCondMin × the largest, an axis is unobservable (e.g. z under
    // pure yaw) -> report the prior rather than a blown-up solve.
    Eigen::SelfAdjointEigenSolver<Mat3> es(AtA);
    if (es.info() != Eigen::Success) return prior.t;
    const Vec3 ev = es.eigenvalues();                // ascending
    const Scalar lo = ev(0), hi = ev(2);
    // Conditioning floor is INCLUSIVE: an exactly-at-floor smallest eigenvalue
    // (lo == kCondMin·hi) is accepted as observable. Below the floor -> rank-deficient
    // (e.g. z under pure yaw) -> report the prior rather than a blown-up solve.
    if (!(hi > Scalar(0)) || lo < kCondMin * hi) {
        return prior.t;                              // ill-conditioned / rank-deficient
    }
    const Vec3 tx = AtA.ldlt().solve(atb_[s]);
    if (!tx.allFinite()) return prior.t;
    if (observable) *observable = true;
    return tx;
}

Scalar Phase2Calibrator::roll(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0 || roll_[s].empty()) return Scalar(0);    // prior roll = 0
    return roll_[s].mode();
}

Vec3 Phase2Calibrator::lever_arm(SourceId id) const {
    if (!configured_) return Vec3::Zero();
    const int s = slot_for(id);
    if (s < 0) return Vec3::Zero();
    const SE3 prior = prior_[s];
    // Per-channel histogram mode is the robust committed estimate; an empty channel
    // (unobservable / unvoted axis) falls back to the prior component.
    Vec3 t = prior.t;
    if (!xyz_[3 * s + 0].empty()) t.x() = xyz_[3 * s + 0].mode();
    if (!xyz_[3 * s + 1].empty()) t.y() = xyz_[3 * s + 1].mode();
    if (!xyz_[3 * s + 2].empty()) t.z() = xyz_[3 * s + 2].mode();
    return t;
}

SE3 Phase2Calibrator::extrinsic(SourceId id) const {
    if (!configured_) return SE3{};
    const int s = slot_for(id);
    if (s < 0) return SE3{};
    SE3 X = prior_[s];
    const Mat3 R_yp = ryp_set_[s] ? ryp_[s] : prior_[s].R;
    X.R = R_yp * roll_about_forward(roll(id));    // yaw/pitch ∘ recovered roll
    X.t = lever_arm(id);
    return X;
}

Scalar Phase2Calibrator::extrinsic_confidence(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0) return Scalar(0);
    return roll_[s].confidence();
}

Scalar Phase2Calibrator::translation_confidence(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0) return Scalar(0);
    const Scalar cx = xyz_[3 * s + 0].confidence();
    const Scalar cy = xyz_[3 * s + 1].confidence();
    const Scalar cz = xyz_[3 * s + 2].confidence();
    return std::min(cx, std::min(cy, cz));
}

Scalar Phase2Calibrator::roll_vote_count(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0) return Scalar(0);
    return roll_[s].total();
}

Scalar Phase2Calibrator::xyz_vote_count(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0) return Scalar(0);
    return rows_[s];
}

} // namespace ofc
