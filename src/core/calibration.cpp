// ofc/core/calibration.cpp — Phase-1 (straight-regime) calibration. See calibration.hpp
// for the full direction->yaw/pitch mapping, the scale derivation, the gate, and the
// reverse-fold; this file is the bounded, no-heap implementation.
//
// STRICT CORE: all storage is fixed-capacity member arrays bound in configure();
// set_prior()/observe()/readouts allocate nothing and run loops bounded by the live
// source count or the channel count. No exceptions; status-code returns; double math.
#include "ofc/core/calibration.hpp"

#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
constexpr Scalar kUnitEps = Scalar(1e-12);   // guard ||t|| ~ 0 before normalizing
constexpr Scalar kPi      = Scalar(3.14159265358979323846);

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

Status Phase1Calibrator::observe(int n, const SourceId* ids, const SE3* reported,
                                 const Vec3& fused_omega, const Vec3& fused_trans) {
    if (!configured_)                       return Status::NotInitialized;
    if (n <= 0 || ids == nullptr || reported == nullptr) return Status::NoData;

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

        // Candidate rotation δR taking e_x -> g_unit; vote its so(3) log (3 channels).
        // The guard above keeps θ < 90°, well clear of the so3::log singularity at π.
        const Mat3 dR  = rotation_between(kFwd, g_unit);
        const Vec3 phi = so3::log(dR);
        so3_[3 * slot + 0].add(phi.x(), Scalar(1));
        so3_[3 * slot + 1].add(phi.y(), Scalar(1));
        so3_[3 * slot + 2].add(phi.z(), Scalar(1));

        // Scale vote = magnitude ratio vs the reference's reported magnitude.
        if (have_ref && scale_calib_[slot] && id != reference_id_) {
            scale_[slot].add(bn / ref_mag, Scalar(1));
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

} // namespace ofc
