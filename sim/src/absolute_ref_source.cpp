// ofc_sim/absolute_ref_source.cpp — see absolute_ref_source.hpp for the model.
#include "ofc_sim/absolute_ref_source.hpp"

#include "ofc/core/lie.hpp"

#include <cassert>
#include <cmath>
#include <random>

namespace ofc {
namespace sim {

namespace {
constexpr Scalar kNanosPerSec = Scalar(1e9);

// Deterministic per-measurement seed: a pure function of (source seed, frontier stamp), so
// the noise on a given fix does not depend on call history (replay-stable, mirrors
// SyntheticSource::mix_seed).
std::uint64_t mix_seed(std::uint32_t base, Timestamp t) {
    std::uint64_t z = static_cast<std::uint64_t>(base);
    z ^= static_cast<std::uint64_t>(t) + 0x9E3779B97F4A7C15ULL + (z << 6) + (z >> 2);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z =  z ^ (z >> 31);
    return z;
}
} // namespace

bool SyntheticAbsoluteRef::evaluate(const State& x, Measurement& out) const {
    if (traj_ == nullptr) return false;
    const Timestamp stamp = x.stamp;

    // --- rate gate: emit at most once per period_s of frontier time --------------------
    if (p_.period_s > Scalar(0) && have_last_emit_) {
        const Scalar dt_s = static_cast<Scalar>(stamp - last_emit_) / kNanosPerSec;
        if (dt_s < p_.period_s) return false;
    }

    // --- GT anchor (lazy, mirrors rig.cpp) ---------------------------------------------
    // The estimator's odom origin sits at the GT pose at (first_frontier - window_s): it
    // starts the odom pose at identity at the first fused tick and composes the first
    // bootstrap window. Establish anchor_inv on the FIRST evaluate (== the first frontier).
    if (!have_anchor_) {
        // FOOTGUN GUARD (relaxed-edge sim): window_s MUST mirror the rig's Config::window_s,
        // else this anchor sits at a different GT pose than the estimator's odom origin and
        // EVERY residual is biased by a constant offset (a silent constant-offset bias). We
        // cannot read the rig's Config from here (evaluate() only gets State x), so we assert
        // the value is at least non-negative and finite as cheap insurance; callers are
        // responsible for setting rp.window_s == cfg.window_s (see the header field comment).
        assert(p_.window_s >= Scalar(0) && std::isfinite(p_.window_s) &&
               "AbsoluteRefParams::window_s must mirror the rig's Config::window_s "
               "(non-negative); a mismatch silently biases every correction residual");
        const Timestamp window_ns =
            static_cast<Timestamp>(std::llround(p_.window_s * kNanosPerSec));
        anchor_inv_  = se3::inverse(traj_->pose(stamp - window_ns));
        have_anchor_ = true;
    }

    // --- true GT position in the odom frame --------------------------------------------
    const SE3 gt_odom = se3::compose(anchor_inv_, traj_->pose(stamp));
    Vec3 z = gt_odom.t;

    // Outlier injection (frontier time): a gross-wrong position so the gate must reject it.
    const Scalar t_s = static_cast<Scalar>(stamp) / kNanosPerSec;
    bool is_outlier = false;
    for (const RefWindow& w : p_.outlier_windows) {
        if (w.contains(t_s)) { is_outlier = true; break; }
    }
    if (is_outlier) {
        z += p_.outlier_offset_m;
    } else if (p_.sigma_pos > Scalar(0)) {
        // Seeded zero-mean Gaussian noise consistent with R = sigma_pos^2 I3.
        std::mt19937_64 gen(mix_seed(p_.seed, stamp));
        std::normal_distribution<Scalar> nz(Scalar(0), p_.sigma_pos);
        const Scalar e0 = nz(gen);
        const Scalar e1 = nz(gen);
        const Scalar e2 = nz(gen);
        z += Vec3(e0, e1, e2);
    }

    // --- measurement: residual = z - h(x), h(x) = x.pose.t -----------------------------
    out.dim = 3;
    out.residual.setZero();
    out.residual.head<3>() = z - x.pose.t;
    // H = [ R | 0_3x3 | 0_3x6 ] (right-error tangent: t perturbs as t + R*rho to 1st order).
    out.H.setZero();
    out.H.block<3, 3>(0, 0) = x.pose.R;
    out.R.setZero();
    const Scalar var = (p_.sigma_pos > Scalar(0)) ? (p_.sigma_pos * p_.sigma_pos)
                                                  : Scalar(1e-6);   // floor for a clean fix
    out.R.block<3, 3>(0, 0) = var * Mat3::Identity();
    out.stamp = stamp;

    last_emit_      = stamp;
    have_last_emit_ = true;
    return true;
}

} // namespace sim
} // namespace ofc
