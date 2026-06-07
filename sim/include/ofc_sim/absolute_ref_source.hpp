// ofc_sim/absolute_ref_source.hpp — a synthetic ICorrection: a GPS-like POSITION fix
// against the GT trajectory (the oracle), expressed in the estimator's odom frame.
//
// RELAXED EDGE (sim/): std containers / std::mt19937 / exceptions / mutable lazy state are
// fine here. This source is the test fixture that exercises the Slice-11 correction path:
// it returns a linearized position measurement that the ESKF Mahalanobis-gates + applies,
// so an absolute reference removes the predict-only odometry drift.
//
// MEASUREMENT MODEL (position fix, dim=3). h(x) = x.pose.t (the fused base position in the
// odom frame). The reference reports the TRUE GT position in that SAME odom frame:
//       z = ( anchor_inv o GT_pose(stamp) ).t
// where anchor_inv = GT_pose(first_frontier - window_s)^{-1} mirrors the rig's GT anchor
// (rig.cpp): the estimator starts its odom pose at identity at the first fused tick and
// composes the first bootstrap window, so its pose origin sits at the GT pose at
// (first_frontier - window_s). We establish anchor_inv LAZILY on the first evaluate() (the
// first frontier), so the ref needs no external wiring — it self-anchors exactly as the
// rig does. residual = z (-) h(x) = z - x.pose.t.
//
// JACOBIAN (H, 3x12). In the ESKF's RIGHT-error tangent T_true = T o exp(eta),
// eta = [rho; theta] (trans;rot), the odom-frame translation perturbs as
//       t(eta) = t + R * J_l(theta) * rho  ~=  t + R * rho   (first order, theta -> 0)
// so dh/d eta = [ R | 0_3x3 ] on the pose block and 0 on the twist block:
//       H = [ R | 0_3x3 | 0_3x6 ]   (3x12),   R = x.pose.R.
// (Derived/verified by the drift-removal test: a wrong H would not reduce the drift.)
//
// NOISE (R). R_noise = sigma_pos^2 * I3. Seeded zero-mean Gaussian noise is added to z
// consistent with R, deterministically per measurement (mirrors SyntheticSource: the draw
// is a pure function of (seed, stamp), so replay is stable regardless of tick alignment).
//
// RATE / OUTLIERS. A fix is emitted only every `period_ns` of frontier time (a configurable
// measurement rate, not every step). Outlier windows (frontier time) replace z with a
// gross-wrong position so the Mahalanobis gate is exercised.
//
// SCOPE: position-only (dim=3). A full pose fix (dim=6) or orientation fix (dim=3) is a
// noted extension — same machinery, a different h/H — not required for Slice 11.
#ifndef OFC_SIM_ABSOLUTE_REF_SOURCE_HPP
#define OFC_SIM_ABSOLUTE_REF_SOURCE_HPP

#include "ofc/core/absolute_ref.hpp"
#include "ofc/core/types.hpp"

#include "ofc_sim/trajectory.hpp"

#include <cstdint>
#include <vector>

namespace ofc {
namespace sim {

// A half-open frontier-time window [start_s, end_s) for injecting outlier fixes.
struct RefWindow {
    Scalar start_s = Scalar(0);
    Scalar end_s   = Scalar(0);
    bool contains(Scalar t_s) const { return t_s >= start_s && t_s < end_s; }
};

struct AbsoluteRefParams {
    // Measurement rate: emit a fix once per this many seconds of frontier time. The first
    // eligible frontier emits; subsequent emits are gated to >= this spacing. <= 0 => every
    // frontier emits.
    Scalar period_s = Scalar(0.2);

    // window_s mirrors the rig/estimator GT anchor: the estimator's odom origin sits at the
    // GT pose at (first_frontier - window_s). FOOTGUN: this MUST equal the Config::window_s
    // the rig was init'd with. If it does not, anchor_inv shifts the reference origin and EVERY
    // residual is biased by a constant offset (a silent constant-offset bias — the drift test
    // would still partially improve, masking the misconfig). Callers must set
    // `rp.window_s = cfg.window_s`. The first evaluate() asserts this value is sane (see .cpp).
    Scalar window_s = Scalar(0.10);

    // Position noise (m). R = sigma_pos^2 * I3; the added measurement noise is drawn with
    // this sigma per axis. 0 => noise-free (a perfect fix).
    Scalar sigma_pos = Scalar(0.05);

    std::uint32_t seed = 0;   // fixed seed (reproducible)

    // Outlier-fix windows (frontier time). A fix whose frontier time falls in one of these
    // is replaced by `outlier_offset_m` added to the true position (a gross wrong fix), so
    // the Mahalanobis gate must reject it.
    std::vector<RefWindow> outlier_windows;
    Vec3 outlier_offset_m = Vec3(20.0, 0.0, 0.0);
};

class SyntheticAbsoluteRef : public ICorrection {
public:
    SyntheticAbsoluteRef(const Trajectory& traj, const AbsoluteRefParams& p)
        : traj_(&traj), p_(p) {}

    // ICorrection: build a position measurement (residual + H + R) at the current frontier
    // state x. Returns false when no fix is due this step (rate gate). x.stamp is the real
    // frontier time the estimator stamps in (see estimator.cpp's correction loop).
    bool evaluate(const State& x, Measurement& out) const override;

    const AbsoluteRefParams& params() const { return p_; }

private:
    const Trajectory* traj_;
    AbsoluteRefParams p_;

    // Lazy GT anchor (mirrors rig.cpp) + rate state. mutable: evaluate() is const but the
    // anchor/last-emit are established on the fly (relaxed-edge sim convenience).
    mutable bool      have_anchor_ = false;
    mutable SE3       anchor_inv_;             // GT_pose(first_frontier - window_s)^{-1}
    mutable bool      have_last_emit_ = false;
    mutable Timestamp last_emit_ = 0;
};

} // namespace sim
} // namespace ofc
#endif // OFC_SIM_ABSOLUTE_REF_SOURCE_HPP
