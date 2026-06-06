// ofc/core/eskf.cpp — predict-only error-state KF integrator (Slice 2).
//
// Error-state ordering (matches the codebase [trans; rot] convention, DECISIONS D21):
//   0..2  pose translation error   3..5  pose rotation error
//   6..8  twist linear (v)         9..11 twist angular (omega)
//
// Right-error propagation. With the pose updated by a RIGHT composition
// T+ = T o delta, a right-invariant error eta defined by T_true = T o exp(eta)
// transports as eta+ = Ad(delta^-1) eta. So the 12x12 propagation Jacobian is
//   F = [[ Ad(delta^-1), 0 ],
//        [ 0,            0 ]]
// the twist block is zero because the twist is RE-READ each step (log(delta)/dt),
// not integrated, so its prior error does not carry into the next step. Covariance:
//   P+ = F P F^T + Qmap,   Qmap = blkdiag( q_pose, q_pose / dt^2 ).
// q_pose is the 6x6 process noise on the pose increment ([trans; rot]); dividing by
// dt^2 maps that delta uncertainty to the velocity (twist) uncertainty. The result
// is symmetrized to kill round-off asymmetry and stays PSD (F P F^T is PSD for any
// F, and Qmap is PSD).
#include "ofc/core/eskf.hpp"

#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
constexpr Scalar kNanosPerSec = Scalar(1e9);
constexpr Scalar kMinDt       = Scalar(1e-9);   // guard against dt <= 0

Timestamp secs_to_ns(Scalar s) {
    return static_cast<Timestamp>(std::llround(s * kNanosPerSec));
}

Mat12 symmetrize(const Mat12& M) {
    return Scalar(0.5) * (M + M.transpose());
}
} // namespace

void Eskf::init(const SE3& pose0, const Mat12& cov0) {
    state_ = State{};
    state_.pose      = pose0;
    state_.twist.xi  = Vec6::Zero();
    state_.twist.cov = Mat6::Identity();
    state_.cov       = symmetrize(cov0);   // ensure clean symmetry from the start
    state_.stamp     = 0;
}

void Eskf::predict(const SE3& delta, Scalar dt, const Mat6& q_pose) {
    const Scalar dt_eff = (dt > kMinDt) ? dt : kMinDt;

    // --- mean propagation ---------------------------------------------------
    state_.pose     = se3::compose(state_.pose, delta);   // T <- T o delta
    state_.twist.xi = se3::log(delta) / dt_eff;           // const-velocity readout

    // --- covariance propagation: P <- F P F^T + Qmap ------------------------
    // Pose block of F is the inverse adjoint (right-error of the right composition).
    const Mat6 Ad_inv = se3::adjoint(se3::inverse(delta));

    Mat12 F = Mat12::Zero();
    F.block<6, 6>(0, 0) = Ad_inv;
    // twist block stays zero: twist error is reset by the fresh log(delta)/dt read.

    Mat12 Q = Mat12::Zero();
    Q.block<6, 6>(0, 0) = q_pose;                         // pose increment noise
    Q.block<6, 6>(6, 6) = q_pose / (dt_eff * dt_eff);     // -> velocity noise

    state_.cov = symmetrize(F * state_.cov * F.transpose() + Q);

    // Keep the published per-twist covariance consistent with the 12x12 block.
    state_.twist.cov = state_.cov.block<6, 6>(6, 6);

    state_.stamp += secs_to_ns(dt_eff);
}

void Eskf::predict_tip(Scalar dt_ahead, Scalar inflation, State& tip_out) const {
    const Scalar infl = std::max(Scalar(1), inflation);
    const Scalar dta  = std::max(Scalar(0), dt_ahead);

    // Const-velocity extrapolation: advance by exp(twist * dt_ahead).
    const SE3 ext = se3::exp(state_.twist.xi * dta);
    tip_out.pose      = se3::compose(state_.pose, ext);
    tip_out.twist     = state_.twist;
    tip_out.cov       = symmetrize(infl * state_.cov);   // inflation >= 1 keeps PSD
    tip_out.twist.cov = tip_out.cov.block<6, 6>(6, 6);
    tip_out.stamp     = state_.stamp + secs_to_ns(dta);
}

Mat6 Eskf::adaptive_q(Scalar spread, Scalar q_scale, const Scalar* q_floor) {
    Mat6 Q = (q_scale * spread * spread) * Mat6::Identity();
    if (q_floor != nullptr) {
        for (int i = 0; i < 6; ++i) Q(i, i) += q_floor[i];
    }
    return Q;
}

} // namespace ofc
