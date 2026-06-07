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

#include <Eigen/Cholesky>

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
    last_nis_        = Scalar(0);          // no update() yet
}

void Eskf::predict(const SE3& delta, Scalar dt, const Mat6& q_pose) {
    const Scalar dt_eff = (dt > kMinDt) ? dt : kMinDt;

    // --- mean propagation ---------------------------------------------------
    state_.pose     = se3::compose(state_.pose, delta);   // T <- T o delta
    state_.twist.xi = se3::log(delta) / dt_eff;           // const-velocity readout
    // ACCEPTED SLICE-2 APPROXIMATION: the published twist.cov below is q_pose/dt^2
    // only (the fresh-readout noise). It ignores that delta's own posterior pose
    // covariance also feeds the readout. This is consistent with the "re-read the
    // twist each step" model; a tighter coupling lands with the later reliability
    // work — recorded as a modeling approximation, not a bug.

    // --- covariance propagation: P <- F P F^T + Qmap ------------------------
    // Pose block of F is the inverse adjoint (right-error of the right composition).
    const Mat6 Ad_inv = se3::adjoint(se3::inverse(delta));

    Mat12 F = Mat12::Zero();
    F.block<6, 6>(0, 0) = Ad_inv;
    // twist block stays zero: twist error is reset by the fresh log(delta)/dt read.
    // ACCEPTED SLICE-2 APPROXIMATION: F and Qmap are block-diagonal, so the
    // pose<->twist cross-covariance is forced to 0 every predict even though
    // twist = log(delta)/dt is a deterministic function of the same delta that moves
    // the pose (they ARE correlated). Tip extrapolation (pose+twist) is therefore
    // mildly over/under-confident. A future correct model couples them via the
    // readout Jacobian (reliability/coupling improve in a later slice).

    Mat12 Q = Mat12::Zero();
    Q.block<6, 6>(0, 0) = q_pose;                         // pose increment noise
    Q.block<6, 6>(6, 6) = q_pose / (dt_eff * dt_eff);     // -> velocity noise

    state_.cov = symmetrize(F * state_.cov * F.transpose() + Q);

    // Keep the published per-twist covariance consistent with the 12x12 block.
    state_.twist.cov = state_.cov.block<6, 6>(6, 6);

    // NOTE: the filter does NOT maintain an absolute timeline. The estimator owns
    // the frontier clock and stamps the published state with the real frontier t1
    // (and the tip with `now`); accumulating dt here from 0 would diverge from that
    // real timeline, so we deliberately leave state_.stamp untouched in predict()
    // (review fix — was `state_.stamp += dt_eff`, which was masked but inconsistent).
}

bool Eskf::update(const Measurement& m, Scalar chi2_threshold) {
    const int n = m.dim;
    if (n <= 0 || n > 6) return false;       // empty / out-of-range measurement => no-op

    // Work in FULL fixed-size 6/12 capacity (strict core: no heap). The measurement is
    // ACTIVE only on its top-left n block; we PAD the unused rows so every fixed-size
    // solve below is well-posed yet the padding contributes nothing:
    //   * H rows  [n..5] -> 0  => P H^T has zero columns there => zero gain columns.
    //   * residual[n..5] -> 0  => zero NIS contribution, zero error injection.
    //   * R / S  padded with 1 on the unused diagonal => S stays invertible (the padded
    //     1x1 blocks are independent identities; with a zero residual + zero H they do not
    //     couple into the active n-block result).
    Eigen::Matrix<Scalar, 6, 12> H = Eigen::Matrix<Scalar, 6, 12>::Zero();
    Vec6 r = Vec6::Zero();
    Mat6 R = Mat6::Identity();               // unused diagonal already 1 (padding)
    H.topRows(n)    = m.H.topRows(n);
    r.head(n)       = m.residual.head(n);
    R.topLeftCorner(n, n) = m.R.topLeftCorner(n, n);

    const Mat12& P = state_.cov;

    // Innovation covariance S = H P H^T + R (active n-block; padded rows -> identity).
    const Mat6 S = H * P * H.transpose() + R;

    // NIS d2 = r^T S^-1 r. Stable solve (no explicit inverse). Padding rows add 0 (r=0).
    const Eigen::LDLT<Mat6> Sldlt = S.ldlt();
    const Vec6  Sinv_r = Sldlt.solve(r);
    const Scalar d2    = r.dot(Sinv_r);
    last_nis_ = d2;                           // exposed even when the gate rejects

    // Mahalanobis gate: reject (state unchanged) when NIS exceeds the threshold.
    if (d2 > chi2_threshold) return false;

    // Kalman gain K = P H^T S^-1  (12 x 6 capacity; unused columns are 0 via padded H).
    // Solve S K^T = (P H^T)^T column-wise: K^T = S^-1 (H P) -> K = (S^-1 (H P))^T.
    const Eigen::Matrix<Scalar, 6, 12> HP = H * P;            // 6 x 12
    const Eigen::Matrix<Scalar, 6, 12> Kt = Sldlt.solve(HP);  // S^-1 (H P) = K^T
    const Eigen::Matrix<Scalar, 12, 6> K  = Kt.transpose();   // 12 x 6

    // Error state dx = K r (12-vector).
    const Eigen::Matrix<Scalar, 12, 1> dx = K * r;

    // Inject into the nominal state in the right-error tangent (matches predict()):
    //   pose  <- pose o exp(dx[0..5])   (full SE(3) exp of the [trans;rot] tangent — the
    //                                    same coupled-SE(3) tangent P is propagated in)
    //   twist <- twist + dx[6..11]      (additive)
    const Vec6 dx_pose  = dx.head<6>();
    const Vec6 dx_twist = dx.tail<6>();
    state_.pose     = se3::compose(state_.pose, se3::exp(dx_pose));
    state_.twist.xi = state_.twist.xi + dx_twist;

    // Covariance: Joseph form (I - K H) P (I - K H)^T + K R K^T -> guaranteed PSD even
    // with round-off / a slightly-off K. This SHRINKS P toward the measurement — the one
    // step that mitigates the Slice-14 covariance-pessimism finding when an absolute ref
    // is present (the predict-only no-ref path is unchanged).
    const Mat12 IKH = Mat12::Identity() - K * H;
    state_.cov = symmetrize(IKH * P * IKH.transpose() + K * R * K.transpose());
    state_.twist.cov = state_.cov.block<6, 6>(6, 6);

    // NOTE: stamp untouched — the estimator owns the frontier clock (as in predict()).
    return true;
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
