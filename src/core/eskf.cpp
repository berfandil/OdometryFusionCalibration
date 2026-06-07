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

Mat18 symmetrize18(const Mat18& M) {
    return Scalar(0.5) * (M + M.transpose());
}

// chi-square quantiles at the 0.95 confidence level, indexed by DOF n = 1..6 (index 0 unused, so
// the table can be indexed directly by n). Used by Eskf::chi2_gate to scale the n=3-tuned base
// threshold to any DOF at the SAME confidence. File-scope constexpr (strict core: no heap).
constexpr Scalar kChi2Q95[7] = {
    Scalar(0),            // [0] unused (n is 1-based)
    Scalar(3.841459),     // n=1
    Scalar(5.991465),     // n=2
    Scalar(7.814728),     // n=3  (the base-tuning DOF -> ratio 1.0 -> identity)
    Scalar(9.487729),     // n=4
    Scalar(11.070498),    // n=5
    Scalar(12.591587),    // n=6
};
} // namespace

void Eskf::init(const SE3& pose0, const Mat12& cov0) {
    state_ = State{};
    state_.pose      = pose0;
    state_.twist.xi  = Vec6::Zero();
    state_.twist.cov = Mat6::Identity();
    state_.cov       = symmetrize(cov0);   // ensure clean symmetry from the start
    state_.stamp     = 0;
    last_nis_        = Scalar(0);          // no update() yet
    // Augmented (bias) mode OFF on init — the 12-DOF path is the default (Slice 11b).
    bias_active_      = false;
    bias_             = Vec6::Zero();
    cov18_            = Mat18::Identity();
    bias_prior_trace_ = Scalar(0);
}

void Eskf::enable_bias(Scalar bias_cov0) {
    // Switch into augmented mode at the CURRENT 12-DOF state. Bias starts at zero; the
    // augmented covariance copies the live pose+twist 12x12 into its top-left block, seeds the
    // bias block to bias_cov0*I6 (the bias prior — must be > 0 to be observable), and zeroes all
    // cross-blocks (the predict builds the pose<->bias coupling from here).
    bias_active_ = true;
    bias_        = Vec6::Zero();
    cov18_       = Mat18::Zero();
    cov18_.block<12, 12>(0, 0) = state_.cov;
    const Scalar bc = (bias_cov0 > Scalar(0)) ? bias_cov0 : Scalar(1e-6);
    cov18_.block<6, 6>(12, 12) = bc * Mat6::Identity();
    cov18_ = symmetrize18(cov18_);
    bias_prior_trace_ = cov18_.block<6, 6>(12, 12).trace();   // for bias_confidence()
}

Scalar Eskf::bias_confidence() const {
    if (!bias_active_ || !(bias_prior_trace_ > Scalar(0))) return Scalar(0);
    const Scalar cur = cov18_.block<6, 6>(12, 12).trace();
    const Scalar reduced = Scalar(1) - cur / bias_prior_trace_;   // 1 - P_bias/P_bias_prior
    return std::max(Scalar(0), std::min(Scalar(1), reduced));
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

Scalar Eskf::chi2_gate(Scalar base_chi2_n3, int n) {
    // base <= 0 is the validate()-forbidden "disabled gate" path: return it untouched so the
    // caller's existing reject-everything behavior is preserved (never NaN/garbage).
    if (!(base_chi2_n3 > Scalar(0))) return base_chi2_n3;

    // Clamp n to the table's supported DOF range [1, 6].
    const int nc = (n < 1) ? 1 : (n > 6 ? 6 : n);

    // Scale the n=3-tuned base by the χ²-quantile ratio q[n]/q[3]. n==3 -> ratio 1.0 -> identity
    // (back-compat). The shared 0.95 confidence cancels in the ratio, so every DOF gates at the
    // SAME confidence with a single config knob.
    return base_chi2_n3 * (kChi2Q95[nc] / kChi2Q95[3]);
}

void Eskf::predict_aug(const SE3& delta, Scalar dt, const Mat6& q_pose, Scalar bias_pn) {
    const Scalar dt_eff = (dt > kMinDt) ? dt : kMinDt;

    // --- de-bias the driving delta -----------------------------------------
    // The bias is a body-twist RATE applied over the window; integrating -bias over dt gives the
    // se(3) 6-vector -bias*dt, so the de-biased delta is Delta_db = delta o exp(-bias*dt).
    const SE3 bias_corr = se3::exp(-bias_ * dt_eff);
    const SE3 delta_db  = se3::compose(delta, bias_corr);

    // --- mean propagation --------------------------------------------------
    state_.pose     = se3::compose(state_.pose, delta_db);   // T <- T o Delta_db
    state_.twist.xi = se3::log(delta_db) / dt_eff;           // de-biased const-velocity readout
    // bias mean is unchanged (random walk; the update is the only thing that moves it).

    // --- covariance propagation: P <- F P F^T + Q (18x18) ------------------
    const Mat6 Ad_inv = se3::adjoint(se3::inverse(delta_db));

    Mat18 F = Mat18::Zero();
    F.block<6, 6>(0, 0)  = Ad_inv;                       // pose <- Ad(Delta_db^-1) * pose
    // pose <- J_pb * bias.  J_pb = -dt*I6: the right-perturbation Delta_db,true =
    // Delta_db o exp(-delta_bias*dt) injects -dt*delta_bias straight into the new pose tangent
    // (exact to first order; the sim pins this sign/scale by recovering the planted bias).
    F.block<6, 6>(0, 12) = -dt_eff * Mat6::Identity();   // THE bias->pose coupling
    // twist row stays zero (re-read each step). bias row = identity (random walk):
    F.block<6, 6>(12, 12) = Mat6::Identity();

    Mat18 Q = Mat18::Zero();
    Q.block<6, 6>(0, 0)   = q_pose;                      // pose increment noise
    Q.block<6, 6>(6, 6)   = q_pose / (dt_eff * dt_eff);  // -> velocity noise
    const Scalar bpn = (bias_pn > Scalar(0)) ? bias_pn : Scalar(0);
    Q.block<6, 6>(12, 12) = (bpn * dt_eff) * Mat6::Identity();   // bias random-walk noise

    cov18_ = symmetrize18(F * cov18_ * F.transpose() + Q);

    // Mirror the published 12x12 (pose+twist) into state_.cov so latest()/frontier read the
    // augmented filter's pose/twist covariance unchanged.
    state_.cov       = cov18_.block<12, 12>(0, 0);
    state_.twist.cov = state_.cov.block<6, 6>(6, 6);
    // stamp untouched — the estimator owns the frontier clock (as in predict()).
}

void Eskf::predict_aug_frozen(const SE3& delta, Scalar dt, const Mat6& q_pose, Scalar bias_pn) {
    const Scalar dt_eff = (dt > kMinDt) ? dt : kMinDt;

    // --- de-bias with the FROZEN learned bias (keep helping the trajectory) -----------
    // delta is the multi-source CONSENSUS here (out of the drives-alone regime), so this is an
    // approximation — but it still removes the bias source's systematic contribution rather than
    // silently dropping the de-bias. The bias MEAN is held (the coupling that would move it is
    // invalid for the consensus, so we do not let an update touch it — see the zeroed cross-cov).
    const SE3 bias_corr = se3::exp(-bias_ * dt_eff);
    const SE3 delta_db  = se3::compose(delta, bias_corr);

    // --- mean propagation (same as predict(), on the de-biased delta) -----------------
    state_.pose     = se3::compose(state_.pose, delta_db);
    state_.twist.xi = se3::log(delta_db) / dt_eff;
    // bias mean unchanged (frozen while out of regime).

    // --- covariance: propagate the 12x12 pose/twist EXACTLY as predict() does, embedded in
    // the 18x18 with NO bias->pose coupling (J_pb = 0: the -dt*I coupling is invalid for the
    // consensus delta), then ZERO the pose<->bias cross-cov so a later predict_aug resume starts
    // consistent. Keeping the block in cov18_ (not on state_.cov via predict()) is the fix: the
    // top-left 12x12 stays propagated across the out-of-regime gap, so resume is never stale.
    const Mat6 Ad_inv = se3::adjoint(se3::inverse(delta_db));

    Mat18 F = Mat18::Zero();
    F.block<6, 6>(0, 0)   = Ad_inv;             // pose <- Ad(Delta_db^-1) * pose (12-DOF F)
    // pose<-bias coupling DELIBERATELY ZERO (invalid for the consensus delta).
    // twist row stays zero (re-read each step). bias row = identity (random walk):
    F.block<6, 6>(12, 12) = Mat6::Identity();

    Mat18 Q = Mat18::Zero();
    Q.block<6, 6>(0, 0)   = q_pose;                      // pose increment noise
    Q.block<6, 6>(6, 6)   = q_pose / (dt_eff * dt_eff);  // -> velocity noise
    const Scalar bpn = (bias_pn > Scalar(0)) ? bias_pn : Scalar(0);
    Q.block<6, 6>(12, 12) = (bpn * dt_eff) * Mat6::Identity();   // bias random-walk noise

    cov18_ = symmetrize18(F * cov18_ * F.transpose() + Q);

    // Zero the (now-meaningless) pose<->bias cross-cov in BOTH off-diagonal blocks so the next
    // drives-alone predict_aug rebuilds the coupling from a clean prior (as enable_bias seeds it).
    cov18_.block<6, 6>(0, 12).setZero();
    cov18_.block<6, 6>(12, 0).setZero();

    // Mirror the published 12x12 (pose+twist), as predict_aug() does.
    state_.cov       = cov18_.block<12, 12>(0, 0);
    state_.twist.cov = state_.cov.block<6, 6>(6, 6);
    // stamp untouched — the estimator owns the frontier clock (as in predict()).
}

bool Eskf::update_aug(const Measurement& m, Scalar chi2_threshold) {
    if (!bias_active_) return false;        // caller should use update() when bias is off
    const int n = m.dim;
    if (n <= 0 || n > 6) return false;      // empty / out-of-range measurement => no-op

    // Pad to full fixed-size capacity (strict core, mirrors update()). H is 6x18: the active
    // measurement on the top-left n x 12 (pose+twist) block; the bias columns (12..17) are
    // ZERO (GPS observes the pose, not the bias directly). Padded rows [n..5] -> 0 with R's
    // unused diagonal at 1 so S stays invertible without coupling into the active n-block.
    Eigen::Matrix<Scalar, 6, 18> H = Eigen::Matrix<Scalar, 6, 18>::Zero();
    Vec6 r = Vec6::Zero();
    Mat6 R = Mat6::Identity();
    H.block(0, 0, n, 12) = m.H.topLeftCorner(n, 12);    // pose+twist part; bias cols stay 0
    r.head(n)            = m.residual.head(n);
    R.topLeftCorner(n, n) = m.R.topLeftCorner(n, n);

    const Mat18& P = cov18_;

    const Mat6 S = H * P * H.transpose() + R;            // innovation covariance (n-block + pad)
    const Eigen::LDLT<Mat6> Sldlt = S.ldlt();
    const Vec6  Sinv_r = Sldlt.solve(r);
    const Scalar d2    = r.dot(Sinv_r);
    last_nis_ = d2;

    if (d2 > chi2_threshold) return false;               // gate: reject (state unchanged)

    // Kalman gain K = P H^T S^-1 (18 x 6 capacity). The bias rows of K are NONZERO via the
    // P_pose,bias cross-covariance the predict built -> the update injects the bias.
    const Eigen::Matrix<Scalar, 6, 18> HP = H * P;       // 6 x 18
    const Eigen::Matrix<Scalar, 6, 18> Kt = Sldlt.solve(HP);
    const Eigen::Matrix<Scalar, 18, 6> K  = Kt.transpose();

    const Eigen::Matrix<Scalar, 18, 1> dx = K * r;       // 18-vector error state

    // Inject (same right-error tangent as update()):
    //   pose  <- pose o exp(dx[0:6]) ;  twist += dx[6:12] ;  bias += dx[12:18]
    state_.pose     = se3::compose(state_.pose, se3::exp(dx.segment<6>(0)));
    state_.twist.xi = state_.twist.xi + dx.segment<6>(6);
    bias_           = bias_ + dx.segment<6>(12);          // <-- the bias removal

    // Joseph-form 18x18 covariance (guaranteed PSD), then symmetrize.
    const Mat18 IKH = Mat18::Identity() - K * H;
    cov18_ = symmetrize18(IKH * P * IKH.transpose() + K * R * K.transpose());

    // Mirror the published 12x12.
    state_.cov       = cov18_.block<12, 12>(0, 0);
    state_.twist.cov = state_.cov.block<6, 6>(6, 6);
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

Mat6 Eskf::adaptive_q(Scalar spread, Scalar q_scale, const Scalar* q_floor, int n_eff) {
    // Divide ONLY the spread term by max(1, n_eff): the geometric median of n_eff agreeing
    // sources has ~1/n_eff the per-window variance of a single source, so the disagreement-
    // derived noise must be scaled down to the FUSED accuracy (Slice 14, approach A). n_eff <= 0
    // and n_eff == 1 both leave the term unchanged (clamp to 1 -> back-compat). The additive
    // floor is the per-axis minimum and is NOT divided. spread == 0 -> the whole term vanishes,
    // so the result is exactly the floor for ANY n_eff (noise-free / single-source invariant).
    const Scalar inv_neff = Scalar(1) / static_cast<Scalar>(n_eff > 1 ? n_eff : 1);
    Mat6 Q = (q_scale * spread * spread * inv_neff) * Mat6::Identity();
    if (q_floor != nullptr) {
        for (int i = 0; i < 6; ++i) Q(i, i) += q_floor[i];
    }
    return Q;
}

} // namespace ofc
