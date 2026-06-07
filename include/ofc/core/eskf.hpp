// ofc/core/eskf.hpp — predict-only error-state KF "robust integrator" (Slice 2).
//
// State (types.hpp State): pose T in SE(3) + body twist (v,omega) in R^6, with a
// dense 12x12 covariance. Error-state ordering matches the rest of the codebase's
// [trans; rot] convention (DECISIONS D21, D8):
//   indices 0..2  -> pose translation error (R^3)
//           3..5  -> pose rotation error    (so(3))
//           6..8  -> twist linear error     (v)
//           9..11 -> twist angular error    (omega)
//
// This slice has NO correction step (no absolute refs yet) — a KF with only a
// predict is an honest integrator (DECISIONS D4). predict() composes the median
// delta onto the pose and propagates covariance with a right-error model; the twist
// is a direct readout log(delta)/dt (constant-velocity over the window). The tip is
// a const-velocity extrapolation to `now` with an inflated covariance.
//
// STRICT CORE: fixed-size Eigen only (no heap); no exceptions; covariance kept
// symmetric PSD by construction (F P F^T + Q, then symmetrized).
#ifndef OFC_CORE_ESKF_HPP
#define OFC_CORE_ESKF_HPP

#include "ofc/core/absolute_ref.hpp"
#include "ofc/core/types.hpp"

namespace ofc {

class Eskf {
public:
    Eskf() = default;

    // Initialize the filter at pose0 with covariance cov0 (12x12, [pose; twist]).
    // Twist starts at zero; stamp at 0. Resets the filter.
    void init(const SE3& pose0, const Mat12& cov0);

    // Predict-only step over a window of length dt (seconds), consuming the fused
    // median delta:
    //   pose  <- pose o delta                          (right composition)
    //   twist <- log(delta)/dt                          (constant-velocity readout)
    //   P     <- F P F^T + Qmap   then symmetrize       (right-error propagation)
    // where F is block-diagonal: the pose block is Ad(delta^-1) (the right-error of
    // a right composition transports through the inverse adjoint), and the twist
    // block is 0 (the twist is re-read each step, not integrated, so its prior error
    // does not carry over). Qmap places the 6x6 process noise `q_pose` into the pose
    // block and `q_pose / dt^2` into the twist block (delta uncertainty over dt maps
    // to a velocity uncertainty). `q_pose` is in [trans; rot] order.
    // Does NOT touch `stamp`: the caller (estimator) owns the absolute frontier clock
    // and stamps the published state itself. dt <= 0 is treated as a no-op-ish guard
    // (uses a tiny positive dt to avoid division by zero).
    void predict(const SE3& delta, Scalar dt, const Mat6& q_pose);

    // Mahalanobis-gated error-state measurement update (Slice 11, D22). Applies a
    // linearized absolute-reference measurement `m` (residual = z (-) h(x), Jacobian H,
    // noise R) at the current frontier, in the SAME right-error tangent predict() lives
    // in (T_true = T o exp(eta), eta = [trans; rot], pose 0..5, twist 6..11). Steps, all
    // sliced to the measurement's ACTIVE dimension n = m.dim (1..6; top-left n rows of H,
    // first n entries of residual, top-left n x n of R):
    //   S  = H P H^T + R                  (innovation covariance, n x n)
    //   d2 = residual^T S^-1 residual     (NIS — Normalized Innovation Squared)
    //   GATE: d2 > chi2_threshold  =>  REJECT (return false, state UNCHANGED). d2 is still
    //         recorded (last_nis()) for diagnostics, accepted OR rejected.
    //   K  = P H^T S^-1                   (Kalman gain, 12 x n; stable ldlt solve, no inverse)
    //   dx = K * residual                 (12-vector error state)
    //   INJECT: pose  <- pose o se3::exp(dx[0..5])   (full SE(3) exp of the [trans;rot]
    //                                                  tangent — the same coupled-SE(3)
    //                                                  tangent the covariance lives in)
    //           twist <- twist + dx[6..11]            (additive)
    //   P  <- (I - K H) P (I - K H)^T + K R K^T       (Joseph form -> guaranteed PSD), then
    //         symmetrize. twist.cov <- P block(6,6).
    // Returns true iff the update was APPLIED (gate passed AND m.dim in 1..6); false if
    // gated out or empty (m.dim <= 0). Does NOT touch `stamp` (estimator owns the clock,
    // like predict()). This is the ONLY step that SHRINKS P toward a measurement — it
    // directly mitigates the Slice-14 covariance-pessimism finding WHEN an absolute ref is
    // present (the predict-only no-ref path is unchanged).
    //
    // `chi2_threshold` is a single scalar applied REGARDLESS of n. Ideally it would be the
    // chi-square quantile for n DOF (so a 6-DOF pose fix and a 3-DOF position fix gate at
    // different magnitudes); the core currently uses one cfg.mahalanobis_chi2 for all n —
    // documented limitation, not over-engineered here.
    bool update(const Measurement& m, Scalar chi2_threshold);

    // NIS (d2) of the LAST update() call — set whether the gate accepted or rejected; 0
    // before any update(). Surfaced by the estimator into Result diagnostics.
    Scalar last_nis() const { return last_nis_; }

    // Const-velocity extrapolation `dt_ahead` seconds past the frontier:
    //   tip.pose  = pose o exp(twist * dt_ahead)
    //   tip.twist = twist (carried)
    //   tip.cov   = inflation * P            (inflation >= 1 keeps it PSD)
    // Writes into `tip_out`; does not mutate the filter. tip_out.stamp = stamp +
    // dt_ahead (ns).
    void predict_tip(Scalar dt_ahead, Scalar inflation, State& tip_out) const;

    // Build the 6x6 adaptive process noise from the inter-source spread:
    //   q_pose = (q_scale * spread^2) * I6 + diag(q_floor)
    // tight agreement (small spread) -> near the floor; disagreement -> grows
    // quadratically (DESIGN §4, DECISIONS D4). q_floor is a 6-vector in [trans; rot]
    // order; pass nullptr for no floor.
    static Mat6 adaptive_q(Scalar spread, Scalar q_scale, const Scalar* q_floor);

    const State& state() const { return state_; }
    Mat12        cov()   const { return state_.cov; }

private:
    State  state_{};
    Scalar last_nis_ = Scalar(0);   // NIS of the last update() (accepted or rejected)
};

} // namespace ofc
#endif // OFC_CORE_ESKF_HPP
