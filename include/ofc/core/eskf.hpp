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
    // Twist starts at zero; stamp at 0. Resets the filter. Also resets the augmented
    // (bias) mode to OFF — the 12-DOF path is the default (Slice 11b).
    void init(const SE3& pose0, const Mat12& cov0);

    // ------------------------------------------------------------------------------------
    // OPTIONAL AUGMENTED (bias) MODE — Slice 11b, Option A (clean single-driving-source).
    //
    // When a SINGLE source has SensorConfig::bias_states and drives the predict, the filter
    // can carry that source's CONSTANT body-twist bias b = [v_bias; omega_bias] in R^6 as an
    // extra state and let an absolute-ref update OBSERVE + REMOVE it through the pose<->bias
    // cross-covariance (the classic loosely-coupled GPS/INS drift correction, D22). The
    // augmented error-state is 18-DOF [delta_pose(6); delta_twist(6); delta_bias(6)] with a
    // dense Mat18 covariance. Enabling this DOES NOT change the 12-DOF path: predict()/update()
    // (the no-bias methods above/below) stay byte-identical; the augmented behavior lives in
    // predict_aug()/update_aug(). The estimator routes to the aug methods only when a bias
    // source is active, so the DEFAULT (no bias source) is unaffected.
    //
    // enable_bias() switches the filter into augmented mode at the current 12-DOF state:
    // bias starts at zero, the bias block of the 18x18 covariance is seeded to bias_cov0 * I6
    // (a prior on how large the bias might be — must be > 0 for the bias to be observable), the
    // pose+twist blocks copy the current 12x12, and all cross-blocks start at 0. Idempotent-safe
    // (re-enabling re-seeds). Call AFTER init().
    void enable_bias(Scalar bias_cov0);
    bool bias_active() const { return bias_active_; }
    Vec6 bias() const { return bias_; }   // current bias estimate (zero if inactive)

    // Bias OBSERVABILITY CONFIDENCE in [0, 1]: how much the bias-block uncertainty has been
    // REDUCED below its enable_bias() prior, = clamp(1 - trace(P_bias)/trace(P_bias_prior), 0, 1).
    // 0 with NO absolute ref (the bias variance only grows under the random walk -> the bias is
    // never DETERMINED), rising toward 1 as absolute-ref updates shrink it through the pose<->bias
    // cross-covariance. This is the bounded "is the bias observed yet" signal surfaced as
    // CalibSnapshot::bias_observable; the observability self-test reads it staying ~0 with no ref.
    Scalar bias_confidence() const;

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
    // `chi2_threshold` is a RAW (already-per-n) threshold: the gate compares d2 against it
    // directly, with no DOF scaling inside update(). The per-`n` χ² scaling is the CALLER's job
    // (the estimator runs cfg.mahalanobis_chi2 through chi2_gate(base, m.dim) before calling
    // update()) — so every measurement DOF gates at the SAME confidence (Slice 11b). Passing a
    // raw threshold here keeps this method a thin, testable gate (test_eskf.cpp calls
    // update(m, 9.0) etc. directly with raw thresholds).
    bool update(const Measurement& m, Scalar chi2_threshold);

    // Per-DOF Mahalanobis χ² gate threshold (Slice 11b). Scales a base threshold that is
    // tuned at the n=3 position-fix DOF (= cfg.mahalanobis_chi2) by the χ²-quantile ratio so
    // every measurement DOF n in 1..6 gates at the SAME confidence level. n is clamped to 1..6;
    // n==3 returns base unchanged (back-compat). base <= 0 returns base (validate() forbids it).
    //
    // Rationale: the gate's statistical confidence should be DOF-invariant — a 6-DOF pose fix and
    // a 3-DOF position fix should be rejected at the same false-reject rate, NOT at the same raw
    // NIS magnitude. Scaling ONE base by the quantile ratio q[n]/q[3] keeps a single config knob
    // (mahalanobis_chi2, tuned at n=3) while making the effective gate per-`n`. The absolute
    // confidence the quantile table is built at cancels in the ratio to first order; 0.95 is the
    // standard choice. This is a pure value function (no state); the const quantile table lives in
    // eskf.cpp (file-scope constexpr, strict-core: no heap).
    static Scalar chi2_gate(Scalar base_chi2_n3, int n);

    // ---- Augmented (bias) predict/update — Slice 11b, Option A --------------------------
    // The bias-aware counterparts of predict()/update(). The estimator calls THESE (not the
    // 12-DOF versions) for the single driving bias source; the 12-DOF path is untouched.
    //
    // predict_aug(): de-biases the driving delta before propagating, then propagates the 18x18
    // covariance with the bias->pose coupling that makes the bias observable.
    //   de-bias:  Delta_db = delta o se3::exp(-bias * dt)   (bias is a RATE over the window;
    //             integrating -bias over dt gives the se(3) 6-vector -bias*dt)
    //   mean:     pose <- pose o Delta_db ;  twist <- log(Delta_db)/dt   (de-biased readout)
    //   F (18x18, right-error, block rows pose/twist/bias):
    //       pose row : [ Ad(Delta_db^-1) | 0 |  J_pb ]   J_pb = -dt * I6   (THE coupling)
    //       twist row: [ 0 | 0 | 0 ]   (twist re-read each step, as in predict())
    //       bias row : [ 0 | 0 | I6 ]  (constant-bias random walk)
    //   P <- F P F^T + Q,  Q = blkdiag( q_pose, q_pose/dt^2, bias_pn*dt*I6 ), symmetrized.
    // The J_pb = -dt*I6 sign/scale is EXACT for the right-perturbation Delta_db,true =
    //   Delta_db o exp(-delta_bias*dt): the bias error enters the new pose tangent directly as
    //   -dt*delta_bias (verified by the sim recovering the planted bias). Without this block the
    //   bias has ZERO cross-covariance with the pose and a GPS update cannot observe it.
    // `bias_pn` is SensorConfig::bias_process_noise (>= 0). dt <= 0 is guarded as in predict().
    void predict_aug(const SE3& delta, Scalar dt, const Mat6& q_pose, Scalar bias_pn);

    // predict_aug_frozen(): the OUT-OF-REGIME augmented predict. Used when the bias filter is
    // active (bias_active_) but the bias source NO LONGER drives the predict ALONE (another
    // source joined the fusion median, n > 1). In that regime `delta` is the multi-source
    // CONSENSUS, not the bias source's own frame-aligned delta, so the exact -dt*I bias->pose
    // coupling no longer holds. We therefore:
    //   * KEEP applying the learned bias to the de-bias (Delta_db = delta o exp(-bias*dt)) so the
    //     trajectory STILL benefits from the bias the drives-alone phase learned (it is not
    //     silently dropped) — approximate for the consensus but in the right direction;
    //   * propagate the pose/twist blocks of cov18_ EXACTLY as the 12-DOF predict() would
    //     (F_pose = Ad(Delta_db^-1), Q_pose/twist), keeping cov18_'s top-left block CONSISTENT
    //     (no stale/un-propagated 12x12 the way a plain predict() on state_.cov would leave it);
    //   * HOLD the bias mean (no coupling injected) and random-walk the bias block
    //     (Q_bias = bias_pn*dt*I6, so bias_confidence correctly DECAYS — the bias is frozen, not
    //     determined, while out of regime);
    //   * ZERO the now-meaningless pose<->bias cross-covariance, so when the bias source later
    //     drives alone again the next predict_aug RESUMES from a CONSISTENT prior (the coupling
    //     is rebuilt fresh from zero, exactly as enable_bias() seeds it) — never a corrupt one.
    // state_.cov mirrors cov18_'s top-left 12x12 as in predict_aug(). bias_active_ stays true so
    // the caller keeps routing absolute-ref updates through update_aug() (which now sees zero
    // bias cross-cov, so it corrects pose/twist but cannot move the frozen bias). dt <= 0 guarded.
    void predict_aug_frozen(const SE3& delta, Scalar dt, const Mat6& q_pose, Scalar bias_pn);

    // update_aug(): the 18-DOF Mahalanobis-gated measurement update. The measurement H is the
    // SAME [H_pose(n x6) | 0(n x6)] the 12-DOF update uses; the augmented H pads a third zero
    // block (n x6) — GPS observes the pose, NOT the bias directly. But the Kalman gain K = P H^T
    // S^-1 (18 x n) has NONZERO bias rows via the P_pose,bias cross-covariance the predict built,
    // so dx[12:18] != 0 and the update INJECTS the bias: bias += dx[12:18] (pose/twist injected as
    // in update()). Joseph-form 18x18 covariance; gates on the RAW (already-per-n) chi2_threshold
    // exactly as update() does (the estimator pre-scales via chi2_gate); last_nis() updated.
    // Returns true iff applied. No-op (false) when bias is inactive (caller should use update()).
    bool update_aug(const Measurement& m, Scalar chi2_threshold);

    // NIS (d2) of the LAST update() / update_aug() call — set whether the gate accepted or
    // rejected; 0 before any update. Surfaced by the estimator into Result diagnostics.
    Scalar last_nis() const { return last_nis_; }

    // Const-velocity extrapolation `dt_ahead` seconds past the frontier:
    //   tip.pose  = pose o exp(twist * dt_ahead)
    //   tip.twist = twist (carried)
    //   tip.cov   = inflation * P            (inflation >= 1 keeps it PSD)
    // Writes into `tip_out`; does not mutate the filter. tip_out.stamp = stamp +
    // dt_ahead (ns).
    void predict_tip(Scalar dt_ahead, Scalar inflation, State& tip_out) const;

    // Build the 6x6 adaptive process noise from the inter-source spread:
    //   q_pose = (q_scale * spread^2 / max(1, n_eff)) * I6 + diag(q_floor)
    // tight agreement (small spread) -> near the floor; disagreement -> grows
    // quadratically (DESIGN §4, DECISIONS D4). q_floor is a 6-vector in [trans; rot]
    // order; pass nullptr for no floor.
    //
    // n_eff = participating fusion-median source count (Slice 14, approach A). The geometric
    // median of n_eff agreeing sources has ~1/n_eff the per-window variance of a single source,
    // so the spread-derived process noise (spread = inter-source DISAGREEMENT) is divided by
    // n_eff to reflect the FUSED median's accuracy, not the raw disagreement (which over-states
    // the fused error by ~n_eff and inflates the predict-only covariance). Only the SPREAD term
    // is reduced — the floor is an additive per-axis minimum unaffected by the source count.
    // n_eff defaults to 1 (no reduction) for back-compat with direct callers (test_eskf) and the
    // n==1 single-source / bias-driving regime; clamped to >= 1 so it never divides by zero.
    static Mat6 adaptive_q(Scalar spread, Scalar q_scale, const Scalar* q_floor, int n_eff = 1);

    const State& state() const { return state_; }
    Mat12        cov()   const { return state_.cov; }

private:
    State  state_{};
    Scalar last_nis_ = Scalar(0);   // NIS of the last update() (accepted or rejected)

    // Augmented (bias) mode — Slice 11b, Option A. Inactive by default; enable_bias() turns it
    // on. When active, cov18_ is the live 18x18 covariance (state_.cov mirrors its top-left
    // 12x12 for the published frontier) and bias_ is the current body-twist bias estimate.
    bool   bias_active_ = false;
    Vec6   bias_   = Vec6::Zero();
    Mat18  cov18_  = Mat18::Identity();
    Scalar bias_prior_trace_ = Scalar(0);   // trace(P_bias) at enable_bias() (for confidence)
};

} // namespace ofc
#endif // OFC_CORE_ESKF_HPP
