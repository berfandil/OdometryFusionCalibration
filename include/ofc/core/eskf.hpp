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

// Generalized multi-bias augmented error-state capacity (Slice 18, 11b Option B):
// [pose(6); twist(6); bias_1(6) .. bias_k(6)] with k <= Eskf::kMaxBiasSources = 4.
// Fixed worst-case size (strict core: allocated once, no heap); the ACTIVE dimension is
// 12 + 6k, the padding rows/cols are kept exactly zero by construction.
using Mat36 = Eigen::Matrix<Scalar, 36, 36>;

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
    //
    // `robust_kappa` (Slice 15): Huber gain down-weighting. 0 (default) = DISABLED, bit-identical
    // to the non-robust update. When > 0, an outlier innovation (per-DOF RMS Mahalanobis
    // dbar = sqrt(d2/n) > kappa) inflates the active R block by dbar/kappa, recomputing S so the
    // gain and the Joseph covariance shrink consistently — bounding the injected correction. The
    // GATE still tests the TRUE (non-robust) NIS; robustness only attenuates accepted fixes.
    //
    // `rot_suppress_kappa` (Slice 15b, lever C4): 0 (default) = DISABLED. When > 0 AND the
    // measurement does NOT observe rotation (H rotation-error columns ~ 0, e.g. a position fix)
    // AND dbar = sqrt(d2/n) > rot_suppress_kappa, the ROTATION rows of the gain are scaled by
    // rot_suppress_kappa/dbar — bounding the heading the fix injects through the pose trans-rot
    // cross-covariance while leaving the translation correction intact. Orthogonal to robust_kappa.
    bool update(const Measurement& m, Scalar chi2_threshold, Scalar robust_kappa = Scalar(0),
                Scalar rot_suppress_kappa = Scalar(0));

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
    bool update_aug(const Measurement& m, Scalar chi2_threshold, Scalar robust_kappa = Scalar(0),
                    Scalar rot_suppress_kappa = Scalar(0));

    // ---- Generalized multi-bias mode — Slice 18 (11b Option B, median-coupled) ----------
    //
    // The k-source generalization of the Option-A machinery above. The estimator de-biases
    // each biased source's SOURCE-FRAME delta B_i' = B_i o exp(-b_i*dt) BEFORE the frame-align
    // and the median; the consensus med is therefore already de-biased, and the bias error of
    // source i enters the consensus tangent through the MEDIAN'S INPUT SENSITIVITY. The exact
    // first-order coupling block (FD-verified, tests/test_multi_bias.cpp — the spike-5a scalar
    // form -dt*omega_i*Ad(X_i) is FALSIFIED there) is
    //
    //     J_i = blkdiag(R_m^T, I) * Omega_i * blkdiag(R_{A_i}, I) * (-dt * Ad(X_i))
    //     Omega_i = M^-1 u_i P_i,  P_i = I6 - xi_i (W xi_i)^T / d_i^2,  M = sum_j u_j P_j,
    //     u_i = w_i/d_i,  xi_i = [t_i - t_m; log(R_m^T R_i)],  W = diag(lambda*I3, I3),
    //
    // with sum_i Omega_i = I exactly. Regimes (all FD-pinned): n == 1 sole participant ->
    // Omega = I (Option A exact); n == 2 geodesic midpoint -> Omega_i = (w_i/sum w)*I (the
    // solver interpolates with FIXED weights there — NOT the d-based form); a source absent
    // from the window -> u_i = 0 -> Omega_i = 0 (no coupling, its bias random-walks frozen);
    // a (weight-dominant / coincident) vertex d_i <= eps -> the median sits ON the vertex set
    // -> the coincident sources split the identity by weight and every other Omega_j -> 0
    // (the documented non-smooth-vertex limit; no 1/d blow-up).
    //
    // median_influence() computes Omega_i for all n median participants from EXACTLY the
    // inputs/weights the median consumed (the estimator passes its post-reliability fusion
    // weights + the solver's lambda/eps conventions); bias_coupling() assembles J_i. Both are
    // pure static value functions (strict core: fixed-size, no heap) so the FD pin can assert
    // the PRODUCTION blocks against its reference construction directly.
    static void median_influence(const SE3& med, const SE3* A, const Scalar* w, int n,
                                 Scalar lambda, Scalar eps, Mat6* Omega);
    static Mat6 bias_coupling(const SE3& med, const SE3& A_i, const Mat6& Omega_i,
                              const SE3& X_i, Scalar dt);

    // ---- Split-median influence (Slice 19, D3 amendment) --------------------------------
    // The per-channel counterpart of median_influence() for the split solver
    // (median::solve_split): the rotation and translation medians are INDEPENDENT
    // Weiszfeld fixed points, so the median influence is BLOCK-DIAGONAL on the
    // [trans; rot] tangent — Omega_i = blkdiag(Omega_trans_i, Omega_rot_i), each 3x3
    // from that channel's own IFT at its fixed point with that channel's FINAL weights
    // (the veto-scaled effective weights solve_split reports via w_*_final — the
    // influence must describe the median ACTUALLY solved):
    //     channel fixed point:  sum_i u_i x_i = 0,  u_i = w_i / d_i
    //       trans: x_i = t_i - t_m,            d_i = |x_i|        (Euclidean)
    //       rot:   x_i = log(R_m^T R_i),       d_i = |x_i|        (geodesic)
    //     Omega_chan_i = M^-1 u_i P_i,  P_i = I3 - x_i x_i^T / d_i^2,  M = sum_j u_j P_j.
    // Same per-channel regimes as the coupled influence (all FD-pinned in
    // tests/test_multi_bias.cpp): n == 1 -> I3; n == 2 -> the FIXED interpolation weight
    // (w_i / sum w) * I3; absent (w = 0) -> 0; coincident vertex (d <= eps) -> the
    // coincident set splits I3 by weight, others -> 0 (no 1/d blow-up). Channels may take
    // different branches independently. w_of clamp/uniform-fallback mirrors the solver per
    // channel. bias_coupling() above consumes the assembled 6x6 unchanged (its transports
    // are themselves [trans; rot] block-diagonal).
    static void median_influence_split(const SE3& med, const SE3* A, const Scalar* w_rot,
                                       const Scalar* w_trans, int n, Scalar eps,
                                       Mat6* Omega);

    // Maximum number of simultaneously-tracked bias sources (compile-time; augmented error
    // dim <= 12 + 6*4 = 36).
    static constexpr int kMaxBiasSources = 4;

    // Maximum number of median participants median_influence() supports — the shared
    // constant the estimator's source cap must not exceed (review MINOR: was a hard-coded
    // 32 duplicating estimator.cpp's kMaxSourcesCap, with the tail UNINITIALIZED if the
    // caps ever diverged; the estimator static_asserts against this and the tail is now
    // zero-filled defensively).
    static constexpr int kMaxMedianInputs = 32;

    // Switch into the generalized k-bias augmented mode at the current 12-DOF state: all k
    // biases start at zero, each bias block of the 36x36 covariance is seeded to bias_cov0*I6
    // (> 0 required for observability), pose+twist copy the live 12x12, cross-blocks zero,
    // padding zero. Mutually exclusive with the Option-A enable_bias() mode (the estimator
    // routes one or the other; init() resets both). k is clamped to [1, kMaxBiasSources].
    //
    // PER-DOF PINNING (Slice 18 review/B2): when `bias_pn` is non-null it points at k
    // per-source 6-vectors of random-walk rates (the same array predict_multi consumes); a
    // DOF whose rate is NOT > 0 is PINNED at zero — its prior variance is seeded to 0
    // instead of bias_cov0 (and predict_multi adds no walk and zeroes its coupling column),
    // so the DOF is excluded from estimation entirely: bias stays exactly 0 with exactly
    // zero variance. The per-source confidence prior counts only the UNPINNED DOFs (a
    // fully-pinned source reads multi_bias_confidence == 0 forever). nullptr = all DOFs
    // free (the pre-B2 behavior).
    void enable_multi_bias(int k, Scalar bias_cov0, const Vec6* bias_pn = nullptr);
    bool multi_bias_active() const { return mb_count_ > 0; }
    int  multi_bias_count()  const { return mb_count_; }
    // Current bias estimate of bias source j (zero when inactive / j out of range).
    Vec6 multi_bias(int j) const;
    // Per-bias observability confidence in [0,1] (the multi-bias analogue of
    // bias_confidence(): 1 - trace(P_bias_j)/trace(prior)). 0 when inactive / out of range.
    Scalar multi_bias_confidence(int j) const;

    // predict_multi(): the generalized augmented predict. `delta` is the ALREADY-DE-BIASED
    // consensus (the estimator de-biased each source before the median — no further de-bias
    // here, unlike Option A's predict_aug). J_bias[j] (j < multi_bias_count()) is the
    // pose<-bias_j coupling block for THIS step (bias_coupling() above; ZERO for a source
    // absent from the window). bias_pn[j] is the per-source PER-DOF random-walk rate
    // 6-vector (Slice 18 review/B2; a uniform vector reproduces the pre-B2 scalar numerics
    // exactly):
    //   F (36x36 capacity, active 12+6k): pose row [Ad(delta^-1) | 0 | J_1 .. J_k];
    //       twist row 0 (re-read each step); bias rows identity (random walk).
    //   Q = blkdiag(q_pose, q_pose/dt^2, diag(pn_1)*dt, .., diag(pn_k)*dt).
    // A DOF whose rate is NOT > 0 is PINNED (see enable_multi_bias): zero walk AND its
    // column of J_bias[j] is zeroed before entering F, so the pinned bias DOF never couples
    // into the pose and never accrues cross-covariance — it stays exactly 0/0-variance.
    // Mean: pose <- pose o delta; twist <- log(delta)/dt. dt <= 0 guarded as in predict().
    void predict_multi(const SE3& delta, Scalar dt, const Mat6& q_pose,
                       const Mat6* J_bias, const Vec6* bias_pn);

    // update_multi(): the generalized augmented measurement update — identical structure to
    // update_aug() at the 36 capacity (H pads zero bias columns; the gain's bias rows are
    // nonzero via the pose<->bias cross-covariance; per-n gate / Joseph / robust kappa /
    // C4 rot-suppression act on the pose rows 0..5 / 3..5 exactly as before — the pose block
    // is ALWAYS the first 6 rows regardless of k). Injects pose/twist/all k biases.
    // Returns true iff applied; no-op (false) when multi-bias is inactive.
    bool update_multi(const Measurement& m, Scalar chi2_threshold,
                      Scalar robust_kappa = Scalar(0), Scalar rot_suppress_kappa = Scalar(0));

    // NIS (d2) of the LAST update() / update_aug() / update_multi() call — set whether the
    // gate accepted or rejected; 0 before any update. Surfaced into Result diagnostics.
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
    //
    // spread is the weighted-RMS split-distance of the fused inputs to the TRUE interior
    // robust median (median.cpp), i.e. the distance-to-centroid of the participating sources —
    // already the correctly-sized per-window uncertainty of the fused delta. (The former
    // /n_eff "median variance reduction" — Slice 14 approach A — was REMOVED in D4: it was a
    // fudge that compensated for the OLD pinning median's INFLATED spread; with the fixed
    // median it over-divides Q and makes the filter overconfident.)
    static Mat6 adaptive_q(Scalar spread, Scalar q_scale, const Scalar* q_floor);

    // Per-channel adaptive process noise for the SPLIT median path (Slice 19, §1.3):
    //   q_pose[0:3] (trans) = (q_scale * spread_trans^2) * I3 + diag(q_floor[0:3])
    //   q_pose[3:6] (rot)   = (q_scale * spread_rot^2)   * I3 + diag(q_floor[3:6])
    // One SHARED scale knob (the split-ON cov-cal sweep showed no per-channel scale is
    // needed to meet the band) — the estimator passes Config::q_scale_split here, NOT the
    // coupled q_scale (review MAJOR-2: the split path has its own calibration); the two
    // spreads are the split solver's per-channel RMS distances (m / rad — no lambda
    // unit-mixing). [trans; rot] order as adaptive_q().
    static Mat6 adaptive_q_split(Scalar spread_trans, Scalar spread_rot, Scalar q_scale,
                                 const Scalar* q_floor);

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

    // Generalized multi-bias mode — Slice 18 (Option B). Inactive (mb_count_ == 0) by
    // default; enable_multi_bias() turns it on. cov36_'s active top-left (12+6k) block is the
    // live covariance (state_.cov mirrors its top-left 12x12); padding stays exactly zero.
    int    mb_count_ = 0;
    Vec6   mb_bias_[kMaxBiasSources];
    Mat36  cov36_ = Mat36::Zero();
    // Per-source bias-block prior trace at enable (for confidence). Per-source since B2's
    // per-DOF pinning: only the UNPINNED DOFs contribute (0 for a fully-pinned source).
    Scalar mb_prior_trace_[kMaxBiasSources] = {};
};

} // namespace ofc
#endif // OFC_CORE_ESKF_HPP
