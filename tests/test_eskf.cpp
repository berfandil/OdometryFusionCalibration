// Slice 2 unit tests: predict-only ESKF integrator (pose compose, twist readout,
// covariance growth + symmetry/PSD, const-velocity tip).
#include <doctest/doctest.h>

#include "ofc/core/absolute_ref.hpp"
#include "ofc/core/eskf.hpp"
#include "ofc/core/lie.hpp"

#include <Eigen/Eigenvalues>

#include <cmath>

using namespace ofc;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

bool is_symmetric(const Mat12& M, Scalar tol = 1e-12) {
    return (M - M.transpose()).cwiseAbs().maxCoeff() < tol;
}

bool is_psd(const Mat12& M, Scalar tol = -1e-9) {
    // Symmetrize for the eigensolver, then check the minimum eigenvalue.
    const Mat12 S = Scalar(0.5) * (M + M.transpose());
    Eigen::SelfAdjointEigenSolver<Mat12> es(S);
    return es.eigenvalues().minCoeff() >= tol;
}

Mat6 q_pose_simple(Scalar v) {
    return v * Mat6::Identity();
}
} // namespace

// ---------------------------------------------------------------------------
// Pose composition + twist readout
// ---------------------------------------------------------------------------
TEST_CASE("predict composes the delta onto the pose") {
    Eskf f;
    f.init(SE3{}, Mat12::Identity());

    SE3 delta;
    delta.R = so3::exp(Vec3(0, 0, kPi / 4));
    delta.t = Vec3(1.0, 0.0, 0.0);

    f.predict(delta, 0.1, q_pose_simple(1e-6));
    CHECK(close(f.state().pose.R, delta.R, 1e-12));
    CHECK(close(f.state().pose.t, delta.t, 1e-12));

    // Compose a second delta -> pose = delta o delta.
    f.predict(delta, 0.1, q_pose_simple(1e-6));
    const SE3 expected = se3::compose(delta, delta);
    CHECK(close(f.state().pose.R, expected.R, 1e-12));
    CHECK(close(f.state().pose.t, expected.t, 1e-12));
}

TEST_CASE("predict sets twist = log(delta)/dt") {
    Eskf f;
    f.init(SE3{}, Mat12::Identity());

    SE3 delta;
    delta.R = so3::exp(Vec3(0, 0, 0.2));
    delta.t = Vec3(0.5, 0.1, 0.0);
    const Scalar dt = 0.1;

    f.predict(delta, dt, q_pose_simple(1e-6));
    const Vec6 expected = se3::log(delta) / dt;
    CHECK(close(f.state().twist.xi, expected, 1e-9));
}

TEST_CASE("predict accumulates a multi-step straight trajectory") {
    Eskf f;
    f.init(SE3{}, Mat12::Identity());

    SE3 step;                          // 1 m forward, no rotation
    step.t = Vec3(1.0, 0.0, 0.0);
    for (int k = 0; k < 5; ++k) f.predict(step, 0.1, q_pose_simple(1e-6));

    CHECK(f.state().pose.t.x() == doctest::Approx(5.0).epsilon(1e-9));
    CHECK(std::abs(f.state().pose.t.y()) < 1e-9);
    CHECK(close(f.state().pose.R, Mat3::Identity(), 1e-12));
}

// ---------------------------------------------------------------------------
// Covariance growth + symmetry + PSD
// ---------------------------------------------------------------------------
TEST_CASE("one predict pins covariance to F P0 F^T + Qmap (analytic)") {
    // Hand-compute the closed-form propagation for a CONCRETE step so a wrong F
    // (e.g. Ad(delta) instead of Ad(delta^-1), or an identity twist block) FAILS.
    // F = blkdiag( Ad(delta^-1), 0 ),  Qmap = blkdiag( q_pose, q_pose/dt^2 ).
    Eskf f;

    // A non-symmetric, non-isotropic P0 so an asymmetric error in F is exposed.
    Mat12 P0 = Mat12::Zero();
    for (int i = 0; i < 12; ++i) P0(i, i) = 0.1 * (i + 1);
    P0(0, 3) = P0(3, 0) = 0.02;     // pose trans<->rot coupling
    P0(6, 9) = P0(9, 6) = 0.03;     // twist v<->omega coupling
    f.init(SE3{}, P0);

    // Delta with BOTH a rotation and an offset translation -> Ad(delta^-1) is neither
    // identity nor equal to Ad(delta) (catches the inverse-adjoint sign of the bug).
    SE3 delta;
    delta.R = so3::exp(Vec3(0.10, -0.20, 0.30));
    delta.t = Vec3(0.7, -0.4, 0.2);
    const Scalar dt = 0.05;

    Mat6 q_pose;
    q_pose.setZero();
    for (int i = 0; i < 6; ++i) q_pose(i, i) = 0.001 * (i + 1);

    f.predict(delta, dt, q_pose);

    // Build the expected propagation independently from the same primitives.
    const Mat6 Ad_inv = se3::adjoint(se3::inverse(delta));
    Mat12 F = Mat12::Zero();
    F.block<6, 6>(0, 0) = Ad_inv;           // twist block stays zero
    Mat12 Q = Mat12::Zero();
    Q.block<6, 6>(0, 0) = q_pose;
    Q.block<6, 6>(6, 6) = q_pose / (dt * dt);
    Mat12 P_expect = F * P0 * F.transpose() + Q;
    P_expect = 0.5 * (P_expect + P_expect.transpose());   // matches eskf symmetrize

    CHECK(close(f.cov(), P_expect, 1e-12));

    // Guard: the wrong (forward-adjoint) F would give a materially different P, so
    // the analytic check above is genuinely discriminating.
    const Mat6 Ad_fwd = se3::adjoint(delta);
    Mat12 F_wrong = Mat12::Zero();
    F_wrong.block<6, 6>(0, 0) = Ad_fwd;
    Mat12 P_wrong = F_wrong * P0 * F_wrong.transpose() + Q;
    P_wrong = 0.5 * (P_wrong + P_wrong.transpose());
    CHECK_FALSE(close(f.cov(), P_wrong, 1e-6));

    // The published twist covariance mirrors the 6x6 twist block.
    CHECK(close(f.state().twist.cov, f.cov().block<6, 6>(6, 6), 1e-12));
    CHECK(is_symmetric(f.cov()));
    CHECK(is_psd(f.cov()));
}

TEST_CASE("covariance grows under predict and stays symmetric PSD") {
    Eskf f;
    f.init(SE3{}, 0.01 * Mat12::Identity());
    const Scalar tr0 = f.cov().trace();

    SE3 delta;
    delta.R = so3::exp(Vec3(0.05, -0.03, 0.2));
    delta.t = Vec3(1.0, 0.2, -0.1);

    for (int k = 0; k < 10; ++k) {
        f.predict(delta, 0.1, q_pose_simple(0.01));
        CHECK(is_symmetric(f.cov()));
        CHECK(is_psd(f.cov()));
    }
    CHECK(f.cov().trace() > tr0);   // process noise + propagation inflate it
}

TEST_CASE("adaptive_q grows with spread and respects the floor") {
    const Scalar floor[6] = {0.001, 0.001, 0.001, 0.002, 0.002, 0.002};

    const Mat6 q_tight = Eskf::adaptive_q(0.0, 1.0, floor);
    // Zero spread -> exactly the floor on the diagonal.
    CHECK(q_tight(0, 0) == doctest::Approx(0.001));
    CHECK(q_tight(3, 3) == doctest::Approx(0.002));

    const Mat6 q_loose = Eskf::adaptive_q(2.0, 1.0, floor);
    // spread^2 * q_scale = 4 added on top of the floor.
    CHECK(q_loose(0, 0) == doctest::Approx(4.0 + 0.001));
    CHECK(q_loose(0, 0) > q_tight(0, 0));

    // No floor pointer -> just the spread term.
    const Mat6 q_nofloor = Eskf::adaptive_q(1.0, 2.0, nullptr);
    CHECK(q_nofloor(0, 0) == doctest::Approx(2.0));
}

TEST_CASE("adaptive_q is exactly q_scale*spread^2 + floor (no /n_eff reduction; D4)") {
    // The former /n_eff "median-variance reduction" (Slice 14 approach A) was REMOVED (D4): it
    // was a fudge for the OLD pinning median's INFLATED spread, and over-divided the now-
    // correctly-sized spread of the TRUE interior median (overconfident). adaptive_q is now the
    // plain q_scale*spread^2 on the diagonal + the additive per-axis floor.
    const Scalar floor[6] = {0.001, 0.001, 0.001, 0.002, 0.002, 0.002};
    const Scalar spread = 2.0, q_scale = 1.0;

    const Mat6 q = Eskf::adaptive_q(spread, q_scale, floor);
    // spread^2 * q_scale = 4 on top of the floor, on EVERY diagonal axis (trans + rot).
    CHECK(q(0, 0) == doctest::Approx(4.0 + 0.001));   // trans axis: spread term + its floor
    CHECK(q(3, 3) == doctest::Approx(4.0 + 0.002));   // rot axis:   same spread term + its floor
    // Off-diagonal is zero (isotropic spread term, diagonal floor).
    CHECK(q(0, 3) == doctest::Approx(0.0));

    // spread == 0 -> the spread term vanishes, so the result is EXACTLY the floor (the noise-free
    // / single-source invariant that keeps the golden committed values unchanged).
    const Mat6 q_zero = Eskf::adaptive_q(0.0, q_scale, floor);
    CHECK(q_zero(0, 0) == doctest::Approx(0.001));
    CHECK(q_zero(3, 3) == doctest::Approx(0.002));

    // q_scale scales the spread term linearly (spread^2 * q_scale).
    const Mat6 q_half = Eskf::adaptive_q(spread, 0.5, floor);
    CHECK((q_half(0, 0) - floor[0]) == doctest::Approx((q(0, 0) - floor[0]) * 0.5));
}

// ---------------------------------------------------------------------------
// Const-velocity tip extrapolation
// ---------------------------------------------------------------------------
TEST_CASE("predict_tip extrapolates const-velocity ahead of the frontier") {
    Eskf f;
    f.init(SE3{}, 0.01 * Mat12::Identity());

    SE3 delta;                          // 1 m forward over 0.1 s -> 10 m/s
    delta.t = Vec3(1.0, 0.0, 0.0);
    f.predict(delta, 0.1, q_pose_simple(0.001));

    State tip;
    f.predict_tip(0.05, 1.5, tip);
    // 10 m/s for 0.05 s -> +0.5 m beyond the frontier x (which is at 1.0).
    CHECK(tip.pose.t.x() == doctest::Approx(1.5).epsilon(1e-6));
    // Tip leads the frontier.
    CHECK(tip.pose.t.x() > f.state().pose.t.x());
    // Inflated covariance: trace strictly larger than the frontier's.
    CHECK(tip.cov.trace() == doctest::Approx(1.5 * f.cov().trace()).epsilon(1e-9));
    CHECK(is_symmetric(tip.cov));
    CHECK(is_psd(tip.cov));
}

TEST_CASE("predict_tip with a turning twist follows the rotation") {
    Eskf f;
    f.init(SE3{}, 0.01 * Mat12::Identity());

    SE3 delta;
    delta.R = so3::exp(Vec3(0, 0, 0.1));   // yaw 0.1 rad over 0.1 s -> 1 rad/s
    delta.t = Vec3(0.5, 0.0, 0.0);
    f.predict(delta, 0.1, q_pose_simple(0.001));

    State tip;
    f.predict_tip(0.1, 1.5, tip);
    // One more window ahead -> pose ~ delta o delta.
    const SE3 expected = se3::compose(delta, delta);
    CHECK(close(tip.pose.R, expected.R, 1e-6));
    CHECK(close(tip.pose.t, expected.t, 1e-6));
}

// ---------------------------------------------------------------------------
// Mahalanobis-gated measurement update (Slice 11)
// ---------------------------------------------------------------------------
namespace {
// Build a POSITION-fix measurement (dim=3) against a target world position z, at the
// CURRENT pose. residual = z - h(x) = z - pose.t; H_pos = [ R | 0_3x3 | 0_3x6 ]
// (3x12, right-error [trans;rot] tangent — the odom-frame translation perturbs as
// t + R*rho to first order under T o exp(eta)); R_noise = sigma^2 I3.
Measurement position_fix(const SE3& pose, const Vec3& z, Scalar sigma) {
    Measurement m;
    m.dim = 3;
    m.residual.setZero();
    m.residual.head<3>() = z - pose.t;
    m.H.setZero();
    m.H.block<3, 3>(0, 0) = pose.R;
    m.R.setZero();
    m.R.block<3, 3>(0, 0) = (sigma * sigma) * Mat3::Identity();
    return m;
}
} // namespace

TEST_CASE("update: empty / out-of-range measurement is a no-op (returns false)") {
    Eskf f;
    f.init(SE3{}, Mat12::Identity());
    const SE3 pose0 = f.state().pose;

    Measurement m;          // dim defaults to 0
    CHECK_FALSE(f.update(m, 9.0));
    CHECK(close(f.state().pose.t, pose0.t, 1e-15));

    m.dim = 7;              // out of range
    CHECK_FALSE(f.update(m, 9.0));
}

TEST_CASE("update: a position fix pulls the pose toward it and shrinks the pose cov") {
    Eskf f;
    f.init(SE3{}, Mat12::Identity());

    // Move the filter to a non-trivial pose first (rotation + translation).
    SE3 delta;
    delta.R = so3::exp(Vec3(0, 0, 0.3));
    delta.t = Vec3(1.0, 0.5, -0.2);
    f.predict(delta, 0.1, q_pose_simple(0.01));

    const Vec3 t_before  = f.state().pose.t;
    const Scalar trP_before = f.cov().block<6, 6>(0, 0).trace();

    // A position measurement 0.2 m off in +x from the current pose (small -> linearization
    // valid). It should pull the pose toward z and reduce the pose-block covariance.
    const Vec3 z = t_before + Vec3(0.2, -0.1, 0.05);
    Measurement m = position_fix(f.state().pose, z, /*sigma=*/0.05);
    CHECK(f.update(m, 1e6));           // huge threshold -> always accepted

    const Vec3 t_after = f.state().pose.t;
    // Moved toward z (closer than before on the residual norm).
    CHECK((t_after - z).norm() < (t_before - z).norm());
    // Pose-block covariance shrank (the one step that reduces P).
    CHECK(f.cov().block<6, 6>(0, 0).trace() < trP_before);
    CHECK(is_symmetric(f.cov()));
    CHECK(is_psd(f.cov()));
    // NIS recorded for an accepted update.
    CHECK(f.last_nis() > 0.0);
}

TEST_CASE("update: a tight, near-zero-residual position fix re-pulls pose almost exactly") {
    Eskf f;
    f.init(SE3{}, 0.5 * Mat12::Identity());
    SE3 delta; delta.t = Vec3(2.0, 0.0, 0.0);
    f.predict(delta, 0.1, q_pose_simple(0.001));

    // A VERY tight measurement (tiny sigma) at a target slightly off -> the posterior pose
    // should land very close to z (gain ~ 1 on the position).
    const Vec3 z = f.state().pose.t + Vec3(0.10, 0.0, 0.0);
    Measurement m = position_fix(f.state().pose, z, /*sigma=*/1e-3);
    CHECK(f.update(m, 1e9));
    CHECK((f.state().pose.t - z).norm() < 1e-2);
}

TEST_CASE("update: a wildly-off measurement is gated out (state unchanged, NIS recorded)") {
    Eskf f;
    f.init(SE3{}, Mat12::Identity());
    SE3 delta; delta.t = Vec3(1.0, 0.0, 0.0);
    f.predict(delta, 0.1, q_pose_simple(0.01));

    const SE3   pose_before = f.state().pose;
    const Mat12 cov_before  = f.cov();
    const Vec6  twist_before = f.state().twist.xi;

    // A gross outlier: 50 m off with a tight sigma -> NIS huge -> rejected by the gate.
    const Vec3 z = f.state().pose.t + Vec3(50.0, 0.0, 0.0);
    Measurement m = position_fix(f.state().pose, z, /*sigma=*/0.1);
    CHECK_FALSE(f.update(m, /*chi2=*/9.0));

    // State + covariance + twist are byte-for-byte untouched.
    CHECK(close(f.state().pose.R, pose_before.R, 1e-15));
    CHECK(close(f.state().pose.t, pose_before.t, 1e-15));
    CHECK(close(f.cov(), cov_before, 1e-15));
    CHECK(close(f.state().twist.xi, twist_before, 1e-15));
    // NIS is still recorded (above threshold) even though the update was rejected.
    CHECK(f.last_nis() > 9.0);
}

// Slice 15: Huber-robust gain down-weighting bounds the heading kick a large position innovation
// injects through the pose trans-rot cross-covariance (the urban12 over-rotation mode). A lateral
// (y) position fix, with a y<->yaw(z) cross-cov in P, rotates the heading via K's rotation row.
TEST_CASE("update: robust kappa attenuates the yaw kick a large position fix injects via cross-cov") {
    // Inject a +y position residual into a filter whose P couples trans-y with yaw-z; return the
    // magnitude of the yaw the update injects. sigma chosen so R is comparable to HPH^T (the
    // realistic-R regime where Huber bites); huge gate so the fix is always accepted (we test the
    // injection magnitude, not the gate). All injected angles stay < pi (so3::log unambiguous).
    auto inj_yaw = [](Scalar resid, Scalar kappa) {
        Eskf f;
        Mat12 P0 = Mat12::Identity();
        P0(1, 5) = P0(5, 1) = Scalar(0.2);             // trans-y <-> yaw-z cross-cov (PSD: 1,0.2;0.2,1)
        f.init(SE3{}, P0);
        const SE3 p0 = f.state().pose;
        const Vec3 z = p0.t + Vec3(Scalar(0), resid, Scalar(0));
        Measurement m = position_fix(p0, z, /*sigma=*/Scalar(0.5));
        const Mat3 R_before = f.state().pose.R;
        const bool ok = f.update(m, /*chi2=*/Scalar(1e12), kappa);   // huge gate -> accepted
        REQUIRE(ok);
        return so3::log(R_before.transpose() * f.state().pose.R).norm();
    };

    // The mechanism: a pure-translation fix injects yaw through the cross-cov.
    CHECK(inj_yaw(Scalar(3.0), Scalar(0)) > Scalar(0));

    // Below threshold (dbar <= kappa): robust is bit-identical to non-robust (w == 1).
    CHECK(inj_yaw(Scalar(0.3), Scalar(0.3)) ==
          doctest::Approx(inj_yaw(Scalar(0.3), Scalar(0))).epsilon(1e-12));

    // Large innovation: robust cuts the injected yaw by more than half vs non-robust.
    const Scalar big0 = inj_yaw(Scalar(9.0), Scalar(0));
    const Scalar bigR = inj_yaw(Scalar(9.0), Scalar(0.3));
    CHECK(bigR < Scalar(0.5) * big0);

    // Huber compresses outliers: the non-robust kick grows ~linearly with the residual, the robust
    // kick grows much slower (sub-linear) -> the robust/non-robust ratio shrinks as residual grows.
    const Scalar ratio0 = inj_yaw(Scalar(9.0), Scalar(0))   / inj_yaw(Scalar(3.0), Scalar(0));
    const Scalar ratioR = inj_yaw(Scalar(9.0), Scalar(0.3)) / inj_yaw(Scalar(3.0), Scalar(0.3));
    CHECK(ratio0 > Scalar(2.5));    // non-robust ~linear (3x residual -> ~3x kick)
    CHECK(ratioR < ratio0);         // robust grows sub-linearly (bounded)
}

// Slice 15b (lever C4): residual-gated rotation-row suppression bounds the heading a large
// position fix injects via the trans-rot cross-cov WHILE KEEPING the translation pull (the property
// Huber lacks). Only fires for rotation-UNOBSERVING fixes above the residual threshold.
TEST_CASE("update: rot-suppress (C4) bounds the yaw kick but keeps the translation pull") {
    // Run a +y position fix into a y<->yaw-coupled P; return (injected yaw magnitude, +y pull).
    auto run_pos = [](Scalar resid, Scalar rot_kappa) {
        Eskf f;
        Mat12 P0 = Mat12::Identity();
        P0(1, 5) = P0(5, 1) = Scalar(0.5);             // trans-y <-> yaw-z cross-cov
        f.init(SE3{}, P0);
        const SE3 p0 = f.state().pose;
        const Vec3 z = p0.t + Vec3(Scalar(0), resid, Scalar(0));
        Measurement m = position_fix(p0, z, /*sigma=*/Scalar(0.5));
        const Mat3 R_before = f.state().pose.R;
        const Vec3 t_before = f.state().pose.t;
        REQUIRE(f.update(m, /*chi2=*/Scalar(1e12), /*robust_kappa=*/Scalar(0), rot_kappa));
        const Scalar yaw  = so3::log(R_before.transpose() * f.state().pose.R).norm();
        const Scalar pull = f.state().pose.t.y() - t_before.y();   // moved toward z (+y)
        return std::make_pair(yaw, pull);
    };

    // Large residual: C4 bounds the heading kick, and (because the un-suppressed over-rotation
    // corrupts even the translation through the SE(3) exp coupling) C4 also CLEANS the position pull.
    const auto off = run_pos(Scalar(8.0), Scalar(0.0));     // (yaw, pull) without C4
    const auto on  = run_pos(Scalar(8.0), Scalar(0.5));     // (yaw, pull) with C4
    CHECK(off.first  > Scalar(0.5));                // a position fix injected a LARGE yaw (the bug)
    CHECK(on.first   < Scalar(0.3) * off.first);    // C4 bounds the heading kick (>70% cut)
    CHECK(on.second  > Scalar(0));                  // C4: translation pulled the RIGHT way (toward z)
    CHECK(on.second  > off.second);                 // ... vs the over-rotation corrupting it without C4

    // Below the residual threshold (dbar <= kappa): C4 is inert -> rotation still corrected.
    CHECK(run_pos(Scalar(0.2), Scalar(0.5)).first ==
          doctest::Approx(run_pos(Scalar(0.2), Scalar(0.0)).first).epsilon(1e-12));

    // A rotation-OBSERVING fix (6-DOF, H observes the rotation error) is NOT suppressed by C4.
    auto run_pose6 = [](Scalar rot_kappa) {
        Eskf f;
        Mat12 P0 = Mat12::Identity();
        P0(1, 5) = P0(5, 1) = Scalar(0.5);
        f.init(SE3{}, P0);
        const Mat3 R_before = f.state().pose.R;
        Measurement m;
        m.dim = 6;
        m.H.setZero();  m.H.block<6, 6>(0, 0) = Mat6::Identity();      // observe pose error directly
        m.residual.setZero();  m.residual(1) = Scalar(8.0);           // large y-translation residual
        m.R.setZero();  m.R.block<6, 6>(0, 0) = Scalar(0.25) * Mat6::Identity();
        REQUIRE(f.update(m, /*chi2=*/Scalar(1e12), /*robust_kappa=*/Scalar(0), rot_kappa));
        return so3::log(R_before.transpose() * f.state().pose.R).norm();
    };
    CHECK(run_pose6(Scalar(0.5)) == doctest::Approx(run_pose6(Scalar(0.0))).epsilon(1e-12));
}

TEST_CASE("update: Joseph form keeps P symmetric PSD across many accepted updates") {
    Eskf f;
    f.init(SE3{}, Mat12::Identity());

    SE3 delta;
    delta.R = so3::exp(Vec3(0.02, -0.01, 0.05));
    delta.t = Vec3(0.5, 0.1, 0.0);

    for (int k = 0; k < 20; ++k) {
        f.predict(delta, 0.1, q_pose_simple(0.005));
        const Vec3 z = f.state().pose.t + Vec3(0.05, -0.02, 0.01);
        Measurement m = position_fix(f.state().pose, z, /*sigma=*/0.05);
        f.update(m, 1e6);
        CHECK(is_symmetric(f.cov()));
        CHECK(is_psd(f.cov()));
    }
}

TEST_CASE("update: the position H = [R|0|0] convention is correct (residual shrinks)") {
    // A discriminating check that H uses pose.R (not identity, not R^T). With the WRONG H
    // the posterior would move the WRONG way / not reduce the residual for a rotated pose.
    Eskf f;
    f.init(SE3{}, Mat12::Identity());
    SE3 delta;
    delta.R = so3::exp(Vec3(0, 0, 1.2));     // a large yaw so R differs from I markedly
    delta.t = Vec3(3.0, -1.0, 0.4);
    f.predict(delta, 0.1, q_pose_simple(0.01));

    const Vec3 z = f.state().pose.t + Vec3(0.3, 0.2, -0.1);
    const Scalar r_before = (z - f.state().pose.t).norm();
    Measurement m = position_fix(f.state().pose, z, /*sigma=*/0.02);
    CHECK(f.update(m, 1e9));
    const Scalar r_after = (z - f.state().pose.t).norm();
    CHECK(r_after < 0.25 * r_before);        // tight fix -> residual collapses
}
