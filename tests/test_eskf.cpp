// Slice 2 unit tests: predict-only ESKF integrator (pose compose, twist readout,
// covariance growth + symmetry/PSD, const-velocity tip).
#include <doctest/doctest.h>

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
