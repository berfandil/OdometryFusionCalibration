// Slice 2 unit tests: weighted split-metric geometric median (Weiszfeld IRLS).
#include <doctest/doctest.h>

#include "ofc/core/lie.hpp"
#include "ofc/core/median.hpp"

#include <cmath>

using namespace ofc;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

median::Params params() {
    median::Params p;
    p.max_iters = 50;
    p.tol       = 1e-8;     // sub-micro tangent step (Weiszfeld near a vertex
                            // converges linearly; this is reached well within 50 iters)
    p.eps       = 1e-9;
    p.lambda    = 1.0;
    return p;
}

SE3 make_se3(const Vec3& w, const Vec3& t) {
    SE3 T;
    T.R = so3::exp(w);
    T.t = t;
    return T;
}
} // namespace

// ---------------------------------------------------------------------------
// Degenerate counts
// ---------------------------------------------------------------------------
TEST_CASE("median of one input is a passthrough") {
    const SE3 x = make_se3(Vec3(0.1, -0.2, 0.3), Vec3(1.0, 2.0, -1.0));
    const Scalar w = 1.0;
    const median::Result r = median::solve(&x, &w, 1, params());
    CHECK(close(r.value.R, x.R, 1e-12));
    CHECK(close(r.value.t, x.t, 1e-12));
    CHECK(r.spread == doctest::Approx(0.0));
    CHECK(r.converged);
}

TEST_CASE("median of two equal-weight inputs is the geodesic midpoint") {
    const SE3 a = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    const SE3 b = make_se3(Vec3(0, 0, kPi / 2), Vec3(2.0, 0.0, 0.0));
    const SE3 xs[2] = {a, b};
    const Scalar ws[2] = {1.0, 1.0};

    const median::Result r = median::solve(xs, ws, 2, params());
    // Equal weights -> midpoint: 45 deg about z, translation (1,0,0).
    const Mat3 expected_R = so3::exp(Vec3(0, 0, kPi / 4));
    CHECK(close(r.value.R, expected_R, 1e-9));
    CHECK(close(r.value.t, Vec3(1.0, 0.0, 0.0), 1e-9));
}

TEST_CASE("median of two inputs honors the weight (biased midpoint)") {
    const SE3 a = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    const SE3 b = make_se3(Vec3(0, 0, 0.0), Vec3(4.0, 0.0, 0.0));
    const SE3 xs[2] = {a, b};
    // Weight b three times a -> the consensus sits at u = 3/4 toward b: x = 3.0.
    const Scalar ws[2] = {1.0, 3.0};
    const median::Result r = median::solve(xs, ws, 2, params());
    CHECK(r.value.t.x() == doctest::Approx(3.0).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// Outlier rejection (the load-bearing property: needs >= 3 inputs)
// ---------------------------------------------------------------------------
TEST_CASE("median of 3+ rejects a single gross outlier (translation)") {
    // Four near-identical inliers clustered around t=(1,0,0); one gross outlier.
    SE3 xs[5];
    xs[0] = make_se3(Vec3(0, 0, 0.00), Vec3(1.00, 0.00, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.01), Vec3(1.02, 0.01, 0.0));
    xs[2] = make_se3(Vec3(0, 0, -0.01), Vec3(0.98, -0.01, 0.0));
    xs[3] = make_se3(Vec3(0, 0, 0.00), Vec3(1.01, 0.00, 0.0));
    xs[4] = make_se3(Vec3(0, 0, 1.50), Vec3(50.0, -40.0, 7.0));   // OUTLIER
    Scalar ws[5] = {1, 1, 1, 1, 1};

    const median::Result r = median::solve(xs, ws, 5, params());
    // The consensus must sit with the inliers, not be dragged to the outlier.
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(0.05));
    CHECK(std::abs(r.value.t.y()) < 0.1);
    CHECK(std::abs(r.value.t.z()) < 0.1);
    CHECK(so3::log(r.value.R).norm() < 0.1);
    CHECK(r.converged);
}

TEST_CASE("median of 3+ rejects a single gross outlier (rotation)") {
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.00, 0.00, 0.30), Vec3(0, 0, 0));
    xs[1] = make_se3(Vec3(0.01, 0.00, 0.31), Vec3(0, 0, 0));
    xs[2] = make_se3(Vec3(0.00, 0.01, 0.29), Vec3(0, 0, 0));
    xs[3] = make_se3(Vec3(2.50, -1.50, -2.0), Vec3(0, 0, 0));   // OUTLIER rotation
    Scalar ws[4] = {1, 1, 1, 1};

    const median::Result r = median::solve(xs, ws, 4, params());
    const Vec3 log_r = so3::log(r.value.R);
    // Consensus rotation near the inlier cluster ~ (0,0,0.3).
    CHECK(std::abs(log_r.z() - 0.30) < 0.05);
    CHECK(std::abs(log_r.x()) < 0.05);
    CHECK(std::abs(log_r.y()) < 0.05);
}

// ---------------------------------------------------------------------------
// Convergence + spread
// ---------------------------------------------------------------------------
TEST_CASE("Weiszfeld converges within the iteration cap") {
    SE3 xs[3];
    xs[0] = make_se3(Vec3(0, 0, 0.1), Vec3(1.0, 0.0, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.2), Vec3(1.1, 0.1, 0.0));
    xs[2] = make_se3(Vec3(0, 0, 0.15), Vec3(0.9, -0.1, 0.0));
    Scalar ws[3] = {1, 1, 1};

    median::Params p = params();
    p.max_iters = 20;
    const median::Result r = median::solve(xs, ws, 3, p);
    CHECK(r.converged);
    CHECK(r.iters <= 20);
}

TEST_CASE("spread is zero for identical inputs and grows with disagreement") {
    SE3 same[3];
    same[0] = make_se3(Vec3(0, 0, 0.2), Vec3(1.0, 0.0, 0.0));
    same[1] = same[0];
    same[2] = same[0];
    Scalar ws[3] = {1, 1, 1};
    const median::Result r_same = median::solve(same, ws, 3, params());
    CHECK(r_same.spread == doctest::Approx(0.0).epsilon(1e-6));

    SE3 spread[3];
    spread[0] = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    spread[1] = make_se3(Vec3(0, 0, 0.0), Vec3(1.0, 0.0, 0.0));
    spread[2] = make_se3(Vec3(0, 0, 0.0), Vec3(-1.0, 0.0, 0.0));
    const median::Result r_spread = median::solve(spread, ws, 3, params());
    CHECK(r_spread.spread > 0.1);
}

TEST_CASE("median is deterministic: identical inputs give identical output") {
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.0, 0.1, 0.3), Vec3(1.0, 0.2, -0.3));
    xs[1] = make_se3(Vec3(0.1, 0.0, 0.2), Vec3(1.1, 0.1, -0.2));
    xs[2] = make_se3(Vec3(-0.1, 0.1, 0.4), Vec3(0.9, 0.3, -0.4));
    xs[3] = make_se3(Vec3(0.0, -0.1, 0.3), Vec3(1.0, 0.2, -0.3));
    Scalar ws[4] = {1.0, 0.7, 1.3, 0.9};

    const median::Result a = median::solve(xs, ws, 4, params());
    const median::Result b = median::solve(xs, ws, 4, params());
    CHECK((a.value.R.array() == b.value.R.array()).all());   // bit-identical
    CHECK((a.value.t.array() == b.value.t.array()).all());
    CHECK(a.spread == b.spread);
}
