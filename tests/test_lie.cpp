// Slice 0 unit tests: SO(3)/SE(3) operations vs analytic identities.
#include <doctest/doctest.h>

#include "ofc/core/lie.hpp"

#include <cmath>

using namespace ofc;

namespace {
constexpr Scalar kEps = 1e-9;
constexpr Scalar kPi  = 3.14159265358979323846;

// Single template avoids ambiguity from Eigen expression templates, which can
// construct any fixed-size Matrix at overload-resolution time.
template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}
} // namespace

TEST_CASE("so3 hat/vee are inverse") {
    const Vec3 w(0.3, -0.7, 1.1);
    CHECK(close(so3::vee(so3::hat(w)), w));
}

TEST_CASE("so3 exp produces a rotation matrix") {
    const Vec3 w(0.4, 0.2, -0.9);
    const Mat3 R = so3::exp(w);
    CHECK(close(R * R.transpose(), Mat3::Identity()));
    CHECK(std::abs(R.determinant() - Scalar(1)) < 1e-12);
}

TEST_CASE("so3 exp/log round-trip") {
    for (const Vec3& w : {Vec3(0, 0, 0), Vec3(1e-10, 0, 0),
                          Vec3(0.1, -0.2, 0.3), Vec3(2.0, -1.0, 0.5)}) {
        CHECK(close(so3::log(so3::exp(w)), w, 1e-8));
    }
}

TEST_CASE("so3 exp matches known 90deg about z") {
    const Mat3 R = so3::exp(Vec3(0, 0, kPi / 2));
    Mat3 expected;
    expected << 0, -1, 0,
                1,  0, 0,
                0,  0, 1;
    CHECK(close(R, expected, 1e-12));
}

TEST_CASE("se3 compose/inverse give identity") {
    SE3 T;
    T.R = so3::exp(Vec3(0.2, -0.5, 0.7));
    T.t = Vec3(1.0, -2.0, 3.0);
    const SE3 I = se3::compose(T, se3::inverse(T));
    CHECK(close(I.R, Mat3::Identity()));
    CHECK(close(I.t, Vec3::Zero()));
}

TEST_CASE("se3 exp/log round-trip") {
    Vec6 xi;
    xi << 1.0, -2.0, 0.5, 0.3, -0.4, 0.9;   // [v; omega]
    CHECK(close(se3::log(se3::exp(xi)), xi, 1e-8));

    Vec6 small;
    small << 1e-10, 0, 0, 1e-10, 0, 0;
    CHECK(close(se3::log(se3::exp(small)), small, 1e-8));
}

TEST_CASE("se3 adjoint matches conjugation: Ad_T * xi == log(T exp(xi) T^-1)") {
    SE3 T;
    T.R = so3::exp(Vec3(0.1, 0.2, -0.3));
    T.t = Vec3(0.5, -1.0, 2.0);
    Vec6 xi;
    xi << 0.05, -0.02, 0.03, 0.01, 0.04, -0.06;   // small twist
    const Vec6 lhs = se3::adjoint(T) * xi;
    const SE3  conj = se3::compose(se3::compose(T, se3::exp(xi)), se3::inverse(T));
    const Vec6 rhs = se3::log(conj);
    CHECK(close(lhs, rhs, 1e-6));
}

TEST_CASE("split_distance is zero for equal poses, positive otherwise") {
    SE3 a;
    a.R = so3::exp(Vec3(0.2, 0, 0));
    a.t = Vec3(1, 2, 3);
    CHECK(se3::split_distance(a, a, 1.0) < kEps);

    SE3 b = a;
    b.t.x() += 1.0;
    CHECK(se3::split_distance(a, b, 1.0) == doctest::Approx(1.0));
    CHECK(se3::split_distance(a, b, 4.0) == doctest::Approx(2.0));   // sqrt(4*1)
}
