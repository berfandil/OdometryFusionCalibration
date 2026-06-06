// ofc/core/lie.cpp — in-house SO(3)/SE(3) operations (Slice 0).
// Conventions: Vec6 tangent is [v; omega] (translation part first).
#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {
namespace {
constexpr Scalar kSmall = 1e-8;   // small-angle threshold
} // namespace

namespace so3 {

Mat3 hat(const Vec3& w) {
    Mat3 W;
    W <<     0, -w.z(),  w.y(),
         w.z(),      0, -w.x(),
        -w.y(),  w.x(),      0;
    return W;
}

Vec3 vee(const Mat3& W) {
    return Vec3(W(2, 1) - W(1, 2),
                W(0, 2) - W(2, 0),
                W(1, 0) - W(0, 1)) * Scalar(0.5);
}

Mat3 exp(const Vec3& w) {
    const Scalar theta2 = w.squaredNorm();
    const Mat3   W = hat(w);
    if (theta2 < kSmall * kSmall) {
        // 2nd-order series; exact enough near zero and avoids division.
        return Mat3::Identity() + W + Scalar(0.5) * W * W;
    }
    const Scalar theta = std::sqrt(theta2);
    const Scalar a = std::sin(theta) / theta;
    const Scalar b = (Scalar(1) - std::cos(theta)) / theta2;
    return Mat3::Identity() + a * W + b * W * W;
}

Vec3 log(const Mat3& R) {
    const Scalar cos_theta = (R.trace() - Scalar(1)) * Scalar(0.5);
    const Scalar c = std::max(Scalar(-1), std::min(Scalar(1), cos_theta));
    const Scalar theta = std::acos(c);
    if (theta < kSmall) {
        return vee(R - R.transpose());   // ~ (R - R^T)/2 unskewed
    }
    const Scalar s = theta / (Scalar(2) * std::sin(theta));
    return s * Vec3(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
}

// Left Jacobian of SO(3), used for SE(3) exp/log.
Mat3 left_jacobian(const Vec3& w) {
    const Scalar theta2 = w.squaredNorm();
    const Mat3   W = hat(w);
    if (theta2 < kSmall * kSmall) {
        return Mat3::Identity() + Scalar(0.5) * W;
    }
    const Scalar theta = std::sqrt(theta2);
    const Scalar a = (Scalar(1) - std::cos(theta)) / theta2;
    const Scalar b = (theta - std::sin(theta)) / (theta2 * theta);
    return Mat3::Identity() + a * W + b * W * W;
}

Mat3 inv_left_jacobian(const Vec3& w) {
    const Scalar theta2 = w.squaredNorm();
    const Mat3   W = hat(w);
    if (theta2 < kSmall * kSmall) {
        return Mat3::Identity() - Scalar(0.5) * W;
    }
    const Scalar theta = std::sqrt(theta2);
    const Scalar half = theta * Scalar(0.5);
    const Scalar cot  = std::cos(half) / std::sin(half);
    const Scalar b = (Scalar(1) / theta2) - cot / (Scalar(2) * theta);
    return Mat3::Identity() - Scalar(0.5) * W + b * W * W;
}

} // namespace so3

namespace se3 {

SE3 compose(const SE3& a, const SE3& b) {
    SE3 out;
    out.R = a.R * b.R;
    out.t = a.R * b.t + a.t;
    return out;
}

SE3 inverse(const SE3& a) {
    SE3 out;
    out.R = a.R.transpose();
    out.t = -out.R * a.t;
    return out;
}

SE3 exp(const Vec6& xi) {
    const Vec3 rho = xi.head<3>();   // translation tangent
    const Vec3 phi = xi.tail<3>();   // rotation tangent
    SE3 out;
    out.R = so3::exp(phi);
    out.t = so3::left_jacobian(phi) * rho;
    return out;
}

Vec6 log(const SE3& T) {
    const Vec3 phi = so3::log(T.R);
    const Vec3 rho = so3::inv_left_jacobian(phi) * T.t;
    Vec6 xi;
    xi.head<3>() = rho;
    xi.tail<3>() = phi;
    return xi;
}

Mat6 adjoint(const SE3& T) {
    // For twist ordering [v; omega]:  Ad = [[R, [t]x R], [0, R]]
    Mat6 A = Mat6::Zero();
    A.block<3, 3>(0, 0) = T.R;
    A.block<3, 3>(0, 3) = so3::hat(T.t) * T.R;
    A.block<3, 3>(3, 3) = T.R;
    return A;
}

SE3 interpolate(const SE3& a, const SE3& b, Scalar u) {
    const Scalar uc = std::max(Scalar(0), std::min(Scalar(1), u));
    SE3 out;
    out.R = a.R * so3::exp(uc * so3::log(a.R.transpose() * b.R));
    out.t = (Scalar(1) - uc) * a.t + uc * b.t;
    return out;
}

Scalar split_distance(const SE3& a, const SE3& b, Scalar lambda) {
    const Vec3 dr = so3::log(a.R.transpose() * b.R);
    const Vec3 dt = a.t - b.t;
    return std::sqrt(dr.squaredNorm() + lambda * dt.squaredNorm());
}

} // namespace se3
} // namespace ofc
