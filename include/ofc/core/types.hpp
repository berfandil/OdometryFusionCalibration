// ofc/core/types.hpp — fixed-size numeric + geometric types for the core.
#ifndef OFC_CORE_TYPES_HPP
#define OFC_CORE_TYPES_HPP

#include <cstdint>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace ofc {

using Scalar    = double;                 // double everywhere; no fast-math.
using Timestamp = std::int64_t;           // nanoseconds since an arbitrary epoch.
using SourceId  = std::uint8_t;

using Vec3  = Eigen::Matrix<Scalar, 3, 1>;
using Vec6  = Eigen::Matrix<Scalar, 6, 1>;   // twist [v; omega] or se(3) error.
using Mat3  = Eigen::Matrix<Scalar, 3, 3>;
using Mat6  = Eigen::Matrix<Scalar, 6, 6>;
using Mat12 = Eigen::Matrix<Scalar, 12, 12>; // pose(6) + twist(6) covariance.

// Rigid transform on SE(3). Stored as R (orthonormal) + t.
struct SE3 {
    Mat3 R = Mat3::Identity();
    Vec3 t = Vec3::Zero();
};

// A relative motion over a time window, with its covariance.
struct Delta {
    SE3       motion;                 // pose change over [t0, t1]
    Mat6      cov = Mat6::Identity(); // 6x6 in [trans; rot] order
    Timestamp t0  = 0;
    Timestamp t1  = 0;
};

// Body twist (v, omega) with covariance — per-sensor KF / tip extrapolation.
struct Twist {
    Vec6 xi  = Vec6::Zero();          // [v; omega]
    Mat6 cov = Mat6::Identity();
};

// The fused base state.
struct State {
    SE3   pose;                       // in the odom frame (anchored at init)
    Twist twist;
    Mat12 cov = Mat12::Identity();    // pose(SO(3)xR^3) + twist
    Timestamp stamp = 0;
};

} // namespace ofc
#endif // OFC_CORE_TYPES_HPP
