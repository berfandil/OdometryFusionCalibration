// ofc/core/lie.hpp — in-house SO(3)/SE(3) operations (no Sophus/manif dependency).
// Only the handful of ops the framework needs. Implemented in src/core (Slice 0).
#ifndef OFC_CORE_LIE_HPP
#define OFC_CORE_LIE_HPP

#include "ofc/core/types.hpp"

namespace ofc {
namespace so3 {

Mat3 hat(const Vec3& w);        // skew-symmetric generator
Vec3 vee(const Mat3& W);
Mat3 exp(const Vec3& w);        // exponential map  so(3) -> SO(3)
Vec3 log(const Mat3& R);        // logarithm map    SO(3) -> so(3)

} // namespace so3

namespace se3 {

SE3  compose(const SE3& a, const SE3& b);   // a * b
SE3  inverse(const SE3& a);
SE3  exp(const Vec6& xi);                   // [v; omega] -> SE3
Vec6 log(const SE3& T);                     // SE3 -> [v; omega]
Mat6 adjoint(const SE3& T);                 // Ad_T

// Geodesic interpolation between two poses at fraction u in [0, 1]:
//   R = Ra * exp(u * log(Ra^T Rb))   (SO(3) slerp via the rotation log)
//   t = (1 - u) * ta + u * tb        (straight-line lerp in R^3)
// u is clamped to [0, 1]. interpolate(a,b,0)=a, interpolate(a,b,1)=b.
SE3  interpolate(const SE3& a, const SE3& b, Scalar u);

// Split-metric geodesic distance used by the geometric median (D3):
//   d^2 = ||log(Ra^T Rb)||^2 + lambda * ||ta - tb||^2
Scalar split_distance(const SE3& a, const SE3& b, Scalar lambda);

} // namespace se3
} // namespace ofc
#endif // OFC_CORE_LIE_HPP
