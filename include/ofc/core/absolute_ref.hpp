// ofc/core/absolute_ref.hpp — optional absolute-reference correction plugin (D22).
// Generic ESKF measurement: residual + Jacobian + noise. Mahalanobis-gated by core.
// Loop-closure is intentionally NOT supported (would need retained past states).
#ifndef OFC_CORE_ABSOLUTE_REF_HPP
#define OFC_CORE_ABSOLUTE_REF_HPP

#include "ofc/core/types.hpp"

#include <Eigen/Core>

namespace ofc {

// A linearized measurement of dimension N (compile-time bounded; here a dynamic
// view kept small). residual = z (-) h(x); H = dh/dx; R = measurement noise.
struct Measurement {
    // Up to 6-dim measurements (position=3, pose=6, orientation=3).
    int  dim = 0;
    Eigen::Matrix<Scalar, 6, 1>  residual = Eigen::Matrix<Scalar, 6, 1>::Zero();
    Eigen::Matrix<Scalar, 6, 12> H        = Eigen::Matrix<Scalar, 6, 12>::Zero();
    Eigen::Matrix<Scalar, 6, 6>  R        = Eigen::Matrix<Scalar, 6, 6>::Identity();
    Timestamp stamp = 0;
};

class ICorrection {
public:
    virtual ~ICorrection() = default;

    // Evaluate the measurement model at the current state. Return false if no
    // measurement is available for this step.
    virtual bool evaluate(const State& x, Measurement& out) const = 0;
};

} // namespace ofc
#endif // OFC_CORE_ABSOLUTE_REF_HPP
