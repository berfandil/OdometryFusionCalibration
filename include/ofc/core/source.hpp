// ofc/core/source.hpp — the odometry source plugin contract (D7).
// One uniform query; adapters convert native forms (twist / increment / absolute).
#ifndef OFC_CORE_SOURCE_HPP
#define OFC_CORE_SOURCE_HPP

#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

namespace ofc {

// Implemented by adapters (relaxed edge). The core only ever calls query().
// Confidence (Delta.cov) is the source's native covariance if provided, else a
// modeled one; the two are combined per Config::conf_combine.
class ISource {
public:
    virtual ~ISource() = default;

    virtual SourceId id() const = 0;

    // Relative SE(3) motion over [t0, t1], with covariance.
    // Returns Status::NoData if the window is not fully covered by the buffer.
    virtual Expected<Delta> query(Timestamp t0, Timestamp t1) const = 0;
};

} // namespace ofc
#endif // OFC_CORE_SOURCE_HPP
