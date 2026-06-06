// ofc/core/estimator.hpp — the caller-pumped facade (D14).
// Single-threaded, lock-free, deterministic. The consumer owns the thread and
// pumps step(); each step does fusion + a bounded slice of calibration work.
#ifndef OFC_CORE_ESTIMATOR_HPP
#define OFC_CORE_ESTIMATOR_HPP

#include "ofc/core/absolute_ref.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/source.hpp"
#include "ofc/core/status.hpp"

namespace ofc {

class Estimator {
public:
    Estimator();
    ~Estimator();

    Estimator(const Estimator&) = delete;
    Estimator& operator=(const Estimator&) = delete;

    // Validate config and preallocate everything. No allocation after this.
    Status init(const Config& cfg);

    // Register sources / corrections (pointers owned by caller; must outlive us).
    Status add_source(const ISource* src);
    Status add_correction(const ICorrection* corr);

    // Advance to wall-clock 'now' (nanoseconds). Fusion runs at the causal
    // frontier; a bounded calibration slice runs deeper. Returns NotReady
    // during WARMUP.
    Status step(Timestamp now);

    // Latest result (valid after a successful step).
    const Result& latest() const;

    // Warm-restart persistence (D23): serialize/deserialize into a fixed buffer.
    // File I/O + double-buffering live in an adapter.
    Expected<int> serialize(unsigned char* buf, int cap) const;  // bytes written
    Status        deserialize(const unsigned char* buf, int len);

private:
    struct Impl;        // pimpl over a preallocated arena (no post-init heap)
    Impl* impl_;
};

} // namespace ofc
#endif // OFC_CORE_ESTIMATOR_HPP
