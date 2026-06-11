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

    // Diagnostic readout (Slice 17): whether the source's turn-regime FULL rotation
    // extrinsic (rot3d) is currently committed — i.e. prior_extrinsic.R is being driven by
    // the axis-correspondence Wahba solve, superseding the Phase-1 yaw/pitch ∘ roll
    // composition. Always false when Config::rot3d_enabled is off / id unknown. (Kept a
    // facade accessor rather than a CalibSnapshot field — the snapshot's extrinsic already
    // carries the rot3d rotation when committed.)
    bool rot3d_committed(SourceId id) const;

    // Diagnostic readout (Slice 17b): whether the source's turn-regime JOINT scale
    // (the κ axis of the 4-unknown lever+scale hand-eye LS, voted as s_res = 1/κ̂) is
    // currently committed — i.e. its residual has been folded into prior_scale by the
    // rising-edge re-anchor and the latch is holding. Always false when
    // Config::joint_lever_scale is off / id unknown. (Facade accessor like
    // rot3d_committed; the snapshot's scale_committed OR-s this flag in.)
    bool scale2_committed(SourceId id) const;

    // Diagnostic readout (Slice 17b, review MAJOR-2): count of joint-scale votes the
    // vote-site RANGE GUARD withheld (residual s_res outside scale_hist's range, incl.
    // a degenerate κ̂ ≤ 0) — cumulative since init. Distinguishes "out-of-regime scale,
    // permanently uncorrectable by the turn path" (votes 0, skipped HIGH — widen the
    // scale_hist range or fix the prior scale) from "κ unexcited" (votes 0, skipped 0).
    // Always 0 when Config::joint_lever_scale is off / id unknown.
    Scalar scale2_skipped(SourceId id) const;

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
