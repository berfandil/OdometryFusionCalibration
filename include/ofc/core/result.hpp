// ofc/core/result.hpp — the per-step output struct (D13).
#ifndef OFC_CORE_RESULT_HPP
#define OFC_CORE_RESULT_HPP

#include "ofc/core/types.hpp"

namespace ofc {

enum class Phase { Init, Warmup, Degraded, Nominal };

// Calibration snapshot for one source (relative to the pinned reference).
//
// Per-DOF confidences are reported SEPARATELY because each calibrated quantity is
// observable in a different motion regime and converges independently (DESIGN §3):
//   * confidence            — the TIME-OFFSET confidence (Slice 5; histogram concentration
//                             of the ‖ω‖-xcorr offset vote). Kept as the primary field for
//                             backward compatibility.
//   * extrinsic_confidence  — the yaw/pitch (so(3) direction) confidence (Slice 6 Phase 1).
//   * scale_confidence      — the per-source scale confidence (Slice 6 Phase 1).
// All are peak concentrations in [0, 1]. `committed` currently reflects the time-offset
// commit gate (Slice 5); per-DOF commit flags arrive with the feedback loop (Slice 8).
struct CalibSnapshot {
    SourceId id            = 0;
    SE3      extrinsic;                 // estimated mount (yaw/pitch in Slice 6; roll Slice 7)
    Scalar   scale         = 1.0;
    Scalar   time_offset_s = 0.0;
    Scalar   confidence    = 0.0;       // TIME-OFFSET peak concentration in [0,1]
    Scalar   extrinsic_confidence = 0.0; // yaw/pitch (so(3)) concentration in [0,1]
    Scalar   scale_confidence     = 0.0; // scale concentration in [0,1]
    bool     committed     = false;
};

// Per-source diagnostics (populated when Config::emit_diagnostics).
struct SourceHealth {
    SourceId id          = 0;
    Scalar   weight      = 0.0;
    Scalar   residual    = 0.0;         // to consensus
    bool     in_window   = false;
    bool     straight    = false;
    bool     turning     = false;
};

struct Result {
    Phase  phase = Phase::Init;

    State  frontier;                    // trustworthy, causal (now - delay)
    State  tip;                         // predict-only extrapolation to now
    bool   tip_valid = false;

    // Calibration + diagnostics (sizes <= max_sources; filled count below).
    CalibSnapshot calib[32];
    SourceHealth  health[32];
    int           source_count = 0;
};

} // namespace ofc
#endif // OFC_CORE_RESULT_HPP
