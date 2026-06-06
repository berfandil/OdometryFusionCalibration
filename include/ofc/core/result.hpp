// ofc/core/result.hpp — the per-step output struct (D13).
#ifndef OFC_CORE_RESULT_HPP
#define OFC_CORE_RESULT_HPP

#include "ofc/core/types.hpp"

namespace ofc {

enum class Phase { Init, Warmup, Degraded, Nominal };

// Calibration snapshot for one source (relative to the pinned reference).
struct CalibSnapshot {
    SourceId id            = 0;
    SE3      extrinsic;                 // estimated mount
    Scalar   scale         = 1.0;
    Scalar   time_offset_s = 0.0;
    Scalar   confidence    = 0.0;       // peak concentration in [0,1]
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
