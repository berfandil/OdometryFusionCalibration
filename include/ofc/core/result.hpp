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
//   * extrinsic_confidence  — the EXTRINSIC ROTATION confidence: yaw/pitch (so(3)
//                             direction, Slice 6 Phase 1) AND, once Phase 2 has votes,
//                             the roll (S¹) concentration. Combined as the MIN of the
//                             two (the weakest rotational DOF bounds the joint rotation
//                             reliability); falls back to the Phase-1 value alone before
//                             Phase 2 has any roll votes.
//   * scale_confidence      — the per-source scale confidence (Slice 6 Phase 1).
//   * translation_confidence— the xyz LEVER-ARM confidence (Slice 7 Phase 2): the MIN
//                             over the 3 xyz-channel histogram concentrations. An
//                             unobservable axis (e.g. z under pure-yaw turning) never
//                             concentrates, so this stays low until the lever arm is
//                             fully observable.
// All are peak concentrations in [0, 1]. Per-DOF commit flags (Slice 8 feedback loop)
// mark which DOF has cleared the commit gate (τ_commit ∧ N_min, hysteresis) and is now
// DRIVING fusion (swapped into the source's effective prior) rather than the config prior:
//   * committed              — the TIME-OFFSET commit (Slice 5; kept as the primary flag
//                              for backward compatibility).
//   * extrinsic_committed    — the yaw/pitch (so(3) direction) ROTATION commit. The roll is
//                              a separate refinement that does NOT gate this flag (combining
//                              them would make the flag oscillate during the window where
//                              roll is observed but not yet committed). Once set, the
//                              recovered rotation feeds back into fusion's prior_extrinsic.R.
//   * scale_committed        — the per-source scale commit (feeds fusion's prior_scale).
//   * translation_committed  — the xyz lever-arm commit (feeds prior_extrinsic.t).
struct CalibSnapshot {
    SourceId id            = 0;
    SE3      extrinsic;                 // estimated mount (yaw/pitch Slice 6; + roll/t Slice 7)
    Scalar   scale         = 1.0;
    Scalar   time_offset_s = 0.0;
    Scalar   confidence    = 0.0;       // TIME-OFFSET peak concentration in [0,1]
    Scalar   extrinsic_confidence = 0.0; // rotation (yaw/pitch ∧ roll) concentration in [0,1]
    Scalar   scale_confidence     = 0.0; // scale concentration in [0,1]
    Scalar   translation_confidence = 0.0; // xyz lever-arm concentration in [0,1] (Slice 7)
    bool     committed     = false;     // TIME-OFFSET commit (Slice 5)
    bool     extrinsic_committed     = false; // yaw/pitch ∘ roll commit (Slice 8)
    bool     scale_committed         = false; // scale commit (Slice 8)
    bool     translation_committed   = false; // xyz lever-arm commit (Slice 8)
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
