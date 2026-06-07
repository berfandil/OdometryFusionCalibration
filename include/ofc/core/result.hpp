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
//
// reliability / bias are the Slice-9 variance-EMA weight refinement diagnostics (D17):
//   * reliability — the current reliability multiplier in [reliability_floor,
//                   reliability_cap], applied to the fusion weight (w = prior ×
//                   reliability × Σ-confidence). 1.0 before warmup / for non-participants.
//                   A NOISY source (large residual scatter) reads < 1; a quiet source > 1.
//   * bias        — the EMA mean of the per-source residual-to-consensus. The systematic
//                   component is left for the calibrator to absorb via the existing
//                   Slice-6/7/8 calibrator observe path; this field is exposed for
//                   DIAGNOSTICS ONLY and is deliberately NOT folded into the fusion weight
//                   (the D17 "bias → calibrator, not weight" contract: a biased-but-
//                   consistent source keeps a high reliability). 0.0 for non-participants.
//                   SEMANTICS: this is the EMA mean of d = split_distance(...), an UNSIGNED
//                   residual MAGNITUDE (>= 0) — "mean residual distance to consensus", NOT a
//                   signed per-DOF offset. It conflates magnitude with direction (a source
//                   biased high or low both read positive; even a zero-mean noisy source
//                   accrues positive bias). The VARIANCE (resid_var), not bias, is what the
//                   weight uses to distinguish noise from a systematic offset.
struct SourceHealth {
    SourceId id          = 0;
    Scalar   weight      = 0.0;
    Scalar   residual    = 0.0;         // to consensus
    Scalar   reliability = 1.0;         // variance-EMA reliability multiplier (Slice 9)
    Scalar   bias        = 0.0;         // EMA-mean residual = systematic component (Slice 9)
    bool     in_window   = false;
    bool     straight    = false;
    bool     turning     = false;
};

// Absolute-reference correction diagnostics for one step (Slice 11, D22). Summarizes the
// registered ICorrection plugins evaluated this step and the Mahalanobis-gated ESKF
// updates they produced. All zero on a step with no corrections registered / none
// available (so corrections-off reads exactly the Slice-2..9 defaults).
//   * corr_evaluated — ICorrection::evaluate() returned a measurement this step (a fix was
//                      available). <= the number of registered corrections.
//   * corr_applied   — gate PASSED and the ESKF update was applied (the fix corrected the
//                      pose).
//   * corr_rejected  — evaluate() yielded a fix but the Mahalanobis gate REJECTED it (NIS >
//                      mahalanobis_chi2; state left unchanged). corr_applied + corr_rejected
//                      == corr_evaluated.
//   * last_nis       — the NIS (Normalized Innovation Squared) of the LAST update() this
//                      step, accepted OR rejected (0.0 if none evaluated). The quantity the
//                      Slice-14 NIS-consistency deferral needs to become computable.
struct CorrectionDiag {
    int    corr_evaluated = 0;
    int    corr_applied   = 0;
    int    corr_rejected  = 0;
    Scalar last_nis       = 0.0;
};

struct Result {
    Phase  phase = Phase::Init;
    // Coarse [0,1] readiness exposed alongside `phase` (DESIGN §7/§8). A monotone
    // map of the phase: Init 0.0 / Warmup 0.25 / Degraded 0.6 / Nominal 1.0.
    Scalar readiness = 0.0;

    State  frontier;                    // trustworthy, causal (now - delay)
    State  tip;                         // predict-only extrapolation to now
    bool   tip_valid = false;

    // Absolute-reference correction summary for this step (Slice 11). Zeroed each step.
    CorrectionDiag correction;

    // Calibration + diagnostics (sizes <= max_sources; filled count below).
    CalibSnapshot calib[32];
    SourceHealth  health[32];
    int           source_count = 0;
};

} // namespace ofc
#endif // OFC_CORE_RESULT_HPP
