# Slice-14 q_scale Calibration Sweep Review (Commit 983fa65)

## Findings

include/ofc/core/config.hpp:88: ✓ CORRECT: q_scale default is 0.5 (calibrated), comment accurately documents the safety-first sweep scope and rationale, and q_floor remains unchanged at 1e-6.

tests/test_cov_calibration.cpp: ✓ GENUINE MULTI-TRAJECTORY GUARD, NON-VACUOUS. The test is a permanent, non-trivial regression guard that correctly asserts both (a) NEVER OVERCONFIDENT (nees < 5.5 on all 3 trajectories) and (b) PESSIMISM REDUCED (nees > 1.6, cutting from ~1.0 at q_scale=1 to ~2.07 at q_scale=0.5). Monte-Carlo N adequate: each trajectory runs M=30 seeds with warmup=20 dropped, yielding N>1000 samples. Regression catch verified: (1) q_scale reverting to 1.0 will trip the >1.6 lower bound (nees_old measured ~1.04); (2) q_scale dropping to 0.15 or 0.2 (both noted as overconfident in commit message) will trip the <5.5 upper bound. The test uses distinct trajectories (nees_traj mixed straight+turns, turning worst-case radial, straight axial) that expose the sweep trade-off.

tests/test_cov_calibration.cpp:23: Commit message grid summary (0.2→7.09 OVERCONFIDENT, 0.3→4.81 at 1x, 0.3→5.8 at 2x): the stated worst-case of 0.2 hitting ~7 on `turning` at 1x is correctly identified as breaking the DOF<6 constraint. The choice of 0.5 (worst-case ~2.9 at 1x, ~3.5 at 2x) is the largest pessimism cut staying inside [2,4] and never exceeding 6, confirming the safety-first logic.

tests/test_validation.cpp:362: ✓ SOUND DESIGN DECISION: the NEES/NIS test cases now read `Config{}.q_scale` (0.5) instead of the old pinned q_scale=1. This is correct because: (1) the test documents that it is now exercising the calibrated default, not an arbitrary placeholder; (2) the /n_eff reduction is proven INDEPENDENT via the knob-OFF guard (line 436–447: n_eff=1 at q_scale=0.5 drops NEES to ~0.70, proving the reduction is ~3x and not conflated with the q_scale change). The two effects are cleanly decoupled: old flow was (q_scale=1.0, adaptive_q_source_reduction=ON/OFF)→mean~[0.35 OFF, ~1.04 ON]; new flow is (q_scale=0.5, adaptive_q_source_reduction=ON/OFF)→mean~[0.70 OFF, ~2.07 ON]. The ratio ON/OFF stays ~3x in both eras, proving independence.

tests/test_validation.cpp:423–424 & 444–445: ✓ BANDS CORRECTLY RETUNED. NEES band changed [0.70,1.50]→[1.50,2.70] to bracket the new calibrated-default ensemble-mean ~2.07 (measured value ~2.067 per commit message). Lower bound 1.50 > 0.70 old + requires q_scale NOT to regress to 1.0 (which would drop to ~1.04). Knob-OFF band changed [0.22,0.50]→[0.45,1.00] to bracket the new reduced-OFF value ~0.70 at q_scale=0.5 (was ~0.35 at q_scale=1). NIS band kept [1.5,3.5]; new measured value ~2.835 sits comfortably inside with headroom. All bands absorb Monte-Carlo seed drift yet break on material parameter drift (band ~1.3x tolerance per prose).

tests/test_validation.cpp:402: ✓ CHECK_FALSE(truly_consistent) LOGIC STILL CORRECT. Explanation is sound: NEES ~2.07 < chi2 99% interval ~[5.4,6.6] for consistency, so verdict is intentionally FALSE because the pessimistic margin is deliberate, not a failure to calibrate. This is NOT a weakened assertion; it is honest about the safety-first trade-off.

include/ofc/core/config.hpp:80–87 & 89–95: ✓ COMMENT & GUIDANCE ACCURATE. Comment precisely documents the calibration scope (4 trajectories × 2 noise levels, M=30, worst-case ~2.9/~3.5), the never-overconfident safety constraint, the lower-value breaking point (q_scale=0.2→~7 on turning at 1x), and the model-mismatch rationale. Deployment guidance for q_floor (raise translation floor for zero-spread absolute-ref rigs) is sound and documented.

tests/test_cov_calibration.cpp:232: ✓ WARM-UP CORRECTLY DROPPED. Code line `if (fused_seen <= warmup) continue;` with warmup=20 is applied, so transient frames are excluded before accumulating steady-state pose NEES. Correct practice.

tests/CMakeLists.txt:30: ✓ CMake WIRING CORRECT. test_cov_calibration.cpp is listed in add_executable and properly linked to ofc_core and ofc_sim.

include/ofc/core/config.hpp lines 40–100: ✓ NO THROWAWAY GRID FILES. Commit touches only 4 files; no grid harness or temporary calibration data left in repo. Confirmed by git show 983fa65 --name-only.

tests/test_validation.cpp & test_cov_calibration.cpp: ✓ ASCII-ONLY CONTENT. No non-ASCII or encoding issues detected in the diff. All Lie group / Eigen operations use standard types (Scalar, Vec6, Mat6, SE3).

Trajectory tests in test_cov_calibration.cpp: ✓ DISTINCT & ADEQUATE. Three trajectories (nees_traj with mixed straight+turns, turning worst-case radial slew, straight axial motion) probe different error-source profiles. These are the same rigs used in the offline grid sweep (mentioned in commit message as {nees_traj, mixed, turning, straight}; the permanent test uses a subset for CI speed).

---

## Summary

**Status: CLEAN.** The calibration is methodologically sound: (1) safety-first selection rule (q_scale=0.5 is the largest pessimism cut that never exceeds DOF=6 on any trajectory); (2) the permanent test_cov_calibration guard is non-vacuous and will catch both regress-to-1.0 and over-aggressive-to-0.15 failures; (3) the NEES/NIS bands in test_validation are correctly retuned to pin the new default and remain regression-sensitive (will break on ~1.3x covariance drift); (4) the /n_eff reduction is proven independent from the q_scale change via the knob-OFF guard; (5) no grid artifacts left behind; (6) no unintended config changes. The two highest-risk areas—the new guard's non-vacuity and the q_scale=0.5 vs grid consistency—both pass: grid reports 0.5 worst-case ~2.9 (1x) / ~3.5 (2x), test commits a <5.5 cap (5.5 ≈ 2×2.9 safety margin), and would catch both regressions and overconfident side-effects.
