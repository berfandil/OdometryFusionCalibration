# Slice 19 — Per-channel split median (D3 amendment: rotation/translation consensus with per-channel weights)

**Goal**: let fusion express per-channel source quality — the fused heading should track the best heading source instead of the middle drifter. Evidence (Slice-18 investigation, `tools/diag_heading_drift.py`): KAIST urban12 fused heading wanders 5–30° (61° spike) while the FOG source holds 4.6° full-drive; all three sources share wheel translation, so the single scalar weight per source cannot prefer FOG's rotation without distorting translation. USER-approved forks (2026-06-12): split the median per channel; weight policy layer (a) config priors only this slice; cross-channel outlier veto opt-in default-ON.

---

## 1. Design

### 1.1 Solver (`median.hpp/.cpp`) — per-channel Weiszfeld

New split solve alongside the existing coupled one (which stays bit-untouched — the default path):

```
SplitResult solve_split(const SE3* deltas, const Scalar* w_rot, const Scalar* w_trans,
                        int n, const Params& p);
```

- **Rotation median**: Weiszfeld IRLS on SO(3) under the geodesic distance `‖log(R_mᵀ R_i)‖`, weights `w_rot`. Same safeguards as the fixed coupled solver (D3): OFF-vertex init (one Karcher mean step about the highest-weight R), Vardi-Zhang coincident-vertex guard, ε-regularized 1/d, `max_iters` cap, n≤2 closed forms (n=2 = weighted geodesic interpolation).
- **Translation median**: Weiszfeld on ℝ³, weighted-mean init, VZ guard, same caps.
- **Per-channel spreads**: weighted RMS geodesic distance (rad) and weighted RMS Euclidean distance (m) — no λ unit-mixing (resolves the D21 scalar-spread wart for the split path).
- `SplitResult = { SE3 value; Scalar spread_rot, spread_trans; int iters_rot, iters_trans; bool converged_rot, converged_trans; }` with `value = {R_rotmedian, t_transmedian}`.
- **Cross-channel outlier veto** (default ON for the split path; `Config::split_veto = true`): after both solves, a source whose normalized channel distance `d_i > kVetoNormDist (3.0) × channel spread` in EITHER channel gets its weight in the OTHER channel scaled by `kVetoWeightScale (0.1)`, and that other channel is re-solved ONCE (bounded WCET: at most one extra IRLS per channel per step). Rationale: hard sensor faults usually corrupt both channels; the veto recovers the coupled solver's whole-source rejection for gross outliers while leaving graceful per-channel weighting for quality differences. Constants fixed (namespace, documented), only the bool is a knob.

### 1.2 Config / weights

- `Config::split_median` (bool, default **false**) — OFF ⇒ the coupled solver runs byte-identically (golden + every existing test untouched). The split solver is mathematically different even at equal weights (coupled IRLS couples the channels), hence opt-in, flip-default-later after real-data soak.
- `SensorConfig::rot_weight_prior` (double, default 1.0, ≥0) — multiplies the source's effective weight in the ROTATION channel only: `w_rot_i = weight_prior_i · reliability_i · rot_weight_prior_i`, `w_trans_i = weight_prior_i · reliability_i`. Policy layer (a): declare the FOG heading-grade in config (datasheet fact). Layers (b) per-channel scatter reliability and (c) GPS-course drift monitor are FOLLOW-UPS, deliberately not in this slice.
- Loader keys: `[global] split_median`, `split_veto`; per-sensor `rot_weight_prior`. Config-hash covers all three. No persistence payload change (no new committed state) — hash only.

### 1.3 Per-channel spread → Q

Split path: `q_trans = q_scale·spread_trans²·I₃ + q_floor[0:3]`, `q_rot = q_scale·spread_rot²·I₃ + q_floor[3:6]` (one `q_scale` knob shared; the sim sweep decides whether a per-channel scale is warranted — add only if the band cannot be met with one). Coupled path untouched.

**Covariance recalibration (the riskiest part)**: the cov-cal NEES band guard [4.0, 5.6] applies to the DEFAULT (coupled) config UNCHANGED. ADD a split-ON cov-cal test: run the existing sweep methodology with `split_median=true`; calibrate the split path's effective q_scale (sweep {0.5, 0.7, 1.0, 1.5}); the split test pins never-overconfident (<5.6) + near-consistency for the calibrated value. If NO swept value lands in band, STOP and report (do not weaken the band, do not ship a knob value outside it).

### 1.4 Estimator wiring

- `split_median` routes the fusion solve to `solve_split`; consensus motion/omega/trans, calibrator feeds, timesync, lifecycle — all consume `SplitResult.value` exactly where `Result.value` flows today.
- Adaptive-Q consumes the two spreads (1.3).
- **Slice-18 Ω_i under split**: the per-channel solvers are independent ⇒ the median influence is block-diagonal `Ω_i = blkdiag(Ω_trans_i, Ω_rot_i)`, each 3×3 from the per-channel Weiszfeld fixed point (same IFT structure per channel with that channel's FINAL weights — including any veto scaling). Extend the FD pin (`tests/test_multi_bias.cpp` harness exists) for the split path; `multi_bias_enabled + split_median` together must be FD-consistent. Coupled path's 6×6 Ω_i untouched.
- Slice-9 reliability stays scalar and multiplies both channels (per-channel reliability = layer (b), out of scope).

## 2. Acceptance

Unit (TDD):
1. Split solver: rotation median rejects a rotation outlier with full translation participation from the same source (the NEW capability — impossible for the coupled solver; pin it as the headline test) and vice versa. n=1/n=2 closed forms; VZ/off-vertex guards mirrored (port the coupled solver's pinning tests incl. the high-weight-outlier guard, per channel).
2. Veto: a gross BOTH-channel faulty source is fully rejected with veto ON (parity with the coupled solver's whole-source rejection); veto OFF leaves the clean channel's contribution; mutation: veto removal fails the parity test.
3. `rot_weight_prior`: a 10× rotation prior on one source pulls the rotation median to it under disagreement while the translation median is unaffected.
4. Default-off: `split_median=false` byte-identical (exact-equality pin + golden untouched).
5. Split-ON cov-cal: calibrated q_scale lands in the band, never-overconfident across the trajectory×noise grid (§1.3; STOP if impossible).
6. Split-ON observability self-tests: every calibration DOF still converges in its regime / not in others (the calibrators consume the split consensus).
7. Ω_i FD pin extended to the split path (block-diagonal), incl. veto-scaled weights.
8. Loader keys; config-hash flips.

Gate: `scripts/dev.ps1 -Task test` fully green.

Real-data (orchestrator, post-merge — KAIST urban07/12/17, recommended urban config + `split_median` + FOG `rot_weight_prior = 10`):
- Fused full-drive heading error **< 10°** (vs the current 5–30° band with a 61° spike; FOG floor is 4.6°).
- urban12 mid-drive transient ≤ the new 129 m floor (expect improvement — truer heading through coasts), tail ≤ 1.4 m, applied fixes ≥ 1494.
- urban07/17 local metrics no regression.

## 3. Status

- [ ] Implemented (TDD, gate green, committed)
- [ ] Reviewed (`reviews/slice-19-findings.md`) + findings fixed
- [ ] Real-data validation table filled in
- [ ] Docs updated (DESIGN §4 / DECISIONS D3-amendment + D29 / CONFIG / ISSUES) — orchestrator
