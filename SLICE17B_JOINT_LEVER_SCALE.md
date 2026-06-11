# Slice 17b — Turn-regime joint lever + scale (4-unknown hand-eye LS)

**Goal**: recover per-source SCALE without the straight regime — the last calibration DOF unobservable on a drone — and immunize the lever solve against an unrecovered scale (Slice-17 finding: scale 1.08 → lever x off by 3.5 cm via the raw-magnitude `t_B` in the rows).

---

## 1. Math + prototype evidence (2026-06-11, `tools/proto_lever_scale.py`)

Hand-eye translation identity for `A = X∘B_true∘X⁻¹` with the observed `t_B = s_res·t_true` (s_res = the residual scale vs the current fusion `prior_scale` — Phase-2 receives DE-SCALED deltas, estimator.cpp `run_calibration_observe`):

```
(R_A − I)·t_X − κ·(R_X t_B) = −t_A,     κ = 1/s_res
```

LINEAR in `u = [t_X; κ]` — 4 unknowns, 3 equations per turn-gated window. The current production solve is the κ≡1 special case; an unmodeled scale lands as a lever error along the dominant motion axis (the Slice-17 3.5 cm).

**Prototype results** (windowed, Wahba R_X from the rot3d pass; truth in brackets):

| Data | lever err | scale err | 3-unknown control |
|---|---|---|---|
| EuRoC V1_01 drone, scale 1.08 inj | **1.8e-8 mm** | **2.7e-12** (1.08) | 88.3 mm lever err |
| KAIST urban07 planar, scale 1.1 inj | **2.1e-4 mm** (xy; z = exactly the null direction, prior-pinned) | **1.4e-9** (1.1) | singular without ridge |

KAIST AtA eigvals [0, 18.8, 18.8, 7435]: exactly ONE null direction (lever z under yaw-only — the existing known limit) — **the κ column is well-conditioned even on planar ground**. Turn-regime scale works everywhere translation accompanies rotation; the straight-regime Phase-1 scale path remains as an independent estimator.

## 2. Design

**Phase2Calibrator** — the lever LS becomes the joint solve:
- Per-slot normal equations grow 3×3→**4×4** (`ata_`), 3→4 (`atb_`). Row block `J = [(R_A−I) | −(R_X t_B)]`, rhs `−t_A`, accumulated exactly where today's rows accumulate (same turn gate, same `kRotRowMin` guards, same R_X source — the running Kabsch R̂ when the rot3d gate is open, else `R_yp·Rx(roll)`).
- Per-window running **ridge solve** centered on the prior `[t_prior; 1]` (κ prior = 1: residual convention). Per-axis info gating extends to the 4th diagonal: an under-informed κ (e.g. near-zero translation) stays at 1 → scale conf 0 (the spine, per-DOF).
- Observable lever axes vote into the existing `xyz_` histograms (unchanged); the κ axis votes `s_res = 1/κ̂` into a NEW per-slot Phase-2 scale histogram configured from `scale_hist` (same range semantics — a residual ratio clustering at 1).
- The de-scale step inside the row: `t_B` is already prior-de-scaled by the estimator; κ̂ is the residual — matching Phase-1's `scale()` convention exactly.
- Readouts: `scale2(id)` (hist mode), `scale2_confidence(id)`, `scale2_vote_count(id)`; `lever_arm()`/`solve_lever_arm()` keep their shape (the full-solve conditioning guard now reads the 4×4 — `kCondMin` applies to the lever 3×3 sub-conditioning as before, implementer pins the exact guard with tests). `reset_lever` clears the 4-unknown state.

**Estimator** — scale commit becomes two-path:
- New `scale2_committed[slot]` via `commit_gate_reanchor` on `scale2_confidence`/votes (same thresholds).
- RISING EDGE of EITHER scale path folds its residual into `prior_scale` and resets BOTH scale histograms (both are residual-convention — a re-anchor stales both). The existing Phase-1 fold logic is the template; no simultaneous-edge double-fold (process the edges sequentially, second sees the reset histogram).
- `scale_calib`/reference exclusions identical to Phase-1's.

**Config**: `Config::joint_lever_scale` (bool, **default false** — the 3-unknown path stays byte-identical; this changes lever numerics even without scale activity, so it cannot ride an existing knob). Loader `[global] joint_lever_scale`. Config-hash. Recipe for full 3D calibration becomes `rot3d_enabled + joint_lever_scale + subbin_centroid + vote_weight=one`.

**Persistence**: `scale2_committed` persisted → **format v3** (old blobs cold-start; established precedent).

**Strict core**: fixed 4×4/4×1 per-slot arrays; `Eigen::LDLT<Mat4>`-equivalent bounded solve (4×4 fixed-size); no heap.

**REJECTED**:
- Voting Phase-2 scale into Phase-1's `scale_` histogram (shared commit) — couples the calibrators' internals; two independent estimators + a shared fold-target (`prior_scale`) is the looser coupling.
- Replacing Phase-1 scale — straight-regime ratio votes remain the denser/simpler source on ground; the two cross-check.
- Always-on joint solve — silently changes every existing lever number; opt-in keeps the trust apparatus untouched.

## 3. Acceptance

Unit (TDD):
1. Joint recovery, multi-axis conjugated stream (clean, planted lever + scale): lever <1e-3 m AND scale <1e-3 of truth; κ prior-pinned (scale conf 0) when translation is negligible.
2. Planar stream: lever xy + scale recovered; lever z stays at prior (per-axis gate) — extends, never weakens, the existing observability self-tests.
3. Scale feedback: scale2 commit folds into `prior_scale` (two-rate re-anchor), BOTH scale hists reset, no commit thrash; Phase-1 scale path's own fold untouched and still green.
4. `joint_lever_scale=false` (default): byte-identical (exact-equality pin).
5. Slice-17 regression: the 17 lever tests still pass with the flag ON (the joint solve must not degrade the unit-scale cases).
6. Loader key; config-hash flip; persistence v3 round-trip + v2-blob reject.

Gate: `scripts/dev.ps1 -Task test` green; cov-cal NEES untouched (default off).

Real-data (orchestrator):
- EuRoC scale-1.08 injection + full recipe: scale ≈1.08 recovered + COMMITTED from turn regime; lever back to ≪1 cm (target: mm).
- KAIST: turn-path scale agrees with the straight-path scale (~1.10); yaw/offset/lever unchanged.

## 4. Status

- [ ] Implemented (TDD, gate green, committed)
- [ ] Reviewed (`reviews/slice-17b-findings.md`) + findings fixed
- [ ] Real-data validation table filled in
- [ ] Docs updated (CONFIG/DECISIONS/DESIGN/ISSUES) — orchestrator
