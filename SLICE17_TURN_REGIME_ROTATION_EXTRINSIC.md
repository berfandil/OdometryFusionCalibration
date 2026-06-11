# Slice 17 — Turn-regime full rotation-extrinsic (axis-correspondence hand-eye)

**Goal**: recover the FULL rotation extrinsic on non-ground (3D, multi-axis-rotation) motion — the gap EuRoC exposed (rotation extrinsic conf 0 on a drone; that error also corrupts the lever hand-eye rows at ~|t|·θ, blocking the <1 cm lever target on 3D motion).

---

## 1. Math + prototype evidence (2026-06-10, `tools/proto_rot_handeye.py`)

A sensor at extrinsic `X` reports `B = X⁻¹∘A∘X`, so rotations satisfy the hand-eye identity `R_A R_X = R_X R_B`, equivalently the **rotation-axis correspondence**

```
a = R_X · b,   a = so3::log(R_A),  b = so3::log(R_B)     (per windowed delta)
```

Accumulate `M = Σ w·a bᵀ` over turn-gated windows; the Wahba/Kabsch solve `R_X = U·diag(1,1,det(UVᵀ))·Vᵀ` (SVD of M) is the LS rotation. Observability lives in `BB = Σ b bᵀ`:

- **rank ≥ 2** (two non-parallel rotation axes) → `R_X` fully determined (Kabsch needs only two correspondences; gate on λ_mid/λ_max, NOT λ_min).
- **rank 1** (ground, yaw-only) → only the 2 axis-tilt DOFs observable; the rotation ABOUT the common axis is blind — and that blind DOF is exactly what Phase-1's straight-regime direction method observes. Complementary by construction (the spine).

**Prototype results** (row-aligned injected-vs-ref increments, composed into windows):

| Data | windows | BB eigvals | recovery |
|---|---|---|---|
| EuRoC V1_01 (drone), clean, win=0.5 s | 204 | 0.30 / 0.42 / 5.3 (rank 3) | **exact** (yaw/pitch/roll 8/5/4° recovered to 4 decimals) |
| EuRoC + synth rot noise 2–5 mrad/step | ~200 | rank 3 | ~1.7° (raw batch Wahba, no robustness layer / no excitation weighting tuning) |
| KAIST urban07 (ground, yaw-only), clean | 253 | **0 / 0 / 19.05 (rank 1)** | err vec = [0, 0, −0.1745]: the injected 10° yaw is EXACTLY the unobservable about-axis DOF; observable tilt DOFs read 0 = truth |

Per-step increments are useless (rotation ≈ 1.5 mrad/step at 200 Hz, gate starves) — **windowed deltas are required**, which is what the calibrator already consumes.

## 2. Design

**Placement**: extend `Phase2Calibrator` (it already owns the turn gate, slots, hand-eye machinery, and the per-window `R_A`/`R_B`). New per-slot state:
- `Mw_[slot]` (3×3 Wahba accumulator), `BBw_[slot]` (3×3 axis Gram), both **decay-aged** by `so3_hist.decay_gamma` per accepted window (independent of the histogram aging mode — a SlidingK ring cannot evict from a rank accumulator; document).
- 3 so(3)-residual histograms `rot3d_[3·slot+c]` configured from `so3_hist` (same shape as Phase-1's channels), basepoint = the slot prior (re-anchorable).

**Per turn-gated window** (inside the existing `observe()` loop, after the `kRotRowMin` guard): `a = so3::log(R_A)`, `b = so3::log(R_B)`; require `‖b‖ ≥ kRotRowMin` too; `Mw += w·a bᵀ`, `BBw += w·b bᵀ` (w = the existing vote weight; |angle| magnitude is baked into the unnormalized a, b — stronger turns weigh more). Then if observable (below), solve the **running** Wahba `R̂` and vote the residual `δφ = so3::log(R̂ · R_bpᵀ)` (basepoint-relative, the contractive parametrization mirroring Phase-1) into the 3 channels.

**Observability gate** (the spine, fixed namespace constant like `kCondMin`): eigendecompose `BBw`; full solve is votable iff `λ_mid ≥ kAxisPairMin · λ_max` (two distinct axes). Rank-1 (ground) stays below the floor forever → channels stay empty → conf 0 → never commits → **planar behavior unchanged**. No partial-subspace voting in this slice (the rank-1 observable DOFs are already covered by Phase-1 + roll; arbitration complexity not worth it — REJECTED alternative below).

**Readout + commit + feedback** (mirror the yaw/pitch DOF):
- `Phase2Calibrator::rot3d(SourceId)` → `exp(φ_mode)·R_bp`-convention extrinsic rotation (implementer pins exact composition with tests: at prior == truth, votes ≈ 0; on clean conjugated data the recovered `R_X` satisfies `a = R_X b` exactly); `rot3d_confidence` = min channel confidence; `rot3d_vote_count`.
- Estimator: new `rot3d_committed[slot]` flag through `commit_gate_reanchor` (same thresholds); **rising edge**: re-anchor basepoint to the recovered rotation + reset the 3 channels (Mw/BBw keep accumulating — they are basepoint-independent).
- **Publish precedence**: when `rot3d_committed`, `prior_extrinsic[i].R ← rot3d(id)` (full 3-DOF — supersedes the `R_yp∘roll` composition); else the existing Phase-1/roll path publishes as today. Rationale: rot3d only ever commits under multi-axis excitation where it strictly dominates (full R vs partial); on ground it never commits and the existing path rules. Phase-1's own commit/re-anchor stays untouched.
- Lever coupling win: once rot3d feeds back, Phase-2's translation rows use a corrected `R_X` → the EuRoC ~3–4 cm lever bias should shrink toward the ridge-only residual.

**Config**: `Config::rot3d_enabled` (bool, **default false** — sim + every existing config byte-identical; turn-regime sim presets must be audited before any default flip). ConfigLoader `[global] rot3d_enabled`. Joins the persistence config-hash. Histograms reuse `so3_hist` (no new HistogramConfig).

**Persistence**: the committed rot3d value + flag must survive warm restart like ext/roll/lever/scale → payload layout changes → **bump `format_version`** (old blobs reject → cold start; consistent with Slice-16's hash-change precedent). Mw/BBw NOT persisted (re-fill, mirroring histogram bins).

**Strict core**: fixed arrays (`3×3·kMaxSources` doubles), Eigen 3×3 `JacobiSVD`/`SelfAdjointEigenSolver` (bounded, already precedented by the lever's eigensolver), no heap/exceptions.

**REJECTED alternatives**:
- Partial-subspace (rank-1) voting on ground — double-estimates pitch against Phase-1 (feedback fight risk) for zero new observability; the planar DOFs are already covered.
- Replacing `best_roll` — the roll path stays as the rank-1 turn-regime estimator; rot3d is additive.
- Batch SVD at readout only (no histograms) — loses the outlier robustness + commit machinery every other DOF has.

**Known limitation**: raw axis-correspondence noise sensitivity (~1.7° at 2–5 mrad/step synthetic) — acceptable because (a) the histogram mode + decay add robustness the prototype lacks, (b) excitation weighting (|a||b| in M) favors high-SNR windows, (c) the EuRoC validation regime (GT-derived odom) is clean; a real-noisy-3D-odometry precision study needs a dataset we don't have wired. Revisit if real 3D sources land.

## 3. Acceptance

Unit (TDD):
1. Conjugated synthetic multi-axis stream (clean): rot3d recovers a planted `R_X` (yaw+pitch+roll) to <1e-3 rad; votes ≈ 0 at prior == truth.
2. **Observability self-test** (load-bearing): yaw-only stream → BBw rank 1 → channels empty, conf 0, never commits; multi-axis stream → converges + commits. Never weaken.
3. Commit + contractive re-anchor: from a large wrong prior, basepoint walks to truth (mirror the Slice-8 ext tests).
4. Precedence: rot3d committed → publishes full R; uncommitted → existing path byte-identical.
5. `rot3d_enabled=false` (default): byte-identical everything (pin by exact equality on a golden-style run).
6. Noise sanity: 2 mrad/step synthetic noise → recovery within a documented loose bound; gate stays honest.
7. Loader key; config-hash flip; persistence round-trip of the committed flag/value + format-version reject of old blobs.

Gate: `scripts/dev.ps1 -Task test` green; cov-cal NEES untouched (default off).

Real-data (orchestrator, post-merge):
- EuRoC `euroc_calib_cen1.ini` + `rot3d_enabled=true`: rotation extrinsic (8/5/4°) recovered + COMMITTED (target err < 0.1° clean-data); **lever re-measured** — expect ≪ 3–4 cm once R_X feeds the translation rows.
- KAIST `calib_vw1_commit_cen.ini` + flag: yaw/scale/offset results UNCHANGED (rot3d inert on planar — never commits).

## 4. Status

- [x] Implemented (TDD, gate green, committed) — `f1eee59`; unit 232 cases
- [x] Reviewed (`reviews/slice-17-findings.md`: APPROVE WITH FIXES, 2 MAJOR + 3 MINOR + 2 NIT) + all fixed — `2d22c28` (+ `f013320` lever-accumulator reset on the rot3d rising edge); unit 239 cases / 14591 asserts
- [x] Real-data validation — table below
- [x] Docs updated (CONFIG/DECISIONS D26/DESIGN/ISSUES) — orchestrator

### Real-data validation (2026-06-11, `rot3d_enabled=true` + centroid + one)

| Probe | result |
|---|---|
| EuRoC rotation extrinsic (truth yaw 8° pitch 5° roll 4°) | log read [0.0635633, 0.0919568, 0.136435] vs truth [0.0635633, 0.0919567, 0.1364349] → err **1.6e-7 rad = 9.1e-6°**, COMMITTED |
| EuRoC lever, unit-scale injection (truth [0.3, −0.2, 0.15]) | **[0.299988, −0.200001, 0.150000]** — µm-level, committed |
| EuRoC lever, scale-1.08 injection | [0.265, −0.201, 0.163] — the 3.5 cm x-residual is ENTIRELY the unrecovered source scale inflating `t_B` in the rows (scale unobservable on a drone: no straight windows) |
| KAIST urban07 planar, flag ON | byte-identical to flag OFF (yaw 0.174711 / scale 1.10008 / toff −0.0500026) — gate honest, inert on ground |

Implementer note kept §2's wording honest: the lever-row coupling is realized calibrator-internally (rows use the running R̂ once the BBw gate opens), and "decay per accepted window" means per rot3d-row-accepted (depositing) window.

**Follow-up identified (next slice candidate)**: the hand-eye translation row is linear in `(t_X, 1/s)` → a turn-regime **joint lever+scale** 4-unknown LS would recover scale WITHOUT the straight regime, closing the last calibration DOF on 3D motion (and the EuRoC lever-under-scale case above).
