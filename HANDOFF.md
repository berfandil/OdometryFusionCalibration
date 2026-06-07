# HANDOFF — pick up the autonomous build cold

Single entry point for another agent (or a fresh session) to continue this project without missing anything. Read this top to bottom, then the source-of-truth docs it points to.

---

## 1. What this is

`OdometryFusionCalibration` — a lightweight **C++14 / AUTOSAR** library that fuses N odometry sources into one robust motion estimate (geometric-median-driven error-state integrator) **and** self-calibrates their relative extrinsics, scale, and time-offsets **online**, with calibration closing back into fusion. Organizing principle = the **observability spine** (each calibration DOF is observable in exactly one motion regime).

**Read the docs in this order (they are the source of truth — keep them current):**
1. `DESIGN.md` — architecture spec.
2. `DECISIONS.md` — every decision (chosen / rejected / why), incl. per-slice impl notes.
3. `CONFIG.md` — every config knob.
4. `ISSUES.md` — the slice roadmap + status checkboxes.
5. `WORKFLOW.md` — the operating model + standing decisions + the sub-agent brief template.
6. `reviews/` — the per-slice code-review records (one file per slice).

---

## 2. Build / test gate (the one command)

```
powershell -ExecutionPolicy Bypass -File scripts/dev.ps1 -Task test
```
Configures (if needed) → builds → runs CTest. **Green = build succeeds AND every doctest passes.** This is the bar for every commit.

Toolchain facts (this Windows box):
- No global `cmake`/compiler on PATH. `scripts/dev.ps1` auto-discovers **Visual Studio 2022 Community** via `vswhere`, sources `vcvars64`, and uses VS-bundled `cmake` + `ninja` + `ctest`. MSVC 19.44, C++14.
- Deps fetched by CMake `FetchContent`: **Eigen 3.4.0** + **doctest 2.4.11**. In-house SO(3)/SE(3) (no Sophus).
- A cosmetic `'vswhere.exe' is not recognized` line prints from inside `vcvars64.bat` — **harmless**, build still runs.
- `LF will be replaced by CRLF` git warnings — **harmless** (Windows line endings).
- **Windows PowerShell 5.1 reads `.ps1` as ANSI → keep scripts pure ASCII** (a stray em-dash breaks parsing).

---

## 3. Current state (as of this handoff)

- **Gate green: ctest 2/2 — `unit` 206 cases / 8352 assertions + `adapters` 20 cases / 319 assertions.** ~62 commits on `main`. Feature work through **11b Option A** + the **per-`n` χ² gate** (`e8491dd`) is committed locally; the per-`n` gate commit + this doc-refresh commit may be ahead of the remote (`0c19628`) until pushed. (The gate also builds the relaxed-edge adapters — `dev.ps1` sets `OFC_BUILD_ADAPTERS=ON`, so CTest runs both `unit` and `adapters`.)
- **Done (all green):**

| Unit | What |
|---|---|
| Slice 0 | in-house SO(3)/SE(3) Lie ops, build, doctest harness |
| Slice 1 | `SourceBuffer` ring buffer + uniform `delta(t0,t1)→(SE3,Σ)`, native⊕modeled Σ |
| Slice 2 | geometric-median fusion + predict-only ESKF integrator + estimator wiring (first tracer bullet) |
| Slice 3 | lifecycle `INIT/WARMUP/DEGRADED/NOMINAL` ladder + `readiness` scalar + degrade-don't-block (reference-only dead-reckon, graceful downgrade on source loss); `min_sources_warn` NOMINAL threshold |
| Slice 4 | `Histogram1D` primitive (fixed bins, decay/sliding-K aging, linear-split, sub-bin, circular, concentration confidence) |
| Sim rig | `sim/` ground-truth oracle: trajectory presets + planted `SyntheticSource`s + rig driver |
| Slice 5 | time-sync (‖ω‖ xcorr, pluggable metric, excitation gate, offset histogram, commit N_min + hysteresis) |
| Slice 6 | Phase-1 calibration: straight-gated yaw/pitch (3-ch so(3)@prior) + per-source scale |
| Slice 7 | Phase-2 calibration: turn-gated roll (S¹) + xyz lever-arm (hand-eye LS), both strategies |
| Slice 8 | commit + feedback loop: per-DOF commit (mass + hysteresis), atomic swap, **contractive** re-anchor, cold-start |
| Slice 9 | weight refinement: variance-EMA reliability (bias/variance split, D17) — noisy source downweighted, biased source kept (bias → calibrator); `reliability_floor`/`reliability_cap`, `SourceHealth.reliability`/`bias` |
| Slice 11 | absolute-ref correction path: `Eskf::update` (Mahalanobis-gated, Joseph, right-error full-SE(3) injection) + `add_correction` wired into `step()` + `Result::CorrectionDiag`; sim drift removed 0.58→0.20 m, outlier gated (NIS ~3e5), NIS now computable |
| Slice 10 | per-sensor fixed-lag RTS twist smoother (`TwistSmoother`, D18): CV forward + backward RTS, deeper calibration frontier `now−delay−L`; variance ↓~0.3×, zero-phase, peaks ~4× sharper, no bias; `per_sensor_kf` OFF byte-identical; dropout time-alignment fixed (push-seq stamping) |
| Slice 12 | warm-restart persistence (D23): core `serialize`/`deserialize` into fixed buffers (`persistence.hpp`, "OFCP" v1, FNV-1a config-hash + checksum + orthonormality guard, `CorruptData`/`VersionMismatch`); re-anchor-and-refill restore resumes near-NOMINAL (warm ~0.001 m vs cold ~0.256 m); crash-mid-write/config/checksum/version all reject. Bins NOT persisted; file double-buffer is relaxed-edge (production adapter → 13) |
| Slice 13 (subset) | relaxed-edge `adapters/` (`ofc_adapters`, core PUBLIC API only, no new deps): file-persistence double-buffer (validity-based overwrite target — torn higher-seq can't clobber last-good), threading wrapper (mutex-guarded snapshot), dep-free config loader (subset → Config). ROS node + true fsync → 13b |
| Slice 11b (Option A) | per-source bias states (augmented 18-DOF ESKF): a SINGLE `bias_states` source driving the predict alone → predict de-biases (`Δ∘exp(-b·dt)`) + builds the pose↔bias cross-cov (`J_pb=-dt·I₆`); absolute-ref update removes the bias (sim: planted recovered, drift 16 m → 0.06 m); no-ref observability self-test; multi-bias guard; `predict_aug_frozen` out-of-regime; default-OFF byte-identical. `bias_process_noise` knob; `CalibSnapshot.bias`/`bias_observable`. Option B + per-n gate + GPS adapter → 11b residual |

The **calibration spine (5–8) is complete** — calibration closes back into fusion and bootstraps from arbitrary priors. **All numbered roadmap slices 0–14 are now addressed** (13 + 14 are partial — see below).

- **Remaining work** (any order):
  - Slice 11b residual — **Option A DONE** (single-driving-source augmented bias filter) + **per-`n` χ² Mahalanobis gate DONE** (`e8491dd`: `Eskf::chi2_gate` scales the n=3-tuned base by `q[n]/q[3]`). Deferred: **Option B** (median-coupled multi-source bias — design in the ISSUES Slice-11b DESIGN NOTE) and the concrete **GPS adapter**.
  - Slice 13b — real ROS node + recorded-bag round-trip; replace the persistence adapter's `flush`/`close` with a real `fsync` (durability). Deferred from Slice 13 (no ROS on the dev box).
  - Slice 14 — `[~]` partial: NEES + golden + **NIS** DONE; **init-P covariance fix DONE** (`70c7d38`, NEES ~0.13→~0.35). Remaining: the CONFIG "tuned"-placeholder sweep, and — for strict no-ref NEES consistency — a **distance-aware covariance model** (counter the predict-only translation Ad-inflation).

---

## 4. The autonomous workflow (how to proceed) — from WORKFLOW.md

Work advances **one slice at a time**, sequentially. Per slice, the **orchestrator** (main thread):
1. **Briefs** an implementer sub-agent (use the brief template in `WORKFLOW.md`; point it at the source docs; give exact interfaces + math + done-criteria).
2. Implementer **plans first, works TDD**, verifies the green gate, **commits**.
3. **Review (orchestrator-driven, file hand-off):** launch a reviewer sub-agent that writes `reviews/slice-<n>-findings.md`; then launch a **fix agent** (fresh `general-purpose`) that reads that file, fixes **all** findings, re-verifies green, commits. *(Sub-agents here cannot spawn sub-agents or be resumed — `SendMessage` is unavailable — so the orchestrator runs review and the findings file is the hand-off.)*
4. Orchestrator **independently re-runs the gate**, **updates the source-of-truth docs** (orchestrator owns DESIGN/DECISIONS/CONFIG/WORKFLOW/ISSUES — implementers/fixers must NOT edit them; they report doc-affecting items back), then **summarizes and waits** for the user.

Use a full-capability agent (`general-purpose`/`claude`) for implement + fix; `caveman:cavecrew-reviewer` or `general-purpose` for review (the latter can write the findings file). Surface any user-judgement question rather than guessing.

---

## 5. Conventions & gotchas a new agent MUST NOT violate

- **Strict core / relaxed edges**: `include/ofc/core` + `src/core` are strict (no heap after `init()`, no exceptions, bounded WCET, fixed-capacity, Status-code returns, `double`). `adapters/`, `sim/`, `tests/` are relaxed (std/exceptions/heap fine).
- **Frame-align**: a source-frame delta `B` maps to base as `A = X∘B∘X⁻¹`, `X = SensorConfig::prior_extrinsic` (sensor→base). Fusion also **de-scales first**: `B_corr = {B.R, B.t/prior_scale}`.
- **ESKF**: state = pose `SE(3)` + twist `ℝ⁶`, error `[trans;rot]` (pose 0–5, twist 6–11), dense 12×12. Predict `F = blkdiag(Ad(delta⁻¹), 0)`, `P ← F P Fᵀ + blkdiag(Q, Q/dt²)`. The pose-block `Ad` is the **full SE(3) adjoint** → the covariance is **coupled SE(3)**, not block-diagonal; NEES/consistency must use the full `se3::log` tangent to match (Slice-14). Predict interval = `[last_frontier, frontier]` (gap/overlap-free); `window_s` is bootstrap/lookback only.
- **Time-offset sign**: positive `prior_time_offset_s` ⇒ source clock **ahead** of base (reads `[t0+off, t1+off]`).
- **Phase-1 direction**: 3-channel so(3) histogram @ the per-sensor **prior** basepoint; reverse-fold by the **consensus (fused) sign**; skip a vote ≥90° off prior (avoids the so(3)-log π singularity).
- **Extrinsic recovery is contractive** via the **inverse** minimal rotation: `extrinsic.R = δRᵀ·R_basepoint` (δR = `rotation_between(e_x, g_obs)`). `forward_axis/yaw/pitch` read `δR·e_x` (unchanged). Do NOT revert this to `δR·R_basepoint` — it breaks the Slice-8 bootstrap.
- **Phase-2**: fix yaw/pitch, recover roll (1-D circular) + xyz (`(R_A−I)t_X = R_X t_B − t_A`, 3×3 normal-eq LS + conditioning floor + prior ridge). Pure **yaw-only** turning leaves lever-arm **z unobservable** — needs multi-axis rotation.
- **`vote_weight`** is honored in both phases: `one`/`rotation`(‖ω‖)/`confidence`(Σ)/`combo`(default = rotation×confidence). Under non-`one`, **`commit_min_votes` is a vote-MASS threshold**, not a count, and saturates at `sliding_k` under SlidingK aging.
- **Observability self-tests are load-bearing** — every calibration slice asserts the DOF converges in its regime and does NOT in the others. Never weaken these.
- **Docs are owned by the orchestrator.** Sub-agents report doc changes; the orchestrator edits DESIGN/DECISIONS/CONFIG/WORKFLOW/ISSUES and marks slice checkboxes.

---

## 6. Known limitations / open items

- Extrinsic bootstrap converges from large priors but has a **~0.08 rad realistic floor** on mixed straight+turn trajectories (windows straddling a regime boundary spread the so(3) mode). Sub-0.04 recovery would need revisiting the canonical-rotation-vs-consensus coupling (see Slice-8 fix report).
- `validate()` still has a `TODO: per-sensor + histogram range checks` — nested `HistogramConfig`s are validated at `Histogram1D::configure()` time, not in the top-level `validate()`.
- Several thresholds are tuned placeholders (CONFIG marks them) pending the Slice-14 sweep.
- **Slice-3 lifecycle scope**: NOMINAL is source-count-driven (`n ≥ min_sources_warn`), not directly calibration-convergence-gated; under `ReferenceOnly` cold-start the DEGRADED→NOMINAL upgrade tracks convergence only *indirectly* (a source joins the median once its extrinsic commits). `min_sources_warn` is validated lower-bound only (`≥1`); a value above `max_sources` is legitimate (NOMINAL never reached). If a future slice wants readiness to encode calibration convergence directly, revisit the ladder.
- **Slice-9 weight scope**: reliability is the variance-EMA quality factor; `sigma_confidence()`'s D21 unit-mixing (mean of trans m² + rot rad² in one scalar) is **left intact** — reliability was added multiplicatively, not as a unit-separation rewrite, so that caveat stays open. `SourceHealth.bias` is an unsigned residual *magnitude* (mean `split_distance`), not a signed per-DOF offset — it cannot itself distinguish direction; the weight uses `resid_var` (scatter), not `bias`.
- **Covariance pessimism — init-P RESOLVED; residual is the predict-only Ad-inflation**: the Monte-Carlo NEES (6-DOF pose, full-SE(3) `se3::log` tangent) was ensemble-mean ≈ **0.13** vs DOF **6** (~46×) because the ESKF initialized `P = I₁₂`. **Fixed** (`70c7d38`): init `P = 0` at the gauge-anchored first fuse and let the first `predict()` establish the one-window covariance `blkdiag(q_pose, q_pose/dt²)` (no seed, no double-count) → NEES now ≈ **0.35** (~17×). The REMAINING pessimism is the predict-only right-error `Ad(delta⁻¹)` **translation-block inflation** as the pose accumulates motion (rotation block is well-calibrated, ~0.14 rad²) — NOT `q_scale`/`q_floor`-tunable. Strict no-ref consistency would need a **distance-aware covariance model** and/or a no-ref correction. Slice 11 mitigates it WHEN an absolute ref is present (the correction shrinks P → NIS ≈ 2.4 vs DOF 3). The NEES test pins the current value (band [0.22,0.50]) + `CHECK_FALSE(truly_consistent)` so a regression (e.g. reverting to `I₁₂`) re-trips.
- **Slice-11 correction-gate — per-`n` now (RESOLVED, `e8491dd`)**: the Mahalanobis gate WAS a single scalar `mahalanobis_chi2` regardless of measurement DOF `n` (~97% quantile for n=3, ~80% for n=6). Now `Eskf::chi2_gate(base, n)` scales the n=3-tuned base by the χ²-quantile ratio `q[n]/q[3]` (const `kChi2Q95`, 0.95 confidence) so every DOF `n∈1..6` gates at the same confidence; the estimator passes `chi2_gate(cfg.mahalanobis_chi2, m.dim)` to `update()`/`update_aug()` (signatures unchanged — they still take a raw per-n threshold). n=3 returns base unchanged, so the dim=3 position fixes shipping now are behaviorally IDENTICAL; the per-`n` scaling becomes load-bearing once a 6-DOF/mixed plugin lands. The drift-removal test still uses a loose `chi2=100` — a **test artifact** driven by the covariance pessimism above (legit drift residuals exceed `chi2=9` because P never shrinks on predict-only stretches), NOT a production value (and unrelated to the per-`n` gate).
- **Covariance tangent (doc-vs-code reconciled)**: the ESKF covariance is **full coupled SE(3)** (propagated by the full `Ad(delta⁻¹)`), NOT block-diagonal "decoupled SO(3)×ℝ³" — that phrase only describes the error ordering `[trans;rot]` + the median's split metric. NEES uses the full `se3::log` to match. (DESIGN §5, D21 now corrected.)
- **Slice-10 smoother scope**: a single shared `TwistSmoother` uses the MAX `kf_process_noise` over enabled sources (not per-source `q`/`r`; `r` fixed 1.0). The refined RTS covariance is computed + exposed but NOT wired into the calibrator vote weight (the deeper path still feeds the raw Σ-confidence) — D18's "refined Σ" is half-wired. Fixed ~7.5 MB footprint (compile-time caps 32×65, paid even when off). The lag `L` is a step count via the nominal tick rate → off-cadence pumping gives an effective time-lag ≠ `calib_lag_s` (rings stay step-aligned regardless). Per-source dropout degrades to a variance loss (no bias) after the push-seq alignment fix.
- **Slice-12 persistence scope**: histogram BINS are not persisted — restore re-anchors the calibrators to the committed values + holds them via hysteresis (an `offset_restored` latch does the same for the time-offset DOF), so a restored DOF's `committed` flag + value are correct immediately but its confidence reads low until the histogram re-fills. `config_hash` covers rig shape + per-sensor priors/flags + histogram configs + commit/gate thresholds + `phase2_strategy`/`calib_lag_s` (runtime knobs excluded).
- **Slice-13 adapters scope**: relaxed-edge only (std/heap/exceptions/threads/file IO). The ROS node is a compile-guarded header sketch (not built — no ROS on the dev box) → **13b**. The file-persistence adapter uses `flush`/`close`, NOT `fsync` — an OS-page-cache power-loss window remains even with the crash-safe target selection (durable `fsync` → **13b**). The config loader covers a documented knob SUBSET (not every `Config` field). Adapters are off by default (`OFC_BUILD_ADAPTERS=OFF`); the gate turns them on.
- **Slice-11b bias scope**: **Option A done** — augmented 18-DOF filter for a SINGLE source that drives the predict ALONE (1-source / `ReferenceOnly` dead-reckon). Out-of-regime (a 2nd source joins the median) → `predict_aug_frozen` holds the learned bias + keeps `cov18_` consistent, but de-biases the *consensus* (an approximation; exact per-contribution de-bias is Option-B territory). **Per-`n` χ² gate also done** (`e8491dd`). Deferred: Option B (median-coupled multi-source), the GPS adapter.
- Not yet built: 11b residual (Option B + GPS adapter), ROS node + durable fsync (13b). Partial: validation (14) — NEES/NIS/golden + init-P fix DONE (`70c7d38`); remaining = the CONFIG "tuned"-placeholder sweep + (for strict NEES consistency) a distance-aware covariance model / no-ref correction.

---

## 7. Resume in one line

All numbered roadmap slices (0–14) are addressed; init-P fix + 11b Option A + the per-`n` χ² gate have landed. Remaining work is carve-outs + polish:
- **11b residual** — Option A (single-driving-source bias) DONE + per-`n` χ² gate DONE (`e8491dd`). Deferred: **Option B** (median-coupled multi-source bias — design in `ISSUES.md` Slice-11b DESIGN NOTE) and the concrete **GPS adapter**.
- **13b** (real ROS node + bag round-trip + durable `fsync`) — env-blocked (no ROS on the dev box).
- **Slice-14 CONFIG "tuned"-placeholder sweep**.
- **distance-aware covariance model** (or a no-ref correction) for strict no-ref NEES consistency (counters the predict-only translation Ad-inflation).
Pick one, author the brief from `WORKFLOW.md`'s template, run the cycle in §4. Verify with §2. Done.
