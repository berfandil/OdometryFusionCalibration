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

- **Gate green: 187 doctest cases / 8084 assertions.** 47 commits on `main`, working tree clean. Remote synced through Slice 11; the Slice-10 commits are local (not pushed).
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

The **calibration spine (5–8) is complete** — calibration closes back into fusion and bootstraps from arbitrary priors.

- **Remaining slices** (any order; recommended: **12 → 13**, plus **11b**):
  - Slice 11b — per-source bias states (augment the core state) + GPS adapter; deferred from Slice 11. The path to the classic GPS/INS bias-removal and to making the gate per-DOF.
  - Slice 12 — warm-restart persistence (double-buffer + config-hash guard).
  - Slice 13 — adapters (YAML/ROS/threading/file-persistence).
  - Slice 14 — `[~]` partial: NEES consistency + golden + **NIS** DONE (Slice 11 made NIS computable). Remaining: the CONFIG "tuned"-placeholder sweep (and the init-P covariance fix is the way to make NEES/NIS strictly consistent).

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
- **Slice-14 finding — published covariance is PESSIMISTIC (~46×) on the predict-only path**: Monte-Carlo NEES (6-DOF pose, full-SE(3) `se3::log` tangent) is ensemble-mean ≈ **0.13** vs DOF **6** — over-conservative, not over-confident. Root cause: ESKF inits `P = I₁₂` (≈100× the steady-state error) and a predict-only stretch has no correction to shrink P; the right-error `Ad(delta⁻¹)` propagation further inflates the translation block over distance. **Not** fixable via `q_scale`/`q_floor` (those only add to P). Fix path: seed init-P from the first-window/median uncertainty (≪ 1.0 m²/rad²). **Slice 11 partially mitigates this when an absolute ref is present** (the correction step shrinks P → NIS ≈ 2.4 vs DOF 3, near-consistent); the no-ref NEES test still trips `CHECK_FALSE(truly_consistent)` so a future init-P fix re-trips it.
- **Slice-11 correction-gate limitation**: the Mahalanobis gate is a **single scalar `mahalanobis_chi2` regardless of measurement DOF `n`** (~97% quantile for n=3, ~80% for n=6). Fine for the dim=3 position fixes shipping now; make it a per-`n` χ² quantile when a 6-DOF/mixed plugin lands (Slice 11b). The drift-removal test uses a loose `chi2=100` — a **test artifact** driven by the covariance pessimism above (legit drift residuals exceed `chi2=9` because P never shrinks on predict-only stretches), NOT a production value.
- **Covariance tangent (doc-vs-code reconciled)**: the ESKF covariance is **full coupled SE(3)** (propagated by the full `Ad(delta⁻¹)`), NOT block-diagonal "decoupled SO(3)×ℝ³" — that phrase only describes the error ordering `[trans;rot]` + the median's split metric. NEES uses the full `se3::log` to match. (DESIGN §5, D21 now corrected.)
- **Slice-10 smoother scope**: a single shared `TwistSmoother` uses the MAX `kf_process_noise` over enabled sources (not per-source `q`/`r`; `r` fixed 1.0). The refined RTS covariance is computed + exposed but NOT wired into the calibrator vote weight (the deeper path still feeds the raw Σ-confidence) — D18's "refined Σ" is half-wired. Fixed ~7.5 MB footprint (compile-time caps 32×65, paid even when off). The lag `L` is a step count via the nominal tick rate → off-cadence pumping gives an effective time-lag ≠ `calib_lag_s` (rings stay step-aligned regardless). Per-source dropout degrades to a variance loss (no bias) after the push-seq alignment fix.
- Not yet built: per-source bias states + GPS adapter (11b), persistence (12), adapters (13). Partial: validation (14) — NEES/NIS/golden done, CONFIG "tuned"-placeholder sweep + the init-P covariance fix pending.

---

## 7. Resume in one line

Pick a slice (suggest **12 persistence**, or **13 adapters**, or **11b bias states**, or a small **init-P covariance fix** to make NEES/NIS strictly consistent), author the brief from `WORKFLOW.md`'s template, run the cycle in §4. Verify with §2. Done.
