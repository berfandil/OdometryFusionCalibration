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

- **Gate green: 156 doctest cases / 5247 assertions.** 35 commits on `main`, **not pushed**, working tree clean.
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

The **calibration spine (5–8) is complete** — calibration closes back into fusion and bootstraps from arbitrary priors.

- **Remaining slices** (any order; recommended: **9 → 14 → 10 → 11 → 12 → 13**):
  - Slice 9 — weight refinement: variance-EMA reliability, bias→calibration (D17).
  - Slice 10 — per-sensor fixed-lag RTS smoother (two-sided, deeper frontier).
  - Slice 11 — absolute-ref plugin (Mahalanobis-gated) + optional per-source GPS/INS bias states.
  - Slice 12 — warm-restart persistence (double-buffer + config-hash guard).
  - Slice 13 — adapters (YAML/ROS/threading/file-persistence).
  - Slice 14 — finish validation: NEES/NIS consistency + recorded-data golden regression (per-slice observability self-tests already exist).

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
- **ESKF**: state = pose `SE(3)` + twist `ℝ⁶`, error `[trans;rot]` (pose 0–5, twist 6–11), dense 12×12. Predict `F = blkdiag(Ad(delta⁻¹), 0)`, `P ← F P Fᵀ + blkdiag(Q, Q/dt²)`. Predict interval = `[last_frontier, frontier]` (gap/overlap-free); `window_s` is bootstrap/lookback only.
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
- Not yet built: weight-refine (9), RTS smoother (10), absolute-ref/bias (11), persistence (12), adapters (13), NEES/NIS+golden (14).

---

## 7. Resume in one line

Pick a slice (suggest **9 weight-refinement**), author the brief from `WORKFLOW.md`'s template, run the cycle in §4. Verify with §2. Done.
