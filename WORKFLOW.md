# Autonomous Workflow

How implementation proceeds for this project. This is the operating contract for the orchestrator and every sub-agent. Keep it current — if a decision changes, update this file in the same step.

## Operating model

Work advances **one step at a time** (a step = one roadmap slice from [`ISSUES.md`](./ISSUES.md), or a setup/data task). Each step is executed by **one dedicated sub-agent**. Steps are **sequential** — never start the next step until the current one is committed and the user has been updated.

### Per-step lifecycle
1. **Brief** — the orchestrator (main thread) authors a detailed task brief for the step's sub-agent: goal, scope, exact interfaces/signatures, math, done-criteria, and the verification command. The brief points the agent at the source-of-truth docs.
2. **Plan-first (TDD)** — the sub-agent writes a short plan *before* coding, then works test-first (red → green → refactor).
3. **Review (orchestrator-driven, file hand-off)** — after the implementer commits, the orchestrator launches a separate **reviewer sub-agent** that writes its findings to `reviews/slice-<n>-findings.md`; the orchestrator then launches a **fix agent** (implementer role) that reads that file and fixes **all** findings, re-verifies green, and commits. (Sub-agents here can't spawn sub-agents or be resumed, so the findings file is the hand-off between the reviewer and the implementer role.)
4. **Verify** — the sub-agent must get a green gate: `powershell -File scripts/dev.ps1 -Task test` (build + all tests pass). No commit on red.
5. **Commit** — the sub-agent commits its own work with a clear message (ending in the standard `Co-Authored-By` line).
6. **Report** — the orchestrator briefly summarizes what happened, proposes the next step, and **waits** for the user's instruction.

### Escalation
If a sub-agent hits a question that needs the user's judgement, it surfaces it; the orchestrator relays it to the user and steers the agent with the answer. Do not guess on user-facing or irreversible decisions.

### Review mechanics
The orchestrator drives the review (step 3) because sub-agents here cannot spawn sub-agents or be resumed:
1. Launch a **reviewer sub-agent** (`caveman:cavecrew-reviewer`) on the slice diff; it **writes `reviews/slice-<n>-findings.md`** (severity-tagged, one finding per line).
2. Launch a **fix agent** (`general-purpose`, implementer role) given that findings-file path; it reads it, fixes every finding, re-runs the green gate, and commits.
3. The orchestrator independently re-runs the gate, updates the source-of-truth docs, summarizes, and waits.
Findings files are kept in `reviews/` as the per-slice review record.

### Sub-agent type
Use a full-capability agent (`claude` / `general-purpose`) for implementation steps — they can spawn their own reviewer sub-agent and run the build. (The `cavecrew-builder` agent refuses 3+ file scope and cannot spawn sub-agents — not suitable for a full slice.)

## Verification gate
- One command: `powershell -File scripts/dev.ps1 -Task test` (configures if needed, builds, runs CTest).
- Green = build succeeds **and** every test passes. This is the bar for every commit.

## Standing decisions (authoritative; mirror of the design)
- **Language/target**: C++14, AUTOSAR profile. Strict core (no heap post-init, no exceptions, bounded WCET, fixed-capacity), relaxed edges (adapters, sim, tests).
- **Dev toolchain (this box)**: Visual Studio 2022 Community (MSVC 19.44). The dev build keeps exceptions on for Eigen/MSVC ergonomics; true `-fno-exceptions`/`-fno-rtti` is the cross-compile target's job.
- **Deps**: Eigen 3.4.0 + doctest 2.4.11, fetched via CMake FetchContent. In-house SO(3)/SE(3) (no Sophus).
- **Test framework**: doctest.
- **Datasets**: automotive (wheel + IMU + GPS). Used via the *real-trajectory + published-extrinsics → synthesized odometry* path to retain calibration ground truth; raw odometry streams for golden/time-sync tests.
- **Build loop**: encapsulated in `scripts/dev.ps1` (auto-discovers VS via vswhere; vcvars → VS-bundled cmake + ninja).

## Source-of-truth documents
- [`DESIGN.md`](./DESIGN.md) — architecture spec.
- [`DECISIONS.md`](./DECISIONS.md) — decision log (chosen / rejected / why).
- [`CONFIG.md`](./CONFIG.md) — every config knob.
- [`ISSUES.md`](./ISSUES.md) — slice roadmap + status.
- `WORKFLOW.md` (this file) — operating model + standing decisions.

## Progress
- [x] Slice 0 — Lie ops, build, doctest harness (green: 14 cases / 25 assertions).
- [x] Slice 1 — Source buffers & uniform delta query (green: 31 cases / 287 assertions).
- [x] Slice 2 — Median fusion + ESKF integrator (first tracer bullet) (green: 56 cases / 491 assertions).
- [x] Slice 4 — Histogram primitive (green: 76 cases / 2651 assertions).
- [x] Sim rig (the calibration oracle) — trajectory + synthetic sources + rig driver (green: 98 cases / 4259 assertions). Exposed + fixed a Slice-2 gap (estimator now applies `prior_scale`).
- [x] Slice 5 — Time-sync (‖ω‖ xcorr → offset histogram, commit N_min + hysteresis) (green: 115 cases / 4370 assertions).
- [x] Slice 6 — Phase-1 calibration: straight-gated yaw/pitch (3-ch so(3) @ prior) + per-source scale; observability self-test (green: 124 cases / 4565 assertions).
- [x] Slice 7 — Phase-2 calibration: turn-gated roll (S¹) + xyz lever-arm (hand-eye LS), both strategies, vote_weight honored, observability self-test (green: 138 cases / 5031 assertions).
- [x] Slice 8 — Commit + feedback loop: per-DOF commit gate (mass + hysteresis), atomic swap, **contractive** re-anchor (extrinsic/scale/offset converge from large priors), cold-start switch (green: 148 cases / 5160 assertions). **Calibration spine complete.**
- [x] Slice 3 — Lifecycle & degrade-don't-block: `INIT→WARMUP→DEGRADED→NOMINAL` ladder driven by `n` vs `min_sources_warn`, `readiness` scalar, reference-only dead-reckon emitted as early as the reference has span, graceful downgrade on source loss (green: 156 cases / 5247 assertions).
- [x] Slice 9 — Weight refinement: variance-EMA reliability (bias/variance split, D17) — noisy source downweighted, biased source kept (bias left to the calibrator); `reliability_floor`/`reliability_cap`; `SourceHealth.reliability`/`bias` diagnostics (green: 161 cases / 5317 assertions).
- [~] Slice 14 — Validation harness (partial): NEES Monte-Carlo consistency + golden regression added; **finding: published Σ is ~46x PESSIMISTIC** (ensemble NEES ~0.13 vs DOF 6 — init `P=I12` + right-error Ad inflation, predict-only has no correction to shrink). NIS deferred to Slice 11; CONFIG placeholder sweep remains (green: 165 cases / 5809 assertions).
- [x] Slice 11 — Absolute-ref correction path: `Eskf::update` (Mahalanobis-gated, Joseph, right-error full-SE(3) injection) + `add_correction` wired into `step()` + `Result::CorrectionDiag`; sim `SyntheticAbsoluteRef` removes drift (0.58→0.20 m), gate rejects outlier (NIS ~3e5), **NIS now computable** (~2.4 vs DOF 3). Bias states + GPS adapter → Slice 11b (green: 175 cases / 7705 assertions).
- [x] Slice 10 — Per-sensor fixed-lag RTS twist smoother (D18): CV forward + backward RTS, deeper calibration frontier `now − delay − L`; variance ↓~0.3×, **zero-phase** verified; calibration peaks ~4× sharper (conf 0.038→0.156) with no estimate bias; `per_sensor_kf` OFF byte-identical; dropout time-alignment fixed via push-seq stamping (green: 187 cases / 8084 assertions).
- [x] Slice 12 — Warm-restart persistence (D23): core `serialize`/`deserialize` into fixed buffers (`persistence.hpp`, "OFCP" v1, FNV-1a config-hash + checksum + orthonormality guard); re-anchor-and-refill restore resumes near-NOMINAL (warm ~0.001 m vs cold ~0.256 m); config-hash/checksum/version/crash-mid-write all reject; histogram bins NOT persisted (re-fill via hysteresis); file double-buffer is relaxed-edge (test); production adapter → Slice 13 (green: 195 cases / 8274 assertions).
- [~] Slice 13 — Adapters (buildable subset; ROS → 13b): `adapters/` tree + `ofc_adapters` (relaxed, core PUBLIC API only, no new deps; gate builds+tests via `dev.ps1 -DOFC_BUILD_ADAPTERS=ON`). File-persistence double-buffer (validity-based overwrite target — torn higher-seq can't clobber last-good; CRITICAL caught+fixed in review), threading wrapper (mutex-guarded snapshot, determinism-tested), dep-free config loader (subset → Config, dup-id/reference cross-check). ROS node + true fsync deferred to 13b (green: ctest 2/2 — unit 195/8274 + adapters 20/319).
- [x] Init-P covariance fix (closes the Slice-14 finding's primary cause): init `P = 0` at the gauge-anchored first fuse, let the first `predict()` establish the one-window covariance `blkdiag(q_pose, q_pose/dt²)` (was `I₁₂`); NEES ~0.13 → ~0.35 (~46× → ~17× pessimistic). Residual = predict-only translation Ad-inflation (needs a distance-aware cov model for strict consistency). Reviewed; no new deps/knobs (green: unit 195/8274 + adapters 20/319).
- [~] Slice 11b Option A — per-source bias states (augmented ESKF): a single `bias_states` source driving the predict alone augments to 18-DOF `[pose;twist;bias]`; predict de-biases (`Δ∘exp(-b·dt)`) + builds the pose↔bias cross-cov (`J_pb=-dt·I₆`); the absolute-ref update removes the bias (sim: planted bias recovered, drift 16 m → 0.06 m); no-ref observability self-test; multi-bias guard; `predict_aug_frozen` out-of-regime; default-OFF byte-identical; new knob `bias_process_noise`. Option B (median-coupled) + per-n χ² gate + GPS adapter deferred (green: unit 200/8322 + adapters 20/319).
- [x] Slice 11b — per-`n` χ² Mahalanobis gate (carve-out, `e8491dd`): `Eskf::chi2_gate(base, n)` scales the n=3-tuned `mahalanobis_chi2` base by the χ²-quantile ratio `q[n]/q[3]` (const `kChi2Q95`, 0.95 confidence) so every measurement DOF n∈1..6 gates at the same confidence; estimator passes `chi2_gate(cfg.mahalanobis_chi2, m.dim)` to both `update()`/`update_aug()` (raw-threshold signatures unchanged). n=3 byte-identical → dim=3 fixes unchanged; load-bearing once a 6-DOF/mixed plugin lands. Resolves the Slice-11 single-scalar-gate limitation (green: unit 206/8352 + adapters 20/319).
- [x] Slice 11b — concrete GPS adapter (carve-out, `4bce8d1`): `ofc_adapters::GpsCorrection` (relaxed edge) turns geodetic GPS fixes (lat/lon/alt + ENU cov) into the dim=3 position correction — WGS-84 geodetic→ECEF→ENU about a datum (configured or lazy-latched) → odom via a configurable `odom_from_enu` SE3; lever-arm `h=t+R·l`, `H=[R|−R·[l]×|0]`; latest-wins emit-once gate. Self-contained adapters test: drift 0.133→0.0075 m + 25 m outlier rejected at the dim=3 per-`n` gate. Surfaced the GPS↔`q_floor` gain coupling (CONFIG §3). 11b residual now = Option B only (green: unit 206/8352 + adapters 27/372).
- [x] Slice 14 — median-variance adaptive-Q reduction (approach A, `730bcfa`): a read-only spike attributed the ~17x predict-only no-ref NEES pessimism to the spread-derived Q over-stating the FUSED median's accuracy (spread = inter-source disagreement; median of n agreeing sources has ~1/n the variance), amplified by the Ad translation transport. Fix = `adaptive_q_source_reduction` (default on) divides the spread term by the participating-median source count `n_eff` -> NEES 0.35 -> ~1.0 (~17x -> ~6x), principled + safe (never overconfident); knob-OFF reproduces the old value. Residual ~6x = the un-calibrated q_scale coefficient (the CONFIG sweep). NIS 2.4 -> 2.7. Golden unchanged (spread=0 noise-free path) (green: unit 207/8537 + adapters 27/372).
- [x] Slice 14 — `q_scale` covariance calibration sweep (`983fa65`): swept `q_scale` {1.0..0.1} × {nees_traj, mixed, turning, straight} × {1×,2× noise} against the NEES harness; `q_scale 1.0 -> 0.5` chosen SAFETY-FIRST (max pessimism cut subject to NEVER overconfident, conservative margin) -> no-ref NEES ~1.0 -> ~2.07 (worst-case ~2.9/~3.5, never >6). Deliberately NOT pushed to NEES=6 (overfits one trajectory; sim under-states real model mismatch). `q_floor` default kept small (GPS-gain trade is a documented deployment choice). Permanent multi-trajectory `test_cov_calibration` guard (never-overconfident + pessimism-reduced); bands re-tuned; golden unchanged (green: unit 208/9154 + adapters 27/372).
- (remaining slices per `ISSUES.md`).

## Sub-agent task-brief template
```
ROLE: Implement <Slice N: title> for OdometryFusionCalibration.
READ FIRST: DESIGN.md §<x>, DECISIONS.md D<n>, CONFIG.md §<x>, WORKFLOW.md.
GOAL: <one sentence>.
SCOPE: files you may touch = <list>. Do not change public headers unless the brief says so.
INTERFACES: <exact signatures / data structures to implement>.
MATH/ALGORITHM: <formulas, edge cases, references>.
METHOD: plan first; TDD (write failing tests, then implement); keep the strict-core profile.
DONE-CRITERIA: <observable, testable conditions> AND `scripts/dev.ps1 -Task test` is green.
REVIEW: after implementing, spawn a reviewer sub-agent; fix ALL findings; re-verify.
COMMIT: commit your work (message + Co-Authored-By line). Report what you did, any deviations, and open questions.
```
