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
