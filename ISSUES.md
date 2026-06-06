# Implementation Roadmap

Independently-grabbable work items as **tracer-bullet vertical slices** — thinnest end-to-end thing first, then layer capability. Each slice ships something runnable + tested. References: [`DESIGN.md`](./DESIGN.md), [`DECISIONS.md`](./DECISIONS.md), [`CONFIG.md`](./CONFIG.md).

Status legend: `[ ]` todo · `[~]` in progress · `[x]` done.

---

## Slice 0 — Skeleton & primitives  `[ ]`
**Goal**: the strict-core foundation compiles and the build splits core / adapters / sim / tests.
- Core types: fixed-size aliases, `Status`/`Expected`/`Optional` (no-except, C++14), `double`-based vectors/matrices.
- In-house **SO(3)/SE(3)** ops: `exp`/`log`/`Adj`/compose/inverse (+ unit tests vs analytic).
- `Config` struct tree + `validate()` (bounds/capacity → `Status`).
- CMake: `ofc_core` static lib (strict flags), optional `adapters`/`sim`/`tests`.
**Done when**: core lib builds with strict flags; Lie-op + `validate()` unit tests pass.
**Deps**: none.

## Slice 1 — Source buffers & uniform delta query  `[x]`
**Goal**: ingest a source and answer `delta(t0,t1) → (SE3, Σ)`.
- Fixed-capacity ring buffer per source (preallocated).
- `ISource` interface + adapters for native forms: twist→integrate, increments→compose, absolute→difference; in-window interpolation.
- Σ: native ⊕ modeled combine (D7); missing native → identity.
**Done when**: integration/interpolation/Σ-combine unit tests pass against planted streams.
**Deps**: Slice 0.

## Slice 2 — Median fusion (FIRST tracer bullet)  `[x]`
**Goal**: end-to-end fused output from N sources, no calibration yet.
- Weiszfeld split-metric geometric median (bounded iters, ε-reg).
- ESKF predict on median delta; **adaptive Q** from spread; pose+twist state, SO(3)×ℝ³ error, dense 12×12.
- Caller-pumped `step()`; `Result` struct (frontier state + predicted tip).
- Frame-align deltas using **prior** extrinsics (calibration still fixed).
- Add `confidence_blend` to the `Config` struct (CONFIG.md §4 knob; from Slice 1 it exists only as a `SourceBuffer::configure()` parameter).
**Done when**: sim rig (known motion, N sources) → fused trajectory tracks GT within noise; NEES in-bounds; deterministic replay byte-stable.
**Deps**: Slice 1. **This is the vertical slice everything else hangs off.**

## Slice 3 — Lifecycle & degrade-don't-block  `[ ]`
**Goal**: `INIT → WARMUP → DEGRADED → NOMINAL`, odom anchored at first tick, readiness in `Result`.
- Reference-only output during DEGRADED; cold-start switch (D6).
- Graceful downgrade on source loss.
**Done when**: lifecycle transitions unit-tested; output appears as early as the reference source allows.
**Deps**: Slice 2.

## Slice 4 — Histogram primitive  `[x]`
**Goal**: the shared robustness estimator (reused by every calibrated quantity).
- Fixed bins, linear-split voting, parabolic sub-bin, peak-concentration confidence.
- Aging: `decay` and `sliding_k` (config). Circular variant (for roll).
**Done when**: mode-recovery + aging + circular-wrap unit tests pass.
**Deps**: Slice 0. (Can proceed in parallel with 2–3.)

## Slice 5 — Time-sync  `[x]`
**Goal**: per-source clock-offset estimation.
- ‖ω‖ resample + pluggable-metric cross-correlation, parabolic sub-sample, excitation gate.
- Vote into per-source offset histogram; feed offset back into the delta query.
**Done when**: planted offset recovered in sim (within sub-sample tol); flat-‖ω‖ correctly skips.
**Deps**: Slices 1, 4.

## Slice 6 — Calibration Phase 1 (yaw/pitch + scale)  `[x]`
**Goal**: straight-regime recovery.
- Straight gate (+ reverse-fold, optional ref cross-check).
- 3-channel **so(3) histogram @ per-sensor prior** → mode → exp → yaw,pitch.
- Scale = magnitude ratio vs pinned reference → scale histogram.
**Done when**: planted yaw/pitch/scale recovered in sim; **observability self-test**: values stay at prior unless straight motion is present.
**Deps**: Slices 4, 2.

## Slice 7 — Calibration Phase 2 (roll + xyz)  `[x]`
**Goal**: turn-regime recovery via hand-eye.
- Turn gate; strategy A (`vs_fused_base`): roll 1-D residual + xyz linear-LS.
- Then strategy B (`pairwise_pinned_ref`); compare empirically.
**Done when**: planted roll/xyz recovered under turning; **observability self-test**: frozen under pure-straight; LS singularity at `R_A=I` handled.
**Deps**: Slice 6.

## Slice 8 — Commit, hysteresis & feedback loop  `[x]`
**Goal**: close the calibration→fusion loop.
- Commit on `τ_commit` ∧ `N_min` with hysteresis; atomic swap of extrinsics/scale/offset into fusion between steps.
- Bootstrap convergence (prior → refined) verified end-to-end.
**Done when**: from wrong-ish priors, sim converges extrinsics to GT and fused accuracy improves; no commit thrash.
**Deps**: Slices 5, 6, 7.

## Slice 9 — Weight refinement  `[ ]`
**Goal**: variance-EMA reliability with bias routed to calibration (D17).
- Floored/capped reliability; systematic residual → calibrator, not weight.
**Done when**: a noisy source is downweighted; a *biased* source is calibrated (not just downweighted) in sim.
**Deps**: Slices 2, 8.

## Slice 10 — Per-sensor fixed-lag RTS smoother  `[ ]`
**Goal**: two-sided pre-calibration smoothing at the deeper frontier (D18).
- CV twist ESKF forward + backward RTS pass over lag `L`.
**Done when**: calibration histogram peaks sharpen vs raw input; no peak shift (zero-phase) in sim.
**Deps**: Slices 1, 6.

## Slice 11 — Absolute-ref plugin + optional bias states  `[ ]`
**Goal**: correction path + classic GPS/INS drift removal.
- `ICorrection` interface, Mahalanobis-gated ESKF update.
- Optional per-source bias states; GPS adapter (in `adapters/`).
**Done when**: sim with GPS fixes removes pose drift; with `bias_states` on, a raw-IMU source's bias is observed and removed.
**Deps**: Slice 2.

## Slice 12 — Persistence (warm restart)  `[ ]`
**Goal**: serialize/deserialize calibration state; survive crashes.
- Core `serialize/deserialize` into fixed buffers; file double-buffer + version + checksum in adapter; config-hash guard.
**Done when**: restart resumes near-NOMINAL; crash mid-write keeps last good state; config change invalidates.
**Deps**: Slice 8.

## Slice 13 — Adapters  `[ ]`
**Goal**: relaxed-edge integrations. YAML/JSON config loader, ROS node (Odometry/TF/calib msgs), threading wrapper, file-persistence backend.
**Done when**: each adapter builds against the core public API; ROS node round-trips on a recorded bag.
**Deps**: Slices 2, 12.

## Slice 14 — Validation harness  `[~]` (sim rig built; observability self-tests + NEES/NIS + golden remain)
**Goal**: the trust apparatus (D24).
- Sim rig (parameterized GT trajectory + sources), observability self-tests (per regime), NEES/NIS Monte-Carlo, recorded-data golden regression.
**Done when**: CI runs unit + observability + consistency + golden; tuned defaults replace the "tuned" placeholders in `CONFIG.md`.
**Deps**: grows alongside Slices 2–11.

---

## Suggested order
`0 → 1 → 2` (tracer bullet) → `4` (parallel) → `3` → `5 → 6 → 7 → 8` (calibration spine) → `9, 10, 11` (refinements) → `12, 13` → `14` runs throughout.
