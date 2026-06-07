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

## Slice 3 — Lifecycle & degrade-don't-block  `[x]`
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

## Slice 9 — Weight refinement  `[x]`
**Goal**: variance-EMA reliability with bias routed to calibration (D17).
- Floored/capped reliability; systematic residual → calibrator, not weight.
**Done when**: a noisy source is downweighted; a *biased* source is calibrated (not just downweighted) in sim.
**Deps**: Slices 2, 8.

## Slice 10 — Per-sensor fixed-lag RTS smoother  `[x]`
**Goal**: two-sided pre-calibration smoothing at the deeper frontier (D18).
- CV twist ESKF forward + backward RTS pass over lag `L`.
**Done when**: calibration histogram peaks sharpen vs raw input; no peak shift (zero-phase) in sim.
**Deps**: Slices 1, 6.

## Slice 11 — Absolute-ref correction path  `[x]`
**Goal**: correction path — absolute refs remove fused pose drift.
- `ICorrection` interface, **Mahalanobis-gated ESKF update** (`Eskf::update`, Joseph form, right-error full-SE(3) injection).
- `Estimator::add_correction` wired into `step()` (after predict, before publish); `Result::CorrectionDiag` (evaluated/applied/rejected/last NIS); `validate()` `mahalanobis_chi2 > 0`.
- Sim `SyntheticAbsoluteRef` (position fix) for tests. **NIS now computable** (closes the Slice-14 NIS deferral; ensemble ~2.4 vs DOF 3, mildly conservative).
**Done when**: sim with absolute fixes removes pose drift (✓ 0.58→0.20 m tail error); a gross outlier fix is Mahalanobis-rejected (✓ NIS ~3e5).
**Deps**: Slice 2.

## Slice 11b — Per-source bias states + GPS adapter  `[ ]`
**Goal**: classic loosely-coupled GPS/INS drift removal via online bias estimation.
- Optional per-source bias states (augment the fixed core state; `SensorConfig::bias_states` is a no-op today); absolute-ref updates observe/remove the bias through filter cross-covariance (D22).
- Concrete GPS adapter (in `adapters/`, Slice 13 territory).
**Done when**: with `bias_states` on, a raw-IMU source's bias is observed and removed in sim.
**Deps**: Slice 11. **Note**: gate single-scalar `mahalanobis_chi2` should become a per-DOF chi² quantile when a 6-DOF/mixed plugin lands.

## Slice 12 — Persistence (warm restart)  `[x]`
**Goal**: serialize/deserialize calibration state; survive crashes.
- Core `serialize/deserialize` into fixed buffers (`persistence.hpp` format primitives: "OFCP" v1, explicit-LE/padding-free, FNV-1a-64 config-hash + FNV-1a-32 checksum + orthonormality guard); config-hash guard rejects a changed rig; new Status `CorruptData`/`VersionMismatch`.
- **Re-anchor-and-refill restore** (committed values + flags + reliability + lifecycle restored, calibrators re-anchored to the restored priors; histogram BINS not persisted — they re-fill, held committed by the hysteresis guard). File double-buffer ping-pong is relaxed-edge (test-side here; production adapter is Slice 13).
**Done when**: restart resumes near-NOMINAL (warm ~0.001 m vs cold ~0.256 m); crash mid-write keeps last good state; config change invalidates. ✓
**Deps**: Slice 8.

## Slice 13 — Adapters  `[~]` (buildable subset done; ROS → 13b)
**Goal**: relaxed-edge integrations. YAML/JSON config loader, ROS node (Odometry/TF/calib msgs), threading wrapper, file-persistence backend.
**Done (subset)**: `adapters/` tree + `ofc_adapters` target (relaxed flags, links the core PUBLIC API only, no new deps; gate builds+tests it via `OFC_BUILD_ADAPTERS=ON` in `dev.ps1`).
- **File-persistence**: production double-buffer ping-pong on the Slice-12 core serialize/deserialize; **validity-based** overwrite-target selection (preserves the highest-seq blob `load()` accepts; a torn higher-seq file can't clobber the last-good) → crash-mid-write keeps last good.
- **Threading wrapper** (`ThreadedEstimator`): a worker pumps `step()`, mutex-guarded `Result` snapshot; determinism vs single-thread reference.
- **Config loader**: dependency-free key=value/INI → `Config` (documented knob subset; owns `vector<SensorConfig>`; duplicate-id + reference cross-check; runs core `validate()`).
**Deferred → Slice 13b**: real ROS node + recorded-bag round-trip (no ROS on the dev box; a compile-guarded header sketch ships); a true `fsync` (the relaxed-edge adapter uses `flush`/`close`, leaving an OS-page-cache power-loss window).
**Deps**: Slices 2, 12.

## Slice 13b — ROS node + durable fsync  `[ ]`
**Goal**: the platform integrations that the dev box can't build/test.
- Real ROS node (Odometry/TF subs → `ISource`; GPS/map-match → `ICorrection`; `Result` → calib/odom/TF msgs); round-trip on a recorded bag.
- Replace the persistence adapter's `flush`/`close` with a real `fsync(fd)` (platform layer) to close the power-loss durability window.
**Deps**: Slice 13.

## Slice 14 — Validation harness  `[~]` (sim rig + observability self-tests + NEES + NIS + golden DONE; init-P covariance fix DONE; CONFIG placeholder sweep + a distance-aware covariance model for strict no-ref NEES consistency remain)
**Goal**: the trust apparatus (D24).
- Sim rig (parameterized GT trajectory + sources), observability self-tests (per regime), NEES/NIS Monte-Carlo, recorded-data golden regression.
**Done when**: CI runs unit + observability + consistency + golden; tuned defaults replace the "tuned" placeholders in `CONFIG.md`.
**Deps**: grows alongside Slices 2–11.

---

## Suggested order
`0 → 1 → 2` (tracer bullet) → `4` (parallel) → `3` → `5 → 6 → 7 → 8` (calibration spine) → `9, 10, 11` (refinements) → `12, 13` → `14` runs throughout.
