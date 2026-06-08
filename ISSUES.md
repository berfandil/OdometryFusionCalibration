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

## Slice 11b — Per-source bias states + GPS adapter  `[~]` (Option A + per-n gate + GPS adapter done; Option B deferred)
**Goal**: classic loosely-coupled GPS/INS drift removal via online bias estimation.
- **Option A DONE** (`48ece29` + `f47b29a`): a SINGLE `bias_states=true` source that drives the predict alone augments the ESKF to **18-DOF** `[pose;twist;bias(6)]`; predict de-biases (`Δ∘exp(-b·dt)`) + builds the pose↔bias cross-cov (`J_pb = -dt·I₆`); the absolute-ref update removes the bias via that cross-cov (sim: planted bias recovered, drift 16 m → 0.06 m; no-ref observability self-test; multi-bias guard → `InvalidConfig`; out-of-regime `predict_aug_frozen` keeps the filter consistent). New knob `SensorConfig::bias_process_noise`; `CalibSnapshot.bias` + `bias_observable`. Default-OFF byte-identical.
- **Per-`n` χ² Mahalanobis gate DONE** (`e8491dd`): `Eskf::chi2_gate(base, n)` scales the n=3-tuned `mahalanobis_chi2` base by the χ²-quantile ratio `q[n]/q[3]` (const `kChi2Q95[7]`, 0.95 confidence, DOF 1..6) so every measurement DOF n∈1..6 gates at the same confidence; the estimator passes `chi2_gate(cfg.mahalanobis_chi2, m.dim)` to both `update()`/`update_aug()`. n=3 returns base unchanged (the dim=3 fixes shipping now are behaviorally identical — load-bearing once a 6-DOF/mixed `ICorrection` lands; n=6 ≈ 1.61×base). Signatures unchanged (still take a raw per-n threshold).
- **GPS adapter DONE** (`4bce8d1`): `ofc_adapters::GpsCorrection` (relaxed edge) — a production `ICorrection` turning geodetic GPS fixes (lat/lon/alt + ENU cov) into the dim=3 position correction. WGS-84 geodetic→ECEF→ENU about a datum (configured or lazy-latched on first fix) → odom via a configurable `odom_from_enu` SE3; lever arm `h(x)=t+R·l`, `H=[R|−R·[l]×|0]`; `R_odom=R_align·(cov_enu+floor·I)·R_alignᵀ`; latest-wins emit-once pending-fix gate. End-to-end (self-contained adapters test): tail drift 0.133 m → 0.0075 m (~94% removed, 89 fixes), gross 25 m outlier Mahalanobis-rejected at the dim=3 per-`n` gate (NIS ~2.5e5). `GpsConfig` is adapter-local (documented in `adapters/README.md`); no new core knobs. **CAVEAT**: GPS Kalman gain couples to `q_floor` — with two near-identical sources the adaptive-q spread is ~0 so `q_pose`→`q_floor`; a too-small `q_floor` keeps P from re-inflating between fixes and fixes barely pull (see CONFIG §3 note).
- **Deferred (11b residual)**: Option B only (median-coupled multi-source bias — see the DESIGN NOTE).
**Done when**: with `bias_states` on, a raw-IMU source's bias is observed and removed in sim (✓ Option A); a concrete GPS adapter removes fused drift end-to-end (✓ `4bce8d1`).
**Deps**: Slice 11.

### DESIGN NOTE (deferred 2026-06-07 — decide A vs B when revisiting)
**The hard part / why deferred.** D22 motivates bias states by *classic loosely-coupled GPS/INS*: the biased sensor **drives the predict**, so an absolute-ref (GPS) correction of the pose observes the accumulated bias through the **pose↔bias cross-covariance** that the predict builds up. But this architecture's predict is driven by the **geometric median of N frame-aligned source deltas**, not by one sensor — so a single source's bias only reaches the pose through its (robust, nonlinear, weight-dependent) contribution to the median. The cross-covariance that makes the bias observable therefore does **not** form the classic way; how the bias couples into the median-driven predict is the open modeling question. Two candidate approaches were identified (both promising — pick one when revisiting):

- **Option A — clean-regime augmented filter (lower risk).** Augment the ESKF state to `[pose(6); twist(6); bias_i(d_i)...]` with a dense augmented covariance (fixed-size, compile-time cap on total bias DOF). Implement the classic mechanism — predict random-walks each bias and builds the pose↔bias cross-covariance; the absolute-ref `update()` removes it — **only for the regime where it is exact**: a *single driving source* (e.g. `ColdStart::ReferenceOnly` dead-reckon, or a 1–2 source rig), where that source's delta IS the predict so the GPS/INS structure holds. Deliver a sim test (a raw-IMU-like biased source + GPS fixes → bias observed & removed) and **document** that N-source-median bias coupling is deferred. Correct + testable + bounded; marks 11b `[~]`.
- **Option B — full median-coupled bias (higher risk).** Model the bias→pose coupling through each source's **median weight** in the augmented predict `F` (the Jacobian block from bias_i to the pose-error is ∝ that source's effective weight in the Weiszfeld solution), so the bias is observable with N sources in the median. Most faithful to D22's general case, but the coupling Jacobian through the robust median is subtle and is exactly where an implementation would silently go wrong — needs careful derivation + observability self-tests; not a one-cycle drop-in.

**State sizing.** Augmented dim = `12 + Σ d_i` over sources with `bias_states=true`. Worst case (all 32 sources, 6-DOF IMU bias) = 12+192 = 204 → a fixed `Mat<204>` is large but strict-core-legal (allocate once); realistically `bias_states` is on for ≤1–2 sources, so a compile-time cap on total bias DOF (e.g. `kMaxBiasDof`) sized for the expected use keeps the augmented `P` small. Bias is off by default (core stays bias-free + agnostic).

**Recommendation when revisiting**: land **A** first (the genuine GPS/INS win in the clean regime, low risk), then evaluate whether **B** (the N-source coupling) is worth the modeling effort for the target automotive rigs.

**UPDATE (`1142e41`) — Option B RECONSIDERABLE.** Spikes 5a/5b had ruled B a NO-GO because the per-source median-influence weights `ωᵢ` were degenerate `(1,0,…,0)` — but that was a symptom of the D3 **median-pinning bug** (the solver returned one vertex, so only one source had influence). With the median fix (true interior blend), the `ωᵢ` are now REAL distributed weights, so B's coupling Jacobian `J_{pose,bᵢ} = −dt·ωᵢ·Ad(Xᵢ)` (derived + FD-verified in spike 5a) now carries genuine multi-source observability. The NO-GO reason is gone; B is a legitimate (still subtle, still higher-risk) future option. The `ωᵢ` are extractable from the estimator's existing per-source `split_distance` to `med.value` + weights. Note the related calibrator-consensus-contamination item (Slice 14) — a true blend means a source's bias/calibration now reaches the consensus, which B's modeling must account for.

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

## Slice 14 — Validation harness  `[~]` (sim rig + observability + NEES + NIS + golden + init-P + median-variance Q reduction + q_scale calibration DONE; remaining placeholders + distance-aware cov SHAPE / NHC for strict no-ref NEES)
**Goal**: the trust apparatus (D24).
- Sim rig (parameterized GT trajectory + sources), observability self-tests (per regime), NEES/NIS Monte-Carlo, recorded-data golden regression.
- **No-ref covariance — NEAR-CONSISTENT (`1142e41`)**: ensemble-mean pose NEES ≈ **4.82** vs DOF 6 (worst-case ~4.85 across trajectories, never overconfident). Got here via: init-P fix (`70c7d38`, 0.13→0.35); then the D3 **median-pinning fix** (the solver was returning ONE source verbatim — see DECISIONS D3) which made the spread correctly sized; then **`q_scale` recalibrated 0.5→0.7** against the true median (the safety-first close-to-6-from-below value). `test_cov_calibration` guards never-overconfident + near-consistency across trajectories.
- **Approach A (`/n_eff`) REMOVED** (`730bcfa` … reverted in `1142e41`): it was calibrating around the D3 bug (false "median variance reduction" rationale — the pinning median averaged nothing); with the true median it over-divided Q → overconfident. `q_scale` alone now calibrates the covariance. (DECISIONS D4.)
- **Remaining**: the non-covariance placeholders (`excitation_min_var`, `kf_process_noise`, `match_metric`, straight/turn calibration gates — separate pass, observability self-tests pin them functionally); and the `Ad` distance *shape* model and/or an NHC no-ref correction for strict per-trajectory consistency to exactly DOF=6.
- **Latent bugs surfaced by the median fix (were MASKED by the pinning bug)**: (1) **cold-start scale-calib spurious commit — FIXED (`3bd91e2`)**: root cause was `scale_hist` inheriting the generic `[-1,1]` default → a unit scale ratio (1.0) hit the half-open upper boundary → last-bin clamp → `mode()` returned the boundary-bin center 0.984375 (the "63/64" = bin 63 of 64), so a true unit scale committed as 0.984. Fix = default `scale_hist` to `[0.5,1.5]` (1.0 interior) + a `validate()` guard; test_sim re-enabled `scale_calib` (dropped the workaround). The noise-free golden `c1->scale` corrected 1.08281→1.1 (1.08281 = 1.1×0.984375 — independent confirmation). (2) **partial scale recovery (~1.07 vs 1.2) — NOT a bug, documented known-limitation**: the scale vote is vs the REFERENCE (not the median), and recovers 1.2 fully when votes are dense; the ~1.07 is a sparse-vote (27/1305 gated) + ±20%-per-window-noise jitter in the SlidingK ring (the median fix only changed *which* steps gate-fire). `test_weights`'s relaxed `>1.05` directional assertion stays; levers if ever needed = denser straight-regime data or `ReferenceOnly` cold-start. No calibrator change.
**Done when**: CI runs unit + observability + consistency + golden; tuned defaults replace the "tuned" placeholders in `CONFIG.md` (the load-bearing `q_scale` is done; the rest are a non-covariance follow-on pass).
**Deps**: grows alongside Slices 2–11.

---

## Suggested order
`0 → 1 → 2` (tracer bullet) → `4` (parallel) → `3` → `5 → 6 → 7 → 8` (calibration spine) → `9, 10, 11` (refinements) → `12, 13` → `14` runs throughout.
