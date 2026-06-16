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

**UPDATE (`1142e41`) — Option B RECONSIDERABLE.** Spikes 5a/5b had ruled B a NO-GO because the per-source median-influence weights `ωᵢ` were degenerate `(1,0,…,0)` — but that was a symptom of the D3 **median-pinning bug** (the solver returned one vertex, so only one source had influence). With the median fix (true interior blend), the `ωᵢ` are now REAL distributed weights, so B's coupling Jacobian `J_{pose,bᵢ} = −dt·ωᵢ·Ad(Xᵢ)` (derived + FD-verified in spike 5a) now carries genuine multi-source observability. [CORRECTION 2026-06-11, Slice 18: spike 5a's FD verification was against the PINNING median; the scalar form is falsified against the true interior median — the exact matrix influence Ω_i replaces it, see SLICE18_MEDIAN_COUPLED_BIAS.md §2 + tests/test_multi_bias.cpp.] The NO-GO reason is gone; B is a legitimate (still subtle, still higher-risk) future option. The `ωᵢ` are extractable from the estimator's existing per-source `split_distance` to `med.value` + weights. Note the related calibrator-consensus-contamination item (Slice 14) — a true blend means a source's bias/calibration now reaches the consensus, which B's modeling must account for.

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
- **Config loader**: dependency-free key=value/INI → `Config` (documented knob subset; owns `vector<SensorConfig>`; duplicate-id + reference cross-check; runs core `validate()`). Now also exposes the **process-noise knobs** `q_scale`/`q_floor`/`adaptive_q` (`e03b237`, the real-data covariance enabler — `q_floor` accepts 1 number → all 6 axes, or 6 → [trans;rot]). **Known gap**: the `ofc_replay` manifest reader rejects full-line `#`/`;` comments (only inline) — minor CLI robustness.
- **GPS correction adapter** (`GpsCorrection`, `4bce8d1`): geodetic GPS → WGS-84 ECEF → ENU → odom dim=3 position `ICorrection` (lever arm; datum/alignment config).
- **Real-data CSV ingestion + replay** (`ddcf436`): `CsvSource : ISource` (generic CSV, 3 forms: absolute/increment/twist, hand-rolled dep-free parser, backed by `SourceBuffer`), `CsvGtTrack` (interpolated GT lookup), `ReplayHarness` (pumps `Estimator::step()` over a tick timeline; computes drift + 6-DOF NEES + GPS-fix NIS vs GT — the real-data analogue of the sim Rig), and the `ofc_replay` CLI (manifest = extended config-loader INI). Equivalence-tested: CSV ingestion == in-memory to 1e-6. **Unblocks real-dataset testing** (KITTI/KAIST/NCLT/EuRoC → CSV via a one-time conversion script). Schema + manifest in `adapters/README.md`.
**Deferred → Slice 13b**: real ROS node + recorded-bag round-trip (no ROS on the dev box; a compile-guarded header sketch ships); a true `fsync` (the relaxed-edge adapter uses `flush`/`close`, leaving an OS-page-cache power-loss window).
**Deps**: Slices 2, 12.

## Slice 13b — ROS node + durable fsync  `[ ]`
**Goal**: the platform integrations that the dev box can't build/test.
- Real ROS node (Odometry/TF subs → `ISource`; GPS/map-match → `ICorrection`; `Result` → calib/odom/TF msgs); round-trip on a recorded bag.
- Replace the persistence adapter's `flush`/`close` with a real `fsync(fd)` (platform layer) to close the power-loss durability window.
**Deps**: Slice 13.

## Slice 14 — Validation harness  `[~]` (sim rig + observability + NEES + NIS + golden + init-P + median-variance Q reduction + q_scale calibration DONE; real-data tested KITTI+KAIST; remaining = non-cov placeholders + a real-data covariance recalibration (q_floor↓ predict + GPS R↑) — see below)
**Goal**: the trust apparatus (D24).
- Sim rig (parameterized GT trajectory + sources), observability self-tests (per regime), NEES/NIS Monte-Carlo, recorded-data golden regression.
- **No-ref covariance — NEAR-CONSISTENT (`1142e41`)**: ensemble-mean pose NEES ≈ **4.82** vs DOF 6 (worst-case ~4.85 across trajectories, never overconfident). Got here via: init-P fix (`70c7d38`, 0.13→0.35); then the D3 **median-pinning fix** (the solver was returning ONE source verbatim — see DECISIONS D3) which made the spread correctly sized; then **`q_scale` recalibrated 0.5→0.7** against the true median (the safety-first close-to-6-from-below value). `test_cov_calibration` guards never-overconfident + near-consistency across trajectories.
- **Approach A (`/n_eff`) REMOVED** (`730bcfa` … reverted in `1142e41`): it was calibrating around the D3 bug (false "median variance reduction" rationale — the pinning median averaged nothing); with the true median it over-divided Q → overconfident. `q_scale` alone now calibrates the covariance. (DECISIONS D4.)
- **Remaining**: the non-covariance placeholders (`excitation_min_var`, `kf_process_noise`, `match_metric`, straight/turn calibration gates — separate pass, observability self-tests pin them functionally); and a **real-data covariance recalibration** [PREDICT SIDE DONE 2026-06-12: q_floor swept on KAIST full drives, recommended urban value 1e-5/1e-6 (urban12 transient 409->129 m, all drives net-better, CONFIG #3); the per-drive consistent floor differs 11x -> scalar floor cannot be cross-drive consistent, Ad-shape remains; R side (GPS cov_floor vs urban multipath) still open]. Real-data findings (KITTI 28 drives + KAIST urban07/12/17) reframed this: the covariance error is a MAGNITUDE problem on EXISTING knobs, not a missing shape term. (a) **GPS-corrected overconfidence is R-dominated** — the GPS measurement noise R (~1 m²) is far tighter than the consumer-GPS-vs-SLAM-GT urban error (>10 m); bigger R (`GpsConfig.cov_floor_m2` 0→100) drops NEES 687→87, and a too-tight R also makes the Mahalanobis gate reject valid GPS (urban12 death-spiral: 16173/16673 rejected → 2 km runaway). (b) **predict-only is `q_floor`-magnitude-tunable** — overconfident at the tiny default (KITTI q_floor=1e-6 → NEES ~400) but PESSIMISTIC at 1e-3 (KAIST → NEES 0.085); the right value is between. A **distance-aware predict-Q SHAPE term (`q_dist`) was tried and REVERTED** (`4fc2be3`→`cc38598`): real data showed it's the wrong lever (q_dist had ~0 effect in the R-dominated GPS regime; made the already-pessimistic predict worse). So the open work = recalibrate `q_floor` (predict, down-ish) + GPS `R` (up) on real data — both EXISTING magnitude knobs; an `Ad`-shape / NHC term is at most a second-order refinement, not the primary lever. See memory `real-dataset-testing`.
  - **Motion-aware predict-Q DE-RISKED → NO-GO (2026-06-16, `tools/motion_stats.py`)**: tested whether a motion-proportional Q collapses the per-drive predict-only NEES spread (urban07 62.6 vs urban17 10.4 = 6× at q_floor=1e-7, split config). **Stage-1**: the per-drive honest-q_floor tracks **yaw_rate²** (urban07/17 ratio 8.4×) NOT speed² (anti-correlated) — predict error is ROTATION-driven (the `Ad` transport). **Stage-2** (re-added the reverted diagonal `add_distance_q` as an env-gated prototype, reverted after): a single per-step `(k·|Δθ|)²` rotation-Q coefficient does **NOT** collapse the spread — it OVER-corrects the lower-turning drive (qdist_rot=25 → urban07 1.28 vs urban17 0.082 = 15×, WORSE than 6× baseline; urban17 paradoxically more sensitive). The aggregate yaw_rate² correlation does not carry to a per-step DIAGONAL term because the `Ad` transport accumulates non-linearly. **Verdict**: the naive diagonal motion-Q is confirmed-wrong (matching the earlier revert); a true `Ad`-transport-SHAPED Q (the coupled covariance transport, not diagonal) is the only principled form, a heavy derivation for a SECOND-ORDER refinement on an already-SAFE (never-overconfident, pessimistic) predict-only NEES → DEFERRED, not worth it now. The real-data covariance frontier is the R-side (GPS-vs-GT frame offset, `GPS_R_NEES_SWEEP.md`), not the predict-Q shape.
- **Latent bugs surfaced by the median fix (were MASKED by the pinning bug)**: (1) **cold-start scale-calib spurious commit — FIXED (`3bd91e2`)**: root cause was `scale_hist` inheriting the generic `[-1,1]` default → a unit scale ratio (1.0) hit the half-open upper boundary → last-bin clamp → `mode()` returned the boundary-bin center 0.984375 (the "63/64" = bin 63 of 64), so a true unit scale committed as 0.984. Fix = default `scale_hist` to `[0.5,1.5]` (1.0 interior) + a `validate()` guard; test_sim re-enabled `scale_calib` (dropped the workaround). The noise-free golden `c1->scale` corrected 1.08281→1.1 (1.08281 = 1.1×0.984375 — independent confirmation). (2) **partial scale recovery (~1.07 vs 1.2) — NOT a bug, documented known-limitation**: the scale vote is vs the REFERENCE (not the median), and recovers 1.2 fully when votes are dense; the ~1.07 is a sparse-vote (27/1305 gated) + ±20%-per-window-noise jitter in the SlidingK ring (the median fix only changed *which* steps gate-fire). `test_weights`'s relaxed `>1.05` directional assertion stays; levers if ever needed = denser straight-regime data or `ReferenceOnly` cold-start. No calibrator change.
**Done when**: CI runs unit + observability + consistency + golden; tuned defaults replace the "tuned" placeholders in `CONFIG.md` (the load-bearing `q_scale` is done; the rest are a non-covariance follow-on pass).
**Deps**: grows alongside Slices 2–11.

## Slice 15 — Robust correction update (large-innovation safety)  `[~]` A shipped opt-in; D validated; urban12 NOT solved
**Goal**: bound the damage a single large-innovation GPS correction can do. The KAIST `urban12` divergence (4214 m global) is NOT pervasive — the local metric (`7c71b63`) shows its median 10 s window is 0.17 m (best of the three drives); the km-scale tail comes from a HANDFUL of catastrophic correction windows (`local max` 856 m). Mechanism: with lever=0, `H` rotation cols are zero but `K`'s rotation rows (fed by the `P` trans↔rot cross-cov) still inject a large heading kick `exp(dx[3..5])` on a large position residual → the continuous filter carries the corruption forever.
**Shipped** (`44e0711`): Huber-robust gain (A) in core `update()`/`update_aug()` — down-weight gain ∝ `kappa/dbar` for outlier NIS by inflating active `R`; `correction_robust_kappa=0` default → bit-identical to non-robust (opt-in). ConfigLoader exposes `correction_robust_kappa` + `mahalanobis_chi2`. Unit-tested (yaw kick cut >50% on a large innovation, sub-linear growth, inert below threshold).
**OUTCOME — A did NOT fix urban12; D (realistic GPS `R` via `cov_floor_m2`) is the keeper.** Sweep (urban12 full): **D alone** (cov_floor_m2=25, tight gate) NEES 59907→**225**, local max 856→**476**, tail 4214→3075, accepted-NIS healthy 0.88 — a real CONSISTENCY win (config-only). **A (Huber) is the wrong lever here**: with the tight gate its trigger band (`d2>kappa²·n`) sits *inside* the gate's reject band (`d2>9`) → inert; opening the gate so it can fire admits genuine multipath (admitted-NIS 38≫3) → net WORSE (tail 4605, worst window 1951), Huber claws back NEES (66k→7.5k) but divergence persists. **The gate is correctly rejecting bad fixes.** A kept as a legit opt-in robustness primitive (loose-gate / non-GPS corrections), NOT as the urban12 fix.
**urban12 DIAGNOSED** (`tools/diagnose_window.py` on `urban12_k0f25_out.csv`): the divergence is pinned to ONE applied position fix at **t=1929.2s** — `rot_err` jumps 0.975→2.162 rad (+68°) in one step while translation is pulled 51 m closer; NIS=**8.45** (just under the gate 9). Trigger: a **522 s GPS-denied coast** drifts DR to ~650 m / ~1 rad, then the returning fix over-rotates heading via the `P` trans-rot cross-cov → death spiral (later fixes are genuine 800 m outliers, correctly gated). A (Huber) can't fix: culprit dbar=1.68, so any practical kappa misses it AND A scales the *good* translation pull too. → **Slice 15b (lever C).**
**Deps**: Slices 13 (local metric), 14 (cov calibration).

## Slice 15b — Bounded heading injection from position-only corrections (lever C)  `[x]` DONE — fixes urban12
**Goal**: stop a position fix from over-rotating heading via the `P` trans-rot cross-cov when its residual is large (the pinned urban12 t=1929 s 68° kick). **Shipped** (`ba895b6`): **C4 — residual-gated rotation-row suppression**. For a rotation-unobserving correction (`H` rotation cols ≈ 0) with `dbar > correction_rot_suppress_kappa`, scale ONLY `K`'s pose-rotation rows (3..5) by `kappa/dbar` (translation/twist/bias untouched, Joseph reuses the modified K → consistent + PSD). `correction_rot_suppress_kappa=0` default → opt-in/inert; orthogonal to Slice-15 Huber. ConfigLoader key added; unit-tested (yaw cut >70%, translation kept, inert below threshold + for rotation-observing fixes). **OUTCOME — first lever that FIXES urban12**: D+C4(κ=0.8) tail **4214→1.99 m**, death spiral averted (gps_applied 96→1375), t=1929 kick eliminated (rot 2.0→0.08), arc recovers to ~2 m. NO regression on urban07/17 (C4 inert there — their dbar 0.18–0.55 < κ=0.8); sim bit-identical at κ=0; validated on FULL drives via the local metric. Residual: a mid-drive transient (max ~300–400 m) that now RECOVERS (full flattening needs the upstream 522 s GPS-coast heading-drift, out of scope). Recommended: `cov_floor_m2=25` + `correction_rot_suppress_kappa=0.8`. Full sweep in `SLICE15B_BOUNDED_HEADING_INJECTION.md`. **Deps**: Slices 13, 15.

## Slice 16 — Centroid sub-bin readout (calib precision floor)  `[x]` DONE — committed yaw 0.229°→0.0102°
**Goal**: the user precision target (rotation < 0.1°) on the committed calibration. Investigation (2026-06-10) discriminated the Slice-15c committed-yaw gap: NOT vote weighting (one ≈ combo to 3e-5 rad), NOT commit dynamics — `Histogram1D::mode()`'s parabolic sub-bin interpolation pulls the read ~70–80% toward the peak-bin center (probes: truth at a bin center reads +0.0045°, at quarter-bin reads −0.264°; same pull visible in scale + time-offset of the same run). The historical "0.10°" was luck (10° ≈ a bin center); real floor ~0.3–0.45° at the so3 bin width 1.79°.
**Shipped** (`c237bfb` + review fixes `00978c8`): `HistogramConfig::subbin_centroid` (default false) — mass-weighted centroid over peak±1, EXACT for split votes at any sub-bin position; wrap-aware, one-sided at boundaries; in the persistence config-hash (pre-16 blobs cold-start, intentional); loader `[global] subbin_centroid` fans out to all five calib histograms. Default-off bit-identical (pinned by exact equality); sim/trust apparatus untouched.
**Validated (KAIST urban07 calib_test, flag ON)**: yaw err sub-bin-phase-INDEPENDENT +0.015° across all probes; **committed yaw 0.0102°** (target <0.1°), committed scale err ~8e-5, time-offset err **2.6 µs** (was 2.2 ms), committed lever ~2.5 mm. Design doc `SLICE16_CENTROID_SUBBIN_READOUT.md`; review `reviews/slice-16-findings.md`; DECISIONS D25.
**Deps**: Slices 8 (commit/re-anchor), 15c (`vote_weight=one` commit enabler).

## Slice 17 — Turn-regime full rotation-extrinsic (rot3d axis-correspondence hand-eye)  `[x]` DONE — 9.1e-6° on drone, committed
**Goal**: recover the FULL rotation extrinsic on 3D/multi-axis motion (the EuRoC gap: conf 0 on drone; the wrong rotation also corrupted the lever rows ~|t|·θ). **Math validated by prototype first** (`tools/proto_rot_handeye.py`): hand-eye axis correspondence `a = R_X·b` over windowed deltas + Wahba/Kabsch; KAIST yaw-only → `Σbbᵀ` rank 1 with the missing yaw EXACTLY the about-axis blind DOF → complementary to Phase-1 by construction (the spine); gate = two distinct axes (`λ_mid ≥ 1e-2·λ_max`).
**Shipped** (`f1eee59` + review fixes `2d22c28` + lever-reset `f013320`): Phase-2 per-slot decay-aged `Mw`/`BBw`, running Kabsch vote into 3 so(3) channels (reuse `so3_hist`, contractive rising-edge re-anchor), `rot3d_committed` commit/feedback (supersedes ext/roll publish; joins the median under `ReferenceOnly`), lever rows driven by the running R̂ once the gate opens + lever-accumulator reset on the rising edge, `rot3d_enabled` default OFF (byte-identical), persistence format v2. Review: 2 MAJOR (lever-row coupling, ReferenceOnly join) + 3 MINOR + 2 NIT — all fixed, mutation-pinned. Unit 239/14591 + adapters 40/687 green.
**Validated**: EuRoC rotation extrinsic err **9.1e-6° COMMITTED**; lever **µm-exact** with unit-scale sources; KAIST planar byte-identical with flag ON (gate honest). **Finding → next-slice candidate**: with an unrecovered source scale the lever reads ~3.5 cm off (raw `t_B` in the rows; scale unobservable on drone — no straight windows); the row is linear in `(t_X, 1/s)` → turn-regime JOINT lever+scale LS would close the last calib DOF on 3D motion. D26; `SLICE17_TURN_REGIME_ROTATION_EXTRINSIC.md`.
**Deps**: Slices 7/8 (Phase-2, commit machinery), 16 (centroid readout for the precision claims).

## Slice 17b — Turn-regime joint lever+scale (4-unknown hand-eye LS)  `[x]` DONE — last calib DOF closed; ALL precision targets met
**Goal**: per-source scale without the straight regime (drone-unobservable; an unrecovered scale degraded the Slice-17 lever to ~3.5 cm). **Math**: the hand-eye translation row is LINEAR in `[t_X; κ=1/s_res]` — prototype `tools/proto_lever_scale.py` (EuRoC lever 1.8e-8 mm + scale 2.7e-12 vs the 3-unknown control's 88 mm; KAIST planar κ well-conditioned).
**Shipped** (`4138c3d` + review fixes `0ea764a`): Phase-2 4×4 joint ridge LS behind `joint_lever_scale` (default OFF, byte-identical), marginal-information (Schur) vote gating (`kJointMarginMin`, ≥3 distinct excitations or conf 0), scale2 histogram + out-of-range skip guard + `scale2_skipped` diagnostic, either-path scale fold into `prior_scale` (both hists + joint accumulator reset, no double-fold), persistence v3. Review: 3 MAJOR (mixed near-null vote bias; silent out-of-range; guard untested) + 4 MINOR + 1 NIT — all fixed, mutation-pinned. Unit 255/16899 + adapters 41/693 green.
**Validated (full recipe)**: EuRoC **scale 1.07962 committed from the turn regime** + **lever sub-mm with scale unknown**; KAIST turn-scale agrees with straight-scale (1.10008), lever ~10× tighter. **ALL calib DOFs now meet the precision targets on both regimes.** D27; `SLICE17B_JOINT_LEVER_SCALE.md`.
**Deps**: Slice 17 (rot3d R̂ in the rows), 16 (centroid), 8 (commit machinery).

## Slice 18 — Median-coupled multi-source bias (11b Option B)  `[x]` landed opt-in + hardening; urban12 heading goal NOT met (observability boundary, documented)
**Goal**: fused heading tracks the best source instead of the middle drifter (urban12: fused 5–30° wander vs FOG 4.6°; diagnostic `tools/diag_heading_drift.py`). **Course-correction en route**: the spike-5a scalar coupling was FD-FALSIFIED against the true interior median (it had been verified against the old pinning median); the exact matrix influence Ω_i = M⁻¹uᵢPᵢ (implicit function theorem at the Weiszfeld fixed point, ΣΩᵢ=I, FD-pinned 1.5% worst) replaced it — `tests/test_multi_bias.cpp` (`def0cee`).
**Shipped** (`2dc6678` + `67b0242` + `57fc874` + harness diag `9756d46`): k≤4 augmented filter (Mat36), single de-bias site (fusion+calib see the same deltas), `multi_bias_enabled` default OFF (Option A byte-identical), per-DOF `bias_process_noise` (0 pins a DOF), `multi_bias_cov0` knob; **out-of-range calib-vote SKIP guards (so3/xyz/rot3d + skip counters)** — closes the measured corruption channel where wandering bias rotated deltas until votes edge-clamped and COMMITTED at ±0.984375 (urban12 tail 1.99→2716 m); guards active regardless of the flag. Review 4 MAJOR + 8 MINOR + 3 NIT all fixed; the MAJOR-1 "scale bug?" challenge falsified by a 5-point noise-scale sweep (linearization-radius signature, not multiplicative). Unit 274/17238 + adapters 43/732 green. Sim: separation, coast 4.4×/14×, NEES never-overconfident.
**Real-data verdict (KAIST, instrumented)**: sim-scale priors → junk sink (one GPS fix kicks yaw bias to +9870°/h, divergence); honest physical priors → INERT (learns ±3°/h of −75°/h in 41 min, bias_observable 0.000, run ≡ flag-off to 5 decimals). Structural: position-only GPS observes yaw-rate bias only through the rotation gain rows C4 deliberately suppresses; the cheap channel (source-vs-source rate differences) is what the median discards by design; source-as-measurement = D4's rejected double-count. **Heading prize needs a D3-level lever (per-DOF consensus weighting) — future decision.** urban07/12 no-regression at physical scale ✓. D28; `SLICE18_MEDIAN_COUPLED_BIAS.md`.
**Deps**: 11b Option A, 15b (C4), median fix (`1142e41`).

## Slice 19 — Per-channel split median (D3 amendment, D29)  `[x]` DONE — fused heading now tracks the best heading source
**Goal**: per-channel source quality in fusion (the urban12 heading prize Slice 18 proved unreachable by de-biasing). USER-approved forks: split per channel; weight policy = config priors only (layer a); cross-channel veto opt-in default-ON.
**Shipped** (`6ebfcc5` + review fixes `6c23894`): `solve_split` — independent SO(3)/ℝ³ Weiszfeld with the FULL D3 safeguard set per channel; per-channel spreads → per-channel Q (`q_scale_split=3.0`, own calibration — the coupled 0.7 is overconfident under split, NEES 21.7 pinned; the coupled mixed-unit spread had padded both Q blocks); LOO-floored cross-channel veto (full-set RMS was dead code, unfloored LOO a hair-trigger on the shared-wheel rig); `rot_weight_prior` per sensor; Ω_i block-diagonal FD-pinned 1.1%; default OFF byte-identical (golden untouched). Review 2 MAJOR (veto floor; q_scale footgun → own knob) + 5 MINOR + 2 NIT all fixed. Unit 305/18911 + adapters 44/752 green.
**Validated (KAIST, urban config + FOG rot_weight_prior=10)**: urban12 local rot p50 0.0133→**0.0041 rad (0.23°)**, rms rot **5.4° ≈ the FOG floor** (was a 5–30° band with a 61° spike); urban07 better everywhere; urban17 neutral. One urban12 worst-window regression 120→211 m logged. **Follow-ups**: ~~layer (b) per-channel scatter reliability~~ **DONE — Slice 19b** (`51b86bf`); ~~layer (c) GPS-course drift monitor~~ **DONE — Slice 19c** (`770f6b2`+`e7bba55`); the one-window regression; veto floors vs ≫1 s windows; flip split default after soak. D29; `SLICE19_SPLIT_CHANNEL_MEDIAN.md`.
**Slice 19b — per-channel scatter reliability (layer b)** `[x]` DONE (`51b86bf`): under split, per-channel residuals vs the split consensus feed two Welford-EMA tracks → `rel_rot`/`rel_trans` per source (`w_trans = clamp(prior·rel_trans·σconf)`, `w_rot = clamp(prior·rel_rot·σconf)·rot_weight_prior`); mixed scalar EMA frozen under split; `SourceHealth` +`reliability_rot/trans`+`bias_rot/trans` (legacy fields = trans channel under split); no new knobs; coupled path byte-identical (pre-change exact pin). Per-channel state NOT persisted (re-warms ~20 fused steps post-restore — format bump only if warm-restart reliability ever matters for split rigs). Known boundary: a source holding ≥ half a channel's mass self-validates (reliability caps) — priors/veto remain the majority-lock-in levers. Review APPROVE no findings (`reviews/slice-19b-findings.md`). Unit 311/19216 + adapters 44/752 green.
**Slice 19c — GPS-course heading-drift monitor (layer c)** `[x]` DONE (`770f6b2` + review fixes `e7bba55`): strict-core `HeadingMonitor` ranks per-source heading drift vs GPS course-over-ground → boosts the heading-grade source's ROTATION-channel weight (`w_rot ·= monitor_boost`), auto-discovering the FOG with no `rot_weight_prior`. Ranking-not-bias ⇒ sidesteps the D28 boundary. Anchor-gated telescoping pairs (bridge turns, never dropped) → median residual split → per-GPS-contiguous-segment 60-s-**median**-block slope → `boost = clip(boost_max·max(min,floor)/max(score,floor),1,boost_max)`, `kScoreFloor=5°/h`, ≥2-scored abstain + rank-flip hysteresis. `heading_monitor` (default false; requires `split_median`; validate() rejects otherwise) + `heading_monitor_boost_max` (10); both hashed; state not persisted; diag in `SourceHealth.{heading_score,heading_boost,heading_scored}`. OFF byte-identical (flag-off AND abstain-ON pinned). Review REQUEST CHANGES (block mean→median CRITICAL; zero-drift boost floor CRITICAL; abstain pin) all fixed. **Validated (KAIST, `rot_weight_prior` REMOVED)**: FOG auto-discovered on all 3 drives (boost 10 vs 1, scores rank-correct); <10° met everywhere; urban17 neutral; cost = discovery latency (urban12 GPS-hostile canyon start → p50 0.23°→0.49° during the ~1300 s abstain) ⇒ layer (c) for UNKNOWN/drifting heading source, config prior (a) when known. Follow-up: warm-start boost from `rot_weight_prior` to erase the latency. Unit 323/37157 + adapters 45/764 green. `SLICE19C_HEADING_MONITOR.md`; review `reviews/slice-19c-findings.md`.
**Deps**: median fix (`1142e41`), Slice 18 (Ω_i harness + the boundary evidence), q_floor recal (`e37256b`), Slice 19b (`51b86bf`, per-channel weights it composes with).

## Slice 20 — Translation-only source lever (velocity/Doppler sensors, D30)  `[x]` DONE (`629d379`)
**Goal**: let a velocity-only source (Doppler radar; `R_B=I`) calibrate its lever + stop it committing a garbage rotation extrinsic. **Motivated by** the nuScenes radar test (5 radars heading-blind → lever conf 0 + garbage rotation). Investigated (the lever LS is already right; only `R_X` falls back to Phase-1 and gets poisoned) + prototyped (`tools/proto_radar_lever.py`: lateral lever recovers ~cm with a trusted `R_X`, along-track is a sensor SNR wall).
**Shipped**: `SensorConfig::translation_only` (default false) pins the rotation extrinsic to the prior (skip Phase-1 direction / rot3d / roll) so the lever row uses `R_X = prior`; keeps scale; along-track withheld by the existing concentration gate (no residual gate). Byte-identical off; `validate()` needs an orthonormal prior rotation; hashed. Review APPROVE clean (`reviews/slice-20-review.md`). 8 unit cases; gate green.
**Real-data (nuScenes radar)**: footgun fix CONFIRMED (all radars hold `extr_rot=prior`, conf 0, no garbage) — the verified win + a real safety fix. But the LEVER does NOT commit on real data (conf 0, readout overshoots lateral): `translation_confidence=min(cx,cy,cz)` pins on the planar-unobservable z + the production per-window/consensus-base estimator is noisier than the proto's idealized GT-base batch. The proto over-promised; the full path revealed the gap (validate-on-the-real-path lesson). **Follow-ups**: `min`-over-OBSERVED-axes lever confidence (drop the unobservable z); lower per-window radar noise / batch lever-only solve at the deeper frontier; along-track stays a sensor SNR wall. D30; `SLICE20_TRANSLATION_ONLY_LEVER.md`.
**Slice 20b — per-axis lever commit** `[x]` DONE (`c85442e`, review APPROVE `reviews/slice-20b-review.md`): for `translation_only` sources the lever commits PER-AXIS (3 independent `commit_gate_reanchor` on per-axis confidence+mass; committed axes fold into `prior.t[c]`, others held; `lever_committed[i]`=OR; `CalibSnapshot.translation_committed_xyz[3]`; persistence v3→v4, pre-v4 cold-start). Non-flagged byte-identical. Sim-validated (cx 0.61 commits where the whole min-over-observed 0.17 could not). **Real radar: `lever_committed_xyz=[0,0,0]` — still no commit** (every axis ~0.55 < 0.6 gate). Chain now complete (footgun + min-over-observed + per-axis commit, all correct/sim-validated); the ISOLATED remaining wall = raw per-window Doppler noise (readout overshoots ~20-28 cm → don't lower the gate). Committing a real radar lever is now a radar-front-end odometry-quality problem, not a calibration one. D30 §4b; `SLICE20_TRANSLATION_ONLY_LEVER.md`.
**Deps**: nuScenes pipeline (`799dcb8`/`5f68121`) + the radar heading-blind finding.

## GPS R-side / with-correction consistency (sweep done 2026-06-15, `GPS_R_NEES_SWEEP.md`)  `[~]`
**Done**: swept `cov_floor_m2` {25..3200} × urban07/12/17 on the recommended config. Findings: `cov_floor_m2=25` is **17–300× overconfident** (with-GPS NEES 103/1889/1198 vs 6); raising R reduces NEES + (urban07) IMPROVES local drift (over-tight R injected GPS jitter); honest-R differs ~7× cross-drive (1500 vs 10000) → scalar R can't be cross-drive consistent; NIS<<3 while NEES>>6 = structural; the honest-R is implausibly large for VRS-RTK → dominated by **GPS-vs-GT frame disagreement** (~8–26 m), not noise; the KAIST per-fix variance is a flat 0.25 m² nominal (no multipath signal).
**Open follow-ups**: (1) raise the default `cov_floor_m2` to ~200–400 (config-only, strictly better — CONFIG updated); (2) a per-fix / multipath-aware **R MODEL** (the flat VRS variance has no signal — needs modelling or a varying-variance GPS); (3) characterize + separate the **GPS-vs-GT frame offset** (alignment/lever/SLAM-drift) from the measurement noise; (4) an **un-gated-NIS** R-calibration methodology (NEES-vs-GT is confounded by the frame offset; gated NIS is biased low). Parallel to the `Ad`-shape `q_floor` (the predict-side cross-drive analogue).

---

## Suggested order
`0 → 1 → 2` (tracer bullet) → `4` (parallel) → `3` → `5 → 6 → 7 → 8` (calibration spine) → `9, 10, 11` (refinements) → `12, 13` → `14` runs throughout.
