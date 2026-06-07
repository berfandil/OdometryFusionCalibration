# OdometryFusionCalibration — Design

A lightweight, AUTOSAR-C++14, dependency-free **library** that:

1. **Fuses** N odometry sources into one robust motion estimate (a geometric-median-driven error-state integrator), and
2. **Self-calibrates** their relative extrinsics, scale, and time-offsets **online** — from odometry disagreement alone, no external reference required.

One robustness primitive runs throughout: **geometric median** for fusion, **histogram-mode** for everything calibrated.

---

## 1. Architecture principles

- **Platform-agnostic core, thin plugins.** Generality lives at the edges (source/correction plugins, adapters), never in the core. This is the resolution of the "agnostic vs lightweight" tension.
- **Error-state, Lie-group.** State on a manifold; errors in the tangent space. ESKF, *not* IEKF or a smoother (deliberately — bounded compute/memory).
- **Separation of concerns.** Fusion (fast, causal, up-to-date) and calibration (slow, latency-tolerant) are distinct subsystems sharing the buffers. Calibration is **not** augmented into the filter state.
- **Determinism.** Caller-pumped, single-threaded, lock-free, bounded WCET. Same input → identical output, every run.

---

## 2. Dataflow

```
per-source ring buffers (fixed capacity, preallocated)
  → TIME-SYNC  [config on/off]
        cross-correlate ‖ω‖(t) per source (extrinsic-invariant → runs before calibration)
        pluggable match metric (L1 / L2 / ratio / NCC ...), excitation-gated
        → per-source offset histogram → mode
  → FUSION  (causal frontier = now − delay_fast)
        integrate each source over window W → SE(3) delta
        scale-correct (B.t / prior_scale) then frame-align (A = X∘B∘X⁻¹)
        weighted geometric median (Weiszfeld, split SO(3)/ℝ³ metric, bounded iters) → u_k
        ESKF: predict on u_k; Q from inter-source spread;
              correct ONLY on optional absolute-ref plugins (Mahalanobis-gated)
        OUTPUT: frontier state (pose_odom, twist, Σ)  +  predict-only tip @ now
  → CALIBRATION  (deeper frontier = now − delay_fast − L, latency-tolerant)
        per-sensor fixed-lag RTS smoother (two-sided)         [config per-sensor on/off]
        Phase 1 — STRAIGHT-gated (‖ω‖<ε_ω ∧ ‖t‖>δ_v, reverse-folded,
                  optional reference-sensor cross-check):
              forward dir → 3-channel so(3) histogram @ per-sensor-prior basepoint
                          → mode → exp → forward axis → yaw,pitch
              translation magnitude ratio vs pinned ref → scale histogram → s_i
        Phase 2 — TURN-gated (‖ω‖>θ_ω): fix yaw,pitch → hand-eye A·X = X·B
              roll = 1-D rotation residual; xyz = linear-LS  (R_A−I)t_X = R_X t_B − t_A
              votes → roll (S¹) + xyz (3×1-D) histograms
              strategy ∈ { vs-fused-base | pairwise-with-pinned-ref }   (A/B compared)
        commit extrinsic/scale/offset when peak-concentration ≥ τ_commit
              AND votes ≥ N_min, with hysteresis (re-open below τ_drop)
        → feeds back into time-sync + fusion (atomic swap between steps)
```

---

## 3. The observability spine (organizing principle)

Each calibration DOF is observable in exactly one motion regime — and the math hands you the gate for free (the lever-arm least-squares goes singular precisely when there is no rotation).

| DOF | Observable when | Recovered by |
|---|---|---|
| yaw, pitch (forward axis) | **straight** (‖ω‖≈0 ∧ ‖t‖>0) | 3-channel so(3) histogram @ prior basepoint |
| **scale** (per-source translation) | **straight** (pure translation → no lever-arm confound) | magnitude-ratio histogram vs pinned reference |
| **roll** (about forward axis) | **turning** (‖ω‖>0) | hand-eye 1-D rotation residual |
| **xyz lever-arm** | **turning** (better with larger ‖ω‖) | hand-eye linear-LS; singular at R_A=I ⇒ re-derives the gate |
| **time-offset** | ‖ω‖ **varies** (distinctive shape) | ‖ω‖ cross-correlation → offset histogram |

Every gate, vote weight, and reliability term derives from this table.

---

## 4. Fusion

- **Geometric median** (robust to outliers — uses the *unsquared* distance; squaring would give the Karcher mean and lose robustness). Solved by **Weiszfeld IRLS**: reweight by `1/d`, re-solve; bounded `max_iters`; ε-regularized weight to avoid the at-a-vertex singularity. **Split metric**: median rotations (SO(3) via log map) and translations (ℝ³) separately, with a rot/trans weight `λ`.
- Requires **≥3 sources** to actually reject an outlier; with 2 it degenerates to a weighted geodesic midpoint (no rejection); with 1 it is pass-through.
- **Weights** `w_i = prior × reliability × Σ-confidence`, bounded.
  - **Reliability** = variance-EMA of the source's **zero-mean** residual scatter to the consensus; slow, **floored and capped** (a source can never collapse to 0 → always recoverable; never dominates).
  - The **systematic** residual component is routed to the **calibrator** (it's a miscalibration), *not* the weight — avoids masking faults and majority lock-in.
  - **Slice-9 impl**: per source, an EMA mean (`resid_mean` = the systematic component) and an EMA variance *around that mean* (`resid_var` = the zero-mean scatter) of the per-step residual-to-consensus distance (`split_distance`), at rate `reliability_ema_alpha`. A constant bias washes into `resid_mean` and so does **not** inflate `resid_var` — that is the bias/variance split. `reliability = clamp(ref_var / resid_var, reliability_floor, reliability_cap)` with `ref_var` = the robust median of warmed-up participants' `resid_var` (so equal-noise sources read ≈ 1, a noisier source < 1, a cleaner source > 1); applied as `w = prior × reliability × Σ-confidence`. Reliability is **causal** (this step's weight uses the prior steps' track; updated after the median) and holds at 1.0 until a per-source warmup count. The bias is **not** folded into the weight (it is left for the calibrator to absorb via the existing Phase-1/2 observe path); `SourceHealth.bias` exposes it for diagnostics only and is an unsigned residual *magnitude*, not a signed offset.
  - The within-solve Weiszfeld `1/d` handles per-step outliers; the EMA is the slow per-source quality track.
- **Adaptive Q**: process noise from the inter-source spread around the median (tight agreement → small Q; disagreement → large Q).

---

## 5. ESKF (the "robust integrator")

- **State**: pose `T ∈ SE(3)` + body twist `(v,ω) ∈ ℝ⁶`. **Error-state**: the pose error is *ordered* `[trans; rot]` and the **median** uses a split `SO(3)×ℝ³` metric — but the **ESKF covariance is NOT block-diagonal**: predict propagates it with the **full SE(3) adjoint** `F = Ad(delta⁻¹)` (the `[t]ₓR` block couples translation and rotation), so the dense 12×12 covariance lives in the **full coupled SE(3) tangent**. "Decoupled" describes the error *parameterization/ordering* (and the median metric), not the covariance — NEES must be computed in the full-SE(3) `se3::log` tangent to match (Slice-14 validation note). **Dense 12×12 covariance** (fixed-size, no-heap friendly).
- **Predict**: `T_k = T_{k-1} ∘ Δ_median`; twist from the windowed median; const-velocity tip extrapolation to `now`.
- **Correct**: only via optional absolute-reference plugins.
- **Calibration params are NOT in the state** — they live in the external histograms and feed back via frame-alignment.
- **Optional per-source bias states** (config per-source): when enabled (e.g. a raw-IMU source), absolute-ref updates observe and remove the bias through filter cross-covariance — recovers classic loosely-coupled GPS/INS drift correction. Off by default (core stays bias-free and agnostic).

---

## 6. Calibration

- **Gauge**: odometry-only ⇒ extrinsics recoverable only *relative to an anchor*. User declares a **base frame**; one **reference sensor** is pinned (extrinsic = prior) as the gauge. All else relative.
- **Bootstrap**: every extrinsic/scale/offset starts at its config prior; fusion runs with current-best values; calibration refines; refined values **commit and swap into the fusion prior atomically** (between steps). On an extrinsic/scale commit the calibration **basepoint re-anchors** (so(3)/scale histograms reset to vote the *residual* around the new basepoint) — making the extrinsic, scale, and time-offset re-anchors **contractive**, so they converge even from a large prior error (verified ~23° extrinsic). Two-rate loop (fusion fast, calibration slow). Realistic convergence floor ~0.08 rad on mixed straight+turn trajectories (windows straddling a regime boundary spread the mode). Commit gate: `confidence ≥ τ_commit ∧ vote-mass ≥ N_min`, hysteresis (re-open below `τ_drop`).
- **Phase 1 (straight)**: forward direction via **3-channel so(3) histogram** (store the full rotation vector — `φ_x≡0` only holds in the basepoint-aligned frame; numerical noise makes it nonzero anyway — take the so(3) **mode**, then `exp → forward axis → yaw,pitch` at the very end). Basepoint = each sensor's **prior** forward (keeps data near the basepoint: small tilt, no pole/antipode issue, isotropic resolution). Reverse motion folded into the consensus hemisphere before voting — fold uses the **fused/consensus travel sign** (not a per-source axis component); a direction ≥90° off the prior basepoint is **skipped** (outside the small-deviation regime, and avoids the so(3)-log π singularity); with fold off, reverse samples are dropped. Optional reference-sensor cross-check for ice-slide / drone niches. Scale recovered here too (magnitude ratio vs the pinned reference, on the de-scaled delta → residual × prior). **Slice-6 impl note:** the straight gate uses `‖ω‖` as a rate but `‖t‖` as a per-step displacement, so `straight_trans_min` is cadence-dependent (tune per tick rate).
- **Phase 2 (turning)**: fix yaw,pitch, solve **hand-eye** `A·X = X·B` against the base motion. **Roll** is the only nonlinear DOF (1-D; bounded scan + parabolic refine, voted to a circular S¹ histogram); **xyz is linear** given rotation (`(R_A−I)t_X = R_X t_B − t_A`, normal-equation least-squares with an eigenvalue conditioning floor + prior-centred ridge). Two strategies compared empirically: vs-fused-base, and pairwise-between-sensors with the reference extrinsic fixed. **Seam**: Phase-2's roll basepoint is Phase-1's yaw/pitch — so Phase-2 presumes a converged (or prior-correct) Phase-1 (the Slice-8 bootstrap loop makes this hold from arbitrary priors; the R_yp fallback strips prior roll to avoid double-counting). **Lever-arm observability**: pure **yaw-only** rotation leaves `t_X.z` unobservable (z is the null space of `R_A−I` under yaw); recovering all 3 lever-arm components needs **multi-axis** rotation (yaw+pitch) — the conditioning floor keeps the unobservable axis at prior with low confidence.
- **Histograms** (the robustness primitive, a cheap outlier-rejecting alternative to pose-graph optimization):
  - Fixed global bins. so(3) (Phase-1 direction) 3-channel; roll = circular S¹; xyz/scale = 1-D.
  - **Aging**: config **exponential decay** *or* **sliding window of K** (forced by bootstrap — early votes are computed with wrong priors and must wash out; also tracks re-mount/thermal drift).
  - Voting: +1, linear-split between nearest bins, or confidence/rotation-magnitude weighted (configurable).
  - Peak: sub-bin via parabolic interpolation. **Confidence** = peak concentration — (peak bin + its two immediate neighbours) ÷ total (a linear-split peak spans ~2 bins, so peak-bin-only under-reports) — drives commit + cold-start.
  - Binning: 1-D bins are half-open `[range_min, range_max)`; non-circular values clamp, circular values wrap. SlidingK aging clamps tiny-negative bin residuals and rebuilds the running total on eviction (float-drift guard).

---

## 7. Timing & runtime

- **Caller-pumped, single-threaded, lock-free.** Each `step()` does fusion (fast) + a **bounded slice** of calibration work. Threading is an external adapter concern.
- **Two frontiers**: fusion at `now − delay_fast` (causal, up-to-date) + predict-only tip; calibration at `now − delay_fast − L` (deeper, two-sided-smoothed).
- **Predict interval** (Slice 2 impl): the fusion predict advances by the motion over `[last_frontier, frontier]` with `dt = frontier − last_frontier`, so integration is gap/overlap-free regardless of tick cadence. `window_s` is the **bootstrap/lookback interval** for the first fuse (and the max lookback), not a fixed per-tick step.
- **Late samples** behind the frontier: dropped or trigger small re-integration (config). No OOSM machinery.
- **Lifecycle**: `INIT → WARMUP → DEGRADED → NOMINAL`. **Degrade-don't-block** — emit best-available output (reference-only dead-reckon) the instant the reference source has span; auto-upgrade as calibration converges; readiness + per-stage confidence exposed; graceful downgrade on fault. **Slice-3 impl**: a single per-step ladder, recomputed every step from the participating-source count `n` (no high-water-mark needed — downgrade is automatic): `INIT` is set at `init()` and observable before the first step; pre-first-fuse steps (degenerate interval `t1≤q0` or `n==0`) are `WARMUP`; a successful fuse latches `ever_fused` and yields `NOMINAL` iff `n ≥ min_sources_warn` else `DEGRADED`; once fused, a non-fusing step (`n==0` / `t1≤q0`) is `DEGRADED`+`NotReady` and **retains the last good frontier** (never zeroed). NOMINAL is thus source-count-driven; under `ReferenceOnly` cold-start it still tracks calibration convergence indirectly (a non-reference source joins the median — raising `n` — only once its extrinsic commits, via `participates()`).
- **Warm-restart persistence** [config on/off, periodic]: persist **calibration state** (committed extrinsics/scale/offsets, histogram bins, reliability weights, lifecycle confidence) → boot near-NOMINAL. **Double-buffer ping-pong** (files A/B + sequence counter, each version-tagged + checksummed; write inactive → fsync → flip) so a crash mid-write never corrupts the last good state. **Invalidate** on setup change: manual delete + **auto config-hash guard** (reject persisted state whose config/prior signature mismatches the current rig). Core does `serialize/deserialize` into fixed buffers; file IO lives in an adapter.

---

## 8. Interfaces

- **Source plugin**: one query `delta(t0,t1) → (SE3 motion, Σ6×6)`. An adapter converts the native form (twist → integrate, increments → compose, absolute → difference) and interpolates in-window. **Confidence** = Σ: provider's native Σ **⊕** a modeled noise Σ (combine config — default **sum**; native missing → identity).
- **Absolute-ref plugin** (optional, 0..N): `evaluate(state) → (residual = z ⊟ h(x), H, R)` + timestamp. Standard ESKF update at the frontier, **Mahalanobis-gated**. Covers position / pose / orientation agnostically. **Loop-closure excluded** (needs retained past states ⇒ smoother territory). Concrete sensors (GPS, map-match) ship as adapters.
- **Output**: transport-agnostic, **callback + poll**. Per-tick result struct: frontier state (pose-in-odom + twist + Σ), predicted tip (+ inflated Σ), per-sensor calibration snapshot (extrinsic, scale, time-offset, and **per-DOF confidences**: `confidence`=time-offset, `extrinsic_confidence`=yaw/pitch so(3), `scale_confidence`=scale — each DOF converges in its own regime), health/diagnostics (per-source residuals, weights, per-source `reliability` + mean-residual `bias` (Slice 9), gating regime, histogram confidences, lifecycle phase, and a coarse `readiness` scalar in [0,1] — Slice 3: Init 0.0 / Warmup 0.25 / Degraded 0.6 / Nominal 1.0). Pose lives in an **odom frame anchored at init** (drifts — inherent to odometry-only). Extrinsics reported in the base frame. ROS/DDS/zmq/logging are adapters outside core.

---

## 9. Platform & safety (AUTOSAR C++14)

- **Strict core**: no heap after `init()` (everything sized from the config struct at construction); no exceptions (status-code + hand-rolled `Optional`/`Expected`); every loop bounded (Weiszfeld iters, calibration slice, linear-LS); fixed-capacity containers; `double`, no fast-math.
- **Relaxed edges**: adapters + unit tests use std containers, exceptions, file IO — they never ship to the ECU.
- **No `std::optional`/`variant`/`string_view`** (C++17) → explicit types / hand-rolled equivalents.
- **Eigen** for linear algebra; **minimal in-house SO(3)/SE(3)** exp/log/adjoint (no Sophus/manif dependency; Sophus is the fallback).
- **Config**: one validated POD struct at construction (`validate()` → status, bounds/capacity checks), preallocate once. Per-sensor sub-configs nest under a global config; every knob documented. File parsing (YAML/JSON/ROS-param) is an adapter that *builds* the struct.

---

## 10. Validation

- **Unit tests** per block (median/Weiszfeld, histogram mode + aging, hand-eye solve, xcorr, ESKF update) with analytic expected values.
- **Simulation rig with known ground truth** (the backbone): synthetic trajectory + N synthetic sources with planted extrinsics/scale/offset/noise/outliers → the only place truth is known, hence the only place calibration *correctness* is checkable. Sweep regimes, noise, dropout.
- **Observability self-tests**: assert each DOF converges *only* in its regime (lever-arm frozen under pure-straight sim; converges under turning) → regression-guards the spine.
- **NEES / NIS consistency** (Monte-Carlo): a filter that publishes Σ must prove it is neither over- nor under-confident. **Slice-14 status**: NEES Monte-Carlo harness implemented (6-DOF pose NEES from the published Σ vs sim GT, ensemble over seeds). **Finding**: the published Σ is **grossly pessimistic** — ensemble-mean NEES ≈ 0.13 vs DOF 6 (~46× over-conservative). Root cause: the ESKF inits `P = I₁₂` (≈100× the steady-state error) and the predict-only integrator (no correction step) never shrinks it; the right-error `Ad` propagation further inflates the translation block over distance. Not tunable via `q_scale`/`q_floor` (those only add to P). Fixing it needs a smaller init-P seeded from the first-window/median uncertainty, and/or the Slice-11 correction step. **NIS** is deferred to Slice 11 (no `ICorrection` updates exist yet → no innovation to normalize). The golden-regression net (deterministic recorded-scenario replay + committed numeric snapshot) is in place.
- **Recorded-data golden regression**: real logs replayed; determinism → byte-stable golden output as a regression net.

---

## 11. Default knobs (all config-exposed)

| Knob | Default / form |
|---|---|
| Confidence combine (native ⊕ modeled Σ) | **sum**; {native-only, modeled-only, sum, max, weighted} |
| Histogram aging | {exponential-decay `γ`, sliding-window `K`} |
| Phase-2 histograms | roll = circular S¹, xyz = 3×1-D ℝ³ |
| Commit gate | peak-concentration ≥ τ_commit ∧ votes ≥ N_min; hysteresis τ_drop |
| Weiszfeld | bounded `max_iters` + tol; ε-regularized `1/d`; split-metric weight `λ` |
| Gates | straight (ε_ω, δ_v), turn (θ_ω), excitation (σ²_ω min) |
| Reverse-fold | flip forward sample into consensus/prior hemisphere (dot<0 → negate) |
| Cold-start | {reference-only-until-confident, median-from-start} |
| Time-sync | on/off; pluggable match metric; bounded lag range |
| Per-sensor KF | on/off; CV twist; smoothing-strength (process noise); fixed-lag L |
| Per-source bias states | on/off (default off) |
| Persistence | on/off; period; {calibration-only, full} |
