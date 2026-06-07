# Design Decision Log

Record of the design-grilling session: each decision, the option chosen, the alternatives rejected, and the rationale (including the catches that grilling surfaced). The resulting specification lives in [`DESIGN.md`](./DESIGN.md); this file is the *why*.

Scope chosen up front: grill the **whole system architecture**, greenfield.

---

## D1 — Motion domain → **platform-agnostic core**
- **Chosen**: abstract manifold + plugins from day one.
- **Rejected**: SE(3) ground robot (my rec), SE(2) planar, aerial 6-DOF.
- **Why / resolution**: agnostic + lightweight only reconcile as **small generic core + thin platform plugins** — generality at the edges, never the center. Load-bearing for everything below.

## D2 — Fusion engine → **ESKF on Lie groups**
- **Chosen**: error-state KF, calibration handling decided separately later.
- **Rejected**: IEKF, fixed-lag smoother, full incremental smoother (GTSAM).
- **Why**: bounded compute/memory honors "lightweight"; error-state keeps manifold math clean. Accept weaker calibration observability vs a smoother; mitigate elsewhere.

## D3 — Propagation input → **weighted geometric median of all odometries** (user's design)
- **Chosen**: buffer each source, integrate over a common window to an SE(3) delta, collapse to one consensus via weighted **geometric median in SE(3)**.
- **Corrections made**: (a) "median" requires the **unsquared** distance (Weiszfeld IRLS) — squaring gives the Karcher *mean* and loses robustness; (b) SE(3) has **no bi-invariant metric** → use a **split SO(3)/ℝ³ metric**; (c) needs **≥3 sources** to reject an outlier.
- **Why**: robust consensus motion instead of trusting one source.

## D4 — ESKF correction step → **robust integrator + optional absolute refs**
- **Chosen**: predict on median; **Q from inter-source spread** (adaptive); corrections only from optional absolute-reference plugins; no absolute ref ⇒ honest dead-reckoning.
- **Rejected**: fixed Q (wastes the spread signal); update from per-sensor KFs (double-counts → inconsistent covariance); const-vel predict with odom as updates (abandons median-drives-predict).
- **Why**: a KF without an update is an integrator — name it honestly; the spread gives adaptive Q for free.

## D5 — Calibration regime gating → **symmetric gates** (+ user refinements)
- **Chosen**: Phase 1 gated straight (`‖ω‖<ε ∧ ‖t‖>δ`); Phase 2 gated turning (`‖ω‖>θ`), votes weighted by `‖ω‖`.
- **Catch (load-bearing)**: roll and the **xyz lever-arm are unobservable in straight motion** — the user's original Phase 2 had no turn-gate. Straight-motion samples give *unconstrained* (not wrong) xyz that smear the histogram.
- **User refinements**: reverse-direction handled explicitly; optional reference-sensor cross-check (ice-slide / drone niche); turn magnitude usable as vote weight.

## D6 — Gauge & bootstrap → **base frame + pinned reference sensor**
- **Chosen**: user base frame; one reference sensor pinned (extrinsic = prior) as the gauge; bootstrap all extrinsics from config priors and iterate; **cold-start switch** (reference-only-until-confident *or* median-from-start).
- **Rejected**: float the gauge (underdetermined, never settles); raw-frame median (breaks for non-aligned mounts).
- **Why**: odometry-only ⇒ **gauge freedom** — extrinsics recoverable only relative to a pinned anchor. The median needs frame-aligned deltas (needs extrinsics) while calibration needs the fused output → circular; bootstrap-from-prior breaks it.

## D7 — Source contract → **uniform delta query + native ⊕ modeled Σ** (+ user refinement)
- **Chosen**: `delta(t0,t1) → (SE3, Σ6×6)`; adapter converts twist/increment/absolute and interpolates in-window. Confidence = Σ.
- **User refinement**: combine the provider's **native** Σ with a **modeled** Σ (configurable); missing native → identity.
- **Why**: one "confidence" definition flows to median weight, adaptive Q, and histogram votes.
- **Impl note (Slice 1)**: missing native Σ is treated as **identity** then run through the configured combine rule (not modeled-only); native-present is decided by the explicit `native_confidence` flag, not by inspecting the matrix. Absolute-pose providers' supplied Σ is taken as the per-step increment covariance (no differencing of absolute covariances — not PSD-safe).

## D8 — ESKF update step (revisited with full picture) → **integrator + optional absolute**
- Confirmed D4 once calibration was settled as a separate subsystem.

## D9 — Histogram vote aging → **configurable decay / sliding-K** (user)
- **Chosen**: config between **exponential decay** and **sliding window of K**.
- **Why (forced)**: bootstrap means early votes are computed with *wrong* priors → systematically biased → must wash out. Also tracks re-mount/thermal drift. Pure accumulation would never forget the bias.

## D10 — Phase-2 cost → **hand-eye, split solve** + **pairwise** (user wants both, compared)
- **Chosen**: implement two strategies and compare: (1) per-sensor hand-eye `A·X=X·B` vs fused base — **roll = 1-D rotation residual, xyz = linear-LS** `(R_A−I)t_X = R_X t_B − t_A`; (2) **pairwise** hand-eye between sensors with the reference extrinsic **fixed**.
- **Key structure**: xyz is *linear* given rotation; roll is the only nonlinear DOF (1-D). The LS singularity at `R_A=I` independently re-derives the turn-gate.

## D11 — Direction histogram (pole problem) → **3-channel so(3)/tangent @ per-sensor prior** (user's idea, refined)
- **Discussion path**: rejected naive 2×1-D (yaw,pitch) — poles. Considered: 3-component Cartesian (anisotropic, ±1 boundary, renormalize); Fibonacci equal-area sphere; tangent-plane @ running estimate (recenter fights aging).
- **Chosen**: the so(3) rotation vector of the minimal rotation, histogrammed per channel, basepoint = **each sensor's prior** forward.
- **Insight**: this *is* the tangent-plane method with a **fixed** basepoint → kills the recenter con; data sits near the basepoint (small tilt, no pole/antipode, isotropic). **User correction**: store **all 3 channels** (φ_x≡0 only in the basepoint-aligned frame; noise makes it nonzero), take the so(3) mode, convert to yaw,pitch **at the end**.

## D12 — Output timing → **lagged frontier + predicted tip**
- **Chosen**: fixed-rate tick; fuse causally at `now − delay` (consistent state + Σ); also expose a predict-only extrapolation to `now` (inflated Σ) for control.
- **Rejected**: lagged-only (no real-time tip); OOSM (complex); fixed-lag re-smoothing (drifts toward rejected smoother).

## D13 — Output contract → **transport-agnostic library + adapters**
- **Chosen**: pure C++ core, no middleware dep; callback + poll; rich result struct (frontier state, tip, per-sensor calib snapshot, diagnostics). ROS/DDS/zmq adapters outside core. Pose in an odom frame anchored at init (drifts); extrinsics in base frame.
- **Rejected**: ROS-first (couples core); pose-only (hides calib/diagnostics); streaming IPC (transport dep).

- **Impl note (Slice 13, adapters subset)**: a relaxed-edge `adapters/` tree (`ofc_adapters`, links the core PUBLIC API only, no new deps; the gate builds+tests it via `OFC_BUILD_ADAPTERS=ON` in `dev.ps1`). Shipped: a **file-persistence** double-buffer ping-pong (production form of the Slice-12 core serialize/deserialize — overwrite target chosen by VALIDITY, not raw seq, so a torn higher-seq file can't clobber the last-good; a review CRITICAL was caught + fixed here), a **threading wrapper** (`ThreadedEstimator` — one worker pumps `step()`, mutex-guarded `Result` snapshot), and a dep-free **config loader** (D19). The **ROS node** is a compile-guarded header sketch only (no ROS on the dev box) and a true `fsync` is deferred — both → **Slice 13b**.

## D14 — Threading → **caller-pumped single-thread** (+ user: **C++14 / AUTOSAR**)
- **Chosen**: core is a state machine the consumer pumps via `step()`; fusion + bounded calibration slice per step; lock-free, deterministic, trivially testable. Threading = external adapter.
- **Constraint added by user**: **C++14** (AUTOSAR). Ripples: no `std::optional`/`variant`/`string_view`, likely no-exceptions, no runtime heap, bounded WCET, `double`/no-fast-math. Mostly *reinforces* prior choices.

## D15 — Memory/error policy → **strict core, relaxed edges**
- **Chosen**: core = no heap post-init, no exceptions, hand-rolled `Optional`/`Expected`, bounded loops, fixed-capacity from config; adapters + tests stay relaxed (off-target).
- **Rejected**: strict-everywhere (painful test ergonomics); pragmatic-harden-later (retrofit touches every file); dual-build (#ifdef drift).

## D16 — Time-sync → **xcorr of ‖ω‖ → offset histogram** (+ user: pluggable metric)
- **Chosen**: cross-correlate the **extrinsic-invariant** ‖ω‖(t) signal (so temporal decouples from spatial and runs first), parabolic sub-sample, excitation-gated, vote into a per-source offset histogram (reuse the robustness primitive). Constant offset + slow track.
- **User note**: make the **match metric pluggable** (L1, L2, ratio, NCC…) — worth sweeping.
- **Rejected**: continuous SSD (no histogram robustness); joint-with-extrinsics (breaks the decoupling); ω-vector correlation (needs extrinsics first).

## D17 — Weight refinement → **variance-EMA reliability, bias → calibration**
- **Chosen**: `w ∝ 1/σ̂²` from **zero-mean** residual scatter (slow, floored/capped); the **systematic** residual routed to the calibrator, not the weight; Weiszfeld `1/d` handles per-step outliers. Final `w = prior × reliability × Σ-confidence`.
- **Why**: downweighting a *systematically* disagreeing source masks a miscalibration the calibrator should fix, and risks majority lock-in suppressing a correct minority. Split bias from variance.
- **Impl note (Slice 9)**: per-source EMA mean (`resid_mean` = bias) + EMA variance around that mean (`resid_var` = zero-mean scatter) of the residual-to-consensus distance, at `reliability_ema_alpha`. `reliability = clamp(ref_var / resid_var, reliability_floor, reliability_cap)`, `ref_var` = robust median of warmed-up participants' variances (equal-noise → ≈1; noisier → <1; cleaner → >1, floored/capped). `w = clamp(prior × reliability × Σ-confidence, weight_floor, weight_cap)`. **Causal** (weight reads the prior-step reliability; the EMA updates after the median, so this step's residual never weights this step's median) and 1.0 until a per-source warmup. The bias is **not** folded into the weight — it is left for the existing Phase-1/2 calibrator observe path to absorb (the routing is architectural, not a new feedback wire); `SourceHealth.bias` is diagnostic-only and is an unsigned residual *magnitude* (mean `split_distance`), not a signed per-DOF offset. `sigma_confidence()`'s unit-mixing (D21) is left intact — reliability is an additive quality factor, not a rewrite of the Σ-confidence; the unit-separation caveat stays open.

## D18 — Per-sensor KF → **CV twist ESKF (ℝ⁶)** (+ user: fixed-lag smoother)
- **Chosen**: per-sensor (config on/off) constant-velocity error-state KF on body twist, random-walk accel; smoothed twist + refined Σ.
- **User insight (load-bearing)**: calibration tolerates latency → run it as a **fixed-lag RTS smoother** (two-sided) at a **deeper frontier** (`now − delay − L`). Two-sided smoothing is near-zero-phase → dissolves the over-smoothing-lag worry. Fusion stays causal; calibration gets nicer odometry.
- **Impl note (Slice 10)**: `ofc::TwistSmoother` (`smoother.hpp`/`.cpp`) — augmented CV model `x = [twist(6); accel(6)]`, `F = [[I, dt·I],[0,I]]`, white-accel `Q = q·[[dt³/3, dt²/2],[dt²/2, dt]]⊗I` (`q = kf_process_noise`), measurement `H=[I,0]`, `R = r·I` (`r` fixed 1.0; only `q/r` shapes smoothing). Forward KF + **backward RTS pass** over a depth-`L+1` window emits the lag-`L` smoothed twist (true fixed-lag two-sided, verified **zero-phase**: emitted err <0.1× the causal-lag floor; variance ↓~0.3×). Wired into the calibration feed at the deeper frontier (gated by `per_sensor_kf`; OFF ⇒ byte-identical). Calibration peaks ~4× sharper (extrinsic_confidence 0.038→0.156) with no estimate bias. **Open**: a single shared smoother uses the MAX `q` over enabled sources (not per-source `q`/`r`); the refined RTS covariance is computed/exposed but NOT yet wired into the calibrator vote weight (the deeper path still feeds the raw Σ-confidence) — D18's "refined Σ" promise is half-wired; ~7.5 MB fixed footprint (compile-time caps); the lag is a step count via the nominal tick rate (off-cadence ⇒ effective time-lag ≠ `calib_lag_s`).

## D19 — Config system → **validated POD struct + adapter parsing**
- **Chosen**: core takes one validated config struct, preallocates once at init; YAML/JSON/ROS loaders are adapters that *build* the struct.
- **Rejected**: built-in file config (parser + heap + IO in core); compile-time config (recompile to retune); hybrid compile-sizes (capacities need recompile).
- **Requirement**: a **config-reference doc** documenting every knob (user-stressed).

## D20 — Scale calibration → **straight-regime magnitude-ratio histogram**
- **Chosen**: per-source translation scale `s_i` = magnitude ratio vs the pinned reference, recovered in the **straight** regime (pure translation → zero lever-arm confound), voted into a scale histogram, fed back **pre-median**. Per-source on/off; metric sources fixed at 1.
- **Rejected**: joint-with-lever-arm in Phase 2 (nonlinear, confounds); none (wrong for wheel-scale drift / mono-VO); single global scale (per-source by nature).
- **Why**: real for automotive (tire pressure/wear drifts wheel-odom scale); slots into existing regime + machinery.

## D21 — ESKF state → **pose + twist, SO(3)×ℝ³ error, calib external**
- **Chosen**: state = `SE(3)` pose + `ℝ⁶` twist; decoupled `SO(3)×ℝ³` pose error (matches split-metric median) + `ℝ⁶` twist; dense 12×12; calibration params **not** in state.
- **Rejected**: pose-only coupled (cruder twist Σ/tip); IEKF-coupled (declined path); augmented-with-calib (contradicts separate-histogram calibration, bloats filter).
- **Impl note (Slice 2)**: error-state ordering is `[trans; rot]` (pose 0–5, twist 6–11). Predict on a delta uses right-error transport `F = blkdiag(Ad(delta⁻¹), 0)` (twist re-read each step as `log(delta)/dt`, so its prior error doesn't carry), `P ← F P Fᵀ + blkdiag(Q, Q/dt²)`. Adaptive `Q = q_scale·spread² · I₆ + diag(q_floor)`. Per source, the delta is **scale-corrected then frame-aligned**: `B_corr = {B.R, B.t / prior_scale}`, then `A = X ∘ B_corr ∘ X⁻¹` with `X = prior_extrinsic` (sensor→base). Fusion applying `prior_scale` is no longer deferred (pre-median scale *calibration* still is). Accepted Slice-2 approximations: zero pose↔twist cross-covariance; `sigma_confidence` unit-mixing (until Slice 9 reliability EMA).
- **Time-offset sign convention (sim/oracle, Slice-5 contract)**: **positive `prior_time_offset_s` shifts the sampled window later** — a source's reported motion for `[t0,t1]` corresponds to true motion over `[t0+off, t1+off]` (source clock ahead of base). Slice-5 time-sync must invert this sign.
- **Impl note (Slice 14, covariance tangent)**: the "decoupled `SO(3)×ℝ³`" wording describes the error *ordering* (`[trans; rot]`) and the **median's** split metric — NOT the covariance. The ESKF propagates `P` with the **full SE(3) adjoint** `F = Ad(delta⁻¹)` (the `[t]ₓR` coupling block), so the dense 12×12 covariance lives in the **full coupled SE(3) tangent**; NEES/consistency must use the full-SE(3) `se3::log` error to match it (verified in Slice-14 review). **Init-P fix (`70c7d38`)**: the filter now inits `P = 0` at the gauge-anchored first fuse and lets the first `predict()` establish the one-window covariance `blkdiag(q_pose, q_pose/dt²)` (was `I₁₂`), cutting the pessimism from ~46× (NEES ≈ 0.13) to ~17× (NEES ≈ 0.35). The residual is the predict-only `Ad(delta⁻¹)` translation-block inflation over accumulated motion (no correction to shrink it; rotation block well-calibrated) — strict no-ref consistency needs a distance-aware covariance model / no-ref correction. See DESIGN §10.

## D22 — Absolute refs → **generic measurement plugin + optional per-source bias states**
- **Chosen**: one interface `evaluate(state) → (residual, H, R)` + stamp, Mahalanobis-gated; position/pose/orientation agnostic; **loop-closure excluded** (needs past states ⇒ smoother). **Plus optional per-source bias states** (config per-source).
- **Why the bias states** (user's GPS/IMU experience): the classic "like a charm" GPS/INS drift removal comes from **online IMU bias estimation via GPS cross-covariance** — the plain plugin (no bias states) corrects fused pose but *can't recalibrate a raw IMU*. Bias is a fast nuisance state observable only via an absolute ref at the fusion level → belongs in the filter (genuinely different from the slow geometric extrinsics in histograms). Off by default; enable for a raw-IMU source.
- **Impl note (Slice 11, correction path)**: the **plugin + Mahalanobis-gated ESKF update half** of D22 shipped (`Eskf::update`: `S=HPHᵀ+R`, NIS gate, gain via one LDLT, Joseph covariance, right-error full-SE(3) injection; `Estimator::add_correction` wired into `step()` after predict; `Result::CorrectionDiag`; `validate()` `mahalanobis_chi2>0`). Sim `SyntheticAbsoluteRef` position fix → drift removed (0.58→0.20 m); gross outlier rejected (NIS ~3e5); NIS now computable (~2.4 vs DOF 3). The gate is a **single scalar regardless of measurement DOF** (fine for the dim=3 fixes; needs a per-`n` χ² quantile for 6-DOF/mixed plugins). Bias states + GPS adapter were deferred to Slice 11b (see below).
- **Impl note (Slice 11b Option A)**: chosen + shipped (`48ece29` + `f47b29a`). The architecture problem: D22's GPS/INS bias-estimation assumes the biased sensor *drives the predict*, but here the **median of N sources** drives it — so the pose↔bias cross-covariance doesn't form classically. **Option A** handles the regime where it IS exact: a SINGLE `bias_states` source driving the predict alone augments the ESKF to **18-DOF** `[pose;twist;bias(6 body-twist rate)]`; predict de-biases `Δ_db = Δ∘exp(-b·dt)` and builds the cross-cov via `J_pb = -dt·I₆` (the exact leading term — reviewer-verified the full Jacobian is `-dt·Ad(exp(b·dt))`, omitted part O(dt²)); the absolute-ref update injects `δb` through the cross-cov (Joseph, 18-DOF). New knob `bias_process_noise`; single-bias guard (`>1 → InvalidConfig`); `predict_aug_frozen` keeps the filter consistent out-of-regime; default-OFF byte-identical. Sim: planted bias recovered, drift 16 m → 0.06 m, no-ref observability self-test. **Still DEFERRED (11b residual)**: **Option B** (median-coupled multi-source bias via the per-source median weight in the augmented `F` — general N-source case, higher risk; design in `ISSUES.md` Slice-11b DESIGN NOTE), the per-`n` χ² Mahalanobis gate, and the concrete GPS adapter.

## D23 — Startup → **lifecycle state machine, degrade-don't-block** (+ user: persistence)
- **Chosen**: `INIT → WARMUP → DEGRADED → NOMINAL`; emit best-available (reference-only) output ASAP, auto-upgrade as calibration converges; readiness + confidence exposed; graceful downgrade on fault; gauge anchored at first fused tick.
- **Impl note (Slice 3)**: the ladder is a single per-step recompute from the participating-source count `n` (no high-water-mark — downgrade falls out for free). `INIT` set at `init()` (observable before the first step); pre-fuse (`t1≤q0` or `n==0`) → `WARMUP`; first successful fuse latches `ever_fused` and gives `NOMINAL` iff `n ≥ min_sources_warn` (new `Config` knob, default 3; `validate()` enforces lower-bound `≥1` only — `> max_sources` is a legitimate "never NOMINAL" config) else `DEGRADED`; once fused, a non-fusing step is `DEGRADED`+`NotReady` and keeps the last good frontier. A coarse `Result::readiness` scalar (Init 0.0 / Warmup 0.25 / Degraded 0.6 / Nominal 1.0) is published alongside `phase` via a single `publish_phase()` sink so the two can never drift. NOMINAL is source-count-driven; calibration-convergence coupling is only indirect (via `ReferenceOnly` cold-start participation) — a deliberate Slice-3 scope boundary (per-DOF convergence is already exposed in `CalibSnapshot`).
- **User addition — warm-restart persistence**: config on/off, periodic; persist **calibration state**; **double-buffer ping-pong** (two files + sequence counter, version-tagged + checksummed) so a crash mid-write keeps the last good state; invalidate on setup change via manual delete **+ auto config-hash guard**; serialize/deserialize in core, file IO in adapter.
- **Impl note (Slice 12)**: core `Estimator::serialize`/`deserialize` into a caller fixed buffer (`persistence.hpp`: "OFCP" magic + `format_version` 1 + FNV-1a-64 `config_hash` + payload + FNV-1a-32 checksum; explicit-little-endian, padding-free, IEEE-754-bit-exact doubles; new Status `CorruptData`/`VersionMismatch`). Validation order magic→version→exact-length→checksum→config-hash (no OOB read before the checksum) + a restored-extrinsic orthonormality guard. `config_hash` covers rig shape + per-sensor priors/flags + all histogram configs + commit/gate thresholds + `phase2_strategy`/`calib_lag_s` (runtime knobs tick_rate/fusion_delay/output are deliberately excluded). **Histogram bins are NOT persisted** (deviates from the literal §7 "histogram bins"): deserialize restores committed values + flags + reliability + lifecycle and **re-anchors** the calibrators; hysteresis (+ an `offset_restored` latch) holds restored-committed DOFs through the empty-histogram refill → near-NOMINAL (warm ~0.001 m vs cold ~0.256 m). The **double-buffer ping-pong file IO is relaxed-edge** (test-side here; the production adapter is Slice 13).

## D24 — Validation → **layered, sim-GT backbone**
- **Chosen**: unit (per block) + **sim rig with known ground truth** (only place calibration *correctness* is checkable) + **observability self-tests** (each DOF converges only in its regime → guards the spine) + **NEES/NIS** consistency (a Σ-publishing filter must prove calibration of its own covariance) + recorded-data **golden regression** (determinism → byte-stable).
- **Impl note (Slice 14)**: NEES Monte-Carlo (6-DOF pose, full-SE(3) tangent, ensemble over seeds) + golden regression (deterministic noise-free committed snapshot + byte-identical replay incl. the calib snapshot) implemented in `tests/test_validation.cpp`. **NIS deferred to Slice 11** (no `ICorrection` updates → no innovation yet). **Finding**: the published Σ is **pessimistic ~46×** (ensemble NEES ≈ 0.13 ≪ DOF 6) — init `P=I₁₂` + predict-only Ad inflation; needs a smaller init-P and/or a Slice-11 correction step (not Q-tunable). The CONFIG "tuned"-placeholder sweep remains open → Slice 14 stays `[~]`.

---

## Leaf defaults (locked)
Confidence-combine = **sum** (configurable); Phase-2 histograms roll = circular S¹, xyz = 3×1-D; commit on peak-concentration ≥ τ ∧ votes ≥ N with hysteresis; Weiszfeld bounded iters + ε-regularized `1/d` + split-metric weight `λ`; gate thresholds documented; reverse-fold into the consensus/prior hemisphere before voting.

## The spine, in one sentence
Each calibration DOF is observable in exactly one motion regime, and the math hands you the gate for free (the lever-arm least-squares goes singular precisely when there is no rotation). Every gate, vote, and weight derives from that.
