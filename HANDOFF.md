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

- **Gate green: ctest 2/2 — `unit` 214 cases + `adapters`.** Pushed through **`4ca0c6c`** (HEAD). Since the last handoff, ALL the real-data validation + a new **Slice 15 series** landed (see the "Real-data + Slice 15" block below). (The gate also builds the relaxed-edge adapters — `dev.ps1` sets `OFC_BUILD_ADAPTERS=ON`, so CTest runs both `unit` and `adapters`. **Use `dev.ps1 -Task clean` before A/B covariance/median measurements — stale ninja artifacts gave misleading NEES numbers.**)
- **USER PRECISION TARGETS (north-star, 2026-06-09):** rotation error **< 0.1°** (1.745e-3 rad), translation error **< 1 cm** — **rotation precision takes priority** over translation when they trade off (heading error is what makes position diverge). Judge drift/calib results against these, not just "better than before." NOT yet met (closest: KAIST yaw calib ~0.1–0.23°, EuRoC lever ~4 cm).

### Real-data validation + Slice 15 series (this session — all pushed)

Real-data path is the `ofc_replay` CLI on a manifest → drift / 6-DOF NEES / NIS + (NEW) a local metric + final per-source calib. Datasets converted by `tools/{kitti,kaist,euroc}_to_csv.py` (read source archives IN MEMORY, pull only the needed CSVs, never unzip; run dirs `*_run/` are gitignored). Data at `C:\workspace\data\{KITTI,KAIST,EuRoC}`.

| What | Outcome |
|---|---|
| **KITTI** (28 drives, smoke test) | pipeline works on real data; single-source → covariance overconfident (q_floor lever); `ConfigLoader` exposes `q_scale`/`q_floor`/`adaptive_q` |
| **KAIST** urban07/12/17 (3-source + GPS + GT) | median-fix robustness CONFIRMED (×3 outlier rejected); urban07/17 bounded; **urban12 was the hard case** |
| **Local metric** (`7c71b63`) | GT-anchored fixed-window relative-pose error (`local_batch_len`): length-FAIR drift (global tail grows with run length; local windows don't). Eval-only, no filter reset. `diagnose_window.py` localizes the worst window. |
| **Slice 15 — robust GPS update** (`44e0711`) | Huber gain down-weight (`correction_robust_kappa`, 0=off). Opt-in primitive; did NOT fix urban12 (gate pre-empts it); kept. Realistic GPS R (`cov_floor_m2`) is the consistency keeper. |
| **Slice 15b — bounded heading injection (lever C4)** (`ba895b6`) | **FIXED urban12**: tail 4214→**1.99 m**. `correction_rot_suppress_kappa` scales ONLY the rotation gain rows of a position fix with a large residual (kills the t=1929 s 68° heading kick via the trans-rot cross-cov; keeps the translation pull). 0=off. Recommended urban cfg: `cov_floor_m2=25`+`correction_rot_suppress_kappa=0.8`. |
| **Online calibration validated on real data** (`a530deb`) | `tools/inject_calib.py` injects a KNOWN extrinsic/scale/time-offset into a source → calibrator recovers it on real motion. KAIST (ground): yaw 0.103° + scale recovered. EuRoC (drone, 3D): **full xyz lever recovered** (KAIST couldn't — planar). Gap: scale/rotation-extrinsic NOT recovered on drone (excitation-regime-sensitive). |
| **Slice 15c — calib-commit on real data** (`4ca0c6c`) | The accurate yaw/scale calib NEVER COMMITTED on real data (conf~1, uncommitted): `commit_min_votes`(200) is a vote-MASS threshold but Combo weight (ω-floored 1e-3 in straight) keeps mass ≪200. **FAILED approach** (reverted): global obs-COUNT commit gate broke the sim's tuned feedback (cov-cal NEES 4.99→0.8 — the sim relies on commit+publish firing). **WORKING fix = config**: exposed `vote_weight` as a manifest key; `vote_weight=one` makes mass==count → commit reachable → yaw/scale/time COMMIT + feed back. Sim untouched (Config{} default=Combo). |

**Slice-15 design docs**: `SLICE15_ROBUST_GPS_UPDATE.md`, `SLICE15B_BOUNDED_HEADING_INJECTION.md` (full sweeps + acceptance). Memory `real-dataset-testing` has the blow-by-blow.
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
| Slice 12 | warm-restart persistence (D23): core `serialize`/`deserialize` into fixed buffers (`persistence.hpp`, "OFCP" v1, FNV-1a config-hash + checksum + orthonormality guard, `CorruptData`/`VersionMismatch`); re-anchor-and-refill restore resumes near-NOMINAL (warm ~0.001 m vs cold ~0.256 m); crash-mid-write/config/checksum/version all reject. Bins NOT persisted; file double-buffer is relaxed-edge (production adapter → 13) |
| Slice 13 (subset) | relaxed-edge `adapters/` (`ofc_adapters`, core PUBLIC API only, no new deps): file-persistence double-buffer (validity-based overwrite target — torn higher-seq can't clobber last-good), threading wrapper (mutex-guarded snapshot), dep-free config loader (subset → Config), **GPS correction adapter** (`GpsCorrection`), and **real-data CSV ingestion** (`ddcf436`): `CsvSource`(ISource, 3 forms) + `CsvGtTrack` + `ReplayHarness` (drift/NEES/NIS vs GT) + the `ofc_replay` CLI. ROS node + true fsync → 13b |
| Slice 11b (Option A) | per-source bias states (augmented 18-DOF ESKF): a SINGLE `bias_states` source driving the predict alone → predict de-biases (`Δ∘exp(-b·dt)`) + builds the pose↔bias cross-cov (`J_pb=-dt·I₆`); absolute-ref update removes the bias (sim: planted recovered, drift 16 m → 0.06 m); no-ref observability self-test; multi-bias guard; `predict_aug_frozen` out-of-regime; default-OFF byte-identical. `bias_process_noise` knob; `CalibSnapshot.bias`/`bias_observable`. Option B + per-n gate + GPS adapter → 11b residual |

The **calibration spine (5–8) is complete** — calibration closes back into fusion and bootstraps from arbitrary priors. **All numbered roadmap slices 0–14 are now addressed** (13 + 14 are partial — see below).

- **Remaining work** (any order):
  - Slice 11b residual — **Option A DONE** (single-driving-source augmented bias filter) + **per-`n` χ² gate DONE** (`e8491dd`) + **GPS adapter DONE** (`4bce8d1`). **Option B** (median-coupled multi-source bias) was a NO-GO under the pinning median but is **reconsiderable** post-fix (ISSUES Slice-11b UPDATE).
  - Slice 13b — real ROS node + recorded-bag round-trip; replace the persistence adapter's `flush`/`close` with a real `fsync` (durability). Deferred from Slice 13 (no ROS on the dev box).
  - Slice 14 — `[~]` partial: NEES + golden + **NIS** DONE; **init-P covariance fix DONE** (`70c7d38`, NEES ~0.13→~0.35); **median-variance Q reduction DONE** (`730bcfa`, approach A — NEES ~0.35→~1.0, see the covariance bullet in §6). Remaining: the CONFIG `q_scale`-coefficient sweep, and — for strict no-ref NEES consistency — the `Ad` distance-*shape* model and/or an NHC no-ref correction (approach B).

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
- **Median solver (`median.cpp`, `1142e41` — was a CORE BUG, now fixed)**: the n≥3 weighted geometric median uses Weiszfeld IRLS under the split metric, **initialized OFF-vertex at the weighted mean** + a **Vardi-Zhang coincident-vertex guard** (`d≤eps` self-term excluded). Do NOT revert it to init-at-the-highest-weight-vertex: that made iter-0's `d_start=0` give a `w/eps≈1e9` self-weight that pinned the result on that vertex in one iteration → fusion returned the highest-weight source verbatim (no blending, a high-weight outlier returned with zero rejection). It is a true interior robust median now; the high-weight-outlier guard in `test_median.cpp` pins this. `q_scale` (covariance) was recalibrated for this true median — see the covariance bullet in §6.
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
- `validate()` now range-checks `scale_hist` (must strictly contain 1.0, `3bd91e2`); the broader `TODO: per-sensor + histogram range checks` for the other nested `HistogramConfig`s remains (they're still validated at `Histogram1D::configure()` time, not top-level).
- Several thresholds are tuned placeholders (CONFIG marks them) pending the Slice-14 sweep. **`q_scale` is now calibrated** (`983fa65`, 1.0→0.5); the rest (`excitation_min_var`, `kf_process_noise`, `match_metric`, straight/turn gates) are a separate non-covariance pass (the observability self-tests pin them functionally).
- **Slice-3 lifecycle scope**: NOMINAL is source-count-driven (`n ≥ min_sources_warn`), not directly calibration-convergence-gated; under `ReferenceOnly` cold-start the DEGRADED→NOMINAL upgrade tracks convergence only *indirectly* (a source joins the median once its extrinsic commits). `min_sources_warn` is validated lower-bound only (`≥1`); a value above `max_sources` is legitimate (NOMINAL never reached). If a future slice wants readiness to encode calibration convergence directly, revisit the ladder.
- **Slice-9 weight scope**: reliability is the variance-EMA quality factor; `sigma_confidence()`'s D21 unit-mixing (mean of trans m² + rot rad² in one scalar) is **left intact** — reliability was added multiplicatively, not as a unit-separation rewrite, so that caveat stays open. `SourceHealth.bias` is an unsigned residual *magnitude* (mean `split_distance`), not a signed per-DOF offset — it cannot itself distinguish direction; the weight uses `resid_var` (scatter), not `bias`.
- **Covariance — NEAR-CONSISTENT now (`1142e41`); the big saga is RESOLVED**: ensemble-mean 6-DOF pose NEES ≈ **4.82** vs DOF **6** (worst-case ~4.85 across {nees_traj, mixed, turning, straight} × {1×,2× noise}; never overconfident, the gap a deliberate safety margin). The road here: (1) ESKF originally seeded `P = I₁₂` → NEES ≈ 0.13 (~46× pessimistic); **init-P fix** (`70c7d38`, seed P=0 + first predict establishes `blkdiag(q_pose, q_pose/dt²)`) → ≈ 0.35. (2) Steps 3–4 added a `/n_eff` divisor (`730bcfa`) + `q_scale=0.5` (`983fa65`) → ~2.07 — but these were **calibrating around a core bug**: the geometric-median solver was **PINNING** (it returned the highest-weight source's delta verbatim, see §5 "median" gotcha + DECISIONS D3) so it averaged NOTHING and the `/n_eff` "median variance reduction" rationale was FALSE. (3) **Median-pinning FIX** (`1142e41`, off-vertex init + Vardi-Zhang guard) made fusion a true interior blend → spread correctly sized → `/n_eff` then OVER-divided Q (overconfident, NEES ~20), so `/n_eff` was **REMOVED** and `q_scale` recalibrated to **0.7** → the near-consistent 4.82. With the true median, smaller `q_scale` → larger NEES (opposite of the pinning sign). Strict per-trajectory DOF=6 would still want an `Ad(δ⁻¹)` distance-*shape* model and/or an NHC no-ref correction (deferred). NIS (with an absolute ref) ≈ **2.7** vs DOF 3. Guards: NEES band [4.0,5.6] + `test_cov_calibration` (never-overconfident <5.5 + near-consistent across trajectories) + `CHECK_FALSE(truly_consistent)`.
- **Slice-11 correction-gate — per-`n` now (RESOLVED, `e8491dd`)**: the Mahalanobis gate WAS a single scalar `mahalanobis_chi2` regardless of measurement DOF `n` (~97% quantile for n=3, ~80% for n=6). Now `Eskf::chi2_gate(base, n)` scales the n=3-tuned base by the χ²-quantile ratio `q[n]/q[3]` (const `kChi2Q95`, 0.95 confidence) so every DOF `n∈1..6` gates at the same confidence; the estimator passes `chi2_gate(cfg.mahalanobis_chi2, m.dim)` to `update()`/`update_aug()` (signatures unchanged — they still take a raw per-n threshold). n=3 returns base unchanged, so the dim=3 position fixes shipping now are behaviorally IDENTICAL; the per-`n` scaling becomes load-bearing once a 6-DOF/mixed plugin lands. The drift-removal test still uses a loose `chi2=100` — a **test artifact** driven by the covariance pessimism above (legit drift residuals exceed `chi2=9` because P never shrinks on predict-only stretches), NOT a production value (and unrelated to the per-`n` gate).
- **Covariance tangent (doc-vs-code reconciled)**: the ESKF covariance is **full coupled SE(3)** (propagated by the full `Ad(delta⁻¹)`), NOT block-diagonal "decoupled SO(3)×ℝ³" — that phrase only describes the error ordering `[trans;rot]` + the median's split metric. NEES uses the full `se3::log` to match. (DESIGN §5, D21 now corrected.)
- **Slice-10 smoother scope**: a single shared `TwistSmoother` uses the MAX `kf_process_noise` over enabled sources (not per-source `q`/`r`; `r` fixed 1.0). The refined RTS covariance is computed + exposed but NOT wired into the calibrator vote weight (the deeper path still feeds the raw Σ-confidence) — D18's "refined Σ" is half-wired. Fixed ~7.5 MB footprint (compile-time caps 32×65, paid even when off). The lag `L` is a step count via the nominal tick rate → off-cadence pumping gives an effective time-lag ≠ `calib_lag_s` (rings stay step-aligned regardless). Per-source dropout degrades to a variance loss (no bias) after the push-seq alignment fix.
- **Slice-12 persistence scope**: histogram BINS are not persisted — restore re-anchors the calibrators to the committed values + holds them via hysteresis (an `offset_restored` latch does the same for the time-offset DOF), so a restored DOF's `committed` flag + value are correct immediately but its confidence reads low until the histogram re-fills. `config_hash` covers rig shape + per-sensor priors/flags + histogram configs + commit/gate thresholds + `phase2_strategy`/`calib_lag_s` (runtime knobs excluded).
- **Slice-13 adapters scope**: relaxed-edge only (std/heap/exceptions/threads/file IO). Adapters now: file-persistence double-buffer, threading wrapper, config loader, and the **GPS correction adapter** (`ofc_adapters::GpsCorrection`, `4bce8d1`). The ROS node is a compile-guarded header sketch (not built — no ROS on the dev box) → **13b**. The file-persistence adapter uses `flush`/`close`, NOT `fsync` — an OS-page-cache power-loss window remains even with the crash-safe target selection (durable `fsync` → **13b**). The config loader covers a documented knob SUBSET (not every `Config` field). Adapters are off by default (`OFC_BUILD_ADAPTERS=OFF`); the gate turns them on.
- **Slice-11b bias scope**: **Option A done** — augmented 18-DOF filter for a SINGLE source that drives the predict ALONE (1-source / `ReferenceOnly` dead-reckon). Out-of-regime (a 2nd source joins the median) → `predict_aug_frozen` holds the learned bias + keeps `cov18_` consistent, but de-biases the *consensus* (an approximation; exact per-contribution de-bias is Option-B territory). **Per-`n` χ² gate done** (`e8491dd`); **GPS adapter done** (`4bce8d1`). Deferred: Option B only (median-coupled multi-source).
- **GPS adapter ↔ `q_floor` coupling (surfaced `4bce8d1`)**: the GPS Kalman gain depends on `q_floor`. With two near-identical sources the adaptive-q spread is ~0 so `q_pose` collapses to `q_floor`; a too-small `q_floor` keeps predict-only `P` from re-inflating between fixes → gain ~0 → fixes barely pull (the end-to-end test raises translation `q_floor` to `1e-3`). Documented in CONFIG §3; the right default is a Slice-14-sweep item. **USER decision open**: revisit the default `q_floor` / call this out more loudly?
- Not yet built: 11b Option B (now RECONSIDERABLE — the median fix makes the per-source `ωᵢ` real; see ISSUES Slice-11b UPDATE), ROS node + durable fsync (13b). Partial: validation (14) — NEES/NIS/golden + init-P fix (`70c7d38`) + the **median-pinning fix + `q_scale=0.7` recalibration** (`1142e41`) DONE → NEES near-consistent ~4.82; `/n_eff` (A) was REMOVED. Remaining = the non-covariance placeholders (excitation_min_var, kf_process_noise, match_metric, calib gates); the `Ad` distance-shape model / NHC no-ref correction for strict DOF=6. The two latent calibration bugs the median fix exposed are RESOLVED: 1(a) scale_hist boundary commit FIXED (`3bd91e2`); 1(b) partial scale recovery is NOT a bug (sparse-vote/noise artifact — documented known-limitation, ISSUES Slice 14).

---

## 7. Resume in one line

Slices 0–14 + the median-pinning fix + the real-data validation + the Slice-15 series are all landed and pushed (`4ca0c6c`). The system fuses + self-calibrates + survives real urban/drone data; urban12 divergence is FIXED. The active frontier is now the **USER PRECISION GOAL** (rotation < 0.1°, rotation-priority — see §3). Open work, roughly priority-ordered:

- **🎯 Committed-yaw under 0.1°** (rotation goal). `vote_weight=one` now makes yaw COMMIT + feed back on real data (`4ca0c6c`), but the *committed* estimate is ~0.23° (vs the 0.10° uncommitted Combo estimate) — committing engages the feedback (the win) but loosened the value. Tighten it: commit-time re-anchor dynamics / vote weighting that keeps both commit AND precision / more excitation.
- **scale + rotation-extrinsic calib on non-ground (3D) motion** — EuRoC exposed these don't recover on a drone (conf 0): scale needs straight/low-ω windows (a drone has few); the rotation-extrinsic estimator is yaw-dominant-ground-tuned. A turn-regime yaw/pitch path would generalize calibration off the ground plane.
- **urban12 residual** — the mid-drive transient (max ~300–400 m, RECOVERS) is the upstream **522 s GPS-coast heading-drift** (DR/FOG-IMU heading-hold over long GPS gaps). Rotation-relevant.
- **q_scale ↔ calib-engagement coupling** — the cov-cal `q_scale=0.7` was tuned with calib effectively frozen; once calib commits + corrects (smaller error), it may want re-tuning. Watch when calib engagement changes the accuracy regime.
- **predict-side `q_floor` recalibration on real data** — D (GPS R↑) done; predict side still sim-tuned. Motivates the deferred distance/motion-aware predict-Q (`Ad` shape).
- **Carve-outs**: 11b Option B (median-coupled multi-source bias — reconsiderable); non-cov CONFIG placeholders (`excitation_min_var`, `kf_process_noise`, `match_metric`, calib gates); 13b ROS node + durable `fsync` (env-blocked).

Pick one; for a core change, design-first (the Slice-15 series shows the pattern: investigate → design doc → implement → validate full-drive via the local metric, and **never weaken the trust apparatus** — verify cov-cal NEES stays ~4.82 + sim green). Verify with §2.
