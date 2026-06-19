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

- **Gate green: ctest 2/2 — `unit` ~343 cases + `adapters` 46/772** (numbers grow per slice). **HEAD = `baf6893`** (2026-06-19). Persistence format is now **v4** (Slice 20b; pre-v4 blobs cold-start). **The SESSION DELTA block below (2026-06-15..19) is the authoritative current state** — the line items in this §3 above/below it are historical. Since the original handoff: Slices 16–20b + the radar/camera real-data threads + the covariance sweeps all landed. (The gate also builds the relaxed-edge adapters — `dev.ps1` sets `OFC_BUILD_ADAPTERS=ON`. **Use `dev.ps1 -Task clean` before A/B covariance/median measurements — stale ninja artifacts gave misleading NEES numbers.**)
- **USER PRECISION TARGETS (north-star, 2026-06-09):** rotation error **< 0.1°** (1.745e-3 rad), translation error **< 1 cm** — **rotation precision takes priority**. **CALIBRATION DOFs: MET on both ground and 3D motion** (Slices 16/17/17b): committed rotation extrinsic 9.1e-6° (EuRoC drone), lever sub-mm with unknown scale, scale err 4e-4, time-offset µs-level. Calib recipe on real-data manifests: `subbin_centroid=true + vote_weight=one + rot3d_enabled=true + joint_lever_scale=true`. The remaining precision gap is fused TRAJECTORY drift on real GPS-corrected data (metres-scale; heading now at the best-source floor — see Slice 19).

### SESSION DELTA 2026-06-15..19 (all pushed; this supersedes the older "current state" numbers above)

**Calibration / fusion slices landed:**
- **Slice 19b** (`51b86bf`) — per-channel scatter reliability on the split path (`rel_rot`/`rel_trans`; layer b).
- **Slice 19c** (`770f6b2`+`e7bba55`) — GPS-course heading-drift MONITOR (auto-discovers the heading-grade source, no `rot_weight_prior`; layer c). **Slice 19d** (`4a77219`) — monitor WARM-START from `rot_weight_prior` (erases discovery-latency: urban12 0.49°→0.23°; fixes a double-count; monitor absorbs the prior when ON).
- **Slice 20** (`629d379`) + **20b** (`c85442e`) — `translation_only` source flag (pins rotation to prior → footgun fix: a velocity/Doppler source no longer commits a garbage rotation extrinsic) + min-over-observed lever conf + PER-AXIS lever commit. Persistence **v4** (pre-v4 cold-start). Real radar: footgun fixed; lever doesn't commit (sensor-noise-bound).

**Real-data covariance (KAIST):**
- **GPS R-side NEES sweep** (`28339aa`) — `cov_floor_m2=25` was 17–300× overconfident; **banked `cov_floor_m2=300`** into the URBAN RECIPE (`d12c419`). Honest-R differs ~7× cross-drive + confounded by GPS-vs-GT FRAME offset (~8–26 m) → a per-fix/multipath R MODEL is the principled follow-up. Also CLOSED the urban12 one-window regression (`127f48b`, 211→110.8 m via cov_floor=300).
- **Motion-aware (Ad-shape) predict-Q DE-RISKED → NO-GO** (`0d4c717`) — diagonal motion-Q doesn't collapse the per-drive predict-only NEES spread; true Ad-transport-shaped Q deferred-second-order on an already-safe NEES.
- **Slice-14 CONFIG placeholder sweep DONE** (`e51788d`) — `match_metric`/`excitation_min_var`/`kf_process_noise`/`straight_omega_max`/`turn_omega_min` swept+justified (KEEP); `straight_trans_min`=0.05 flagged cadence-coupled (use ~0.02 for slow/high-cadence rigs).

**MANY-SENSOR + RADAR/CAMERA real-data threads (`SURROUND_MULTISOURCE_TEST.md`, `RADAR_SCAN_ODOMETRY.md`):**
- Synth 9-source surround (4cam+2radar+3 real KAIST) — all calib DOFs recover. Real nuScenes 7-source fusion (5 radar Doppler + CAN wheel/imu) — bounded NEES 1.6/2.6.
- **Radar scan-matching odometry** (user-contributed descriptor algorithm, `61b79f3`; 3D-ready) — recovers full R+t. On 2D nuScenes ARS408 rotation is below the noise floor (sensor-bound, not algorithm); **synth-4D radar (`345c782`) PROVES density/3D unlocks rotation** (yaw corr 0→0.995). Longer-baseline 2D fix falsified.
- **Real camera surround (`a39521d`+`fc74b72`+`2cf90fb`)** — monocular VO on nuScenes; **5/6 cameras recover their ROTATION extrinsic** (+5°→<3°) on real VO — the first real-data multi-extrinsic rotation recovery (the venue the heading-blind radar couldn't be). Camera LEVERS WARM with 10-scene accumulation (length was the limit, not observability) but don't COMMIT (planar-urban + disjoint-20s-fragment data limit; needs continuous non-planar drives — nuScenes can't provide).

**NEXT (data queued): REAL 4D RADAR — RADIal.** User downloaded RADIal (Valeo, open) to **`C:\workspace\data\RaDIaL\`**: `Ready_to_use/` (123 GB, extracted — the processed/synced data we run on), `scripts-*.zip` (81 MB, NOT extracted — the devkit/loaders to read first), `miscellaneous-*.zip` (4 KB, NOT extracted — calib/extrinsics). PLAN: extract scripts+miscellaneous → recon the `Ready_to_use` radar format (sparse x,y,z+Doppler point cloud? or a range-Doppler SPECTRUM needing a CFAR peak-extraction step — RADIal is a detection dataset, so likely the latter) + the GPS(→GT pose) + CAN-speed/IMU(→reference odom) → build `tools/radial_to_csv.py` (mirror nuscenes_to_csv) → run the PROVEN scan-match 3D path on real imaging radar = confirm the synth-4D rotation unlock on real data. Do NOT start until the user says go.

### Slices 16–19 + recalibration (pre-session delta — pushed)

| What | Outcome |
|---|---|
| **Slice 16 — centroid sub-bin readout** (`c237bfb`/`00978c8`, D25) | The calib precision floor was `mode()`'s parabolic pull-toward-bin-center (~70–80% of the sub-bin offset; the historical "0.10°" was luck). Opt-in `subbin_centroid` (all 5 calib hists): committed yaw 0.229°→**0.0102°**, time-offset 2.2 ms→**2.6 µs**. |
| **Slice 17 — rot3d turn-regime rotation extrinsic** (`f1eee59`/`2d22c28`/`f013320`, D26) | Axis-correspondence Wahba hand-eye in Phase-2 (two-axis BBw gate = the spine; prototype `tools/proto_rot_handeye.py`). EuRoC drone: **9.1e-6° committed** (was unobservable); KAIST planar byte-identical (gate honest). |
| **Slice 17b — joint lever+scale** (`4138c3d`/`0ea764a`, D27) | Hand-eye translation row is linear in `[t_X; 1/s]` → 4-unknown LS: **scale committed FROM TURNS on a drone** (1.07962 vs 1.08), lever sub-mm with scale unknown. Schur-marginal vote gating (≥3 distinct excitations), scale2 out-of-range skip guard. |
| **Slice 18 — median-coupled multi-source bias** (`def0cee`/`2dc6678`/`67b0242`/`57fc874`, D28) | Spike-5a's scalar Jacobian **FD-FALSIFIED** (was verified against the pinning median) → exact matrix influence Ω_i landed. Opt-in; sim-validated. **Real-data boundary documented**: position-only GPS at honest priors cannot observe per-source yaw bias (C4 suppresses the carrying rows; sim-scale priors = junk sink, +9870°/h). **Keeper: out-of-range calib-vote SKIP guards (so3/xyz/rot3d)** — wandering inputs were edge-clamping votes that COMMITTED at ±0.984375; guards active always. |
| **q_floor recalibrated on real data** (`e37256b`) | Predict-only NEES ∝ 1/q_floor exactly; per-drive consistent floor differs 11× (scalar floor can't be cross-drive consistent → `Ad`-shape still the principled fix). **Recommended urban: `1e-5 ×3, 1e-6 ×3`** — urban12 transient 409→129 m, all drives net-better. |
| **Slice 19 — per-channel split median** (`6ebfcc5`/`6c23894`, **D29 = D3 amendment, USER-approved**) | Independent SO(3)/ℝ³ Weiszfeld (full D3 safeguards per channel), per-channel spread→Q (**`q_scale_split=3.0`** — coupled 0.7 is overconfident under split, NEES 21.7: the mixed-unit spread had padded both Q blocks), LOO-floored cross-channel veto, per-sensor `rot_weight_prior`. **KAIST + FOG prior 10: urban12 fused heading rms 5.4° ≈ the FOG floor** (was a 5–30° band + 61° spike), local rot p50 **0.23°**; urban07 better everywhere; urban17 neutral. Default OFF (flip after soak). |

**URBAN TRAJECTORY RECIPE (current best)**: `q_floor = 1e-5 1e-5 1e-5 1e-6 1e-6 1e-6` + **`cov_floor_m2 = 300`** (was 25 — the R-side NEES sweep `28339aa`/`GPS_R_NEES_SWEEP.md` found 25 is 17–300× overconfident; ~300 cuts NEES, IMPROVES urban07 local drift, and shrinks the urban12 worst window 211→~110 m — strictly better, config-only) + `correction_rot_suppress_kappa = 0.8` + `split_median = true` + FOG-grade source `rot_weight_prior = 10`. Slice docs: `SLICE16_…`/`SLICE17_…`/`SLICE17B_…`/`SLICE18_…`/`SLICE19_…` (each has the validation tables + the authoritative amendments).

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
- **Histogram readout precision** (Slice 16): the parabolic `mode()` pulls ~70–80% toward the peak-bin center — any sub-bin precision claim needs `subbin_centroid=true`. **Out-of-range calibration votes are SKIPPED, never edge-clamped** (Slices 17b/18 guards): boundary-bin mass committing at ±0.984375 destroyed urban12 once; do not "simplify" the guards back to clamping.
- **Split median** (Slice 19, D29): `split_median=true` MUST pair with its own `q_scale_split` (default 3.0) — the coupled `q_scale=0.7` is grossly overconfident under split (NEES ~21, pinned). Veto floors (0.01 rad / 0.02 m) are per-WINDOW scales — revisit for ≫1 s windows.
- **Multi-bias** (Slice 18, D28): scale `multi_bias_cov0`/`bias_process_noise` to the SENSOR (sim-scale defaults on real data = junk sink, one GPS fix → +9870°/h); per-DOF pn 0 pins a DOF; do NOT expect yaw-bias learning from position-only GPS (the documented boundary).
- **Persistence format v3** (17b) — pre-v3 blobs cold-start by design; Slice-16's flag is in the config-hash (pre-16 configs hash differently).
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
- 11b Option B is now BUILT (Slice 18, D28 — opt-in, with the spike-5a Jacobian corrected to the exact Ω_i and the real-data observability boundary documented; do not re-attempt GPS-fed yaw de-biasing without new evidence). Still not built: ROS node + durable fsync (13b). Partial: validation (14) — NEES/NIS/golden/init-P/median-fix/q_scale + the predict-side q_floor real-data recalibration (`e37256b`) DONE; remaining = the non-covariance placeholders (excitation_min_var, kf_process_noise, match_metric, calib gates), the `Ad` distance-shape model (now quantified: 11× per-drive floor spread), and the GPS R-side honest-NEES sweep. Historical: the two latent calib bugs the median fix exposed are RESOLVED (`3bd91e2`; sparse-vote artifact documented).

---

## 7. Resume in one line

Slices 0–20b landed + pushed (HEAD `baf6893`). **The calibration story is COMPLETE** (every calib DOF meets the precision targets on ground AND 3D — §3 recipes), **the heading prize is delivered** (Slice 19), and the **MANY-SENSOR + radar/camera real-data threads are done or precisely diagnosed** (SESSION DELTA block in §3). The tractable covariance/calibration backlog is cleared. Remaining work is data-gated or research-y:

- **IMMEDIATE NEXT (data in hand, awaiting user "go"): REAL 4D RADAR — RADIal.** Data at **`C:\workspace\data\RaDIaL\`** (`Ready_to_use/` 123 GB extracted; `scripts-*.zip` 81 MB + `miscellaneous-*.zip` 4 KB NOT extracted). Build `tools/radial_to_csv.py` → run the PROVEN scan-match 3D path on real imaging radar → confirm the synth-4D rotation unlock (yaw corr 0→0.995) on real data. RADIal is a DETECTION dataset → the radar may be a range-Doppler SPECTRUM needing a CFAR peak-extraction step (recon `scripts`/`Ready_to_use` first). The non-planar/elevation content (if any) could also warm the radar LEVER. See §3 SESSION DELTA + `RADAR_SCAN_ODOMETRY.md`.
- **GPS per-fix / multipath-aware R MODEL + GPS-vs-GT frame-offset characterization** — the principled R-side follow-up the sweep surfaced (`GPS_R_NEES_SWEEP.md`); research-y (the flat VRS per-fix variance has no multipath signal; the implausibly-large honest-R is mostly the GPS⊥GT frame offset). An un-gated-NIS calibration methodology is the cleaner R target.
- **Camera LEVER commit** — observable (warms with accumulation) but needs a CONTINUOUS NON-PLANAR driving dataset (nuScenes can't provide; trainval is also 20 s + planar). Data-bound.
- **Noisy `ego_pose`-as-GPS on nuScenes** — activates the heading monitor (inert without GPS) + bounds the 7-source drift; data in hand, small build (offered, not yet built).
- **Deferred/closed**: `Ad`-shape predict-Q (de-risked NO-GO, `0d4c717` — diagonal doesn't collapse the spread, transport-shaped is second-order on an already-safe NEES); 13b ROS node + durable `fsync` (env-blocked); default flips (`split_median`/`heading_monitor` ON) need real-world soak first (re-bases golden/cov-cal). Closed: per-source bias via position-only GPS (D28); scale/rot-extrinsic on 3D (17/17b); the urban12 window (`127f48b`); the non-cov placeholder sweep (`e51788d`).

Pick one; for a core change, design-first (the 15→19 series shows the pattern: investigate → falsifiable prototype → design doc → implement → review → validate full-drive, and **never weaken the trust apparatus** — cov-cal bands, observability self-tests, golden, byte-identical defaults). Verify with §2.
