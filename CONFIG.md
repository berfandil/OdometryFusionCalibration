# Configuration Reference

Every tunable knob, grouped by subsystem. Config is one validated POD struct passed at construction (see [`DESIGN.md`](./DESIGN.md) §9); the core preallocates from it once and never allocates again. File parsing (YAML/JSON/ROS-param) is an **adapter** that builds this struct — the core itself is parser-free.

Conventions: **type** is the logical type; **default** is the shipped value; **range** is the validated bound (`validate()` rejects out-of-range). `[strict]` knobs size preallocation and cannot change after `init()`.

---

## 1. Runtime / timing (`RuntimeConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `tick_rate_hz` | double | 50.0 | (0, 1000] | Fusion step cadence (caller is expected to pump at ~this rate). |
| `fusion_delay_s` | double | 0.05 | [0, 2] | Causal frontier offset `now − delay`; must cover the slowest source latency. |
| `window_s` | double | 0.1 | (0, 5] | Bootstrap/lookback interval for the first fuse (and max lookback). Steady-state predict integrates `[last_frontier, frontier]`, so tick cadence need not equal this. |
| `calib_lag_s` | double | 0.2 | [0, 10] | Extra lag `L` for the calibration frontier (`now − delay − L`); enables two-sided smoothing. |
| `late_sample_policy` | enum | `drop` | {`drop`, `reintegrate`} | Handling of samples arriving behind the frontier. |
| `max_sources` | int | 8 | [1, 32] | **[strict]** Sizes per-source buffers/structures. |
| `buffer_capacity` | int | 512 | [16, 65536] | **[strict]** Ring-buffer depth per source (samples). |

## 2. Gauge & frames (`GaugeConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `base_frame` | string/id | `"base"` | — | Name of the user-defined base frame. |
| `reference_sensor_id` | int | 0 | [0, max_sources) | The pinned gauge anchor; its extrinsic is held at its prior. |
| `cold_start_mode` | enum | `reference_only` | {`reference_only`, `median_from_start`} | Fusion behavior before calibration is confident. |

## 3. Median fusion (`MedianConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `weiszfeld_max_iters` | int | 10 | [1, 100] | Hard iteration cap (bounded WCET). |
| `weiszfeld_tol` | double | 1e-6 | (0, 1) | Convergence tolerance (step norm). |
| `weiszfeld_eps` | double | 1e-9 | (0, 1) | ε-regularization of the `1/d` weight (vertex singularity guard). |
| `metric_lambda` | double | 1.0 | (0, ∞) | Rotation-vs-translation weight in the split SO(3)/ℝ³ metric. |
| `min_sources_warn` | int | 3 | [1, ∞) | Lifecycle NOMINAL threshold (Slice 3): a fused step with `n >= min_sources_warn` participating sources is NOMINAL, else DEGRADED; below 3 outlier rejection also degrades. `validate()` enforces only `>= 1`; a value above `max_sources` is legitimate (NOMINAL is then never reached). |
| `adaptive_q` | bool | true | — | Derive process noise Q from the inter-source spread. |
| `q_scale` | double | **0.7** | (0, ∞) | Multiplier on the spread-derived Q (`q_pose = q_scale·spread²·I6 + q_floor`). **Calibrated against the TRUE interior median (`1142e41`)** — swept over {0.5,0.7,1.0,1.5,2.0} × {nees_traj, mixed, turning, straight} × {1×,2× noise}; `0.7` is the SAFETY-FIRST choice: it gets the worst-case no-ref pose NEES as close to DOF=6 as possible FROM BELOW while NEVER overconfident — worst-case ~4.85 (1×) / ~3.9 (2×) vs DOF 6 (`0.5` hits ~6.8 → overconfident). **Sign note**: with the true median (and `adaptive_q_source_reduction` now removed) SMALLER `q_scale` → LARGER NEES (less Q, tighter P), the OPPOSITE of the historical pinning-median behaviour. The gap to 6 is a deliberate safety margin (the sim under-states real model mismatch). History: was a placeholder `1.0`, then `0.5` (calibrated against the buggy pinning median in `983fa65`, now superseded). |
| `split_median` | bool | false | — | Per-channel split median (D29/Slice 19): independent rotation/translation Weiszfeld consensus with per-channel weights + per-channel spread→Q. OFF = the coupled D3 solver, byte-identical. Enables `rot_weight_prior` to express heading-grade sources (the urban12 heading fix: fused rot p50 0.23°). |
| `split_veto` | bool | true | — | Cross-channel outlier veto on the split path: a source ≥3× the (floored leave-one-out) spread off in one channel is ×0.1-weighted in the other (one bounded re-solve). Floors `0.01 rad / 0.02 m` stop near-coincident hair-triggers; revisit for ≫1 s windows. |
| `q_scale_split` | double | **3.0** | (0, ∞) | The split path's own calibrated Q multiplier (used INSTEAD of `q_scale` when `split_median`) — the coupled 0.7 is grossly overconfident under split (NEES ~21, pinned): the coupled mixed-unit spread had padded both channels' Q. Calibrated like `q_scale` (band [4.0, 5.6], never-overconfident; post-veto-floor worst 4.998). |
| `heading_monitor` | bool | false | — | GPS-course heading-drift monitor (D29/Slice 19c, layer c; **requires `split_median=true`** — else `validate()` rejects). Ranks per-source heading drift against GPS course-over-ground and boosts the heading-grade source's ROTATION-channel weight (`·monitor_boost`), auto-discovering the FOG with no `rot_weight_prior`. OFF byte-identical (both paths AND the abstain state). Validated on KAIST: discovers the FOG on all three drives (boost 10 vs 1); the cost is discovery latency (early-drive abstain) so the config prior (layer a) wins when the heading-grade source is known. Loader: `[global] heading_monitor`. |
| `heading_monitor_boost_max` | double | 10.0 | ≥1 | The winner's rotation-channel boost cap for `heading_monitor` (`boost = clip(boost_max·max(min_score,floor)/max(score_i,floor), 1, boost_max)`; GPS-noise floor `kScoreFloor=2.42e-5 rad/s≈5°/h`). On KAIST 10 reproduces the manual `rot_weight_prior=10`. Loader: `[global] heading_monitor_boost_max`. |
| `q_floor` | double[6] | small | ≥0 | Per-axis minimum process noise. **GPS/absolute-ref coupling (Slice 11b)**: `q_floor` is the floor `q_pose` collapses to when the inter-source spread is ~0 (e.g. two near-identical sources, or a single source). If it is too small the predict-only `P` never re-inflates between absolute fixes, so the Kalman gain on a GPS/absolute correction → ~0 and fixes barely pull the estimate. If absolute fixes seem ignored, raise the translation `q_floor` (the GPS adapter end-to-end test uses `1e-3`). **Slice-14 sweep verdict (`983fa65`)**: the trade is a genuine deployment choice, so the small default is KEPT. On a spread>0 (disagreeing-source) rig a larger floor only ADDS no-ref pessimism (e.g. no-ref NIS 2.94→0.78, drift-removal tail-err 0.020→0.071 m as trans floor 1e-6→1e-2); the floor only HELPS the zero-spread absolute-ref/GPS rigs (identical sources → `q_pose=floor`). **Guidance**: keep small (default) for no-ref / multi-source consistency; raise translation floor (~`1e-3`) for GPS/absolute-ref rigs whose sources agree. **Real-data recalibration (2026-06-12, KAIST full-drive sweep)**: predict-only NEES scales exactly ∝ 1/q_floor (10×/decade; estimate path invariant); the per-drive consistent floor differs 11× between drives (urban07 ≈7e-6 vs urban17 ≈6e-7) — a scalar floor canNOT be cross-drive consistent (the deferred `Ad` distance-shape model remains the principled fix). **Recommended URBAN GPS-rig value: `1e-5 1e-5 1e-5 1e-6 1e-6 1e-6`** (replaces the earlier ad-hoc `1e-3`): urban07 local p50 2.37→1.35 m / rot rms 0.084→0.049; urban17 local max 28.6→11.8 / rot max 0.204→0.072; **urban12 mid-drive transient 409→129 m**, tail 1.99→1.39 m, applied fixes 1375→1494 (truer coast heading keeps returning fixes in-gate) — no GPS-gain starvation at this floor with `cov_floor_m2=25`. With-GPS NEES becomes more overconfident (R-side under-states urban multipath — the open R-side item, separate knob). |

## 4. Source weights (`WeightConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `reliability_ema_alpha` | double | 0.02 | (0, 1] | EMA rate for the variance-based reliability track (smaller = slower). Under `split_median=true` it governs BOTH per-channel tracks (Slice 19b). |
| `reliability_floor` | double | 0.2 | (0, 1] | Minimum reliability multiplier (Slice 9): a noisy source is downweighted but never collapses. Applies per channel under split (Slice 19b). |
| `reliability_cap` | double | 5.0 | [1, ∞) | Maximum reliability multiplier (Slice 9): a clean source is upweighted but never dominates. Applies per channel under split (Slice 19b). |
| `weight_floor` | double | 0.05 | (0, 1) | Minimum source weight (guarantees recovery; never collapses to 0). |
| `weight_cap` | double | 10.0 | [1, ∞) | Maximum source weight (no single source dominates). |
| `confidence_combine` | enum | `sum` | {`native_only`, `modeled_only`, `sum`, `max`, `weighted`} | How native Σ and modeled Σ combine. |
| `confidence_blend` | double | 0.5 | [0, 1] | Blend factor when `confidence_combine = weighted`. |

## 5. Time-sync (`TimeSyncConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `enabled` | bool | true | — | Master on/off (off if inputs are pre-synced). |
| `match_metric` | enum | `l2` | {`l1`, `l2`, `ratio`, `ncc`} | Cross-correlation cost between ‖ω‖ signals (worth sweeping). |
| `max_lag_s` | double | 0.1 | (0, 2] | Bounded search range for the offset. |
| `excitation_min_var` | double | tuned | ≥0 | Minimum variance of ‖ω‖ to accept a window's estimate. |
| `offset_hist` | HistogramConfig | — | — | Per-source offset histogram (see §8). |

## 6. Calibration gates (`GateConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `straight_omega_max` | double (rad/s) | 0.05 | ≥0 | Phase-1 straight gate: `‖ω‖ < ε_ω`. |
| `straight_trans_min` | double (m) | 0.05 | ≥0 | Phase-1 motion gate: `‖t‖ > δ_v`. Per-step **displacement** (not a speed) → cadence-dependent; tune per tick rate. |
| `turn_omega_min` | double (rad/s) | 0.2 | ≥0 | Phase-2 turn gate: `‖ω‖ > θ_ω` (distinct from straight gate). |
| `reverse_fold` | bool | true | — | Fold reverse motion into the consensus hemisphere before voting. |
| `ref_cross_check` | bool | false | — | Require the reference sensor to also read "straight fwd/back" (ice-slide / drone niche). |

## 7. Calibration phases (`CalibConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `phase2_strategy` | enum | `vs_fused_base` | {`vs_fused_base`, `pairwise_pinned_ref`} | Hand-eye strategy (the two are compared empirically). |
| `commit_concentration` | double | 0.6 | (0, 1] | `τ_commit`: peak mass ÷ total required to commit. |
| `commit_min_votes` | int | 200 | [1, ∞) | `N_min` votes required to commit. **Must be consistent with the target histogram's aging**: under `sliding_k` aging, `total()` saturates at ~`sliding_k`, so set `commit_min_votes ≤ sliding_k` (or use decay aging) or commit can never trigger. Also a vote-**mass** threshold (not a count) under any non-`one` `vote_weight` (votes are Σ-confidence/‖ω‖-scaled) — re-tune per `vote_weight`. |
| `commit_drop` | double | 0.4 | (0, commit_concentration) | `τ_drop`: hysteresis — re-open below this. |
| `vote_weight` | enum | `combo` | {`one`, `confidence`, `rotation`, `combo`} | Per-vote weight: `one`=1; `rotation`=‖fused_ω‖ (Phase-1 floors it since straight-gated ω≈0 → effectively no-op there; meaningful in Phase-2's turn regime); `confidence`=per-source Σ-confidence (the inverse-cov fusion weight); `combo`=rotation×confidence. Honored in both calibration phases. |
| `rot3d_enabled` | bool | false | — | Turn-regime FULL rotation-extrinsic (axis-correspondence hand-eye, D26/Slice 17): per-slot Wahba accumulator over windowed rotation-axis pairs `a = R_X·b`, two-axis observability gate (`λ_mid ≥ kAxisPairMin·λ_max` on `Σ b bᵀ` — rank-1/ground rigs never vote/commit), votes into 3 so(3) channels reusing `so3_hist` (basepoint-relative, contractive re-anchor). Committed rot3d supersedes the Phase-1/roll rotation publish + lets the source join the median under `ReferenceOnly`. Accumulator decay uses `so3_hist.decay_gamma` regardless of aging mode. In the persistence config-hash; persistence format v2 (pre-17 blobs cold-start). Loader: `[global] rot3d_enabled`. |
| `joint_lever_scale` | bool | false | — | Turn-regime JOINT lever+scale (D27/Slice 17b): Phase-2's lever LS becomes the 4-unknown `[t_X; κ=1/s_res]` system (the hand-eye translation row is linear in both), recovering the residual scale WITHOUT the straight regime; votes gated on marginal (Schur) information (`kJointMarginMin` — needs ≥3 distinct turn excitations); `s_res` votes into a per-slot scale2 histogram (`scale_hist` shape) with an out-of-range SKIP guard (`scale2_skipped(id)` counts withheld votes; a true residual outside the `scale_hist` range never commits on a turn-only rig — widen the range if that's your regime). Either scale path's commit folds into `prior_scale` + resets both scale hists. Persistence format v3. Loader: `[global] joint_lever_scale`. |
| `so3_hist` | HistogramConfig | — | — | Phase-1 direction (3-channel so(3)); also the Slice-17 rot3d channels. |
| `roll_hist` | HistogramConfig | — | — | Phase-2 roll (circular S¹). |
| `xyz_hist` | HistogramConfig×3 | — | — | Phase-2 translation (3×1-D ℝ³). |
| `scale_hist` | HistogramConfig | **bins=64, range=[0.5,1.5]** | — | Per-source scale (1-D). **Defaults to `[0.5,1.5]`, NOT the generic `[-1,1]`** (`3bd91e2`): a scale is a positive ratio clustering at 1.0; on `[-1,1]` a unit ratio hits the half-open upper boundary → last-bin clamp → `mode()` returns the boundary-bin center (~0.984) → a true unit scale commits as 0.984 (was a real bug — masked by the old pinning median). The range must strictly contain 1.0; `validate()` now enforces this. |

## 8. Histogram primitive (`HistogramConfig`)

Shared by every calibrated quantity (so(3), roll, xyz, scale, time-offset).

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `bins` | int | 64 | [4, 4096] | **[strict]** Bin count per channel. |
| `range_min` / `range_max` | double | — | — | Value range (per quantity; circular for roll). |
| `circular` | bool | false | — | Wrap-around bins (true for roll/yaw-like). |
| `aging_mode` | enum | `decay` | {`decay`, `sliding_k`} | Vote-aging strategy. |
| `decay_gamma` | double | 0.999 | (0, 1) | Per-vote bin decay (when `aging_mode = decay`). |
| `sliding_k` | int | 1000 | [1, 4096] | Window size (when `aging_mode = sliding_k`). Upper bound = the fixed-capacity ring (`Histogram1D::kMaxSlidingK`). |
| `vote_split` | bool | true | — | Linear-split a vote between the two nearest bins. |
| `subbin` | bool | true | — | Parabolic sub-bin peak interpolation. |
| `subbin_centroid` | bool | false | — | Opt-in (honored only when `subbin`): replace the parabola with the mass-weighted centroid over peak±1 — exact for split votes at any sub-bin position, removing the parabola's pull-toward-bin-center bias (~70–80% of the sub-bin offset; D25, Slice 16). In the persistence config-hash (a flip rejects stale restores). Loader: `[global] subbin_centroid` sets it on all five calib histograms. |

## 9. Per-sensor (`SensorConfig`, one per source)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `id` | int | — | [0, max_sources) | Source index. |
| `prior_extrinsic` | SE3 | identity | — | Mount prior (also the histogram basepoint for Phase 1). |
| `prior_scale` | double | 1.0 | (0, ∞) | Translation-scale prior. |
| `prior_time_offset_s` | double | 0.0 | [−max_lag, max_lag] | Clock-offset prior. Sign: **positive ⇒ source clock ahead of base** (its `[t0,t1]` maps to true `[t0+off, t1+off]`). |
| `weight_prior` | double | 1.0 | (0, ∞) | Initial trust weight. |
| `rot_weight_prior` | double | 1.0 | ≥0 | Slice 19 (`split_median` only): multiplies this source's weight in the ROTATION channel — declare a heading-grade source (FOG/fiber gyro) from its datasheet. All-zero rotation weights in a step fall back to uniform. Loader: per-sensor `rot_weight_prior`. |
| `translation_only` | bool | false | — | Slice 20 (D30): declare a VELOCITY/translation-only source (Doppler radar, optical-flow; per-step rotation uninformative, `R_B≈I`). PINS the rotation extrinsic to `prior_extrinsic` (skips Phase-1 direction / rot3d / roll calibration for this source → no garbage rotation commit) so the hand-eye lever LS uses a clean `R_X = prior`; keeps scale calibration. Requires an orthonormal `prior_extrinsic` rotation (`validate()`). Hashed (pre-20 blobs cold-start). Loader: per-sensor `translation_only`. NOTE (real-data): the footgun fix (no garbage rotation) is solid, but committed LEVER recovery on a planar real radar is NOT yet delivered (`min(cx,cy,cz)` pins on the unobservable z + per-window radar noise — see SLICE20 follow-ups). |
| `native_confidence` | bool | true | — | Use the provider's reported Σ if present. |
| `modeled_noise_per_m` | double | tuned | ≥0 | Synthetic translation noise (per metre). |
| `modeled_noise_per_rad` | double | tuned | ≥0 | Synthetic rotation noise (per radian). |
| `per_sensor_kf` | bool | false | — | Enable the per-sensor fixed-lag RTS twist smoother feeding calibration (Slice 10, **wired**). On ⇒ calibration runs at the deeper frontier `now − delay − L` with two-sided (zero-phase) smoothed twists; OFF ⇒ byte-identical to no-smoothing. |
| `kf_process_noise` | double | tuned | ≥0 | Smoothing strength `q` (larger = looser/closer-tracking; only `q/r` shapes smoothing, `r` fixed at 1.0). **Note**: a single shared smoother uses the MAX `q` over enabled sources (not yet per-source). |
| `scale_calib` | bool | true | — | Calibrate scale; if false, fixed at `prior_scale`. |
| `bias_states` | bool | false | — | Add nuisance-bias states (GPS/INS-style drift removal); enable for raw-IMU sources. **Slice-11b Option A**: when set on the SINGLE source that drives the predict alone, the ESKF augments to 18-DOF `[pose;twist;bias]` and an absolute-ref update removes the bias via the pose↔bias cross-covariance. `>1` source set → `InvalidConfig` UNLESS `multi_bias_enabled` (Slice 18, ≤4 sources). |
| `bias_process_noise` | double / 6×double | 1e-4 | ≥0 | Random-walk process noise for the per-source bias state. **Per-DOF since Slice 18**: 1 value = uniform, 6 values = per-twist-DOF `[v;ω]`; a DOF set to **0 is PINNED at zero** (excluded from estimation — e.g. `0 0 0 0 0 1e-10` = yaw-rate-only bias). **SCALE TO YOUR SENSOR** (with `multi_bias_cov0`): real yaw-rate biases are ~1e-4 rad/s; sim-scale values (1e-4..1e-6) on real data make the bias a junk sink (a single GPS fix kicked the estimate to +9870°/h on KAIST — see D28). |
| `multi_bias_enabled` (global) | bool | false | — | Slice 18 (D28): median-coupled multi-source bias via the exact Ω_i influence. Opt-in; Option-A behavior byte-identical when false. **Real-data caveat (D28)**: position-only GPS at honest priors cannot observe per-source yaw-rate bias (C4 suppresses the very gain rows that carry it) — validated INERT on KAIST urban; the win regime is sim-like rich corrections. |
| `multi_bias_cov0` (global) | double | 0.04 | >0 | Initial per-DOF bias prior variance ((rad/s)²/(m/s)²). Default is SIM-scale; for real sensors set ≈ (expected bias 1σ)² — e.g. 1e-6 ⇒ σ ≈ 206°/h yaw. Too large = junk sink (D28). |
| `is_reference` | bool | (id==ref) | — | Marks the gauge anchor. |

## 10. Absolute references (`AbsoluteRefConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `mahalanobis_chi2` | double | 9.0 | (0, ∞) | Innovation χ² (NIS) gate on absolute-ref updates: reject when `dᵀS⁻¹d > threshold` (GPS multipath / outliers). `validate()` enforces `> 0` (Slice 11). **Now per-`n` (Slice 11b, `e8491dd`)**: this value is the **base threshold tuned at the n=3 position-fix DOF**; the estimator scales it by the χ²-quantile ratio `q[n]/q[3]` via `Eskf::chi2_gate(base, m.dim)` (const `kChi2Q95`, 0.95 confidence) so every measurement DOF `n∈1..6` gates at the same confidence — n=3 unchanged, n=6 ≈ 1.61×base. One knob, DOF-invariant confidence. The dim=3 fixes shipping now are behaviorally identical to the old single-scalar gate; this becomes load-bearing once a 6-DOF/mixed `ICorrection` lands. (The `chi2=100` value in the drift-removal test is an unrelated covariance-pessimism test artifact, not a production value.) |
| Plugin-specific extrinsic / noise | — | — | — | Supplied by the concrete adapter (GPS, map-match). |

## 11. Lifecycle & persistence (`PersistenceConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `enabled` | bool | false | — | Warm-restart persistence on/off. |
| `period_s` | double | 5.0 | (0, ∞) | Snapshot frequency. |
| `scope` | enum | `calibration_only` | {`calibration_only`, `full`} | Persist calibration state, or also raw buffers. |
| `path_a` / `path_b` | string | — | — | Double-buffer ping-pong file paths (adapter-side). |
| `config_hash_guard` | bool | true | — | Reject persisted state whose config/prior signature mismatches the current rig. |

## 12. Output (`OutputConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `emit_predicted_tip` | bool | true | — | Include the predict-only extrapolation to `now`. |
| `emit_diagnostics` | bool | true | — | Populate per-source residuals/weights/gating in the result. |
| `tip_cov_inflation` | double | 1.5 | [1, ∞) | Σ inflation factor applied to the predicted tip. |

---

*Defaults marked "tuned" are placeholders pending the simulation-rig sweep (see [`ISSUES.md`](./ISSUES.md), validation slice). `q_scale` was recalibrated against the true interior median (`1142e41`, now `0.7`); the remaining placeholders (`excitation_min_var`, `kf_process_noise`, `match_metric`, the straight/turn calibration gates) are a separate future pass — they are calibration/smoothing knobs the observability self-tests already pin functionally, not covariance-consistency knobs.*
