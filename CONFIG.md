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
| `q_scale` | double | **0.5** | (0, ∞) | Multiplier on the spread-derived Q. **Calibrated (Slice-14 sweep, `983fa65`)** — swept over {1.0,0.5,0.3,0.2,0.15,0.1} × {nees_traj, mixed, turning, straight} × {1×,2× noise}; `0.5` is the SAFETY-FIRST choice (max pessimism cut subject to NEVER overconfident, with margin): worst-case no-ref pose NEES ~2.9 (1×) / ~3.5 (2×) vs DOF 6 — `0.3` reaches ~5.8 at 2× (no headroom), `0.2` is already overconfident (~7.1 on `turning`). Was `1.0` (a placeholder; NEES ~1.0 post-approach-A). Not pushed to NEES≈6 (≈`q_scale 0.17`) because that overfits one trajectory + the sim under-states real model mismatch — the gap to 6 is a deliberate safety margin. |
| `adaptive_q_source_reduction` | bool | true | — | **Slice-14 approach A**: divide the spread-derived adaptive Q by `n_eff` = the participating fusion-median source count (`q_pose = q_scale·spread²/max(1,n_eff)·I + q_floor`). The geometric median of `n` agreeing sources has ~`1/n` the per-window variance of one source, so the spread (which measures DISAGREEMENT, not fused accuracy) over-states Q by ~`n`. Cuts the predict-only no-ref NEES pessimism ~17× → ~6× (ensemble NEES 0.35 → ~1.0, still conservative — never overconfident). `false` recovers the pre-A covariance exactly. Floor is NOT divided. Runtime knob (excluded from the persistence `config_hash`). **Interaction**: ON shrinks predict-only `P` between absolute fixes → smaller Kalman gain (the `q_floor`↔GPS coupling below) for rigs whose sources actually disagree. |
| `q_floor` | double[6] | small | ≥0 | Per-axis minimum process noise. **GPS/absolute-ref coupling (Slice 11b)**: `q_floor` is the floor `q_pose` collapses to when the inter-source spread is ~0 (e.g. two near-identical sources, or a single source). If it is too small the predict-only `P` never re-inflates between absolute fixes, so the Kalman gain on a GPS/absolute correction → ~0 and fixes barely pull the estimate. If absolute fixes seem ignored, raise the translation `q_floor` (the GPS adapter end-to-end test uses `1e-3`). **Slice-14 sweep verdict (`983fa65`)**: the trade is a genuine deployment choice, so the small default is KEPT. On a spread>0 (disagreeing-source) rig a larger floor only ADDS no-ref pessimism (e.g. no-ref NIS 2.94→0.78, drift-removal tail-err 0.020→0.071 m as trans floor 1e-6→1e-2); the floor only HELPS the zero-spread absolute-ref/GPS rigs (identical sources → `q_pose=floor`). **Guidance**: keep small (default) for no-ref / multi-source consistency; raise translation floor (~`1e-3`) for GPS/absolute-ref rigs whose sources agree. |

## 4. Source weights (`WeightConfig`)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `reliability_ema_alpha` | double | 0.02 | (0, 1] | EMA rate for the variance-based reliability track (smaller = slower). |
| `reliability_floor` | double | 0.2 | (0, 1] | Minimum reliability multiplier (Slice 9): a noisy source is downweighted but never collapses. |
| `reliability_cap` | double | 5.0 | [1, ∞) | Maximum reliability multiplier (Slice 9): a clean source is upweighted but never dominates. |
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
| `so3_hist` | HistogramConfig | — | — | Phase-1 direction (3-channel so(3)). |
| `roll_hist` | HistogramConfig | — | — | Phase-2 roll (circular S¹). |
| `xyz_hist` | HistogramConfig×3 | — | — | Phase-2 translation (3×1-D ℝ³). |
| `scale_hist` | HistogramConfig | — | — | Per-source scale (1-D). |

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

## 9. Per-sensor (`SensorConfig`, one per source)

| Knob | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `id` | int | — | [0, max_sources) | Source index. |
| `prior_extrinsic` | SE3 | identity | — | Mount prior (also the histogram basepoint for Phase 1). |
| `prior_scale` | double | 1.0 | (0, ∞) | Translation-scale prior. |
| `prior_time_offset_s` | double | 0.0 | [−max_lag, max_lag] | Clock-offset prior. Sign: **positive ⇒ source clock ahead of base** (its `[t0,t1]` maps to true `[t0+off, t1+off]`). |
| `weight_prior` | double | 1.0 | (0, ∞) | Initial trust weight. |
| `native_confidence` | bool | true | — | Use the provider's reported Σ if present. |
| `modeled_noise_per_m` | double | tuned | ≥0 | Synthetic translation noise (per metre). |
| `modeled_noise_per_rad` | double | tuned | ≥0 | Synthetic rotation noise (per radian). |
| `per_sensor_kf` | bool | false | — | Enable the per-sensor fixed-lag RTS twist smoother feeding calibration (Slice 10, **wired**). On ⇒ calibration runs at the deeper frontier `now − delay − L` with two-sided (zero-phase) smoothed twists; OFF ⇒ byte-identical to no-smoothing. |
| `kf_process_noise` | double | tuned | ≥0 | Smoothing strength `q` (larger = looser/closer-tracking; only `q/r` shapes smoothing, `r` fixed at 1.0). **Note**: a single shared smoother uses the MAX `q` over enabled sources (not yet per-source). |
| `scale_calib` | bool | true | — | Calibrate scale; if false, fixed at `prior_scale`. |
| `bias_states` | bool | false | — | Add nuisance-bias states (GPS/INS-style drift removal); enable for raw-IMU sources. **Slice-11b Option A**: when set on the SINGLE source that drives the predict alone, the ESKF augments to 18-DOF `[pose;twist;bias]` and an absolute-ref update removes the bias via the pose↔bias cross-covariance. `>1` source set → `InvalidConfig`. |
| `bias_process_noise` | double | 1e-4 | ≥0 | Random-walk process noise for the per-source bias state (Slice-11b; only active when `bias_states` + the source drives the predict). Larger = the bias tracks faster / is trusted less. |
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

*Defaults marked "tuned" are placeholders pending the simulation-rig sweep (see [`ISSUES.md`](./ISSUES.md), validation slice). `q_scale` was calibrated in the Slice-14 covariance sweep (`983fa65`, now `0.5`); the remaining placeholders (`excitation_min_var`, `kf_process_noise`, `match_metric`, the straight/turn calibration gates) are a separate future pass — they are calibration/smoothing knobs the observability self-tests already pin functionally, not covariance-consistency knobs.*
