# adapters — relaxed edge

Middleware and platform integrations. **Not** subject to the strict-core profile (these never ship to the ECU): std containers, exceptions, file I/O, and third-party deps are all fine here.

## Planned adapters (ISSUES.md Slice 13)
- `config_yaml/` — parse YAML/JSON/ROS-param into `ofc::Config` (the core stays parser-free, D19).
- `ros/` — node wrapping `ofc::Estimator`: subscribe odometry sources, publish `nav_msgs/Odometry` + TF + a calibration message; implement `ISource` / `ICorrection` over ROS topics.
- `threading/` — optional background-thread wrapper around the caller-pumped core (D14).
- `persistence_file/` — double-buffer ping-pong file backend for `serialize/deserialize` (D23): two files + sequence counter, version tag + checksum, config-hash guard.
- `gps/` — `ICorrection` for absolute position/pose (D22).

Each adapter depends only on the public `ofc/core` headers.

## GPS correction adapter (`gps_correction.hpp`, D22, Slice 11b/13)

A production `ofc::ICorrection` that turns real GPS fixes into the core's **dim=3 position**
correction — the same Mahalanobis-gated ESKF update path the sim ref (`ofc::sim::SyntheticAbsoluteRef`)
exercises — so a live GPS stream removes the fused-odometry drift. Position-only; a full-pose
(dim=6) or orientation fix is a noted extension (same machinery, a different `h`/`H`).

- `GpsFix` — one fix: WGS-84 `lat_deg`/`lon_deg`/`alt_m`, an ENU position covariance `cov_enu`
  (m²), a `stamp` (ns, same clock as `Estimator::step`), and a `valid` no-fix sentinel.
- `GpsConfig` — the datum (`has_datum` + `datum_lat/lon/alt`, or a **lazy** latch of the first
  valid fix when `has_datum` is false), the rigid `odom_from_enu` alignment (ENU → the estimator's
  odom frame; identity = "GPS ENU already aligned to odom" — a wrong rotation biases every
  residual, the same footgun as the sim ref's `window_s`), the antenna `lever_arm_base` (m, in the
  base frame), and an optional `cov_floor_m2` R-inflation floor.
- `GpsCorrection` — construct with a `GpsConfig`, register via `Estimator::add_correction(&gps)`,
  and `submit_fix(...)` from the GPS subscriber. The adapter holds the most-recent **unconsumed**
  valid fix (latest-wins); `evaluate()` (called inside `step()`) emits it **once**, then returns
  false until a newer fix arrives. An invalid fix is ignored.

Pipeline: geodetic → ECEF (WGS-84) → ENU about the datum → odom (`z = odom_from_enu ∘ p_enu`).
Measurement: `h(x) = x.pose.t + x.pose.R·lever`, `residual = z ⊟ h(x)`,
`H = [ R | −R·[lever]× | 0 ]` (3×12), `R_odom = R_align·(cov_enu + cov_floor·I)·R_alignᵀ`.
`out.stamp` is the **fix's** acquisition time (a deviation from the sim ref, which stamps the
frontier time because its "fix" is sampled at the frontier).

> Note on `q_floor`: the GPS gain depends on the per-step process noise. With a near-zero
> `Config::q_floor` (and identical sources giving zero adaptive-q spread) the position covariance
> never re-inflates between fixes, so the Kalman gain collapses and a fix barely pulls. A realistic
> translation `q_floor` (the odometry-is-uncertain floor) is what lets a GPS fix actually correct
> the drift — see the end-to-end test.

## Real-data ingestion: CSV source + GT track + replay harness + `ofc_replay` CLI (Slice 13)

Until now the system was exercised only by the **sim rig** (synthetic sources sampled from ground
truth). This scaffolding lets the sim-validated system run on a **real dataset**: a generic-CSV
`ISource`, a GT-track reader, an offline replay+eval harness, and a thin CLI. The user converts any
public dataset (KITTI/EuRoC/…) to the CSV schema below with a small script. Dependency-free (Eigen
+ doctest only — the CSV is parsed by hand, **no csv library**).

### `CsvSource` (`csv_source.hpp`) — a file-backed `ISource`

Loads a CSV into a `SourceBuffer` (the **same** buffer the sim's `SyntheticSource::build_buffer`
uses) and answers `query(t0,t1)` straight from it — so the CSV ingestion path is **equivalent** to
the in-memory path (proven by the equivalence test). The native form is chosen by a `# form: <tag>`
header line **or** the `CsvSourceConfig::form` ctor arg (`force_form=true` makes the arg win and
rejects a disagreeing tag). Quaternion is **(qw, qx, qy, qz)** w-first, normalized on read.

CSV schema (numbers decimal; `#`/`//`/blank lines skipped; a leading non-numeric header row
tolerated; whitespace trimmed; **comma OR whitespace** delimited):

```
# form: absolute        (or: increment, or: twist)
# absolute / increment:
t_ns, x, y, z, qw, qx, qy, qz [, var_tx, var_ty, var_tz, var_rx, var_ry, var_rz]
# twist:
t_ns, vx, vy, vz, wx, wy, wz [, var_tx, var_ty, var_tz, var_rx, var_ry, var_rz]
```

* `t_ns` — sample time (ns, int64), strictly increasing.
* **absolute** — `x,y,z`+quat is the cumulative pose in the source's own integrated frame.
  **increment** — the per-step relative motion from the previous row (the first row is the identity
  seed). **twist** — body linear `vx,vy,vz` (m/s) + angular `wx,wy,wz` (rad/s).
* The optional trailing **6** numbers are the per-step increment covariance **diagonal** in
  `[trans; rot]` order (for `twist`, the per-**second** twist covariance). Absent → the buffer's
  modeled-cov path is used (combined per `Config::conf_combine`).

A malformed row (bad number/quaternion, wrong column count, non-increasing stamp) returns a
`Status` + a human `error()` string with the line number; `< 2` rows → `NoData`.

### `CsvGtTrack` (`csv_source.hpp`) — a GT absolute-pose track (eval only)

```
t_ns, x, y, z, qw, qx, qy, qz       (no covariance columns; extra columns ignored)
```

`pose_at(t)` returns the **linearly SE(3)-interpolated** (`se3::interpolate`) pose between the two
bracketing samples, clamped at the ends.

### `ReplayHarness` (`replay_harness.hpp`) — offline replay + eval

The real-data analogue of the sim `Rig`. Given a `Config` + a set of `CsvSource`s + an optional
`GpsCorrection` (+ a sorted `std::vector<GpsFix>`) + an optional `CsvGtTrack`, it:
merges the source spans into a regular tick timeline at `tick_rate_hz` (folding in each GPS-fix
frontier tick), submits each GPS fix just before the step whose frontier reaches its stamp, calls
`Estimator::step()` at each tick, and records per published step the frontier pose (t + w-first
quat) + the 12×12 cov diagonal, tip, per-source `CalibSnapshot`, and `CorrectionDiag`. The GT pose
is anchored exactly as the sim Rig (at the GT pose at `first_frontier − window_s`, DESIGN §7) so the
frontier-vs-GT error compares like-for-like.

With a GT track it also computes, per step:
`trans_err_m = ‖frontier.t − gt.t‖`, `rot_err_rad = ‖so3::log(gt.Rᵀ frontier.R)‖`, and the **6-DOF
pose NEES** `e = se3::log(T_estⁱⁿᵛ ∘ T_gt)`, `nees = eᵀ (P_pose)⁻¹ e` (P_pose = `frontier.cov(0..5)`)
— the **exact** convention of `tests/test_validation.cpp pose_nees`. The **NIS** of accepted GPS
fixes is read from `CorrectionDiag.last_nis` on steps with `corr_applied > 0`. `RunSummary`
aggregates max/RMS/tail drift, ensemble-mean pose NEES (target ~DOF=6), and mean accepted-GPS NIS
(target ~DOF=3). `write_results_csv()` emits the per-step CSV + a trailing `# summary:` block.

### `ofc_replay` CLI (`tools/ofc_replay/main.cpp`)

`ofc_replay <manifest.ini> [results_out.csv]` — reads a manifest, loads the CsvSources / GPS fixes /
GT track, runs the harness, writes the results CSV, and prints the summary to stdout. Built only
when `OFC_BUILD_ADAPTERS=ON` (so the default core build is unaffected).

The **manifest** is the `config_loader` INI format **extended** with per-source CSV paths, a GPS
block, a GT block, and replay knobs. Extension keys (`csv`/`form` in `[sensor.N]`; the whole
`[gps]`/`[gt]`/`[replay]` sections) are consumed by the CLI; every other line is passed through to
the existing `ConfigLoader` (so the core `Config` is validated normally). Concrete example:

```ini
[global]
tick_rate_hz = 50
fusion_delay_s = 0.05
window_s = 0.10
max_sources = 2
reference_sensor_id = 0
cold_start = median_from_start

[sensor.0]
id = 0
is_reference = true
csv = data/odom0.csv         ; EXTENSION: the CSV for this source
form = increment             ; EXTENSION: absolute|increment|twist (else the file's # form: tag)

[sensor.1]
id = 1
csv = data/odom1.csv
form = increment

[gps]                        ; EXTENSION (optional): the absolute ref
csv = data/gps.csv           ; rows: t_ns, lat_deg, lon_deg, alt_m [, var_e, var_n, var_u]
datum_lat_deg = 47.0         ; omit the 3 datum_* to lazily latch the first fix
datum_lon_deg = 8.0
datum_alt_m = 400.0
odom_from_enu_yaw = 0.0      ; ENU->odom yaw (rad)
lever_x = 0.0
lever_y = 0.0
lever_z = 0.0                ; antenna lever in the base frame (m)
cov_floor_m2 = 0.0

[gt]                         ; EXTENSION (optional): GT track (eval only)
csv = data/gt.csv            ; rows: t_ns, x, y, z, qw, qx, qy, qz

[replay]                     ; EXTENSION (optional): harness knobs
tail_window_s = 1.0
warmup_steps = 20
out = results.csv
```

GPS fix CSV schema (one fix per row): `t_ns, lat_deg, lon_deg, alt_m [, var_e, var_n, var_u]`
(ENU position variances in m²; absent → 1 m² isotropic). The CLI submits fixes against the harness
timeline so `step()` consumes each at the matching frontier.
