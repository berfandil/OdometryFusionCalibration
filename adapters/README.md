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
