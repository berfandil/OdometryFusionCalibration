# adapters — relaxed edge

Middleware and platform integrations. **Not** subject to the strict-core profile (these never ship to the ECU): std containers, exceptions, file I/O, and third-party deps are all fine here.

## Planned adapters (ISSUES.md Slice 13)
- `config_yaml/` — parse YAML/JSON/ROS-param into `ofc::Config` (the core stays parser-free, D19).
- `ros/` — node wrapping `ofc::Estimator`: subscribe odometry sources, publish `nav_msgs/Odometry` + TF + a calibration message; implement `ISource` / `ICorrection` over ROS topics.
- `threading/` — optional background-thread wrapper around the caller-pumped core (D14).
- `persistence_file/` — double-buffer ping-pong file backend for `serialize/deserialize` (D23): two files + sequence counter, version tag + checksum, config-hash guard.
- `gps/` — `ICorrection` for absolute position/pose (D22).

Each adapter depends only on the public `ofc/core` headers.
