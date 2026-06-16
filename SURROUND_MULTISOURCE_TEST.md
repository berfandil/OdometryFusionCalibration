# Many-source fusion test — surround cameras + radars + real KAIST (9 sources)

**Goal (user, 2026-06-14)**: fuse many heterogeneous sensors at once (a 4-camera surround view + IMU + wheel-odo + radars), not just the 3 co-located KAIST sources. This is the first test of **multiple sources with genuinely different extrinsics** — the 4 surround cameras have 0/-90/180/90 deg yaw mounts + lever arms; KAIST's 3 real sources are co-located (extrinsic ~= I).

The system fuses **odometry** (per-sensor SE3 increments + covariance), not raw images/radar — turning cameras/radar into odometry is an upstream VO / radar-ego-velocity job. So the test feeds 9 odometry streams. Approach (user-chosen): **mix real + synthesized on KAIST** now; real nuScenes per-sensor pipeline as a follow-up.

## Rig (9 sources, KAIST urban drive)

| id | source | origin | extrinsic | role |
|---|---|---|---|---|
| 0 | wheel | REAL KAIST | identity (ref) | reference |
| 1 | wheel+FOG | REAL KAIST | identity | heading-grade |
| 2 | wheel+IMU | REAL KAIST | identity | — |
| 3 | cam_front | SYNTH (GT-derived) | yaw 0, lever [1.5,0,1.2] | surround VO |
| 4 | cam_right | SYNTH | yaw -90, lever [0,-0.9,1.0] | surround VO |
| 5 | cam_rear | SYNTH | yaw 180, lever [-1.5,0,1.2] | surround VO |
| 6 | cam_left | SYNTH | yaw 90, lever [0,0.9,1.0] | surround VO |
| 7 | radar_fl | SYNTH | yaw -45, lever [1.6,0.7,0.5] | ego-velocity (heading-WEAK) |
| 8 | radar_fr | SYNTH | yaw 45, lever [1.6,-0.7,0.5] | ego-velocity (heading-WEAK) |

Synthesizer `tools/kaist_surround_synth.py` (`54834c5`, `8fc2602`): from the KAIST GT trajectory it forms base-frame increments, conjugates each by the planted mount extrinsic `delta_sensor = E^-1 o delta_ref o E` (the `inject_calib` convention), scales, adds seeded se3 noise, and writes increment CSVs with modeled 6-var columns (cameras heading-decent -> small rot var; radars heading-weak -> ~100x rot var). `--noise-scale` multiplies the per-step sigmas. Recipe: urban split recipe + `split_median=true` + `heading_monitor=true`, `max_sources=9`, **no** `rot_weight_prior` (the monitor discovers the grade).

## Results

### Sanity (noise-free synth + TRUE priors) — convention pin
urban12: local rot p50 **0.00165 rad (0.095 deg)**, tail 2.0 m. The noise-free GT-derived cameras fed their true priors re-align to base EXACTLY and pull the fused heading near-perfect -> the frame-align conjugation convention is correct (an inverted E would have made the 6 synth sources garbage and wrecked the median).

### Main (noisy, TRUE priors) — does 9-source fusion work?
At a **realistic** per-step noise level (`--noise-scale 0.1`: cam sigma_rot 0.015 deg/step, radar 0.15 deg/step):

| drive | gps_applied | tail trans | local rot p50 | vs 3-real-source baseline |
|---|---|---|---|---|
| urban12 | 1328/2479 | 2.21 m | 0.0086 rad (0.49 deg) | baseline tail 1.88 m / p50 0.23 deg |
| urban07 | 551/556 | 7.69 m | 0.034 rad (1.9 deg) | baseline tail 8.33 m / p50 0.93 deg |

9-source fusion is **stable and bounded**, comparable to the 3-real-source baseline. It is slightly looser because we *added 6 noisier synthetic sensors* — the robust split median + reliability absorb them gracefully instead of letting them hurt (the correct behaviour: more sensors is not automatically better when the additions are noisier than your best existing source).

**Quality auto-differentiation (the many-sensor machinery, both drives):**
- **Heading monitor** auto-discovered the heading-grade sources among 9: real FOG (src1) **boost 10**, and the clean synthetic cameras boosted too (front cam boost 10 on urban12; 8-9 on urban07); both radars **boost 1** (scores ~100x worse — correctly not boosted). No `rot_weight_prior` needed.
- **19b per-channel reliability** floored both radars to **0.2 / 0.2** (heading-weak, downweighted in both channels) and capped the cleanest camera rotation at 5, while trimming the cameras' noisier translation.

### Calibration recovery (noisy, PERTURBED priors +5 deg yaw +0.1 m lever) — multi-extrinsic headline
All 4 surround cameras, each perturbed +5 deg, recovered **simultaneously** on BOTH drives:

| cam | planted yaw | recovered (urban07) | recovered (urban12) |
|---|---|---|---|
| front (src3) | 0 deg | 0.27 deg | -0.72 deg |
| right (src4) | -90 deg | -89.75 deg | -90.2 deg |
| rear (src5) | 180 deg | -179.5 deg | 179.0 deg |
| left (src6) | 90 deg | 90.7 deg | 89.6 deg |

4 different camera mounts calibrated at once, back to <1 deg of truth, committed. Side-camera **lateral levers recovered** (y: 0.882/0.899 vs 0.9, conf ~0.99 committed); front/rear **along-motion x-levers correctly withheld at conf 0** (unobservable on planar urban — the known z/along-axis lever limitation). Radar rotation poorly recovered (heading-weak by design, conf ~0.1). Drives stayed bounded (urban12 C tail 2.26 m, urban07 C tail 7.45 m).

## Lessons / boundaries

- **Per-step noise must be physically scaled.** iid per-step rotation noise integrates to a heading random walk `sigma_step * sqrt(N_steps)`; on a 41-min 100 Hz drive (~246k steps) the first (over-pessimistic) cam 0.15 deg/step + radar 1.5 deg/step level reached tens of degrees of random-walk heading, overwhelmed the GPS gate, and death-spiralled urban12 (tail 7670 m, gps_applied 84). `--noise-scale 0.1` -> realistic VO -> bounded (tail 2.2 m, gps_applied 1328). This is a modelling lesson, not a fusion limitation.
- **Bootstrap needs a non-poisoned consensus.** At the absurd ns=1 noise, urban12 also death-spiralled the perturbed-prior run; at realistic noise it recovers on both drives. The earlier single-source `inject_calib` recovery worked because the mis-calibrated source was a MINORITY against a good consensus; if excessive noise corrupts the consensus the calibrator has nothing clean to vote against.
- The radars are a good **per-channel split-median** test: strong translation, weak heading -> reliability floors their rotation channel, the monitor keeps them off the heading podium, and their rotation extrinsic stays unobservable — all correct.

---

# Real nuScenes pipeline — 7 REAL sources (5 radar Doppler + CAN wheel + IMU)

**Done (2026-06-15)**: the follow-up to the synthesized rig, with ACTUAL per-sensor odometry instead of GT-synthesis. Converter `tools/nuscenes_to_csv.py` (`799dcb8`); data nuScenes v1.0-mini + CAN-bus expansion at `C:\workspace\data\nuScenes`; run dir `nuscenes_run/` (gitignored). Cameras (real VO) deliberately deferred — radar+CAN first cut (user-chosen).

## Rig (7 real sources, per nuScenes mini scene)

| id | source | from | odometry |
|---|---|---|---|
| 0 | wheel (ref) | CAN `vehicle_monitor` (~2 Hz) | vehicle_speed (km/h) + yaw_rate (deg/s) -> 2D odom |
| 1 | imu | CAN `ms_imu` (~100 Hz) raw gyro + wheel speed | gyro rotation (heading-grade) + forward |
| 2-6 | radar x5 | each radar `.pcd` Doppler (~13 Hz) | seeded-RANSAC ego-velocity (Kellner radial, RAW vx/vy) -> translation; R=I + (0.3 rad)^2 rot var (heading-blind) |
| GT | ego_pose | LiDAR-localization global poses | absolute track |

Each radar's real `calibrated_sensor` mount -> `prior_extrinsic` (a REAL multi-extrinsic rig: 5 radars at genuinely different positions). Raw streams only — the CAN `pose` message (pre-filtered) is NOT used (raw-odometry input contract). No GPS block (ego_pose is GT only -> predict + calibrate, no absolute correction).

## Results (scene-0061 92 m/19 s, scene-0796 237 m/20 s)

**Converter self-check (unit/convention pin):** unit auto-pick resolved vehicle_speed=km/h, yaw_rate=deg/s (the m/s / rad/s readings were 3.6x / 57x off, rejected). Each source's dead-reckoned path tracks ego_pose GT within ~10%: wheel dist-ratio 0.93/0.97, radars 0.96-1.01 — EXCEPT RADAR_FRONT_LEFT on scene-0796 (real Doppler dropout on 147/270 sweeps -> those emit zero-trans + huge cov, fusion ignores them; flagged, not masked). This validates the radar ego-velocity math (raw vx/vy radial RANSAC, not the circular ego-compensated vx_comp).

| metric | scene-0061 | scene-0796 |
|---|---|---|
| steps / fused | 977 / 969 | 987 / 979 |
| tail trans / rot | 5.40 m / 0.112 rad | 9.17 m / 0.059 rad |
| rms trans / rot | 3.88 m / 0.072 rad | 5.45 m / 0.043 rad |
| pose NEES (DOF 6) | **1.59** | **2.57** |

**Verdict: real 7-source fusion works and stays bounded** on genuine nuScenes sensor noise — all 5 radar Doppler streams + CAN wheel + CAN IMU fuse with no filter divergence; NEES 1.6/2.6 (conservative, no blow-up). Drift is pure dead-reckoning (no absolute correction by design) growing with path length (9 m over 237 m).

## Limits (honest)
- **20 s scenes -> no calibration convergence** (predicted): mounts stay near their priors; this run validates real-sensor FUSION + drift/consistency with mounts pinned, NOT extrinsic recovery. Recovery needs the longer trainval logs.
- **Heading monitor inert here**: it scores against a GPS-derived course; with no GPS block there's no anchor (`scored=no`, boost 1 all). Radar-rotation down-weighting is done instead by the split median + the per-delta (0.3 rad)^2 rot covariance (rotation taken from wheel/imu). The 19b reliability layer is quiet on 20 s scenes (only wheel rel_rot=5) — too short to warm the EMA + radar R=I is not a scatter outlier.

## Multi-scene calibration accumulation (mini, 2026-06-15)

To get past the 20 s limit WITHOUT downloading trainval: `tools/nuscenes_concat.py` (`5f68121`) concatenates all 10 mini scenes per source into one ~200 s INCREMENT megastream (re-base each scene to 0, offset by cumulative end-time + one tick; one identity seed at the start; teleport lives only in absolute GT, which calibration never touches) -> ONE persistent calibrator accumulates votes across all scenes. Verified strictly monotonic + row-exact; ~1567 deg total yaw turned / 1115 m path / ~13k RANSAC radar increments = abundant excitation.

**Result: accumulation works mechanically but the radar LEVER does NOT recover** (perturbed +0.20 m -> conf 0, never committed, in EVERY config tested). What accumulates: radar SCALE commits + tightens (1.02-1.03 on one 20 s scene -> ~1.004-1.013 over 200 s, toward 1.0). The +0.20 m lever perturbation passes straight through (one radar's perturbed z = exactly +0.20000 m in the readout).

**Verified across 3 configs** (single 20 s scene; 10-scene rot3d+joint ON; 10-scene plain Phase-2):
- 10-scene vs 1-scene: scale + the (held) rotation prior commit in BOTH -> accumulation is NOT a vote-count problem.
- rot3d+joint ON: the rot3d/joint gate keys off the SOURCE's rotation excitation, which a heading-blind radar (R=I) never has -> the lever path never opens -> conf 0 (rotation stays at the true prior).
- rot3d+joint OFF (plain Phase-2): the lever is STILL conf 0, AND Phase-1 now commits GARBAGE radar rotation (e.g. front radar rx -0.984 rad vs true ~0, at conf 0.95) -- a heading-blind source has no real rotation to calibrate, so Phase-1 fits noise. (Footgun: do NOT run rotation calibration on a translation-only source; pin its rotation prior.)

**Verdict (verified): trainval will NOT fix the radar lever -- it is a CALIBRATOR-DESIGN limit, not a data-quantity gap.** The lever-arm signal physically EXISTS in the radar's turn-velocity (`v_sensor = v_ego + omega x lever`, omega from wheel/imu), so the lever is observable IN PRINCIPLE -- but the hand-eye lever rows are coupled to a source rotation the heading-blind radar cannot provide (rot3d gate never opens; Phase-1 gives garbage R_X). Extracting it needs a PURPOSE-BUILT translation-only-source lever estimator (correlate the radar's measured velocity against the base omega over turns) -- a core feature, not a download. What accumulation DID prove: votes accumulate continuously across concatenated increment streams (a reusable short-log pattern), and radar scale converges.

---

# Real camera surround — 6-cam VO multi-extrinsic recovery (the deferred half, now REAL)

**Done (2026-06-16)**: the camera surround half, with ACTUAL monocular visual odometry (not GT-synth). The proper venue for real multi-extrinsic recovery — cameras self-calibrate rotation, unlike the heading-blind radar. Tools `tools/nuscenes_cam_vo.py` (`a39521d`) + `tools/nuscenes_surround.py` + `tools/nuscenes_surround_score.py` (`fc74b72`); opencv 4.13 (relaxed-edge, tool-only).

**1-cam de-risk (`a39521d`, GO)**: ORB+KLT + Essential-matrix pose on CAM_FRONT (scene-0061) → per-frame rotation median **0.05–0.10°**, yaw correlation **0.87–0.97** on turning scenes (low on straight = correct), translation-DIRECTION error 1–4° (monocular up-to-scale; the calibrator recovers scale). Stretch: the camera fuses as a FULL-DOF source and its rotation extrinsic IS observable (perturbed +5° → recovered to **2.08°**, conf 0.68, committed; scale committed) — the exact contrast with the radar's structural conf-0.

**6-cam surround (`fc74b72`)**: VO on all 6 cameras (each its own K + `ego_from_cam` mount) → 8-source fuse (CAN wheel ref + imu + 6 VO cams), split_median, no GPS, no `translation_only`. RECOVER run (+5° yaw perturbation per camera):

| metric | result |
|---|---|
| cameras recovering +5° → <3° yaw err | **5 of 6** (scene-0061; CAM_BACK_RIGHT the miss, worst VO conf 0.17) |
| fused NEES (CTRL) | 0.11–0.12 (bounded, conservative) |
| tail drift | 1.76–2.07 m |
| 6-cam vs 1-cam cumulative heading (high-fallback scene) | +22.1° → **+2.0°** (redundancy closes the gap) |

**Verdict: the real 6-cam surround fuses + recovers the camera ROTATION extrinsics simultaneously on real VO** — the first real-data MULTI-extrinsic rotation recovery, the venue the heading-blind radar could not be. Honest limits: (1) mount YAW recovers but optical-axis ROLL is weakly observable on planar urban (the camera analogue of the radar lever limit — yaw is the headline metric); (2) levers stay conf-0 on 20 s scenes (too short to warm the lever estimator, same as the synth surround); (3) side/back cameras recover worse than front/front-corner (less forward parallax, higher fallback) — reported per-camera, not hidden; (4) the heading monitor is inert without GPS (same as the 7-source radar run).

## Natural next steps
- **Translation-only-source lever estimator** (core feature) → the radar (velocity-only) lever unlock; trainval data alone will not (the algorithm-bound is the calibrator; the SENSOR-bound radar rotation is unlocked by a 4D radar — `RADAR_SCAN_ODOMETRY.md`).
- **Real 4D radar dataset** (View-of-Delft / K-Radar) → confirm the synth-4D radar rotation unlock on real imaging radar (download-gated; the converter + 3D path are ready).
- **Longer camera scenes (trainval)** → warm the camera LEVER estimator past the 20 s limit (the levers were conf-0 only for scene length, not observability); also closes the optical-axis-roll gap with more motion variety.
- **Noisy ego_pose as a sparse GPS-style correction** → bounds the dead-reckoning drift AND activates the heading monitor on nuScenes (mildly circular since ego_pose is GT; a deliberate stand-in).
