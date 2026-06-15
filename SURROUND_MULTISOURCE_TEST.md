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

## Follow-up (user-requested next step)
**Real nuScenes per-sensor pipeline**: 6 surround cameras + 5 radars + IMU + wheel(CAN) + GT poses, with ACTUAL per-sensor odometry (run VO on the cameras + radar ego-velocity) instead of GT-synthesized streams. Heavier (download + per-sensor odometry pipelines, no odometry-quality GT), deferred to a follow-up session. The synthesized rig here is the rehearsal: it pins the conventions, the 9-source stability, and the quality-differentiation before the real-pipeline cost.
