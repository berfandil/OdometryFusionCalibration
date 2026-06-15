#!/usr/bin/env python3
# nuscenes_to_csv.py - convert a nuScenes (v1.0-mini) scene into a 7-source many-sensor replay
# (see adapters/README.md schema). The FIRST real many-sensor fusion test: 5 radars + wheel + IMU,
# all RAW odometry, fused against ego_pose ground truth. No GPS block (ego_pose is GT only -> no
# absolute-ref correction; predict + calibrate only, like the no-GPS KAIST runs).
#
# SEVEN sources, all in INCREMENT form on their OWN sensor timeline:
#   src0 (reference) = wheel  : CAN vehicle_monitor (vehicle_speed + yaw_rate, ~2 Hz)
#   src1            = imu     : CAN ms_imu gyro (so3 exp of rotation_rate) + wheel forward speed
#   src2..6         = radar_* : per-radar Doppler ego-velocity (RANSAC over detections, ~13 Hz)
#                              RADAR_FRONT, _FRONT_LEFT, _FRONT_RIGHT, _BACK_LEFT, _BACK_RIGHT
# + GT (ego_pose -> absolute-pose CSV, first-relative to the scene start).
#
# WHY radars are heading-WEAK by construction: a single radar's Doppler observes only its 2-D
# ego-VELOCITY (vS_x, vS_y) in its sensor frame; it CANNOT observe its own rotation. So each radar
# increment carries R = I and a LARGE rotation variance, while wheel + IMU carry the heading. The
# split-median + 19b reliability + heading monitor should rank imu/wheel above the radars for
# rotation and floor the radars' rot channel.
#
# Radar mounts (calibrated_sensor ego_from_sensor) are emitted as PRIORS (prior_extrinsic), NOT
# recovered: a ~20 s mini scene is far too short for online calibration to converge -- this run
# validates real multi-sensor FUSION + drift/consistency, with the mounts pinned.
#
# Pure ASCII; uses numpy (lstsq for the RANSAC inlier refit); SEEDED RNG (random.Random(seed)) for
# the RANSAC sampling -- no time-based randomness -> reproducible.
#
# Usage:
#   python tools/nuscenes_to_csv.py <dataroot> <out_dir> <scene-name> [scene-name ...] [--seed N]
#   e.g. python tools/nuscenes_to_csv.py C:/workspace/data/nuScenes nuscenes_run scene-0061
#
# <dataroot> layout (confirmed):
#   <dataroot>/v1.0-mini/v1.0-mini/*.json          (tables)
#   <dataroot>/v1.0-mini/{samples,sweeps}/RADAR_*/  (PCD blobs)
#   <dataroot>/can_bus/can_bus/<scene>_{ms_imu,vehicle_monitor}.json
import json
import math
import os
import re
import struct
import sys

import numpy as np


# ---------- small matrix / quaternion helpers (mirror kaist_to_csv.py; quat is w-first) ----------
def mat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat_t(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def mat_vec(a, v):
    return [sum(a[i][k] * v[k] for k in range(3)) for i in range(3)]


def rot_z(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]


def quat_wxyz_to_mat(qw, qx, qy, qz):
    n = math.sqrt(qw * qw + qx * qx + qy * qy + qz * qz) or 1.0
    qw, qx, qy, qz = qw / n, qx / n, qy / n, qz / n
    return [[1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - qz * qw), 2 * (qx * qz + qy * qw)],
            [2 * (qx * qy + qz * qw), 1 - 2 * (qx * qx + qz * qz), 2 * (qy * qz - qx * qw)],
            [2 * (qx * qz - qy * qw), 2 * (qy * qz + qx * qw), 1 - 2 * (qx * qx + qy * qy)]]


def mat_to_quat(R):  # w-first
    tr = R[0][0] + R[1][1] + R[2][2]
    if tr > 0:
        s = math.sqrt(tr + 1.0) * 2.0
        w, x, y, z = 0.25 * s, (R[2][1] - R[1][2]) / s, (R[0][2] - R[2][0]) / s, (R[1][0] - R[0][1]) / s
    elif R[0][0] > R[1][1] and R[0][0] > R[2][2]:
        s = math.sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2.0
        w, x, y, z = (R[2][1] - R[1][2]) / s, 0.25 * s, (R[0][1] + R[1][0]) / s, (R[0][2] + R[2][0]) / s
    elif R[1][1] > R[2][2]:
        s = math.sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2.0
        w, x, y, z = (R[0][2] - R[2][0]) / s, (R[0][1] + R[1][0]) / s, 0.25 * s, (R[1][2] + R[2][1]) / s
    else:
        s = math.sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2.0
        w, x, y, z = (R[1][0] - R[0][1]) / s, (R[0][2] + R[2][0]) / s, (R[1][2] + R[2][1]) / s, 0.25 * s
    n = math.sqrt(w * w + x * x + y * y + z * z) or 1.0
    return (w / n, x / n, y / n, z / n)


def so3_exp(w):  # Rodrigues; w = rotation vector (rad)
    th = math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
    if th < 1e-12:
        return [[1, -w[2], w[1]], [w[2], 1, -w[0]], [-w[1], w[0], 1]]
    k = [w[0] / th, w[1] / th, w[2] / th]
    c, s = math.cos(th), math.sin(th)
    K = [[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]]
    K2 = mat_mul(K, K)
    return [[(1.0 if i == j else 0.0) + s * K[i][j] + (1 - c) * K2[i][j]
             for j in range(3)] for i in range(3)]


def quat_to_ypr(qw, qx, qy, qz):
    # Recover (yaw, pitch, roll) from a quaternion under the Rz(yaw)Ry(pitch)Rx(roll) convention
    # (config_loader parse_extrinsic). R = quat; yaw=atan2(R10,R00), pitch=asin(-R20),
    # roll=atan2(R21,R22).
    R = quat_wxyz_to_mat(qw, qx, qy, qz)
    yaw = math.atan2(R[1][0], R[0][0])
    pitch = math.asin(max(-1.0, min(1.0, -R[2][0])))
    roll = math.atan2(R[2][1], R[2][2])
    return yaw, pitch, roll


# ---------- nuScenes table + CAN bus loading ----------
RADAR_CHANNELS = ["RADAR_FRONT", "RADAR_FRONT_LEFT", "RADAR_FRONT_RIGHT",
                  "RADAR_BACK_LEFT", "RADAR_BACK_RIGHT"]
# source id assignment: wheel=0, imu=1, then the radars 2..6 in the order above.
RADAR_SRC_BASE = 2

# nuScenes timestamps are MICROSECONDS; the CsvSource / ReplayHarness clock is NANOSECONDS
# (kNanosPerSec = 1e9, and the manifest's tick_rate_hz/window_s/fusion_delay_s are seconds ->
# converted to ns internally). So every emitted t_ns column must be us*1000.
US_TO_NS = 1000


def table_dir(dataroot):
    return os.path.join(dataroot, "v1.0-mini", "v1.0-mini")


def blob_dir(dataroot):
    return os.path.join(dataroot, "v1.0-mini")


def canbus_dir(dataroot):
    return os.path.join(dataroot, "can_bus", "can_bus")


def load_json(path):
    with open(path) as f:
        return json.load(f)


def index_by(records, key):
    return {r[key]: r for r in records}


def load_tables(dataroot):
    td = table_dir(dataroot)
    return {
        "scene": load_json(os.path.join(td, "scene.json")),
        "sample": load_json(os.path.join(td, "sample.json")),
        "sample_data": load_json(os.path.join(td, "sample_data.json")),
        "calibrated_sensor": load_json(os.path.join(td, "calibrated_sensor.json")),
        "ego_pose": load_json(os.path.join(td, "ego_pose.json")),
    }


def scene_sample_tokens(tables, scene_name):
    scene = next((s for s in tables["scene"] if s["name"] == scene_name), None)
    if scene is None:
        raise ValueError("scene %s not found" % scene_name)
    by_tok = index_by(tables["sample"], "token")
    toks, t = [], scene["first_sample_token"]
    while t:
        toks.append(t)
        t = by_tok[t]["next"]
    return scene, set(toks)


def radar_stream(tables, sample_tokens, channel):
    # All sample_data rows for this radar channel within the scene (samples + sweeps = the union ->
    # dense ~13 Hz odometry), sorted by timestamp.
    pat = "/%s/" % channel
    cs = index_by(tables["calibrated_sensor"], "token")
    ep = index_by(tables["ego_pose"], "token")
    rows = []
    for sd in tables["sample_data"]:
        if pat not in sd["filename"]:
            continue
        if sd["sample_token"] not in sample_tokens:
            continue
        rows.append(sd)
    rows.sort(key=lambda r: r["timestamp"])
    if not rows:
        return None, []
    calib = cs[rows[0]["calibrated_sensor_token"]]
    stream = [(r["timestamp"], r["filename"], ep[r["ego_pose_token"]]) for r in rows]
    return calib, stream


def ego_pose_track(tables, sample_tokens):
    # Every ego_pose attached to a sample_data within the scene (any sensor), sorted+deduped by
    # timestamp. This is the GT trajectory (global_from_ego).
    ep = index_by(tables["ego_pose"], "token")
    seen = {}
    for sd in tables["sample_data"]:
        if sd["sample_token"] not in sample_tokens:
            continue
        p = ep[sd["ego_pose_token"]]
        seen[p["timestamp"]] = p
    return [seen[t] for t in sorted(seen)]


def load_can(dataroot, scene_name, kind):
    p = os.path.join(canbus_dir(dataroot), "%s_%s.json" % (scene_name, kind))
    if not os.path.exists(p):
        return []
    return load_json(p)


# ---------- radar PCD (binary little-endian; 43-byte records) ----------
# FIELDS x y z dyn_prop id rcs vx vy vx_comp vy_comp is_quality_valid ambig_state x_rms y_rms
#        invalid_state pdh0 vx_rms vy_rms ; SIZE 4 4 4 1 2 4 4 4 4 4 1 1 1 1 1 1 1 1.
PCD_FMT = "<fffBHfffffBBBBBBBB"
PCD_REC = struct.calcsize(PCD_FMT)  # 43


def read_radar_pcd(path):
    # Return a list of (x, y, vx, vy) detections in the sensor frame.
    with open(path, "rb") as f:
        data = f.read()
    idx = data.find(b"DATA binary\n")
    if idx < 0:
        return []
    header = data[:idx].decode("ascii", errors="replace")
    m = re.search(r"POINTS (\d+)", header)
    n = int(m.group(1)) if m else 0
    body = data[idx + len(b"DATA binary\n"):]
    out = []
    for i in range(n):
        off = i * PCD_REC
        if off + PCD_REC > len(body):
            break
        rec = struct.unpack_from(PCD_FMT, body, off)
        x, y = rec[0], rec[1]
        vx, vy = rec[6], rec[7]
        out.append((x, y, vx, vy))
    return out


# ---------- radar ego-velocity via seeded RANSAC ----------
# Static-world model: for detection i with bearing (cos,sin) = (x/rho, y/rho) and radial speed
# v_r = (vx*x + vy*y)/rho, we have  -v_r = vS_x*cos + vS_y*sin  (vS = sensor velocity in sensor
# frame). LINEAR in (vS_x, vS_y). RANSAC over 2-point samples; final lstsq refit on inliers.
def radar_ego_velocity(dets, rng, ransac_iters=60, thresh=0.1, min_inliers=3):
    pts = []  # (cos, sin, b)  where b = -v_r  (so b = vS.[cos,sin])
    for (x, y, vx, vy) in dets:
        rho = math.hypot(x, y)
        if rho < 0.5:
            continue
        c, s = x / rho, y / rho
        v_r = (vx * x + vy * y) / rho
        pts.append((c, s, -v_r))
    if len(pts) < min_inliers:
        return None, 0, 0.0  # too few -> low confidence
    n = len(pts)
    best_inl, best_v = [], None
    for _ in range(ransac_iters):
        i, j = rng.randrange(n), rng.randrange(n)
        if i == j:
            continue
        c0, s0, b0 = pts[i]
        c1, s1, b1 = pts[j]
        det = c0 * s1 - c1 * s0
        if abs(det) < 1e-6:
            continue
        # solve [c0 s0; c1 s1] [vx;vy] = [b0;b1]
        vx = (b0 * s1 - b1 * s0) / det
        vy = (c0 * b1 - c1 * b0) / det
        inl = [k for k in range(n)
               if abs(pts[k][0] * vx + pts[k][1] * vy - pts[k][2]) <= thresh]
        if len(inl) > len(best_inl):
            best_inl, best_v = inl, (vx, vy)
    if best_v is None or len(best_inl) < min_inliers:
        return None, len(best_inl), 0.0
    # least-squares refit on inliers (numpy lstsq)
    A = np.array([[pts[k][0], pts[k][1]] for k in best_inl], dtype=float)
    bb = np.array([pts[k][2] for k in best_inl], dtype=float)
    sol, _res, _rank, _sv = np.linalg.lstsq(A, bb, rcond=None)
    vS_x, vS_y = float(sol[0]), float(sol[1])
    resid = A.dot(sol) - bb
    rms = float(np.sqrt(np.mean(resid * resid))) if len(resid) else 0.0
    return (vS_x, vS_y), len(best_inl), rms


# ---------- CSV emit ----------
INC_HEADER = ("# form: increment\n"
              "# t_ns, x, y, z, qw, qx, qy, qz, var_x, var_y, var_z, var_rx, var_ry, var_rz\n")


def open_inc(path, comment):
    f = open(path, "w", newline="")
    f.write(INC_HEADER)
    f.write("# %s\n" % comment)
    # identity seed (increment-form contract); modest seed variances
    f.write("0,0,0,0,1,0,0,0,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4\n")
    return f


def emit_inc(f, t_ns, R, t, var6):
    qw, qx, qy, qz = mat_to_quat(R)
    f.write("%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n"
            % (t_ns, t[0], t[1], t[2], qw, qx, qy, qz,
               var6[0], var6[1], var6[2], var6[3], var6[4], var6[5]))


# ---------- dead-reckoning self-check (integrate increments, compare to GT) ----------
def path_length(xy):
    L = 0.0
    for k in range(1, len(xy)):
        L += math.hypot(xy[k][0] - xy[k - 1][0], xy[k][1] - xy[k - 1][1])
    return L


def integrate_se2(increments):
    # increments: list of (R(3x3), t(3)) body increments. Return list of global (x,y) positions.
    R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    p = [0.0, 0.0, 0.0]
    xy = [(0.0, 0.0)]
    for (Rd, td) in increments:
        p = [p[i] + mat_vec(R, td)[i] for i in range(3)]
        R = mat_mul(R, Rd)
        xy.append((p[0], p[1]))
    return xy


# ---------- per-source builders ----------
def build_wheel(vm, t0_ns, speed_scale, yaw_scale):
    # vehicle_monitor: vehicle_speed + yaw_rate at ~2 Hz. forward = v*dt along +x; yaw = yaw_rate*dt.
    # speed_scale/yaw_scale chosen by the unit self-check (km/h->m/s = 1/3.6, deg/s->rad/s = pi/180).
    vm = sorted(vm, key=lambda r: r["utime"])
    rows, increments = [], []
    for k in range(1, len(vm)):
        tp, tc = vm[k - 1]["utime"], vm[k]["utime"]
        dt = (tc - tp) * 1e-6
        if dt <= 0:
            continue
        v = vm[k]["vehicle_speed"] * speed_scale
        w = vm[k]["yaw_rate"] * yaw_scale
        d = v * dt
        dyaw = w * dt
        R = rot_z(dyaw)
        t = [d, 0.0, 0.0]
        rows.append((tc - t0_ns, R, t))
        increments.append((R, t))
    return rows, increments


def speed_at(vm_sorted, vm_ts, ts_us):
    # nearest-prior vehicle_speed (m/s already scaled outside) for IMU translation pairing.
    import bisect
    k = bisect.bisect_right(vm_ts, ts_us) - 1
    k = min(max(k, 0), len(vm_sorted) - 1)
    return vm_sorted[k]


def build_imu(imu, vm, t0_ns, speed_scale):
    # ms_imu gyro rotation_rate [wx,wy,wz] (rad/s) -> R = exp(omega*dt). Translation = wheel forward
    # speed (interpolated to imu time) * dt along +x  (a wheel+IMU pairing; heading from gyro).
    import bisect
    imu = sorted(imu, key=lambda r: r["utime"])
    vm_sorted = sorted(vm, key=lambda r: r["utime"])
    vm_ts = [r["utime"] for r in vm_sorted]
    rows, increments = [], []
    for k in range(1, len(imu)):
        tp, tc = imu[k - 1]["utime"], imu[k]["utime"]
        dt = (tc - tp) * 1e-6
        if dt <= 0:
            continue
        wx, wy, wz = imu[k]["rotation_rate"]
        R = so3_exp([wx * dt, wy * dt, wz * dt])
        j = min(max(bisect.bisect_right(vm_ts, tc) - 1, 0), len(vm_sorted) - 1)
        v = vm_sorted[j]["vehicle_speed"] * speed_scale
        t = [v * dt, 0.0, 0.0]
        rows.append((tc - t0_ns, R, t))
        increments.append((R, t))
    return rows, increments


def build_radar(stream, rng, ransac_iters, blobroot):
    # Per radar: ego-velocity per sweep, integrated over dt to a sensor-frame translation increment.
    # R = I (heading-blind). Returns (rows, increments, stats).
    rows, increments = [], []
    n_ok = n_low = 0
    prev_t = None
    for (ts_us, filename, _ego) in stream:
        path = os.path.join(blobroot, filename)
        dets = read_radar_pcd(path)
        vel, n_inl, rms = radar_ego_velocity(dets, rng, ransac_iters)
        if prev_t is None:
            prev_t = ts_us
            continue
        dt = (ts_us - prev_t) * 1e-6
        prev_t = ts_us
        if dt <= 0:
            continue
        if vel is None:
            # low-confidence sweep: zero translation, HUGE covariance so fusion ignores it.
            R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            t = [0.0, 0.0, 0.0]
            var6 = [100.0, 100.0, 100.0, (0.5) ** 2, (0.5) ** 2, (0.5) ** 2]
            n_low += 1
        else:
            vS_x, vS_y = vel
            t = [vS_x * dt, vS_y * dt, 0.0]
            R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]  # radar cannot observe its own rotation
            # trans var from the RANSAC residual scaled by dt (floor 2 mm); LARGE rot var so the
            # split median + 19b reliability down-weight radar rotation.
            s_t = (rms * dt + 2e-3) ** 2
            s_r = (0.3) ** 2  # (0.3 rad)^2 heading-weak
            var6 = [s_t, s_t, 100.0, s_r, s_r, s_r]  # z translation unobserved -> large
            n_ok += 1
        rows.append((ts_us, R, t, var6))
        increments.append((R, t))
    return rows, increments, {"ok": n_ok, "low": n_low}


# ---------- main convert ----------
def convert(dataroot, out_dir, scene_name, tables, seed=12345, ransac_iters=60):
    scene, sample_tokens = scene_sample_tokens(tables, scene_name)
    gt = ego_pose_track(tables, sample_tokens)
    if len(gt) < 2:
        raise ValueError("scene %s has < 2 ego_poses" % scene_name)
    t0_ns_us = gt[0]["timestamp"]  # microseconds
    t0_ns = t0_ns_us  # we keep timestamps in microseconds-as-int relative to scene start

    # --- GT track (global_from_ego), first-relative to scene start (t0-ego frame) ---
    R0 = quat_wxyz_to_mat(*gt[0]["rotation"])
    R0t = mat_t(R0)
    p0 = gt[0]["translation"]
    gt_xy = []
    gt_R = []
    gt_rows = []
    for p in gt:
        Ri = quat_wxyz_to_mat(*p["rotation"])
        ti = p["translation"]
        d = [ti[0] - p0[0], ti[1] - p0[1], ti[2] - p0[2]]
        trel = mat_vec(R0t, d)
        Rrel = mat_mul(R0t, Ri)
        qw, qx, qy, qz = mat_to_quat(Rrel)
        gt_rows.append((p["timestamp"] - t0_ns, trel, (qw, qx, qy, qz)))
        gt_xy.append((trel[0], trel[1]))
        gt_R.append(Rrel)
    gt_len = path_length(gt_xy)
    # GT heading-vs-time for rotating radar increments in the self-check
    gt_ts = [p["timestamp"] for p in gt]

    # --- CAN sources ---
    vm = load_can(dataroot, scene_name, "vehicle_monitor")
    imu = load_can(dataroot, scene_name, "ms_imu")
    if not vm:
        raise ValueError("no vehicle_monitor CAN for %s" % scene_name)
    if not imu:
        raise ValueError("no ms_imu CAN for %s" % scene_name)

    # === UNIT self-check for the wheel: try km/h vs m/s and deg/s vs rad/s, pick the combo whose
    #     dead-reckoned path best matches GT distance + heading. ===
    speed_choices = [("km/h->m/s", 1.0 / 3.6), ("m/s", 1.0)]
    yaw_choices = [("deg/s->rad/s", math.pi / 180.0), ("rad/s", 1.0)]
    best = None
    for (sname, sval) in speed_choices:
        for (yname, yval) in yaw_choices:
            _rows, incs = build_wheel(vm, t0_ns, sval, yval)
            xy = integrate_se2(incs)
            L = path_length(xy)
            dist_ratio = L / gt_len if gt_len > 1e-6 else 0.0
            # heading check: net yaw turned vs GT net yaw (atan2 of final R relative)
            est_yaw = 0.0
            for (Rd, _t) in incs:
                est_yaw += math.atan2(Rd[1][0], Rd[0][0])
            gt_net_yaw = math.atan2(gt_R[-1][1][0], gt_R[-1][0][0])
            yaw_err = abs(((est_yaw - gt_net_yaw + math.pi) % (2 * math.pi)) - math.pi)
            # score: distance ratio closeness to 1 dominates; heading as a tiebreak
            score = abs(dist_ratio - 1.0) + 0.2 * yaw_err
            cand = (score, sname, sval, yname, yval, dist_ratio, yaw_err)
            if best is None or cand[0] < best[0]:
                best = cand
    _, sname, speed_scale, yname, yaw_scale, wheel_ratio, wheel_yaw_err = best
    print("  [unit-check] wheel speed=%s yaw=%s -> dist_ratio=%.3f yaw_err=%.3f rad"
          % (sname, yname, wheel_ratio, wheel_yaw_err))

    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, scene_name)

    # --- emit wheel (src0) ---
    wheel_rows, wheel_incs = build_wheel(vm, t0_ns, speed_scale, yaw_scale)
    fw = open_inc(base + "_wheel.csv", "wheel src0: vehicle_monitor speed+yaw_rate (%s, %s)" % (sname, yname))
    # modest 6-var: forward trans good, lateral/vert weak, yaw decent, roll/pitch unobserved
    wheel_var = [4e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2]
    for (t, R, tt) in wheel_rows:
        emit_inc(fw, t * US_TO_NS, R, tt, wheel_var)
    fw.close()

    # --- emit imu (src1) ---
    imu_rows, imu_incs = build_imu(imu, vm, t0_ns, speed_scale)
    fi = open_inc(base + "_imu.csv", "imu src1: ms_imu gyro so3 + wheel forward speed")
    imu_var = [4e-4, 1e-2, 1e-2, 1e-5, 1e-5, 1e-5]  # heading-grade rot
    for (t, R, tt) in imu_rows:
        emit_inc(fi, t * US_TO_NS, R, tt, imu_var)
    fi.close()

    # --- emit radars (src2..6) ---
    import random
    blobroot = blob_dir(dataroot)
    radar_stats = {}
    radar_ratios = {}
    radar_priors = {}
    for ri, channel in enumerate(RADAR_CHANNELS):
        calib, stream = radar_stream(tables, sample_tokens, channel)
        sid = RADAR_SRC_BASE + ri
        short = channel.lower().replace("radar_", "radar_")  # radar_front, radar_front_left, ...
        outp = base + "_" + short + ".csv"
        if not stream:
            print("  [warn] %s: no sweeps in scene; skipping" % channel)
            continue
        rng = random.Random(seed * 1000003 + sid)
        rows, incs, stats = build_radar(stream, rng, ransac_iters, blobroot)
        radar_stats[channel] = stats
        # prior_extrinsic = calibrated_sensor (ego_from_sensor) as yaw pitch roll x y z (radians)
        qw, qx, qy, qz = calib["rotation"]
        yaw, pitch, roll = quat_to_ypr(qw, qx, qy, qz)
        lx, ly, lz = calib["translation"]
        radar_priors[channel] = (sid, short, yaw, pitch, roll, lx, ly, lz)
        f = open_inc(outp, "%s src%d: Doppler ego-velocity (RANSAC); R=I heading-blind" % (channel, sid))
        # rows here carry their own per-row var6 (4-tuple form)
        for (t_us, R, tt, var6) in rows:
            emit_inc(f, (t_us - t0_ns) * US_TO_NS, R, tt, var6)
        f.close()
        # self-check: integrate radar translation increments rotated by GT heading at each sweep,
        # compare path LENGTH to GT (radars are heading-blind so we compare length, the robust check).
        # Build a GT-heading-rotated path: for each increment, rotate sensor t by ego_from_sensor R
        # then by the GT ego heading at that time.
        import bisect
        ER = quat_wxyz_to_mat(qw, qx, qy, qz)
        # integrate using GT ego rotation at each sweep time
        p = [0.0, 0.0, 0.0]
        xy = [(0.0, 0.0)]
        sweep_ts = [s[0] for s in stream][1:]  # increments start at 2nd sweep
        for (idx, (Rd, td)) in enumerate(incs):
            ts_us = sweep_ts[idx] if idx < len(sweep_ts) else stream[-1][0]
            gi = min(max(bisect.bisect_right(gt_ts, ts_us) - 1, 0), len(gt_R) - 1)
            # body increment in ego frame = ER * t_sensor; then rotate to global by GT ego R
            t_ego = mat_vec(ER, td)
            t_glob = mat_vec(gt_R[gi], t_ego)
            p = [p[i] + t_glob[i] for i in range(3)]
            xy.append((p[0], p[1]))
        L = path_length(xy)
        radar_ratios[channel] = (L / gt_len if gt_len > 1e-6 else 0.0, stats)

    # --- GT CSV ---
    with open(base + "_gt.csv", "w", newline="") as fgt:
        fgt.write("# t_ns, x, y, z, qw, qx, qy, qz  (ego_pose GT, first-relative)\n")
        for (t, trel, q) in gt_rows:
            fgt.write("%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f,%.9f\n"
                      % (t * US_TO_NS, trel[0], trel[1], trel[2], q[0], q[1], q[2], q[3]))

    # --- manifest ---
    write_manifest(out_dir, scene_name, radar_priors)

    dur_s = (gt[-1]["timestamp"] - gt[0]["timestamp"]) * 1e-6
    return {
        "scene": scene_name, "samples": len(sample_tokens), "gt_n": len(gt), "gt_len_m": gt_len,
        "dur_s": dur_s, "wheel_n": len(wheel_rows), "imu_n": len(imu_rows),
        "wheel_ratio": wheel_ratio, "wheel_yaw_err": wheel_yaw_err,
        "radar_ratios": radar_ratios, "radar_stats": radar_stats,
        "speed_unit": sname, "yaw_unit": yname,
    }


def write_manifest(out_dir, scene_name, radar_priors):
    lines = []
    lines.append("[global]")
    lines.append("tick_rate_hz = 50")
    lines.append("fusion_delay_s = 0.05")
    lines.append("window_s = 0.1")
    lines.append("max_sources = 7")
    lines.append("reference_sensor_id = 0")
    lines.append("cold_start = median_from_start")
    lines.append("split_median = true")
    lines.append("heading_monitor = true")
    lines.append("adaptive_q = true")
    lines.append("q_scale = 0.7")
    lines.append("q_floor = 0.01 0.01 0.01 0.001 0.001 0.001")
    lines.append("")
    lines.append("[sensor.0]")
    lines.append("id = 0")
    lines.append("is_reference = true")
    lines.append("csv = %s_wheel.csv" % scene_name)
    lines.append("form = increment")
    lines.append("")
    lines.append("[sensor.1]")
    lines.append("id = 1")
    lines.append("csv = %s_imu.csv" % scene_name)
    lines.append("form = increment")
    lines.append("")
    # radars in source-id order
    for channel in RADAR_CHANNELS:
        if channel not in radar_priors:
            continue
        sid, short, yaw, pitch, roll, lx, ly, lz = radar_priors[channel]
        lines.append("[sensor.%d]" % sid)
        lines.append("id = %d" % sid)
        lines.append("csv = %s_%s.csv" % (scene_name, short))
        lines.append("form = increment")
        lines.append("prior_extrinsic = %.9f %.9f %.9f %.6f %.6f %.6f"
                     % (yaw, pitch, roll, lx, ly, lz))
        lines.append("")
    lines.append("[gt]")
    lines.append("csv = %s_gt.csv" % scene_name)
    lines.append("")
    lines.append("[replay]")
    lines.append("warmup_steps = 20")
    lines.append("local_batch_len = 500")
    lines.append("out = %s_out.csv" % scene_name)
    lines.append("")
    with open(os.path.join(out_dir, scene_name + ".ini"), "w", newline="") as f:
        f.write("\n".join(lines))


def main():
    args = sys.argv[1:]
    seed = 12345
    pos = []
    i = 0
    while i < len(args):
        if args[i] == "--seed":
            i += 1
            seed = int(args[i])
        else:
            pos.append(args[i])
        i += 1
    if len(pos) < 3:
        print("usage: nuscenes_to_csv.py <dataroot> <out_dir> <scene-name> [scene-name ...] [--seed N]")
        sys.exit(2)
    dataroot, out_dir = pos[0], pos[1]
    scenes = pos[2:]
    print("loading nuScenes tables...")
    tables = load_tables(dataroot)
    for scene_name in scenes:
        print("=== %s ===" % scene_name)
        try:
            info = convert(dataroot, out_dir, scene_name, tables, seed=seed)
        except Exception as e:  # noqa: BLE001
            import traceback
            print("ERR %s: %s" % (scene_name, e))
            traceback.print_exc()
            continue
        print("OK  %s  samples=%d gt_n=%d dur=%.1fs gt_len=%.1f m wheel=%d imu=%d"
              % (info["scene"], info["samples"], info["gt_n"], info["dur_s"],
                 info["gt_len_m"], info["wheel_n"], info["imu_n"]))
        print("    units: speed=%s yaw=%s ; wheel dist_ratio=%.3f yaw_err=%.3f rad"
              % (info["speed_unit"], info["yaw_unit"], info["wheel_ratio"], info["wheel_yaw_err"]))
        for ch, (ratio, stats) in info["radar_ratios"].items():
            print("    %-18s dist_ratio=%.3f  ransac_ok=%d low=%d"
                  % (ch, ratio, stats["ok"], stats["low"]))


if __name__ == "__main__":
    main()
