#!/usr/bin/env python3
# nuscenes_cam_vo.py - MONOCULAR visual odometry on a single nuScenes camera (CAM_FRONT) for one
# scene, the de-risking step BEFORE building the 6-camera surround. Proves whether ONE camera's
# monocular VO produces USABLE odometry (rotation + translation DIRECTION) vs ego_pose GT on real
# nuScenes driving, and emits the project's increment-form CSV so the camera can be fed to the
# calibrator as a real FULL-6DOF rotation+translation source (unlike the heading-blind radar).
#
# WHY a single CAM_FRONT first: the surround build is 6 cameras x this same pipeline. If monocular
# VO is too fragile on nuScenes (dynamic objects + 12 Hz fast urban motion + low texture), the
# surround is not worth building on this data without a heavier pipeline. A clean GO/NO-GO with
# numbers here de-risks that large build.
#
# THE PIPELINE (per consecutive CAM_FRONT frame pair, samples+sweeps union, time-ordered):
#   1. Load both frames grayscale. Detect ORB features in the previous frame; track them into the
#      current frame with KLT (pyramidal Lucas-Kanade) optical flow -- KLT gives sub-pixel, well-
#      distributed correspondences and is more robust on consecutive ~100 ms driving frames than
#      ORB descriptor matching (which struggles with repetitive urban texture). Forward-backward
#      flow check rejects unstable tracks.
#   2. RANSAC-filter the correspondences with the ESSENTIAL matrix (cv2.findEssentialMat with the
#      CAM_FRONT intrinsic K), then cv2.recoverPose -> R (cam-to-cam rotation) + t_unit (UNIT
#      translation DIRECTION; monocular scale is ambiguous). Reject pairs with too few inliers
#      (emit identity + huge covariance, count the fallback).
#   3. Map the OpenCV camera-frame relative pose (x-right, y-down, z-forward) to the project body
#      frame. recoverPose returns (R, t) with X_prev ~= R @ X_curr + t in camera coords, i.e. the
#      camera-frame body delta over [t_prev, t_curr]. Conjugate by the CAM_FRONT mount rotation
#      X_R = R_ego_from_cam (calibrated_sensor) to express the increment in the ego/body convention:
#      R_body = X_R @ R_cam @ X_R^T,  t_body_dir = X_R @ t_cam_unit.
#   4. Emit increment CSV (t_ns,x,y,z,qw..qz,6var). Translation = t_body_dir scaled by a NOMINAL
#      placeholder (the GT speed at that time, for the validation only) -- the REAL per-source scale
#      is the calibrator's Phase-1 job (a monocular camera is a rotation + translation-DIRECTION
#      source; scale is recovered online). Mark high rot/trans covariance where inliers are low.
#
# REUSE, don't reimplement: imports tools/nuscenes_to_csv.py for the table joins, ego_pose GT
# loading, the v1.0-mini layout, the SE3/quat/SO3 helpers, the increment-CSV emit. Only the VO
# front-end + the camera-vs-GT validation is new here.
#
# VALIDATION (the GO/NO-GO -- NUMBERS):
#   1. Per-frame recovered R vs GT camera motion B_cam_gt = X_cam^-1 o A_gt o X_cam: rotation error
#      (deg, median+p90), the per-step YAW CORRELATION recovered-vs-GT (the headline: does monocular
#      VO recover usable rotation on real driving?), and cumulative heading vs GT over the scene.
#   2. Translation DIRECTION: angle between the recovered t_unit and the GT translation direction
#      (deg, median). Magnitude is scale-ambiguous -> NOT scored.
#   3. Robustness: avg ORB features, avg tracked, avg E-matrix inliers, % frames that fell back to
#      identity (low-texture / dynamic / blur failures).
#
# Pure ASCII; numpy + opencv; SEEDED (cv2.setRNGSeed + numpy seed) -> reproducible RANSAC.
#
# Usage:
#   python tools/nuscenes_cam_vo.py <dataroot> <out_dir> <scene> [scene ...]
#       [--max-features N] [--min-inliers N] [--ransac-prob P] [--ransac-thresh PX]
#       [--fb-thresh PX] [--seed N] [--quiet]
#   e.g. python tools/nuscenes_cam_vo.py C:/workspace/data/nuScenes nuscenes_run scene-0061
import math
import os
import sys

import numpy as np

import nuscenes_to_csv as n2c

try:
    import cv2
except ImportError:
    print("ERROR: opencv not installed. Run: pip install opencv-python-headless")
    sys.exit(2)

CAM_CHANNEL = "CAM_FRONT"


# ---------- CAM_FRONT sample_data stream + intrinsics/mount (mirrors radar_stream) ----------
def cam_stream(tables, sample_tokens, channel):
    # All sample_data rows for this camera channel within the scene (samples + sweeps union ->
    # dense ~10-12 Hz frames), time-ordered. Returns (calibrated_sensor, [(ts_us, filename), ...]).
    pat = "/%s/" % channel
    cs = n2c.index_by(tables["calibrated_sensor"], "token")
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
    stream = [(r["timestamp"], r["filename"]) for r in rows]
    return calib, stream


# ---------- VO front-end: ORB detect + KLT track + Essential-matrix pose ----------
def detect_track(prev_gray, curr_gray, orb, fb_thresh):
    # Detect ORB keypoints in prev, track into curr with pyramidal KLT, forward-backward checked.
    # Returns (pts_prev, pts_curr) Nx2 float32 arrays of surviving correspondences, plus n_detected.
    kps = orb.detect(prev_gray, None)
    if not kps:
        return None, None, 0
    p0 = np.array([kp.pt for kp in kps], dtype=np.float32).reshape(-1, 1, 2)
    n_det = len(kps)
    lk = dict(winSize=(21, 21), maxLevel=3,
              criteria=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 30, 0.01))
    p1, st1, _ = cv2.calcOpticalFlowPyrLK(prev_gray, curr_gray, p0, None, **lk)
    if p1 is None:
        return None, None, n_det
    # forward-backward consistency check: track curr->prev, keep points that round-trip within thresh
    p0r, st0, _ = cv2.calcOpticalFlowPyrLK(curr_gray, prev_gray, p1, None, **lk)
    if p0r is None:
        return None, None, n_det
    st = (st1.reshape(-1) == 1) & (st0.reshape(-1) == 1)
    fb = np.linalg.norm((p0 - p0r).reshape(-1, 2), axis=1)
    good = st & (fb < fb_thresh)
    pts_prev = p0.reshape(-1, 2)[good]
    pts_curr = p1.reshape(-1, 2)[good]
    return pts_prev, pts_curr, n_det


def recover_pose(pts_prev, pts_curr, K, prob, thresh, min_inliers):
    # Essential matrix (RANSAC) + recoverPose. Returns (R, t_unit, n_inliers) where the pose
    # satisfies X_prev ~= R @ X_curr + t in CAMERA coordinates (recoverPose maps the SECOND/curr
    # camera points back into the FIRST/prev frame), i.e. the camera body delta prev->curr. t is a
    # UNIT vector (monocular scale ambiguity). None if too few points/inliers.
    if pts_prev is None or len(pts_prev) < min_inliers:
        return None, None, 0
    E, mask = cv2.findEssentialMat(pts_curr, pts_prev, K, method=cv2.RANSAC,
                                   prob=prob, threshold=thresh)
    if E is None or E.shape != (3, 3):
        return None, None, 0
    # recoverPose: points1=curr, points2=prev -> returns R,t with prev = R*curr + t (cheirality-checked)
    n_in, R, t, mask2 = cv2.recoverPose(E, pts_curr, pts_prev, K, mask=mask)
    if n_in < min_inliers:
        return None, None, int(n_in)
    tu = t.reshape(3)
    nrm = np.linalg.norm(tu)
    if nrm < 1e-9:
        return None, None, int(n_in)
    return R, tu / nrm, int(n_in)


# ---------- camera-frame pose -> project body-frame increment ----------
def cam_to_body(R_cam, t_cam_unit, X_R):
    # The recovered (R_cam, t_cam) is the camera body delta prev->curr in OpenCV camera coords
    # (x-right, y-down, z-forward). The project body frame is the ego frame (x-fwd, y-left, z-up).
    # X_R = R_ego_from_cam (calibrated_sensor rotation) maps a camera-frame vector to the ego frame.
    # A body delta transforms by conjugation: R_body = X_R @ R_cam @ X_R^T; a direction by X_R @ t.
    Xt = X_R.T
    R_body = X_R @ R_cam @ Xt
    t_body = X_R @ t_cam_unit
    return R_body, t_body


# ---------- GT camera frame-to-frame motion (mirrors the radar-GT delta) ----------
def gt_cam_delta(gt_ts, gt_R, gt_p, X_R3, X_t3, t_prev_us, t_curr_us):
    # B_cam_gt = X_cam^-1 o A_gt o X_cam, the GT camera body delta prev->curr in the CAMERA frame.
    # A_gt is the GT ego body delta over [t_prev, t_curr] (global poses composed; the X^-1 A X
    # cancels the global frame). Returns (R_gt 3x3 list, t_gt 3 list) in the camera frame.
    import bisect
    gi0 = min(max(bisect.bisect_right(gt_ts, t_prev_us) - 1, 0), len(gt_R) - 1)
    gi1 = min(max(bisect.bisect_right(gt_ts, t_curr_us) - 1, 0), len(gt_R) - 1)
    R0g, p0g = gt_R[gi0], gt_p[gi0]
    R1g, p1g = gt_R[gi1], gt_p[gi1]
    R0t = n2c.mat_t(R0g)
    Agt_R = n2c.mat_mul(R0t, R1g)
    Agt_t = n2c.mat_vec(R0t, [p1g[i] - p0g[i] for i in range(3)])
    # B_cam_gt = X^-1 o A o X  with X = (X_R, X_t) ego_from_cam
    XRt = n2c.mat_t(X_R3)
    BR = n2c.mat_mul(XRt, n2c.mat_mul(Agt_R, X_R3))
    inner = [n2c.mat_vec(Agt_R, X_t3)[i] + Agt_t[i] - X_t3[i] for i in range(3)]
    Bt = n2c.mat_vec(XRt, inner)
    return BR, Bt


# ---------- validation metrics ----------
def rot_err_deg(Re, Rg):
    Rt = n2c.mat_mul(n2c.mat_t(Re), Rg)
    tr = Rt[0][0] + Rt[1][1] + Rt[2][2]
    c = max(-1.0, min(1.0, (tr - 1.0) / 2.0))
    return math.degrees(math.acos(c))


def yaw_of(R):
    return math.atan2(R[1][0], R[0][0])


def angle_between_deg(a, b):
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(x * x for x in b))
    if na < 1e-9 or nb < 1e-9:
        return float("nan")
    c = sum(a[i] * b[i] for i in range(3)) / (na * nb)
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def pct(vals, q):
    if not vals:
        return float("nan")
    s = sorted(vals)
    k = min(len(s) - 1, int(round(q * (len(s) - 1))))
    return s[k]


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = math.sqrt(sxx * syy)
    if den < 1e-12:
        return float("nan")
    return sxy / den


# ---------- per-scene VO ----------
def vo_scene(stream, calib, blobroot, gt_ts, gt_R, gt_p, vm, speed_scale, args):
    # Run monocular VO over the CAM_FRONT stream. Returns (rows, stats, perframe). rows are
    # increment-CSV tuples (t_us, R3, t3, var6); perframe holds (R_body, t_body, Rg_cam, tg_cam,
    # recovered, ts_us, t_prev_us) for the vs-GT validation.
    import bisect
    X_R = np.array(n2c.quat_wxyz_to_mat(*calib["rotation"]))   # R_ego_from_cam, numpy for the VO math
    X_R3 = X_R.tolist()                                        # list form for the n2c helpers
    X_t3 = list(calib["translation"])

    # CAN speed for the NOMINAL translation placeholder (validation only; real scale is the
    # calibrator's). nearest-prior vehicle_speed (SI) held over the frame.
    vm_sorted = sorted(vm, key=lambda r: r["utime"]) if vm else []
    vm_ts = [r["utime"] for r in vm_sorted]

    orb = cv2.ORB_create(nfeatures=args.max_features, fastThreshold=12)

    rows = []
    perframe = []
    stats = {"n_pairs": 0, "n_fallback": 0, "sum_det": 0, "sum_track": 0, "sum_inl": 0,
             "n_recov": 0}

    prev_gray = None
    prev_ts = None
    for (ts_us, filename) in stream:
        path = os.path.join(blobroot, filename)
        gray = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        if gray is None:
            continue
        if prev_gray is None:
            prev_gray, prev_ts = gray, ts_us
            continue

        dt = (ts_us - prev_ts) * 1e-6
        if dt <= 0:
            prev_gray, prev_ts = gray, ts_us
            continue

        pts_prev, pts_curr, n_det = detect_track(prev_gray, gray, orb, args.fb_thresh)
        n_track = 0 if pts_prev is None else len(pts_prev)
        R_cam, t_cam, n_inl = recover_pose(pts_prev, pts_curr, args.K,
                                           args.ransac_prob, args.ransac_thresh, args.min_inliers)
        stats["n_pairs"] += 1
        stats["sum_det"] += n_det
        stats["sum_track"] += n_track
        stats["sum_inl"] += n_inl

        recovered = R_cam is not None
        if recovered:
            R_body, t_body = cam_to_body(R_cam, t_cam, X_R)
            R_body_l = R_body.tolist()
            t_dir = t_body.tolist()
            # the EMITTED increment frame: body (default; the mount already applied -> ego-aligned,
            # extrinsic ~= I) OR raw OpenCV camera optical frame (--raw-cam-frame; leaves the real
            # ~90 deg mount in, so the calibrator has a genuine camera extrinsic to RECOVER -- the
            # stretch calibration test). The vs-GT validation always uses the body comparison below.
            if args.raw_cam_frame:
                R_emit_l = R_cam.tolist()
                t_emit_dir = t_cam.tolist()
            else:
                R_emit_l = R_body_l
                t_emit_dir = t_dir
            # nominal translation placeholder = GT speed * dt (validation only; scale is calibrator's)
            v = 0.0
            if vm_sorted:
                j = min(max(bisect.bisect_right(vm_ts, ts_us) - 1, 0), len(vm_sorted) - 1)
                v = vm_sorted[j]["vehicle_speed"] * speed_scale
            t_emit = [t_emit_dir[i] * v * dt for i in range(3)]
            # confidence: rot var tightens with inliers, translation modest (direction-only)
            s_r = max(0.003, 0.5 / max(1, n_inl)) ** 2     # rad^2; ~0.5/inliers floor 0.003
            s_t = (0.05) ** 2
            var6 = [s_t, s_t, s_t, s_r, s_r, s_r]
            stats["n_recov"] += 1
            R_out, t_out = R_emit_l, t_emit
        else:
            R_body_l = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            t_dir = [0.0, 0.0, 0.0]
            R_out = R_body_l
            t_out = [0.0, 0.0, 0.0]
            var6 = [100.0, 100.0, 100.0, (0.5) ** 2, (0.5) ** 2, (0.5) ** 2]
            stats["n_fallback"] += 1

        rows.append((ts_us, R_out, t_out, var6))

        # GT camera frame-to-frame motion over [t_prev, t_curr], in the CAMERA frame
        Rg_cam, tg_cam = gt_cam_delta(gt_ts, gt_R, gt_p, X_R3, X_t3, prev_ts, ts_us)
        # validation ALWAYS compares in the BODY frame: R_body_l / t_dir are the body-frame recovered
        # increment regardless of the emit mode, and the GT cam delta is mapped to the body frame.
        Rg_body = (X_R @ np.array(Rg_cam) @ X_R.T).tolist()
        tg_body = (X_R @ np.array(tg_cam)).tolist()
        perframe.append((R_body_l, t_dir, Rg_body, tg_body, recovered, ts_us, prev_ts))

        prev_gray, prev_ts = gray, ts_us

    return rows, stats, perframe


def report_cam(channel, perframe, stats):
    # Per-frame rotation error (deg) + translation-DIRECTION error (deg) recovered-vs-GT over
    # recovered frames; the headline per-step yaw correlation; cumulative + per-step heading track.
    rot_errs = []
    dir_errs = []
    est_yaws, gt_yaws = [], []
    cum_est_yaw = 0.0
    cum_gt_yaw = 0.0
    for rec in perframe:
        (R_inc, t_dir, Rg, tg, recovered, ts_us, t_prev_us) = rec
        cum_est_yaw += yaw_of(R_inc)
        cum_gt_yaw += yaw_of(Rg)
        if not recovered:
            continue
        rot_errs.append(rot_err_deg(R_inc, Rg))
        # translation DIRECTION error: recovered unit dir vs GT body-delta direction
        dir_errs.append(angle_between_deg(t_dir, tg))
        est_yaws.append(math.degrees(yaw_of(R_inc)))
        gt_yaws.append(math.degrees(yaw_of(Rg)))
    n = stats["n_pairs"]
    return {
        "channel": channel,
        "n_pairs": n,
        "n_recov": stats["n_recov"],
        "fallback_pct": 100.0 * stats["n_fallback"] / n if n else 0.0,
        "avg_det": stats["sum_det"] / n if n else 0.0,
        "avg_track": stats["sum_track"] / n if n else 0.0,
        "avg_inl": stats["sum_inl"] / n if n else 0.0,
        "rot_med": pct(rot_errs, 0.5),
        "rot_p90": pct(rot_errs, 0.9),
        "dir_med": pct(dir_errs, 0.5),
        "dir_p90": pct(dir_errs, 0.9),
        "yaw_corr": pearson(est_yaws, gt_yaws),
        "cum_est_yaw_deg": math.degrees(cum_est_yaw),
        "cum_gt_yaw_deg": math.degrees(cum_gt_yaw),
        "rot_errs": rot_errs,
    }


# ---------- CSV emit (mirrors radar_scan_odometry) ----------
def emit_vo_csv(path, comment, rows, t0_us):
    f = n2c.open_inc(path, comment)
    for (t_us, R, t, var6) in rows:
        n2c.emit_inc(f, (t_us - t0_us) * n2c.US_TO_NS, R, t, var6)
    f.close()


def scene_gt_and_speed(dataroot, scene_name, tables):
    # Shared per-scene setup: GT ego_pose track (timestamps + R + p) and the CAN-derived speed
    # scale for the nominal translation placeholder. Returns (gt, t0_us, gt_ts, gt_R, gt_p, vm,
    # speed_scale). REUSED by the single-camera convert_scene AND the 6-camera surround driver.
    scene, sample_tokens = n2c.scene_sample_tokens(tables, scene_name)
    gt = n2c.ego_pose_track(tables, sample_tokens)
    if len(gt) < 2:
        raise ValueError("scene %s has < 2 ego_poses" % scene_name)
    t0_us = gt[0]["timestamp"]
    gt_ts = [p["timestamp"] for p in gt]
    gt_R = [n2c.quat_wxyz_to_mat(*p["rotation"]) for p in gt]
    gt_p = [list(p["translation"]) for p in gt]

    vm = n2c.load_can(dataroot, scene_name, "vehicle_monitor")
    speed_scale = 1.0 / 3.6
    if vm:
        try:
            from radar_scan_odometry import wheel_unit_scales
            speed_scale, _yaw_scale = wheel_unit_scales(vm, gt_R)
        except Exception:  # noqa: BLE001
            speed_scale = 1.0 / 3.6
    return sample_tokens, gt, t0_us, gt_ts, gt_R, gt_p, vm, speed_scale


def run_channel_vo(channel, dataroot, scene_name, tables, sample_tokens, t0_us,
                   gt_ts, gt_R, gt_p, vm, speed_scale, args):
    # Run the proven monocular VO on ONE camera channel. Returns (calib, rows, stats, perframe)
    # where rows are increment-CSV tuples in the EMITTED frame (raw camera optical if
    # args.raw_cam_frame, else ego/body) and perframe is the body-frame vs-GT validation. The
    # camera's own intrinsic K + ego_from_cam mount come from its calibrated_sensor. This is the
    # per-camera core the 6-camera surround driver loops over (NO file IO -- the caller emits).
    calib, stream = cam_stream(tables, sample_tokens, channel)
    if not stream:
        return None, None, None, None
    args.K = np.array(calib["camera_intrinsic"], dtype=float)
    blobroot = n2c.blob_dir(dataroot)
    rows, stats, perframe = vo_scene(stream, calib, blobroot, gt_ts, gt_R, gt_p,
                                     vm, speed_scale, args)
    return calib, stream, (rows, stats, perframe)


def convert_scene(dataroot, out_dir, scene_name, tables, args):
    (sample_tokens, gt, t0_us, gt_ts, gt_R, gt_p,
     vm, speed_scale) = scene_gt_and_speed(dataroot, scene_name, tables)

    calib, stream, vo = run_channel_vo(CAM_CHANNEL, dataroot, scene_name, tables, sample_tokens,
                                       t0_us, gt_ts, gt_R, gt_p, vm, speed_scale, args)
    if stream is None:
        raise ValueError("no %s frames in scene %s" % (CAM_CHANNEL, scene_name))
    rows, stats, perframe = vo

    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, scene_name)
    suf = "_camframe" if args.raw_cam_frame else ""
    outp = base + "_camvo_front" + suf + ".csv"
    frame_note = "raw OpenCV camera optical frame (mount left in for calibration)" \
        if args.raw_cam_frame else "ego/body frame (mount applied)"
    emit_vo_csv(outp, "%s monocular VO: ORB+KLT+Essential (R + t-DIRECTION); %s; scale=calibrator"
                % (CAM_CHANNEL, frame_note), rows, t0_us)
    rep = report_cam(CAM_CHANNEL, perframe, stats)
    rep["n_frames"] = len(stream)
    rep["dur_s"] = (stream[-1][0] - stream[0][0]) * 1e-6
    return rep, scene_name


def print_report(scene_name, rep, args):
    print("=" * 100)
    print("MONOCULAR CAMERA VO -- %s  (%s, %d frames, %.1f s)  "
          "[max_feat=%d min_inl=%d ransac_thr=%.2fpx prob=%.4f]"
          % (scene_name, rep["channel"], rep["n_frames"], rep["dur_s"],
             args.max_features, args.min_inliers, args.ransac_thresh, args.ransac_prob))
    print("=" * 100)
    print("ROBUSTNESS  (per-frame averages over %d frame pairs)" % rep["n_pairs"])
    print("  ORB detected/frame : %7.1f" % rep["avg_det"])
    print("  KLT tracked/frame  : %7.1f" % rep["avg_track"])
    print("  E-matrix inliers   : %7.1f" % rep["avg_inl"])
    print("  recovered frames   : %7d / %d" % (rep["n_recov"], rep["n_pairs"]))
    print("  fallback (identity): %7.1f%%   (low-texture / dynamic / blur failures)"
          % rep["fallback_pct"])
    print("-" * 100)
    print("ROTATION vs GT  (over recovered frames; B_cam_gt = X_cam^-1 o A_gt o X_cam)")
    print("  per-frame rot err  : median %.3f deg   p90 %.3f deg" % (rep["rot_med"], rep["rot_p90"]))
    print("  per-step YAW CORR  : %.3f   (recovered-vs-GT; HEADLINE: usable rotation?)"
          % rep["yaw_corr"])
    print("  cumulative heading : recovered %+.2f deg   vs GT %+.2f deg   (err %+.2f deg)"
          % (rep["cum_est_yaw_deg"], rep["cum_gt_yaw_deg"],
             rep["cum_est_yaw_deg"] - rep["cum_gt_yaw_deg"]))
    print("-" * 100)
    print("TRANSLATION DIRECTION vs GT  (magnitude scale-ambiguous -> NOT scored)")
    print("  t-direction err    : median %.3f deg   p90 %.3f deg" % (rep["dir_med"], rep["dir_p90"]))
    print("=" * 100)


class Args:
    pass


def parse_args(argv):
    a = Args()
    a.max_features = 1500
    a.min_inliers = 30
    a.ransac_prob = 0.999
    a.ransac_thresh = 1.0     # px (Essential-matrix Sampson distance)
    a.fb_thresh = 1.0         # px forward-backward KLT consistency
    a.raw_cam_frame = False   # emit in raw camera optical frame (leave the mount in) for calibration
    a.seed = 12345
    a.quiet = False
    pos = []
    i = 0
    while i < len(argv):
        x = argv[i]
        if x == "--max-features":
            i += 1; a.max_features = int(argv[i])
        elif x == "--min-inliers":
            i += 1; a.min_inliers = int(argv[i])
        elif x == "--ransac-prob":
            i += 1; a.ransac_prob = float(argv[i])
        elif x == "--ransac-thresh":
            i += 1; a.ransac_thresh = float(argv[i])
        elif x == "--fb-thresh":
            i += 1; a.fb_thresh = float(argv[i])
        elif x == "--raw-cam-frame":
            a.raw_cam_frame = True
        elif x == "--seed":
            i += 1; a.seed = int(argv[i])
        elif x == "--quiet":
            a.quiet = True
        else:
            pos.append(x)
        i += 1
    a.pos = pos
    return a


def main():
    args = parse_args(sys.argv[1:])
    if len(args.pos) < 3:
        print("usage: nuscenes_cam_vo.py <dataroot> <out_dir> <scene> [scene ...] "
              "[--max-features N] [--min-inliers N] [--ransac-prob P] [--ransac-thresh PX] "
              "[--fb-thresh PX] [--raw-cam-frame] [--seed N] [--quiet]")
        sys.exit(2)
    cv2.setRNGSeed(args.seed)
    np.random.seed(args.seed)
    dataroot, out_dir = args.pos[0], args.pos[1]
    scenes = args.pos[2:]
    print("loading nuScenes tables...")
    tables = n2c.load_tables(dataroot)
    for scene_name in scenes:
        try:
            rep, sc = convert_scene(dataroot, out_dir, scene_name, tables, args)
        except Exception as e:  # noqa: BLE001
            import traceback
            print("ERR %s: %s" % (scene_name, e))
            traceback.print_exc()
            continue
        print_report(sc, rep, args)


if __name__ == "__main__":
    main()
