#!/usr/bin/env python3
# nuscenes_surround.py - the 6-CAMERA SURROUND real multi-extrinsic recovery on nuScenes. Extends
# the proven single-camera VO (tools/nuscenes_cam_vo.py, the GO de-risk) to ALL 6 nuScenes cameras
# and recovers their mount ROTATION extrinsics SIMULTANEOUSLY through the calibrator -- the real-data
# version of the synthesized camera surround in SURROUND_MULTISOURCE_TEST.md (which used planted
# cameras; here the 6 mounts are the REAL calibrated_sensor extrinsics and the odometry is REAL VO).
#
# WHY this is the proper multi-extrinsic-recovery venue (the radar could not be): a monocular camera
# observes FULL-6DOF motion -- it recovers rotation (per-frame ~0.05-0.10 deg, yaw corr 0.87-0.97 on
# turns, de-risk a39521d) -- so each camera's rotation extrinsic is OBSERVABLE and self-calibrates.
# The heading-blind radar (R=I) never opens the rot3d gate, so its rotation extrinsic stays
# unobservable (the SURROUND_MULTISOURCE_TEST radar finding). 6 cameras at genuinely different mounts
# (front/sides/back, yaws spanning -146..+159 deg) = a real surround multi-extrinsic rig.
#
# THE 8 SOURCES (mirrors the synth-surround + the nuScenes 7-source manifests):
#   src0 (reference) = wheel   : CAN vehicle_monitor (speed + yaw_rate)        [nuscenes_to_csv]
#   src1            = imu       : CAN ms_imu gyro so3 + wheel forward            [nuscenes_to_csv]
#   src2..7         = 6 cameras : monocular VO (ORB+KLT+Essential), RAW camera   [nuscenes_cam_vo]
#                                 optical frame (mount LEFT IN -> a real extrinsic to recover)
#   GT             = ego_pose   : LiDAR-localization global track (eval only)
#
# Each camera's prior_extrinsic = its REAL ego_from_cam (calibrated_sensor) as yaw pitch roll x y z
# (radians), via the same quat_to_ypr the radar manifest uses. split_median=true (cameras differ in
# per-frame quality), NO translation_only (cameras HAVE rotation), NO GPS (ego_pose is GT only). The
# calib recipe: subbin_centroid + vote_weight=one (Slice 16/17), rot3d_enabled + joint_lever_scale
# (the rotation extrinsic + lever, Slice 17/17b). max_sources=8.
#
# TWO manifests per scene (mirrors the synth surround Run C):
#   CTRL    : the 6 TRUE camera mounts as priors -> do the cameras self-calibrate, holding the truth?
#   RECOVER : each camera's prior perturbed +5 deg YAW + +0.10 m on ONE lever axis (varied per
#             camera) -> does the multi-extrinsic calibration recover the 6 rotation extrinsics at
#             once? The HEADLINE: 6 real cameras calibrated simultaneously on real VO.
#
# REUSE, don't reimplement: imports nuscenes_cam_vo (the proven per-camera VO front-end + the mount
# conjugation + the vs-GT validation) and nuscenes_to_csv (wheel/imu/gt builders, table joins,
# quat<->ypr, the increment-CSV emit + manifest helpers). Pure ASCII; SEEDED (cv2 + numpy + RANSAC).
#
# Usage:
#   python tools/nuscenes_surround.py <dataroot> <out_dir> <scene> [scene ...]
#       [--prefix NAME] [--min-inliers N] [--max-features N] [--seed N]
#       [--perturb-yaw-deg D] [--perturb-lever-m M] [--quiet]
#   e.g. python tools/nuscenes_surround.py C:/workspace/data/nuScenes nuscenes_run scene-0061
import math
import os
import sys

import numpy as np

import nuscenes_to_csv as n2c
import nuscenes_cam_vo as cvo

# the 6 surround cameras, source ids 2..7 (wheel=0, imu=1 from nuscenes_to_csv). Order: front, the
# two front corners, back, the two back corners -- spanning the full yaw circle.
CAM_CHANNELS = ["CAM_FRONT", "CAM_FRONT_LEFT", "CAM_FRONT_RIGHT",
                "CAM_BACK", "CAM_BACK_LEFT", "CAM_BACK_RIGHT"]
CAM_SRC_BASE = 2
# short csv suffix per channel (lower-cased channel)
def short_of(ch):
    return ch.lower()


# per-camera RECOVER perturbation: +perturb_yaw on the extrinsic yaw, +perturb_lever on ONE lever
# axis, the axis VARIED per camera (mirrors the synth surround's "varied per camera" perturbation
# and probes both the along-motion (often unobservable on planar drives) and lateral lever axes).
PERTURB_LEVER_AXIS = {
    "CAM_FRONT": 0,        # x (along-motion) -- often withheld on planar urban
    "CAM_FRONT_LEFT": 1,   # y (lateral)      -- observable
    "CAM_FRONT_RIGHT": 1,  # y (lateral)
    "CAM_BACK": 0,         # x (along-motion)
    "CAM_BACK_LEFT": 1,    # y (lateral)
    "CAM_BACK_RIGHT": 1,   # y (lateral)
}


def yaw_of(R):
    return math.atan2(R[1][0], R[0][0])


def best_available_heading(per_cam_perframe, gt_ts):
    # 6-cam-vs-1-cam REDUNDANCY: at each GT-anchored time bin take the recovered yaw increment from
    # the camera whose per-frame rotation is the LEAST in error vs GT (cameras degenerate on
    # DIFFERENT frames -> a cross-camera pick should close the cumulative-heading gap a single camera
    # leaves). Returns (cum_best_deg, cum_gt_deg). We bin per camera by frame index against GT yaw
    # increments; the "best" picks, per frame, the camera with the smallest |est-gt| yaw increment.
    # This is the ORACLE redundancy ceiling (uses GT to pick) -- it answers "is the heading RECOVERABLE
    # from the surround if you always trusted the locally-best camera?", the redundancy question.
    # We also compute the front-only cumulative for the contrast.
    # Align by collecting, per camera, the list of (recovered, est_yaw_inc, gt_yaw_inc) in frame order
    # and walking them in lockstep up to the shortest length (cameras have ~equal frame counts).
    seqs = []
    for pf in per_cam_perframe:
        s = []
        for rec in pf:
            (R_inc, t_dir, Rg, tg, recovered, ts_us, t_prev_us) = rec
            s.append((recovered, yaw_of(R_inc), yaw_of(Rg)))
        seqs.append(s)
    if not seqs:
        return 0.0, 0.0
    nmin = min(len(s) for s in seqs)
    cum_best = 0.0
    cum_gt = 0.0
    for k in range(nmin):
        gt_inc = seqs[0][k][2]  # GT yaw increment (same GT for all; front cam's bin is representative)
        cum_gt += gt_inc
        # pick the camera with the smallest |est - gt| this frame among recovered cameras
        best_err = None
        best_est = gt_inc
        for s in seqs:
            recovered, est, gtinc = s[k]
            if not recovered:
                continue
            e = abs(est - gtinc)
            if best_err is None or e < best_err:
                best_err = e
                best_est = est
        cum_best += best_est
    return math.degrees(cum_best), math.degrees(cum_gt)


# ---------- manifest writer (8 sources: wheel, imu, 6 cameras) ----------
def write_surround_manifest(out_dir, scene_name, prefix, cam_priors, recover):
    # cam_priors: list of (sid, short, channel, yaw, pitch, roll, lx, ly, lz) for the 6 cameras.
    # recover=False -> CTRL (true mounts); recover=(yaw_d, lever_m) -> RECOVER (perturbed priors).
    tag = "recover" if recover else "ctrl"
    name = "%s_%s" % (prefix, tag)
    L = []
    L.append("[global]")
    L.append("tick_rate_hz = 50")
    L.append("fusion_delay_s = 0.05")
    L.append("window_s = 0.1")
    L.append("max_sources = 8")
    L.append("reference_sensor_id = 0")
    L.append("cold_start = median_from_start")
    L.append("split_median = true")
    L.append("timesync_enabled = false")
    L.append("vote_weight = one")
    L.append("subbin_centroid = true")
    L.append("rot3d_enabled = true")
    L.append("joint_lever_scale = true")
    L.append("heading_monitor = true")
    L.append("adaptive_q = true")
    L.append("q_scale = 0.7")
    L.append("q_floor = 0.01 0.01 0.01 0.001 0.001 0.001")
    L.append("")
    L.append("[sensor.0]")
    L.append("id = 0")
    L.append("is_reference = true")
    L.append("csv = %s_wheel.csv" % prefix)
    L.append("form = increment")
    L.append("")
    L.append("[sensor.1]")
    L.append("id = 1")
    L.append("csv = %s_imu.csv" % prefix)
    L.append("form = increment")
    L.append("")
    truth_lines = []
    for (sid, short, channel, yaw, pitch, roll, lx, ly, lz) in cam_priors:
        p_yaw, p_pitch, p_roll, p_lx, p_ly, p_lz = yaw, pitch, roll, lx, ly, lz
        if recover:
            yaw_d, lever_m = recover
            p_yaw = yaw + math.radians(yaw_d)
            axis = PERTURB_LEVER_AXIS.get(channel, 1)
            lev = [lx, ly, lz]
            lev[axis] += lever_m
            p_lx, p_ly, p_lz = lev
        L.append("[sensor.%d]" % sid)
        L.append("id = %d" % sid)
        L.append("csv = %s_%s.csv" % (prefix, short))
        L.append("form = increment")
        L.append("prior_extrinsic = %.9f %.9f %.9f %.6f %.6f %.6f"
                 % (p_yaw, p_pitch, p_roll, p_lx, p_ly, p_lz))
        L.append("")
        truth_lines.append("src%d %-16s true_yaw=%+.6f rad (%+.3f deg)  prior_yaw=%+.6f rad (%+.3f deg)  "
                           "true_lever=[%.4f %.4f %.4f] perturb_axis=%s"
                           % (sid, channel, yaw, math.degrees(yaw), p_yaw, math.degrees(p_yaw),
                              lx, ly, lz, "xyz"[PERTURB_LEVER_AXIS.get(channel, 1)] if recover else "-"))
    L.append("[gt]")
    L.append("csv = %s_gt.csv" % prefix)
    L.append("")
    L.append("[replay]")
    L.append("warmup_steps = 20")
    L.append("local_batch_len = 500")
    L.append("out = %s_out.csv" % name)
    L.append("")
    ini_path = os.path.join(out_dir, name + ".ini")
    with open(ini_path, "w", newline="") as f:
        f.write("\n".join(L))
    # truth sidecar (for the recovery scoring readout)
    with open(os.path.join(out_dir, name + "_truth.txt"), "w", newline="") as f:
        f.write("\n".join(truth_lines) + "\n")
    return ini_path


# ---------- wheel/imu/gt emit (reuse nuscenes_to_csv builders) ----------
def emit_can_and_gt(dataroot, out_dir, scene_name, prefix, tables, gt, t0_us, gt_R):
    # Reuse nuscenes_to_csv's wheel/imu/gt builders + the same unit self-check, so src0/src1 and the
    # GT match the proven 7-source run exactly. Returns (wheel_n, imu_n, wheel_ratio).
    vm = n2c.load_can(dataroot, scene_name, "vehicle_monitor")
    imu = n2c.load_can(dataroot, scene_name, "ms_imu")
    if not vm:
        raise ValueError("no vehicle_monitor CAN for %s" % scene_name)
    if not imu:
        raise ValueError("no ms_imu CAN for %s" % scene_name)
    # GT track first-relative, plus the unit self-check (mirrors nuscenes_to_csv.convert)
    R0 = n2c.quat_wxyz_to_mat(*gt[0]["rotation"])
    R0t = n2c.mat_t(R0)
    p0 = gt[0]["translation"]
    gt_rows = []
    gt_xy = []
    for p in gt:
        Ri = n2c.quat_wxyz_to_mat(*p["rotation"])
        ti = p["translation"]
        d = [ti[0] - p0[0], ti[1] - p0[1], ti[2] - p0[2]]
        trel = n2c.mat_vec(R0t, d)
        Rrel = n2c.mat_mul(R0t, Ri)
        qw, qx, qy, qz = n2c.mat_to_quat(Rrel)
        gt_rows.append((p["timestamp"] - t0_us, trel, (qw, qx, qy, qz)))
        gt_xy.append((trel[0], trel[1]))
    gt_len = n2c.path_length(gt_xy)

    speed_choices = [("km/h->m/s", 1.0 / 3.6), ("m/s", 1.0)]
    yaw_choices = [("deg/s->rad/s", math.pi / 180.0), ("rad/s", 1.0)]
    best = None
    gt_Rrel = [n2c.mat_mul(R0t, n2c.quat_wxyz_to_mat(*p["rotation"])) for p in gt]
    for (sname, sval) in speed_choices:
        for (yname, yval) in yaw_choices:
            _rows, incs = n2c.build_wheel(vm, t0_us, sval, yval)
            xy = n2c.integrate_se2(incs)
            L = n2c.path_length(xy)
            dist_ratio = L / gt_len if gt_len > 1e-6 else 0.0
            est_yaw = sum(math.atan2(Rd[1][0], Rd[0][0]) for (Rd, _t) in incs)
            gt_net_yaw = math.atan2(gt_Rrel[-1][1][0], gt_Rrel[-1][0][0])
            yaw_err = abs(((est_yaw - gt_net_yaw + math.pi) % (2 * math.pi)) - math.pi)
            score = abs(dist_ratio - 1.0) + 0.2 * yaw_err
            cand = (score, sname, sval, yname, yval, dist_ratio)
            if best is None or cand[0] < best[0]:
                best = cand
    _, sname, speed_scale, yname, yaw_scale, wheel_ratio = best

    base = os.path.join(out_dir, prefix)
    # wheel
    wheel_rows, _ = n2c.build_wheel(vm, t0_us, speed_scale, yaw_scale)
    fw = n2c.open_inc(base + "_wheel.csv", "wheel src0: vehicle_monitor speed+yaw_rate (%s,%s)" % (sname, yname))
    wheel_var = [4e-4, 1e-2, 1e-2, 1e-4, 1e-2, 1e-2]
    for (t, R, tt) in wheel_rows:
        n2c.emit_inc(fw, t * n2c.US_TO_NS, R, tt, wheel_var)
    fw.close()
    # imu
    imu_rows, _ = n2c.build_imu(imu, vm, t0_us, speed_scale)
    fi = n2c.open_inc(base + "_imu.csv", "imu src1: ms_imu gyro so3 + wheel forward")
    imu_var = [4e-4, 1e-2, 1e-2, 1e-5, 1e-5, 1e-5]
    for (t, R, tt) in imu_rows:
        n2c.emit_inc(fi, t * n2c.US_TO_NS, R, tt, imu_var)
    fi.close()
    # gt
    with open(base + "_gt.csv", "w", newline="") as fgt:
        fgt.write("# t_ns, x, y, z, qw, qx, qy, qz  (ego_pose GT, first-relative)\n")
        for (t, trel, q) in gt_rows:
            fgt.write("%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f,%.9f\n"
                      % (t * n2c.US_TO_NS, trel[0], trel[1], trel[2], q[0], q[1], q[2], q[3]))
    return len(wheel_rows), len(imu_rows), wheel_ratio, gt_len


# ---------- per-scene surround build ----------
def build_scene(dataroot, out_dir, scene_name, tables, args):
    prefix = args.prefix if args.prefix else scene_name
    (sample_tokens, gt, t0_us, gt_ts, gt_R, gt_p,
     vm, speed_scale) = cvo.scene_gt_and_speed(dataroot, scene_name, tables)
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, prefix)

    # FORCE raw-cam-frame for the surround: the mount must be LEFT IN so the calibrator has a real
    # per-camera extrinsic to recover. (The body-frame mode is only for the vs-GT VO validation,
    # which we still run below against the body-frame perframe regardless of emit mode.)
    args.raw_cam_frame = True

    wheel_n, imu_n, wheel_ratio, gt_len = emit_can_and_gt(
        dataroot, out_dir, scene_name, prefix, tables, gt, t0_us, gt_R)

    cam_priors = []
    cam_reports = []
    per_cam_perframe = []
    for ri, channel in enumerate(CAM_CHANNELS):
        sid = CAM_SRC_BASE + ri
        calib, stream, vo = cvo.run_channel_vo(channel, dataroot, scene_name, tables, sample_tokens,
                                               t0_us, gt_ts, gt_R, gt_p, vm, speed_scale, args)
        if stream is None:
            print("  [warn] %s: no frames in %s; skipping" % (channel, scene_name))
            continue
        rows, stats, perframe = vo
        # emit the raw-cam-frame increment CSV for this camera
        short = short_of(channel)
        outp = base + "_" + short + ".csv"
        cvo.emit_vo_csv(outp, "%s monocular VO: ORB+KLT+Essential (R + t-DIRECTION); raw camera "
                        "optical frame (mount left in for calibration); scale=calibrator"
                        % channel, rows, t0_us)
        # prior_extrinsic = the REAL calibrated_sensor mount as yaw pitch roll x y z (radians)
        qw, qx, qy, qz = calib["rotation"]
        yaw, pitch, roll = n2c.quat_to_ypr(qw, qx, qy, qz)
        lx, ly, lz = calib["translation"]
        cam_priors.append((sid, short, channel, yaw, pitch, roll, lx, ly, lz))
        rep = cvo.report_cam(channel, perframe, stats)
        rep["sid"] = sid
        rep["n_frames"] = len(stream)
        cam_reports.append(rep)
        per_cam_perframe.append(perframe)

    # 6-cam-vs-1-cam redundancy heading (oracle best-available vs front-only vs GT)
    cum_best, cum_gt = best_available_heading(per_cam_perframe, gt_ts)
    front_rep = next((r for r in cam_reports if r["channel"] == "CAM_FRONT"), None)
    cum_front = front_rep["cum_est_yaw_deg"] if front_rep else float("nan")

    # write the two manifests
    write_surround_manifest(out_dir, scene_name, prefix, cam_priors, recover=False)
    write_surround_manifest(out_dir, scene_name, prefix, cam_priors,
                            recover=(args.perturb_yaw_deg, args.perturb_lever_m))

    return {
        "scene": scene_name, "prefix": prefix, "wheel_n": wheel_n, "imu_n": imu_n,
        "wheel_ratio": wheel_ratio, "gt_len": gt_len, "dur_s": (gt[-1]["timestamp"] - t0_us) * 1e-6,
        "cam_reports": cam_reports, "cum_best_deg": cum_best, "cum_gt_deg": cum_gt,
        "cum_front_deg": cum_front, "n_cams": len(cam_priors),
    }


def print_scene(info, args):
    print("=" * 108)
    print("6-CAMERA SURROUND VO -- %s  (prefix %s, %d cams, %.1f s, gt_len %.1f m)  "
          "[min_inl=%d max_feat=%d]"
          % (info["scene"], info["prefix"], info["n_cams"], info["dur_s"], info["gt_len"],
             args.min_inliers, args.max_features))
    print("=" * 108)
    print("  wheel rows=%d (dist_ratio %.3f)   imu rows=%d" %
          (info["wheel_n"], info["wheel_ratio"], info["imu_n"]))
    print("-" * 108)
    print("PER-CAMERA VO vs ego_pose GT  (B_cam_gt = X_cam^-1 o A_gt o X_cam, body frame)")
    print("  %-16s %4s  %7s %7s  %6s  %7s %7s  %7s  %8s  %8s" %
          ("camera", "sid", "rotMed", "rotP90", "fb%", "avgInl", "recov", "yawCorr",
           "cumEst", "cumGT"))
    for r in info["cam_reports"]:
        print("  %-16s %4d  %6.3fd %6.3fd  %5.1f  %7.1f %4d/%-3d %7.3f  %+7.1f  %+7.1f" %
              (r["channel"], r["sid"], r["rot_med"], r["rot_p90"], r["fallback_pct"],
               r["avg_inl"], r["n_recov"], r["n_pairs"], r["yaw_corr"],
               r["cum_est_yaw_deg"], r["cum_gt_yaw_deg"]))
    print("-" * 108)
    print("6-CAM vs 1-CAM cumulative heading vs GT  (redundancy: cross-camera best-available)")
    print("  GT net yaw            : %+8.2f deg" % info["cum_gt_deg"])
    print("  CAM_FRONT only        : %+8.2f deg   (err %+7.2f deg)" %
          (info["cum_front_deg"], info["cum_front_deg"] - info["cum_gt_deg"]))
    print("  6-cam best-available  : %+8.2f deg   (err %+7.2f deg)   <- redundancy ceiling" %
          (info["cum_best_deg"], info["cum_best_deg"] - info["cum_gt_deg"]))
    print("=" * 108)


class Args:
    pass


def parse_args(argv):
    a = Args()
    a.prefix = ""
    a.min_inliers = 15
    a.max_features = 1500
    a.ransac_prob = 0.999
    a.ransac_thresh = 1.0
    a.fb_thresh = 1.0
    a.raw_cam_frame = True
    a.perturb_yaw_deg = 5.0
    a.perturb_lever_m = 0.10
    a.seed = 12345
    a.quiet = False
    pos = []
    i = 0
    while i < len(argv):
        x = argv[i]
        if x == "--prefix":
            i += 1; a.prefix = argv[i]
        elif x == "--min-inliers":
            i += 1; a.min_inliers = int(argv[i])
        elif x == "--max-features":
            i += 1; a.max_features = int(argv[i])
        elif x == "--ransac-prob":
            i += 1; a.ransac_prob = float(argv[i])
        elif x == "--ransac-thresh":
            i += 1; a.ransac_thresh = float(argv[i])
        elif x == "--fb-thresh":
            i += 1; a.fb_thresh = float(argv[i])
        elif x == "--perturb-yaw-deg":
            i += 1; a.perturb_yaw_deg = float(argv[i])
        elif x == "--perturb-lever-m":
            i += 1; a.perturb_lever_m = float(argv[i])
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
        print("usage: nuscenes_surround.py <dataroot> <out_dir> <scene> [scene ...] "
              "[--prefix NAME] [--min-inliers N] [--max-features N] "
              "[--perturb-yaw-deg D] [--perturb-lever-m M] [--seed N] [--quiet]")
        sys.exit(2)
    cvo.cv2.setRNGSeed(args.seed)
    np.random.seed(args.seed)
    dataroot, out_dir = args.pos[0], args.pos[1]
    scenes = args.pos[2:]
    print("loading nuScenes tables...")
    tables = n2c.load_tables(dataroot)
    for scene_name in scenes:
        try:
            info = build_scene(dataroot, out_dir, scene_name, tables, args)
        except Exception as e:  # noqa: BLE001
            import traceback
            print("ERR %s: %s" % (scene_name, e))
            traceback.print_exc()
            continue
        print_scene(info, args)


if __name__ == "__main__":
    main()
