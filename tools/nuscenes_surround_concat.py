#!/usr/bin/env python3
# nuscenes_surround_concat.py - WARM the 6-camera LEVER estimators past the 20 s scene limit by
# ACCUMULATING surround VO across ALL 10 nuScenes-mini scenes into ONE persistent calibrator. This
# composes the two proven tools:
#   * nuscenes_surround.py  - the 6-cam VO front-end + the REAL ego_from_cam mounts as priors + the
#                             8-source surround manifest (CAN wheel ref + imu + 6 VO cameras).
#   * nuscenes_concat.py    - the increment-space concatenation across the 10 mini scenes (re-base
#                             each scene's timestamps to 0, offset by cumulative end-time + one tick,
#                             ONE identity seed per source at the very start, per-delta covariance
#                             preserved; teleport lives only in absolute GT, which calib never uses).
#
# WHY (the headline question): the 6-cam surround (nuscenes_surround.py) recovered 5/6 camera ROTATION
# extrinsics on real VO, but the camera LEVERS stayed conf-0 -- diagnosed as SCENE LENGTH (20 s too
# short to warm the lever estimator), NOT observability (cameras HAVE rotation, so the hand-eye lever
# rows are well-posed -- unlike the heading-blind radar whose lever is structurally stuck). The radar
# concat (nuscenes_concat.py) proved increment-space concatenation accumulates votes into ONE
# persistent calibrator over ~200 s. REUSE that pattern for the camera VO streams: the camera MOUNTS
# are CONSTANT across scenes, so the lever/rotation votes accumulate into ONE persistent calibrator.
#
# THE 8 CONCATENATED SOURCES (mirrors the surround manifest, source-id order):
#   src0 (reference) = wheel   : CAN vehicle_monitor                                  [nuscenes_to_csv]
#   src1            = imu       : CAN ms_imu gyro + wheel forward                       [nuscenes_to_csv]
#   src2..7         = 6 cameras : monocular VO (ORB+KLT+Essential), RAW camera optical  [nuscenes_cam_vo]
#                                 frame (real ego_from_cam mount left in -> a real extrinsic to recover)
#   GT             = ego_pose   : LiDAR-localization track (eval only; teleports at scene boundaries)
#
# COMPUTE: step 1 (the heavy step, 6 cams x ~250 frames x 10 scenes) is CACHED per scene -- a scene is
# skipped if its per-scene gt CSV already exists, so reruns are cheap. step 2 (concat) is instant.
#
# MANIFEST: mirrors the surround RECOVER manifest -- each camera prior_extrinsic = its REAL mount
# PERTURBED by +5 deg yaw + +0.10 m on ONE lever axis (varied per camera, SAME as the surround RECOVER
# so it is directly comparable), split_median=true, NO translation_only, NO GPS, FULL calib recipe
# (subbin_centroid + vote_weight=one + rot3d_enabled + joint_lever_scale), GT track, max_sources=8.
#
# Pure ASCII; numpy + opencv (via the imported VO); SEEDED (cv2.setRNGSeed + numpy seed).
#
# Usage:
#   python tools/nuscenes_surround_concat.py <dataroot> <out_dir> [scene ...]
#       [--name BASENAME] [--perturb-yaw-deg D] [--perturb-lever-m M]
#       [--min-inliers N] [--max-features N] [--seed N] [--keep-temp]
#   default scenes = all 10 v1.0-mini scenes. Emits <out_dir>/<BASENAME>_{wheel,imu,cam_*}.csv,
#   <BASENAME>_gt.csv, <BASENAME>_recover.ini, <BASENAME>_recover_truth.txt.
import os
import shutil
import sys

import numpy as np

import nuscenes_to_csv as n2c
import nuscenes_surround as sur
import nuscenes_concat as cc

# all 10 v1.0-mini scenes (scene.json order) -- same list nuscenes_concat uses.
ALL_MINI = list(cc.ALL_MINI)

# per-source file suffixes for the 8-source surround (must match nuscenes_surround naming:
# <prefix>_wheel.csv, <prefix>_imu.csv, <prefix>_<channel.lower()>.csv). Order = source-id order.
SRC_SUFFIXES = ["wheel", "imu"] + [sur.short_of(ch) for ch in sur.CAM_CHANNELS]


def build_per_scene(dataroot, scene_dir, scene, tables, args):
    # Generate (or reuse cached) the per-scene surround CSVs for one scene. The heavy 6-cam VO runs
    # only if the per-scene gt CSV is absent (skip-if-exists). Returns the build_scene info dict
    # (carries the per-camera VO reports for the single-scene-vs-accumulated readout).
    prefix = scene
    gt_path = os.path.join(scene_dir, prefix + "_gt.csv")
    cam_done = all(os.path.exists(os.path.join(scene_dir, "%s_%s.csv" % (prefix, suf)))
                   for suf in SRC_SUFFIXES)
    if os.path.exists(gt_path) and cam_done:
        print("=== %s: cached (skip VO) ===" % scene)
        # still need the per-scene info for the report + the true mounts; re-derive cheaply from the
        # cached per-scene CTRL manifest priors below (mounts constant) -- return None to signal cache.
        return None
    print("=== %s: running 6-cam VO (heavy) ===" % scene)
    args.prefix = prefix
    info = sur.build_scene(dataroot, scene_dir, scene, tables, args)
    return info


def read_cam_priors_from_manifest(ini_path):
    # Read the 6 camera priors (TRUE mounts, source-id order) from a per-scene CTRL surround manifest.
    # Returns list of (sid, short, yaw, pitch, roll, lx, ly, lz). Mounts are CONSTANT across scenes so
    # any scene's CTRL manifest carries the truth. (The CTRL manifest holds the UNPERTURBED priors.)
    priors = []
    with open(ini_path) as f:
        lines = f.read().splitlines()
    sid = csvname = prior = None
    for line in lines + [""]:
        s = line.strip()
        if s.startswith("[sensor.") or s == "":
            if sid is not None and prior is not None and csvname and sid >= sur.CAM_SRC_BASE:
                short = os.path.basename(csvname)[:-4].split("_", 1)[1]  # "<scene>_cam_front" -> "cam_front"
                priors.append((sid, short, prior[0], prior[1], prior[2],
                               prior[3], prior[4], prior[5]))
            sid = csvname = prior = None
        if s.startswith("id ="):
            sid = int(s.split("=", 1)[1])
        elif s.startswith("csv ="):
            csvname = s.split("=", 1)[1].strip()
        elif s.startswith("prior_extrinsic ="):
            prior = [float(v) for v in s.split("=", 1)[1].split()]
    priors.sort(key=lambda p: p[0])
    return priors


def channel_of(short):
    # "cam_front" -> "CAM_FRONT"; matches sur.short_of inverse.
    return short.upper()


def write_recover_manifest(out_dir, basename, cam_priors, perturb_yaw_deg, perturb_lever_m):
    # Mirror sur.write_surround_manifest's RECOVER manifest, but pointing at the CONCATENATED CSVs and
    # using <basename>_<suffix>.csv naming. Each camera prior = REAL mount perturbed +yaw + +lever on
    # one varied axis (sur.PERTURB_LEVER_AXIS), SAME perturbation as the single-scene surround RECOVER.
    import math
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
    L.append("csv = %s_wheel.csv" % basename)
    L.append("form = increment")
    L.append("")
    L.append("[sensor.1]")
    L.append("id = 1")
    L.append("csv = %s_imu.csv" % basename)
    L.append("form = increment")
    L.append("")
    truth_lines = []
    for (sid, short, yaw, pitch, roll, lx, ly, lz) in cam_priors:
        channel = channel_of(short)
        p_yaw = yaw + math.radians(perturb_yaw_deg)
        axis = sur.PERTURB_LEVER_AXIS.get(channel, 1)
        lev = [lx, ly, lz]
        lev[axis] += perturb_lever_m
        p_lx, p_ly, p_lz = lev
        L.append("[sensor.%d]" % sid)
        L.append("id = %d" % sid)
        L.append("csv = %s_%s.csv" % (basename, short))
        L.append("form = increment")
        L.append("prior_extrinsic = %.9f %.9f %.9f %.6f %.6f %.6f"
                 % (p_yaw, pitch, roll, p_lx, p_ly, p_lz))
        L.append("")
        truth_lines.append("src%d %-16s true_yaw=%+.6f rad (%+.3f deg)  prior_yaw=%+.6f rad (%+.3f deg)  "
                           "true_lever=[%.4f %.4f %.4f] perturb_axis=%s"
                           % (sid, channel, yaw, math.degrees(yaw), p_yaw, math.degrees(p_yaw),
                              lx, ly, lz, "xyz"[axis]))
    L.append("[gt]")
    L.append("csv = %s_gt.csv" % basename)
    L.append("")
    L.append("[replay]")
    L.append("warmup_steps = 20")
    L.append("local_batch_len = 500")
    L.append("out = %s_recover_out.csv" % basename)
    L.append("")
    ini_path = os.path.join(out_dir, basename + "_recover.ini")
    with open(ini_path, "w", newline="") as f:
        f.write("\n".join(L))
    with open(os.path.join(out_dir, basename + "_recover_truth.txt"), "w", newline="") as f:
        f.write("\n".join(truth_lines) + "\n")
    return ini_path


def main():
    argv = sys.argv[1:]
    seed = 12345
    basename = "mini_surround_concat"
    keep_temp = False
    perturb_yaw_deg = 5.0
    perturb_lever_m = 0.10
    min_inliers = 15
    max_features = 1500
    pos = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--seed":
            i += 1; seed = int(argv[i])
        elif a == "--name":
            i += 1; basename = argv[i]
        elif a == "--keep-temp":
            keep_temp = True
        elif a == "--perturb-yaw-deg":
            i += 1; perturb_yaw_deg = float(argv[i])
        elif a == "--perturb-lever-m":
            i += 1; perturb_lever_m = float(argv[i])
        elif a == "--min-inliers":
            i += 1; min_inliers = int(argv[i])
        elif a == "--max-features":
            i += 1; max_features = int(argv[i])
        else:
            pos.append(a)
        i += 1
    if len(pos) < 2:
        print("usage: nuscenes_surround_concat.py <dataroot> <out_dir> [scene ...] "
              "[--name BASENAME] [--perturb-yaw-deg D] [--perturb-lever-m M] "
              "[--min-inliers N] [--max-features N] [--seed N] [--keep-temp]")
        sys.exit(2)
    dataroot, out_dir = pos[0], pos[1]
    scenes = pos[2:] if len(pos) > 2 else list(ALL_MINI)

    sur.cvo.cv2.setRNGSeed(seed)
    np.random.seed(seed)

    os.makedirs(out_dir, exist_ok=True)
    scene_dir = os.path.join(out_dir, "_per_scene_surround")
    os.makedirs(scene_dir, exist_ok=True)

    # surround VO args (same defaults as nuscenes_surround.parse_args, raw cam frame forced).
    args = sur.Args()
    args.prefix = ""
    args.min_inliers = min_inliers
    args.max_features = max_features
    args.ransac_prob = 0.999
    args.ransac_thresh = 1.0
    args.fb_thresh = 1.0
    args.raw_cam_frame = True
    args.perturb_yaw_deg = perturb_yaw_deg
    args.perturb_lever_m = perturb_lever_m
    args.seed = seed
    args.quiet = False

    print("loading nuScenes tables...")
    tables = n2c.load_tables(dataroot)

    # --- step 1: per-scene 6-cam VO (cached) ---
    used_scenes = []
    scene_infos = []
    for scene in scenes:
        try:
            info = build_per_scene(dataroot, scene_dir, scene, tables, args)
        except Exception as e:  # noqa: BLE001
            import traceback
            print("ERR %s: %s" % (scene, e))
            traceback.print_exc()
            continue
        # confirm the per-scene gt landed (a scene with no frames is skipped)
        if not os.path.exists(os.path.join(scene_dir, scene + "_gt.csv")):
            print("  [warn] %s produced no per-scene CSVs; dropping from concat" % scene)
            continue
        used_scenes.append(scene)
        scene_infos.append(info)  # may be None for a cached scene
    if not used_scenes:
        print("no scenes generated; abort")
        sys.exit(1)

    # capture the TRUE camera mounts (constant across scenes) from the first used scene's CTRL manifest
    cam_priors = read_cam_priors_from_manifest(os.path.join(scene_dir, used_scenes[0] + "_ctrl.ini"))
    if len(cam_priors) != len(sur.CAM_CHANNELS):
        print("  [warn] expected %d camera priors, got %d (some camera missing frames in scene %s)"
              % (len(sur.CAM_CHANNELS), len(cam_priors), used_scenes[0]))

    # --- step 2: compute per-scene boundary offsets (one shared span per scene) ---
    TICK_GAP_NS = 20_000_000  # 20 ms (one 50 Hz tick) gap between scenes -> strict monotonicity
    offsets = []
    cum = 0
    spans = []
    for scene in used_scenes:
        offsets.append(cum)
        span = cc.scene_span_ns(scene_dir, scene, SRC_SUFFIXES)
        spans.append(span)
        cum += span + TICK_GAP_NS
    total_span_s = cum / 1e9

    # --- step 3: concatenate each source + GT into the megastream ---
    base = os.path.join(out_dir, basename)
    counts = {}
    for suf in SRC_SUFFIXES:
        outp = "%s_%s.csv" % (base, suf)
        counts[suf] = cc.concat_source(outp, used_scenes, scene_dir, suf,
                                       "CONCAT %s over %d scenes" % (suf, len(used_scenes)),
                                       offsets, is_gt=False)
    counts["gt"] = cc.concat_source(base + "_gt.csv", used_scenes, scene_dir, "gt",
                                    "concat gt", offsets, is_gt=True)

    # --- step 4: RECOVER manifest + truth sidecar (perturbed priors) ---
    write_recover_manifest(out_dir, basename, cam_priors, perturb_yaw_deg, perturb_lever_m)

    # --- report concatenation sanity ---
    print("\n=== SURROUND CONCAT SUMMARY (%s) ===" % basename)
    print("scenes used: %d (%s)" % (len(used_scenes), ", ".join(used_scenes)))
    print("total megastream span: %.1f s" % total_span_s)
    print("per-scene spans (s): %s" % ", ".join("%.1f" % (s / 1e9) for s in spans))
    print("perturb: +%.1f deg yaw + +%.2f m lever (axis varied per camera, SAME as surround RECOVER)"
          % (perturb_yaw_deg, perturb_lever_m))
    print("row counts (megastream, excl headers):")
    for suf in SRC_SUFFIXES + ["gt"]:
        print("  %-18s %d" % (suf, counts[suf]))

    # accumulation excitation: total turn + path over GT, summed across scenes
    tot_turn_deg, tot_path_m = accumulation_excitation(scene_dir, used_scenes, scene_infos)
    print("accumulation excitation: total turn %.1f deg / path %.1f m over %d scenes"
          % (tot_turn_deg, tot_path_m, len(used_scenes)))

    # monotonicity check (mirror nuscenes_concat)
    print("monotonicity check:")
    all_ok = True
    for suf in SRC_SUFFIXES + ["gt"]:
        p = "%s_%s.csv" % (base, suf)
        ok, n, first, last = cc.check_monotonic(p)
        all_ok = all_ok and ok
        print("  %-18s %s  rows=%d t[0]=%d t[-1]=%d" % (suf, "OK " if ok else "NONMONO", n, first, last))
    print("ALL MONOTONIC" if all_ok else "*** NON-MONOTONIC ***")

    if not keep_temp:
        # keep the per-scene dir so reruns are cheap (the heavy VO is cached); only purge on request.
        pass
    print("\nmanifest: %s_recover.ini   truth: %s_recover_truth.txt" % (base, base))
    print("(per-scene VO cached in %s -- reruns skip the heavy step)" % scene_dir)


def accumulation_excitation(scene_dir, used_scenes, scene_infos):
    # Total turn (deg) and path length (m) accumulated across the used scenes, read from each scene's
    # GT track (the same teleporting-but-per-scene-valid GT the concat writes). Per-scene net yaw +
    # path are summed -> the abundant-excitation stat (matches the radar concat's report style).
    import math
    tot_turn = 0.0
    tot_path = 0.0
    for scene in used_scenes:
        gt = os.path.join(scene_dir, scene + "_gt.csv")
        if not os.path.exists(gt):
            continue
        xs = []
        ys = []
        yaws = []
        with open(gt) as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split(",")
                x, y = float(parts[1]), float(parts[2])
                qw, qx, qy, qz = (float(parts[4]), float(parts[5]),
                                  float(parts[6]), float(parts[7]))
                # yaw from quat (Rz(yaw)Ry(pitch)Rx(roll)): atan2(2(qw qz + qx qy), 1-2(qy^2+qz^2))
                yaw = math.atan2(2.0 * (qw * qz + qx * qy), 1.0 - 2.0 * (qy * qy + qz * qz))
                xs.append(x); ys.append(y); yaws.append(yaw)
        for k in range(1, len(xs)):
            tot_path += math.hypot(xs[k] - xs[k - 1], ys[k] - ys[k - 1])
        for k in range(1, len(yaws)):
            d = yaws[k] - yaws[k - 1]
            d = (d + math.pi) % (2 * math.pi) - math.pi
            tot_turn += abs(d)
    return math.degrees(tot_turn), tot_path


if __name__ == "__main__":
    main()
