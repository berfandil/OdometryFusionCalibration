#!/usr/bin/env python3
# nuscenes_concat.py - ACCUMULATE calibration across ALL nuScenes-mini scenes into ONE persistent
# calibrator. A single ~20 s mini scene is far too short for online calibration to converge (the
# nuScenes 7-source run pinned the radar mounts as priors and explicitly did NOT recover them). But
# the radar MOUNTS are CONSTANT across every scene of the same vehicle, so calibration VOTES
# ACCUMULATE. The source streams are INCREMENTS (relative motion), so concatenating scenes is
# CONTINUOUS in increment-space: the teleport at a scene boundary lives only in ABSOLUTE poses, and
# calibration never uses absolute poses -- it votes from per-delta source-vs-consensus residuals.
# One persistent calibrator over all 10 mini scenes = ~200 s of accumulated excitation.
#
# REUSE, don't reimplement: this tool IMPORTS tools/nuscenes_to_csv.py and calls its per-scene
# convert() (the radar Doppler RANSAC, the CAN wheel/imu builders, the GT track, the unit self-check,
# the real calibrated_sensor -> prior_extrinsic mapping) UNCHANGED. It then CONCATENATES the per-scene
# per-source increment CSVs into one megastream per source.
#
# CONCATENATION CONTRACT (per source) ----------------------------------------------------------------
#   * Re-base each scene's timestamps to 0, then OFFSET by the running cumulative end-time + ONE tick
#     gap, so the megastream is STRICTLY MONOTONIC across the boundary.
#   * Keep exactly ONE identity-seed row at the very START (scene 0's seed). DROP every subsequent
#     scene's identity-seed row: each later scene's first real increment just CONTINUES the relative-
#     motion stream (an identity increment mid-stream would inject a spurious zero-motion step).
#   * Preserve the per-delta 6-variance columns exactly as the per-scene converter emitted them.
#   * The scene boundary uses ONE shared scene-span (max last-t over all of that scene's sources) so
#     every source is offset by the same boundary -> sources stay time-aligned scene to scene.
#   GT is concatenated the SAME way (re-base + offset). Its ABSOLUTE poses TELEPORT at each boundary
#   -- that is FINE: GT here is only a sanity drift glance, never an input to calibration.
#
# PERTURB MODE (the experiment) ----------------------------------------------------------------------
#   --perturb-lever shifts each radar's prior_extrinsic LEVER by +0.20 m on ONE axis (varied per
#   radar), leaving the rotation prior at truth. Radars are HEADING-BLIND (R=I) so their rotation
#   extrinsic is UNOBSERVABLE by construction -- the recoverable target is the LEVER only. The
#   manifest then carries the perturbed priors; the true mounts are written to a sidecar
#   <out>_truth.txt so the analysis can compare recovered-vs-true-vs-perturbed.
#
# Pure ASCII, numpy ok (via the imported converter), SEEDED (passes --seed through to convert()).
#
# Usage:
#   python tools/nuscenes_concat.py <dataroot> <out_dir> [scene ...] [--seed N]
#                                   [--perturb-lever] [--name BASENAME] [--keep-temp]
#   default scenes = all 10 v1.0-mini scenes. Emits <out_dir>/<BASENAME>_{wheel,imu,radar_*}.csv,
#   <BASENAME>_gt.csv, <BASENAME>.ini, <BASENAME>_truth.txt.
import os
import shutil
import sys

import nuscenes_to_csv as n2c

# all 10 v1.0-mini scenes (scene.json order)
ALL_MINI = ["scene-0061", "scene-0103", "scene-0553", "scene-0655", "scene-0757",
            "scene-0796", "scene-0916", "scene-1077", "scene-1094", "scene-1100"]

# per-source file suffixes (must match nuscenes_to_csv naming). Order = source-id order.
SRC_SUFFIXES = ["wheel", "imu",
                "radar_front", "radar_front_left", "radar_front_right",
                "radar_back_left", "radar_back_right"]

# Per-radar lever perturbation: +0.20 m on a DIFFERENT axis per radar (axis index 0=x,1=y,2=z).
# Varying the axis exercises all three lever DOFs across the 5 radars.
PERTURB_AXIS = {
    "radar_front": 0,        # +x
    "radar_front_left": 1,   # +y
    "radar_front_right": 2,  # +z
    "radar_back_left": 0,    # +x
    "radar_back_right": 1,   # +y
}
PERTURB_M = 0.20


def read_csv_rows(path):
    # Return (comment_lines, seed_line_or_None, data_lines). data_lines EXCLUDES the identity seed.
    comments, data = [], []
    seed = None
    with open(path) as f:
        for line in f:
            s = line.rstrip("\n")
            if not s.strip():
                continue
            if s.startswith("#"):
                comments.append(s)
                continue
            # first non-comment row is the identity seed (t_ns=0, identity quat)
            if seed is None:
                seed = s
                continue
            data.append(s)
    return comments, seed, data


def parse_t(line):
    return int(line.split(",", 1)[0])


def reoffset(line, old_to_new_offset):
    # shift the t_ns (first field) by offset; keep the rest byte-identical
    i = line.index(",")
    t = int(line[:i]) + old_to_new_offset
    return "%d%s" % (t, line[i:])


def scene_span_ns(scene_dir, scene, suffixes):
    # the shared scene boundary = max last real t_ns across all of this scene's sources + GT.
    span = 0
    for suf in suffixes + ["gt"]:
        p = os.path.join(scene_dir, "%s_%s.csv" % (scene, suf))
        if not os.path.exists(p):
            continue
        _c, _seed, data = read_csv_rows(p)
        if data:
            span = max(span, parse_t(data[-1]))
    return span


def concat_source(out_path, scenes, scene_dir, suffix, header_comment, offsets, is_gt=False):
    # Concatenate one source (or GT) across scenes into out_path. ONE seed at the very start (non-GT).
    total = 0
    with open(out_path, "w", newline="") as g:
        if is_gt:
            g.write("# t_ns, x, y, z, qw, qx, qy, qz  (CONCAT GT; abs poses TELEPORT at boundaries)\n")
        else:
            g.write(n2c.INC_HEADER)
            g.write("# %s\n" % header_comment)
        wrote_seed = False
        for si, scene in enumerate(scenes):
            p = os.path.join(scene_dir, "%s_%s.csv" % (scene, suffix))
            if not os.path.exists(p):
                continue
            _c, seed, data = read_csv_rows(p)
            off = offsets[si]
            if not is_gt and not wrote_seed and seed is not None:
                # keep scene-0's identity seed at t=0 (off should be 0 for scene 0)
                g.write(reoffset(seed, off) + "\n")
                wrote_seed = True
                total += 1
            for line in data:
                g.write(reoffset(line, off) + "\n")
                total += 1
    return total


def main():
    args = sys.argv[1:]
    seed = 12345
    perturb = False
    basename = None
    keep_temp = False
    pos = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--seed":
            i += 1
            seed = int(args[i])
        elif a == "--perturb-lever":
            perturb = True
        elif a == "--name":
            i += 1
            basename = args[i]
        elif a == "--keep-temp":
            keep_temp = True
        else:
            pos.append(a)
        i += 1
    if len(pos) < 2:
        print("usage: nuscenes_concat.py <dataroot> <out_dir> [scene ...] "
              "[--seed N] [--perturb-lever] [--name BASENAME] [--keep-temp]")
        sys.exit(2)
    dataroot, out_dir = pos[0], pos[1]
    scenes = pos[2:] if len(pos) > 2 else list(ALL_MINI)
    if basename is None:
        basename = "mini_concat" if len(scenes) > 1 else scenes[0] + "_solo"

    os.makedirs(out_dir, exist_ok=True)
    scene_dir = os.path.join(out_dir, "_per_scene")
    os.makedirs(scene_dir, exist_ok=True)

    print("loading nuScenes tables...")
    tables = n2c.load_tables(dataroot)

    # --- per-scene generation (REUSE convert()) ---
    radar_priors = None  # captured from the LAST scene's write_manifest priors (mounts are constant)
    scene_infos = []
    used_scenes = []
    for scene in scenes:
        print("=== generating %s ===" % scene)
        try:
            info = n2c.convert(dataroot, scene_dir, scene, tables, seed=seed)
        except Exception as e:  # noqa: BLE001
            import traceback
            print("ERR %s: %s" % (scene, e))
            traceback.print_exc()
            continue
        scene_infos.append(info)
        used_scenes.append(scene)
    if not used_scenes:
        print("no scenes generated; abort")
        sys.exit(1)

    # capture the real radar priors (yaw pitch roll x y z) from a per-scene manifest -- the mount is
    # CONSTANT across scenes, so read them once (use the first generated scene's manifest).
    radar_priors = read_priors_from_manifest(os.path.join(scene_dir, used_scenes[0] + ".ini"))

    # --- compute per-scene boundary offsets (one shared span per scene) ---
    TICK_GAP_NS = 20_000_000  # 20 ms (one tick at 50 Hz) gap between scenes -> strict monotonicity
    offsets = []
    cum = 0
    spans = []
    for scene in used_scenes:
        offsets.append(cum)
        span = scene_span_ns(scene_dir, scene, SRC_SUFFIXES)
        spans.append(span)
        cum += span + TICK_GAP_NS
    total_span_s = cum / 1e9

    # --- concatenate each source ---
    base = os.path.join(out_dir, basename)
    counts = {}
    for suf in SRC_SUFFIXES:
        outp = "%s_%s.csv" % (base, suf)
        counts[suf] = concat_source(outp, used_scenes, scene_dir, suf,
                                     "CONCAT %s over %d scenes" % (suf, len(used_scenes)),
                                     offsets, is_gt=False)
    # GT
    counts["gt"] = concat_source(base + "_gt.csv", used_scenes, scene_dir, "gt",
                                 "concat gt", offsets, is_gt=True)

    # --- manifest + truth sidecar ---
    write_concat_manifest(out_dir, basename, radar_priors, perturb)
    write_truth(out_dir, basename, radar_priors, perturb)

    # --- report concatenation sanity ---
    print("\n=== CONCAT SUMMARY (%s) ===" % basename)
    print("scenes used: %d (%s)" % (len(used_scenes), ", ".join(used_scenes)))
    print("total megastream span: %.1f s" % total_span_s)
    print("per-scene spans (s): %s" % ", ".join("%.1f" % (s / 1e9) for s in spans))
    print("perturb-lever: %s" % ("ON (+%.2f m, axis varied per radar)" % PERTURB_M if perturb else "OFF (true priors)"))
    print("row counts (megastream, excl headers):")
    for suf in SRC_SUFFIXES + ["gt"]:
        print("  %-20s %d" % (suf, counts[suf]))

    # monotonicity check
    print("monotonicity check:")
    all_ok = True
    for suf in SRC_SUFFIXES + ["gt"]:
        p = "%s_%s.csv" % (base, suf)
        ok, n, first, last = check_monotonic(p)
        all_ok = all_ok and ok
        print("  %-20s %s  rows=%d t[0]=%d t[-1]=%d" % (suf, "OK " if ok else "NONMONO", n, first, last))
    print("ALL MONOTONIC" if all_ok else "*** NON-MONOTONIC ***")

    if not keep_temp:
        shutil.rmtree(scene_dir, ignore_errors=True)
    print("\nmanifest: %s.ini   truth: %s_truth.txt" % (base, base))


def check_monotonic(path):
    last = None
    n = 0
    first_t = None
    ok = True
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            t = int(s.split(",", 1)[0])
            if first_t is None:
                first_t = t
            if last is not None and t <= last:
                ok = False
            last = t
            n += 1
    return ok, n, first_t if first_t is not None else 0, last if last is not None else 0


def read_priors_from_manifest(ini_path):
    # Return {short_name: (sid, yaw, pitch, roll, x, y, z)} parsed from a per-scene .ini.
    priors = {}
    cur_sid = None
    cur_csv = None
    with open(ini_path) as f:
        lines = f.read().splitlines()
    sid = csvname = prior = None
    for line in lines:
        s = line.strip()
        if s.startswith("[sensor."):
            sid = csvname = prior = None
        elif s.startswith("id ="):
            sid = int(s.split("=", 1)[1])
        elif s.startswith("csv ="):
            csvname = s.split("=", 1)[1].strip()
        elif s.startswith("prior_extrinsic ="):
            vals = [float(v) for v in s.split("=", 1)[1].split()]
            prior = vals
        # flush when we hit a blank line after a sensor with a prior
        if (s == "" or line is lines[-1]) and sid is not None and prior is not None and csvname:
            short = _short_from_csv(csvname)
            priors[short] = (sid, prior[0], prior[1], prior[2], prior[3], prior[4], prior[5])
            sid = csvname = prior = None
    # final flush
    if sid is not None and prior is not None and csvname:
        short = _short_from_csv(csvname)
        priors[short] = (sid, prior[0], prior[1], prior[2], prior[3], prior[4], prior[5])
    return priors


def _short_from_csv(csvname):
    # "scene-0061_radar_front_left.csv" -> "radar_front_left"
    base = os.path.basename(csvname)[:-4]  # strip .csv
    idx = base.find("_radar_")
    if idx >= 0:
        return base[idx + 1:]
    return base.split("_", 1)[1]  # wheel/imu


def write_concat_manifest(out_dir, basename, radar_priors, perturb):
    L = []
    L.append("[global]")
    L.append("tick_rate_hz = 50")
    L.append("fusion_delay_s = 0.05")
    L.append("window_s = 0.1")
    L.append("max_sources = 7")
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
    # radars in source-id order (by sid in radar_priors)
    radars_sorted = sorted(radar_priors.items(), key=lambda kv: kv[1][0])
    for short, (sid, yaw, pitch, roll, lx, ly, lz) in radars_sorted:
        plx, ply, plz = lx, ly, lz
        if perturb:
            ax = PERTURB_AXIS.get(short, 0)
            if ax == 0:
                plx += PERTURB_M
            elif ax == 1:
                ply += PERTURB_M
            else:
                plz += PERTURB_M
        L.append("[sensor.%d]" % sid)
        L.append("id = %d" % sid)
        L.append("csv = %s_%s.csv" % (basename, short))
        L.append("form = increment")
        L.append("prior_extrinsic = %.9f %.9f %.9f %.6f %.6f %.6f"
                 % (yaw, pitch, roll, plx, ply, plz))
        L.append("")
    L.append("[gt]")
    L.append("csv = %s_gt.csv" % basename)
    L.append("")
    L.append("[replay]")
    L.append("warmup_steps = 20")
    L.append("local_batch_len = 500")
    L.append("out = %s_out.csv" % basename)
    L.append("")
    with open(os.path.join(out_dir, basename + ".ini"), "w", newline="") as f:
        f.write("\n".join(L))


def write_truth(out_dir, basename, radar_priors, perturb):
    # sidecar listing true vs perturbed lever per radar (the analysis reads this)
    radars_sorted = sorted(radar_priors.items(), key=lambda kv: kv[1][0])
    with open(os.path.join(out_dir, basename + "_truth.txt"), "w", newline="") as f:
        f.write("# radar lever truth/perturbed sidecar for %s\n" % basename)
        f.write("# short, sid, true_x, true_y, true_z, pert_x, pert_y, pert_z, axis, "
                "true_yaw, true_pitch, true_roll\n")
        for short, (sid, yaw, pitch, roll, lx, ly, lz) in radars_sorted:
            plx, ply, plz = lx, ly, lz
            ax = PERTURB_AXIS.get(short, 0)
            if perturb:
                if ax == 0:
                    plx += PERTURB_M
                elif ax == 1:
                    ply += PERTURB_M
                else:
                    plz += PERTURB_M
            f.write("%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.9f,%.9f,%.9f\n"
                    % (short, sid, lx, ly, lz, plx, ply, plz, ax, yaw, pitch, roll))


if __name__ == "__main__":
    main()
