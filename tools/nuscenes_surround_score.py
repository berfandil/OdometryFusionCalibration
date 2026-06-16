#!/usr/bin/env python3
# nuscenes_surround_score.py - score the 6-camera surround multi-extrinsic recovery. Parses an
# ofc_replay summary (the "# final calib" block of a CTRL or RECOVER run) + the scene's true camera
# mounts + the run's _truth.txt sidecar (the perturbed priors), and prints the per-camera
# rotation-extrinsic recovery table: TRUE yaw / PERTURBED-prior yaw / RECOVERED yaw (deg error),
# the full-rotation geodesic error, conf + committed, and the lever recovery.
#
# extr_rot in the summary is so3::log(c.extrinsic.R) -- the ABSOLUTE recovered ego_from_cam rotation
# (the calibrator refines the prior into this), so we exp() it back to R, extract its mount yaw under
# the same Rz(yaw)Ry(pitch)Rx(roll) convention quat_to_ypr uses, and compare to the true mount yaw.
#
# Usage: python tools/nuscenes_surround_score.py <dataroot> <scene> <run_out_summary.txt> [--label L]
#        where run_out_summary.txt is the stdout-captured summary (or the *_out.csv has the calib in
#        a comment block -- we accept either: we scan for the "extr_rot" lines).
import math
import os
import re
import sys

import numpy as np

import nuscenes_to_csv as n2c
import nuscenes_cam_vo as cvo

CAMS = ["CAM_FRONT", "CAM_FRONT_LEFT", "CAM_FRONT_RIGHT",
        "CAM_BACK", "CAM_BACK_LEFT", "CAM_BACK_RIGHT"]
CAM_SRC_BASE = 2


def so3_exp(w):
    return np.array(n2c.so3_exp([w[0], w[1], w[2]]))


def mount_yaw_deg(R):
    # mount yaw under Rz(yaw)Ry(pitch)Rx(roll): yaw = atan2(R10, R00)
    return math.degrees(math.atan2(R[1][0], R[0][0]))


def geodesic_deg(Ra, Rb):
    Rt = np.array(Ra).T @ np.array(Rb)
    c = max(-1.0, min(1.0, (np.trace(Rt) - 1.0) / 2.0))
    return math.degrees(math.acos(c))


def parse_calib(path):
    # Returns {sid: dict(extr_rot=[3], extr_conf, extr_committed, lever=[3], lever_conf,
    #          lever_committed_xyz=[3], scale, scale_conf)}
    out = {}
    with open(path) as f:
        for line in f:
            m = re.search(r"src(\d+)\s+scale=([-\d.eE+]+)\s+\(conf\s+([-\d.eE+]+)", line)
            if not m or "extr_rot" not in line:
                continue
            sid = int(m.group(1))
            scale = float(m.group(2))
            scale_conf = float(m.group(3))
            er = re.search(r"extr_rot\[rx,ry,rz\]=\[([-\d.eE+]+),([-\d.eE+]+),([-\d.eE+]+)\]\s+\(conf\s+([-\d.eE+]+)(,committed)?\)", line)
            lv = re.search(r"lever\[x,y,z\]=\[([-\d.eE+]+),([-\d.eE+]+),([-\d.eE+]+)\]\s+\(conf\s+([-\d.eE+]+)(,committed)?\)", line)
            lc = re.search(r"lever_committed_xyz=\[([01]),([01]),([01])\]", line)
            rec = {"scale": scale, "scale_conf": scale_conf}
            if er:
                rec["extr_rot"] = [float(er.group(1)), float(er.group(2)), float(er.group(3))]
                rec["extr_conf"] = float(er.group(4))
                rec["extr_committed"] = er.group(5) is not None
            if lv:
                rec["lever"] = [float(lv.group(1)), float(lv.group(2)), float(lv.group(3))]
                rec["lever_conf"] = float(lv.group(4))
                rec["lever_committed"] = lv.group(5) is not None
            if lc:
                rec["lever_committed_xyz"] = [int(lc.group(1)), int(lc.group(2)), int(lc.group(3))]
            out[sid] = rec
    return out


def true_mounts(dataroot, scene):
    tables = n2c.load_tables(dataroot)
    _scene, toks = n2c.scene_sample_tokens(tables, scene)
    out = {}
    for ri, ch in enumerate(CAMS):
        sid = CAM_SRC_BASE + ri
        calib, stream = cvo.cam_stream(tables, toks, ch)
        if not stream:
            continue
        qw, qx, qy, qz = calib["rotation"]
        R = n2c.quat_wxyz_to_mat(qw, qx, qy, qz)
        yaw, pitch, roll = n2c.quat_to_ypr(qw, qx, qy, qz)
        lx, ly, lz = calib["translation"]
        out[sid] = {"channel": ch, "R": R, "yaw_deg": math.degrees(yaw),
                    "lever": [lx, ly, lz]}
    return out


def parse_truth_sidecar(path):
    # src2 CAM_FRONT ... prior_yaw=+X rad (+Y deg) ... perturb_axis=z
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            m = re.search(r"src(\d+).*prior_yaw=([-\d.eE+]+)\s+rad", line)
            ax = re.search(r"perturb_axis=(\w)", line)
            if m:
                out[int(m.group(1))] = {"prior_yaw_deg": math.degrees(float(m.group(2))),
                                        "perturb_axis": ax.group(1) if ax else "-"}
    return out


def main():
    if len(sys.argv) < 4:
        print("usage: nuscenes_surround_score.py <dataroot> <scene> <summary.txt> [--label L] [--truth T]")
        sys.exit(2)
    dataroot, scene, summary = sys.argv[1], sys.argv[2], sys.argv[3]
    label = "RUN"
    truth_path = None
    i = 4
    while i < len(sys.argv):
        if sys.argv[i] == "--label":
            i += 1; label = sys.argv[i]
        elif sys.argv[i] == "--truth":
            i += 1; truth_path = sys.argv[i]
        i += 1

    tm = true_mounts(dataroot, scene)
    calib = parse_calib(summary)
    truth = parse_truth_sidecar(truth_path) if truth_path else {}

    print("=" * 116)
    print("CAMERA ROTATION-EXTRINSIC RECOVERY  --  %s  (scene %s)" % (label, scene))
    print("=" * 116)
    print("  %-16s %4s | %9s %9s %9s %9s | %6s %6s %4s | lever recovery" %
          ("camera", "sid", "trueYaw", "priorYaw", "recovYaw", "yawErr", "geoErr", "conf", "cmt"))
    print("  %-16s %4s | %9s %9s %9s %9s | %6s %6s %4s |" %
          ("", "", "(deg)", "(deg)", "(deg)", "(deg)", "(deg)", "", ""))
    print("-" * 116)
    for sid in sorted(tm):
        t = tm[sid]
        c = calib.get(sid, {})
        if "extr_rot" not in c:
            print("  %-16s %4d | (no calib)" % (t["channel"], sid))
            continue
        R_rec = so3_exp(c["extr_rot"])
        recov_yaw = mount_yaw_deg(R_rec.tolist())
        geo = geodesic_deg(t["R"], R_rec.tolist())
        prior_yaw = truth.get(sid, {}).get("prior_yaw_deg", float("nan"))
        ax = truth.get(sid, {}).get("perturb_axis", "-")
        yaw_err = recov_yaw - t["yaw_deg"]
        # normalize yaw_err to [-180,180]
        yaw_err = ((yaw_err + 180.0) % 360.0) - 180.0
        cmt = "Y" if c.get("extr_committed") else "."
        # lever recovery: recovered lever vs true on the perturbed axis
        lev_str = ""
        if "lever" in c:
            axidx = {"x": 0, "y": 1, "z": 2}.get(ax, None)
            lc = c.get("lever_committed_xyz", [0, 0, 0])
            if axidx is not None:
                lev_str = "axis %s: true %.3f recov %.3f (conf %.2f, cmt %s)" % (
                    ax, t["lever"][axidx], c["lever"][axidx], c.get("lever_conf", 0.0),
                    "Y" if lc[axidx] else ".")
            else:
                lev_str = "conf %.2f cmt[%d%d%d]" % (c.get("lever_conf", 0.0), lc[0], lc[1], lc[2])
        print("  %-16s %4d | %+8.2f  %+8.2f  %+8.2f  %+8.2f | %6.2f %6.2f  %-3s | %s" %
              (t["channel"], sid, t["yaw_deg"], prior_yaw, recov_yaw, yaw_err, geo,
               c.get("extr_conf", 0.0), cmt, lev_str))
    print("=" * 116)


if __name__ == "__main__":
    main()
