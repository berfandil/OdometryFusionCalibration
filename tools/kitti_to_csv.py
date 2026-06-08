#!/usr/bin/env python3
# kitti_to_csv.py - convert KITTI raw OXTS (GPS/IMU/INS) into the OdometryFusionCalibration
# replay CSV schema (see adapters/README.md) + an ofc_replay manifest, reading the *_sync.zip
# files IN MEMORY (nothing is unzipped to disk).
#
# Per KITTI drive (a folder under <kitti_root> containing <date>_drive_NNNN_sync.zip) it emits,
# into <out_dir>:
#   <drive>_imu.csv  - twist source: t_ns, vf, vl, vu, wf, wl, wu   (vehicle/body frame, x-fwd y-left z-up)
#   <drive>_gps.csv  - absolute ref: t_ns, lat_deg, lon_deg, alt_m, var_e, var_n, var_u
#   <drive>_gt.csv   - GT track    : t_ns, x, y, z, qw, qx, qy, qz  (convertOxtsToPose, FIRST-relative)
#   <drive>_gps.ini  - manifest WITH the GPS correction
#   <drive>_nogps.ini- manifest WITHOUT GPS (predict-only: isolates velocity-integration drift)
#
# Frames: KITTI world is ENU (Mercator x=East, y=North, z=Up). GT is made FIRST-relative
# (T0_i = inv(T0).Ti), so it lives in the t0-vehicle frame -- the SAME frame the estimator
# integrates body twist in (body starts at identity). odom_from_enu_yaw = -yaw0 rotates the
# GPS ENU displacement into that t0-vehicle frame. (Flat-ground pitch/roll-at-t0 ignored; the
# manifest only exposes a yaw alignment.)
#
# OXTS twist (vf,vl,vu,wf,wl,wu) is the body velocity/rate -> a ready-made twist source. NOTE:
# OXTS GPS, IMU and the convertOxtsToPose GT all come from the SAME RT3003 INS, so this is a
# single-odometry-source SMOKE TEST (predict + GPS-correction + NEES vs an absolute-position GT),
# NOT a multi-source fusion test.
#
# Usage: python kitti_to_csv.py <kitti_root> <out_dir> [drive_glob]
#   e.g. python tools/kitti_to_csv.py C:/workspace/data/KITTI kitti_run
import math
import os
import sys
import zipfile

A_EARTH = 6378137.0  # KITTI Mercator earth radius (m)
DEG = math.pi / 180.0

# OXTS column indices (0-based) from oxts/dataformat.txt
C_LAT, C_LON, C_ALT, C_ROLL, C_PITCH, C_YAW = 0, 1, 2, 3, 4, 5
C_VF, C_VL, C_VU = 8, 9, 10
C_WF, C_WL, C_WU = 20, 21, 22
C_POS_ACC = 23


def tod_ns(line):
    # "2011-09-26 13:02:25.964389445" -> nanoseconds-of-day (monotonic within a drive)
    tod = line.strip().split(" ")[1]
    hh, mm, ss = tod.split(":")
    s_int, s_frac = (ss.split(".") + ["0"])[:2]
    return (((int(hh) * 3600 + int(mm) * 60 + int(s_int)) * 1_000_000_000)
            + int(s_frac.ljust(9, "0")[:9]))


def mat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat_t(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def mat_vec(a, v):
    return [sum(a[i][k] * v[k] for k in range(3)) for i in range(3)]


def rot_from_rpy(roll, pitch, yaw):
    # convertOxtsToPose: R = Rz(yaw) @ Ry(pitch) @ Rx(roll)
    cr, sr = math.cos(roll), math.sin(roll)
    cp, sp = math.cos(pitch), math.sin(pitch)
    cy, sy = math.cos(yaw), math.sin(yaw)
    Rx = [[1, 0, 0], [0, cr, -sr], [0, sr, cr]]
    Ry = [[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]]
    Rz = [[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]]
    return mat_mul(Rz, mat_mul(Ry, Rx))


def mat_to_quat(R):
    # w-first quaternion from a rotation matrix (Shepperd's method, stable)
    tr = R[0][0] + R[1][1] + R[2][2]
    if tr > 0:
        s = math.sqrt(tr + 1.0) * 2.0
        w = 0.25 * s
        x = (R[2][1] - R[1][2]) / s
        y = (R[0][2] - R[2][0]) / s
        z = (R[1][0] - R[0][1]) / s
    elif R[0][0] > R[1][1] and R[0][0] > R[2][2]:
        s = math.sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2.0
        w = (R[2][1] - R[1][2]) / s
        x = 0.25 * s
        y = (R[0][1] + R[1][0]) / s
        z = (R[0][2] + R[2][0]) / s
    elif R[1][1] > R[2][2]:
        s = math.sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2.0
        w = (R[0][2] - R[2][0]) / s
        x = (R[0][1] + R[1][0]) / s
        y = 0.25 * s
        z = (R[1][2] + R[2][1]) / s
    else:
        s = math.sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2.0
        w = (R[1][0] - R[0][1]) / s
        x = (R[0][2] + R[2][0]) / s
        y = (R[1][2] + R[2][1]) / s
        z = 0.25 * s
    n = math.sqrt(w * w + x * x + y * y + z * z) or 1.0
    return (w / n, x / n, y / n, z / n)


def load_oxts(zip_path):
    """Read oxts/{timestamps.txt, data/*.txt} from a *_sync.zip IN MEMORY. Returns list of
    (t_ns, row[floats]) sorted by frame index, t_ns made relative to the first frame."""
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        ts_name = next(n for n in names if n.endswith("/oxts/timestamps.txt"))
        data_names = sorted(n for n in names if "/oxts/data/" in n and n.endswith(".txt"))
        stamps = z.read(ts_name).decode().splitlines()
        rows = []
        for i, dn in enumerate(data_names):
            vals = [float(x) for x in z.read(dn).decode().split()]
            rows.append((tod_ns(stamps[i]), vals))
    t0 = rows[0][0]
    return [(t - t0, v) for (t, v) in rows]


def convert_drive(zip_path, drive, out_dir):
    rows = load_oxts(zip_path)
    if len(rows) < 2:
        return None
    lat0 = rows[0][1][C_LAT]
    scale = math.cos(lat0 * DEG)
    yaw0 = rows[0][1][C_YAW]

    def merc(lat, lon):
        mx = scale * lon * math.pi * A_EARTH / 180.0
        my = scale * A_EARTH * math.log(math.tan((90.0 + lat) * math.pi / 360.0))
        return mx, my

    # World (ENU/Mercator) pose of frame 0, for the first-relative transform.
    mx0, my0 = merc(rows[0][1][C_LAT], rows[0][1][C_LON])
    t0_w = [mx0, my0, rows[0][1][C_ALT]]
    R0 = rot_from_rpy(rows[0][1][C_ROLL], rows[0][1][C_PITCH], rows[0][1][C_YAW])
    R0t = mat_t(R0)

    base = os.path.join(out_dir, drive)
    with open(base + "_imu.csv", "w") as f_imu, \
         open(base + "_gps.csv", "w") as f_gps, \
         open(base + "_gt.csv", "w") as f_gt:
        f_imu.write("# form: twist\n# t_ns, vf, vl, vu, wf, wl, wu\n")
        f_gps.write("# t_ns, lat_deg, lon_deg, alt_m, var_e, var_n, var_u\n")
        f_gt.write("# t_ns, x, y, z, qw, qx, qy, qz\n")
        for t_ns, v in rows:
            f_imu.write("%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n"
                        % (t_ns, v[C_VF], v[C_VL], v[C_VU], v[C_WF], v[C_WL], v[C_WU]))
            pa = max(v[C_POS_ACC], 0.05)
            f_gps.write("%d,%.10f,%.10f,%.6f,%.6g,%.6g,%.6g\n"
                        % (t_ns, v[C_LAT], v[C_LON], v[C_ALT], pa * pa, pa * pa, (2 * pa) ** 2))
            # GT: world pose -> first-relative (t0-vehicle frame)
            mx, my = merc(v[C_LAT], v[C_LON])
            Ri = rot_from_rpy(v[C_ROLL], v[C_PITCH], v[C_YAW])
            Rrel = mat_mul(R0t, Ri)
            dt = [mx - t0_w[0], my - t0_w[1], v[C_ALT] - t0_w[2]]
            trel = mat_vec(R0t, dt)
            qw, qx, qy, qz = mat_to_quat(Rrel)
            f_gt.write("%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f,%.9f\n"
                       % (t_ns, trel[0], trel[1], trel[2], qw, qx, qy, qz))

    dur_s = rows[-1][0] / 1e9
    write_manifests(out_dir, drive, yaw0)
    return {"frames": len(rows), "dur_s": dur_s, "yaw0": yaw0}


def write_manifests(out_dir, drive, yaw0):
    glob_block = (
        "[global]\n"
        "tick_rate_hz = 10\n"
        "fusion_delay_s = 0.10\n"
        "window_s = 0.10\n"
        "max_sources = 1\n"
        "reference_sensor_id = 0\n"
        "cold_start = median_from_start\n"
        "timesync_enabled = false\n\n"
        "[sensor.0]\n"
        "id = 0\n"
        "is_reference = true\n"
        "csv = %s_imu.csv\n"
        "form = twist\n\n" % drive
    )
    gt_block = "[gt]\ncsv = %s_gt.csv\n\n" % drive
    gps_block = (
        "[gps]\n"
        "csv = %s_gps.csv\n"
        "odom_from_enu_yaw = %.9f\n"   # -yaw0: ENU -> t0-vehicle frame
        "cov_floor_m2 = 0.0\n\n" % (drive, -yaw0)
    )
    with open(os.path.join(out_dir, drive + "_nogps.ini"), "w") as f:
        f.write(glob_block + gt_block + "[replay]\nwarmup_steps = 5\n")
    with open(os.path.join(out_dir, drive + "_gps.ini"), "w") as f:
        f.write(glob_block + gps_block + gt_block + "[replay]\nwarmup_steps = 5\n")


def main():
    if len(sys.argv) < 3:
        print("usage: kitti_to_csv.py <kitti_root> <out_dir> [drive_substring]")
        sys.exit(2)
    root, out_dir = sys.argv[1], sys.argv[2]
    filt = sys.argv[3] if len(sys.argv) > 3 else ""
    os.makedirs(out_dir, exist_ok=True)
    zips = []
    for dirpath, _, files in os.walk(root):
        for fn in files:
            if fn.endswith("_sync.zip") and (filt in fn):
                zips.append(os.path.join(dirpath, fn))
    zips.sort()
    print("found %d sync zips" % len(zips))
    for zp in zips:
        drive = os.path.basename(zp)[: -len("_sync.zip")]  # 2011_09_26_drive_0001
        try:
            info = convert_drive(zp, drive, out_dir)
            print("OK  %-28s frames=%-5d dur=%5.1fs yaw0=%+.3f"
                  % (drive, info["frames"], info["dur_s"], info["yaw0"]))
        except Exception as e:  # noqa: BLE001 - report + continue
            print("ERR %-28s %s" % (drive, e))


if __name__ == "__main__":
    main()
