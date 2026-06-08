#!/usr/bin/env python3
# kaist_to_csv.py - convert KAIST Complex Urban raw navigation data into a 3-source replay
# (see adapters/README.md schema). Builds THREE odometry sources from one drive that share the
# wheel forward speed but differ in their ROTATION source -> a realistic multi-source disagreement
# that exercises the median fusion + reliability + GPS correction on REAL motion:
#   src0 (reference) = wheel differential-drive  (forward + yaw from L/R wheel counts)
#   src1            = wheel forward + FOG 3-axis rotation (KVH FOG incremental angles)
#   src2            = wheel forward + Xsens IMU orientation delta (quaternion)
# + GPS (gps.csv -> GpsCorrection) and GT (global_pose.csv, first-relative).
#
# All sources are emitted in INCREMENT form on the ENCODER timeline (FOG summed / IMU sampled per
# encoder window). Sensors are co-located on the vehicle frame (identity extrinsic) -> this tests
# fusion/robustness/correction, NOT extrinsic/scale calibration.
#
# Reads the per-sequence tarballs (urban<NN>_{data,calibration,pose}.tar.gz) selectively into the
# gitignored kaist_run/raw/<seq>/ (only the small nav CSVs + calib + global_pose; never the LiDAR).
#
# Usage: python tools/kaist_to_csv.py <kaist_root> <out_dir> <seq> [seq ...]
#   e.g. python tools/kaist_to_csv.py C:/workspace/data/KAIST kaist_run urban07 urban12 urban17
import bisect
import math
import os
import sys
import tarfile


# ---------- small matrix / quaternion helpers (no numpy) ----------
def mat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat_t(a):
    return [[a[j][i] for j in range(3)] for i in range(3)]


def rot_xyz(rr, rp, ry):
    cr, sr = math.cos(rr), math.sin(rr)
    cp, sp = math.cos(rp), math.sin(rp)
    cy, sy = math.cos(ry), math.sin(ry)
    Rx = [[1, 0, 0], [0, cr, -sr], [0, sr, cr]]
    Ry = [[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]]
    Rz = [[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]]
    return mat_mul(Rz, mat_mul(Ry, Rx))


def quat_xyzw_to_mat(qx, qy, qz, qw):
    n = math.sqrt(qx * qx + qy * qy + qz * qz + qw * qw) or 1.0
    qx, qy, qz, qw = qx / n, qy / n, qz / n, qw / n
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


# ---------- selective tar extraction (in-memory stream, small files only) ----------
def extract_nav(kaist_root, seq, raw_dir):
    os.makedirs(raw_dir, exist_ok=True)
    need = ["encoder.csv", "fog.csv", "gps.csv", "xsens_imu.csv", "global_pose.csv", "EncoderParameter.txt"]
    if all(os.path.exists(os.path.join(raw_dir, n)) for n in need):
        return   # already extracted (skip the slow 2.6 GB re-stream)
    sd = os.path.join(kaist_root, seq)
    for arch in (seq + "_calibration.tar.gz", seq + "_pose.tar.gz"):
        with tarfile.open(os.path.join(sd, arch)) as t:
            for m in t:
                if m.isfile() and (m.name.endswith(".txt") or m.name.endswith("global_pose.csv")):
                    open(os.path.join(raw_dir, m.name.split("/")[-1]), "wb").write(t.extractfile(m).read())
    want = {"encoder.csv", "fog.csv", "gps.csv", "vrs_gps.csv", "xsens_imu.csv"}
    got = set()
    with tarfile.open(os.path.join(sd, seq + "_data.tar.gz")) as t:
        for m in t:
            bn = m.name.split("/")[-1]
            if m.isfile() and bn in want and bn not in got:
                open(os.path.join(raw_dir, bn), "wb").write(t.extractfile(m).read())
                got.add(bn)
            if got >= want:
                break


def read_rows(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in "#/":
                continue
            out.append([x for x in line.replace(",", " ").split()])
    return out


def enc_params(raw_dir):
    res = dl = dr = base = None
    for line in open(os.path.join(raw_dir, "EncoderParameter.txt")):
        low = line.lower()
        if "resolution" in low:
            res = float(line.split(":")[1])
        elif "left wheel diameter" in low:
            dl = float(line.split(":")[1])
        elif "right wheel diameter" in low:
            dr = float(line.split(":")[1])
        elif "wheel base" in low:
            base = float(line.split(":")[1])
    return res, dl, dr, base


def convert(kaist_root, out_dir, seq):
    raw = os.path.join(out_dir, "raw", seq)
    extract_nav(kaist_root, seq, raw)
    res, dl, dr, base = enc_params(raw)

    enc = [(int(r[0]), float(r[1]), float(r[2])) for r in read_rows(os.path.join(raw, "encoder.csv"))]
    fog = [(int(r[0]), float(r[1]), float(r[2]), float(r[3])) for r in read_rows(os.path.join(raw, "fog.csv"))]
    imu = [(int(r[0]), float(r[1]), float(r[2]), float(r[3]), float(r[4]))
           for r in read_rows(os.path.join(raw, "xsens_imu.csv"))]   # ts, qx,qy,qz,qw
    gps = [(int(r[0]), float(r[1]), float(r[2]), float(r[3])) for r in read_rows(os.path.join(raw, "gps.csv"))]
    gp = read_rows(os.path.join(raw, "global_pose.csv"))             # ts, 12 (3x4 row-major)

    # Align the whole run to the GT (global_pose) span: GT can start well after the encoder
    # (SLAM init delay, ~15 s in urban07). Use the GT start as t0 and clip sources/GPS to
    # [t0, tend] so the harness's GT anchor + NEES compare like-for-like.
    t0 = int(gp[0][0])
    tend = int(gp[-1][0])
    enc = [e for e in enc if t0 <= e[0] <= tend]
    fog_ts = [f[0] for f in fog]
    imu_ts = [i[0] for i in imu]

    def imu_R_at(ts):
        k = bisect.bisect_left(imu_ts, ts)
        k = min(max(k, 0), len(imu) - 1)
        _, qx, qy, qz, qw = imu[k]
        return quat_xyzw_to_mat(qx, qy, qz, qw)

    os.makedirs(out_dir, exist_ok=True)
    base_out = os.path.join(out_dir, seq)
    f0 = open(base_out + "_wheel.csv", "w")     # src0
    f1 = open(base_out + "_wheelfog.csv", "w")  # src1
    f2 = open(base_out + "_wheelimu.csv", "w")  # src2
    for f in (f0, f1, f2):
        f.write("# form: increment\n# t_ns, x, y, z, qw, qx, qy, qz\n")

    def emit(f, t_ns, R, tx):
        qw, qx, qy, qz = mat_to_quat(R)
        f.write("%d,%.6f,0,0,%.9f,%.9f,%.9f,%.9f\n" % (t_ns, tx, qw, qx, qy, qz))

    # First row = identity seed (per the increment-form contract).
    for f in (f0, f1, f2):
        f.write("%d,0,0,0,1,0,0,0\n" % 0)

    prev_imu_R = imu_R_at(enc[0][0])
    for k in range(1, len(enc)):
        tp, lp, rp = enc[k - 1]
        tc, lc, rc = enc[k]
        if tc <= tp:
            continue
        rel = tc - t0
        dleft = ((lc - lp) / res) * math.pi * dl
        dright = ((rc - rp) / res) * math.pi * dr
        d_center = 0.5 * (dleft + dright)
        # src0: wheel differential-drive yaw
        d_yaw = (dright - dleft) / base
        emit(f0, rel, rot_xyz(0, 0, d_yaw), d_center)
        # src1: FOG 3-axis incremental rotation summed over (tp, tc]
        i0 = bisect.bisect_right(fog_ts, tp)
        i1 = bisect.bisect_right(fog_ts, tc)
        sr = sp = sy = 0.0
        for j in range(i0, i1):
            sr += fog[j][1]; sp += fog[j][2]; sy += fog[j][3]
        emit(f1, rel, rot_xyz(sr, sp, sy), d_center)
        # src2: Xsens IMU orientation delta (body frame) R_prev^T R_cur
        cur_imu_R = imu_R_at(tc)
        dR = mat_mul(mat_t(prev_imu_R), cur_imu_R)
        emit(f2, rel, dR, d_center)
        prev_imu_R = cur_imu_R
    for f in (f0, f1, f2):
        f.close()

    # GPS CSV (geodetic -> GpsCorrection). pos var unknown here; use a modest 1 m^2 / 4 m^2 (up).
    with open(base_out + "_gps.csv", "w") as fg:
        fg.write("# t_ns, lat_deg, lon_deg, alt_m, var_e, var_n, var_u\n")
        for ts, lat, lon, alt in gps:
            if t0 <= ts <= tend:
                fg.write("%d,%.10f,%.10f,%.6f,1,1,4\n" % (ts - t0, lat, lon, alt))

    # GT: global_pose 3x4 [R|t] (UTM) -> first-relative (t0-vehicle frame).
    R0 = [[float(gp[0][1 + 4 * i + j]) for j in range(3)] for i in range(3)]
    t0w = [float(gp[0][4]), float(gp[0][8]), float(gp[0][12])]
    R0t = mat_t(R0)
    yaw0 = math.atan2(R0[1][0], R0[0][0])
    with open(base_out + "_gt.csv", "w") as fgt:
        fgt.write("# t_ns, x, y, z, qw, qx, qy, qz\n")
        for r in gp:
            Ri = [[float(r[1 + 4 * i + j]) for j in range(3)] for i in range(3)]
            ti = [float(r[4]), float(r[8]), float(r[12])]
            Rrel = mat_mul(R0t, Ri)
            d = [ti[0] - t0w[0], ti[1] - t0w[1], ti[2] - t0w[2]]
            trel = [sum(R0t[a][b] * d[b] for b in range(3)) for a in range(3)]
            qw, qx, qy, qz = mat_to_quat(Rrel)
            fgt.write("%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f,%.9f\n"
                      % (int(r[0]) - t0, trel[0], trel[1], trel[2], qw, qx, qy, qz))

    write_manifest(out_dir, seq, yaw0)
    dur = (enc[-1][0] - t0) / 1e9
    return {"enc": len(enc), "fog": len(fog), "imu": len(imu), "gps": len(gps), "dur_s": dur, "yaw0": yaw0}


def write_manifest(out_dir, seq, yaw0):
    txt = (
        "[global]\n"
        "tick_rate_hz = 100\n"
        "fusion_delay_s = 0.05\n"
        "window_s = 0.05\n"
        "max_sources = 3\n"
        "reference_sensor_id = 0\n"
        "cold_start = median_from_start\n"
        "timesync_enabled = false\n"
        "q_floor = 0.001 0.001 0.001 0.0001 0.0001 0.0001\n\n"
        "[sensor.0]\nid = 0\nis_reference = true\ncsv = %s_wheel.csv\nform = increment\n\n"
        "[sensor.1]\nid = 1\ncsv = %s_wheelfog.csv\nform = increment\n\n"
        "[sensor.2]\nid = 2\ncsv = %s_wheelimu.csv\nform = increment\n\n"
        "[gps]\ncsv = %s_gps.csv\nodom_from_enu_yaw = %.9f\nlever_x = -0.32\nlever_z = 1.7\ncov_floor_m2 = 0.0\n\n"
        "[gt]\ncsv = %s_gt.csv\n\n"
        "[replay]\nwarmup_steps = 20\n"
        % (seq, seq, seq, seq, -yaw0, seq)
    )
    open(os.path.join(out_dir, seq + ".ini"), "w").write(txt)


def main():
    if len(sys.argv) < 4:
        print("usage: kaist_to_csv.py <kaist_root> <out_dir> <seq> [seq ...]")
        sys.exit(2)
    root, out_dir = sys.argv[1], sys.argv[2]
    for seq in sys.argv[3:]:
        try:
            info = convert(root, out_dir, seq)
            print("OK  %-8s enc=%-6d fog=%-7d imu=%-6d gps=%-5d dur=%6.1fs yaw0=%+.3f"
                  % (seq, info["enc"], info["fog"], info["imu"], info["gps"], info["dur_s"], info["yaw0"]))
        except Exception as e:  # noqa: BLE001
            import traceback
            print("ERR %-8s %s" % (seq, e)); traceback.print_exc()


if __name__ == "__main__":
    main()
