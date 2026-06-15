#!/usr/bin/env python3
# kaist_surround_synth.py - synthesize SURROUND automotive odometry sources (4 cameras + 2 radars)
# from a KAIST ground-truth absolute-pose CSV, for the MANY-SOURCE fusion test (9 sources at once:
# 3 REAL KAIST streams + 6 SYNTHESIZED surround streams). See adapters/README.md for the CSV/manifest
# schema and tools/inject_calib.py for the extrinsic-conjugation convention this tool reuses EXACTLY.
#
# WHAT IT DOES -----------------------------------------------------------------------------------
# Given a GT abs-pose CSV (t_ns, x,y,z, qw,qx,qy,qz), it:
#   1. forms base-frame per-step increments  delta_k = T_{k-1}^-1 o T_k  from consecutive GT poses;
#   2. for each synthetic sensor in the hardcoded SENSORS table, CONJUGATES each delta by the
#      sensor's planted extrinsic E (E = reference_from_sensor, i.e. sensor->base):
#         delta_sensor = E^-1 o delta_ref o E             (the inject_calib.py adjoint/conjugation)
#      then SCALES its translation by the planted scale, then (unless --noise-free) adds zero-mean
#      Gaussian noise on the se3 tangent (per-sensor sigma_trans, sigma_rot) via a right-perturbation
#         delta_noisy = delta_sensor o exp(xi),  xi = [n_t(3); n_r(3)] ~ N(0, diag(sigma^2));
#   3. writes an INCREMENT-form CSV per sensor WITH the 6 per-delta variance columns
#         t_ns, x,y,z, qw,qx,qy,qz, var_x,var_y,var_z, var_rx,var_ry,var_rz
#      set to the sensor's MODELED per-step variances (so confidence weighting + NEES are meaningful:
#      cameras = heading-decent -> small rot variance; radars = heading-weak -> LARGE rot variance).
#
# CONVENTION (must match config_loader parse_extrinsic + inject_calib) ----------------------------
#   E.R = Rz(yaw) Ry(pitch) Rx(roll);  E.t = lever[x,y,z]  (metres).
#   The fusion frame-aligns a source back via  A = X o B o X^-1,  X = prior_extrinsic, de-scaling
#   B.t/prior_scale first. So feeding prior_extrinsic = the planted E (+ prior_scale = planted scale)
#   re-aligns a noise-free synth source to the base EXACTLY (zero residual) -> the correctness pin.
#
# ANGLE UNITS: this tool's SENSORS table is in DEGREES (human-readable, matches inject_calib's CLI).
#   config_loader parse_extrinsic reads RADIANS, so --emit-manifest prints the prior_extrinsic line
#   in RADIANS (yaw pitch roll x y z). Keep the two consistent.
#
# NOISE NUMBERS (principled; per-step, 100 Hz encoder cadence so a "step" ~ 1 cm-1 m of motion) -----
#   cam_*  (metric surround visual odometry, full 6DOF, good heading):
#       sigma_trans = 0.5% of the per-step translation magnitude + 1e-3 m floor   (good VO scale)
#       sigma_rot   = 0.15 deg/step                                               (heading-decent)
#   radar_* (automotive ego-velocity radar: strong translation, weak heading):
#       sigma_trans = 1.0% of the per-step translation magnitude + 2e-3 m floor   (decent forward)
#       sigma_rot   = 1.5 deg/step                                                (heading-WEAK)
#   The variance COLUMNS use a representative per-step translation (a fixed nominal so the columns
#   are constant + comparable): var_trans = (pct*NOMINAL_STEP_M + floor)^2 ; var_rot = (sigma_rot)^2
#   in rad^2. NOMINAL_STEP_M defaults to 0.05 m (~5 m/s at 100 Hz) and is overridable via --nominal.
#   This makes radars' rot variance ~100x the cameras' -> the split-channel median down-weights radar
#   heading while keeping its translation, and the heading monitor should NOT pick a radar as winner.
#
# Pure ASCII, stdlib only, SEEDED RNG (Mersenne via random.Random(seed)); timestamps + seeds passed
# explicitly, no time-based randomness -> reproducible.
#
# Usage:
#   python tools/kaist_surround_synth.py GT.csv OUT_PREFIX [--noise-free] [--seed N] [--nominal M]
#                                        [--emit-manifest]
#   e.g. python tools/kaist_surround_synth.py kaist_run/urban12_gt.csv kaist_run/urban12 --seed 12
# Emits OUT_PREFIX_<sensorname>.csv for each sensor (e.g. urban12_cam_front.csv ...). With
# --emit-manifest, also prints the [sensor.N] prior_extrinsic/prior_scale lines (radians) to stdout.
import math
import random
import sys


# ---------- minimal SO3 / SE3 in plain lists (mirrors inject_calib.py) ----------
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def matvec(A, v):
    return [sum(A[i][k] * v[k] for k in range(3)) for i in range(3)]


def transpose(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]


def Rx(a):
    c, s = math.cos(a), math.sin(a)
    return [[1, 0, 0], [0, c, -s], [0, s, c]]


def Ry(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, 0, s], [0, 1, 0], [-s, 0, c]]


def Rz(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]


def quat_to_R(qw, qx, qy, qz):
    n = math.sqrt(qw * qw + qx * qx + qy * qy + qz * qz) or 1.0
    qw, qx, qy, qz = qw / n, qx / n, qy / n, qz / n
    return [[1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - qz * qw), 2 * (qx * qz + qy * qw)],
            [2 * (qx * qy + qz * qw), 1 - 2 * (qx * qx + qz * qz), 2 * (qy * qz - qx * qw)],
            [2 * (qx * qz - qy * qw), 2 * (qy * qz + qx * qw), 1 - 2 * (qx * qx + qy * qy)]]


def R_to_quat(R):
    tr = R[0][0] + R[1][1] + R[2][2]
    if tr > 0:
        s = math.sqrt(tr + 1.0) * 2
        qw = 0.25 * s
        qx = (R[2][1] - R[1][2]) / s
        qy = (R[0][2] - R[2][0]) / s
        qz = (R[1][0] - R[0][1]) / s
    elif R[0][0] > R[1][1] and R[0][0] > R[2][2]:
        s = math.sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2
        qw = (R[2][1] - R[1][2]) / s
        qx = 0.25 * s
        qy = (R[0][1] + R[1][0]) / s
        qz = (R[0][2] + R[2][0]) / s
    elif R[1][1] > R[2][2]:
        s = math.sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2
        qw = (R[0][2] - R[2][0]) / s
        qx = (R[0][1] + R[1][0]) / s
        qy = 0.25 * s
        qz = (R[1][2] + R[2][1]) / s
    else:
        s = math.sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2
        qw = (R[1][0] - R[0][1]) / s
        qx = (R[0][2] + R[2][0]) / s
        qy = (R[1][2] + R[2][1]) / s
        qz = 0.25 * s
    n = math.sqrt(qw * qw + qx * qx + qy * qy + qz * qz) or 1.0
    return qw / n, qx / n, qy / n, qz / n


# SE3 compose: (R3, t3) = (R1,t1) o (R2,t2)  ->  R3 = R1 R2, t3 = R1 t2 + t1.
def se3_mul(R1, t1, R2, t2):
    return matmul(R1, R2), [matvec(R1, t2)[i] + t1[i] for i in range(3)]


def se3_inv(R, t):
    Rt = transpose(R)
    return Rt, [-x for x in matvec(Rt, t)]


# SO3 exp of a small rotation vector w (rad) -> R (Rodrigues).
def so3_exp(w):
    th = math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
    if th < 1e-12:
        # first-order: I + [w]x
        return [[1, -w[2], w[1]], [w[2], 1, -w[0]], [-w[1], w[0], 1]]
    k = [w[0] / th, w[1] / th, w[2] / th]
    c, s = math.cos(th), math.sin(th)
    # K = [k]x ; R = I + sin(th) K + (1-cos(th)) K^2
    K = [[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]]
    K2 = matmul(K, K)
    return [[(1.0 if i == j else 0.0) + s * K[i][j] + (1 - c) * K2[i][j]
             for j in range(3)] for i in range(3)]


# se3 exp of a tangent xi = [trans(3); rot(3)] -> (R, t). Right-perturbation of an increment:
# we use the standard se3 exp so the rotation noise couples into translation via the left-Jacobian.
def se3_exp(xi):
    rho = [xi[0], xi[1], xi[2]]
    phi = [xi[3], xi[4], xi[5]]
    th = math.sqrt(phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2])
    R = so3_exp(phi)
    if th < 1e-12:
        V = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    else:
        k = [phi[0] / th, phi[1] / th, phi[2] / th]
        K = [[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]]
        K2 = matmul(K, K)
        a = (1 - math.cos(th)) / th
        b = (th - math.sin(th)) / th
        V = [[(1.0 if i == j else 0.0) + a * K[i][j] + b * K2[i][j]
              for j in range(3)] for i in range(3)]
    t = matvec(V, rho)
    return R, t


# ---------- SENSOR TABLE (surround automotive; E = reference_from_sensor) ----------
# name, yaw_deg, pitch_deg, roll_deg, lever[x,y,z] (m), scale,
#   sigma_trans_pct (fraction of per-step |t|), sigma_trans_floor_m, sigma_rot_deg (per step)
SENSORS = [
    ("cam_front",  0.0,   0.0, 0.0, [1.5,  0.0, 1.2], 1.0, 0.005, 1e-3, 0.15),
    ("cam_right", -90.0,  0.0, 0.0, [0.0, -0.9, 1.0], 1.0, 0.005, 1e-3, 0.15),
    ("cam_rear",   180.0, 0.0, 0.0, [-1.5, 0.0, 1.2], 1.0, 0.005, 1e-3, 0.15),
    ("cam_left",   90.0,  0.0, 0.0, [0.0,  0.9, 1.0], 1.0, 0.005, 1e-3, 0.15),
    ("radar_fl",  -45.0,  0.0, 0.0, [1.6,  0.7, 0.5], 1.0, 0.010, 2e-3, 1.50),
    ("radar_fr",   45.0,  0.0, 0.0, [1.6, -0.7, 0.5], 1.0, 0.010, 2e-3, 1.50),
]


def extrinsic_E(yaw_deg, pitch_deg, roll_deg, lever):
    y, p, r = math.radians(yaw_deg), math.radians(pitch_deg), math.radians(roll_deg)
    ER = matmul(matmul(Rz(y), Ry(p)), Rx(r))
    return ER, list(lever)


def read_gt(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#/":
                continue
            p = s.replace(",", " ").split()
            t = int(p[0])
            x, y, z = float(p[1]), float(p[2]), float(p[3])
            qw, qx, qy, qz = float(p[4]), float(p[5]), float(p[6]), float(p[7])
            rows.append((t, x, y, z, qw, qx, qy, qz))
    return rows


def main():
    args = sys.argv[1:]
    flags = {"noise_free": False, "seed": 12345, "nominal": 0.05, "emit_manifest": False,
             "noise_scale": 1.0}
    pos = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--noise-free":
            flags["noise_free"] = True
        elif a == "--emit-manifest":
            flags["emit_manifest"] = True
        elif a == "--seed":
            i += 1
            flags["seed"] = int(args[i])
        elif a == "--nominal":
            i += 1
            flags["nominal"] = float(args[i])
        elif a == "--noise-scale":
            i += 1
            flags["noise_scale"] = float(args[i])
        else:
            pos.append(a)
        i += 1
    if len(pos) < 2:
        print("usage: kaist_surround_synth.py GT.csv OUT_PREFIX "
              "[--noise-free] [--seed N] [--nominal M] [--emit-manifest]")
        sys.exit(2)
    gt_path, out_prefix = pos[0], pos[1]

    gt = read_gt(gt_path)
    if len(gt) < 2:
        print("error: GT has < 2 rows")
        sys.exit(1)

    # Precompute base-frame increments delta_k = T_{k-1}^-1 o T_k (k=1..N-1), keyed to t_k.
    deltas = []  # (t_ns, Rd, td)
    Rprev = quat_to_R(gt[0][4], gt[0][5], gt[0][6], gt[0][7])
    tprev = [gt[0][1], gt[0][2], gt[0][3]]
    for k in range(1, len(gt)):
        t = gt[k][0]
        Rc = quat_to_R(gt[k][4], gt[k][5], gt[k][6], gt[k][7])
        tc = [gt[k][1], gt[k][2], gt[k][3]]
        Rpi, tpi = se3_inv(Rprev, tprev)
        Rd, td = se3_mul(Rpi, tpi, Rc, tc)
        deltas.append((t, Rd, td))
        Rprev, tprev = Rc, tc

    nominal = flags["nominal"]
    ns = flags["noise_scale"]
    for (name, yaw, pitch, roll, lever, scale, s_pct, s_floor, s_rot_deg) in SENSORS:
        ER, Et = extrinsic_E(yaw, pitch, roll, lever)
        ERt, Einv_t = se3_inv(ER, Et)
        # --noise-scale uniformly multiplies all per-step sigmas (and so the variance columns
        # by ns^2). The TABLE values are a representative-but-pessimistic per-step level; the
        # physical heading random walk is sigma_step * sqrt(N_steps), so a long 100 Hz drive
        # needs small per-step sigma to stay GPS-correctable. ns < 1 dials toward realistic VO.
        s_pct *= ns
        s_floor *= ns
        s_rot = math.radians(s_rot_deg) * ns
        # Modeled per-step variance columns (constant; representative nominal step).
        var_t = (s_pct * nominal + s_floor) ** 2
        var_r = s_rot ** 2
        # Per-sensor seeded RNG (derive a distinct, stable stream per sensor name for
        # independence; sum-of-ordinals keeps it int + reproducible across runs/platforms).
        name_salt = sum((idx + 1) * ord(ch) for idx, ch in enumerate(name))
        rng = random.Random(flags["seed"] * 1000003 + name_salt)
        out_path = "%s_%s.csv" % (out_prefix, name)
        n_out = 0
        with open(out_path, "w", newline="") as g:
            g.write("# form: increment\n")
            g.write("# t_ns, x, y, z, qw, qx, qy, qz, var_x, var_y, var_z, var_rx, var_ry, var_rz\n")
            g.write("# synth surround sensor '%s'  E=yaw%.4f pitch%.4f roll%.4f (deg) "
                    "lever%s scale=%.4f\n" % (name, yaw, pitch, roll, lever, scale))
            g.write("# noise: sigma_trans=%.3f%%*|t|+%.4g m  sigma_rot=%.4g deg/step%s\n"
                    % (s_pct * 100.0, s_floor, s_rot_deg,
                       "  (NOISE-FREE)" if flags["noise_free"] else ""))
            # First row = identity seed (increment-form contract).
            g.write("0,0,0,0,1,0,0,0,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n"
                    % (var_t, var_t, var_t, var_r, var_r, var_r))
            n_out += 1
            for (t, Rd, td) in deltas:
                # delta_sensor = E^-1 o delta_ref o E  (conjugation; the inject_calib math).
                R1, t1 = se3_mul(ERt, Einv_t, Rd, td)
                Rs, ts = se3_mul(R1, t1, ER, Et)
                ts = [c * scale for c in ts]
                if not flags["noise_free"]:
                    # per-step translation magnitude for the percent term
                    mag = math.sqrt(ts[0] * ts[0] + ts[1] * ts[1] + ts[2] * ts[2])
                    st = s_pct * mag + s_floor
                    xi = [rng.gauss(0.0, st), rng.gauss(0.0, st), rng.gauss(0.0, st),
                          rng.gauss(0.0, s_rot), rng.gauss(0.0, s_rot), rng.gauss(0.0, s_rot)]
                    Rn, tn = se3_exp(xi)
                    Rs, ts = se3_mul(Rs, ts, Rn, tn)
                ow, ox, oy, oz = R_to_quat(Rs)
                g.write("%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
                        "%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n"
                        % (t, ts[0], ts[1], ts[2], ow, ox, oy, oz,
                           var_t, var_t, var_t, var_r, var_r, var_r))
                n_out += 1
        print("wrote %s  (%d rows)  E=yaw%.3f deg lever%s scale=%.3f var_t=%.3g var_r=%.3g%s"
              % (out_path, n_out, yaw, lever, scale, var_t, var_r,
                 "  noise-free" if flags["noise_free"] else ""))

    if flags["emit_manifest"]:
        print("\n# --- prior_extrinsic / prior_scale lines (RADIANS, yaw pitch roll x y z) ---")
        for (name, yaw, pitch, roll, lever, scale, _sp, _sf, _sr) in SENSORS:
            print("# %s:" % name)
            print("prior_extrinsic = %.9f %.9f %.9f %.6f %.6f %.6f"
                  % (math.radians(yaw), math.radians(pitch), math.radians(roll),
                     lever[0], lever[1], lever[2]))
            print("prior_scale = %.6f" % scale)


if __name__ == "__main__":
    main()
