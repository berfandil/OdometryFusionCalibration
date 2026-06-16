#!/usr/bin/env python
# Per-drive MOTION STATISTICS from a KAIST GT track (absolute pose CSV).
# Pure ASCII, numpy, deterministic (no randomness).
#
# Reads "# t_ns, x, y, z, qw, qx, qy, qz" rows, computes per-step
# (consecutive GT samples) translation magnitude and rotation-angle
# magnitude, plus aggregate motion statistics used to test the
# motion-aware-Q hypothesis (does the honest predict-only q_floor track
# a motion stat?).
import sys
import numpy as np


def load_gt(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if len(parts) < 8:
                continue
            rows.append([float(p) for p in parts[:8]])
    a = np.array(rows, dtype=float)
    t = a[:, 0] * 1e-9  # ns -> s
    p = a[:, 1:4]
    q = a[:, 4:8]  # qw,qx,qy,qz
    return t, p, q


def quat_to_R(q):
    # q = [qw,qx,qy,qz], normalized
    w, x, y, z = q
    n = np.sqrt(w * w + x * x + y * y + z * z)
    if n == 0:
        return np.eye(3)
    w, x, y, z = w / n, x / n, y / n, z / n
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def rot_angle(Ra, Rb):
    # geodesic angle between two rotations = ||log(Ra^T Rb)||
    Rrel = Ra.T @ Rb
    c = (np.trace(Rrel) - 1.0) * 0.5
    c = max(-1.0, min(1.0, c))
    return np.arccos(c)


def yaw_of(R):
    # yaw about world z (atan2 of R[1,0], R[0,0])
    return np.arctan2(R[1, 0], R[0, 0])


def main():
    path = sys.argv[1]
    name = sys.argv[2] if len(sys.argv) > 2 else path
    t, p, q = load_gt(path)
    N = len(t)
    dt = np.diff(t)
    dp = np.diff(p, axis=0)
    step_trans = np.linalg.norm(dp, axis=1)  # m per GT step

    # per-step rotation angle + yaw delta
    Rs = [quat_to_R(q[i]) for i in range(N)]
    step_rot = np.array([rot_angle(Rs[i], Rs[i + 1]) for i in range(N - 1)])
    yaws = np.array([yaw_of(Rs[i]) for i in range(N)])
    dyaw = np.diff(yaws)
    dyaw = (dyaw + np.pi) % (2 * np.pi) - np.pi  # wrap

    # speeds / rates (per second), guard tiny dt
    good = dt > 1e-6
    speed = step_trans[good] / dt[good]
    yaw_rate = np.abs(dyaw[good]) / dt[good]
    rot_rate = step_rot[good] / dt[good]

    total_time = t[-1] - t[0]
    total_dist = float(np.sum(step_trans))
    total_yaw = float(np.sum(np.abs(dyaw)))  # total |yaw turned| (rad)
    total_rot = float(np.sum(step_rot))

    # fraction of TIME turning (yaw rate above a small threshold)
    turn_thresh = np.deg2rad(2.0)  # 2 deg/s
    frac_turn = float(np.mean(yaw_rate > turn_thresh))

    mean_dt = float(np.mean(dt[good]))

    print("=== %s ===" % name)
    print("  samples           : %d" % N)
    print("  total_time_s       : %.1f" % total_time)
    print("  mean_GT_dt_s       : %.4f  (%.1f Hz)" % (mean_dt, 1.0 / mean_dt))
    print("  total_dist_m       : %.1f" % total_dist)
    print("  total_yaw_turned_rad: %.2f  (%.0f deg)" % (total_yaw, np.rad2deg(total_yaw)))
    print("  total_rot_angle_rad : %.2f" % total_rot)
    print("  frac_time_turning  : %.3f  (>2 deg/s)" % frac_turn)
    print("  --- speed (m/s) ---")
    print("    mean=%.3f median=%.3f p90=%.3f max=%.3f" % (
        np.mean(speed), np.median(speed), np.percentile(speed, 90), np.max(speed)))
    print("  --- yaw rate (rad/s) ---")
    print("    mean=%.4f median=%.4f p90=%.4f max=%.4f" % (
        np.mean(yaw_rate), np.median(yaw_rate), np.percentile(yaw_rate, 90), np.max(yaw_rate)))
    print("  --- rot rate (rad/s) ---")
    print("    mean=%.4f median=%.4f p90=%.4f max=%.4f" % (
        np.mean(rot_rate), np.median(rot_rate), np.percentile(rot_rate, 90), np.max(rot_rate)))
    print("  --- key squared stats (the motion-Q hypothesis levers) ---")
    print("    mean_speed^2        = %.4f  (m/s)^2" % (np.mean(speed) ** 2))
    print("    mean(speed^2)       = %.4f" % np.mean(speed ** 2))
    print("    rms_speed           = %.4f" % np.sqrt(np.mean(speed ** 2)))
    print("    mean_yaw_rate^2     = %.6f  (rad/s)^2" % (np.mean(yaw_rate) ** 2))
    print("    mean(yaw_rate^2)    = %.6f" % np.mean(yaw_rate ** 2))
    print("    rms_yaw_rate        = %.6f" % np.sqrt(np.mean(yaw_rate ** 2)))


if __name__ == "__main__":
    main()
