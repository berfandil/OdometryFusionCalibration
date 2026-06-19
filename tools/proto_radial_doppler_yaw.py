#!/usr/bin/env python3
# proto_radial_doppler_yaw.py - Doppler-NATIVE ego-motion from a single RADIal radar frame, and the
# falsifiable test of whether YAW RATE is recoverable from Doppler (the direction flagged at the end of
# RADAR_SCAN_ODOMETRY.md §5 after the descriptor scan-match failed to recover rotation on real radar).
#
# OBSERVABILITY (why this is the right question):
#   A static world point at sensor-frame position r (bearing unit u=r/|r|) shows radial velocity
#       v_radial = -(v + omega x r) . u = -v.u           because (omega x r).u = 0  (omega x r _|_ r).
#   So a single radar's INSTANTANEOUS Doppler gives ONLY the linear velocity v of the SENSOR -- omega
#   does NOT appear per-point. omega couples in ONLY through the lever arm: the radar's own linear
#   velocity is v_radar = v_cog + omega x L. For a forward lever Lx under (near-)non-holonomic motion
#   (v_cog ~ (speed, 0) in body frame), the LATERAL component is  v_y = omega_z * Lx. Hence
#       omega_z = v_y / Lx .
#   This is the SAME lever<->rotation coupling the calibrator uses (lever from known omega); here we run
#   it FORWARD (known/assumed lever -> omega). The signal is small: an ~0.08 rad/s urban turn with
#   Lx~3.5 m gives v_y ~ 0.28 m/s, near the Doppler-fit noise floor -> the test is whether the per-frame
#   lateral Doppler velocity v_y CORRELATES with the lidar yaw rate at all.
#
# METHOD per frame:
#   RANSAC fit of  v_radial = a*cos(az) + b*sin(az)  over the detections (a=-vx, b=-vy), rejecting
#   movers as the residual minority -> per-frame (vx, vy) = sensor velocity in body frame. Then the
#   forward speed = vx (compare to lidar |disp|/dt -> validates the LINEAR recovery + calibrates hz),
#   and the lateral velocity vy -> omega_z = vy/Lx (the rotation test). HEADLINE = corr(vy, lidar
#   per-frame yaw); if that is ~0, single-radar Doppler yaw is buried (the honest negative); if high,
#   yaw is recoverable up to the lever scale Lx (fit Lx by regressing vy on lidar yaw-rate).
#
# Reference + GT = lidar ICP odometry (proto_radial_lidar_odom). Pure ASCII; numpy; seeded.
#
# Usage: python tools/proto_radial_doppler_yaw.py <run_dir> [--hz HZ] [--lx M] [--rmin M] [--rmax M]
#        [--ransac-iters N] [--tau M_S] [--quiet]
import math
import os
import sys
import io
import contextlib

import numpy as np

import proto_radial_lidar_odom as plo


def ransac_ego(arr, rmin, rmax, az_max, iters, tau, rng):
    # RANSAC fit v_radial = a cos(az) + b sin(az) over FOV detections; reject movers. Return
    # ((a,b), n_inl, n_pts) with a=-vx, b=-vy (sensor body-frame velocity). 2-point minimal sample.
    x, y, v = arr[5], arr[6], arr[8]
    r = np.sqrt(x * x + y * y + np.asarray(arr[7]) ** 2)
    az = np.arctan2(y, x)
    keep = (r >= rmin) & (r <= rmax) & (np.abs(np.degrees(az)) <= az_max) & (x > 0)
    az = az[keep]; vv = np.asarray(v)[keep]
    n = az.shape[0]
    if n < 6:
        return None, 0, n
    A = np.stack([np.cos(az), np.sin(az)], 1)
    best_inl = None
    for _ in range(iters):
        s = rng.sample(range(n), 2)
        try:
            sol = np.linalg.solve(A[s], vv[s])
        except np.linalg.LinAlgError:
            continue
        res = np.abs(vv - A @ sol)
        inl = np.where(res < tau)[0]
        if best_inl is None or len(inl) > len(best_inl):
            best_inl = inl
    if best_inl is None or len(best_inl) < 4:
        return None, 0, n
    sol, *_ = np.linalg.lstsq(A[best_inl], vv[best_inl], rcond=None)
    return (float(sol[0]), float(sol[1])), len(best_inl), n


def run(run_dir, args):
    import random
    rng = random.Random(args.seed)
    radar_dir = os.path.join(run_dir, 'radar')
    files = sorted(f for f in os.listdir(radar_dir) if f.endswith('.npy'))

    # lidar GT: per-frame yaw + displacement
    la = plo.parse_args([run_dir, '--voxel', '1.0', '--target-n', '1200'])
    with contextlib.redirect_stdout(io.StringIO()):
        poses_R, poses_p, incs = plo.run(run_dir, run_dir, la)
    gt_yaw = {}; gt_disp = {}
    for (k, R_inc, t_inc, n_inl, rms, npts) in incs:
        gt_yaw[k] = math.degrees(math.atan2(R_inc[1][0], R_inc[0][0]))
        gt_disp[k] = math.sqrt(t_inc[0] ** 2 + t_inc[1] ** 2 + t_inc[2] ** 2)

    dt = 1.0 / args.hz
    vxs, vys, gyaws, gspeeds, fk = [], [], [], [], []
    for k, fn in enumerate(files):
        arr = np.load(os.path.join(radar_dir, fn), allow_pickle=True)
        sol, n_inl, n = ransac_ego(arr, args.rmin, args.rmax, args.az_max,
                                   args.ransac_iters, args.tau, rng)
        if sol is None or k not in gt_yaw:
            continue
        a, b = sol
        vx, vy = -a, -b
        vxs.append(vx); vys.append(vy)
        gyaws.append(gt_yaw[k]); gspeeds.append(gt_disp[k] / dt); fk.append(k)

    vxs = np.array(vxs); vys = np.array(vys)
    gyaws = np.array(gyaws); gspeeds = np.array(gspeeds)

    # LINEAR check: radar forward speed |vx| vs lidar speed (validates Doppler + calibrates hz)
    speed_corr = float(np.corrcoef(np.abs(vxs), gspeeds)[0, 1]) if len(vxs) > 2 else float('nan')
    # implied hz so that |vx| matches lidar disp/frame: dt_imp = median(disp/|vx|)
    mask = (np.abs(vxs) > 0.5)
    dt_imp = float(np.median(np.array([gt_disp[k] for k in fk])[mask] / np.abs(vxs)[mask])) if mask.any() else float('nan')

    # ROTATION test: does lateral velocity vy correlate with lidar per-frame yaw?
    yaw_corr = float(np.corrcoef(vys, gyaws)[0, 1]) if len(vys) > 2 else float('nan')
    # fit lever Lx from vy = omega_z * Lx ; lidar omega_z (rad/s) = radians(gt_yaw)/dt ;
    # so vy = (radians(gt_yaw)/dt) * Lx -> slope of vy vs (radians(gt_yaw)/dt) = Lx
    omega = np.radians(gyaws) / dt
    if np.var(omega) > 1e-9:
        Lx_fit = float(np.dot(vys, omega) / np.dot(omega, omega))
    else:
        Lx_fit = float('nan')

    print('=' * 96)
    print('DOPPLER-NATIVE YAW TEST  %s   hz=%.2f Lx=%.2f' %
          (os.path.basename(run_dir.rstrip('/\\')), args.hz, args.lx))
    print('  frames fit: %d / %d' % (len(vxs), len(files)))
    print('  LINEAR (sanity): forward speed |vx| vs lidar speed corr = %.3f   implied hz = %.2f'
          % (speed_corr, 1.0 / dt_imp if dt_imp == dt_imp and dt_imp > 0 else float('nan')))
    print('     radar |vx| med %.2f m/s   lidar speed med %.2f m/s'
          % (float(np.median(np.abs(vxs))), float(np.median(gspeeds))))
    print('  ROTATION (headline): lateral vy vs lidar per-frame yaw   corr = %.3f' % yaw_corr)
    print('     vy std %.3f m/s   (expected lever signal omega*Lx ~ %.3f m/s at peak)'
          % (float(np.std(vys)), float(np.max(np.abs(omega))) * args.lx))
    print('     lever Lx fit (vy = omega_z*Lx) = %.2f m   (physical front lever ~ 2-4 m if real)' % Lx_fit)
    # integrated yaw via omega_z = vy/Lx
    yaw_from_vy = float(np.degrees(np.sum(vys / args.lx * dt)))
    gt_net = float(np.sum(gyaws))
    print('  integrated yaw (omega=vy/Lx): recovered %.1f deg   lidar GT %.1f deg' % (yaw_from_vy, gt_net))
    print('=' * 96)


class Args:
    pass


def parse_args(argv):
    a = Args()
    a.hz = 4.1
    a.lx = 3.5
    a.rmin = 1.0
    a.rmax = 100.0
    a.az_max = 70.0
    a.ransac_iters = 200
    a.tau = 0.5
    a.seed = 12345
    a.quiet = False
    pos = []
    i = 0
    while i < len(argv):
        x = argv[i]
        if x == '--hz': i += 1; a.hz = float(argv[i])
        elif x == '--lx': i += 1; a.lx = float(argv[i])
        elif x == '--rmin': i += 1; a.rmin = float(argv[i])
        elif x == '--rmax': i += 1; a.rmax = float(argv[i])
        elif x == '--az-max': i += 1; a.az_max = float(argv[i])
        elif x == '--ransac-iters': i += 1; a.ransac_iters = int(argv[i])
        elif x == '--tau': i += 1; a.tau = float(argv[i])
        elif x == '--seed': i += 1; a.seed = int(argv[i])
        elif x == '--quiet': a.quiet = True
        else: pos.append(x)
        i += 1
    a.pos = pos
    return a


def main():
    args = parse_args(sys.argv[1:])
    if len(args.pos) < 1:
        print('usage: proto_radial_doppler_yaw.py <run_dir> [--hz HZ] [--lx M] [--tau M_S] [--quiet]')
        sys.exit(2)
    run(args.pos[0], args)


if __name__ == '__main__':
    main()
