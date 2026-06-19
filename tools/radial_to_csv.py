#!/usr/bin/env python3
# radial_to_csv.py - run the PROVEN descriptor scan-matching odometry (RADAR_SCAN_ODOMETRY.md, the 3D
# path) on REAL 4D imaging-radar point clouds from the RADIal dataset, with LIDAR scan-match odometry as
# the reference + de-facto GT. THE GOAL (HANDOFF "REAL 4D RADAR - RADIal"): confirm the synth-4D rotation
# unlock (synth_4d_radar.py: yaw correlation 0 -> 0.995 as detection density rises) on REAL imaging radar.
#
# WHY LIDAR IS THE REFERENCE: RADIal Ready_to_use ships NO GPS / CAN / IMU / ego-pose (only camera, the
# Valeo Scala lidar, the TI-cascade radar, and freespace labels). So the dense lidar (laser_PCL, ~14k
# pts/frame) is the only ego-motion truth: tools/proto_radial_lidar_odom.py recovers the per-frame ego
# SE(3) by ICP. That ego delta is (a) the static-filter PRIOR for the radar (replaces nuScenes CAN) and
# (b) the GT to score the radar-recovered motion against. The radar SCAN-MATCH itself is INDEPENDENT of
# this prior (it recovers motion from radar geometry alone); the prior only rejects movers + gates
# RANSAC blunders, exactly as the CAN prior did in radar_scan_odometry.py.
#
# REAL 4D RADAR FORMAT (radar_PCL/<frame>.npy, shape 9 x Npts; ~500 detections/frame):
#   row0 range(m)  row1 azimuth(deg)  row2 elevation(deg, 1-deg quantized)  row3 power
#   row4 doppler_bin  row5 x  row6 y  row7 z  row8 v (radial velocity, m/s, signed)
# We hand the 3D static filter (vx,vy,vz) = v * unit_bearing (a real radar measures only the radial
# component; static_filter_3d projects onto the bearing anyway -- same convention as synth_4d_radar).
#
# Frames carry NO real timestamp; both sensors are stamped at a nominal fixed --hz (RADIal radar ~ low
# rate). Pose recovery (Kabsch) is dt-INDEPENDENT; --hz only scales the static-filter twist + the emitted
# CSV timestamps. Radar frame k is time-aligned with lidar frame k (RADIal synced them; same numSample).
#
# Usage:
#   python tools/radial_to_csv.py <run_dir> [<out_dir>] [--hz HZ] [--doppler-sign S] [--tau-v M_S]
#       [--num-points N] [--baseline N] [--mount yaw,pitch,roll,lx,ly,lz] [--emit] [--quiet]
#   run_dir = radial_extract/<run>/ (has radar/NNNN.npy + laser/NNNN.npy). Default out_dir = run_dir.
import math
import os
import sys

import numpy as np

import nuscenes_to_csv as n2c
import radar_scan_odometry as rso
import synth_4d_radar as s4d
import proto_radial_lidar_odom as plo


def read_radar_dets(path, rmin, rmax, az_max, doppler_sign, num_points=0):
    # radar_PCL .npy (9 x Npts). Return list of (x,y,z, vx,vy,vz) detections; velocity vector is the
    # signed radial speed * unit bearing. FOV/range gate (forward x>0, |az|<=az_max, range in [rmin,rmax]).
    # num_points>0 caps to the N NEAREST returns (closest scatterers dominate + the nearest-N set is
    # near-stable frame-to-frame -> descriptor persistence; mirrors synth_4d_radar) and bounds the
    # O(N^2) descriptor/match cost.
    arr = np.load(path, allow_pickle=True)
    x, y, z, v = arr[5], arr[6], arr[7], arr[8]
    rng = np.sqrt(x * x + y * y + z * z)
    az = np.degrees(np.arctan2(y, x))
    keep = (rng >= rmin) & (rng <= rmax) & (np.abs(az) <= az_max) & (x > 0.0)
    idx = np.where(keep)[0]
    if num_points > 0 and idx.shape[0] > num_points:
        near = np.argpartition(rng[idx], num_points - 1)[:num_points]
        idx = idx[near]
    dets = []
    for i in idx:
        r = float(rng[i])
        if r < 1e-6:
            continue
        ux, uy, uz = float(x[i]) / r, float(y[i]) / r, float(z[i]) / r
        vr = doppler_sign * float(v[i])
        dets.append((float(x[i]), float(y[i]), float(z[i]), vr * ux, vr * uy, vr * uz))
    return dets


def lidar_reference(run_dir, args):
    # Run the lidar ICP odometry -> per-frame ego deltas (the static-filter prior + GT) + abs poses.
    # Returns (deltas, poses_R, poses_p) where deltas[k] (k>=1) = (R_inc 3x3 list, t_inc 3 list) is the
    # ego body delta over [k-1, k]. deltas[0] = identity (no predecessor).
    la = plo.parse_args([run_dir, '--voxel', str(args.lidar_voxel),
                         '--target-n', str(args.lidar_n)])
    poses_R, poses_p, incs = plo.run(run_dir, run_dir, la)
    n = len(poses_R)
    deltas = [([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]], [0.0, 0.0, 0.0])] * n
    for (k, R_inc, t_inc, n_inl, rms, npts) in incs:
        deltas[k] = (R_inc, t_inc)
    return deltas, poses_R, poses_p


def run(run_dir, out_dir, args):
    import random
    radar_dir = os.path.join(run_dir, 'radar')
    files = sorted(f for f in os.listdir(radar_dir) if f.endswith('.npy'))
    nframe = len(files)

    # --- reference + GT: lidar ICP odometry (ego deltas, same frame indexing as radar) ---
    deltas, poses_R, poses_p = lidar_reference(run_dir, args)

    # radar mount prior X_radar = ego(lidar)_from_radar; default identity (unknown). Yaw is preserved
    # under any lever; a small mount rotation barely shifts the radar yaw from the ego yaw -> the lidar
    # ego yaw is a valid GT for the radar yaw-correlation headline.
    R_x, t_x = s4d.mount_from_spec(args.mount)

    period_us = int(round(1e6 / args.hz))
    dt = 1.0 / args.hz
    N = max(1, int(args.baseline))
    ransac_rng = random.Random(args.seed * 1000003 + 777)

    stats = {"n_pairs": 0, "n_fallback": 0, "sum_det": 0, "sum_static": 0,
             "sum_corr": 0, "sum_inl": 0, "n_recov": 0, "n_gated": 0}
    perframe = []
    rows = []
    hist = []   # (k, dets)
    for k, fn in enumerate(files):
        dets = read_radar_dets(os.path.join(radar_dir, fn), args.rmin, args.rmax,
                               args.az_max, args.doppler_sign, args.num_points)
        hist.append((k, dets))
        if len(hist) > N + 1:
            hist.pop(0)
        if k < N:
            continue
        k0, prev_dets = hist[0]

        # --- ego delta over [k-N, k] from lidar (compose the per-frame deltas) ---
        A_R = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
        A_t = [0.0, 0.0, 0.0]
        for j in range(k0 + 1, k + 1):
            Rd, td = deltas[j]
            A_t = [A_t[i] + n2c.mat_vec(A_R, td)[i] for i in range(3)]
            A_R = n2c.mat_mul(A_R, Rd)
        # expected radar motion B = X_radar^-1 o A o X_radar (the GT radar motion + static-filter prior)
        B_R, B_t = rso.expected_radar_motion(A_R, A_t, R_x, t_x)
        dt_win = dt * N

        # --- A. 3D static filter both frames with the ego prior ---
        prev_kept = rso.static_filter_3d(prev_dets, B_R, B_t, dt_win, args.tau_v)
        curr_kept = rso.static_filter_3d(dets, B_R, B_t, dt_win, args.tau_v)
        stats["sum_det"] += len(dets)
        stats["sum_static"] += len(curr_kept)
        stats["n_pairs"] += 1

        # --- B. descriptors (3D) ---
        descs_prev = rso.build_descriptors(prev_kept)
        descs_curr = rso.build_descriptors(curr_kept)

        # --- C. match (fast exact matcher) ---
        pairs = s4d.match_frames_fast(descs_curr, descs_prev, args.tau_d, args.rho)
        stats["sum_corr"] += len(pairs)

        # --- D. recover (3D Kabsch + RANSAC, minimal sample 3) ---
        R3, t3, n_inl, n_corr = rso.ransac_kabsch(curr_kept, prev_kept, pairs, ransac_rng,
                                                  args.tau_inlier, args.ransac_iters,
                                                  args.min_corr, dim=3)
        stats["sum_inl"] += n_inl

        recovered = R3 is not None
        if recovered:
            (vp, omp) = rso.twist_from_se3_3d(B_R, B_t, dt_win)
            rot_pred = math.sqrt(omp[0] ** 2 + omp[1] ** 2 + omp[2] ** 2) * dt_win
            spd_pred = math.sqrt(vp[0] ** 2 + vp[1] ** 2 + vp[2] ** 2) * dt_win
            rot_inc = math.radians(rso.rot_err_deg(R3.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
            t_mag = float(np.linalg.norm(t3))
            if rot_inc > rot_pred + math.radians(args.yaw_margin_deg) or t_mag > spd_pred + args.trans_margin_m:
                stats["n_gated"] += 1
                recovered = False

        if not recovered:
            R_inc = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            t_inc = [0.0, 0.0, 0.0]
            var6 = [100.0, 100.0, 100.0, 0.25, 0.25, 0.25]
            stats["n_fallback"] += 1
        else:
            R_inc = [[float(R3[a, b]) for b in range(3)] for a in range(3)]
            t_inc = [float(t3[0]), float(t3[1]), float(t3[2])]
            stats["n_recov"] += 1
            s_t = (0.10) ** 2
            s_r = max(0.01, 0.3 / max(1, n_inl)) ** 2
            var6 = [s_t, s_t, s_t, s_r, s_r, s_r]
        rows.append((k * period_us, R_inc, t_inc, var6))

        # GT radar motion (= B from the lidar ego delta) for the vs-GT validation
        perframe.append((R_inc, t_inc, B_R, B_t, recovered, k * period_us, k0 * period_us, k))

    rep = rso.report_radar(os.path.basename(run_dir.rstrip('/\\')), perframe, stats, args.baseline)

    if args.emit:
        os.makedirs(out_dir, exist_ok=True)
        stem = os.path.basename(run_dir.rstrip('/\\'))
        nsuf = "" if N == 1 else ("_n%d" % N)
        outp = os.path.join(out_dir, stem + '_radar_scan' + nsuf + '.csv')
        f = n2c.open_inc(outp, 'RADIal 4D radar scan-match odometry (3D R+t), baseline N=%d' % N)
        for (t_us, R, t, var6) in rows:
            n2c.emit_inc(f, t_us * n2c.US_TO_NS, R, t, var6)
        f.close()
        rep['_csv'] = outp
    return rep, stats, perframe


def print_report(run_dir, rep, stats, args):
    print('=' * 104)
    print('REAL 4D RADAR SCAN-MATCH (RADIal)  %s   N=%d hz=%.1f dop_sign=%+d tau_v=%.2f'
          % (os.path.basename(run_dir.rstrip('/\\')), max(1, args.baseline), args.hz,
             args.doppler_sign, args.tau_v))
    print('=' * 104)
    print('  pairs %d   avg det %.0f  static %.0f  corr %.1f  inlier %.1f   fallbk %.0f%% gated %.0f%%'
          % (rep['n_pairs'], rep['avg_det'], rep['avg_static'], rep['avg_corr'], rep['avg_inl'],
             rep['fallback_pct'], rep['gated_pct']))
    print('  RECOVERED vs LIDAR-GT (over recovered frames):')
    print('    trans  med %.3f m  p90 %.3f m     rot med %.3f deg  p90 %.3f deg'
          % (rep['trans_med'], rep['trans_p90'], rep['rot_med'], rep['rot_p90']))
    print('  ROTATION HEADLINE:')
    print('    per-step YAW CORRELATION (recovered vs GT)  = %.3f   (synth-4D: 0 -> 0.995 with density)'
          % rep['yaw_corr'])
    print('    cumulative heading  recovered %+.2f deg   GT %+.2f deg   err %+.2f deg'
          % (rep['cum_est_yaw_deg'], rep['cum_gt_yaw_deg'],
             rep['cum_est_yaw_deg'] - rep['cum_gt_yaw_deg']))
    print('    non-overlap heading recovered %+.2f deg   GT %+.2f deg   err %+.2f deg'
          % (rep['nov_est_yaw_deg'], rep['nov_gt_yaw_deg'], rep['nov_heading_err_deg']))
    if rep.get('_csv'):
        print('  wrote %s' % rep['_csv'])
    print('=' * 104)


class Args:
    pass


def parse_args(argv):
    a = Args()
    a.hz = 17.0
    a.doppler_sign = 1.0
    a.tau_v = 1.5
    a.rmin = 1.0
    a.rmax = 100.0
    a.az_max = 70.0
    a.num_points = 300        # cap to N nearest FOV detections (0 = all; bounds O(N^2) + aids persistence)
    a.baseline = 1
    a.mount = "0,0,0,0,0,0"    # X_radar = lidar_from_radar prior; identity (unknown) by default
    a.tau_d = 0.5
    a.rho = 0.5
    a.tau_inlier = 0.5
    a.ransac_iters = 200
    a.min_corr = 6
    a.yaw_margin_deg = 8.0
    a.trans_margin_m = 1.5
    a.lidar_voxel = 1.0
    a.lidar_n = 1200
    a.seed = 12345
    a.emit = False
    a.quiet = False
    pos = []
    i = 0
    while i < len(argv):
        x = argv[i]
        if x == '--hz': i += 1; a.hz = float(argv[i])
        elif x == '--doppler-sign': i += 1; a.doppler_sign = float(argv[i])
        elif x == '--tau-v': i += 1; a.tau_v = float(argv[i])
        elif x == '--rmin': i += 1; a.rmin = float(argv[i])
        elif x == '--rmax': i += 1; a.rmax = float(argv[i])
        elif x == '--az-max': i += 1; a.az_max = float(argv[i])
        elif x == '--num-points': i += 1; a.num_points = int(argv[i])
        elif x == '--baseline': i += 1; a.baseline = int(argv[i])
        elif x == '--mount': i += 1; a.mount = argv[i]
        elif x == '--tau-d': i += 1; a.tau_d = float(argv[i])
        elif x == '--rho': i += 1; a.rho = float(argv[i])
        elif x == '--tau-inlier': i += 1; a.tau_inlier = float(argv[i])
        elif x == '--ransac-iters': i += 1; a.ransac_iters = int(argv[i])
        elif x == '--min-corr': i += 1; a.min_corr = int(argv[i])
        elif x == '--yaw-margin-deg': i += 1; a.yaw_margin_deg = float(argv[i])
        elif x == '--trans-margin-m': i += 1; a.trans_margin_m = float(argv[i])
        elif x == '--tau-inlier': i += 1; a.tau_inlier = float(argv[i])
        elif x == '--lidar-voxel': i += 1; a.lidar_voxel = float(argv[i])
        elif x == '--lidar-n': i += 1; a.lidar_n = int(argv[i])
        elif x == '--seed': i += 1; a.seed = int(argv[i])
        elif x == '--emit': a.emit = True
        elif x == '--quiet': a.quiet = True
        else: pos.append(x)
        i += 1
    a.pos = pos
    return a


def main():
    args = parse_args(sys.argv[1:])
    if len(args.pos) < 1:
        print('usage: radial_to_csv.py <run_dir> [<out_dir>] [--hz HZ] [--doppler-sign S] [--tau-v M_S] '
              '[--baseline N] [--mount y,p,r,lx,ly,lz] [--emit] [--quiet]')
        sys.exit(2)
    run_dir = args.pos[0]
    out_dir = args.pos[1] if len(args.pos) > 1 else run_dir
    rep, stats, perframe = run(run_dir, out_dir, args)
    print_report(run_dir, rep, stats, args)


if __name__ == '__main__':
    main()
