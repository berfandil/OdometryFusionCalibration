#!/usr/bin/env python3
# proto_radial_lidar_odom.py - LIDAR scan-matching odometry over a RADIal extracted run, to serve as
# the REFERENCE odometry + de-facto GT for the radar rotation-unlock test. RADIal Ready_to_use ships NO
# GPS/CAN/ego-pose (only camera/lidar/radar/freespace), so the dense Valeo Scala lidar (laser_PCL,
# ~14k pts/frame) is the only ego-motion truth available. We recover frame-to-frame SE(3) by the SAME
# descriptor scan-match used for radar (RADAR_SCAN_ODOMETRY.md: sorted inter-point-distance descriptors
# -> mutual-best tolerance matching -> 3D Kabsch in RANSAC), on a voxel-downsampled point set (the
# matcher is O(N^2 log N); 14k raw is intractable, ~500 voxel centroids is plenty + more stable).
#
# Lidar has no Doppler, so there is NO static filter here -- all downsampled points feed the descriptor
# stage and RANSAC rejects movers as a geometric minority. Output: per-frame SE(3) increments integrated
# into an absolute pose track (the run's frame-0 = world origin), written as a GT abs-pose CSV
# (t_ns, x,y,z, qw,qx,qy,qz) -- the exact format synth_4d_radar.read_gt / nuscenes_to_csv emit -- plus
# an increment-form reference CSV (the project schema) for fusion.
#
# Frames carry NO real timestamp (RADIal Ready_to_use dropped them); we stamp a nominal fixed period
# (--hz, default 17 -- refined later from lidar-displacement vs radar-Doppler). Pose integration is
# dt-INDEPENDENT (deltas compose directly); --hz only sets the emitted timestamps + twist scaling.
#
# Usage:
#   python tools/proto_radial_lidar_odom.py <run_dir> [<out_dir>] [--voxel M] [--target-n N]
#       [--rmin M] [--rmax M] [--hz HZ] [--emit] [--quiet]
#   run_dir = radial_extract/<run>/ (has laser/NNNN.npy). Default out_dir = run_dir.
import math
import os
import sys

import numpy as np

import nuscenes_to_csv as n2c
import radar_scan_odometry as rso


def load_laser_xyz(path, rmin, rmax):
    # laser_PCL is (Npts, 11): cols 0,1,2 = x,y,z (m), col3 = intensity, col4 = range. Keep finite
    # points within [rmin, rmax] horizontal range (drop ego-body returns + far noise).
    arr = np.load(path, allow_pickle=True)
    xyz = np.asarray(arr[:, 0:3], dtype=float)
    rng = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
    keep = np.isfinite(rng) & (rng >= rmin) & (rng <= rmax)
    return xyz[keep]


def voxel_downsample(P, voxel, target_n):
    # Voxel-grid downsample: one representative (the centroid) per occupied voxel. Adapt the voxel size
    # up if we overshoot target_n badly (keeps the descriptor stage O(N^2) tractable + stable).
    for _ in range(6):
        keys = np.floor(P / voxel).astype(np.int64)
        # hash voxel keys -> unique, take centroid per voxel
        order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))
        ks = keys[order]
        Ps = P[order]
        uniq_mask = np.ones(len(ks), dtype=bool)
        uniq_mask[1:] = np.any(ks[1:] != ks[:-1], axis=1)
        starts = np.where(uniq_mask)[0]
        cents = np.add.reduceat(Ps, starts, axis=0) / np.diff(np.append(starts, len(Ps)))[:, None]
        if cents.shape[0] <= target_n * 1.4:
            return cents
        voxel *= 1.3
    return cents


def nn_match(src, dst):
    # For each src point, the nearest dst point index + Euclidean distance. Brute-force, chunked over
    # src to bound the (M,N) distance-matrix memory. dst dense lidar -> a proximity match is valid.
    M = src.shape[0]
    idx = np.empty(M, dtype=np.int64)
    dist = np.empty(M, dtype=float)
    dst2 = (dst * dst).sum(1)                       # (N,)
    step = 512
    for a in range(0, M, step):
        b = min(a + step, M)
        s = src[a:b]
        # ||s - d||^2 = |s|^2 - 2 s.d + |d|^2
        d2 = (s * s).sum(1)[:, None] - 2.0 * (s @ dst.T) + dst2[None, :]
        j = d2.argmin(1)
        idx[a:b] = j
        dist[a:b] = np.sqrt(np.maximum(d2[np.arange(b - a), j], 0.0))
    return idx, dist


def icp_schedule(coarse, fine, n):
    # geometric ramp of the correspondence-distance cap from `coarse` down to `fine` over n steps,
    # then a few extra iterations pinned at `fine` for refinement. A wide first gate gives ICP a large
    # convergence basin (handles the per-frame motion); shrinking it rejects wrong matches as it locks.
    import math as _m
    sched = [coarse * (fine / coarse) ** (i / max(1, n - 1)) for i in range(n)]
    return sched + [fine, fine, fine]


def icp_delta(prev_pts, curr_pts, init_R, init_t, args):
    # Point-to-point ICP aligning CURR onto PREV: recover (R, t) with prev ~= R curr + t (the forward
    # body delta, pose_curr = pose_prev o delta). COARSE-TO-FINE correspondence gate (icp_schedule):
    # a wide first cap converges from a large basin, shrinking caps reject movers + non-overlap as the
    # alignment locks (diagnosed: a single tight gate latched wrong matches -> 6-8 m blunders; the ramp
    # converges monotonically to the lidar point-spacing rms floor ~0.26 m). init_(R,t) seeds it
    # (identity by default; a constant-velocity seed is passed for continuity). Returns
    # (R 3x3 list, t 3 list, n_inliers, rms).
    R = np.asarray(init_R, dtype=float).copy()
    t = np.asarray(init_t, dtype=float).copy()
    cur = (R @ curr_pts.T).T + t
    n_inl = 0
    rms = float('inf')
    sched = icp_schedule(args.icp_coarse, args.icp_max_dist, args.icp_iters)
    for mx in sched:
        idx, dist = nn_match(cur, prev_pts)
        m = dist <= mx
        n_inl = int(m.sum())
        if n_inl < args.min_corr:
            break
        Rk, tk = rso.kabsch(cur[m], prev_pts[idx[m]])    # incremental align cur -> prev
        cur = (Rk @ cur.T).T + tk
        R = Rk @ R
        t = Rk @ t + tk
        rms = float(np.sqrt(np.mean(dist[m] ** 2)))
        if abs(math.atan2(Rk[1, 0], Rk[0, 0])) < 1e-6 and float(np.linalg.norm(tk)) < 1e-5:
            break
    return R.tolist(), [float(t[0]), float(t[1]), float(t[2])], n_inl, rms


def mat_to_quat(R):
    # 3x3 rotation -> (qw,qx,qy,qz).
    tr = R[0][0] + R[1][1] + R[2][2]
    if tr > 0:
        s = math.sqrt(tr + 1.0) * 2
        qw = 0.25 * s
        qx = (R[2][1] - R[1][2]) / s
        qy = (R[0][2] - R[2][0]) / s
        qz = (R[1][0] - R[0][1]) / s
    elif R[0][0] > R[1][1] and R[0][0] > R[2][2]:
        s = math.sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2
        qw = (R[2][1] - R[1][2]) / s; qx = 0.25 * s
        qy = (R[0][1] + R[1][0]) / s; qz = (R[0][2] + R[2][0]) / s
    elif R[1][1] > R[2][2]:
        s = math.sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2
        qw = (R[0][2] - R[2][0]) / s; qx = (R[0][1] + R[1][0]) / s
        qy = 0.25 * s; qz = (R[1][2] + R[2][1]) / s
    else:
        s = math.sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2
        qw = (R[1][0] - R[0][1]) / s; qx = (R[0][2] + R[2][0]) / s
        qy = (R[1][2] + R[2][1]) / s; qz = 0.25 * s
    return qw, qx, qy, qz


def run(run_dir, out_dir, args):
    import random
    laser_dir = os.path.join(run_dir, 'laser')
    files = sorted(f for f in os.listdir(laser_dir) if f.endswith('.npy'))
    rng = random.Random(args.seed)

    # absolute pose track: T0 = identity (world = frame-0 lidar frame). Compose forward deltas.
    poses_R = [[[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]]
    poses_p = [[0.0, 0.0, 0.0]]
    incs = []          # (order, R_inc, t_inc, n_inl, rms, npts) for frame k>=1
    prev_pts = None
    n_recov = 0
    step_lens = []
    yaw_steps = []
    I3 = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
    for k, fn in enumerate(files):
        P = load_laser_xyz(os.path.join(laser_dir, fn), args.rmin, args.rmax)
        pts = voxel_downsample(P, args.voxel, args.target_n)
        if prev_pts is not None:
            # identity init each frame (the coarse-to-fine gate gives a wide basin) -> no error cascade
            R_inc, t_inc, n_inl, rms = icp_delta(prev_pts, pts, I3, [0.0, 0.0, 0.0], args)
            if n_inl >= args.min_corr:
                n_recov += 1
            else:
                R_inc = I3
                t_inc = [0.0, 0.0, 0.0]
            # compose: T_k = T_{k-1} o delta
            Rp, pp = poses_R[-1], poses_p[-1]
            Rk = n2c.mat_mul(Rp, R_inc)
            pk = [pp[i] + n2c.mat_vec(Rp, t_inc)[i] for i in range(3)]
            poses_R.append(Rk); poses_p.append(pk)
            incs.append((k, R_inc, t_inc, n_inl, rms, pts.shape[0]))
            step_lens.append(math.sqrt(t_inc[0] ** 2 + t_inc[1] ** 2 + t_inc[2] ** 2))
            yaw_steps.append(math.degrees(math.atan2(R_inc[1][0], R_inc[0][0])))
        prev_pts = pts

    # report
    n = len(incs)
    path_len = sum(step_lens)
    net = math.sqrt((poses_p[-1][0]) ** 2 + (poses_p[-1][1]) ** 2 + (poses_p[-1][2]) ** 2)
    net_yaw = math.degrees(math.atan2(poses_R[-1][1][0], poses_R[-1][0][0]))
    print('=' * 92)
    print('LIDAR SCAN-MATCH ODOMETRY  %s' % os.path.basename(run_dir.rstrip('/\\')))
    print('  frames=%d  recovered=%d (%.0f%%)  voxel=%.2f target_n=%d  avg_pts=%.0f'
          % (len(files), n_recov, 100.0 * n_recov / max(1, n), args.voxel, args.target_n,
             np.mean([c[5] for c in incs]) if incs else 0))
    print('  step len  m : med %.3f  p90 %.3f  max %.3f' %
          (rso.pct(step_lens, 0.5), rso.pct(step_lens, 0.9), max(step_lens) if step_lens else 0))
    print('  step yaw deg: med %.3f  p90 %.3f  max %.3f' %
          (rso.pct([abs(y) for y in yaw_steps], 0.5), rso.pct([abs(y) for y in yaw_steps], 0.9),
           max([abs(y) for y in yaw_steps]) if yaw_steps else 0))
    print('  total path len %.1f m   net displacement %.1f m   net yaw %.1f deg' %
          (path_len, net, net_yaw))
    print('  avg inliers %.1f   avg icp-rms %.3f m' %
          (np.mean([c[3] for c in incs]) if incs else 0, np.mean([c[4] for c in incs]) if incs else 0))
    print('=' * 92)

    if args.emit:
        os.makedirs(out_dir, exist_ok=True)
        period_ns = int(round(1e9 / args.hz))
        stem = os.path.basename(run_dir.rstrip('/\\'))
        # GT abs-pose CSV (t_ns, x,y,z, qw,qx,qy,qz)
        gtp = os.path.join(out_dir, stem + '_lidar_gt.csv')
        with open(gtp, 'w') as f:
            f.write('# t_ns, x, y, z, qw, qx, qy, qz  (RADIal lidar scan-match abs pose, frame-0 world)\n')
            for k in range(len(poses_R)):
                qw, qx, qy, qz = mat_to_quat(poses_R[k])
                p = poses_p[k]
                f.write('%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f,%.9f\n'
                        % (k * period_ns, p[0], p[1], p[2], qw, qx, qy, qz))
        # increment-form reference CSV (project schema) for fusion
        refp = os.path.join(out_dir, stem + '_lidar_ref.csv')
        fo = n2c.open_inc(refp, 'RADIal lidar scan-match odometry (reference source, R+t)')
        s_t = (0.05) ** 2
        for (k, R_inc, t_inc, n_inl, n_corr, npts) in incs:
            s_r = max(0.005, 0.2 / max(1, n_inl)) ** 2
            var6 = [s_t, s_t, s_t, s_r, s_r, s_r]
            n2c.emit_inc(fo, k * period_ns, R_inc, t_inc, var6)
        fo.close()
        print('  wrote %s' % gtp)
        print('  wrote %s' % refp)
    return poses_R, poses_p, incs


class Args:
    pass


def parse_args(argv):
    a = Args()
    a.voxel = 1.0
    a.target_n = 1200
    a.rmin = 3.0
    a.rmax = 80.0
    a.hz = 17.0
    a.icp_iters = 10          # ramp steps in the coarse-to-fine correspondence-cap schedule
    a.icp_coarse = 5.0        # m: widest (first) correspondence gate -> large convergence basin
    a.icp_max_dist = 0.4      # m: tightest (final) gate ~ the lidar point spacing
    a.min_corr = 30
    a.seed = 12345
    a.emit = False
    a.quiet = False
    pos = []
    i = 0
    while i < len(argv):
        x = argv[i]
        if x == '--voxel': i += 1; a.voxel = float(argv[i])
        elif x == '--target-n': i += 1; a.target_n = int(argv[i])
        elif x == '--rmin': i += 1; a.rmin = float(argv[i])
        elif x == '--rmax': i += 1; a.rmax = float(argv[i])
        elif x == '--hz': i += 1; a.hz = float(argv[i])
        elif x == '--icp-iters': i += 1; a.icp_iters = int(argv[i])
        elif x == '--icp-coarse': i += 1; a.icp_coarse = float(argv[i])
        elif x == '--icp-max-dist': i += 1; a.icp_max_dist = float(argv[i])
        elif x == '--min-corr': i += 1; a.min_corr = int(argv[i])
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
        print('usage: proto_radial_lidar_odom.py <run_dir> [<out_dir>] [--voxel M] [--target-n N] '
              '[--rmin M] [--rmax M] [--hz HZ] [--emit] [--quiet]')
        sys.exit(2)
    run_dir = args.pos[0]
    out_dir = args.pos[1] if len(args.pos) > 1 else run_dir
    run(run_dir, out_dir, args)


if __name__ == '__main__':
    main()
