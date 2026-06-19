"""Per-pair radar scan-match diagnostic: for each consecutive frame pair, report matched corr,
RANSAC inliers, and the recovered per-step yaw in 3D vs 2D (z-zeroed) vs the lidar GT yaw. Separates
the two failure suspects: poor scatterer PERSISTENCE (low inliers) vs coarse ELEVATION (3D ~ 2D)."""
import os, sys, math, io, contextlib
import numpy as np
import radar_scan_odometry as rso
import synth_4d_radar as s4d
import proto_radial_lidar_odom as plo
import radial_to_csv as r2c

run = sys.argv[1]
NP = int(sys.argv[2]) if len(sys.argv) > 2 else 300
radar_dir = os.path.join(run, 'radar')
files = sorted(f for f in os.listdir(radar_dir) if f.endswith('.npy'))

# lidar GT deltas
la = plo.parse_args([run, '--voxel', '1.0', '--target-n', '1200'])
with contextlib.redirect_stdout(io.StringIO()):
    poses_R, poses_p, incs = plo.run(run, run, la)
gt_yaw = {}
for (k, R_inc, t_inc, n_inl, rms, npts) in incs:
    gt_yaw[k] = math.degrees(math.atan2(R_inc[1][0], R_inc[0][0]))

import random
rng = random.Random(7)

def dets_of(fn):
    return r2c.read_radar_dets(os.path.join(radar_dir, fn), 1.0, 100.0, 70.0, 1.0, NP)

def yaw_recover(prev_xyz, curr_xyz, use3d):
    if not use3d:
        prev_xyz = [(x, y, 0.0) for (x, y, z) in prev_xyz]
        curr_xyz = [(x, y, 0.0) for (x, y, z) in curr_xyz]
    dp = rso.build_descriptors(prev_xyz)
    dc = rso.build_descriptors(curr_xyz)
    pairs = s4d.match_frames_fast(dc, dp, 0.5, 0.5)
    R, t, ni, nc = rso.ransac_kabsch(curr_xyz, prev_xyz, pairs, rng, 1.0, 200, 6, dim=3)
    if R is None:
        return None, 0, len(pairs)
    return math.degrees(math.atan2(R[1, 0], R[0, 0])), ni, len(pairs)

prev = None
rows = []
print('pair |  corr inl | yaw3d  yaw2d   GTyaw')
e3, e2, g = [], [], []
for k, fn in enumerate(files):
    d = dets_of(fn)
    xyz = [(x, y, z) for (x, y, z, vx, vy, vz) in d]
    if prev is not None:
        y3, ni3, nc = yaw_recover(prev, xyz, True)
        y2, ni2, _ = yaw_recover(prev, xyz, False)
        gy = gt_yaw.get(k, float('nan'))
        if y3 is not None and k % 2 == 0:
            print('%4d | %4d %4d | %+6.2f %+6.2f  %+6.2f' % (k, nc, ni3, y3, y2 if y2 else 0, gy))
        if y3 is not None and not math.isnan(gy):
            e3.append(y3); e2.append(y2 if y2 is not None else 0); g.append(gy)
    prev = xyz

print('\nyaw correlation vs GT:  3D = %.3f   2D = %.3f   (n=%d)'
      % (rso.pearson(e3, g), rso.pearson(e2, g), len(g)))
print('GT yaw range: [%.2f, %.2f]  std %.2f' % (min(g), max(g), float(np.std(g))))
print('3D recovered yaw std %.2f   2D %.2f' % (float(np.std(e3)), float(np.std(e2))))
