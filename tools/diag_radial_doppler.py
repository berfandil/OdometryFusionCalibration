"""Calibrate the radar Doppler sign + the true frame rate (hz) by fitting the radar ego-velocity
(measured radial velocity vs azimuth) and comparing to the lidar per-frame displacement.

For a STATIC world point at azimuth az, the measured radial velocity is
    v_radial = -(vx_ego*cos(az) + vy_ego*sin(az))          (2D, sensor frame)
Least-squares fit [cos,sin] -> -(vx,vy)_ego from each radar frame's detections (RANSAC-ish via median
trim). Then the implied ego SPEED (m/frame = speed*dt) should match the lidar |displacement|/frame, which
fixes dt = lidar_disp / radar_speed_per_sec; hz = 1/dt. The fit sign vs lidar forward gives doppler_sign.
"""
import os, sys, math
import numpy as np
import proto_radial_lidar_odom as plo

run = sys.argv[1]
radar_dir = os.path.join(run, 'radar')
files = sorted(f for f in os.listdir(radar_dir) if f.endswith('.npy'))

# lidar per-frame displacement (m/frame) from ICP
la = plo.parse_args([run, '--voxel', '1.0', '--target-n', '1200'])
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    poses_R, poses_p, incs = plo.run(run, run, la)
disp = {}
for (k, R_inc, t_inc, n_inl, rms, npts) in incs:
    disp[k] = math.sqrt(t_inc[0]**2 + t_inc[1]**2 + t_inc[2]**2)

def fit_ego(arr):
    # robust LS fit of v_radial = a*cos(az)+b*sin(az); return (a,b)=-(vx,vy), inlier frac.
    x, y, v = arr[5], arr[6], arr[8]
    az = np.arctan2(y, x)
    A = np.stack([np.cos(az), np.sin(az)], 1)
    sol, *_ = np.linalg.lstsq(A, v, rcond=None)
    res = v - A @ sol
    inl = np.abs(res) < (2.0 * np.median(np.abs(res)) + 0.5)
    if inl.sum() > 10:
        sol, *_ = np.linalg.lstsq(A[inl], v[inl], rcond=None)
    return sol, inl.mean()

print('frame | lidar_disp(m) | radar fit (a=-vx, b=-vy) m/s | inl% | implied_vx(=-a) ')
ratios = []
for k in range(2, len(files), 8):
    arr = np.load(os.path.join(radar_dir, files[k]), allow_pickle=True)
    (a, b), inl = fit_ego(arr)
    vx = -a   # ego forward velocity (m/s) if sign convention v_radial=-(v.bearing)
    print('%5d | %12.3f | a=%+7.2f b=%+7.2f | %4.0f%% | vx=%+7.2f' %
          (k, disp.get(k, float('nan')), a, b, 100*inl, vx))
    if disp.get(k, 0) > 0.2 and abs(vx) > 0.5:
        ratios.append(disp[k] / abs(vx))   # = dt (s/frame)

if ratios:
    dt = np.median(ratios)
    print('\nmedian dt = lidar_disp / radar_vx = %.4f s/frame  -> hz = %.2f' % (dt, 1.0/dt))
    print('(sign: radar vx should be POSITIVE forward; if a>0 consistently, doppler_sign=+1 with'
          ' v_radial=-(v.bearing). If implied vx is negative, flip doppler_sign.)')
