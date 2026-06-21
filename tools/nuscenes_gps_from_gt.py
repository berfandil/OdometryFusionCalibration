#!/usr/bin/env python3
# nuscenes_gps_from_gt.py - synthesize a NOISY GPS source from the nuScenes ego_pose GT track and wire
# it into the 7-source replay manifest as the absolute-position reference. RADIal aside done; this is
# the "noisy ego_pose-as-GPS on nuScenes" item (HANDOFF): the nuScenes 7-source fusion (CAN wheel+imu +
# 5 radar Doppler) has split_median + heading_monitor configured but the monitor is INERT (no GPS course
# to score against) and the dead-reckoned drift is unbounded (~5 m tail). A GPS-like absolute ref (a)
# BOUNDS the drift and (b) ACTIVATES the GPS-course heading monitor (Slice 19c/d) so it can discover the
# heading-grade source among the 7.
#
# WHY FAKE-GEODETIC: nuScenes ego_pose is a LOCAL metric pose (the t0-ego frame = the fusion's odom
# frame, the frame <scene>_gt.csv is in). The ofc_adapters GpsCorrection takes GEODETIC fixes
# (lat/lon/alt) and runs geodetic -> WGS-84 ECEF -> ENU(datum) -> odom. We INVERT that exactly: pick a
# datum, treat the (noised) ego x,y,z as the ENU vector, map ENU -> ECEF -> geodetic with the SAME
# WGS-84 constants the adapter uses, and emit lat/lon/alt. With the manifest's odom_from_enu = identity
# the adapter round-trips lat/lon/alt back to the noised ego x,y,z -> a position fix in the odom frame.
# (The ego frame is declared == ENU here; physical compass meaning is irrelevant since we both generate
# and consume in the same convention -- the datum just keeps the WGS-84 math non-degenerate.)
#
# Noise: zero-mean Gaussian, sigma_h horizontal (x,y) + sigma_v vertical (z), at a GPS-like rate (--hz,
# subsampled from the ~169 Hz GT). The manifest's [gps] cov_floor_m2 = sigma_h^2 gives the filter the
# honest measurement noise (the emitted fixes carry zero explicit cov -> the floor IS the R).
#
# Usage:
#   python tools/nuscenes_gps_from_gt.py <scene_dir> <scene_name> [--hz HZ] [--sigma-h M] [--sigma-v M]
#       [--datum lat,lon,alt] [--seed N] [--suffix S]
#   e.g. python tools/nuscenes_gps_from_gt.py nuscenes_run/_per_scene scene-0061 --hz 2 --sigma-h 1.5
#   Writes <scene>_gps.csv + a new manifest <scene><suffix>.ini (copy of <scene>.ini + [gps], new out=).
import math
import os
import sys

# WGS-84 (must match adapters/src/gps_correction.cpp)
A = 6378137.0
F = 1.0 / 298.257223563
E2 = F * (2.0 - F)
B = A * (1.0 - F)
EP2 = E2 / (1.0 - E2)
D2R = math.pi / 180.0


def geodetic_to_ecef(lat, lon, h):
    sl, cl = math.sin(lat), math.cos(lat)
    so, co = math.sin(lon), math.cos(lon)
    N = A / math.sqrt(1.0 - E2 * sl * sl)
    return ((N + h) * cl * co, (N + h) * cl * so, (N * (1.0 - E2) + h) * sl)


def ecef_to_enu_rot(lat0, lon0):
    # rows = East, North, Up in ECEF (same as the adapter)
    sp, cp = math.sin(lat0), math.cos(lat0)
    sl, cl = math.sin(lon0), math.cos(lon0)
    return [[-sl, cl, 0.0],
            [-sp * cl, -sp * sl, cp],
            [cp * cl, cp * sl, sp]]


def ecef_to_geodetic(x, y, z):
    # Bowring closed-form inverse.
    p = math.hypot(x, y)
    th = math.atan2(z * A, p * B)
    st, ct = math.sin(th), math.cos(th)
    lat = math.atan2(z + EP2 * B * st ** 3, p - E2 * A * ct ** 3)
    lon = math.atan2(y, x)
    N = A / math.sqrt(1.0 - E2 * math.sin(lat) ** 2)
    alt = p / math.cos(lat) - N
    return lat, lon, alt


def enu_to_geodetic(enu, datum_ecef, R):
    # ecef = datum_ecef + R^T @ enu  (R = ecef->enu so R^T = enu->ecef)
    ex = datum_ecef[0] + R[0][0] * enu[0] + R[1][0] * enu[1] + R[2][0] * enu[2]
    ey = datum_ecef[1] + R[0][1] * enu[0] + R[1][1] * enu[1] + R[2][1] * enu[2]
    ez = datum_ecef[2] + R[0][2] * enu[0] + R[1][2] * enu[1] + R[2][2] * enu[2]
    return ecef_to_geodetic(ex, ey, ez)


def gauss(rng):
    return rng.gauss(0.0, 1.0)


def main():
    argv = sys.argv[1:]
    if len(argv) < 2:
        print('usage: nuscenes_gps_from_gt.py <scene_dir> <scene_name> [--hz HZ] [--sigma-h M] '
              '[--sigma-v M] [--datum lat,lon,alt] [--seed N] [--suffix S]')
        sys.exit(2)
    scene_dir, scene = argv[0], argv[1]
    hz = 2.0
    sigma_h = 1.5
    sigma_v = 3.0
    datum = (42.336, -71.057, 0.0)   # Boston (a nuScenes city); any non-polar datum works
    seed = 12345
    suffix = '_gps'
    suppress_kappa = 0.8   # correction_rot_suppress_kappa added to [global] (0=off) -> kills the
                           # heading kick a large-residual position fix injects via the trans-rot
                           # cross-cov (the urban recipe; without it a noisy GPS wrecks the heading)
    i = 2
    while i < len(argv):
        x = argv[i]
        if x == '--hz': i += 1; hz = float(argv[i])
        elif x == '--sigma-h': i += 1; sigma_h = float(argv[i])
        elif x == '--sigma-v': i += 1; sigma_v = float(argv[i])
        elif x == '--datum': i += 1; datum = tuple(float(v) for v in argv[i].split(','))
        elif x == '--seed': i += 1; seed = int(argv[i])
        elif x == '--suffix': i += 1; suffix = argv[i]
        elif x == '--suppress-kappa': i += 1; suppress_kappa = float(argv[i])
        i += 1

    import random
    rng = random.Random(seed)

    lat0, lon0, alt0 = datum
    datum_ecef = geodetic_to_ecef(lat0 * D2R, lon0 * D2R, alt0)
    R = ecef_to_enu_rot(lat0 * D2R, lon0 * D2R)

    gt_path = os.path.join(scene_dir, scene + '_gt.csv')
    gps_path = os.path.join(scene_dir, scene + '_gps.csv')
    rows = []
    with open(gt_path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in '#/':
                continue
            p = s.split(',')
            rows.append((int(p[0]), float(p[1]), float(p[2]), float(p[3])))

    period_ns = int(round(1e9 / hz))
    last_t = None
    n = 0
    with open(gps_path, 'w', newline='') as g:
        g.write('# t_ns, lat_deg, lon_deg, alt_m  (synthetic noisy GPS from nuScenes ego_pose; '
                'sigma_h=%.2f m sigma_v=%.2f m hz=%.1f datum=%.4f,%.4f)\n'
                % (sigma_h, sigma_v, hz, lat0, lon0))
        for (t_ns, x, y, z) in rows:
            if last_t is not None and (t_ns - last_t) < period_ns:
                continue
            last_t = t_ns
            enu = (x + sigma_h * gauss(rng), y + sigma_h * gauss(rng), z + sigma_v * gauss(rng))
            lat, lon, alt = enu_to_geodetic(enu, datum_ecef, R)
            g.write('%d,%.12f,%.12f,%.6f\n' % (t_ns, lat / D2R, lon / D2R, alt))
            n += 1
    print('wrote %s  (%d fixes @ %.1f Hz, sigma_h=%.2f sigma_v=%.2f)' % (gps_path, n, hz, sigma_h, sigma_v))

    # --- patch the manifest: insert [gps] before [gt], new out= ---
    ini_in = os.path.join(scene_dir, scene + '.ini')
    ini_out = os.path.join(scene_dir, scene + suffix + '.ini')
    gps_block = [
        '[gps]',
        'csv = %s_gps.csv' % scene,
        'datum_lat_deg = %.6f' % lat0,
        'datum_lon_deg = %.6f' % lon0,
        'datum_alt_m = %.3f' % alt0,
        'odom_from_enu_yaw = 0.0',
        'lever_x = 0.0', 'lever_y = 0.0', 'lever_z = 0.0',
        'cov_floor_m2 = %.4f' % (sigma_h * sigma_h),
        '',
    ]
    out_lines = []
    inserted = False
    have_suppress = False
    in_global = False
    for line in open(ini_in):
        ls = line.rstrip('\n')
        st = ls.strip()
        if st == '[gt]' and not inserted:
            out_lines.extend(gps_block)
            inserted = True
        if st.startswith('out =') or st.startswith('out='):
            ls = 'out = %s%s_out.csv' % (scene, suffix)
        if st.startswith('correction_rot_suppress_kappa'):
            have_suppress = True
            if suppress_kappa > 0:
                ls = 'correction_rot_suppress_kappa = %.4f' % suppress_kappa
        # inject suppress_kappa at the end of [global] if the source manifest lacks it
        if in_global and (st.startswith('[') and st != '[global]'):
            if not have_suppress and suppress_kappa > 0:
                out_lines.append('correction_rot_suppress_kappa = %.4f' % suppress_kappa)
                have_suppress = True
            in_global = False
        if st == '[global]':
            in_global = True
        out_lines.append(ls)
    if not inserted:
        out_lines.extend([''] + gps_block)
    with open(ini_out, 'w', newline='') as f:
        f.write('\n'.join(out_lines) + '\n')
    print('wrote %s' % ini_out)


if __name__ == '__main__':
    main()
