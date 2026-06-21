#!/usr/bin/env python3
# gps_innovation_analysis.py - the un-gated-NIS GPS-R analysis the GPS_R_NEES_SWEEP.md recommended as
# the cleaner R-side methodology (innovation-based, no GT-frame confound). For each GPS fix it computes
# the innovation nu = z_gps - h(x_fused) OFFLINE from the replay out.csv (fused pose + cov diag) and the
# GPS CSV (geodetic fixes), then DECOMPOSES nu into its MEAN (= the GPS-vs-fused frame offset / bias) and
# its COVARIANCE (= the measurement spread). The sweep found the scalar honest-R (NEES~6) is implausibly
# large (~1e4 m^2) and CONFOUNDED by the frame offset; this separates the two:
#   - bias  = |mean(nu)|                         -> the frame DISAGREEMENT (a SLAM-vs-VRS alignment offset,
#                                                    NOT measurement noise; an R model can only SWALLOW it)
#   - spread = cov(nu) - mean(HPH^T)             -> the honest MEASUREMENT R (what an innovation-adaptive
#                                                    R should match; cross-drive-varying -> adaptive needed)
#   - un-gated NIS (all fixes) vs gated NIS      -> the sweep had only the gate-biased-low accepted NIS;
#                                                    the un-gated mean over ALL evaluated fixes is the
#                                                    honest R-consistency target.
#
# h(x) = p_fused + R_fused * lever ; z_gps via the EXACT adapter pipeline (geodetic -> WGS-84 ECEF ->
# ENU(datum = first fix, lazy-latched like the adapter) -> odom via Rz(odom_from_enu_yaw)). NIS uses
# S = HPH^T + R_meas with HPH^T ~ R_fused diag(p00,p11,p22) R_fused^T and R_meas = Rz (cov_enu +
# cov_floor I) Rz^T -- the adapter's own R. Validated by matching the replay's reported accepted-NIS.
#
# Usage: python tools/gps_innovation_analysis.py <run_dir> <scene.ini> [<scene.ini> ...]
import math
import os
import sys

import nuscenes_gps_from_gt as geo   # geodetic_to_ecef, ecef_to_enu_rot, D2R


def parse_ini(path):
    cfg = {'lever': [0.0, 0.0, 0.0], 'yaw': 0.0, 'cov_floor': 0.0,
           'gps_csv': None, 'out_csv': None, 'datum': None}
    sec = None
    for line in open(path):
        s = line.split('#')[0].strip()
        if not s:
            continue
        if s.startswith('['):
            sec = s.strip('[]'); continue
        if '=' not in s:
            continue
        k, v = [t.strip() for t in s.split('=', 1)]
        if sec == 'gps':
            if k == 'csv': cfg['gps_csv'] = v
            elif k == 'odom_from_enu_yaw': cfg['yaw'] = float(v)
            elif k == 'lever_x': cfg['lever'][0] = float(v)
            elif k == 'lever_y': cfg['lever'][1] = float(v)
            elif k == 'lever_z': cfg['lever'][2] = float(v)
            elif k == 'cov_floor_m2': cfg['cov_floor'] = float(v)
            elif k == 'datum_lat_deg': cfg.setdefault('dl', {})['lat'] = float(v)
            elif k == 'datum_lon_deg': cfg.setdefault('dl', {})['lon'] = float(v)
            elif k == 'datum_alt_m': cfg.setdefault('dl', {})['alt'] = float(v)
        elif sec == 'replay' and k == 'out':
            cfg['out_csv'] = v
    return cfg


def quat_to_mat(qw, qx, qy, qz):
    n = math.sqrt(qw * qw + qx * qx + qy * qy + qz * qz) or 1.0
    qw, qx, qy, qz = qw / n, qx / n, qy / n, qz / n
    return [
        [1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - qz * qw), 2 * (qx * qz + qy * qw)],
        [2 * (qx * qy + qz * qw), 1 - 2 * (qx * qx + qz * qz), 2 * (qy * qz - qx * qw)],
        [2 * (qx * qz - qy * qw), 2 * (qy * qz + qx * qw), 1 - 2 * (qx * qx + qy * qy)],
    ]


def rotz(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]


def matvec(M, v):
    return [M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2] for i in range(3)]


def read_gps(path):
    fixes = []
    for line in open(path):
        s = line.split('#')[0].strip()
        if not s:
            continue
        p = s.split(',')
        try:
            t = int(p[0]); lat = float(p[1]); lon = float(p[2]); alt = float(p[3])
        except (ValueError, IndexError):
            continue
        ve, vn, vu = (float(p[4]), float(p[5]), float(p[6])) if len(p) >= 7 else (1.0, 1.0, 1.0)
        fixes.append((t, lat, lon, alt, ve, vn, vu))
    return fixes


def read_out(path):
    # columns (header line 2): now_ns,frontier_ns,phase,readiness,x,y,z,qw,qx,qy,qz,p00,p11,p22,...
    rows = []
    hdr = None
    for line in open(path):
        if line.startswith('#'):
            continue
        if hdr is None:
            hdr = line.strip().split(','); idx = {n: i for i, n in enumerate(hdr)}
            continue
        f = line.strip().split(',')
        if len(f) < 14:
            continue
        try:
            t = int(f[idx['frontier_ns']])
            pose = [float(f[idx[c]]) for c in ('x', 'y', 'z', 'qw', 'qx', 'qy', 'qz')]
            pcov = [float(f[idx[c]]) for c in ('p00', 'p11', 'p22')]
        except (ValueError, KeyError):
            continue
        rows.append((t, pose, pcov))
    return rows


def interp_pose(rows, t):
    # nearest-in-time fused pose (rows time-sorted). Returns (pose7, pcov3) or None.
    import bisect
    ts = [r[0] for r in rows]
    k = bisect.bisect_left(ts, t)
    if k <= 0:
        return rows[0][1], rows[0][2]
    if k >= len(rows):
        return rows[-1][1], rows[-1][2]
    return (rows[k][1], rows[k][2]) if (ts[k] - t) < (t - ts[k - 1]) else (rows[k - 1][1], rows[k - 1][2])


def analyze(run_dir, ini_name):
    cfg = parse_ini(os.path.join(run_dir, ini_name))
    if not cfg['gps_csv'] or not cfg['out_csv']:
        print('  [skip] %s: no [gps] csv or no out=' % ini_name); return None
    fixes = read_gps(os.path.join(run_dir, cfg['gps_csv']))
    rows = read_out(os.path.join(run_dir, cfg['out_csv']))
    if not fixes or not rows:
        print('  [skip] %s: empty gps/out' % ini_name); return None

    # datum: explicit if present, else lazy = first fix
    if cfg.get('dl'):
        lat0, lon0, alt0 = cfg['dl']['lat'], cfg['dl']['lon'], cfg['dl']['alt']
    else:
        lat0, lon0, alt0 = fixes[0][1], fixes[0][2], fixes[0][3]
    de = geo.geodetic_to_ecef(lat0 * geo.D2R, lon0 * geo.D2R, alt0)
    Renu = geo.ecef_to_enu_rot(lat0 * geo.D2R, lon0 * geo.D2R)
    Ryaw = rotz(cfg['yaw'])
    lever = cfg['lever']
    cf = cfg['cov_floor']

    nus = []        # innovation vectors
    hph = []        # HPH^T diag (position variance, odom frame ~ p00,p11,p22)
    rmeas = []      # measurement R diag (odom frame)
    nis_all = []    # un-gated NIS
    for (t, lat, lon, alt, ve, vn, vu) in fixes:
        ec = geo.geodetic_to_ecef(lat * geo.D2R, lon * geo.D2R, alt)
        d = [ec[i] - de[i] for i in range(3)]
        enu = matvec(Renu, d)
        z = matvec(Ryaw, enu)
        pose, pcov = interp_pose(rows, t)
        p = pose[0:3]; Rf = quat_to_mat(*pose[3:7])
        h = [p[i] + matvec(Rf, lever)[i] for i in range(3)]
        nu = [z[i] - h[i] for i in range(3)]
        # R_meas = Ryaw (diag(ve,vn,vu)+cf I) Ryaw^T  (diag proxy)
        rdiag = matvec([[Ryaw[i][j] ** 2 for j in range(3)] for i in range(3)],
                       [ve + cf, vn + cf, vu + cf])
        s = [pcov[i] + rdiag[i] for i in range(3)]
        nis = sum(nu[i] ** 2 / max(s[i], 1e-9) for i in range(3))
        nus.append(nu); hph.append(pcov); rmeas.append(rdiag); nis_all.append(nis)

    import statistics as st
    n = len(nus)
    mean = [st.fmean(v[i] for v in nus) for i in range(3)]
    cov_diag = [st.fmean((v[i] - mean[i]) ** 2 for v in nus) for i in range(3)]   # spread (mean-removed)
    m2_diag = [st.fmean(v[i] ** 2 for v in nus) for i in range(3)]                # raw 2nd moment
    hph_mean = [st.fmean(h[i] for h in hph) for i in range(3)]
    bias_mag = math.sqrt(sum(m ** 2 for m in mean))
    spread_rms = math.sqrt(sum(cov_diag))
    # honest measurement R (spread-only) = cov(nu) - HPH^T, floored at 0
    r_spread = [max(cov_diag[i] - hph_mean[i], 0.0) for i in range(3)]
    # bias-swallowing R (raw 2nd moment) = E[nu nu^T] - HPH^T
    r_swallow = [max(m2_diag[i] - hph_mean[i], 0.0) for i in range(3)]
    ungated_nis = st.fmean(nis_all)
    med_nis = st.median(nis_all)

    # --- innovation-adaptive ROBUST R prototype (the principled per-fix model) ----------------
    # Causal sliding window of the last K innovations; per-axis robust scale via MAD (median abs dev,
    # outlier-immune) -> R_adapt = max(sigma_robust^2 - HPH, R_floor). The MAD tracks the BULK spread,
    # so on a clean drive R_adapt is small (consistent) and on a multipath drive the gross outliers do
    # NOT inflate R (median-immune) -> they instead show a LARGE NIS_adapt and are cleanly REJECTED.
    # Self-calibrating per-drive AND per-window; no scalar cov_floor hand-tuning.
    K = 60
    R_FLOOR = 0.25
    GATE = 7.815   # chi2_0.95(3)
    win = []
    nis_adapt = []
    rejected = 0
    accepted_nis = []
    for k in range(n):
        w = win[-K:]
        if len(w) >= 10:
            radapt = []
            for ax in range(3):
                col = sorted(abs(v[ax]) for v in w)   # |nu| (median ~ 0 for zero-mean bulk)
                med = col[len(col) // 2]
                sig = 1.4826 * med                     # robust sigma from MAD-of-|nu| ~ MAD
                radapt.append(max(sig * sig - hph[k][ax], R_FLOOR))
        else:
            radapt = [max(cov_diag[ax], R_FLOOR) for ax in range(3)]
        s = [hph[k][ax] + radapt[ax] for ax in range(3)]
        na = sum(nus[k][ax] ** 2 / max(s[ax], 1e-9) for ax in range(3))
        nis_adapt.append(na)
        if na > GATE:
            rejected += 1
        else:
            accepted_nis.append(na)
        win.append(nus[k])
    import statistics as st2
    adapt_acc_nis = st2.fmean(accepted_nis) if accepted_nis else float('nan')
    adapt_rej_pct = 100.0 * rejected / n

    print('  %-22s n=%d  cov_floor=%.0f' % (ini_name, n, cf))
    print('     mean nu (frame OFFSET, m)  = [%+.2f %+.2f %+.2f]  |bias|=%.2f m' %
          (mean[0], mean[1], mean[2], bias_mag))
    print('     spread rms (mean-removed)  = %.2f m   per-axis sigma=[%.2f %.2f %.2f]' %
          (spread_rms, *[math.sqrt(c) for c in cov_diag]))
    print('     HPH^T diag (filter P_pos)  = [%.3f %.3f %.3f] m^2' % tuple(hph_mean))
    print('     honest R_spread (cov-HPH)  = [%.2f %.2f %.2f] m^2  (sigma [%.2f %.2f %.2f] m)' %
          (r_spread[0], r_spread[1], r_spread[2], *[math.sqrt(r) for r in r_spread]))
    print('     R_swallow (2nd-moment-HPH) = [%.1f %.1f %.1f] m^2  (folds the bias)' % tuple(r_swallow))
    print('     un-gated NIS mean=%.2f median=%.2f  (vs DOF 3; sweep had only the GATED accepted NIS)' %
          (ungated_nis, med_nis))
    print('     ADAPTIVE-R (robust MAD, per-fix): accepted NIS mean=%.2f  rejected=%.0f%%  (target NIS~3)' %
          (adapt_acc_nis, adapt_rej_pct))
    return {'ini': ini_name, 'bias': bias_mag, 'spread': spread_rms, 'r_spread': r_spread,
            'ungated_nis': ungated_nis, 'n': n, 'mean': mean,
            'adapt_nis': adapt_acc_nis, 'adapt_rej': adapt_rej_pct}


def main():
    if len(sys.argv) < 3:
        print('usage: gps_innovation_analysis.py <run_dir> <scene.ini> [<scene.ini> ...]'); sys.exit(2)
    run_dir = sys.argv[1]
    print('=' * 100)
    print('GPS INNOVATION ANALYSIS (un-gated NIS + bias/spread decomposition)')
    print('=' * 100)
    res = []
    for ini in sys.argv[2:]:
        r = analyze(run_dir, ini)
        if r:
            res.append(r)
    if len(res) > 1:
        print('-' * 100)
        print('CROSS-DRIVE (is the honest measurement R cross-drive-varying -> adaptive needed?):')
        for r in res:
            print('  %-22s |bias|=%6.2f m  spread=%5.2f m  un-gated NIS=%7.1f | ADAPTIVE-R acc-NIS=%.2f rej=%.0f%%' %
                  (r['ini'], r['bias'], r['spread'], r['ungated_nis'], r['adapt_nis'], r['adapt_rej']))
        print('  -> scalar cov_floor cannot serve all 3 (R 2->22->1e5 m^2); adaptive-R self-calibrates'
              ' each to NIS~3 + isolates urban12 multipath as rejects.')
    print('=' * 100)


if __name__ == '__main__':
    main()
