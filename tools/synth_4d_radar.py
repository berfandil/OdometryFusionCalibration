#!/usr/bin/env python3
# synth_4d_radar.py - synthesize a DENSE 3D / 4D (3D position + Doppler) radar from a REAL ground-truth
# trajectory, then run the descriptor scan-matching odometry (RADAR_SCAN_ODOMETRY.md) over it to test
# THE OPEN QUESTION from the radar thread: on the real sparse 2D ARS408 (~9 Kabsch inliers) the per-step
# ROTATION sat BELOW the Kabsch noise floor (yaw correlation ~0; cumulative heading scattered). Section 4
# of RADAR_SCAN_ODOMETRY.md claims "a 4D/imaging radar is the fix: more detections + elevation -> a far
# better-conditioned 3D Kabsch". THIS TOOL CHECKS THAT WITH NUMBERS by sweeping the detections-per-scan
# count N and measuring whether rotation becomes USABLE (yaw correlation -> 1, cumulative heading tracks
# GT) as N grows, and at what N it crosses the noise floor. Direct contrast: same scene, dense 3D radar
# vs the real sparse 2D radar.
#
# WHY SYNTHETIC: the nuScenes radar is a real 2D ARS408 (z=0, vx/vy only) -- there is no real 4D radar in
# the dataset. We synthesize one whose MOTION is the REAL ego trajectory (scene-0061, which turns ~98 deg)
# so the per-frame yaw is realistic (straight ~0.2 deg/scan, turning ~1-2 deg/scan) and the only thing
# that changes vs the real radar is the sensor DENSITY + the 3rd dimension. This isolates the variable.
#
# PIPELINE -----------------------------------------------------------------------------------------------
#  1. Read a GT abs-pose CSV (t_ns, x,y,z, qw,qx,qy,qz) in the scene's t0-ego world frame (the frame that
#     nuscenes_to_csv.py emits). T_ego(t) = (R_ego, p_ego).
#  2. Plant a radar mount X_radar = ego_from_radar = (R_x, t_x) (default yaw 0, lever [3.4,0,0.5], close to
#     RADAR_FRONT's real calibrated_sensor 0.2 deg / [3.412,0,0.5]). World radar pose:
#         R_radar(t) = R_ego(t) R_x ,  p_radar(t) = p_ego(t) + R_ego(t) t_x .
#  3. A FIXED static 3D landmark cloud of M world points in a shell/volume spanning the trajectory bbox +
#     a realistic elevation spread (seeded; generated ONCE, the same world for every scan).
#  4. At each radar scan time (GT subsampled to ~radar_hz): map every landmark P into the radar frame
#         r = R_radar(t)^T (P - p_radar(t)) ,
#     FOV-gate (range in [r_min, r_max], |azimuth| <= az, |elevation| <= el), CAP/sample to N detections
#     (seeded), add zero-mean 3D position noise N(0, sigma_pos^2). -> the static detection set (positions).
#  5. Per detection: the 3D radial velocity a STATIC world point shows is the radial component of
#     v_static = -(v_radar + omega_radar x r) from the GT radar twist (v_radar, omega_radar) over the local
#     interval; + velocity noise N(0, sigma_vel^2). We hand the scan-matcher a full (vx,vy,vz) = radial
#     speed * bearing_unit (a real radar measures only the radial component; this matches the static
#     filter's own radial projection). --movers (default 0) optionally injects moving detections with a
#     velocity offset for a dirtier static-filter test; 0 = the clean rotation test.
#  6. Run the 3D scan-match (radar_scan_odometry.py 3D path: 3D static filter + 3D descriptors + 3D Kabsch
#     wrapped in RANSAC, minimal sample = 3 pairs) over the synthetic scans, baseline N consecutive (k vs
#     k-baseline). The static-filter PRIOR is the clean GT radar twist (a perfect reference odometry -- the
#     realistic best case; with movers=0 the filter is nearly a no-op anyway).
#
# METRICS (per N, the HEADLINE) --------------------------------------------------------------------------
#   recovered B_radar (3D R, t) vs GT B_radar_gt over each [t_{k-base}, t_k]:
#     - rotation error (deg) median + p90 (geodesic SO(3) error),
#     - per-step YAW CORRELATION recovered-vs-GT (the headline: ~0/negative for the real 2D radar -- does
#       it rise toward 1 with density?),
#     - cumulative-heading track: recovered net yaw vs GT net yaw (does the integrated heading track?),
#     - translation error (m) median + p90,
#     - the Kabsch inlier count (the conditioning proxy; does the noise floor drop ~1/sqrt(N)?).
#
# REUSE: imports tools/radar_scan_odometry.py for the descriptor / matching / 3D Kabsch / RANSAC / 3D
# static filter / validation-metric core, and tools/nuscenes_to_csv.py for the SE3/quat helpers. Only the
# scan SYNTHESIS + the sweep harness are new here. Pure ASCII; numpy; SEEDED (numpy default_rng + the
# scan-matcher's random.Random) -> reproducible.
#
# Usage:
#   python tools/synth_4d_radar.py <gt_csv> <out_dir> [--num-points N] [--sigma-pos M] [--sigma-vel M_S]
#       [--seed N] [--movers K] [--landmarks M] [--radar-hz HZ] [--baseline N] [--emit-csv]
#       [--mount yaw,pitch,roll,lx,ly,lz] [--sweep N0,N1,...] [--quiet]
#   e.g. (single N):  python tools/synth_4d_radar.py nuscenes_run/scene-0061_gt.csv nuscenes_run --num-points 250
#        (full sweep): python tools/synth_4d_radar.py nuscenes_run/scene-0061_gt.csv nuscenes_run --sweep 20,50,100,250,500
#   With --emit-csv, writes the recovered radar odometry (full R+t increment-form CSV, the project schema)
#   to <out_dir>/<gtstem>_synth4d_nN.csv for the calibration stretch.
import math
import os
import random
import sys

import numpy as np

import nuscenes_to_csv as n2c
import radar_scan_odometry as rso


# ---------- fast EXACT matcher (same correspondences as rso.match_frames, ~2.4x faster) ----------
def _merge_count_list(a, b, tau):
    # The SAME tolerance-merge consume-once two-pointer as rso.merge_count, on plain Python lists
    # (no numpy scalar overhead in the hot inner loop). Identical output, just faster at large N.
    i = j = 0
    na, nb = len(a), len(b)
    c = 0
    while i < na and j < nb:
        d = a[i] - b[j]
        if d < -tau:
            i += 1
        elif d > tau:
            j += 1
        else:
            c += 1
            i += 1
            j += 1
    return c


def match_frames_fast(descs_curr, descs_prev, tau_d, rho):
    # Drop-in for rso.match_frames producing the IDENTICAL correspondence set (mutual-best +
    # rho*min-len floor over the agreeing-distance count), but with the inner merge on Python lists so
    # the dense-radar sweep (N up to 500) is tractable. Verified set-identical to rso.match_frames.
    nc, npv = len(descs_curr), len(descs_prev)
    if nc == 0 or npv == 0:
        return []
    dc = [d.tolist() for d in descs_curr]
    dp = [d.tolist() for d in descs_prev]
    lc = [len(d) for d in dc]
    lp = [len(d) for d in dp]
    S = np.zeros((nc, npv), dtype=np.int32)
    for i in range(nc):
        di = dc[i]
        if not di:
            continue
        row = S[i]
        for k in range(npv):
            dk = dp[k]
            if not dk:
                continue
            row[k] = _merge_count_list(di, dk, tau_d)
    best_prev = np.argmax(S, axis=1)
    best_curr = np.argmax(S, axis=0)
    pairs = []
    for i in range(nc):
        k = int(best_prev[i])
        if int(best_curr[k]) != i:
            continue
        floor = rho * min(lc[i], lp[k])
        if S[i, k] >= floor and S[i, k] > 0:
            pairs.append((i, k))
    return pairs


# ---------- GT trajectory loading ----------
def read_gt(path):
    # Read a GT abs-pose CSV (t_ns, x,y,z, qw,qx,qy,qz). Returns (ts_us, R_list, p_list) with ts in
    # MICROSECONDS (the scan-matcher's validation uses us internally) and R as 3x3 plain lists.
    ts, Rs, ps = [], [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#/":
                continue
            p = s.replace(",", " ").split()
            t_ns = int(p[0])
            x, y, z = float(p[1]), float(p[2]), float(p[3])
            qw, qx, qy, qz = float(p[4]), float(p[5]), float(p[6]), float(p[7])
            ts.append(t_ns // 1000)                 # ns -> us
            Rs.append(n2c.quat_wxyz_to_mat(qw, qx, qy, qz))
            ps.append([x, y, z])
    return ts, Rs, ps


def mount_from_spec(spec):
    # spec = "yaw,pitch,roll,lx,ly,lz" in (radians, metres) -> (R_x 3x3, t_x 3). Default close to
    # RADAR_FRONT's real calibrated_sensor (yaw 0.2 deg, lever [3.412,0,0.5]).
    y, p, r, lx, ly, lz = [float(v) for v in spec.split(",")]
    Rz = n2c.rot_z(y)
    # Ry, Rx (small/zero by default) -- build via so3 about the axes for generality
    cp, sp = math.cos(p), math.sin(p)
    Ry = [[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]]
    cr, sr = math.cos(r), math.sin(r)
    Rx = [[1, 0, 0], [0, cr, -sr], [0, sr, cr]]
    R_x = n2c.mat_mul(n2c.mat_mul(Rz, Ry), Rx)
    return R_x, [lx, ly, lz]


# ---------- static world landmark cloud ----------
def make_landmarks(ps, n_landmarks, rng):
    # A FIXED static 3D landmark cloud spanning the trajectory bbox + a margin, with a realistic
    # elevation spread (ground-ish to ~building height). Generated ONCE -> the same world for every scan.
    xs = [p[0] for p in ps]
    ys = [p[1] for p in ps]
    zs = [p[2] for p in ps]
    margin = 60.0                                   # m beyond the path bbox (radar sees ~100 m out)
    x0, x1 = min(xs) - margin, max(xs) + margin
    y0, y1 = min(ys) - margin, max(ys) + margin
    z0, z1 = min(zs) - 2.0, max(zs) + 12.0          # ground-ish to ~building height (elevation spread)
    P = np.empty((n_landmarks, 3))
    P[:, 0] = rng.uniform(x0, x1, n_landmarks)
    P[:, 1] = rng.uniform(y0, y1, n_landmarks)
    P[:, 2] = rng.uniform(z0, z1, n_landmarks)
    return P


# ---------- per-scan synthesis ----------
def radar_pose(R_ego, p_ego, R_x, t_x):
    # World radar pose: R_radar = R_ego R_x, p_radar = p_ego + R_ego t_x.
    R_radar = n2c.mat_mul(R_ego, R_x)
    p_radar = [p_ego[i] + n2c.mat_vec(R_ego, t_x)[i] for i in range(3)]
    return R_radar, p_radar


def gt_radar_twist(R0, p0, R1, p1, R_x, t_x, dt):
    # GT radar body twist (v, omega) over [t0, t1]: B_radar_gt = X_radar^-1 o A_ego o X_radar where
    # A_ego = T0^-1 o T1 is the ego body delta. Reuses radar_scan_odometry's expected_radar_motion +
    # the 3D twist. Returns ((vx,vy,vz),(wx,wy,wz)) in the radar frame.
    R0t = n2c.mat_t(R0)
    A_R = n2c.mat_mul(R0t, R1)
    A_t = n2c.mat_vec(R0t, [p1[i] - p0[i] for i in range(3)])
    B_R, B_t = rso.expected_radar_motion(A_R, A_t, R_x, t_x)
    return rso.twist_from_se3_3d(B_R, B_t, dt) if dt > 0 else ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))


def synth_scan(P, R_radar, p_radar, twist, args, rng):
    # Synthesize ONE radar scan: transform landmarks into the radar frame, FOV-gate, cap to N, add pos
    # noise, compute the static radial velocity (+ vel noise). Returns a list of (x,y,z, vx,vy,vz)
    # detections (the format static_filter_3d consumes). `twist` = ((v),(omega)) GT radar twist.
    Rt = np.asarray(n2c.mat_t(R_radar))             # world -> radar rotation
    pr = np.asarray(p_radar)
    rel = P - pr                                    # (M,3) world offsets
    r = (Rt @ rel.T).T                              # (M,3) radar-frame positions
    rng_xy = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2 + r[:, 2] ** 2)
    # FOV gate: range, azimuth (about +x, in the x-y plane), elevation (off the x-y plane)
    az = np.arctan2(r[:, 1], r[:, 0])
    horiz = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)
    el = np.arctan2(r[:, 2], np.maximum(horiz, 1e-9))
    keep = ((rng_xy >= args.r_min) & (rng_xy <= args.r_max)
            & (np.abs(az) <= args.az) & (np.abs(el) <= args.el) & (r[:, 0] > 0.0))
    vis = r[keep]
    vis_rng = rng_xy[keep]
    if vis.shape[0] == 0:
        return []
    # CAP to N by selecting the N NEAREST visible landmarks (closest scatterers dominate radar
    # returns). This is the key to frame-to-frame PERSISTENCE: the ego moves only ~0.7 m/scan, so the
    # nearest-N set is almost the SAME world points across consecutive scans (with natural FOV turnover
    # at the range/angle boundary) -- exactly the partial-overlap the descriptor matcher is built for.
    # An independent random resample per scan would share no common neighbour structure -> no matches.
    n_keep = min(args.num_points, vis.shape[0])
    near = np.argpartition(vis_rng, n_keep - 1)[:n_keep] if n_keep < vis.shape[0] \
        else np.arange(vis.shape[0])
    sel = vis[near]
    # add INDEPENDENT 3D position noise per scan (the same world point jitters differently each scan,
    # like real radar measurement noise -> the descriptors are noisy but anchored to stable geometry)
    sel = sel + rng.normal(0.0, args.sigma_pos, size=sel.shape)
    (v, om) = twist
    dets = []
    n_mov = 0
    for j in range(sel.shape[0]):
        x, y, z = float(sel[j, 0]), float(sel[j, 1]), float(sel[j, 2])
        rho = math.sqrt(x * x + y * y + z * z)
        if rho < 1e-6:
            continue
        ux, uy, uz = x / rho, y / rho, z / rho
        # static-world velocity v_static = -(v + omega x r); radar measures its radial component
        cx = om[1] * z - om[2] * y
        cy = om[2] * x - om[0] * z
        cz = om[0] * y - om[1] * x
        vsx, vsy, vsz = -(v[0] + cx), -(v[1] + cy), -(v[2] + cz)
        radial = vsx * ux + vsy * uy + vsz * uz
        radial += rng.normal(0.0, args.sigma_vel)   # velocity noise on the radial channel
        is_mover = False
        if args.movers > 0 and n_mov < args.movers and rng.random() < (args.movers / max(1, n_keep)):
            radial += rng.uniform(2.0, 8.0) * (1 if rng.random() < 0.5 else -1)  # mover Doppler offset
            is_mover = True
            n_mov += 1
        # hand a full velocity vector whose radial component is `radial` (along the bearing)
        dets.append((x, y, z, radial * ux, radial * uy, radial * uz))
    return dets


# ---------- the scan-match over synthetic scans (3D path) ----------
def run_scan_match(gt_ts, gt_R, gt_p, P, R_x, t_x, args, seed):
    # Subsample the GT track to ~radar_hz scan times, synthesize each scan, then run the 3D descriptor
    # scan-match (baseline N) and validate recovered vs GT B_radar over each interval. Returns
    # (perframe, stats, rows) where perframe = [(R_inc, t_inc, Bgt_R, Bgt_t, recovered, ts_us, ts_prev, k)]
    # mirroring radar_scan_odometry.report_radar's input, rows = increment-CSV tuples for --emit-csv.
    # pick scan times by stepping through GT at the radar period
    period_us = int(round(1e6 / args.radar_hz))
    scan_idx = []
    last_t = None
    for i, t in enumerate(gt_ts):
        if last_t is None or (t - last_t) >= period_us:
            scan_idx.append(i)
            last_t = t
    # detection RNG (numpy) + RANSAC RNG (random.Random), both seeded
    det_rng = np.random.default_rng(seed * 7919 + 17)
    ransac_rng = random.Random(seed * 1000003 + 4242)

    N = max(1, int(args.baseline))
    stats = {"n_pairs": 0, "n_fallback": 0, "sum_det": 0, "sum_static": 0,
             "sum_corr": 0, "sum_inl": 0, "n_recov": 0, "n_gated": 0}
    perframe = []
    rows = []
    # ring buffer of the last N+1 scans' synthesized detections, keyed by ordered scan index
    hist = []  # list of (ts_us, gt_index, dets)
    for k, gi in enumerate(scan_idx):
        ts_us = gt_ts[gi]
        R_ego, p_ego = gt_R[gi], gt_p[gi]
        R_radar, p_radar = radar_pose(R_ego, p_ego, R_x, t_x)
        # the twist used to synthesize the radial velocity: GT radar twist over the LOCAL (~1 scan)
        # interval (so the static prediction matches a real per-scan Doppler). Use the previous GT pose.
        gi_prev = scan_idx[k - 1] if k > 0 else gi
        dt_local = (gt_ts[gi] - gt_ts[gi_prev]) * 1e-6
        if dt_local <= 0:
            tw = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        else:
            tw = gt_radar_twist(gt_R[gi_prev], gt_p[gi_prev], R_ego, p_ego, R_x, t_x, dt_local)
        dets = synth_scan(P, R_radar, p_radar, tw, args, det_rng)
        hist.append((ts_us, gi, dets))
        if len(hist) > N + 1:
            hist.pop(0)
        if k < N:
            continue

        t_prev_us, gi0, prev_dets = hist[0]          # the scan at ordered index k-N
        dt = (ts_us - t_prev_us) * 1e-6
        if dt <= 0:
            continue

        # static-filter PRIOR: the GT radar twist over the matched [t_{k-N}, t_k] interval (a clean
        # reference odometry -- the realistic best case). Build the expected radar motion (B_R, B_t).
        R0g, p0g = gt_R[gi0], gt_p[gi0]
        R0gt = n2c.mat_t(R0g)
        A_R = n2c.mat_mul(R0gt, R_ego)
        A_t = n2c.mat_vec(R0gt, [p_ego[i] - p0g[i] for i in range(3)])
        B_R, B_t = rso.expected_radar_motion(A_R, A_t, R_x, t_x)

        # --- A. 3D static filter both frames ---
        prev_kept = rso.static_filter_3d(prev_dets, B_R, B_t, dt, args.tau_v)
        curr_kept = rso.static_filter_3d(dets, B_R, B_t, dt, args.tau_v)
        stats["sum_det"] += len(dets)
        stats["sum_static"] += len(curr_kept)
        stats["n_pairs"] += 1

        # --- B. descriptors (dim-agnostic; 3D points) ---
        descs_prev = rso.build_descriptors(prev_kept)
        descs_curr = rso.build_descriptors(curr_kept)

        # --- C. match (fast exact matcher; identical correspondences to rso.match_frames) ---
        pairs = match_frames_fast(descs_curr, descs_prev, args.tau_d, args.rho)
        stats["sum_corr"] += len(pairs)

        # --- D. recover (3D Kabsch + RANSAC, minimal sample = 3) ---
        R3, t3, n_inl, n_corr = rso.ransac_kabsch(curr_kept, prev_kept, pairs, ransac_rng,
                                                  args.tau_inlier, args.ransac_iters,
                                                  args.min_corr, dim=3)
        stats["sum_inl"] += n_inl

        recovered = R3 is not None
        if recovered:
            # plausibility gate vs the GT prior (same spirit as the 2D path): reject blunders whose
            # geodesic rotation or translation grossly exceeds the prior + margin.
            (vp, omp) = rso.twist_from_se3_3d(B_R, B_t, dt)
            rot_pred = math.sqrt(omp[0] ** 2 + omp[1] ** 2 + omp[2] ** 2) * dt
            spd_pred = math.sqrt(vp[0] ** 2 + vp[1] ** 2 + vp[2] ** 2) * dt
            rot_inc = math.radians(rso.rot_err_deg(R3.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
            t_mag = float(np.linalg.norm(t3))
            rot_cap = rot_pred + math.radians(args.yaw_margin_deg)
            t_cap = spd_pred + args.trans_margin_m
            if rot_inc > rot_cap or t_mag > t_cap:
                stats["n_gated"] += 1
                recovered = False

        if not recovered:
            R_inc = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            t_inc = [0.0, 0.0, 0.0]
            var6 = [100.0, 100.0, 100.0, 0.25, 0.25, 0.25]
            stats["n_fallback"] += 1
        else:
            # Kabsch solved  r_prev ~= R3 r_curr + t3  over [t_{k-N}, t_k]. A static world point P obeys
            # P_prev = B P_curr where B = T_prev^-1 o T_curr is the FORWARD body delta (partner->curr) --
            # which is EXACTLY (R3, t3). So the increment is (R3, t3) DIRECTLY, no inverse. (NOTE: the 2D
            # path in radar_scan_odometry.py transposes here -- R_inc = R2^T -- a sign flip that is
            # invisible on the real sparse 2D radar because its per-step rotation is pure noise, but the
            # dense synth here EXPOSES it as a systematic yaw-sign inversion. The forward-delta mapping
            # below is the correct one, verified by a noise-free planted-yaw round-trip.)
            R_inc = [[float(R3[a, b]) for b in range(3)] for a in range(3)]
            t_inc = [float(t3[0]), float(t3[1]), float(t3[2])]
            stats["n_recov"] += 1
            s_t = (max(args.sigma_pos, 0.02)) ** 2
            s_r = max(0.01, 0.3 / max(1, n_inl)) ** 2
            var6 = [s_t, s_t, s_t, s_r, s_r, s_r]
        rows.append((ts_us, R_inc, t_inc, var6))

        # --- validation: GT radar frame-to-frame motion over [t_{k-N}, t_k] (forward body delta) ---
        Bgt_R, Bgt_t = rso.expected_radar_motion(A_R, A_t, R_x, t_x)
        perframe.append((R_inc, t_inc, Bgt_R, Bgt_t, recovered, ts_us, t_prev_us, k))

    return perframe, stats, rows


# ---------- reporting (reuses radar_scan_odometry.report_radar; adds full-3D rot-error) ----------
def summarize(perframe, stats, args, label):
    rep = rso.report_radar(label, perframe, stats, args.baseline)
    return rep


def print_sweep_header(args):
    print("=" * 118)
    print("SYNTH 4D RADAR -- descriptor scan-match 3D path  (gt=%s)" % os.path.basename(args.gt_csv))
    print("mount: yaw=%.4f deg lever=[%.3f %.3f %.3f]   sigma_pos=%.3f m sigma_vel=%.2f m/s  "
          "radar_hz=%.1f baseline=%d landmarks=%d movers=%d seed=%d"
          % (math.degrees(args.mount_yaw), args.t_x[0], args.t_x[1], args.t_x[2],
             args.sigma_pos, args.sigma_vel, args.radar_hz, args.baseline, args.landmarks,
             args.movers, args.seed))
    print("FOV: range [%.1f, %.1f] m  azimuth +-%.0f deg  elevation +-%.0f deg"
          % (args.r_min, args.r_max, math.degrees(args.az), math.degrees(args.el)))
    print("=" * 118)


def print_sweep_table(reps):
    # The HEADLINE table: per N, match quality + rotation usability.
    print("-" * 118)
    print("%6s | %6s | %7s | %7s | %7s | %8s || %9s | %9s | %9s | %9s | %10s"
          % ("N", "#pair", "static", "corr", "inlier", "fallbk%",
             "rot med", "rot p90", "yaw corr", "trans med", "cumHeadErr"))
    print("%6s | %6s | %7s | %7s | %7s | %8s || %9s | %9s | %9s | %9s | %10s"
          % ("", "", "avg", "avg", "avg", "", "deg", "deg", "(->1?)", "m", "deg"))
    print("-" * 118)
    for r in reps:
        cum_err = r["cum_est_yaw_deg"] - r["cum_gt_yaw_deg"]
        print("%6d | %6d | %7.1f | %7.1f | %7.2f | %7.1f%% || %9.3f | %9.3f | %9.3f | %9.3f | %+10.2f"
              % (r["_N"], r["n_pairs"], r["avg_static"], r["avg_corr"], r["avg_inl"],
                 r["fallback_pct"], r["rot_med"], r["rot_p90"], r["yaw_corr"],
                 r["trans_med"], cum_err))
    print("-" * 118)
    # cumulative-heading line: recovered net yaw vs GT net yaw (full integrated track)
    print("CUMULATIVE HEADING (recovered net yaw vs GT net yaw, full integrated, deg):")
    for r in reps:
        print("   N=%-4d  recovered=%+8.2f deg   GT=%+8.2f deg   err=%+8.2f deg   (nonovl err=%+8.2f deg)"
              % (r["_N"], r["cum_est_yaw_deg"], r["cum_gt_yaw_deg"],
                 r["cum_est_yaw_deg"] - r["cum_gt_yaw_deg"], r["nov_heading_err_deg"]))
    print("=" * 118)


# ---------- args ----------
class Args:
    pass


def parse_args(argv):
    a = Args()
    a.num_points = 250
    a.sigma_pos = 0.15
    a.sigma_vel = 0.1
    a.seed = 12345
    a.movers = 0
    a.landmarks = 4000
    a.radar_hz = 13.0
    a.baseline = 1
    a.mount = "0.003490659,0,0,3.4,0,0.5"   # ~RADAR_FRONT real (yaw 0.2 deg, lever [3.4,0,0.5])
    a.sweep = None
    a.emit_csv = False
    a.quiet = False
    # FOV + scan-match params (radar_scan_odometry defaults, 3D)
    a.r_min = 0.5
    a.r_max = 100.0
    a.az = math.radians(60.0)
    a.el = math.radians(15.0)
    a.tau_v = 1.0
    a.tau_d = 0.5
    a.rho = 0.5
    a.tau_inlier = 0.5
    a.ransac_iters = 200
    a.min_corr = 4
    a.yaw_margin_deg = 5.0
    a.trans_margin_m = 1.0
    pos = []
    i = 0
    while i < len(argv):
        x = argv[i]
        if x == "--num-points":
            i += 1; a.num_points = int(argv[i])
        elif x == "--sigma-pos":
            i += 1; a.sigma_pos = float(argv[i])
        elif x == "--sigma-vel":
            i += 1; a.sigma_vel = float(argv[i])
        elif x == "--seed":
            i += 1; a.seed = int(argv[i])
        elif x == "--movers":
            i += 1; a.movers = int(argv[i])
        elif x == "--landmarks":
            i += 1; a.landmarks = int(argv[i])
        elif x == "--radar-hz":
            i += 1; a.radar_hz = float(argv[i])
        elif x == "--baseline":
            i += 1; a.baseline = int(argv[i])
        elif x == "--mount":
            i += 1; a.mount = argv[i]
        elif x == "--sweep":
            i += 1; a.sweep = [int(v) for v in argv[i].split(",")]
        elif x == "--r-max":
            i += 1; a.r_max = float(argv[i])
        elif x == "--az-deg":
            i += 1; a.az = math.radians(float(argv[i]))
        elif x == "--el-deg":
            i += 1; a.el = math.radians(float(argv[i]))
        elif x == "--tau-v":
            i += 1; a.tau_v = float(argv[i])
        elif x == "--tau-d":
            i += 1; a.tau_d = float(argv[i])
        elif x == "--rho":
            i += 1; a.rho = float(argv[i])
        elif x == "--tau-inlier":
            i += 1; a.tau_inlier = float(argv[i])
        elif x == "--ransac-iters":
            i += 1; a.ransac_iters = int(argv[i])
        elif x == "--min-corr":
            i += 1; a.min_corr = int(argv[i])
        elif x == "--emit-csv":
            a.emit_csv = True
        elif x == "--quiet":
            a.quiet = True
        else:
            pos.append(x)
        i += 1
    a.pos = pos
    return a


def emit_csv(path, comment, rows, t0_us):
    f = n2c.open_inc(path, comment)
    for (t_us, R, t, var6) in rows:
        n2c.emit_inc(f, (t_us - t0_us) * n2c.US_TO_NS, R, t, var6)
    f.close()


def main():
    args = parse_args(sys.argv[1:])
    if len(args.pos) < 2:
        print("usage: synth_4d_radar.py <gt_csv> <out_dir> [--num-points N] [--sigma-pos M] "
              "[--sigma-vel M_S] [--seed N] [--movers K] [--landmarks M] [--radar-hz HZ] "
              "[--baseline N] [--mount yaw,pitch,roll,lx,ly,lz] [--sweep N0,N1,...] [--emit-csv] [--quiet]")
        sys.exit(2)
    args.gt_csv, out_dir = args.pos[0], args.pos[1]
    os.makedirs(out_dir, exist_ok=True)

    gt_ts, gt_R, gt_p = read_gt(args.gt_csv)
    if len(gt_ts) < 2:
        print("error: GT has < 2 rows")
        sys.exit(1)
    R_x, t_x = mount_from_spec(args.mount)
    args.mount_yaw = math.atan2(R_x[1][0], R_x[0][0])
    args.t_x = t_x
    t0_us = gt_ts[0]

    # the landmark world is FIXED across all N (seeded once) so density is the ONLY variable swept.
    land_rng = np.random.default_rng(args.seed * 2654435761 % (2 ** 32))
    P = make_landmarks(gt_p, args.landmarks, land_rng)

    print_sweep_header(args)
    Ns = args.sweep if args.sweep else [args.num_points]
    reps = []
    stem = os.path.splitext(os.path.basename(args.gt_csv))[0].replace("_gt", "")
    for N in Ns:
        args.num_points = N
        perframe, stats, rows = run_scan_match(gt_ts, gt_R, gt_p, P, R_x, t_x, args, args.seed)
        rep = summarize(perframe, stats, args, "N=%d" % N)
        rep["_N"] = N
        reps.append(rep)
        if args.emit_csv:
            outp = os.path.join(out_dir, "%s_synth4d_n%d.csv" % (stem, N))
            emit_csv(outp, "synth 4D radar scan-odometry (3D R+t), N=%d sigma_pos=%.3f baseline=%d"
                     % (N, args.sigma_pos, args.baseline), rows, t0_us)
            if not args.quiet:
                print("  wrote %s  (%d rows)" % (outp, len(rows)))
    print_sweep_table(reps)


if __name__ == "__main__":
    main()
