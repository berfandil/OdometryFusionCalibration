#!/usr/bin/env python3
# radar_scan_odometry.py - descriptor-based radar scan-matching odometry front-end that recovers
# FULL rigid motion (rotation + translation) from sparse nuScenes radar detections. The principled
# fix flagged at the end of Slice 20b: the Doppler ego-velocity front-end (nuscenes_to_csv.py) gives
# a radar TRANSLATION only (R_B = I, heading-blind), so the radar's rotation extrinsic + hand-eye
# lever stay unobservable. This tool gives the radar a real frame-to-frame R AND t -> an ordinary
# full-3DOF (planar) source, no translation_only flag. Implements RADAR_SCAN_ODOMETRY.md sections 1-3.
#
# THE ALGORITHM (per RADAR_SCAN_ODOMETRY.md sec 1; 2D, point dimension parameterized for a 3D lift):
#   For each consecutive scan pair [prev, curr] of one radar (samples+sweeps union, time-ordered):
#   A. STATIC FILTER - from the reference CAN odometry over [t_prev,t_curr] (base-frame delta A_base)
#      + the radar's calibrated_sensor extrinsic X_radar, form the radar's expected ego-motion
#      B_radar_est = X_radar^-1 o A_base o X_radar, twist (v, omega) = log(B_radar_est)/dt. For each
#      detection i at r_i with measured radial velocity, the static-world prediction is the radial
#      component of -(v + omega x r_i). Keep i iff |measured_radial - predicted_static| < tau_v.
#   B. DESCRIPTOR - per kept detection i: desc(i) = sorted([ ||r_i - r_j|| : j != i kept ]).
#   C. MATCH - for curr i vs prev k: similarity = tolerance-merge COUNT of agreeing distances
#      (|a-b| < tau_d) via a single linear two-pointer walk of the two sorted lists. Accept i<->k iff
#      MUTUAL best match AND count >= rho*min(len_i, len_k).
#   D. RECOVER - 2D Umeyama/Kabsch (rotation + translation, NO scale) on the matched (r_curr, r_prev)
#      pairs, wrapped in seeded RANSAC (minimal sample = 2 pairs in 2D, inliers by residual, refit).
#      Output B_radar = (R 2D, t 2D) over [t_prev,t_curr] WITH a real rotation. Too few corr/inliers
#      -> identity + low-confidence flag (counted).
#
# REUSE, don't reimplement: imports tools/nuscenes_to_csv.py for the radar PCD binary parser, the
# table joins (radar_stream / ego_pose_track), the CAN loading, the SE3/quat/SO3 helpers, the
# calibrated_sensor extrinsic, the increment-CSV emit. Only the scan-matching core is new here.
#
# 3D READINESS (RADAR_SCAN_ODOMETRY.md sec 2): the point dimension is a parameter DIM (2 here). The
# descriptor (B) and matching (C) are dimension-agnostic. Kabsch (D) handles DIM via SVD. The static
# filter (A) uses a planar twist (scalar yaw rate) because the ARS408 is a 2D radar (z=0, vx/vy only,
# verified) -- a 4D/imaging radar drops in with a 3D twist + 3D measured velocity.
#
# Pure ASCII; numpy ok; SEEDED RNG (random.Random) for the RANSAC sampling -> reproducible.
#
# LONGER MATCHING BASELINE (--baseline N, default 1): match scan k to scan k-N instead of k-1, so the
# true accumulated yaw over the matched window is ~N x larger while the Kabsch rotation-noise floor
# stays ~constant -> rotation SNR rises ~N x. Trade-off: larger N -> radar moved more -> less FOV
# overlap -> fewer correspondences. Sweep N to find the sweet-spot before overlap collapse. The
# per-frame static filter is UNCHANGED; only the matched PAIR (and its dt / CAN prior / GT delta /
# plausibility cap) changes. HEADLINE metrics: per-step yaw correlation (recovered vs GT over the same
# interval) + a NON-OVERLAPPING cumulative-heading track (sum over k=N,2N,3N,... windows).
#
# Usage:
#   python tools/radar_scan_odometry.py <dataroot> <out_dir> <scene> [scene ...]
#       [--radar NAME] [--baseline N] [--tau-v M_S] [--tau-d M] [--rho R] [--tau-inlier M]
#       [--ransac-iters N] [--min-corr N] [--seed N] [--quiet]
#   e.g. python tools/radar_scan_odometry.py C:/workspace/data/nuScenes radar_scan_run scene-0061
#        python tools/radar_scan_odometry.py C:/workspace/data/nuScenes run scene-0061 --baseline 8
#
#   --radar NAME (repeatable) restricts to a subset of the 5 radars; default = all 5.
#   --baseline N matches scan k with scan k-N (default 1 = consecutive). Sweep to test radar rotation.
#   Emits <out_dir>/<scene>_scan_<radar>.csv (increment form, the project schema, WITH rotation),
#   a drop-in alternative to the Doppler front-end, plus a per-radar vs-GT validation report.
import math
import os
import sys

import numpy as np

import nuscenes_to_csv as n2c


# ---------- SE(2) / planar helpers (the static-filter twist + Kabsch live in 2D) ----------
def planar_twist_from_se3(R, t, dt):
    # The radar's expected ego-motion B_radar_est = (R 3x3, t 3) over dt. Return the PLANAR
    # instantaneous twist (v=(vx,vy), omega=yaw rate). yaw = atan2(R10,R00); v = log-translation/dt
    # using the SE(2) left-Jacobian inverse so (v,omega) integrate back to (R,t) exactly. For the
    # tiny per-sweep yaw here V^-1 ~ I, but we use the exact form for correctness.
    yaw = math.atan2(R[1][0], R[0][0])
    tx, ty = t[0], t[1]
    if abs(yaw) < 1e-9:
        vx, vy = tx, ty
    else:
        # SE(2) V matrix: t = V * (v_body * 1)  with theta = yaw.
        # V = (1/theta) [ sin th , -(1-cos th) ; (1-cos th), sin th ]
        s, c = math.sin(yaw), math.cos(yaw)
        a = s / yaw
        b = (1.0 - c) / yaw
        det = a * a + b * b
        # V^-1 = (1/det) [ a, b ; -b, a ]
        vx = (a * tx + b * ty) / det
        vy = (-b * tx + a * ty) / det
    return (vx / dt, vy / dt), (yaw / dt)


def static_radial_prediction(r, v, omega):
    # Velocity a STATIC world point at sensor position r would show: v_static = -(v + omega x r).
    # In 2D, omega x r = omega * (-r_y, r_x). Return the (vx,vy) full vector; the radial projection
    # is done by the caller against the detection's own bearing.
    cross = (-omega * r[1], omega * r[0])
    return (-(v[0] + cross[0]), -(v[1] + cross[1]))


# ---------- A. static-object filter (reference-odometry prior, not self-fit) ----------
def static_filter(dets, B_R, B_t, dt, tau_v, min_rho=0.5):
    # dets: list of (x, y, vx, vy) raw detections. B_R/B_t: the radar's expected ego-motion over dt
    # (already mapped into the radar frame by the caller). Keep i iff its measured radial velocity
    # matches the static-world prediction within tau_v. Return list of (x, y) kept positions.
    if dt <= 0:
        return [(x, y) for (x, y, _vx, _vy) in dets]
    (v, omega) = planar_twist_from_se3(B_R, B_t, dt)
    kept = []
    for (x, y, vx, vy) in dets:
        rho = math.hypot(x, y)
        if rho < min_rho:
            continue
        c, s = x / rho, y / rho
        meas_radial = vx * c + vy * s             # radial component of the measured velocity
        vs = static_radial_prediction((x, y), v, omega)
        pred_radial = vs[0] * c + vs[1] * s       # radial component of the static prediction
        if abs(meas_radial - pred_radial) < tau_v:
            kept.append((x, y))
    return kept


# ---------- A (3D). full-3D static-object filter (4D / imaging radar; RADAR_SCAN_ODOMETRY.md sec 2) --
def so3_log(R):
    # Rotation vector (axis*angle, rad) of a 3x3 rotation. R is a plain list-of-lists.
    tr = R[0][0] + R[1][1] + R[2][2]
    c = max(-1.0, min(1.0, (tr - 1.0) / 2.0))
    th = math.acos(c)
    if th < 1e-9:
        # first-order: w = [R21-R12, R02-R20, R10-R01]/2
        return [(R[2][1] - R[1][2]) * 0.5, (R[0][2] - R[2][0]) * 0.5, (R[1][0] - R[0][1]) * 0.5]
    k = th / (2.0 * math.sin(th))
    return [(R[2][1] - R[1][2]) * k, (R[0][2] - R[2][0]) * k, (R[1][0] - R[0][1]) * k]


def twist_from_se3_3d(R, t, dt):
    # Full-3D instantaneous twist (v in R^3, omega in R^3) from the radar's expected ego-motion
    # (R 3x3, t 3) over dt. omega = log(R)/dt; for the tiny per-scan rotation the SE(3) left-Jacobian
    # V^-1 ~ I, so v = t/dt is accurate to second order. Returns ((vx,vy,vz),(wx,wy,wz)).
    w = so3_log(R)
    return (t[0] / dt, t[1] / dt, t[2] / dt), (w[0] / dt, w[1] / dt, w[2] / dt)


def static_filter_3d(dets, B_R, B_t, dt, tau_v, min_rho=0.5):
    # The 3D analogue of static_filter for a 4D/imaging radar. dets: list of (x,y,z, vx,vy,vz) raw
    # detections (3D position + 3D measured velocity). The velocity a STATIC world point at sensor
    # position r would show is v_static = -(v + omega x r) with r,v,omega in R^3 (sec 2). Keep i iff
    # the radial component of the measured velocity matches the radial component of v_static within
    # tau_v. Return list of (x,y,z) kept positions.
    if dt <= 0:
        return [(x, y, z) for (x, y, z, _vx, _vy, _vz) in dets]
    (v, omega) = twist_from_se3_3d(B_R, B_t, dt)
    kept = []
    for (x, y, z, vx, vy, vz) in dets:
        rho = math.sqrt(x * x + y * y + z * z)
        if rho < min_rho:
            continue
        ux, uy, uz = x / rho, y / rho, z / rho                 # unit bearing (3D)
        meas_radial = vx * ux + vy * uy + vz * uz
        # omega x r
        cx = omega[1] * z - omega[2] * y
        cy = omega[2] * x - omega[0] * z
        cz = omega[0] * y - omega[1] * x
        vs = (-(v[0] + cx), -(v[1] + cy), -(v[2] + cz))        # v_static = -(v + omega x r)
        pred_radial = vs[0] * ux + vs[1] * uy + vs[2] * uz
        if abs(meas_radial - pred_radial) < tau_v:
            kept.append((x, y, z))
    return kept


# ---------- B. per-detection descriptor (sorted inter-point distance multiset) ----------
def build_descriptors(pts):
    # pts: list of (x, y) (DIM-agnostic; here DIM=2). Return a list of sorted numpy distance arrays,
    # one per point: desc(i) = sort({ ||r_i - r_j|| : j != i }). Vectorized pairwise distances.
    n = len(pts)
    if n < 2:
        return [np.zeros(0) for _ in range(n)]
    P = np.asarray(pts, dtype=float)              # (n, DIM)
    # full pairwise distance matrix
    diff = P[:, None, :] - P[None, :, :]
    D = np.sqrt(np.sum(diff * diff, axis=2))      # (n, n), zero diagonal
    descs = []
    for i in range(n):
        row = np.delete(D[i], i)                  # drop self-distance
        row.sort()
        descs.append(row)
    return descs


# ---------- C. cross-frame matching (linear tolerance-merge of two sorted lists) ----------
def merge_count(a, b, tau_d):
    # a, b: sorted numpy arrays. Count agreeing pairs (|a_i - b_j| < tau_d) by a single linear
    # two-pointer walk: a matched pair advances BOTH pointers (each distance used once); otherwise
    # advance the smaller. This is the partial-overlap-robust similarity score.
    i = j = 0
    na, nb = len(a), len(b)
    count = 0
    while i < na and j < nb:
        d = a[i] - b[j]
        if abs(d) < tau_d:
            count += 1
            i += 1
            j += 1
        elif d < 0:
            i += 1
        else:
            j += 1
    return count


def match_frames(descs_curr, descs_prev, tau_d, rho):
    # Build the agreeing-count matrix, then accept i<->k as a correspondence iff it is a MUTUAL best
    # match AND count >= rho*min(len_i,len_k). Returns list of (i_curr, k_prev) index pairs.
    nc, npv = len(descs_curr), len(descs_prev)
    if nc == 0 or npv == 0:
        return []
    S = np.zeros((nc, npv), dtype=int)
    for i in range(nc):
        di = descs_curr[i]
        for k in range(npv):
            S[i, k] = merge_count(di, descs_prev[k], tau_d)
    # best previous for each current, best current for each previous
    best_prev = np.argmax(S, axis=1)              # for curr i -> prev k
    best_curr = np.argmax(S, axis=0)              # for prev k -> curr i
    pairs = []
    for i in range(nc):
        k = int(best_prev[i])
        if int(best_curr[k]) != i:                # mutual-best test
            continue
        floor = rho * min(len(descs_curr[i]), len(descs_prev[k]))
        if S[i, k] >= floor and S[i, k] > 0:
            pairs.append((i, k))
    return pairs


# ---------- D. motion recovery (2D/3D Kabsch, no scale) wrapped in seeded RANSAC ----------
def kabsch(src, dst):
    # Solve dst ~= R @ src + t (rotation + translation, NO scale) by Umeyama/Kabsch. src/dst are
    # (m, DIM) arrays. Returns (R DIMxDIM, t DIM). DIM-agnostic (SVD on the DIMxDIM covariance).
    cs = src.mean(axis=0)
    cd = dst.mean(axis=0)
    S = src - cs
    Dd = dst - cd
    H = S.T @ Dd                                   # (DIM, DIM)
    U, _sv, Vt = np.linalg.svd(H)
    Dsign = np.eye(H.shape[0])
    Dsign[-1, -1] = np.sign(np.linalg.det(Vt.T @ U.T))  # reflection fix -> proper rotation
    R = Vt.T @ Dsign @ U.T
    t = cd - R @ cs
    return R, t


def ransac_kabsch(curr_pts, prev_pts, pairs, rng, tau_inlier, iters, min_corr, dim=2):
    # curr_pts/prev_pts: kept point lists for the two frames. pairs: (i_curr,k_prev) correspondences.
    # Recover R,t with r_prev ~= R @ r_curr + t (the radar increment maps CURR -> PREV; the increment
    # CSV contract is "motion from prev to curr" expressed as the body delta -- see emit note below).
    # Seeded RANSAC: sample DIM pairs, fit Kabsch, score inliers by ||r_prev - (R r_curr + t)|| <
    # tau_inlier, keep the best, refit on its inliers. Returns (R, t, n_inliers, n_corr).
    m = len(pairs)
    if m < max(min_corr, dim + 1):
        return None, None, 0, m
    src = np.array([curr_pts[i] for (i, _k) in pairs], dtype=float)   # CURR
    dst = np.array([prev_pts[k] for (_i, k) in pairs], dtype=float)   # PREV
    best_inl = None
    sample_n = dim                                  # minimal sample: 2 in 2D, 3 in 3D
    for _ in range(iters):
        idx = rng.sample(range(m), sample_n)
        try:
            R, t = kabsch(src[idx], dst[idx])
        except np.linalg.LinAlgError:
            continue
        pred = (R @ src.T).T + t
        res = np.linalg.norm(dst - pred, axis=1)
        inl = np.where(res < tau_inlier)[0]
        if best_inl is None or len(inl) > len(best_inl):
            best_inl = inl
    if best_inl is None or len(best_inl) < max(min_corr, dim + 1):
        return None, None, 0 if best_inl is None else len(best_inl), m
    R, t = kabsch(src[best_inl], dst[best_inl])     # final refit on inliers
    return R, t, len(best_inl), m


# ---------- per-radar scan-matching odometry over a scene ----------
def base_delta_from_can(vm_sorted, vm_ts, t_prev_us, t_curr_us, speed_scale, yaw_scale):
    # Reference CAN base-frame delta A_base over [t_prev, t_curr] (ego frame). The radar sweep window
    # (~75 ms) is SHORTER than the CAN vehicle_monitor period (~500 ms), so integrating the discrete
    # samples that fall strictly INSIDE the window usually integrates ZERO samples -> a spurious
    # identity delta -> zero predicted ego velocity -> the static filter rejects every moving-scene
    # detection. Instead evaluate the INSTANTANEOUS CAN twist at the window MIDPOINT (nearest-in-time
    # speed + yaw_rate) and hold it constant over dt. dt is small so a constant-twist delta is exact
    # enough for the static gate. Returns (R 3x3, t 3).
    import bisect
    dt = (t_curr_us - t_prev_us) * 1e-6
    if dt <= 0:
        return [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0]
    t_mid = (t_prev_us + t_curr_us) // 2
    k = min(max(bisect.bisect_right(vm_ts, t_mid) - 1, 0), len(vm_sorted) - 1)
    v = vm_sorted[k]["vehicle_speed"] * speed_scale
    w = vm_sorted[k]["yaw_rate"] * yaw_scale
    R = n2c.rot_z(w * dt)
    t = [v * dt, 0.0, 0.0]
    return R, t


def expected_radar_motion(A_R, A_t, X_R, X_t):
    # B_radar_est = X_radar^-1 o A_base o X_radar  (the base ego-motion seen in the radar frame).
    # X = (X_R, X_t) is ego_from_radar (calibrated_sensor). Compose SE(3): apply X, then A, then X^-1.
    XRt = n2c.mat_t(X_R)
    # B_R = X_R^T A_R X_R
    BR = n2c.mat_mul(XRt, n2c.mat_mul(A_R, X_R))
    # B_t = X_R^T ( A_R X_t + A_t - X_t )
    inner = [n2c.mat_vec(A_R, X_t)[i] + A_t[i] - X_t[i] for i in range(3)]
    Bt = n2c.mat_vec(XRt, inner)
    return BR, Bt


def wheel_unit_scales(vm, gt_R):
    # Mirror nuscenes_to_csv's unit self-check (km/h vs m/s, deg/s vs rad/s) so the static filter's
    # reference twist is in SI. Pick the combo whose dead-reckoned net yaw + distance best matches GT.
    vm_s = sorted(vm, key=lambda r: r["utime"])
    speed_choices = [("km/h", 1.0 / 3.6), ("m/s", 1.0)]
    yaw_choices = [("deg/s", math.pi / 180.0), ("rad/s", 1.0)]
    gt_net_yaw = math.atan2(gt_R[-1][1][0], gt_R[-1][0][0])
    best = None
    for (_sn, sv) in speed_choices:
        for (_yn, yv) in yaw_choices:
            est_yaw = 0.0
            L = 0.0
            for k in range(1, len(vm_s)):
                dt = (vm_s[k]["utime"] - vm_s[k - 1]["utime"]) * 1e-6
                if dt <= 0:
                    continue
                est_yaw += vm_s[k]["yaw_rate"] * yv * dt
                L += vm_s[k]["vehicle_speed"] * sv * dt
            yaw_err = abs(((est_yaw - gt_net_yaw + math.pi) % (2 * math.pi)) - math.pi)
            score = yaw_err + 0.001 * abs(L)
            if best is None or score < best[0]:
                best = (score, sv, yv)
    return best[1], best[2]


def radar_scan_odometry(stream, calib, vm, blobroot, gt_ts, gt_R, gt_p,
                        speed_scale, yaw_scale, rng, args):
    # The full per-radar pipeline. stream = [(ts_us, filename, ego_pose), ...] time-ordered. calib =
    # calibrated_sensor (ego_from_radar). Returns (rows, stats, perframe) where rows are increment-CSV
    # tuples (t_us, R3, t3, var6) and perframe holds the per-pair record (R_est, t_est, R_gt, t_gt,
    # recovered, ts_curr_us, ts_part_us, k_idx) for the vs-GT validation.
    #
    # LONGER MATCHING BASELINE (--baseline N): the scan at ordered index k is matched with the scan
    # at index k-N (default N=1 = consecutive, the original behaviour). Both frames are STILL static-
    # filtered the same way per-frame (the per-frame filter is UNCHANGED); only the matched PAIR
    # changes. The recovered B_radar then spans [t_{k-N}, t_k] (~N x the inter-scan time), so the true
    # accumulated yaw is ~N x larger while the Kabsch rotation-noise floor stays ~constant -> the
    # rotation SNR should rise. dt, the CAN prior, the plausibility cap and the GT delta ALL span the
    # longer [t_{k-N}, t_k] interval. Trade-off: larger N -> the radar moved more -> less FOV overlap
    # -> fewer correspondences -> matching degrades, so there is a sweet-spot N.
    import bisect
    vm_sorted = sorted(vm, key=lambda r: r["utime"])
    vm_ts = [r["utime"] for r in vm_sorted]
    # radar extrinsic ego_from_radar (X_radar)
    X_R = n2c.quat_wxyz_to_mat(*calib["rotation"])
    X_t = list(calib["translation"])
    X_Rt = n2c.mat_t(X_R)

    N = max(1, int(args.baseline))

    rows = []
    perframe = []
    stats = {"n_pairs": 0, "n_fallback": 0, "sum_det": 0, "sum_static": 0,
             "sum_corr": 0, "sum_inl": 0, "n_recov": 0, "n_gated": 0}

    # ring buffer of the last N+1 scans' RAW detections keyed by ordered index k. The partner for the
    # scan at index k is the one at index k-N (so we only need to remember the last N+1 frames).
    hist = []  # list of (ts_us, dets), most-recent last; trimmed to length N+1
    for k, (ts_us, filename, _ego) in enumerate(stream):
        dets = n2c.read_radar_pcd(os.path.join(blobroot, filename))
        hist.append((ts_us, dets))
        if len(hist) > N + 1:
            hist.pop(0)
        if k < N:
            # not enough history yet to reach back N scans (the first N scans are skipped, harmless)
            continue

        t_prev_us, prev_dets = hist[0]   # the scan at index k-N
        dt = (ts_us - t_prev_us) * 1e-6
        if dt <= 0:
            continue

        # --- reference CAN base delta over [t_{k-N}, t_k] -> expected radar motion ---
        A_R, A_t = base_delta_from_can(vm_sorted, vm_ts, t_prev_us, ts_us, speed_scale, yaw_scale)
        B_R, B_t = expected_radar_motion(A_R, A_t, X_R, X_t)

        # --- A. static filter BOTH frames with the same expected-motion prior (per-frame, UNCHANGED) ---
        prev_kept = static_filter(prev_dets, B_R, B_t, dt, args.tau_v)
        curr_kept = static_filter(dets, B_R, B_t, dt, args.tau_v)
        stats["sum_det"] += len(dets)
        stats["sum_static"] += len(curr_kept)
        stats["n_pairs"] += 1

        # --- B. descriptors ---
        descs_prev = build_descriptors(prev_kept)
        descs_curr = build_descriptors(curr_kept)

        # --- C. match ---
        pairs = match_frames(descs_curr, descs_prev, args.tau_d, args.rho)
        stats["sum_corr"] += len(pairs)

        # --- D. recover ---
        R2, t2, n_inl, n_corr = ransac_kabsch(curr_kept, prev_kept, pairs, rng,
                                              args.tau_inlier, args.ransac_iters,
                                              args.min_corr, dim=2)
        stats["sum_inl"] += n_inl

        # PHYSICAL-PLAUSIBILITY GATE: sparse radar descriptors occasionally produce a confident but
        # GROSSLY wrong Kabsch fit (tens of degrees of yaw) that survives RANSAC when the geometry is
        # near-degenerate. The median frame is fine but these blunders destroy the integrated heading.
        # Reject any recovered increment whose yaw or translation exceeds the CAN-predicted radar
        # motion over the (longer) [t_{k-N}, t_k] interval by a generous margin -> treat as fallback.
        # This uses the SAME reference prior the static filter already trusts (not GT).
        recovered = R2 is not None
        if recovered:
            yaw_inc = abs(math.atan2(R2[1, 0], R2[0, 0]))   # |yaw| of the recovered map (sign-free)
            (vp, omp) = planar_twist_from_se3(B_R, B_t, dt)
            yaw_pred = abs(omp * dt)
            spd_pred = math.hypot(vp[0], vp[1]) * dt
            yaw_cap = yaw_pred + math.radians(args.yaw_margin_deg)   # prior yaw + slack
            t_cap = spd_pred + args.trans_margin_m                   # prior step + slack
            t_mag = math.hypot(float(t2[0]), float(t2[1]))
            if yaw_inc > yaw_cap or t_mag > t_cap:
                stats["n_gated"] += 1
                recovered = False

        if not recovered:
            # fallback: identity increment, low confidence (huge covariance)
            R3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            t3 = [0.0, 0.0, 0.0]
            var6 = [100.0, 100.0, 100.0, (0.5) ** 2, (0.5) ** 2, (0.5) ** 2]
            stats["n_fallback"] += 1
            R_inc, t_inc = R3, t3
        else:
            # Kabsch gave r_prev ~= R2 r_curr + t2 = M_prev^-1 M_curr (M = radar pose), which IS the
            # forward body delta over [t_{k-N}, t_k] (pose_curr = pose_prev o delta) -> the increment
            # is the Kabsch map DIRECTLY: R_inc = R2, t_inc = t2 (planar; z row = identity). The earlier
            # R2^T / -R2^T t2 was a yaw-SIGN FLIP (invisible on noise-dominated real 2D radar where the
            # rotation is unusable either way, but a clean -97 vs +98 deg inversion on clean dense data
            # -- caught + verified by the synth_4d_radar planted-yaw round-trip vs B_radar_gt).
            R_inc = [[float(R2[0, 0]), float(R2[0, 1]), 0.0],
                     [float(R2[1, 0]), float(R2[1, 1]), 0.0],
                     [0.0, 0.0, 1.0]]
            t_inc = [float(t2[0]), float(t2[1]), 0.0]
            stats["n_recov"] += 1
            # confidence variances: scale rot var by inlier count (more inliers -> tighter)
            s_t = (0.05) ** 2                          # ~5 cm/scan translation floor
            s_r = max(0.02, 0.3 / max(1, n_inl)) ** 2  # rot tightens with inliers
            var6 = [s_t, s_t, 100.0, s_r, s_r, s_r]

        rows.append((ts_us, R_inc, t_inc, var6))

        # --- validation: GT radar frame-to-frame motion over [t_{k-N}, t_k] ---
        # B_radar_gt = X_radar^-1 o A_gt o X_radar  (forward body delta partner->curr in radar frame)
        gi0 = min(max(bisect.bisect_right(gt_ts, t_prev_us) - 1, 0), len(gt_R) - 1)
        gi1 = min(max(bisect.bisect_right(gt_ts, ts_us) - 1, 0), len(gt_R) - 1)
        R0g, p0g = gt_R[gi0], gt_p[gi0]
        R1g, p1g = gt_R[gi1], gt_p[gi1]
        R0gt = n2c.mat_t(R0g)
        Agt_R = n2c.mat_mul(R0gt, R1g)
        Agt_t = n2c.mat_vec(R0gt, [p1g[i] - p0g[i] for i in range(3)])
        Bgt_R, Bgt_t = expected_radar_motion(Agt_R, Agt_t, X_R, X_t)
        perframe.append((R_inc, t_inc, Bgt_R, Bgt_t, recovered, ts_us, t_prev_us, k))

    return rows, stats, perframe


# ---------- validation metrics ----------
def rot_err_deg(Re, Rg):
    # geodesic rotation error between two 3x3 rotations, in degrees.
    Rt = n2c.mat_mul(n2c.mat_t(Re), Rg)
    tr = Rt[0][0] + Rt[1][1] + Rt[2][2]
    c = max(-1.0, min(1.0, (tr - 1.0) / 2.0))
    return math.degrees(math.acos(c))


def yaw_of(R):
    return math.atan2(R[1][0], R[0][0])


def pct(vals, q):
    if not vals:
        return float("nan")
    s = sorted(vals)
    k = min(len(s) - 1, int(round(q * (len(s) - 1))))
    return s[k]


def pearson(xs, ys):
    # Pearson correlation of two equal-length sequences; nan if degenerate (no variance / <2 points).
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = math.sqrt(sxx * syy)
    if den < 1e-12:
        return float("nan")
    return sxy / den


def report_radar(channel, perframe, stats, baseline):
    # Per-radar: translation error (m) + rotation error (deg), recovered vs GT, median + p90, over
    # frames where a real motion was recovered (fallback frames excluded from the error stats but
    # counted in the fallback rate). Plus the HEADLINE rotation metrics: the per-step yaw correlation
    # (recovered vs GT yaw of B_radar over the SAME [t_{k-N}, t_k] interval -- was ~0 at N=1; does it
    # rise toward 1 as N grows?) and a NON-OVERLAPPING cumulative-heading track.
    N = max(1, int(baseline))
    trans_errs, rot_errs = [], []
    # per-step yaw arrays for correlation (recovered-motion frames only -- fallback yaw is a forced 0
    # and would spuriously inflate/deflate the correlation; correlate where we actually claimed a R).
    est_yaws, gt_yaws = [], []
    # cumulative heading over ALL frames (fallback contributes 0 recovered yaw, like the Doppler R=I)
    cum_est_yaw = 0.0
    cum_gt_yaw = 0.0
    # NON-OVERLAPPING cumulative heading: sum the recovered (and GT) yaw over the pairs whose curr
    # ordered-index k is a multiple of N (k=N,2N,3N,...). Those windows tile [t_0,t_N,t_2N,...]
    # end-to-end with NO overlap, so summing them is a true heading track vs the GT heading at the
    # same break times. Fallback frames keep their 0 recovered yaw (honest: no heading recovered).
    nov_est_yaw = 0.0
    nov_gt_yaw = 0.0
    n_nov = 0
    for rec in perframe:
        (R_inc, t_inc, Rg, tg, recovered, ts_us, t_prev_us, k) = rec
        cum_est_yaw += yaw_of(R_inc)
        cum_gt_yaw += yaw_of(Rg)
        if k % N == 0:
            nov_est_yaw += yaw_of(R_inc)
            nov_gt_yaw += yaw_of(Rg)
            n_nov += 1
        if not recovered:
            continue
        te = math.hypot(t_inc[0] - tg[0], t_inc[1] - tg[1])
        re = rot_err_deg(R_inc, Rg)
        trans_errs.append(te)
        rot_errs.append(re)
        est_yaws.append(math.degrees(yaw_of(R_inc)))
        gt_yaws.append(math.degrees(yaw_of(Rg)))
    n = stats["n_pairs"]
    nov_est_deg = math.degrees(nov_est_yaw)
    nov_gt_deg = math.degrees(nov_gt_yaw)
    out = {
        "channel": channel,
        "baseline": N,
        "n_pairs": n,
        "n_recov": stats["n_recov"],
        "fallback_pct": 100.0 * stats["n_fallback"] / n if n else 0.0,
        "gated_pct": 100.0 * stats["n_gated"] / n if n else 0.0,
        "avg_det": stats["sum_det"] / n if n else 0.0,
        "avg_static": stats["sum_static"] / n if n else 0.0,
        "avg_corr": stats["sum_corr"] / n if n else 0.0,
        "avg_inl": stats["sum_inl"] / n if n else 0.0,
        "trans_med": pct(trans_errs, 0.5),
        "trans_p90": pct(trans_errs, 0.9),
        "rot_med": pct(rot_errs, 0.5),
        "rot_p90": pct(rot_errs, 0.9),
        "yaw_corr": pearson(est_yaws, gt_yaws),   # HEADLINE: per-step recovered-vs-GT yaw correlation
        "cum_est_yaw_deg": math.degrees(cum_est_yaw),
        "cum_gt_yaw_deg": math.degrees(cum_gt_yaw),
        # non-overlapping cumulative-heading track + final error vs GT net yaw over those windows
        "nov_est_yaw_deg": nov_est_deg,
        "nov_gt_yaw_deg": nov_gt_deg,
        "nov_heading_err_deg": nov_est_deg - nov_gt_deg,
        "n_nov": n_nov,
    }
    return out


# ---------- main ----------
class Args:
    pass


def parse_args(argv):
    a = Args()
    a.tau_v = 1.0         # m/s static-velocity gate (radar Doppler is noisy; ~35% kept on mini)
    a.tau_d = 0.5         # m descriptor distance-agreement tolerance
    a.rho = 0.5           # min agreeing-fraction floor for a correspondence
    a.tau_inlier = 0.5    # m RANSAC inlier residual threshold
    a.ransac_iters = 100
    a.min_corr = 3        # min correspondences/inliers to recover a real motion
    a.yaw_margin_deg = 5.0   # plausibility gate: max |recovered yaw| above the CAN-predicted yaw
    a.trans_margin_m = 1.0   # plausibility gate: max |recovered trans| above the CAN-predicted step
    a.baseline = 1        # match scan k with scan k-N (N=1 = consecutive, the original behaviour)
    a.seed = 12345
    a.radars = None       # None -> all 5
    a.quiet = False
    pos = []
    i = 0
    while i < len(argv):
        x = argv[i]
        if x == "--tau-v":
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
        elif x == "--yaw-margin-deg":
            i += 1; a.yaw_margin_deg = float(argv[i])
        elif x == "--trans-margin-m":
            i += 1; a.trans_margin_m = float(argv[i])
        elif x == "--baseline":
            i += 1; a.baseline = int(argv[i])
        elif x == "--seed":
            i += 1; a.seed = int(argv[i])
        elif x == "--radar":
            i += 1
            a.radars = (a.radars or []) + [argv[i].upper()]
        elif x == "--quiet":
            a.quiet = True
        else:
            pos.append(x)
        i += 1
    a.pos = pos
    return a


def emit_scan_csv(path, comment, rows, t0_us):
    f = n2c.open_inc(path, comment)
    for (t_us, R, t, var6) in rows:
        n2c.emit_inc(f, (t_us - t0_us) * n2c.US_TO_NS, R, t, var6)
    f.close()


def convert_scene(dataroot, out_dir, scene_name, tables, args):
    import random
    scene, sample_tokens = n2c.scene_sample_tokens(tables, scene_name)
    gt = n2c.ego_pose_track(tables, sample_tokens)
    if len(gt) < 2:
        raise ValueError("scene %s has < 2 ego_poses" % scene_name)
    t0_us = gt[0]["timestamp"]

    # GT track in the GLOBAL (sensor-absolute) frame for the radar-GT delta (we compose deltas, so
    # using global poses directly is fine; the X_radar^-1 A X cancels the global frame).
    gt_ts = [p["timestamp"] for p in gt]
    gt_R = [n2c.quat_wxyz_to_mat(*p["rotation"]) for p in gt]
    gt_p = [list(p["translation"]) for p in gt]

    vm = n2c.load_can(dataroot, scene_name, "vehicle_monitor")
    if not vm:
        raise ValueError("no vehicle_monitor CAN for %s" % scene_name)
    speed_scale, yaw_scale = wheel_unit_scales(vm, gt_R)

    blobroot = n2c.blob_dir(dataroot)
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, scene_name)

    channels = n2c.RADAR_CHANNELS if args.radars is None else \
        [c for c in n2c.RADAR_CHANNELS if c in args.radars]

    reports = []
    for ri, channel in enumerate(n2c.RADAR_CHANNELS):
        if channel not in channels:
            continue
        calib, stream = n2c.radar_stream(tables, sample_tokens, channel)
        if not stream:
            if not args.quiet:
                print("  [warn] %s: no sweeps in scene; skipping" % channel)
            continue
        sid = n2c.RADAR_SRC_BASE + ri
        rng = random.Random(args.seed * 1000003 + sid)
        rows, stats, perframe = radar_scan_odometry(
            stream, calib, vm, blobroot, gt_ts, gt_R, gt_p,
            speed_scale, yaw_scale, rng, args)
        short = channel.lower()
        # tag the CSV with the baseline so sweep runs (different N) do not overwrite each other
        nbase = max(1, int(args.baseline))
        nsuf = "" if nbase == 1 else ("_n%d" % nbase)
        outp = base + "_scan_" + short.replace("radar_", "radar_") + nsuf + ".csv"
        emit_scan_csv(outp, "%s scan-odometry: descriptor scan-match (R+t), baseline N=%d"
                      % (channel, nbase), rows, t0_us)
        rep = report_radar(channel, perframe, stats, args.baseline)
        reports.append(rep)
    return reports, scene_name


def print_report(scene_name, reports, args):
    nbase = max(1, int(args.baseline))
    print("=" * 110)
    print("RADAR SCAN-MATCHING ODOMETRY -- %s   N=%d  (tau_v=%.2f tau_d=%.2f rho=%.2f tau_inl=%.2f "
          "iters=%d)" % (scene_name, nbase, args.tau_v, args.tau_d, args.rho, args.tau_inlier,
                         args.ransac_iters))
    print("=" * 110)
    # match-quality block
    print("-" * 110)
    print("MATCH QUALITY (per frame averages)  -- overlap collapse shows here as corr/inlier -> 0, fallbk -> 100%")
    print("%-20s | %5s | %7s | %7s | %7s | %7s | %8s | %7s"
          % ("radar", "#pair", "det", "static", "corr", "inlier", "fallbk%", "gated%"))
    print("-" * 90)
    for r in reports:
        print("%-20s | %5d | %7.1f | %7.1f | %7.2f | %7.2f | %7.1f%% | %6.1f%%"
              % (r["channel"], r["n_pairs"], r["avg_det"], r["avg_static"],
                 r["avg_corr"], r["avg_inl"], r["fallback_pct"], r["gated_pct"]))
    # accuracy block
    print("-" * 110)
    print("RECOVERED vs GT  (over frames with a real recovered motion; errors are over the [t_{k-N}, t_k] interval)")
    print("%-20s | %12s | %12s | %12s | %12s"
          % ("radar", "trans med m", "trans p90 m", "rot med deg", "rot p90 deg"))
    print("-" * 78)
    for r in reports:
        print("%-20s | %12.3f | %12.3f | %12.3f | %12.3f"
              % (r["channel"], r["trans_med"], r["trans_p90"], r["rot_med"], r["rot_p90"]))
    # HEADLINE rotation block: per-step yaw correlation + non-overlapping cumulative heading track
    print("-" * 110)
    print("ROTATION HEADLINE  (per-step yaw CORRELATION recovered-vs-GT; was ~0 at N=1 -> rises toward 1?)")
    print("%-20s | %9s | %16s | %16s | %16s"
          % ("radar", "yaw corr", "nonovl recov hdg", "nonovl GT hdg", "nonovl hdg err"))
    print("-" * 88)
    for r in reports:
        print("%-20s | %9.3f | %12.2f deg | %12.2f deg | %12.2f deg"
              % (r["channel"], r["yaw_corr"], r["nov_est_yaw_deg"], r["nov_gt_yaw_deg"],
                 r["nov_heading_err_deg"]))
    print("=" * 110)


def main():
    args = parse_args(sys.argv[1:])
    if len(args.pos) < 3:
        print("usage: radar_scan_odometry.py <dataroot> <out_dir> <scene> [scene ...] "
              "[--radar NAME] [--baseline N] [--tau-v M_S] [--tau-d M] [--rho R] [--tau-inlier M] "
              "[--ransac-iters N] [--min-corr N] [--yaw-margin-deg D] [--trans-margin-m M] "
              "[--seed N] [--quiet]")
        sys.exit(2)
    dataroot, out_dir = args.pos[0], args.pos[1]
    scenes = args.pos[2:]
    print("loading nuScenes tables...")
    tables = n2c.load_tables(dataroot)
    for scene_name in scenes:
        try:
            reports, sc = convert_scene(dataroot, out_dir, scene_name, tables, args)
        except Exception as e:  # noqa: BLE001
            import traceback
            print("ERR %s: %s" % (scene_name, e))
            traceback.print_exc()
            continue
        print_report(sc, reports, args)


if __name__ == "__main__":
    main()
