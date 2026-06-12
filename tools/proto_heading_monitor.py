# Prototype (Slice 19c): GPS-course heading-drift monitor.
# Hypothesis: GPS course-over-ground can RANK the heading quality of N odometry
# sources online (auto-discover the heading-grade source, e.g. the FOG) without
# config. Weighting needs only a RANKING, not an absolute bias estimate, which
# sidesteps the Slice-18/D28 observability boundary: position-only GPS cannot
# observe per-source yaw bias absolutely, but source-vs-course drift RATES are
# directly comparable across sources.
#
# Pipeline (all online-implementable: running accumulators, one block buffer,
# one slope-sample pool per source):
#   1. GPS lat/lon -> local EN (equirectangular about the first fix). Course
#      sample k = atan2(dN, dE) over fix interval k..k+1 (ENU/CCW-from-east so
#      delta-course matches delta-yaw sign), stamped at the interval midpoint.
#      ANCHOR validity gates (course samples are only trusted as pair anchors
#      when the chord direction is a clean heading proxy):
#        - fix spacing dt in (0.05, dtmax]
#        - vmin <= v_gps <= vmax            (course is meaningless when slow)
#        - |v_gps - v_odom| < max(3, 0.5 v) (cross-validation vs the sources'
#          own translation increments; kills multipath position jumps)
#        - net odometry forward motion > 0  (reverse flips course by 180 deg)
#        - |yaw rate| < omega_max over the fix interval (during a turn the
#          chord direction != heading; sub-fix-rate turn transients produced
#          -30..-45 deg chord errors on KAIST)
#   2. Pairs of CONSECUTIVE valid anchors with midpoint gap <= pairgap. Long
#      pairs deliberately BRIDGE low-speed turn stretches: wheel-yaw drift is
#      turn-correlated, so gating turns out of the pair INTERIOR would hide
#      exactly the junk source. Per source i: r_i = wrap(dcourse - dyaw_i),
#      dyaw from composed increment quaternions at the anchor midtimes. A
#      fixed mounting/frame yaw offset cancels in the deltas (no
#      odom_from_enu_yaw needed). NOTHING is hard-dropped after this point:
#      every hard drop of a pair breaks the telescoping of the cumulative sum
#      and injects the dropped chunk as a permanent offset.
#   3. Cross-source split: course error is COMMON to all sources, so
#      m = median_i(r_i) estimates it (on KAIST it matches across sources to
#      <0.1 deg). dev_i = r_i - m is the per-source channel:
#        rate channel  rc_i = m + clamp(dev_i, +-event)  -> cumR_i = sum(rc_i)
#        event channel P_i += |dev_i| - event excess (wheel slips / encoder
#        glitches appear as a step in ONE source only; a magnitude gate would
#        reject exactly the bridge pair that catches the step, the clamp+P
#        split keeps it and attributes it to the right source)
#   4. Drift-rate estimate: cumR_i(t) telescopes to (course error at current
#      anchor) - (course error at segment-start anchor) - yaw drift, so its
#      SLOPE is the drift rate and anchor noise enters only at the evaluation
#      points. Estimator: split each GPS-contiguous SEGMENT (no anchor gap >
#      pairgap) into blocks of block_s seconds, take block medians of
#      (t, cumR), pool pairwise block slopes with baseline >= min_base WITHIN
#      segments (never across: a GPS-denied / multipath stretch steps the
#      staircase between segments), estimate = baseline-weighted median of the
#      pool (Theil-Sen flavored; an isolated multipath block is outvoted).
#   5. Score_i (deg/h) = |slope_i| + P_i / T_valid. Ranking by score;
#      discovery latency = earliest time after which the score ranking equals
#      the final ranking at every later evaluation.
#   6. Candidate weight rule: rot_weight_i = clip(wcap * min_j score_j /
#      score_i, 1, wcap) -- the winner gets the manual prior (wcap = 10).
#
# Measured floor on KAIST (VRS-RTK, 1 Hz): the absolute drift estimate is good
# to ~5-15 deg/h; sources below that floor are mutually indistinguishable.
# With < 3 sources the median split degrades (m = 0 fallback, no event
# channel consensus) -- the monitor is designed for N >= 3.
#
# Falsification probes: --deny T0 T1 masks GPS fixes (the monitor must freeze,
# not emit garbage); --sweep grids vmin x omega_max x pairgap and reports
# whether any reasonable setting flips the final ranking; --gt scores the
# ranking against the ground-truth reference (integrated source yaw vs GT yaw).
#
# Usage:
#   proto_heading_monitor.py GPS.csv SRC.csv [SRC2.csv ...] [--gt GT.csv]
#       [--vmin 3] [--vmax 30] [--dtmax 3] [--omega 3] [--pairgap 60]
#       [--event 5] [--block 60] [--minbase 120] [--win 600]
#       [--deny T0 T1] [--sweep] [--series OUT_PREFIX]
import sys, csv, math
from collections import deque
import numpy as np

R_EARTH = 6378137.0
SIGMA_POS = 0.5  # m, 1-sigma horizontal (VRS-RTK var 0.25 m^2 in the CSVs)


def wrap_pi(a):
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def load_gps(path):
    t, lat, lon = [], [], []
    with open(path) as f:
        for r in csv.reader(f):
            if not r or r[0].startswith('#'):
                continue
            t.append(int(r[0]) * 1e-9)
            lat.append(float(r[1]))
            lon.append(float(r[2]))
    t = np.asarray(t)
    lat = np.radians(np.asarray(lat))
    lon = np.radians(np.asarray(lon))
    E = (lon - lon[0]) * math.cos(lat[0]) * R_EARTH
    N = (lat - lat[0]) * R_EARTH
    return t, E, N


def quat_mul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw)


def quat_to_yaw(qw, qx, qy, qz):
    return math.atan2(2 * (qw * qz + qx * qy), 1 - 2 * (qy * qy + qz * qz))


def load_source(path):
    """Increment CSV -> dict(t, yaw (unwrapped, composed), fwd (cum x))."""
    ts, ys, xs = [], [], []
    q = (1.0, 0.0, 0.0, 0.0)
    prev = None
    with open(path) as f:
        for r in csv.reader(f):
            if not r or r[0].startswith('#'):
                continue
            q = quat_mul(q, (float(r[4]), float(r[5]), float(r[6]), float(r[7])))
            n = math.sqrt(sum(c * c for c in q))
            q = tuple(c / n for c in q)
            y = quat_to_yaw(*q)
            if prev is not None:
                y = prev + wrap_pi(y - prev)
            prev = y
            ts.append(int(r[0]) * 1e-9)
            ys.append(y)
            xs.append(float(r[1]))
    return dict(t=np.asarray(ts), yaw=np.asarray(ys), fwd=np.cumsum(xs))


def load_yaw_absolute(path):
    ts, ys = [], []
    prev = None
    with open(path) as f:
        for r in csv.reader(f):
            if not r or r[0].startswith('#'):
                continue
            y = quat_to_yaw(float(r[4]), float(r[5]), float(r[6]), float(r[7]))
            if prev is not None:
                y = prev + wrap_pi(y - prev)
            prev = y
            ts.append(int(r[0]) * 1e-9)
            ys.append(y)
    return np.asarray(ts), np.asarray(ys)


def course_samples(t, E, N, srcs, p, deny=None):
    """Valid anchor course samples (t_mid, course, speed) + rejection counts."""
    if deny is not None:
        keep = (t < deny[0]) | (t > deny[1])
        t, E, N = t[keep], E[keep], N[keep]
    dt = np.diff(t)
    v = np.hypot(np.diff(E), np.diff(N)) / np.maximum(dt, 1e-9)
    dxs = [np.interp(t[1:], s['t'], s['fwd']) - np.interp(t[:-1], s['t'], s['fwd'])
           for s in srcs]
    dx = np.median(np.vstack(dxs), axis=0)
    v_odo = np.abs(dx) / np.maximum(dt, 1e-9)
    dyaws = [np.diff(np.interp(t, s['t'], s['yaw'])) for s in srcs]
    om = np.abs(np.median(np.vstack(dyaws), axis=0)) / np.maximum(dt, 1e-9)
    g_dt = (dt > 0.05) & (dt <= p['dtmax'])
    g_v = (v >= p['vmin']) & (v <= p['vmax'])
    g_xval = np.abs(v - v_odo) < np.maximum(3.0, 0.5 * v)
    g_fwd = dx > 0
    g_om = om < p['omega_max']
    ok = g_dt & g_v & g_xval & g_fwd & g_om
    rej = dict(dt=int((~g_dt).sum()), speed=int((g_dt & ~g_v).sum()),
               xval=int((g_dt & g_v & ~g_xval).sum()),
               reverse=int((g_dt & g_v & g_xval & ~g_fwd).sum()),
               turn=int((g_dt & g_v & g_xval & g_fwd & ~g_om).sum()))
    tm = 0.5 * (t[:-1] + t[1:])
    return tm[ok], np.arctan2(np.diff(N), np.diff(E))[ok], v[ok], rej


def build_pairs(tm, course, srcs, pairgap):
    """Anchor-to-anchor pair stream: (t_pair, dt_pair, r[nsrc, npair])."""
    yaw_mid = [np.interp(tm, s['t'], s['yaw']) for s in srcs]
    tp, dtp, rows = [], [], []
    for k in range(len(tm) - 1):
        gap = tm[k + 1] - tm[k]
        if gap > pairgap:
            continue
        dc = course[k + 1] - course[k]
        rows.append([wrap_pi(dc - (ym[k + 1] - ym[k])) for ym in yaw_mid])
        tp.append(tm[k + 1])
        dtp.append(gap)
    r = np.asarray(rows).T if rows else np.zeros((len(srcs), 0))
    return np.asarray(tp), np.asarray(dtp), r


def wmedian(vals, wts):
    if len(vals) == 0:
        return None
    i = np.argsort(vals)
    v, w = np.asarray(vals)[i], np.asarray(wts)[i]
    cw = np.cumsum(w)
    return v[min(np.searchsorted(cw, 0.5 * cw[-1]), len(v) - 1)]


class SlopeTracker:
    """Online per-source drift-rate tracker: segment/block bookkeeping +
    Theil-Sen-flavored slope pool. Feed (t, rc) pairs; segment breaks are
    signaled by the caller. Evaluate cumulative + windowed slope."""

    def __init__(self, block_s, min_base, win):
        self.block_s = block_s
        self.min_base = min_base
        self.win = win
        self.cum = 0.0
        self.blocks = []      # finalized blocks of CURRENT segment: (t_med, c_med)
        self.cur = []         # open block: list of (t, cum)
        self.pool_s = []      # cumulative slope samples
        self.pool_w = []      # baselines (weights)
        self.wpool = deque()  # (t, slope, baseline) for the windowed estimate

    def _finalize_block(self):
        if not self.cur:
            return
        ts = np.array([t for t, _ in self.cur])
        cs = np.array([c for _, c in self.cur])
        blk = (float(np.median(ts)), float(np.median(cs)))
        for (tb, cb) in self.blocks:
            base = blk[0] - tb
            if base >= self.min_base:
                s = (blk[1] - cb) / base
                self.pool_s.append(s)
                self.pool_w.append(base)
                self.wpool.append((blk[0], s, base))
        self.blocks.append(blk)
        self.cur = []

    def segment_break(self):
        self._finalize_block()
        self.blocks = []

    def add(self, t, rc):
        if self.cur and t - self.cur[0][0] >= self.block_s:
            self._finalize_block()
        self.cum += rc
        self.cur.append((t, self.cum))

    def slope_cum(self):
        return wmedian(self.pool_s, self.pool_w)

    def slope_win(self, now):
        while self.wpool and self.wpool[0][0] < now - self.win:
            self.wpool.popleft()
        if not self.wpool:
            return None
        return wmedian([s for _, s, _ in self.wpool],
                       [b for _, _, b in self.wpool])


def run_monitor(tp, dtp, r, p):
    """Returns per-source series [(t, score_cum, score_win, frozen)], stats."""
    nsrc = r.shape[0]
    if nsrc >= 3:
        m = np.median(r, axis=0)
    else:
        m = np.zeros(r.shape[1])
    dev = r - m
    ev = p['event']
    trackers = [SlopeTracker(p['block_s'], p['min_base'], p['win'])
                for _ in range(nsrc)]
    P = np.zeros(nsrc)
    T = 0.0
    stats = [dict(events=0, event_sum_deg=0.0) for _ in range(nsrc)]
    series = [[] for _ in range(nsrc)]
    prev_t = None
    for k in range(r.shape[1]):
        t = tp[k]
        if prev_t is not None and t - prev_t > p['pairgap']:
            for tr in trackers:
                tr.segment_break()
        prev_t = t
        T += dtp[k]
        for i in range(nsrc):
            d = dev[i, k]
            excess = max(abs(d) - ev, 0.0)
            if excess > 0:
                P[i] += excess
                stats[i]['events'] += 1
                stats[i]['event_sum_deg'] += math.degrees(excess)
            trackers[i].add(t, m[k] + max(-ev, min(ev, d)))
            sc = trackers[i].slope_cum()
            sw = trackers[i].slope_win(t)
            cum_score = (abs(sc) + P[i] / T) if sc is not None else None
            win_score = (abs(sw) + P[i] / T) if sw is not None else None
            series[i].append((t, cum_score, win_score, sw is None))
    return series, stats


def deg_h(rad_per_s):
    return math.degrees(rad_per_s) * 3600.0


def last_value(series, idx):
    return next((e[idx] for e in reversed(series) if e[idx] is not None), None)


def rank_of(scores):
    return tuple(sorted(range(len(scores)), key=lambda i: scores[i]))


def discovery_latency(series):
    """Earliest t after which the cumulative-score ranking equals the final
    ranking at every later eval point. Returns (t, final_rank)."""
    if not series or not series[0]:
        return None, None
    finals = [last_value(s, 1) for s in series]
    if any(f is None for f in finals):
        return None, None
    final = rank_of(finals)
    stable_from = None
    for k in range(len(series[0])):
        cs = [s[k][1] for s in series]
        if any(c is None for c in cs):
            continue
        if rank_of(cs) == final:
            if stable_from is None:
                stable_from = series[0][k][0]
        else:
            stable_from = None
    return stable_from, final


def gt_reference(gt_path, srcs):
    """Final drift rate (rad/s) of each source's integrated yaw vs GT yaw."""
    gt_t, gt_y = load_yaw_absolute(gt_path)
    rates = []
    for s in srcs:
        t0 = max(s['t'][0], gt_t[0])
        t1 = min(s['t'][-1], gt_t[-1])
        off = np.interp(t0, gt_t, gt_y) - np.interp(t0, s['t'], s['yaw'])
        e1 = (np.interp(t1, s['t'], s['yaw']) + off) - np.interp(t1, gt_t, gt_y)
        rates.append(e1 / (t1 - t0))
    return rates


def weight_rule(scores_deg_h, wcap=10.0, floor_deg_h=5.0):
    """rot_weight_i = clip(wcap * min_j score_j / score_i, 1, wcap)."""
    mags = [max(s, floor_deg_h) for s in scores_deg_h]
    best = min(mags)
    return [min(wcap, max(1.0, wcap * best / m)) for m in mags]


def analyze(tg, E, N, srcs, p, deny=None):
    tm, course, v, rej = course_samples(tg, E, N, srcs, p, deny)
    tp, dtp, r = build_pairs(tm, course, srcs, p['pairgap'])
    series, stats = run_monitor(tp, dtp, r, p)
    lat, final = discovery_latency(series)
    return dict(tm=tm, v=v, rej=rej, tp=tp, dtp=dtp, r=r,
                series=series, stats=stats, latency=lat, final=final)


def main():
    argv = sys.argv[1:]
    p = dict(vmin=3.0, vmax=30.0, dtmax=3.0, omega_max=math.radians(3.0),
             pairgap=60.0, event=math.radians(5.0), block_s=60.0,
             min_base=120.0, win=600.0)
    opt = dict(gt=None, deny=None, sweep=False, series=None)
    paths = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--gt':
            opt['gt'] = argv[i + 1]; i += 2
        elif a == '--vmin':
            p['vmin'] = float(argv[i + 1]); i += 2
        elif a == '--vmax':
            p['vmax'] = float(argv[i + 1]); i += 2
        elif a == '--dtmax':
            p['dtmax'] = float(argv[i + 1]); i += 2
        elif a == '--omega':
            p['omega_max'] = math.radians(float(argv[i + 1])); i += 2
        elif a == '--pairgap':
            p['pairgap'] = float(argv[i + 1]); i += 2
        elif a == '--event':
            p['event'] = math.radians(float(argv[i + 1])); i += 2
        elif a == '--block':
            p['block_s'] = float(argv[i + 1]); i += 2
        elif a == '--minbase':
            p['min_base'] = float(argv[i + 1]); i += 2
        elif a == '--win':
            p['win'] = float(argv[i + 1]); i += 2
        elif a == '--deny':
            opt['deny'] = (float(argv[i + 1]), float(argv[i + 2])); i += 3
        elif a == '--sweep':
            opt['sweep'] = True; i += 1
        elif a == '--series':
            opt['series'] = argv[i + 1]; i += 2
        else:
            paths.append(a); i += 1
    if len(paths) < 2:
        sys.exit('usage: proto_heading_monitor.py GPS.csv SRC.csv [SRC2.csv ...]'
                 ' [--gt GT.csv] [--vmin M] [--omega DPS] [--pairgap S]'
                 ' [--win S] [--deny T0 T1] [--sweep] [--series PREFIX]')
    gps_path, src_paths = paths[0], paths[1:]
    names = [q.split('/')[-1].split('\\')[-1].rsplit('.', 1)[0] for q in src_paths]

    tg, E, N = load_gps(gps_path)
    srcs = [load_source(q) for q in src_paths]
    if len(srcs) < 3:
        print('WARNING: < 3 sources, median consensus degenerates (m = 0)')

    print('== monitor (vmin=%.1f m/s, omega_max=%.1f deg/s, pairgap=%.0f s,'
          ' event=%.0f deg, block=%.0f s, win=%.0f s)'
          % (p['vmin'], math.degrees(p['omega_max']), p['pairgap'],
             math.degrees(p['event']), p['block_s'], p['win']))
    a = analyze(tg, E, N, srcs, p, opt['deny'])
    med_v = np.median(a['v']) if len(a['v']) else 0.0
    floor = math.degrees(math.sqrt(2) * SIGMA_POS / max(med_v, 1e-9))
    print('  GPS: %d fixes -> %d anchors (rej dt/speed/xval/reverse/turn'
          ' %d/%d/%d/%d/%d), %d pairs, med speed %.1f m/s,'
          ' per-anchor course noise ~%.1f deg'
          % (len(tg), len(a['tm']), a['rej']['dt'], a['rej']['speed'],
             a['rej']['xval'], a['rej']['reverse'], a['rej']['turn'],
             len(a['tp']), med_v, floor))
    cums = []
    for i, nm in enumerate(names):
        c = last_value(a['series'][i], 1)
        w = last_value(a['series'][i], 2)
        cums.append(deg_h(c) if c is not None else float('nan'))
        print('  %-12s score cum %8.1f deg/h  win %8.1f deg/h  events %3d'
              ' (%.1f deg)' % (nm, cums[-1],
                               deg_h(w) if w is not None else float('nan'),
                               a['stats'][i]['events'],
                               a['stats'][i]['event_sum_deg']))
    if a['final'] is not None:
        print('  ranking (best->worst): %s   discovery latency: %s'
              % (' > '.join(names[i] for i in a['final']),
                 '%.0f s' % a['latency'] if a['latency'] is not None
                 else 'never stable'))
        w = weight_rule(cums)
        print('  candidate rot_weights (wcap=10): %s'
              % ', '.join('%s=%.1f' % (n, x) for n, x in zip(names, w)))

    gt_order = None
    if opt['gt']:
        gt_rates = gt_reference(opt['gt'], srcs)
        gt_order = rank_of([abs(x) for x in gt_rates])
        print('== GT reference (integrated source yaw vs GT yaw)')
        for nm, rr in zip(names, gt_rates):
            print('  %-12s %+9.1f deg/h' % (nm, deg_h(rr)))
        print('  GT ranking: %s   monitor full order %s, top-1 %s'
              % (' > '.join(names[i] for i in gt_order),
                 'MATCHES' if a['final'] == gt_order else 'MISMATCH',
                 'MATCHES' if a['final'] and a['final'][0] == gt_order[0]
                 else 'MISMATCH'))

    if opt['deny']:
        d0, d1 = opt['deny']
        print('== denial probe: GPS masked over [%g, %g] s' % (d0, d1))
        inside = int(np.sum((a['tp'] >= d0) & (a['tp'] <= d1)))
        print('  pairs inside masked window: %d (expect 0)' % inside)
        for i, nm in enumerate(names):
            def cum_at(tq, i=i):
                ks = [e[1] for e in a['series'][i]
                      if e[0] <= tq and e[1] is not None]
                return deg_h(ks[-1]) if ks else float('nan')
            frz = [e[3] for e in a['series'][i] if d0 + p['win'] < e[0] < d1]
            print('  %-12s cum before %8.1f / at-end %8.1f / final %8.1f deg/h'
                  '  win-frozen inside: %s'
                  % (nm, cum_at(d0), cum_at(d1), cum_at(1e12),
                     'yes' if all(frz) else ('n/a' if not frz else 'NO')))

    if opt['series']:
        for i, nm in enumerate(names):
            out = opt['series'] + '_' + nm + '_mon.csv'
            with open(out, 'w') as f:
                f.write('# t_s,cum_score_deg_h,win_score_deg_h,win_frozen\n')
                for (t, c, w, fz) in a['series'][i]:
                    f.write('%.3f,%s,%s,%d\n'
                            % (t, '%.2f' % deg_h(c) if c is not None else '',
                               '%.2f' % deg_h(w) if w is not None else '',
                               int(fz)))
            print('  wrote %s' % out)

    if opt['sweep']:
        print('== sweep: final ranking + latency per setting%s'
              % (' ([order/top1] vs GT)' if gt_order is not None else ''))
        for vmin in (2.0, 3.0, 5.0):
            for omax in (1.5, 3.0, 10.0):
                for pg in (30.0, 60.0, 120.0):
                    p2 = dict(p)
                    p2['vmin'] = vmin
                    p2['omega_max'] = math.radians(omax)
                    p2['pairgap'] = pg
                    a2 = analyze(tg, E, N, srcs, p2, opt['deny'])
                    rk = (' > '.join(names[i] for i in a2['final'])
                          if a2['final'] else 'n/a')
                    ok = ''
                    if gt_order is not None and a2['final'] is not None:
                        ok = '  [%s/%s]' % (
                            'OK' if a2['final'] == gt_order else 'flip',
                            'OK' if a2['final'][0] == gt_order[0] else 'MISS')
                    print('  vmin=%.0f omax=%4.1f pg=%3.0f : %s  latency %s%s'
                          % (vmin, omax, pg, rk,
                             '%.0f s' % a2['latency']
                             if a2['latency'] is not None else 'never', ok))


if __name__ == '__main__':
    main()
