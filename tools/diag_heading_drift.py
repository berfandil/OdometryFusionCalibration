# Diagnose per-source heading (yaw) drift vs GT on a KAIST run — the urban12 522 s
# GPS-coast residual question: which odometry source's heading drifts, how fast, and
# would the fused median have had a better heading available?
# Usage: diag_heading_drift.py GT.csv SRC.csv [SRC2.csv ...] [--coast t0 t1]
import sys, csv, math

def load_rows(path):
    out = []
    with open(path) as f:
        for r in csv.reader(f):
            if not r or r[0].startswith('#'): continue
            out.append((int(r[0]), [float(x) for x in r[1:8]]))
    return out

def quat_to_yaw(qw,qx,qy,qz):
    # ZYX yaw of the rotation
    return math.atan2(2*(qw*qz+qx*qy), 1-2*(qy*qy+qz*qz))

def quat_mul(a,b):
    aw,ax,ay,az = a; bw,bx,by,bz = b
    return (aw*bw-ax*bx-ay*by-az*bz, aw*bx+ax*bw+ay*bz-az*by,
            aw*by-ax*bz+ay*bw+az*bx, aw*bz+ax*by-ay*bx+az*bw)

def unwrap(prev, y):
    while y - prev > math.pi:  y -= 2*math.pi
    while y - prev < -math.pi: y += 2*math.pi
    return y

args = sys.argv[1:]
coast = None
if '--coast' in args:
    i = args.index('--coast'); coast = (float(args[i+1]), float(args[i+2])); args = args[:i]
gt_path, src_paths = args[0], args[1:]

gt = load_rows(gt_path)
t0_ns = gt[0][0]
gt_yaw = []
prev = None
for t, v in gt:
    y = quat_to_yaw(v[3],v[4],v[5],v[6])
    if prev is not None: y = unwrap(prev, y)
    prev = y
    gt_yaw.append(((t - t0_ns)*1e-9, y))

def gt_yaw_at(ts):
    # linear scan ok for diagnostics; gt sorted
    import bisect
    times = [a for a,_ in gt_yaw]
    i = bisect.bisect_left(times, ts)
    if i <= 0: return gt_yaw[0][1]
    if i >= len(gt_yaw): return gt_yaw[-1][1]
    (ta,ya),(tb,yb) = gt_yaw[i-1], gt_yaw[i]
    f = (ts-ta)/(tb-ta) if tb>ta else 0
    return ya + f*(yb-ya)

print(f"GT span: {gt_yaw[0][0]:.1f} .. {gt_yaw[-1][0]:.1f} s")
for sp in src_paths:
    rows = load_rows(sp)
    q = (1.0,0.0,0.0,0.0)
    yaw_traj = []
    prev = None
    for t, v in rows:
        q = quat_mul(q, (v[3],v[4],v[5],v[6]))
        n = math.sqrt(sum(c*c for c in q)); q = tuple(c/n for c in q)
        ts = (t - t0_ns)*1e-9
        y = quat_to_yaw(*q)
        if prev is not None: y = unwrap(prev, y)
        prev = y
        yaw_traj.append((ts, y))
    # align integrated yaw to GT at the first shared time
    ts0 = max(yaw_traj[0][0], gt_yaw[0][0])
    import bisect
    times = [a for a,_ in yaw_traj]
    def src_yaw_at(ts):
        i = bisect.bisect_left(times, ts)
        if i <= 0: return yaw_traj[0][1]
        if i >= len(yaw_traj): return yaw_traj[-1][1]
        (ta,ya),(tb,yb) = yaw_traj[i-1], yaw_traj[i]
        f = (ts-ta)/(tb-ta) if tb>ta else 0
        return ya + f*(yb-ya)
    off = gt_yaw_at(ts0) - src_yaw_at(ts0)
    # error trajectory at 10 s samples
    t_end = min(yaw_traj[-1][0], gt_yaw[-1][0])
    errs = []
    ts = ts0
    while ts <= t_end:
        e = (src_yaw_at(ts) + off) - gt_yaw_at(ts)
        errs.append((ts, e))
        ts += 10.0
    emax = max(abs(e) for _,e in errs)
    efin = errs[-1][1]
    line = f"{sp:32s} full-drive: final {math.degrees(efin):+8.2f} deg, max {math.degrees(emax):7.2f} deg, rate {math.degrees(efin)/( (t_end-ts0)/3600):+8.1f} deg/h"
    if coast:
        e0 = (src_yaw_at(coast[0]) + off) - gt_yaw_at(coast[0])
        e1 = (src_yaw_at(coast[1]) + off) - gt_yaw_at(coast[1])
        line += f" | coast[{coast[0]:.0f},{coast[1]:.0f}]s drift: {math.degrees(e1-e0):+7.2f} deg"
    print(line)
