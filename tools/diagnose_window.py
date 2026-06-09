# Localize urban12's divergence event from an ofc_replay results CSV (Slice 15 #2).
# Eval-only: finds the step-window where the global error leaps the most, then dumps that
# window's GPS activity (applied/rejected fixes + NIS) and heading/translation arc so we can
# tell a GPS-denied dead-reckoning COAST from an admitted-fix YANK. No core code involved.
import csv, sys, math

path = sys.argv[1] if len(sys.argv) > 1 else "kaist_run/urban12_k0f25_out.csv"
WIN  = int(sys.argv[2]) if len(sys.argv) > 2 else 1000   # window length (steps), matches local metric
WARM = 20

# --- load the per-fused-step rows -----------------------------------------------------------
rows = []
with open(path, newline="") as f:
    rd = csv.reader(f)
    hdr = None
    for r in rd:
        if not r or r[0].startswith("#"):
            continue
        if hdr is None and r[0] == "now_ns":
            hdr = {name: i for i, name in enumerate(r)}
            continue
        if hdr is None:
            continue
        rows.append(r)
col = lambda r, k: r[hdr[k]]
t0_ns = int(rows[0][hdr["now_ns"]])
def t_s(r): return (int(col(r, "now_ns")) - t0_ns) / 1e9

# corr_applied/rejected: per-step (0/1) or cumulative? detect.
ca = [int(col(r, "corr_applied")) for r in rows]
cumulative = ca[-1] > 5 and all(ca[i] <= ca[i+1] for i in range(0, len(ca)-1, max(1, len(ca)//50)))
def applied_at(i):  # 1 if a fix was APPLIED at step i
    return (ca[i] - ca[i-1] > 0) if (cumulative and i > 0) else (ca[i] == 1 if not cumulative else ca[i] > 0)
cr = [int(col(r, "corr_rejected")) for r in rows]
def rejected_at(i):
    return (cr[i] - cr[i-1] > 0) if (cumulative and i > 0) else (cr[i] >= 1 and not cumulative)

n = len(rows)
print(f"# rows={n}  span={t_s(rows[-1]):.1f}s  corr_counts={'cumulative' if cumulative else 'per-step'}")

# --- coarse arc of the whole drive ----------------------------------------------------------
print("\n# coarse trans_err arc (every ~5%):")
for k in range(0, n, max(1, n // 20)):
    r = rows[k]
    print(f"  t={t_s(r):7.1f}s  trans_err={float(col(r,'trans_err_m')):10.2f}m  rot_err={float(col(r,'rot_err_rad')):6.3f}rad")

# --- per-window GLOBAL-error growth; find the steepest --------------------------------------
def te(i): return float(col(rows[i], "trans_err_m"))
def re_(i): return float(col(rows[i], "rot_err_rad"))
wins = []
for s in range(WARM, n - WIN + 1, WIN):
    e = s + WIN - 1
    d_trans = te(e) - te(s)
    d_rot   = re_(e) - re_(s)
    napp = sum(1 for i in range(s, e+1) if applied_at(i))
    nrej = sum(1 for i in range(s, e+1) if rejected_at(i))
    wins.append((s, e, d_trans, d_rot, napp, nrej))

wins_by_trans = sorted(wins, key=lambda w: w[2], reverse=True)
print("\n# top-5 windows by GLOBAL trans-error GROWTH (the divergence events):")
print("  win_start_s  win_end_s   d_trans_m   d_rot_rad  applied  rejected  trans_end_m")
for (s, e, dt, dr, na, nr) in wins_by_trans[:5]:
    print(f"  {t_s(rows[s]):10.1f}  {t_s(rows[e]):8.1f}  {dt:11.1f}  {dr:9.3f}  {na:7d}  {nr:8d}  {te(e):11.1f}")

# --- top windows by HEADING-error growth (the ROOT, since trans follows heading) -----------
wins_by_rot = sorted(wins, key=lambda w: abs(w[3]), reverse=True)
print("\n# top-5 windows by |HEADING-error growth| (the root cause -- trans follows heading):")
print("  win_start_s  win_end_s   d_rot_rad   d_trans_m  applied  rejected  rot_end_rad")
for (s, e, dt, dr, na, nr) in wins_by_rot[:5]:
    print(f"  {t_s(rows[s]):10.1f}  {t_s(rows[e]):8.1f}  {dr:9.3f}  {dt:10.1f}  {na:7d}  {nr:8d}  {re_(e):9.3f}")

def zoom(s, e, label):
    dt_, dr_ = te(e)-te(s), re_(e)-re_(s)
    na_ = sum(1 for i in range(s, e+1) if applied_at(i))
    nr_ = sum(1 for i in range(s, e+1) if rejected_at(i))
    print(f"\n# {label}: t=[{t_s(rows[s]):.1f}, {t_s(rows[e]):.1f}]s  d_trans={dt_:.1f}m  d_rot={dr_:.3f}rad  applied={na_} rejected={nr_}")
    last_app = None
    for i in range(s, max(0, s-60000), -1):
        if applied_at(i): last_app = i; break
    print(("  entering-coast: last applied fix at t=%.1fs (%.1fs before)" % (t_s(rows[last_app]), t_s(rows[s])-t_s(rows[last_app]))) if last_app is not None else "  entering-coast: NO prior applied fix")
    step = max(1, (e-s) // 30)
    for i in range(s, e+1, step):
        mark = "  <-- APPLIED" if applied_at(i) else ("  (rej)" if rejected_at(i) else "")
        print(f"    t={t_s(rows[i]):7.1f}s  trans_err={te(i):9.2f}m  rot_err={re_(i):6.3f}rad  NIS={float(col(rows[i],'last_nis')):9.2f}{mark}")
    apps = [(i, float(col(rows[i],'last_nis'))) for i in range(s, e+1) if applied_at(i)]
    print(f"  applied-fix count in window: {len(apps)}" + (f"  NIS range [{min(a[1] for a in apps):.2f},{max(a[1] for a in apps):.2f}]" if apps else "  (pure DR coast)"))

# Zoom the worst HEADING window (root) and the worst TRANS window (symptom tail).
zoom(*wins_by_rot[0][:2], "WORST HEADING window (root)")

s, e, dt, dr, na, nr = wins_by_trans[0]
print(f"\n# WORST window: t=[{t_s(rows[s]):.1f}, {t_s(rows[e]):.1f}]s  d_trans={dt:.1f}m  d_rot={dr:.3f}rad  applied={na} rejected={nr}")
# coasting gap entering the window: steps since the last applied fix before s
last_app = None
for i in range(s, max(0, s-50000), -1):
    if applied_at(i): last_app = i; break
if last_app is not None:
    print(f"  entering-coast: last applied fix at t={t_s(rows[last_app]):.1f}s ({t_s(rows[s])-t_s(rows[last_app]):.1f}s before window)")
else:
    print("  entering-coast: NO applied fix before this window")

print("  trace (downsampled ~40 pts) + GPS events:")
step = max(1, WIN // 40)
for i in range(s, e+1, step):
    r = rows[i]
    print(f"    t={t_s(r):7.1f}s  trans_err={te(i):9.2f}m  rot_err={re_(i):6.3f}rad")
# list every applied/rejected fix in the window with its NIS
print("  GPS fixes in window (applied & rejected, with NIS):")
shown = 0
for i in range(s, e+1):
    a, rj = applied_at(i), rejected_at(i)
    if a or rj:
        print(f"    t={t_s(rows[i]):7.1f}s  {'APPLIED ' if a else 'rejected'}  NIS={float(col(rows[i],'last_nis')):10.2f}  trans_err={te(i):9.2f}m rot_err={re_(i):6.3f}")
        shown += 1
        if shown > 40:
            print("    ... (truncated)")
            break
if shown == 0:
    print("    (none -- pure dead-reckoning coast)")
