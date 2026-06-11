# Prototype: turn-regime rotation-extrinsic via rotation-axis correspondence (Wahba).
# Math under test:  delta_s = E^-1 o delta_ref o E  =>  R_B = R_E^T R_A R_E
#                   =>  log(R_A) = R_E log(R_B)  (axis-vector correspondence)
# Solve R_E by SVD of M = sum w_i * a_i b_i^T (Wahba / Kabsch), gate on rotation magnitude.
# Also: rank/conditioning analysis of sum(b b^T) -> per-DOF observability (yaw-only ground case).
import sys, math, csv
import numpy as np

def load_increments(path):
    rows = []
    with open(path) as f:
        rd = csv.reader(f)
        for r in rd:
            if not r or r[0].startswith('#'):
                continue
            t = int(r[0]); q = [float(x) for x in r[4:8]]  # qw qx qy qz
            rows.append((t, q))
    return rows

def quat_to_R(qw, qx, qy, qz):
    n = math.sqrt(qw*qw+qx*qx+qy*qy+qz*qz) or 1.0
    qw, qx, qy, qz = qw/n, qx/n, qy/n, qz/n
    return np.array([
        [1-2*(qy*qy+qz*qz), 2*(qx*qy-qz*qw),   2*(qx*qz+qy*qw)],
        [2*(qx*qy+qz*qw),   1-2*(qx*qx+qz*qz), 2*(qy*qz-qx*qw)],
        [2*(qx*qz-qy*qw),   2*(qy*qz+qx*qw),   1-2*(qx*qx+qy*qy)]])

def so3_log(R):
    c = max(-1.0, min(1.0, (np.trace(R)-1)/2))
    th = math.acos(c)
    if th < 1e-12:
        return np.zeros(3)
    v = np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]])
    return v * (th/(2*math.sin(th)))

def Rz(a): c,s=math.cos(a),math.sin(a); return np.array([[c,-s,0],[s,c,0],[0,0,1]])
def Ry(a): c,s=math.cos(a),math.sin(a); return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def Rx(a): c,s=math.cos(a),math.sin(a); return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def so3_exp(v):
    th = np.linalg.norm(v)
    if th < 1e-12:
        return np.eye(3)
    k = v/th
    K = np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[0],k[1]*0+0,0]])
    K = np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]])
    return np.eye(3) + math.sin(th)*K + (1-math.cos(th))*(K@K)

ref_csv, inj_csv = sys.argv[1], sys.argv[2]
yaw, pitch, roll = (math.radians(float(sys.argv[i])) for i in (3,4,5))
win   = int(sys.argv[6])   if len(sys.argv) > 6 else 100     # steps per window (0.5 s @200Hz)
gate  = float(sys.argv[7]) if len(sys.argv) > 7 else 0.05    # rad per WINDOW
noise = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0     # synth rot noise rad/STEP
E = Rz(yaw) @ Ry(pitch) @ Rx(roll)

ref = load_increments(ref_csv)
inj = load_increments(inj_csv)
nrow = min(len(ref), len(inj))
rng = np.random.default_rng(42)

M = np.zeros((3,3)); BB = np.zeros((3,3)); used = 0; n = 0
for w0 in range(0, nrow - win, win):
    Ra = np.eye(3); Rb = np.eye(3)
    for i in range(w0, w0 + win):
        Ria = quat_to_R(*ref[i][1])
        Rib = quat_to_R(*inj[i][1])
        if noise > 0:
            Ria = Ria @ so3_exp(rng.normal(0, noise, 3))
            Rib = Rib @ so3_exp(rng.normal(0, noise, 3))
        Ra = Ra @ Ria
        Rb = Rb @ Rib
    n += 1
    a = so3_log(Ra); b = so3_log(Rb)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < gate or nb < gate:
        continue
    M += np.outer(a, b)          # weight = |angle| baked in (a,b unnormalized)
    BB += np.outer(b, b)
    used += 1

U, S, Vt = np.linalg.svd(M)
D = np.diag([1, 1, np.linalg.det(U @ Vt)])
RX = U @ D @ Vt

err = so3_log(RX @ E.T)
ev = np.linalg.eigvalsh(BB)
print(f"pairs used        : {used}/{n} (gate {gate} rad/step, noise {noise})")
print(f"BB eigvals        : {ev[0]:.4g} {ev[1]:.4g} {ev[2]:.4g}  (cond lo/hi {ev[0]/ev[2]:.3g})")
print(f"R_X recovered err : {np.linalg.norm(err):.6g} rad = {math.degrees(np.linalg.norm(err)):.4g} deg")
print(f"  err vec         : {err}")
# per-Euler readout for intuition
def euler_zyx(R):
    return (math.atan2(R[1,0], R[0,0]),
            math.atan2(-R[2,0], math.hypot(R[0,0], R[1,0])),
            math.atan2(R[2,1], R[2,2]))
ey, ep, er = euler_zyx(RX)
print(f"  euler zyx (deg) : yaw {math.degrees(ey):.4f} pitch {math.degrees(ep):.4f} roll {math.degrees(er):.4f}"
      f"  (truth {math.degrees(yaw):.4f} {math.degrees(pitch):.4f} {math.degrees(roll):.4f})")
