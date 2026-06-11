# Prototype: turn-regime JOINT lever + scale (Slice 17b).
# Hand-eye translation identity for A = X o B_true o X^-1 with reported t_rep = s * t_true:
#   (R_A - I) t_X - kappa * (R_X t_rep) = -t_A,   kappa = 1/s
# Linear in u = [t_X; kappa] (4 unknowns, 3 equations per window). Accumulate 4x4 normal
# equations over turn-gated windows; recover lever t_X and scale s = 1/kappa.
# R_X comes from the rot3d Wahba solve (run jointly here, mirroring Slice 17).
# Usage: proto_lever_scale.py REF.csv INJ.csv yaw_deg pitch_deg roll_deg lx ly lz scale [win] [gate]
import sys, math, csv
import numpy as np

def load(path):
    rows = []
    with open(path) as f:
        for r in csv.reader(f):
            if not r or r[0].startswith('#'): continue
            rows.append(([float(x) for x in r[1:4]], [float(x) for x in r[4:8]]))
    return rows

def quat_to_R(qw,qx,qy,qz):
    n = math.sqrt(qw*qw+qx*qx+qy*qy+qz*qz) or 1.0
    qw,qx,qy,qz = qw/n,qx/n,qy/n,qz/n
    return np.array([[1-2*(qy*qy+qz*qz),2*(qx*qy-qz*qw),2*(qx*qz+qy*qw)],
                     [2*(qx*qy+qz*qw),1-2*(qx*qx+qz*qz),2*(qy*qz-qx*qw)],
                     [2*(qx*qz-qy*qw),2*(qy*qz+qx*qw),1-2*(qx*qx+qy*qy)]])

def so3_log(R):
    c = max(-1.0, min(1.0, (np.trace(R)-1)/2)); th = math.acos(c)
    if th < 1e-12: return np.zeros(3)
    return np.array([R[2,1]-R[1,2],R[0,2]-R[2,0],R[1,0]-R[0,1]])*(th/(2*math.sin(th)))

def Rz(a): c,s=math.cos(a),math.sin(a); return np.array([[c,-s,0],[s,c,0],[0,0,1]])
def Ry(a): c,s=math.cos(a),math.sin(a); return np.array([[c,0,s],[0,1,0],[-s,0,c]])
def Rx(a): c,s=math.cos(a),math.sin(a); return np.array([[1,0,0],[0,c,-s],[0,s,c]])

ref_csv, inj_csv = sys.argv[1], sys.argv[2]
yaw,pitch,roll = (math.radians(float(sys.argv[i])) for i in (3,4,5))
lever_t = np.array([float(sys.argv[i]) for i in (6,7,8)])
s_t = float(sys.argv[9])
win  = int(sys.argv[10])  if len(sys.argv) > 10 else 100
gate = float(sys.argv[11]) if len(sys.argv) > 11 else 0.05
E = Rz(yaw)@Ry(pitch)@Rx(roll)

ref = load(ref_csv); inj = load(inj_csv)
nrow = min(len(ref), len(inj))

# Pass 1: rot3d Wahba (Slice-17 mirror) for R_X.
M = np.zeros((3,3))
wins = []
for w0 in range(0, nrow - win, win):
    Ra = np.eye(3); ta = np.zeros(3); Rb = np.eye(3); tb = np.zeros(3)
    for i in range(w0, w0+win):
        Ria = quat_to_R(*ref[i][1]); tia = np.array(ref[i][0])
        Rib = quat_to_R(*inj[i][1]); tib = np.array(inj[i][0])
        ta = Ra@tia + ta; Ra = Ra@Ria
        tb = Rb@tib + tb; Rb = Rb@Rib
    a = so3_log(Ra); b = so3_log(Rb)
    if np.linalg.norm(a) < gate or np.linalg.norm(b) < gate: continue
    M += np.outer(a,b)
    wins.append((Ra, ta, tb))
U,S,Vt = np.linalg.svd(M)
RX = U@np.diag([1,1,np.linalg.det(U@Vt)])@Vt
print(f"R_X err: {math.degrees(np.linalg.norm(so3_log(RX@E.T))):.3g} deg   windows: {len(wins)}")

# Pass 2: joint 4-unknown LS  [(R_A - I) | -R_X t_rep] u = -t_A
AtA = np.zeros((4,4)); Atb = np.zeros(4)
for Ra, ta, tb in wins:
    J = np.hstack([Ra - np.eye(3), (-(RX@tb)).reshape(3,1)])
    AtA += J.T@J; Atb += J.T@(-ta)
ev = np.linalg.eigvalsh(AtA)
# Planar/rank-deficient guard (mirrors production ridge + Phase-1-supplied R_X on ground):
if ev[0] < 1e-9 * ev[3]:
    print("  [rank-deficient: re-solving with truth R_X (Phase-1 ground path) + ridge]")
    AtA = np.zeros((4,4)); Atb = np.zeros(4)
    for Ra, ta, tb in wins:
        J = np.hstack([Ra - np.eye(3), (-(E@tb)).reshape(3,1)])
        AtA += J.T@J; Atb += J.T@(-ta)
    ev = np.linalg.eigvalsh(AtA)
    AtA = AtA + 1e-9 * ev[3] * np.eye(4)
u = np.linalg.solve(AtA, Atb)
tX, kappa = u[:3], u[3]
print(f"AtA eigvals       : {ev[0]:.4g} {ev[1]:.4g} {ev[2]:.4g} {ev[3]:.4g}  (lo/hi {ev[0]/ev[3]:.3g})")
print(f"lever recovered   : [{tX[0]:.6f}, {tX[1]:.6f}, {tX[2]:.6f}]  truth {lever_t}  err {np.linalg.norm(tX-lever_t)*1000:.3g} mm")
print(f"scale recovered   : {1.0/kappa:.6f}  truth {s_t}  err {abs(1.0/kappa - s_t):.2e}")
# Control: 3-unknown solve with raw t_rep (the Slice-17 status quo) for comparison
AtA3 = np.zeros((3,3)); Atb3 = np.zeros(3)
for Ra, ta, tb in wins:
    A3 = Ra - np.eye(3)
    AtA3 += A3.T@A3; Atb3 += A3.T@(RX@tb - ta)
t3 = np.linalg.solve(AtA3, Atb3)
print(f"3-unknown control : [{t3[0]:.4f}, {t3[1]:.4f}, {t3[2]:.4f}]  err {np.linalg.norm(t3-lever_t)*1000:.3g} mm")
