# Inject a KNOWN extrinsic + scale + time-offset into an increment-form source CSV, so online
# calibration can be tested on REAL motion with a known ground-truth answer (Slice 13 #1).
#
# A sensor mounted at extrinsic E (reference_from_sensor) measures its body motion as
#   delta_sensor = E^-1 o delta_ref o E        (adjoint / conjugation)
# so we transform every per-step increment that way, then scale its translation and shift its
# timestamps. The calibrator should recover (E, scale, time_offset) as that source's CalibSnapshot.
#
# Usage:
#   python inject_calib.py IN.csv OUT.csv  yaw_deg pitch_deg roll_deg  lx ly lz  scale  toff_s
import sys, math

inp, out = sys.argv[1], sys.argv[2]
yaw, pitch, roll = (math.radians(float(sys.argv[i])) for i in (3, 4, 5))
lever = [float(sys.argv[i]) for i in (6, 7, 8)]
scale = float(sys.argv[9])
toff_ns = int(round(float(sys.argv[10]) * 1e9))

# --- minimal SO3/SE3 in plain lists ---------------------------------------------------------
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def matvec(A, v):
    return [sum(A[i][k] * v[k] for k in range(3)) for i in range(3)]
def transpose(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]
def Rx(a): c,s=math.cos(a),math.sin(a); return [[1,0,0],[0,c,-s],[0,s,c]]
def Ry(a): c,s=math.cos(a),math.sin(a); return [[c,0,s],[0,1,0],[-s,0,c]]
def Rz(a): c,s=math.cos(a),math.sin(a); return [[c,-s,0],[s,c,0],[0,0,1]]

# E.R = Rz(yaw) Ry(pitch) Rx(roll)  (matches config_loader parse_extrinsic "yaw pitch roll x y z").
ER = matmul(matmul(Rz(yaw), Ry(pitch)), Rx(roll))
Et = lever
ERt = transpose(ER)                      # E^-1 rotation
Einv_t = [-x for x in matvec(ERt, Et)]   # E^-1 translation = -R^T t

def quat_to_R(qw, qx, qy, qz):
    n = math.sqrt(qw*qw+qx*qx+qy*qy+qz*qz) or 1.0
    qw,qx,qy,qz = qw/n,qx/n,qy/n,qz/n
    return [[1-2*(qy*qy+qz*qz), 2*(qx*qy-qz*qw),   2*(qx*qz+qy*qw)],
            [2*(qx*qy+qz*qw),   1-2*(qx*qx+qz*qz), 2*(qy*qz-qx*qw)],
            [2*(qx*qz-qy*qw),   2*(qy*qz+qx*qw),   1-2*(qx*qx+qy*qy)]]
def R_to_quat(R):
    tr = R[0][0]+R[1][1]+R[2][2]
    if tr > 0:
        s = math.sqrt(tr+1.0)*2; qw=0.25*s
        qx=(R[2][1]-R[1][2])/s; qy=(R[0][2]-R[2][0])/s; qz=(R[1][0]-R[0][1])/s
    elif R[0][0]>R[1][1] and R[0][0]>R[2][2]:
        s=math.sqrt(1.0+R[0][0]-R[1][1]-R[2][2])*2; qw=(R[2][1]-R[1][2])/s
        qx=0.25*s; qy=(R[0][1]+R[1][0])/s; qz=(R[0][2]+R[2][0])/s
    elif R[1][1]>R[2][2]:
        s=math.sqrt(1.0+R[1][1]-R[0][0]-R[2][2])*2; qw=(R[0][2]-R[2][0])/s
        qx=(R[0][1]+R[1][0])/s; qy=0.25*s; qz=(R[1][2]+R[2][1])/s
    else:
        s=math.sqrt(1.0+R[2][2]-R[0][0]-R[1][1])*2; qw=(R[1][0]-R[0][1])/s
        qx=(R[0][2]+R[2][0])/s; qy=(R[1][2]+R[2][1])/s; qz=0.25*s
    return qw,qx,qy,qz

# SE3 compose (R3 = R1 R2, t3 = R1 t2 + t1)
def se3_mul(R1,t1, R2,t2):
    return matmul(R1,R2), [matvec(R1,t2)[i]+t1[i] for i in range(3)]

n_out = 0
with open(inp) as f, open(out, "w", newline="") as g:
    g.write("# form: increment\n# t_ns, x, y, z, qw, qx, qy, qz  (calib-injected)\n")
    for line in f:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split(",")
        t = int(p[0]); x,y,z = float(p[1]),float(p[2]),float(p[3])
        qw,qx,qy,qz = float(p[4]),float(p[5]),float(p[6]),float(p[7])
        Rd = quat_to_R(qw,qx,qy,qz); td = [x,y,z]
        # delta_s = E^-1 o delta_ref o E
        R1,t1 = se3_mul(ERt, Einv_t, Rd, td)     # E^-1 o delta
        Rs,ts = se3_mul(R1, t1, ER, Et)          # (E^-1 o delta) o E
        ts = [c*scale for c in ts]               # scale the translation magnitude
        ow,ox,oy,oz = R_to_quat(Rs)
        g.write("%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n" % (t+toff_ns, ts[0],ts[1],ts[2], ow,ox,oy,oz))
        n_out += 1
print("wrote %s  (%d increment rows)  E=yaw%.3f pitch%.3f roll%.3f lever%s scale=%.4f toff=%.3fs"
      % (out, n_out, yaw, pitch, roll, lever, scale, toff_ns/1e9))
