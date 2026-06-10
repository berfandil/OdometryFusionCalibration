# Convert a EuRoC MAV ASL-format sequence to our replay CSVs (Slice 13 #1 / EuRoC).
#
# EuRoC is a vis-inertial DRONE dataset: full 6-DOF aggressive 3D motion + sub-mm 6-DOF GT, but NO
# wheel/GPS odometry. Its value for us is the 3-D excitation that planar KAIST lacked -> it can test
# the FULL extrinsic (yaw/pitch/roll + xyz lever) that KAIST left weakly observable. So we derive a
# CLEAN increment-form odometry from the GT trajectory (per-step pose increments) -> a "perfect"
# source following real 3-D motion. The calibration test then makes 3 sources from it (ref, clean,
# inject_calib'd) and checks the calibrator recovers the injected transform on real 3-D motion.
#
# Reads the nested zip IN MEMORY (vicon_room1.zip -> V1_01_easy.zip -> mav0/...), pulling only the
# GT csv (~5 MB) -- never unzips, never touches the camera images, leaves the source zip untouched.
#
# Usage: python euroc_to_csv.py  C:/workspace/data/EuRoC/vicon_room1.zip  V1_01_easy  euroc_run
import zipfile, io, os, sys, math

src_zip = sys.argv[1] if len(sys.argv) > 1 else r"C:/workspace/data/EuRoC/vicon_room1.zip"
seq     = sys.argv[2] if len(sys.argv) > 2 else "V1_01_easy"
outdir  = sys.argv[3] if len(sys.argv) > 3 else "euroc_run"
os.makedirs(outdir, exist_ok=True)

# --- locate + read the GT csv from the nested zip (in memory) -------------------------------
outer = zipfile.ZipFile(src_zip)
inner_name = next(e for e in outer.namelist() if e.endswith("%s.zip" % seq))
inner = zipfile.ZipFile(io.BytesIO(outer.read(inner_name)))
gt_entry = next(e for e in inner.namelist()
                if e.endswith("state_groundtruth_estimate0/data.csv") and not e.startswith("__"))
rows = inner.read(gt_entry).decode().splitlines()

# --- minimal SO3/SE3 (lists) ----------------------------------------------------------------
def quat_to_R(qw,qx,qy,qz):
    n=math.sqrt(qw*qw+qx*qx+qy*qy+qz*qz) or 1.0; qw,qx,qy,qz=qw/n,qx/n,qy/n,qz/n
    return [[1-2*(qy*qy+qz*qz),2*(qx*qy-qz*qw),2*(qx*qz+qy*qw)],
            [2*(qx*qy+qz*qw),1-2*(qx*qx+qz*qz),2*(qy*qz-qx*qw)],
            [2*(qx*qz-qy*qw),2*(qy*qz+qx*qw),1-2*(qx*qx+qy*qy)]]
def R_to_quat(R):
    tr=R[0][0]+R[1][1]+R[2][2]
    if tr>0:
        s=math.sqrt(tr+1.0)*2; qw=0.25*s
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
def mT(A): return [[A[j][i] for j in range(3)] for i in range(3)]
def mv(A,v): return [sum(A[i][k]*v[k] for k in range(3)) for i in range(3)]
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def se3_inv(R,t): Ri=mT(R); return Ri,[-x for x in mv(Ri,t)]
def se3_mul(R1,t1,R2,t2): return mm(R1,R2),[mv(R1,t2)[i]+t1[i] for i in range(3)]

# --- parse GT poses -------------------------------------------------------------------------
poses=[]
for ln in rows:
    ln=ln.strip()
    if not ln or ln.startswith("#"): continue
    p=ln.split(",")
    t=int(p[0]); px,py,pz=float(p[1]),float(p[2]),float(p[3])
    qw,qx,qy,qz=float(p[4]),float(p[5]),float(p[6]),float(p[7])
    poses.append((t,quat_to_R(qw,qx,qy,qz),[px,py,pz]))
print("GT poses: %d  span %.1fs" % (len(poses), (poses[-1][0]-poses[0][0])/1e9))

# --- write GT track (absolute; harness anchors at first frontier) ---------------------------
gt_path=os.path.join(outdir, "%s_gt.csv"%seq)
with open(gt_path,"w",newline="") as g:
    g.write("# EuRoC GT (state_groundtruth_estimate0)  t_ns,x,y,z,qw,qx,qy,qz\n")
    for t,R,tt in poses:
        qw,qx,qy,qz=R_to_quat(R)
        g.write("%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n"%(t,tt[0],tt[1],tt[2],qw,qx,qy,qz))

# --- write CLEAN increment-form odometry derived from GT pose increments --------------------
# increment[k] = pose[k-1]^-1 o pose[k]  (the true body motion between GT samples).
odom_path=os.path.join(outdir, "%s_odom.csv"%seq)
with open(odom_path,"w",newline="") as o:
    o.write("# form: increment  (GT-derived clean odometry)\n# t_ns,x,y,z,qw,qx,qy,qz\n")
    o.write("%d,0,0,0,1,0,0,0\n"%poses[0][0])          # identity seed at t0
    for k in range(1,len(poses)):
        Rp,tp=poses[k-1][1],poses[k-1][2]
        Rc,tc=poses[k][1],poses[k][2]
        Ri,ti=se3_inv(Rp,tp)
        Rd,td=se3_mul(Ri,ti,Rc,tc)
        qw,qx,qy,qz=R_to_quat(Rd)
        o.write("%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n"%(poses[k][0],td[0],td[1],td[2],qw,qx,qy,qz))
print("wrote %s  and  %s"%(gt_path,odom_path))
