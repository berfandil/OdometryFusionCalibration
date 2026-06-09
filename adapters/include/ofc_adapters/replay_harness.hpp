// ofc_adapters/replay_harness.hpp — OFFLINE REPLAY + EVAL harness (Slice 13, real-data
// scaffolding). The real-data analogue of the sim Rig (sim/src/rig.cpp): it pumps
// ofc::Estimator::step() over a timeline merged from a set of CsvSources (+ an optional GPS
// stream + an optional GT track), records the per-published-step output, and — with a GT track —
// computes drift / 6-DOF pose NEES / accepted-GPS NIS exactly as the sim validation harness does
// (tests/test_validation.cpp pose_nees + the NIS pattern).
//
// RELAXED EDGE. Depends only on public ofc/core headers + the adapter's CsvSource / GpsCorrection.
//
// TIMELINE: the harness merges every source stamp + every GPS-fix stamp into a sorted, de-duped
// tick timeline, then steps the estimator at each tick. (Stepping at the source stamps — not a
// synthetic fixed rate — keeps real, irregular sensor cadences faithful; the estimator's frontier
// is `now - fusion_delay`, so a tick must be at least one fusion_delay past the data start before
// a fuse succeeds.) A GPS fix is submitted to the GpsCorrection just before the step whose
// frontier first reaches the fix's stamp, so step() consumes it at the matching frontier.
//
// GT ANCHORING: the estimator anchors its odom pose at the GT pose at (first_frontier - window_s)
// (DESIGN §7, gauge anchored at init), the SAME convention the sim Rig uses. The harness records
// the GT pose expressed RELATIVE to that origin so the frontier-vs-GT error compares like-for-like.
//
// EVAL (per published step, WITH a GT track):
//   trans_err_m  = ||frontier.t - gt.t||
//   rot_err_rad  = ||so3::log(gt.R^T frontier.R)||
//   pose_nees    = e^T (P_pose)^-1 e,  e = se3::log(T_est^-1 o T_gt), P_pose = frontier.cov(0..5)
//                  (the EXACT convention of tests/test_validation.cpp pose_nees).
//   nis          = CorrectionDiag.last_nis on steps where a GPS fix was APPLIED (corr_applied>0).
// AGGREGATE (RunSummary): max/RMS/tail drift, ensemble-mean pose NEES, mean accepted-GPS NIS.
#ifndef OFC_ADAPTERS_REPLAY_HARNESS_HPP
#define OFC_ADAPTERS_REPLAY_HARNESS_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include "ofc_adapters/csv_source.hpp"
#include "ofc_adapters/gps_correction.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace ofc {
namespace adapters {

// One recorded published step.
struct ReplayRecord {
    Timestamp now       = 0;     // wall-clock 'now' passed to step()
    bool      fused     = false; // step() returned Ok
    Result    result{};          // the fused Result this step (frontier/tip/calib/correction)

    // GT eval (valid only when has_gt). gt_pose is in the estimator's odom frame (relative to the
    // anchor at first_frontier - window_s).
    bool   has_gt        = false;
    SE3    gt_pose;
    Scalar trans_err_m   = Scalar(0);
    Scalar rot_err_rad   = Scalar(0);
    Scalar pose_nees     = Scalar(0);   // 6-DOF, valid when has_gt
};

// Aggregate verdict over a run (the summary the user reads to judge real-data consistency).
struct RunSummary {
    int    steps          = 0;   // ticks stepped
    int    fused_steps    = 0;   // steps that fused
    int    gps_applied    = 0;   // accepted GPS fixes (corr_applied accumulated)
    int    gps_rejected   = 0;   // Mahalanobis-rejected GPS fixes
    int    gps_evaluated  = 0;   // fixes evaluated

    // Drift (only meaningful with a GT track; 0 otherwise).
    bool   has_gt         = false;
    Scalar max_trans_m    = Scalar(0);
    Scalar max_rot_rad    = Scalar(0);
    Scalar rms_trans_m    = Scalar(0);
    Scalar rms_rot_rad    = Scalar(0);
    Scalar tail_trans_m   = Scalar(0);  // mean translation error over the tail window
    Scalar tail_rot_rad   = Scalar(0);

    // Consistency (only when the relevant samples exist).
    long   nees_count     = 0;
    Scalar mean_pose_nees = Scalar(0);  // ensemble-mean 6-DOF pose NEES (target ~ DOF=6)
    long   nis_count      = 0;
    Scalar mean_nis       = Scalar(0);  // mean accepted-GPS NIS (target ~ DOF=3 position fix)

    // Local (GT-anchored fixed-window) relative-pose error — the KITTI-style segment metric. Each
    // window of `local_batch_len` consecutive fused-with-GT records is re-anchored to GT at its
    // start (E = rel_gt^-1 o rel_est), so ONLY intra-window drift is measured. This makes drift
    // comparable across recordings of different total length (the global tail grows with run
    // length; the local error does not). Populated only when local_batch_len > 0 and enough
    // post-warmup fused-with-GT records exist; counts are 0 otherwise.
    int    local_batch_len    = 0;          // window length (num fused-with-GT records); 0 = off
    int    local_batch_count  = 0;          // complete windows measured
    int    local_dropped      = 0;          // trailing fused-gt records not forming a full window
    Scalar local_mean_trans_m = Scalar(0);
    Scalar local_med_trans_m  = Scalar(0);  // median window translation error (outlier-robust)
    Scalar local_max_trans_m  = Scalar(0);  // worst window (localizes a corruption event)
    Scalar local_mean_rot_rad = Scalar(0);
    Scalar local_med_rot_rad  = Scalar(0);
    Scalar local_max_rot_rad  = Scalar(0);
};

// Inputs to a replay run. The Config + the CsvSources + the optional GPS + GT are all caller-owned
// and must outlive run_replay(). The Config's sensor storage follows the usual outlive contract.
struct ReplayInputs {
    const Config*               cfg = nullptr;     // validated config (init's the estimator)
    std::vector<CsvSource*>     sources;           // odometry sources (>= 1)
    GpsCorrection*              gps = nullptr;      // optional absolute ref
    std::vector<GpsFix>         gps_fixes;          // fixes to submit (sorted by stamp); empty -> none
    const CsvGtTrack*           gt  = nullptr;      // optional GT track (eval only)

    // Tail window (seconds, trailing the run) for the tail-mean drift metric. Default 1 s.
    Scalar tail_window_s = Scalar(1.0);
    // Warmup steps to skip before accumulating NEES (the calibration/cov transient). Default 20,
    // matching the sim validation harness. The local-window metric uses the SAME warmup gate.
    int    warmup_steps  = 20;
    // Local relative-pose-error window length (num post-warmup fused-with-GT records per GT-anchored
    // window). 0 disables the local metric. Same value across recordings -> drift comparable across
    // runs of different total length (each window re-anchors to GT, so only intra-window drift
    // counts). See RunSummary::local_* and write_summary().
    int    local_batch_len = 0;
};

// The harness. run() owns the Estimator; it is re-init'd on each run().
class ReplayHarness {
public:
    ReplayHarness() = default;

    // Pump the estimator over the merged timeline, record each published step, and compute the
    // GT eval + the aggregate summary. Returns the run status:
    //   Ok            — at least one step fused.
    //   InvalidConfig — null cfg / no sources / estimator init or registration failed.
    //   NoData        — no usable timeline (no source stamps) or zero fuses.
    Status run(const ReplayInputs& in);

    const std::vector<ReplayRecord>& records() const { return records_; }
    const RunSummary&                summary() const { return summary_; }
    const std::string&               error()   const { return error_; }

    // Write the per-step results CSV (one row per FUSED step) + a trailing `# summary:` block.
    // Columns:
    //   now_ns, frontier_ns, phase, readiness,
    //   x,y,z, qw,qx,qy,qz,                              (frontier pose, t + w-first quat)
    //   p00..p05,                                        (frontier pose-cov DIAGONAL, [trans;rot])
    //   tip_valid, tip_x,tip_y,tip_z,
    //   corr_evaluated, corr_applied, corr_rejected, last_nis,
    //   has_gt, gt_x,gt_y,gt_z, trans_err_m, rot_err_rad, pose_nees,
    //   src0_id,src0_scale,src0_off_s,src0_conf, ... (one block per source, up to source_count)
    // Returns Ok or InvalidConfig (no records).
    Status write_results_csv(std::ostream& os) const;
    Status write_results_csv_file(const std::string& path) const;

    // Write just the human-readable summary block (also used by the CLI's stdout report).
    void write_summary(std::ostream& os) const;

private:
    Estimator                 est_;
    std::vector<ReplayRecord>  records_;
    RunSummary                 summary_{};
    std::string                error_;
};

} // namespace adapters
} // namespace ofc
#endif // OFC_ADAPTERS_REPLAY_HARNESS_HPP
