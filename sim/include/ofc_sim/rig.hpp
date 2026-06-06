// ofc_sim/rig.hpp — drives the Estimator over a GT trajectory and records the run.
//
// RELAXED EDGE (sim/). The rig owns a Trajectory + a set of SyntheticSources + an
// ofc::Estimator. run() steps the estimator across a time range at a fixed tick rate
// and records, per successful tick, { now, the fused Result, the GT pose at the fused
// frontier }. Tests assert the fused frontier tracks GT and that replay is byte-stable.
//
// The GT comparison is anchored the same way the estimator anchors its odom frame: the
// estimator starts its pose at identity at the FIRST fused tick (DESIGN §7, gauge
// anchored at init), so the recorded GT pose is expressed RELATIVE to the GT pose at
// that first fused frontier — pose_error() then compares like-for-like.
#ifndef OFC_SIM_RIG_HPP
#define OFC_SIM_RIG_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/types.hpp"

#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <vector>

namespace ofc {
namespace sim {

// One recorded tick.
struct Record {
    Timestamp now      = 0;     // wall-clock 'now' passed to step()
    Result    result{};         // the fused Result at this tick
    SE3       gt_frontier;      // GT base pose at the fused frontier, in the odom frame
                                // (relative to the first fused frontier's GT pose)
    bool      fused = false;    // step() returned Ok (a fuse happened)
};

class Rig {
public:
    Rig() = default;

    // Bind the trajectory used as the GT oracle (must outlive the rig).
    void set_trajectory(const Trajectory& traj) { traj_ = &traj; }

    // Configure + init the estimator. The SensorConfig array is owned by the caller
    // and must outlive the rig (mirrors the core's config contract). Returns the
    // Estimator's init status.
    Status init(const Config& cfg);

    // Register a synthetic source (pointer owned by the caller; must outlive the rig).
    Status add_source(const SyntheticSource* src);

    // Step from `from_s` to `to_s` (trajectory seconds) at `tick_rate_hz`. Records each
    // tick. Clears any previous recording. Returns the number of SUCCESSFUL fuses.
    int run(Scalar from_s, Scalar to_s, Scalar tick_rate_hz);

    const std::vector<Record>& records() const { return records_; }
    const Estimator&           estimator() const { return est_; }

    // Fused-vs-GT pose error at a recorded tick: split into translation (m) and
    // rotation (rad) magnitude. Both are zero when the fused frontier matches GT.
    static void pose_error(const Record& r, Scalar& trans_err_m, Scalar& rot_err_rad);

    // Max translation / rotation error over all FUSED records (0 if none fused).
    void max_error(Scalar& max_trans_m, Scalar& max_rot_rad) const;

private:
    const Trajectory*                    traj_ = nullptr;
    Estimator                            est_;
    std::vector<const SyntheticSource*>  sources_;
    std::vector<Record>                  records_;
    Scalar                               window_s_ = Scalar(0.10);  // from Config (GT anchor)
};

} // namespace sim
} // namespace ofc
#endif // OFC_SIM_RIG_HPP
