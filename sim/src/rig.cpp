// ofc_sim/rig.cpp — see rig.hpp.
#include "ofc_sim/rig.hpp"

#include "ofc/core/lie.hpp"

#include <cmath>

namespace ofc {
namespace sim {

namespace {
constexpr Scalar kNanosPerSec = Scalar(1e9);
Timestamp secs_to_ns(Scalar s) {
    return static_cast<Timestamp>(std::llround(s * kNanosPerSec));
}
} // namespace

Status Rig::init(const Config& cfg) {
    records_.clear();
    sources_.clear();
    window_s_ = cfg.window_s;     // bootstrap/lookback of the first fuse (DESIGN §7)
    return est_.init(cfg);
}

Status Rig::add_source(const SyntheticSource* src) {
    if (src == nullptr) return Status::InvalidConfig;
    const Status s = est_.add_source(src);
    if (ok(s)) sources_.push_back(src);
    return s;
}

Status Rig::add_correction(const ICorrection* corr) {
    return est_.add_correction(corr);
}

int Rig::run(Scalar from_s, Scalar to_s, Scalar tick_rate_hz) {
    records_.clear();
    if (traj_ == nullptr || tick_rate_hz <= Scalar(0) || to_s <= from_s) return 0;

    const Timestamp tick = secs_to_ns(Scalar(1) / tick_rate_hz);
    const Timestamp t0   = secs_to_ns(from_s);
    const Timestamp tend = secs_to_ns(to_s);

    // The estimator anchors its odom pose at IDENTITY before the first predict, then
    // composes the first window's delta (DESIGN §7, gauge anchored at init). So the
    // pose published at the first fused frontier is the integrated motion over the
    // bootstrap window [first_frontier - window_s, first_frontier], NOT identity —
    // i.e. the estimator's pose origin sits at the GT pose at (first_frontier -
    // window_s). We anchor GT there so pose_error() compares like-for-like: both the
    // fused pose and the recorded GT are expressed relative to that origin and the
    // estimator integrates gap-free from it onward.
    const Timestamp window_ns = secs_to_ns(window_s_);
    bool      have_anchor = false;
    SE3       gt_anchor_inv;       // inverse of the GT pose at the estimator's origin

    int fuses = 0;
    for (Timestamp now = t0; now <= tend; now += tick) {
        const Status st = est_.step(now);
        Record rec;
        rec.now    = now;
        rec.fused  = ok(st);
        if (ok(st)) {
            rec.result = est_.latest();
            // GT base pose at the fused frontier (the published frontier stamp).
            const Timestamp frontier = rec.result.frontier.stamp;
            const SE3 gt_abs = traj_->pose(frontier);
            if (!have_anchor) {
                // Origin = GT pose at the START of the first bootstrap window.
                gt_anchor_inv = se3::inverse(traj_->pose(frontier - window_ns));
                have_anchor   = true;
            }
            // Express GT in the estimator's odom frame (anchored at that origin).
            rec.gt_frontier = se3::compose(gt_anchor_inv, gt_abs);
            ++fuses;
        }
        records_.push_back(rec);
    }
    return fuses;
}

void Rig::pose_error(const Record& r, Scalar& trans_err_m, Scalar& rot_err_rad) {
    const SE3& fused = r.result.frontier.pose;
    const SE3& gt    = r.gt_frontier;
    trans_err_m = (fused.t - gt.t).norm();
    rot_err_rad = so3::log(gt.R.transpose() * fused.R).norm();
}

void Rig::max_error(Scalar& max_trans_m, Scalar& max_rot_rad) const {
    max_trans_m = Scalar(0);
    max_rot_rad = Scalar(0);
    for (const Record& r : records_) {
        if (!r.fused) continue;
        Scalar te, re;
        pose_error(r, te, re);
        if (te > max_trans_m) max_trans_m = te;
        if (re > max_rot_rad) max_rot_rad = re;
    }
}

} // namespace sim
} // namespace ofc
