// ofc_sim/trajectory.cpp — see trajectory.hpp for the time->pose model.
#include "ofc_sim/trajectory.hpp"

#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {
namespace sim {

namespace {
constexpr Scalar kNanosPerSec = Scalar(1e9);
Scalar ns_to_s(Timestamp t) { return static_cast<Scalar>(t) / kNanosPerSec; }
} // namespace

Trajectory& Trajectory::add_segment(const Vec6& body_twist, Scalar duration_s) {
    Segment s;
    s.twist      = body_twist;
    s.duration_s = duration_s;
    return add_segment(s);
}

Trajectory& Trajectory::add_segment(const Segment& seg) {
    segs_.push_back(seg);
    rebuild_();
    return *this;
}

void Trajectory::rebuild_() {
    total_s_ = Scalar(0);
    cum_t_.clear();
    cum_pose_.clear();
    cum_t_.reserve(segs_.size() + 1);
    cum_pose_.reserve(segs_.size() + 1);

    SE3    pose = start_;
    Scalar t    = Scalar(0);
    cum_t_.push_back(t);
    cum_pose_.push_back(pose);
    for (const Segment& s : segs_) {
        // pose advances by the constant-twist flow over the whole segment.
        pose = se3::compose(pose, se3::exp(s.twist * s.duration_s));
        t   += s.duration_s;
        cum_t_.push_back(t);
        cum_pose_.push_back(pose);
    }
    total_s_ = t;
}

// Find the active segment for trajectory-relative seconds `rel` (rel in [0, total]).
// Returns the segment index and sets tau_local to the time into that segment. For
// rel at/after the end returns the last segment with tau_local == its duration (so the
// integrated pose equals the final boundary pose); for an empty trajectory returns -1.
int Trajectory::segment_at_(Scalar rel, Scalar& tau_local) const {
    tau_local = Scalar(0);
    if (segs_.empty()) return -1;
    if (rel <= Scalar(0)) { tau_local = Scalar(0); return 0; }
    if (rel >= total_s_) {
        const int last = static_cast<int>(segs_.size()) - 1;
        tau_local = segs_[last].duration_s;
        return last;
    }
    // cum_t_[k] is the start time of segment k. Find the segment whose [start,end)
    // contains rel. Linear scan (sim is off-target; trajectories are short).
    for (int k = 0; k < static_cast<int>(segs_.size()); ++k) {
        const Scalar seg_start = cum_t_[k];
        const Scalar seg_end   = cum_t_[k + 1];
        if (rel < seg_end) {
            tau_local = rel - seg_start;
            return k;
        }
    }
    const int last = static_cast<int>(segs_.size()) - 1;
    tau_local = segs_[last].duration_s;
    return last;
}

SE3 Trajectory::pose(Timestamp t) const {
    return pose_s(ns_to_s(t));
}

SE3 Trajectory::pose_s(Scalar t_s) const {
    if (segs_.empty()) return start_;
    const Scalar rel = t_s - t0_s_;
    Scalar tau = Scalar(0);
    const int k = segment_at_(rel, tau);
    if (k < 0) return start_;
    // pose = boundary pose at the start of segment k, advanced by tau under its twist.
    return se3::compose(cum_pose_[k], se3::exp(segs_[k].twist * tau));
}

Vec6 Trajectory::twist(Timestamp t) const {
    return twist_s(ns_to_s(t));
}

Vec6 Trajectory::twist_s(Scalar t_s) const {
    if (segs_.empty()) return Vec6::Zero();
    const Scalar rel = t_s - t0_s_;
    // At rest outside the trajectory's span. HALF-OPEN convention: the span is [0, total),
    // so twist drops to zero AT exactly rel == total (the `>=`), one tick before pose_s
    // stops moving (pose_s clamps to the final boundary pose only for rel > total via
    // segment_at_'s `rel >= total_s_` branch, which still returns the last segment's twist
    // for the pose integral). Net: twist(end) == 0 while pose(end) is still the integrand's
    // endpoint — harmless for queries, but twist(end) is NOT the instantaneous integrand of
    // pose near end. Callers needing the end-of-motion twist should query just before end.
    if (rel < Scalar(0) || rel >= total_s_) return Vec6::Zero();
    Scalar tau = Scalar(0);
    const int k = segment_at_(rel, tau);
    if (k < 0) return Vec6::Zero();
    return segs_[k].twist;
}

// --- Presets ---------------------------------------------------------------

Trajectory Trajectory::straight(Scalar v, Scalar duration_s) {
    Trajectory tr;
    Vec6 xi; xi << v, Scalar(0), Scalar(0), Scalar(0), Scalar(0), Scalar(0);
    tr.add_segment(xi, duration_s);
    return tr;
}

Trajectory Trajectory::turning(Scalar v, Scalar wz, Scalar duration_s) {
    Trajectory tr;
    Vec6 xi; xi << v, Scalar(0), Scalar(0), Scalar(0), Scalar(0), wz;
    tr.add_segment(xi, duration_s);
    return tr;
}

Trajectory Trajectory::omega_varying(Scalar v, Scalar wz, Scalar seg_s) {
    // Alternating turn rates separated by straight gaps -> ||omega||(t) has a
    // distinctive shape for time-sync cross-correlation (DESIGN §6 time-offset row).
    Trajectory tr;
    Vec6 straight;  straight << v, 0, 0, 0, 0, 0;
    Vec6 left;      left     << v, 0, 0, 0, 0,  wz;
    Vec6 right;     right    << v, 0, 0, 0, 0, -wz;
    Vec6 hard;      hard     << v, 0, 0, 0, 0,  Scalar(2) * wz;
    tr.add_segment(straight, seg_s);
    tr.add_segment(left,     seg_s);
    tr.add_segment(straight, seg_s * Scalar(0.5));
    tr.add_segment(hard,     seg_s);
    tr.add_segment(right,    seg_s);
    tr.add_segment(straight, seg_s);
    return tr;
}

Trajectory Trajectory::mixed() {
    // A representative driving scenario: every regime gets excitation somewhere.
    Trajectory tr;
    Vec6 straight;  straight << Scalar(2.0), 0, 0, 0, 0, 0;                 // cruise
    Vec6 arc_l;     arc_l    << Scalar(2.0), 0, 0, 0, 0, Scalar(0.5);       // left arc
    Vec6 arc_r;     arc_r    << Scalar(1.5), 0, 0, 0, 0, Scalar(-0.6);      // right arc
    Vec6 spin;      spin     << 0, 0, 0, 0, 0, Scalar(0.8);                 // turn in place
    Vec6 climb;     climb    << Scalar(1.5), 0, Scalar(0.3), 0, Scalar(0.2), 0; // pitch + rise
    tr.add_segment(straight, Scalar(1.0));
    tr.add_segment(arc_l,    Scalar(1.5));
    tr.add_segment(straight, Scalar(0.8));
    tr.add_segment(arc_r,    Scalar(1.2));
    tr.add_segment(spin,     Scalar(0.7));
    tr.add_segment(climb,    Scalar(1.0));
    tr.add_segment(straight, Scalar(1.0));
    return tr;
}

} // namespace sim
} // namespace ofc
