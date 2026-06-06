// ofc_sim/trajectory.hpp — parameterized ground-truth SE(3) trajectory (the oracle).
//
// RELAXED EDGE (sim/): lives under sim/, never ships to the ECU, so std containers
// / exceptions are fine. Depends only on the public ofc/core headers.
//
// Time -> pose model: the trajectory is a sequence of constant-twist SEGMENTS. Each
// segment has a body twist xi = [v; omega] (Vec6, the ofc convention — translation
// first) held constant for a duration. Within a segment the pose evolves by the
// exact constant-twist flow
//      pose(t) = pose(seg_start) o se3::exp(xi * tau)        (tau = t - seg_start)
// integrated with the in-house se3::exp, so the model matches exactly what the
// SourceBuffer / Estimator integrate with (a fair oracle). Segment boundary poses
// are precomputed cumulatively, so pose continuity across boundaries is exact by
// construction. The body twist twist(t) is piecewise constant (the active segment's
// xi); it is the analytic derivative of pose within a segment and matches a finite
// difference of pose(t) away from boundaries.
//
// Regime presets implement the observability spine (DESIGN §6): each calibration DOF
// is observable in exactly one motion regime, so the validation harness needs a
// trajectory for each regime in isolation plus a mixed driving scenario.
#ifndef OFC_SIM_TRAJECTORY_HPP
#define OFC_SIM_TRAJECTORY_HPP

#include "ofc/core/types.hpp"

#include <vector>

namespace ofc {
namespace sim {

// One constant-twist arc. `twist` is the body twist [v; omega] held constant for
// `duration_s` seconds. Examples: straight = [v,0,0, 0,0,0]; turn-in-place =
// [0,0,0, 0,0,wz]; curved (arc) = [v,0,0, 0,0,wz].
struct Segment {
    Vec6   twist    = Vec6::Zero();   // [v; omega] body twist (ofc convention)
    Scalar duration_s = Scalar(0);
};

class Trajectory {
public:
    // Empty trajectory (pose == identity everywhere, twist == zero). Build with
    // add_segment() or use a preset.
    Trajectory() = default;

    // The trajectory's time origin (seconds). pose(t<=t0_s) == start_pose; the first
    // segment begins at t0_s. Defaults to 0.
    explicit Trajectory(Scalar t0_s) : t0_s_(t0_s) {}

    // Append a constant-twist segment. Returns *this for fluent construction.
    Trajectory& add_segment(const Vec6& body_twist, Scalar duration_s);
    Trajectory& add_segment(const Segment& seg);

    // Set the pose at t0 (the trajectory start). Default identity.
    Trajectory& set_start(const SE3& start) { start_ = start; rebuild_(); return *this; }

    // --- Queries (the oracle) ------------------------------------------------
    // GT pose at time t (nanoseconds). Clamped: t before t0 -> start pose; t after the
    // last segment -> the final pose (extrapolation does NOT continue past the end).
    SE3 pose(Timestamp t) const;

    // GT body twist [v; omega] at time t (nanoseconds). Piecewise constant: the active
    // segment's twist. Before t0 or after the end -> zero (the trajectory is at rest).
    Vec6 twist(Timestamp t) const;

    // Convenience seconds overloads.
    SE3  pose_s(Scalar t_s) const;
    Vec6 twist_s(Scalar t_s) const;

    Scalar t0_s() const { return t0_s_; }
    Scalar duration_s() const { return total_s_; }      // sum of segment durations
    Scalar end_s() const { return t0_s_ + total_s_; }
    int    segment_count() const { return static_cast<int>(segs_.size()); }

    // --- Regime presets (DESIGN §6 observability spine) ----------------------
    // straight: ||omega|| == 0, ||v|| > 0 (recovers yaw/pitch forward axis + scale).
    static Trajectory straight(Scalar v = Scalar(2.0), Scalar duration_s = Scalar(5.0));

    // turning: ||omega|| > 0 (recovers roll + xyz lever-arm via hand-eye). Forward
    // motion + a steady yaw rate (a curved arc) by default.
    static Trajectory turning(Scalar v = Scalar(2.0), Scalar wz = Scalar(0.5),
                              Scalar duration_s = Scalar(5.0));

    // omega_varying: ||omega|| changes over time (a distinctive shape for time-sync
    // cross-correlation). Alternating turn rates separated by straight gaps.
    static Trajectory omega_varying(Scalar v = Scalar(2.0), Scalar wz = Scalar(0.6),
                                    Scalar seg_s = Scalar(1.0));

    // mixed: a representative driving scenario exercising every regime — straight,
    // curved arcs (both turn directions), a turn-in-place, and pitch (so all DOF see
    // excitation somewhere). The default scenario the rig end-to-end test drives.
    static Trajectory mixed();

private:
    void rebuild_();   // recompute cumulative segment-boundary poses + times
    int  segment_at_(Scalar t_s, Scalar& tau_local) const;  // active seg + local time

    Scalar              t0_s_   = Scalar(0);
    Scalar              total_s_ = Scalar(0);
    SE3                 start_;                 // pose at t0
    std::vector<Segment> segs_;
    // Precomputed boundary state (size = segs_.size() + 1):
    //   cum_t_[k]    = cumulative seconds at the START of segment k (relative to t0)
    //   cum_pose_[k] = GT pose at the start of segment k
    std::vector<Scalar> cum_t_;
    std::vector<SE3>    cum_pose_;
};

} // namespace sim
} // namespace ofc
#endif // OFC_SIM_TRAJECTORY_HPP
