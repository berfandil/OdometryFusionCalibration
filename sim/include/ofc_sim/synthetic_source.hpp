// ofc_sim/synthetic_source.hpp — a synthetic ISource with PLANTED calibration error.
//
// RELAXED EDGE (sim/): std containers / std::mt19937 / exceptions are fine here.
//
// Measurement model (the inverse of the estimator's frame-align, D21). Over a query
// window [t0, t1] the source:
//   1. shifts the window by its planted time offset:   ts = t0 + off, te = t1 + off
//   2. reads the TRUE base motion from the trajectory:
//          A = pose(ts)^{-1} o pose(te)                 (base-frame delta)
//      (left-invariant body delta: the motion expressed at the body frame at ts).
//   3. reports it in its OWN sensor frame by undoing the mount:
//          B = X^{-1} o A o X                            (hand-eye inverse)
//      so that the estimator's frame-align  A' = X o B o X^{-1}  recovers A exactly
//      (DECISIONS D21: the estimator applies A = X o B o X^{-1}, X = sensor->base).
//   4. scales the translation part of B by `scale` (per-source translation scale, D20).
//   5. adds zero-mean Gaussian noise (seeded, deterministic) in the body tangent.
//   6. attaches a covariance consistent with the injected noise.
//
// Identity X, scale 1, no noise  =>  B == A == the GT base delta (verified in tests).
//
// Injection: an outlier window replaces B with a gross-wrong delta; a dropout window
// makes query() return Status::NoData (a sensor gap).
//
// Determinism: the noise PRNG is a std::mt19937 seeded once from `seed` (a fixed seed
// per source). Because query() is const but must advance the PRNG, the generator is
// mutable. To keep the stream reproducible *independent of the call sequence*, each
// window draws from a sub-generator seeded by hashing (seed, t0, t1) — so the noise on
// a given window is a pure function of that window, not of how many queries preceded
// it. This makes the rig deterministic regardless of tick alignment.
#ifndef OFC_SIM_SYNTHETIC_SOURCE_HPP
#define OFC_SIM_SYNTHETIC_SOURCE_HPP

#include "ofc/core/buffer.hpp"
#include "ofc/core/source.hpp"
#include "ofc/core/types.hpp"

#include "ofc_sim/trajectory.hpp"

#include <vector>

namespace ofc {
namespace sim {

// A half-open time window [start_s, end_s) for injecting faults, in trajectory time.
struct Window {
    Scalar start_s = Scalar(0);
    Scalar end_s   = Scalar(0);
    bool contains(Scalar t_s) const { return t_s >= start_s && t_s < end_s; }
};

// Planted parameters for one synthetic source.
struct SourceParams {
    SourceId id = 0;

    SE3    X;                       // planted extrinsic (sensor->base). Identity = aligned.
    Scalar scale         = Scalar(1.0);   // planted translation scale (D20).
    Scalar time_offset_s = Scalar(0.0);   // planted clock offset (this source leads/lags).

    // Noise model (deterministic, seeded). Per-distance translational sigma and
    // per-angle rotational sigma, applied in the body tangent of the reported delta.
    Scalar noise_trans_per_m   = Scalar(0.0);  // [m sigma per m travelled]
    Scalar noise_rot_per_rad   = Scalar(0.0);  // [rad sigma per rad rotated]
    Scalar noise_trans_floor   = Scalar(0.0);  // [m] minimum translational sigma
    Scalar noise_rot_floor     = Scalar(0.0);  // [rad] minimum rotational sigma
    std::uint32_t seed         = 0;            // fixed per-source seed (reproducible)

    // Whether to report a native covariance (Delta.cov). When true the cov is
    // synthesized from the noise model; when false an identity cov is reported.
    bool report_native_cov = true;

    // Injection windows (trajectory time). Empty by default.
    std::vector<Window> outlier_windows;       // query() returns a gross-wrong delta
    std::vector<Window> dropout_windows;       // query() returns Status::NoData

    // The gross-wrong delta applied during an outlier window (a fixed body twist
    // integrated over the window). Defaults to a large lateral + yaw error.
    Vec6 outlier_twist = (Vec6() << Scalar(5.0), Scalar(3.0), Scalar(0.0),
                                    Scalar(0.0), Scalar(0.0), Scalar(2.0)).finished();
};

class SyntheticSource : public ISource {
public:
    SyntheticSource(const Trajectory& traj, const SourceParams& p);

    SourceId id() const override { return p_.id; }

    // ANALYTIC delta (priority path): the measurement model above, computed directly
    // from the trajectory. NoData on a dropout window (or a non-positive window).
    Expected<Delta> query(Timestamp t0, Timestamp t1) const override;

    const SourceParams& params() const { return p_; }

    // --- Buffer-backed mode (optional, exercises the SourceBuffer path) -------
    // Populate an internal SourceBuffer by sampling the source's reported motion as
    // INCREMENTS at `rate_hz` over [from_s, to_s] (trajectory time, BEFORE the time
    // offset — the offset is baked into the samples' stamps so query stamps line up
    // with the analytic path). After this, use query_buffered() to answer from the
    // buffer instead of analytically. Returns Ok or a buffer status.
    Status build_buffer(Scalar from_s, Scalar to_s, Scalar rate_hz,
                        OdomForm form = OdomForm::Increment, int capacity = 0);

    // Answer from the internal buffer (must have called build_buffer first). Falls back
    // to NoData if not built or the window is uncovered.
    Expected<Delta> query_buffered(Timestamp t0, Timestamp t1) const;

    // The analytic TRUE base delta A over [t0+off, t1+off] (no X, no scale, no noise) —
    // exposed for tests / the rig's GT comparisons.
    SE3 base_delta(Timestamp t0, Timestamp t1) const;

private:
    bool in_any_(const std::vector<Window>& ws, Timestamp t0, Timestamp t1) const;
    // The clean reported sensor-frame delta B = scale_t( X^{-1} o A o X ) (no noise).
    SE3  clean_reported_(Timestamp t0, Timestamp t1) const;
    // Synthesize the reported covariance from the noise model + motion magnitude.
    Mat6 modeled_cov_(const SE3& reported) const;

    const Trajectory* traj_;
    SourceParams      p_;

    SourceBuffer      buf_;
    bool              buffered_ = false;
};

} // namespace sim
} // namespace ofc
#endif // OFC_SIM_SYNTHETIC_SOURCE_HPP
