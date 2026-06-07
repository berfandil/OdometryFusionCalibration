// ofc_sim/synthetic_source.hpp — a synthetic ISource with PLANTED calibration error.
//
// RELAXED EDGE (sim/): std containers / std::mt19937 / exceptions are fine here.
//
// Measurement model (the inverse of the estimator's frame-align, D21). Over a query
// window [t0, t1] the source:
//   1. shifts the SAMPLED window LATER by its planted time offset (canonical sign —
//      positive off => the source clock is ahead of base time):
//          ts = t0 + off, te = t1 + off       (reads a later trajectory slice for off>0)
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
// Determinism: each window draws from a std::mt19937_64 seeded by hashing (seed, t0, t1)
// — so the noise on a given window is a pure function of that window, not of how many
// queries preceded it. This makes the rig deterministic regardless of tick alignment.
// The six body-tangent noise components are drawn into locals in explicit statement order
// before assembling the Vec6, so the draw order is fixed across C++14 toolchains (an Eigen
// comma-initializer would leave it unspecified pre-C++17).
//
// GOLDEN-REGRESSION CAVEAT (DESIGN §10): byte-stable golden output requires a FIXED stdlib.
// std::mt19937_64 is portable, but std::normal_distribution's transform of the engine
// stream is NOT specified by the standard — different stdlib implementations (libstdc++ /
// libc++ / MSVC) emit different normal variates from the same engine. So replay/golden
// regression is bit-exact only on a single stdlib; cross-stdlib comparison must use a
// tolerance, not byte equality.
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
    // Planted clock offset. CANONICAL SIGN CONVENTION (the contract Slice-5 time-sync
    // must invert): POSITIVE time_offset_s shifts the SAMPLED window LATER. A query for
    // the window [t0, t1] reads the TRUE trajectory over [t0 + off, t1 + off] — i.e. the
    // source's clock is AHEAD of base time by `off`, so at base time t it reports the
    // motion the base body actually had at t + off. Negative off reads an earlier slice.
    Scalar time_offset_s = Scalar(0.0);

    // Planted CONSTANT body-twist bias (Slice 11b, Option A). A rate b = [v_bias; omega_bias]
    // in R^6 added to the reported body motion over the window: the reported delta becomes
    // B <- B o exp(b * dt) (dt = window length in seconds). This is exactly the raw-IMU-style
    // rate offset the augmented ESKF de-biases with Delta o exp(-b*dt); planting it here lets a
    // bias-states source + an absolute ref recover it. Zero = no bias (every other source). The
    // bias is applied BEFORE the per-window Gaussian noise and is NOT injected on outlier windows
    // (the outlier already replaces the whole delta).
    Vec6 body_twist_bias = Vec6::Zero();

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
    // Synthesize the reported covariance from the noise model + motion magnitude of the
    // CLEAN pre-noise delta (so the reported Sigma == the injected sigma; on an outlier
    // window pass the clean window delta, not the gross outlier delta).
    Mat6 modeled_cov_(const SE3& clean_delta) const;

    const Trajectory* traj_;
    SourceParams      p_;

    SourceBuffer      buf_;
    bool              buffered_ = false;
};

} // namespace sim
} // namespace ofc
#endif // OFC_SIM_SYNTHETIC_SOURCE_HPP
