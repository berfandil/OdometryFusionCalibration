// ofc_sim/synthetic_source.cpp — see synthetic_source.hpp for the measurement model.
#include "ofc_sim/synthetic_source.hpp"

#include "ofc/core/lie.hpp"

#include <cmath>
#include <random>

namespace ofc {
namespace sim {

namespace {

// A deterministic per-window seed: noise on a window is a pure function of the source
// seed and the window endpoints (NOT of the call sequence) -> the rig is deterministic
// regardless of tick alignment. A cheap splitmix-style mix of the three values.
std::uint64_t mix_seed(std::uint32_t base, Timestamp t0, Timestamp t1) {
    std::uint64_t z = static_cast<std::uint64_t>(base);
    z ^= static_cast<std::uint64_t>(t0) + 0x9E3779B97F4A7C15ULL + (z << 6) + (z >> 2);
    z ^= static_cast<std::uint64_t>(t1) + 0x9E3779B97F4A7C15ULL + (z << 6) + (z >> 2);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z =  z ^ (z >> 31);
    return z;
}

} // namespace

SyntheticSource::SyntheticSource(const Trajectory& traj, const SourceParams& p)
    : traj_(&traj), p_(p) {}

SE3 SyntheticSource::base_delta(Timestamp t0, Timestamp t1) const {
    // CANONICAL SIGN: positive time_offset_s shifts the SAMPLED window LATER, so a query
    // for [t0,t1] reads the true trajectory over [t0+off, t1+off]; A = pose(ts)^{-1} o
    // pose(te). (This is the offset Slice-5 time-sync must recover with the right sign.)
    const Timestamp off = static_cast<Timestamp>(
        std::llround(p_.time_offset_s * Scalar(1e9)));
    const SE3 pose_ts = traj_->pose(t0 + off);
    const SE3 pose_te = traj_->pose(t1 + off);
    return se3::compose(se3::inverse(pose_ts), pose_te);
}

SE3 SyntheticSource::clean_reported_(Timestamp t0, Timestamp t1) const {
    const SE3 A = base_delta(t0, t1);
    // B = X^{-1} o A o X  (hand-eye inverse of the estimator's A = X o B o X^{-1}).
    const SE3 Xinv = se3::inverse(p_.X);
    SE3 B = se3::compose(se3::compose(Xinv, A), p_.X);
    // Scale the translation part (per-source translation scale, D20).
    B.t = p_.scale * B.t;
    return B;
}

Mat6 SyntheticSource::modeled_cov_(const SE3& clean_delta) const {
    // Covariance consistent with the injected noise: per-axis variance from the
    // distance/angle-scaled sigma (floored), in [trans; rot] order (matches Delta::cov).
    // Computed from the CLEAN pre-noise delta so the reported Sigma equals the sigma the
    // noise was actually drawn with (review fix #5). On an outlier window the caller
    // passes the CLEAN window delta here, not the gross outlier delta, so the reported
    // Sigma is normal and does not down-weight the outlier (review fix #6).
    const Vec6 xi = se3::log(clean_delta);
    const Scalar dist  = xi.head<3>().norm();
    const Scalar angle = xi.tail<3>().norm();
    const Scalar st = std::max(p_.noise_trans_floor, p_.noise_trans_per_m  * dist);
    const Scalar sr = std::max(p_.noise_rot_floor,   p_.noise_rot_per_rad  * angle);
    // Guard against an all-zero covariance (the estimator inverts a confidence from it).
    const Scalar vt = std::max(st * st, Scalar(1e-12));
    const Scalar vr = std::max(sr * sr, Scalar(1e-12));
    Mat6 cov = Mat6::Zero();
    cov(0, 0) = vt; cov(1, 1) = vt; cov(2, 2) = vt;
    cov(3, 3) = vr; cov(4, 4) = vr; cov(5, 5) = vr;
    return cov;
}

bool SyntheticSource::in_any_(const std::vector<Window>& ws,
                              Timestamp t0, Timestamp t1) const {
    // A window is "active" for [t0,t1] if the (offset-shifted) MIDPOINT falls in it — the
    // fault applies to the motion sampled over this query. Membership is ALL-OR-NOTHING:
    // a query whose [t0,t1] straddles a fault-window edge is faulted (or not) for the WHOLE
    // delta based solely on its midpoint; there is no partial-window fault. So a window
    // boundary landing mid-query silently flips the entire reported delta. Acceptable for
    // an oracle (faults are coarse-grained by design); the midpoint is shifted by the same
    // time_offset_s as the sampled motion so the fault tracks the shifted window.
    const Scalar mid_s =
        (static_cast<Scalar>(t0) + static_cast<Scalar>(t1)) * Scalar(0.5) / Scalar(1e9)
        + p_.time_offset_s;
    for (const Window& w : ws) {
        if (w.contains(mid_s)) return true;
    }
    return false;
}

Expected<Delta> SyntheticSource::query(Timestamp t0, Timestamp t1) const {
    if (t1 <= t0) return Status::OutOfRange;

    // Dropout: report a sensor gap.
    if (in_any_(p_.dropout_windows, t0, t1)) return Status::NoData;

    // The CLEAN window delta (the reported sensor-frame motion absent any fault). Used
    // both as the noise/covariance reference and as the reported delta on clean windows.
    // On an outlier window we still derive Sigma and the noise sigmas from this clean
    // delta so the reported Sigma is normal (review fixes #5/#6).
    const SE3 B_clean = clean_reported_(t0, t1);

    SE3 B;
    if (in_any_(p_.outlier_windows, t0, t1)) {
        // Gross-wrong delta: integrate the outlier body twist over the window.
        const Scalar dt = static_cast<Scalar>(t1 - t0) / Scalar(1e9);
        B = se3::exp(p_.outlier_twist * dt);
    } else {
        B = B_clean;
    }

    // Add zero-mean Gaussian noise in the body tangent (deterministic per window). The
    // noise sigmas scale with the CLEAN delta's magnitude (the motion actually intended),
    // so the reported Sigma below matches the sigma the noise was drawn with (fix #5).
    const bool noisy = (p_.noise_trans_per_m   > Scalar(0)) ||
                       (p_.noise_rot_per_rad   > Scalar(0)) ||
                       (p_.noise_trans_floor   > Scalar(0)) ||
                       (p_.noise_rot_floor     > Scalar(0));
    if (noisy) {
        const Vec6 xi = se3::log(B_clean);
        const Scalar dist  = xi.head<3>().norm();
        const Scalar angle = xi.tail<3>().norm();
        const Scalar st = std::max(p_.noise_trans_floor, p_.noise_trans_per_m * dist);
        const Scalar sr = std::max(p_.noise_rot_floor,   p_.noise_rot_per_rad * angle);

        std::mt19937_64 gen(mix_seed(p_.seed, t0, t1));
        std::normal_distribution<Scalar> nt(Scalar(0), st);
        std::normal_distribution<Scalar> nr(Scalar(0), sr);
        // Draw into locals in EXPLICIT statement order, then assign. An Eigen comma-
        // initializer (eps << nt(gen), ...) leaves the relative evaluation order of the
        // six generator calls unspecified pre-C++17 (this target is C++14), so the noise
        // would not be bit-reproducible across toolchains (review fix #3, DESIGN §10).
        const Scalar e0 = nt(gen);
        const Scalar e1 = nt(gen);
        const Scalar e2 = nt(gen);
        const Scalar e3 = nr(gen);
        const Scalar e4 = nr(gen);
        const Scalar e5 = nr(gen);
        Vec6 eps;
        eps << e0, e1, e2, e3, e4, e5;
        // Right-perturb in the body tangent: B <- B o exp(eps).
        B = se3::compose(B, se3::exp(eps));
    }

    Delta d;
    d.motion = B;
    d.t0     = t0;
    d.t1     = t1;
    // Cov from the CLEAN delta (normal even on outlier windows) — fixes #5/#6.
    d.cov    = p_.report_native_cov ? modeled_cov_(B_clean) : Mat6::Identity();
    return d;
}

// --- Buffer-backed mode ----------------------------------------------------

Status SyntheticSource::build_buffer(Scalar from_s, Scalar to_s, Scalar rate_hz,
                                     OdomForm form, int capacity) {
    if (rate_hz <= Scalar(0) || to_s <= from_s) return Status::OutOfRange;
    const int n = static_cast<int>(std::floor((to_s - from_s) * rate_hz)) + 1;
    const int cap = (capacity > 0) ? capacity : (n + 2);

    SensorConfig sc;
    sc.id = p_.id;
    // The buffer feeds its native cov straight through (the sim already models cov).
    const Status cs = buf_.configure(sc, form, cap, ConfCombine::NativeOnly);
    if (!ok(cs)) return cs;
    buf_.reset();

    // Sample the reported sensor-frame motion as INCREMENTS between consecutive sample
    // times. The buffer integrates increments; query(t0,t1) then yields the same
    // sensor-frame delta the analytic path reports (minus per-window noise, which the
    // buffered path does not reproduce — analytic is the noise oracle). Stamps are the
    // *query* timeline (no offset baked in) so the estimator's windows line up; the
    // offset already lives inside clean_reported_ via base_delta's shift.
    const Timestamp step = static_cast<Timestamp>(std::llround(Scalar(1e9) / rate_hz));
    const Timestamp t0ns = static_cast<Timestamp>(std::llround(from_s * Scalar(1e9)));

    SE3 prev = SE3{};  // identity at the first sample (relative form)
    for (int k = 0; k < n; ++k) {
        const Timestamp ts = t0ns + static_cast<Timestamp>(k) * step;
        if (form == OdomForm::Increment) {
            // Increment from the previous sample time to this one, in the sensor frame.
            const Timestamp tp = t0ns + static_cast<Timestamp>(k - 1) * step;
            const SE3 incr = (k == 0) ? SE3{} : clean_reported_(tp, ts);
            const Status ps = buf_.push_increment(ts, incr, Mat6::Identity());
            if (!ok(ps)) return ps;
        } else if (form == OdomForm::AbsolutePose) {
            // Cumulative reported pose from the first sample (sensor-frame integration).
            prev = (k == 0) ? SE3{}
                            : se3::compose(prev, clean_reported_(t0ns +
                                  static_cast<Timestamp>(k - 1) * step, ts));
            const Status ps = buf_.push_absolute(ts, prev, Mat6::Identity());
            if (!ok(ps)) return ps;
        } else { // Twist
            // Body twist of the reported motion at this sample (sensor frame).
            Vec6 xi = Vec6::Zero();
            if (k != 0) {
                const SE3 incr = clean_reported_(
                    t0ns + static_cast<Timestamp>(k - 1) * step, ts);
                xi = se3::log(incr) * rate_hz;
            }
            const Status ps = buf_.push_twist(ts, xi, Mat6::Identity());
            if (!ok(ps)) return ps;
        }
    }
    buffered_ = true;
    return Status::Ok;
}

Expected<Delta> SyntheticSource::query_buffered(Timestamp t0, Timestamp t1) const {
    if (!buffered_) return Status::NoData;
    if (in_any_(p_.dropout_windows, t0, t1)) return Status::NoData;
    return buf_.query(t0, t1);
}

} // namespace sim
} // namespace ofc
