// test_sim.cpp — the simulation rig (ground-truth oracle, sim/ relaxed edge, D24).
//
// Covers: Trajectory (continuity, twist == finite-diff of pose, regime presets),
// SyntheticSource (identity reproduces the GT base delta, planted X/scale/offset,
// outlier + dropout injection), and the Rig end-to-end (fused frontier tracks GT with
// priors == planted; deterministic replay).
#include <doctest/doctest.h>

#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"

#include "ofc_sim/rig.hpp"
#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <cmath>
#include <memory>
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;
Timestamp secs(double s) { return static_cast<Timestamp>(std::llround(s * 1e9)); }

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

bool se3_close(const SE3& a, const SE3& b, Scalar tol = 1e-9) {
    return close(a.R, b.R, tol) && close(a.t, b.t, tol);
}
} // namespace

// ===========================================================================
// Trajectory
// ===========================================================================

TEST_CASE("trajectory: pose is continuous across segment boundaries") {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    Vec6 turn;     turn     << 2.0, 0, 0, 0, 0, 0.5;
    tr.add_segment(straight, 1.0);
    tr.add_segment(turn,     1.0);

    // Approaching the boundary from either side gives the same pose (continuity).
    const SE3 just_before = tr.pose(secs(1.0 - 1e-6));
    const SE3 at_boundary = tr.pose(secs(1.0));
    const SE3 just_after  = tr.pose(secs(1.0 + 1e-6));
    CHECK(se3_close(just_before, at_boundary, 1e-4));
    CHECK(se3_close(at_boundary, just_after, 1e-4));

    // Boundary pose equals the integrated first segment (pure straight, 2 m over 1 s).
    CHECK(close(at_boundary.t, Vec3(2.0, 0, 0), 1e-9));
    CHECK(close(at_boundary.R, Mat3::Identity(), 1e-9));
}

TEST_CASE("trajectory: twist(t) matches a finite-difference of pose(t)") {
    Trajectory tr = Trajectory::mixed();
    const double h = 1e-4;
    // Sample interior points of several segments, avoiding boundaries.
    for (double t = 0.2; t < tr.duration_s() - 0.2; t += 0.37) {
        const SE3 p0 = tr.pose(secs(t - h));
        const SE3 p1 = tr.pose(secs(t + h));
        // Body twist from the relative motion: log(p0^{-1} p1) / (2h).
        const SE3 rel = se3::compose(se3::inverse(p0), p1);
        const Vec6 fd = se3::log(rel) / (2 * h);
        const Vec6 an = tr.twist(secs(t));
        CHECK(close(fd, an, 1e-3));
    }
}

TEST_CASE("trajectory: straight() preset has ||omega|| == 0 and ||v|| > 0") {
    Trajectory tr = Trajectory::straight(2.0, 5.0);
    for (double t = 0.1; t < 4.9; t += 0.5) {
        const Vec6 xi = tr.twist(secs(t));
        CHECK(xi.tail<3>().norm() == doctest::Approx(0.0));
        CHECK(xi.head<3>().norm() > 0.5);
    }
    // Net motion is a pure forward translation, no rotation.
    const SE3 end = tr.pose(secs(5.0));
    CHECK(close(end.R, Mat3::Identity(), 1e-9));
    CHECK(end.t.x() == doctest::Approx(10.0));
}

TEST_CASE("trajectory: turning() preset has ||omega|| > 0") {
    Trajectory tr = Trajectory::turning(2.0, 0.5, 5.0);
    for (double t = 0.1; t < 4.9; t += 0.5) {
        const Vec6 xi = tr.twist(secs(t));
        CHECK(xi.tail<3>().norm() > 0.1);
    }
    // It actually rotated (yaw accumulates).
    const SE3 end = tr.pose(secs(5.0));
    CHECK(so3::log(end.R).norm() > 0.5);
}

TEST_CASE("trajectory: omega_varying() preset changes ||omega|| over time") {
    Trajectory tr = Trajectory::omega_varying();
    // There exist times with near-zero turn rate AND times with high turn rate.
    Scalar min_w = 1e9, max_w = 0;
    for (double t = 0.05; t < tr.duration_s(); t += 0.05) {
        const Scalar w = tr.twist(secs(t)).tail<3>().norm();
        min_w = std::min(min_w, w);
        max_w = std::max(max_w, w);
    }
    CHECK(min_w < 0.05);     // a straight (zero-omega) stretch exists
    CHECK(max_w > 0.5);      // a strong-turn stretch exists
}

TEST_CASE("trajectory: omega_ramp() preset ramps ||omega|| up then down (smooth)") {
    const Scalar peak = 1.0;
    Trajectory tr = Trajectory::omega_ramp(2.0, peak, 60, 0.02);
    // ‖ω‖(t) climbs from ~0 to the apex over the first half, then falls back to ~0.
    const Scalar mid = tr.duration_s() * 0.5;
    const Scalar early = tr.twist_s(mid * 0.5).tail<3>().norm();
    const Scalar apex  = tr.twist_s(mid).tail<3>().norm();
    const Scalar late  = tr.twist_s(mid + mid * 0.5).tail<3>().norm();
    CHECK(apex > early);                 // rising limb
    CHECK(apex > late);                  // falling limb
    CHECK(apex <= peak + 1e-9);          // capped at the apex yaw rate
    // Adjacent micro-segments differ by only a small step (the ramp is near-smooth),
    // so the per-step ‖ω‖ increment stays well below the apex magnitude.
    const Scalar w0 = tr.twist_s(0.01).tail<3>().norm();
    const Scalar w1 = tr.twist_s(0.03).tail<3>().norm();
    CHECK(std::abs(w1 - w0) < 0.1 * peak);
}

TEST_CASE("trajectory: clamps before t0 and after the end") {
    Trajectory tr = Trajectory::straight(2.0, 1.0);
    // Before start: identity pose, zero twist.
    CHECK(se3_close(tr.pose(secs(-1.0)), SE3{}, 1e-12));
    CHECK(tr.twist(secs(-1.0)).norm() == doctest::Approx(0.0));
    // After end: final pose, zero twist (no extrapolation of motion).
    const SE3 end = tr.pose(secs(1.0));
    CHECK(se3_close(tr.pose(secs(5.0)), end, 1e-12));
    CHECK(tr.twist(secs(5.0)).norm() == doctest::Approx(0.0));
}

// ===========================================================================
// SyntheticSource
// ===========================================================================

TEST_CASE("source: identity X, scale 1, no noise reproduces the GT base delta") {
    Trajectory tr = Trajectory::mixed();
    SourceParams p;
    p.id = 0;                    // X = identity, scale = 1, no noise (defaults)
    SyntheticSource src(tr, p);

    for (double t0 = 0.0; t0 < tr.duration_s() - 0.3; t0 += 0.41) {
        const Timestamp a = secs(t0), b = secs(t0 + 0.2);
        const Expected<Delta> q = src.query(a, b);
        REQUIRE(q.ok());
        // GT base delta over [t0, t0+0.2].
        const SE3 A = se3::compose(se3::inverse(tr.pose(a)), tr.pose(b));
        CHECK(se3_close(q.value().motion, A, 1e-9));
    }
}

TEST_CASE("source: planted X reports the conjugated (hand-eye inverse) motion") {
    Trajectory tr = Trajectory::turning(2.0, 0.5, 5.0);
    SE3 X;
    X.R = so3::exp(Vec3(0.1, -0.2, kPi / 4));   // a real 3-DOF mount rotation
    X.t = Vec3(0.4, -0.3, 0.2);

    SourceParams p;
    p.id = 1;
    p.X  = X;
    SyntheticSource src(tr, p);

    const Timestamp a = secs(1.0), b = secs(1.5);
    const Expected<Delta> q = src.query(a, b);
    REQUIRE(q.ok());

    const SE3 A = se3::compose(se3::inverse(tr.pose(a)), tr.pose(b));
    // Reported B should be X^{-1} A X ...
    const SE3 B_expected =
        se3::compose(se3::compose(se3::inverse(X), A), X);
    CHECK(se3_close(q.value().motion, B_expected, 1e-9));
    // ... and the estimator's frame-align A' = X B X^{-1} recovers A exactly.
    const SE3 A_recovered =
        se3::compose(se3::compose(X, q.value().motion), se3::inverse(X));
    CHECK(se3_close(A_recovered, A, 1e-9));
}

TEST_CASE("source: planted X under STRAIGHT motion (disambiguates X.t handedness)") {
    // Review fix #8. Under pure-straight motion A = {I, d} (no rotation), the hand-eye
    // conjugation B = X^{-1} A X collapses to B = {I, X.R^T d}: the rotation cancels and
    // the reported translation is the base translation expressed in the SENSOR frame via
    // X.R^T. This isolates the rotation handedness of X (a transposed X.R yields X.R d,
    // not X.R^T d) and confirms the round-trip A = X B X^{-1} restores X.t's contribution
    // — a case turning() motion masks because there both R and t move together.
    Trajectory tr = Trajectory::straight(2.0, 5.0);
    SE3 X;
    X.R = so3::exp(Vec3(0.1, -0.2, kPi / 4));   // a real 3-DOF mount rotation
    X.t = Vec3(0.4, -0.3, 0.2);

    SourceParams p;
    p.id = 1;
    p.X  = X;
    SyntheticSource src(tr, p);

    const Timestamp a = secs(1.0), b = secs(1.5);
    const Expected<Delta> q = src.query(a, b);
    REQUIRE(q.ok());

    const SE3 A = se3::compose(se3::inverse(tr.pose(a)), tr.pose(b));
    // A is a pure forward translation (no rotation) over this straight window.
    REQUIRE(close(A.R, Mat3::Identity(), 1e-9));
    REQUIRE(A.t.tail<2>().norm() < 1e-9);     // motion is along +x only

    // Reported B == X^{-1} A X. With A.R == I this equals {I, X.R^T * A.t}.
    const SE3 B_expected =
        se3::compose(se3::compose(se3::inverse(X), A), X);
    CHECK(se3_close(q.value().motion, B_expected, 1e-9));
    CHECK(close(q.value().motion.R, Mat3::Identity(), 1e-9));     // rotation cancels
    CHECK(close(q.value().motion.t, X.R.transpose() * A.t, 1e-9)); // X.R^T d, not X.R d
    // The transposed-X.R alternative (the bug this guards) gives a DIFFERENT translation.
    CHECK_FALSE(close(q.value().motion.t, X.R * A.t, 1e-6));

    // And the estimator's frame-align A' = X B X^{-1} recovers A exactly (X.t restored).
    const SE3 A_recovered =
        se3::compose(se3::compose(X, q.value().motion), se3::inverse(X));
    CHECK(se3_close(A_recovered, A, 1e-9));
}

TEST_CASE("source: scale multiplies the translation part of the reported delta") {
    Trajectory tr = Trajectory::straight(2.0, 5.0);
    SourceParams p1; p1.id = 0; p1.scale = 1.0;
    SourceParams p2; p2.id = 1; p2.scale = 1.3;
    SyntheticSource s1(tr, p1), s2(tr, p2);

    const Timestamp a = secs(0.5), b = secs(1.0);
    const Delta d1 = s1.query(a, b).value();
    const Delta d2 = s2.query(a, b).value();
    // Identity mount + pure-straight motion: rotation identical, translation scaled.
    CHECK(close(d2.motion.R, d1.motion.R, 1e-12));
    CHECK(close(d2.motion.t, 1.3 * d1.motion.t, 1e-9));
}

TEST_CASE("source: a non-zero time offset shifts the sampled interval") {
    Trajectory tr = Trajectory::omega_varying();
    const Scalar off = 0.5;
    SourceParams p; p.id = 0; p.time_offset_s = off;
    SyntheticSource src(tr, p);

    const Timestamp a = secs(0.7), b = secs(1.1);
    const Delta d = src.query(a, b).value();
    // The reported delta corresponds to the window shifted by +off.
    const SE3 A_shifted =
        se3::compose(se3::inverse(tr.pose(a + secs(off))), tr.pose(b + secs(off)));
    CHECK(se3_close(d.motion, A_shifted, 1e-9));
    // And it does NOT match the unshifted window (the offset is observable).
    const SE3 A_unshifted =
        se3::compose(se3::inverse(tr.pose(a)), tr.pose(b));
    CHECK_FALSE(se3_close(d.motion, A_unshifted, 1e-3));
}

TEST_CASE("source: noise is zero-mean, bounded, and deterministic per window") {
    Trajectory tr = Trajectory::straight(2.0, 5.0);
    SourceParams p;
    p.id = 0;
    p.noise_trans_per_m = 0.02;
    p.noise_rot_per_rad = 0.0;
    p.noise_trans_floor = 0.005;
    p.seed = 12345;
    SyntheticSource src(tr, p);

    const Timestamp a = secs(1.0), b = secs(1.2);
    const SE3 clean = se3::compose(se3::inverse(tr.pose(a)), tr.pose(b));

    // Deterministic: same window -> bit-identical reported delta across calls.
    const Delta d1 = src.query(a, b).value();
    const Delta d2 = src.query(a, b).value();
    CHECK((d1.motion.R.array() == d2.motion.R.array()).all());
    CHECK((d1.motion.t.array() == d2.motion.t.array()).all());

    // Noise is small relative to the ~0.4 m motion but non-zero.
    CHECK((d1.motion.t - clean.t).norm() > 0.0);
    CHECK((d1.motion.t - clean.t).norm() < 0.2);

    // Different seed -> different noise.
    SourceParams q = p; q.seed = 999;
    SyntheticSource src2(tr, q);
    const Delta d3 = src2.query(a, b).value();
    CHECK_FALSE(se3_close(d3.motion, d1.motion, 1e-9));

    // Covariance is consistent with the injected noise (positive trans variance).
    CHECK(d1.cov(0, 0) > 0.0);
}

TEST_CASE("source: injected noise is empirically zero-mean over many windows") {
    // Review fix #7: average the recovered noise over many INDEPENDENT windows (the
    // per-window seed is hash(seed,t0,t1), so distinct windows draw independently) and
    // assert the empirical mean is ~0. A planted bias (non-zero-mean noise) would fail.
    Trajectory tr = Trajectory::straight(2.0, 60.0);   // long, constant-magnitude motion
    SourceParams p;
    p.id = 0;
    p.noise_trans_per_m = 0.02;
    p.noise_rot_per_rad = 0.02;
    p.noise_trans_floor = 0.005;
    p.noise_rot_floor   = 0.005;
    p.seed = 7u;
    SyntheticSource src(tr, p);

    // A no-noise twin gives the clean delta to subtract off; recover the injected body-
    // tangent noise as eps = log(B_clean^{-1} o B_noisy).
    SourceParams pc = p;
    pc.noise_trans_per_m = 0; pc.noise_rot_per_rad = 0;
    pc.noise_trans_floor = 0; pc.noise_rot_floor = 0;
    SyntheticSource clean(tr, pc);

    Vec6 sum = Vec6::Zero();
    int  n   = 0;
    // Distinct, non-overlapping 0.2 s windows across the trajectory.
    for (double t0 = 0.5; t0 < 55.0; t0 += 0.2) {
        const Timestamp a = secs(t0), b = secs(t0 + 0.2);
        const SE3 bn = src.query(a, b).value().motion;
        const SE3 bc = clean.query(a, b).value().motion;
        const Vec6 eps = se3::log(se3::compose(se3::inverse(bc), bn));
        sum += eps;
        ++n;
    }
    CHECK(n > 200);
    const Vec6 mean = sum / static_cast<Scalar>(n);
    // Per-axis sigma ~0.005..0.02; with >200 samples the standard error of the mean is
    // ~sigma/sqrt(n) ~ 1e-3, so a true-zero-mean noise lands well under 5e-3. A planted
    // bias of even 0.01 would blow past this.
    CHECK(mean.head<3>().norm() < 5e-3);
    CHECK(mean.tail<3>().norm() < 5e-3);
}

TEST_CASE("source: dropout window returns NoData") {
    Trajectory tr = Trajectory::mixed();
    SourceParams p;
    p.id = 0;
    p.dropout_windows.push_back(Window{1.0, 2.0});
    SyntheticSource src(tr, p);

    // Mid-dropout: NoData.
    CHECK(src.query(secs(1.4), secs(1.6)).status() == Status::NoData);
    // Outside the dropout: Ok.
    CHECK(src.query(secs(0.2), secs(0.4)).ok());
    CHECK(src.query(secs(2.5), secs(2.7)).ok());
}

TEST_CASE("source: outlier window yields a grossly wrong delta") {
    Trajectory tr = Trajectory::straight(2.0, 5.0);
    SourceParams p;
    p.id = 0;
    p.outlier_windows.push_back(Window{2.0, 3.0});
    SyntheticSource src(tr, p);

    const Timestamp a = secs(2.4), b = secs(2.6);
    const SE3 clean = se3::compose(se3::inverse(tr.pose(a)), tr.pose(b));
    const Delta bad = src.query(a, b).value();
    // The outlier delta is far from the clean motion (gross lateral + yaw error).
    CHECK((bad.motion.t - clean.t).norm() > 0.5);

    // A window outside the outlier range is clean.
    const Delta good = src.query(secs(0.4), secs(0.6)).value();
    const SE3 cg = se3::compose(se3::inverse(tr.pose(secs(0.4))), tr.pose(secs(0.6)));
    CHECK(se3_close(good.motion, cg, 1e-9));
}

TEST_CASE("source: outlier window reports a NORMAL covariance (not synthesized from the "
          "gross delta)") {
    // Review fix #6: the outlier Sigma must be derived from the CLEAN window delta, not
    // the gross outlier delta. Otherwise a huge synthesized variance would let
    // Sigma-downweighting (not pure geometric-median rejection) reject the outlier.
    Trajectory tr = Trajectory::straight(2.0, 5.0);
    SourceParams p;
    p.id = 0;
    p.noise_trans_per_m = 0.02;          // a normal noise model
    p.noise_rot_per_rad = 0.01;
    p.outlier_windows.push_back(Window{2.0, 3.0});
    SyntheticSource src(tr, p);

    const Timestamp a = secs(2.4), b = secs(2.6);     // inside the outlier window
    const Timestamp c = secs(0.4), d = secs(0.6);     // a clean window, same duration
    const Delta bad   = src.query(a, b).value();
    const Delta clean = src.query(c, d).value();

    // The outlier delta is gross...
    const SE3 cg = se3::compose(se3::inverse(tr.pose(a)), tr.pose(b));
    CHECK((bad.motion.t - cg.t).norm() > 0.5);
    // ...but its reported Sigma equals the clean window's Sigma (same motion magnitude),
    // i.e. it is NOT inflated by the gross delta. (Straight motion of equal length over
    // both windows -> identical clean dist/angle -> identical modeled cov.)
    CHECK(close(bad.cov, clean.cov, 1e-12));
    // And the trans variance is the small noise-model variance, far below the gross
    // delta's squared magnitude (which would be O(1) m^2 if synthesized from the outlier).
    CHECK(bad.cov(0, 0) < 0.05);
}

TEST_CASE("source: buffer-backed mode reproduces the analytic clean delta") {
    Trajectory tr = Trajectory::turning(2.0, 0.4, 5.0);
    SE3 X; X.R = so3::exp(Vec3(0, 0, kPi / 6)); X.t = Vec3(0.2, 0.1, 0);
    SourceParams p; p.id = 0; p.X = X; p.scale = 1.1;
    SyntheticSource src(tr, p);

    REQUIRE(ok(src.build_buffer(0.0, 5.0, 500.0, OdomForm::Increment)));
    // Query a window well inside the buffered span; the buffered increment-integrated
    // delta should match the analytic clean (no-noise) reported delta closely.
    const Timestamp a = secs(1.0), b = secs(1.4);
    const Expected<Delta> qb = src.query_buffered(a, b);
    REQUIRE(qb.ok());
    SourceParams pc = p;                    // a no-noise twin for the analytic ref
    SyntheticSource clean(tr, pc);
    const SE3 ref = clean.query(a, b).value().motion;
    CHECK(se3_close(qb.value().motion, ref, 2e-3));
}

// ===========================================================================
// Rig end-to-end
// ===========================================================================

namespace {
// Build a rig config with N sources whose priors EQUAL the planted params (so fusion
// needs no calibration). Returns the SensorConfig storage by out-param (must outlive
// the estimator).
Config make_rig_config(const std::vector<SourceParams>& planted,
                       std::vector<SensorConfig>& sensors_out) {
    sensors_out.clear();
    for (const SourceParams& sp : planted) {
        SensorConfig sc;
        sc.id              = sp.id;
        sc.prior_extrinsic = sp.X;          // prior == planted
        sc.prior_scale     = sp.scale;
        sc.weight_prior    = 1.0;
        // scale_calib stays ON (the SensorConfig default). With the scale_hist default range
        // fixed to [0.5, 1.5] (1.0 strictly interior), a source whose true residual scale is a
        // UNIT ratio now COMMITS 1.0 (not the old ~0.984 boundary-bin artifact), so the
        // rising-edge feedback leaves prior_scale unchanged and the fused frontier tracks GT
        // tightly. This re-enables the scale-calibration coverage the old workaround removed.
        sensors_out.push_back(sc);
    }
    Config c;
    c.max_sources    = static_cast<int>(planted.size());
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.sensors        = sensors_out.data();
    c.sensor_count   = static_cast<int>(sensors_out.size());
    // Priors == planted here (no calibration needed), so every source drives the median
    // from the start (the Slice-2 semantics). ReferenceOnly cold-start is exercised in
    // test_calib_feedback.
    c.cold_start     = ColdStart::MedianFromStart;
    return c;
}
} // namespace

TEST_CASE("rig: fused frontier tracks GT with >=3 sources, priors == planted") {
    Trajectory tr = Trajectory::mixed();

    // Three sources: one aligned reference + two with real mounts. Priors == planted,
    // so the estimator's frame-align cancels the mount exactly (no calibration needed).
    std::vector<SourceParams> planted(3);
    planted[0].id = 0;   // reference, identity mount
    planted[1].id = 1;
    planted[1].X.R = so3::exp(Vec3(0, 0, kPi / 5));
    planted[1].X.t = Vec3(0.3, -0.2, 0.1);
    planted[2].id = 2;
    planted[2].X.R = so3::exp(Vec3(0.05, 0.1, -kPi / 6));
    planted[2].X.t = Vec3(-0.25, 0.15, 0.05);

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_config(planted, sensors);

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);

    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    CHECK(fuses > 10);

    Scalar max_te, max_re;
    rig.max_error(max_te, max_re);
    // Clean sources + correct priors -> the fused frontier tracks GT tightly. Tolerance
    // accounts for const-twist integration vs the estimator's windowed composition.
    CHECK(max_te < 5e-2);
    CHECK(max_re < 5e-2);
}

TEST_CASE("rig: fused frontier tracks GT with a planted translation scale != 1") {
    // End-to-end scale gate (review fix #1/#4): each source over-reports its translation
    // by `scale`; the matching SensorConfig::prior_scale equals the planted scale, so the
    // estimator must de-scale B.t before the frame-align to recover the true base motion.
    // This FAILS if the estimator ignores prior_scale (the fused frontier drifts by the
    // scale factor). Pure-straight motion makes the translation error maximally visible.
    Trajectory tr = Trajectory::mixed();

    std::vector<SourceParams> planted(3);
    planted[0].id = 0;   planted[0].scale = 1.2;   // identity mount, over-reports 20%
    planted[1].id = 1;   planted[1].scale = 0.8;   // under-reports 20%
    planted[1].X.R = so3::exp(Vec3(0, 0, kPi / 5));
    planted[1].X.t = Vec3(0.3, -0.2, 0.1);
    planted[2].id = 2;   planted[2].scale = 1.5;   // over-reports 50%
    planted[2].X.R = so3::exp(Vec3(0.05, 0.1, -kPi / 6));
    planted[2].X.t = Vec3(-0.25, 0.15, 0.05);

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_config(planted, sensors);   // prior_scale == planted scale

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);

    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    CHECK(fuses > 10);

    Scalar max_te, max_re;
    rig.max_error(max_te, max_re);
    // With prior_scale applied the de-scaled deltas equal the true base motion -> the
    // fused frontier tracks GT within the same tolerance as the unit-scale rig.
    CHECK(max_te < 5e-2);
    CHECK(max_re < 5e-2);
}

TEST_CASE("rig: fused frontier tracks the OFFSET-shifted GT with a planted time offset") {
    // End-to-end time-offset gate (review fix #4). Every source shares the SAME planted
    // offset `off` (positive => the sampled window is shifted LATER, canonical sign). The
    // Slice-2 estimator does NOT time-sync (that is Slice 5), so with priors == planted it
    // integrates the base motion the sources report — which is the GT motion over the
    // OFFSET-shifted timeline. The fused frontier therefore tracks GT(t + off), not GT(t).
    //
    // GT-ANCHOR PIN (finding rig.cpp:61): the rig's built-in pose_error anchors at the
    // UNSHIFTED GT (pose(first_frontier - window_s)), so it would (correctly) report a
    // mismatch here. Because all sources share one offset they still share the frontier
    // stamp, so the single-anchor assumption holds — but the anchor must be taken on the
    // SHIFTED timeline. We therefore reproduce the rig's anchor logic with the +off shift:
    //   origin   = GT pose at (first_fused_frontier - window_s + off)
    //   GT(rec)  = origin^{-1} o GT pose at (frontier_stamp + off)
    // and compare the fused frontier pose against THAT. This both documents the convention
    // and proves the offset propagates with the correct sign end-to-end.
    Trajectory tr = Trajectory::mixed();
    const Scalar off    = 0.1;
    const Scalar wind_s = 0.10;     // == cfg.window_s below (the bootstrap lookback)

    std::vector<SourceParams> planted(3);
    planted[0].id = 0;   planted[0].time_offset_s = off;   // identity mount
    planted[1].id = 1;   planted[1].time_offset_s = off;
    planted[1].X.R = so3::exp(Vec3(0, 0, kPi / 5));
    planted[1].X.t = Vec3(0.3, -0.2, 0.1);
    planted[2].id = 2;   planted[2].time_offset_s = off;
    planted[2].X.R = so3::exp(Vec3(0.05, 0.1, -kPi / 6));
    planted[2].X.t = Vec3(-0.25, 0.15, 0.05);

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_config(planted, sensors);
    cfg.window_s = wind_s;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);

    const int fuses = rig.run(0.2, tr.duration_s() - 0.3, 50.0);
    CHECK(fuses > 10);

    // Offset-aware GT comparison (the pinned convention above).
    const Timestamp off_ns  = secs(off);
    const Timestamp wind_ns = secs(wind_s);
    bool have_anchor = false;
    SE3  origin_inv;
    Scalar max_te = 0, max_re = 0;
    for (const Record& r : rig.records()) {
        if (!r.fused) continue;
        const Timestamp frontier = r.result.frontier.stamp;
        if (!have_anchor) {
            origin_inv  = se3::inverse(tr.pose(frontier - wind_ns + off_ns));
            have_anchor = true;
        }
        const SE3 gt_shifted = se3::compose(origin_inv, tr.pose(frontier + off_ns));
        const SE3& fused = r.result.frontier.pose;
        max_te = std::max(max_te, (fused.t - gt_shifted.t).norm());
        max_re = std::max(max_re, so3::log(gt_shifted.R.transpose() * fused.R).norm());
    }
    // The fused frontier tracks the offset-shifted GT within the same tolerance as the
    // unshifted rig (same windowed-integration vs const-twist discretization error).
    CHECK(max_te < 5e-2);
    CHECK(max_re < 5e-2);

    // Sanity: the offset is OBSERVABLE end-to-end — the fused frontier does NOT track the
    // UNSHIFTED GT (otherwise a no-op / sign-flipped offset would pass). Use the rig's
    // own (unshifted-anchored) error, which must be materially larger than the tolerance.
    Scalar un_te, un_re;
    rig.max_error(un_te, un_re);
    CHECK(un_te > 5e-2);
}

TEST_CASE("rig: replay is deterministic (identical recorded stream)") {
    Trajectory tr = Trajectory::mixed();

    auto run_once = [&](std::vector<Record>& out) {
        std::vector<SourceParams> planted(3);
        planted[0].id = 0;
        planted[1].id = 1; planted[1].X.t = Vec3(0.3, 0, 0);
        planted[2].id = 2; planted[2].X.t = Vec3(-0.3, 0, 0);
        // Plant noise on every source — determinism must hold WITH noise.
        for (auto& sp : planted) {
            sp.noise_trans_per_m = 0.01;
            sp.noise_rot_per_rad = 0.005;
            sp.noise_trans_floor = 0.002;
            sp.seed = 100u + sp.id;
        }
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
        std::vector<SensorConfig> sensors;
        Config cfg = make_rig_config(planted, sensors);

        Rig rig;
        rig.set_trajectory(tr);
        rig.init(cfg);
        for (auto& s : srcs) rig.add_source(s.get());
        rig.run(0.2, 5.0, 50.0);
        out = rig.records();
    };

    std::vector<Record> a, b;
    run_once(a);
    run_once(b);
    REQUIRE(a.size() == b.size());
    REQUIRE(!a.empty());
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].fused == b[i].fused);
        if (!a[i].fused) continue;
        const State& fa = a[i].result.frontier;
        const State& fb = b[i].result.frontier;
        // Bit-identical fused pose + covariance (seeded PRNG -> reproducible).
        CHECK((fa.pose.R.array() == fb.pose.R.array()).all());
        CHECK((fa.pose.t.array() == fb.pose.t.array()).all());
        CHECK((fa.twist.xi.array() == fb.twist.xi.array()).all());
        CHECK((fa.cov.array() == fb.cov.array()).all());
        CHECK((a[i].gt_frontier.t.array() == b[i].gt_frontier.t.array()).all());
    }
}

TEST_CASE("rig: an injected outlier source is rejected by the median fusion") {
    Trajectory tr = Trajectory::straight(2.0, 5.0);

    // Three clean sources + one source that goes outlier over a window.
    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].outlier_windows.push_back(Window{1.0, 4.0});

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_config(planted, sensors);

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);

    rig.run(0.2, 4.5, 50.0);
    Scalar max_te, max_re;
    rig.max_error(max_te, max_re);
    // Despite one source emitting a gross-wrong delta over [1,4]s, the median of the
    // remaining three keeps the fused frontier near GT. With fix #6 the outlier reports a
    // NORMAL covariance (Sigma derived from the clean window delta, not the gross delta),
    // so this isolates pure geometric-median robustness — the outlier is NOT down-weighted
    // via an inflated Sigma. The sources are weighted equally (sigma_confidence equal), so
    // rejection here is geometric.
    CHECK(max_te < 1e-1);
    CHECK(max_re < 1e-1);
}
