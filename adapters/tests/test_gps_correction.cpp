// Adapters Slice 11b/13: the production GPS ICorrection adapter. A real GPS stream (geodetic
// lat/lon/alt + ENU position covariance) is turned into the core's dim=3 position correction —
// the SAME path the sim ref (SyntheticAbsoluteRef) exercises — so absolute fixes remove the
// fused-odometry drift. These tests are SELF-CONTAINED: the adapters test exe links only
// ofc_core + ofc_adapters + doctest (NOT ofc_sim), so the end-to-end rig uses adapter_test_util's
// TwistSource / make_rig_cfg. No RNG (deterministic).
#include <doctest/doctest.h>

#include "ofc_adapters/gps_correction.hpp"

#include "ofc/core/eskf.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include "adapter_test_util.hpp"

#include <cmath>
#include <map>
#include <vector>

using namespace ofc;
using namespace ofc::adapters;
using namespace adptest;

namespace {

constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;

// A reference WGS-84 datum (mid-latitude, away from poles/dateline so the local-tangent
// approximations the tests assert are well behaved).
GpsConfig datum_cfg(double lat0_deg, double lon0_deg, double alt0_m) {
    GpsConfig c;
    c.has_datum     = true;
    c.datum_lat_deg = lat0_deg;
    c.datum_lon_deg = lon0_deg;
    c.datum_alt_m   = alt0_m;
    return c;
}

GpsFix make_fix(double lat, double lon, double alt) {
    GpsFix f;
    f.lat_deg = lat;
    f.lon_deg = lon;
    f.alt_m   = alt;
    f.stamp   = 0;
    f.valid   = true;
    return f;
}

} // namespace

// ---------------------------------------------------------------------------------------------
// 1. Geodetic -> ENU correctness. A small latitude delta is ~north, a small longitude delta is
//    ~east*cos(lat0), a fix AT the datum is ~0. (1e-5 deg ~ 1.11 m on the ground.)
// ---------------------------------------------------------------------------------------------
TEST_CASE("GPS geodetic->ENU: small lat/lon deltas map to the expected N/E meters") {
    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    GpsCorrection gps(datum_cfg(lat0, lon0, alt0));

    // A fix exactly AT the datum -> ENU ~ 0.
    State x;  // identity pose
    Measurement m;
    gps.submit_fix(make_fix(lat0, lon0, alt0));
    CHECK(gps.evaluate(x, m));
    CHECK(gps.last_enu().norm() < 1e-6);

    // +1e-5 deg latitude -> ~ +1.11 m North (E,U ~ 0).
    gps.submit_fix(make_fix(lat0 + 1e-5, lon0, alt0));
    CHECK(gps.evaluate(x, m));
    Vec3 enu_n = gps.last_enu();
    const double meters_per_deg = 1e-5 * (2.0 * 3.14159265358979323846 * 6378137.0) / 360.0; // ~1.113 m
    CHECK(std::abs(enu_n.x()) < 5e-3);                                 // East ~ 0
    CHECK(enu_n.y() == doctest::Approx(meters_per_deg).epsilon(5e-3)); // North
    CHECK(std::abs(enu_n.z()) < 5e-3);

    // +1e-5 deg longitude -> ~ +1.11 * cos(lat0) m East (N,U ~ 0).
    gps.submit_fix(make_fix(lat0, lon0 + 1e-5, alt0));
    CHECK(gps.evaluate(x, m));
    Vec3 enu_e = gps.last_enu();
    const double east_expected = meters_per_deg * std::cos(lat0 * kDeg2Rad);
    CHECK(enu_e.x() == doctest::Approx(east_expected).epsilon(5e-3)); // East
    CHECK(std::abs(enu_e.y()) < 5e-3);
    CHECK(std::abs(enu_e.z()) < 5e-3);

    // +alt -> straight up.
    gps.submit_fix(make_fix(lat0, lon0, alt0 + 2.0));
    CHECK(gps.evaluate(x, m));
    Vec3 enu_u = gps.last_enu();
    CHECK(enu_u.z() == doctest::Approx(2.0).epsilon(1e-4));
    CHECK(std::abs(enu_u.x()) < 5e-3);
    CHECK(std::abs(enu_u.y()) < 5e-3);
}

// ---------------------------------------------------------------------------------------------
// 2. Datum lazy-latch: with !has_datum the FIRST valid fix latches the datum (its ENU ~ 0) and
//    has_datum() flips true.
// ---------------------------------------------------------------------------------------------
TEST_CASE("GPS datum lazy-latch: first fix sets origin, has_datum flips, its ENU ~ 0") {
    GpsConfig c;  // has_datum = false (default)
    GpsCorrection gps(c);
    CHECK_FALSE(gps.has_datum());

    State x;
    Measurement m;
    gps.submit_fix(make_fix(47.0, 8.0, 400.0));
    CHECK(gps.has_datum());                 // latched on submit
    CHECK(gps.evaluate(x, m));
    CHECK(gps.last_enu().norm() < 1e-6);    // first fix sits at the origin

    // A subsequent fix +1e-5 lat is ~1.11 m north of the (latched) origin.
    gps.submit_fix(make_fix(47.0 + 1e-5, 8.0, 400.0));
    CHECK(gps.evaluate(x, m));
    CHECK(gps.last_enu().y() > 1.0);
}

// ---------------------------------------------------------------------------------------------
// 3. Measurement shape with lever=0: byte-matches the sim-ref model H = [R | 0 | 0],
//    residual = z - x.pose.t, dim == 3.
// ---------------------------------------------------------------------------------------------
TEST_CASE("GPS measurement shape (lever=0): H==[R|0|0], residual==z-x.pose.t, dim==3") {
    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    GpsCorrection gps(datum_cfg(lat0, lon0, alt0));

    // A non-identity frontier pose so R shows up in H.
    State x;
    x.pose.R = so3::exp(Vec3(0.1, -0.2, 0.3));
    x.pose.t = Vec3(5.0, -3.0, 1.0);

    // A fix offset from the datum -> a nonzero ENU = z (identity alignment).
    gps.submit_fix(make_fix(lat0 + 2e-5, lon0 - 1e-5, alt0 + 0.5));
    Measurement m;
    REQUIRE(gps.evaluate(x, m));
    CHECK(m.dim == 3);

    const Vec3 z = gps.last_enu();          // identity odom_from_enu -> z == ENU
    // residual = z - x.pose.t (lever = 0).
    CHECK((m.residual.head<3>() - (z - x.pose.t)).norm() < 1e-9);

    // H pose-translation block == R, rotation block == 0, twist block == 0.
    Eigen::Matrix<Scalar, 3, 12> H = m.H.topRows<3>();
    CHECK((H.block<3, 3>(0, 0) - x.pose.R).norm() < 1e-12);
    CHECK(H.block<3, 3>(0, 3).norm() < 1e-12);
    CHECK(H.block<3, 6>(0, 6).norm() < 1e-12);
}

// ---------------------------------------------------------------------------------------------
// 4. Lever arm: nonzero antenna lever -> H rotation block == -R*[l]x and residual subtracts R*l.
//    A non-identity R makes the difference bite.
// ---------------------------------------------------------------------------------------------
TEST_CASE("GPS lever arm: H rot block == -R*[l]x, residual subtracts R*lever") {
    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    GpsConfig c = datum_cfg(lat0, lon0, alt0);
    const Vec3 lever(0.4, -0.1, 0.2);
    c.lever_arm_base = lever;
    GpsCorrection gps(c);

    State x;
    x.pose.R = so3::exp(Vec3(0.2, 0.3, -0.1));
    x.pose.t = Vec3(1.0, 2.0, 3.0);

    gps.submit_fix(make_fix(lat0 + 1e-5, lon0 + 1e-5, alt0));
    Measurement m;
    REQUIRE(gps.evaluate(x, m));

    const Vec3 z = gps.last_enu();
    const Vec3 h = x.pose.t + x.pose.R * lever;      // predicted antenna position
    CHECK((m.residual.head<3>() - (z - h)).norm() < 1e-9);

    // H = [ R | -R*[l]x | 0 ].
    const Mat3 expect_rot = -x.pose.R * so3::hat(lever);
    Eigen::Matrix<Scalar, 3, 12> H = m.H.topRows<3>();
    CHECK((H.block<3, 3>(0, 0) - x.pose.R).norm() < 1e-12);
    CHECK((H.block<3, 3>(0, 3) - expect_rot).norm() < 1e-12);
    CHECK(H.block<3, 6>(0, 6).norm() < 1e-12);
}

// ---------------------------------------------------------------------------------------------
// 5. Covariance mapping: a non-identity odom_from_enu.R rotates a diagonal ENU cov into the
//    expected R_odom = Ralign * cov_enu * Ralign^T; cov_floor adds to the diagonal.
// ---------------------------------------------------------------------------------------------
TEST_CASE("GPS covariance mapping: rotated ENU cov + floor") {
    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    GpsConfig c = datum_cfg(lat0, lon0, alt0);
    c.odom_from_enu.R = so3::exp(Vec3(0.0, 0.0, 0.5));  // 0.5 rad yaw ENU->odom
    c.cov_floor_m2    = 0.04;
    GpsCorrection gps(c);

    GpsFix f = make_fix(lat0 + 1e-5, lon0, alt0);
    f.cov_enu = Vec3(0.25, 1.0, 4.0).asDiagonal();      // sx^2, sy^2, su^2
    gps.submit_fix(f);

    State x;
    Measurement m;
    REQUIRE(gps.evaluate(x, m));

    const Mat3 Ralign = c.odom_from_enu.R;
    const Mat3 expect = Ralign * (f.cov_enu + c.cov_floor_m2 * Mat3::Identity()) * Ralign.transpose();
    CHECK((m.R.topLeftCorner<3, 3>() - expect).norm() < 1e-9);
}

// ---------------------------------------------------------------------------------------------
// 6. Pending-fix gate: submit one fix -> evaluate emits once, a second evaluate returns false;
//    submit a newer fix -> emits again; an invalid fix is ignored.
// ---------------------------------------------------------------------------------------------
TEST_CASE("GPS pending-fix gate: emit-once, newer replaces, invalid ignored") {
    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    GpsCorrection gps(datum_cfg(lat0, lon0, alt0));

    State x;
    Measurement m;

    // No fix submitted yet -> nothing to emit.
    CHECK_FALSE(gps.evaluate(x, m));

    gps.submit_fix(make_fix(lat0 + 1e-5, lon0, alt0));
    CHECK(gps.evaluate(x, m));     // emits once
    CHECK_FALSE(gps.evaluate(x, m)); // consumed -> false until a newer fix

    // A newer fix re-arms the adapter.
    gps.submit_fix(make_fix(lat0 + 2e-5, lon0, alt0));
    CHECK(gps.evaluate(x, m));
    CHECK_FALSE(gps.evaluate(x, m));

    // An invalid (no-fix) submission is ignored: still nothing to emit.
    GpsFix bad = make_fix(lat0 + 5e-5, lon0, alt0);
    bad.valid = false;
    gps.submit_fix(bad);
    CHECK_FALSE(gps.evaluate(x, m));
}

// ---------------------------------------------------------------------------------------------
// 7. End-to-end drift removal THROUGH the estimator (the payoff). A 2-source rig with a planted
//    translation scale makes the fused odom DRIFT. The drift-free "truth" is a TWIN estimator with
//    scale=1 (same config / window / bootstrap -> the SAME odom frame + anchor, no GT-formula
//    anchor-matching footgun). We GPS-correct the drifting estimator with fixes whose ENU == the
//    twin's frontier position at the matching frontier stamp (datum at the odom origin, identity
//    alignment -> ENU IS the odom position). We assert the tail position error vs the twin is
//    reduced vs a no-GPS baseline, AND that a gross-outlier fix is Mahalanobis-REJECTED at the
//    dim=3 per-n gate (NIS > gate) -> the per-n chi2_gate path is consumed end-to-end.
//    Deterministic (no RNG).
// ---------------------------------------------------------------------------------------------
namespace {

// Inverse local-tangent: turn a target ENU (== odom position, identity alignment) into a geodetic
// fix around (lat0,lon0,alt0). a = WGS-84 semi-major; meters per rad of lat ~ a, per rad of lon ~
// a*cos(lat0). The forward WGS-84 round-trip lands back within a few mm at these ranges, so the
// fix the adapter reconstructs is the intended odom target to mm.
GpsFix enu_to_fix(double lat0, double lon0, double alt0, const Vec3& enu, Timestamp stamp,
                  Scalar sigma_m) {
    const double a = 6378137.0;
    GpsFix f;
    f.lat_deg = lat0 + (enu.y() / a) / kDeg2Rad;                          // North -> dlat
    f.lon_deg = lon0 + (enu.x() / (a * std::cos(lat0 * kDeg2Rad))) / kDeg2Rad; // East -> dlon
    f.alt_m   = alt0 + enu.z();
    f.cov_enu = Mat3::Identity() * (sigma_m * sigma_m);
    f.stamp   = stamp;
    f.valid   = true;
    return f;
}

// Run the drift-free TWIN (scale=1) and record its frontier position at EACH published stamp.
// The drifting estimator shares the config/window, so these positions ARE the drift-free truth in
// the same odom frame; we look them up by stamp to build the GPS target.
std::map<Timestamp, Vec3> truth_positions(const Config& cfg, const Vec6& xi,
                                           const std::vector<Timestamp>& ts) {
    std::vector<TwistSource> srcs;
    srcs.emplace_back(0, xi, Scalar(1.0));
    srcs.emplace_back(1, xi, Scalar(1.0));
    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);
    std::map<Timestamp, Vec3> out;
    for (Timestamp now : ts) {
        if (est.step(now) != Status::Ok) continue;
        const Result& r = est.latest();
        out[r.frontier.stamp] = r.frontier.pose.t;
    }
    return out;
}

} // namespace

TEST_CASE("GPS end-to-end: removes fused-odometry drift through the estimator and rejects an outlier") {
    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_cfg(sensors);
    // chi2=100 IS A TEST ARTIFACT (mirrors tests/test_correction.cpp's drift-removal case), NOT a
    // recommended production value. Root cause: the Slice-14 covariance pessimism — the ESKF inits
    // P=I12 and the predict-only stretches never shrink it, so the accumulated-drift residual's NIS
    // legitimately exceeds the default chi2=9 gate even though the drift IS the signal we want to
    // admit. We widen to 100 ONLY so the mis-calibrated-P case can still demonstrate drift removal;
    // the SAME run independently proves a gross 25 m outlier (NIS ~ (25/0.05)^2 = 2.5e5) is still
    // rejected at this loosened threshold -> the per-n chi2_gate path is exercised end-to-end.
    cfg.mahalanobis_chi2 = 100.0;
    // Process noise: the two identical sources give ZERO median spread, so adaptive_q collapses to
    // q_floor (default 1e-6 -> the position covariance never re-inflates between fixes -> a near-
    // zero Kalman gain -> the fix barely pulls). A realistic per-step odometry process noise
    // (q_floor) lets the position covariance grow each predict, so a GPS fix is appropriately
    // trusted and actually corrects the drift. 1e-3 m^2/step on translation is a modest,
    // physically-sensible "the odometry is uncertain" floor (the rot axes stay tight).
    cfg.q_floor[0] = cfg.q_floor[1] = cfg.q_floor[2] = 1e-3;  // translation [x,y,z]

    Vec6 xi;
    xi << 2.0, 0, 0, 0, 0, 0.3;        // forward + yaw (matches make_rig_sources)
    const Scalar planted_scale = 1.04; // 4% translation over-scale -> growing drift (matches the
                                       // core drift-removal test; gentle enough that the fixes can
                                       // out-pull the per-step drift accumulation)

    // 50 Hz over ~2 s.
    std::vector<Timestamp> ts;
    const Timestamp tick = secs_to_ns(0.02);
    for (Timestamp now = secs_to_ns(0.2); now <= secs_to_ns(2.0); now += tick) ts.push_back(now);

    // The drift-free truth: a scale=1 twin's frontier position per stamp (same odom frame/anchor).
    const std::map<Timestamp, Vec3> truth = truth_positions(cfg, xi, ts);

    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    const Scalar n3_gate = Eskf::chi2_gate(cfg.mahalanobis_chi2, 3);
    const Timestamp frontier_tick = secs_to_ns(0.02); // the published-frontier advance per step
    const Timestamp fix_period = frontier_tick;        // a fix every step (dense -> strong tracking)

    // Tail window: average the fused-vs-truth error over the LAST `tail_window_s` of the run
    // (mirrors tests/test_correction.cpp's tail-mean metric). A single final-frontier sample is
    // noisy under the per-correction sawtooth; the tail mean is the stable steady-state error.
    const Scalar  tail_window_s = 0.5;
    const Timestamp tail_start  = ts.back() - secs_to_ns(tail_window_s);

    // One drifting run over `ts`. If with_gps, feed a fix (== the twin truth at the frontier) every
    // fix_period and inject ONE gross outlier near elapsed~0.8 s. Returns the TAIL-MEAN fused-vs-
    // truth position error + out-params (outlier rejected? fixes applied?). The drifting estimator
    // + GPS adapter share the lifetime of the call (the adapter is registered on `est`).
    auto run = [&](bool with_gps, bool* outlier_rejected, int* applied) -> Scalar {
        std::vector<TwistSource> srcs;
        srcs.emplace_back(0, xi, planted_scale);
        srcs.emplace_back(1, xi, planted_scale);
        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);

        GpsConfig gcfg = datum_cfg(lat0, lon0, alt0);   // identity alignment, zero lever
        GpsCorrection gps(gcfg);
        if (with_gps) REQUIRE(est.add_correction(&gps) == Status::Ok);

        Timestamp next_fix = -1;
        bool outlier_injected = false;
        if (outlier_rejected) *outlier_rejected = false;
        if (applied) *applied = 0;

        Scalar tail_err_sum = 0.0; int tail_n = 0;
        for (Timestamp now : ts) {
            if (est.step(now) != Status::Ok) continue;
            const Result& r = est.latest();
            if (next_fix < 0) next_fix = r.frontier.stamp;  // arm at the first frontier

            // Diagnostics of the step that just ran (it consumed any fix submitted last iteration).
            const CorrectionDiag& cd = r.correction;
            if (applied) *applied += cd.corr_applied;
            if (outlier_rejected && cd.corr_rejected > 0 && cd.last_nis > n3_gate) {
                *outlier_rejected = true;
            }

            // Accumulate the tail error vs the truth twin at this same frontier stamp.
            auto truth_it = truth.find(r.frontier.stamp);
            if (truth_it != truth.end() && r.frontier.stamp >= tail_start) {
                tail_err_sum += (r.frontier.pose.t - truth_it->second).norm();
                ++tail_n;
            }

            // Submit a fix every step (it is consumed by the NEXT step(), since the correction runs
            // inside step()). Dense fixes keep the position covariance pinned near R between
            // corrections, so the steady-state tail error stays low. Target the truth at the
            // consuming frontier (current stamp + one published tick) so the residual reflects the
            // drift at the moment the fix is actually applied.
            if (with_gps && r.frontier.stamp >= next_fix) {
                const Timestamp consume_stamp = r.frontier.stamp + frontier_tick;
                auto it = truth.find(consume_stamp);
                if (it == truth.end()) it = truth.find(r.frontier.stamp); // last step: no look-ahead
                if (it != truth.end()) {
                    Vec3 target = it->second;
                    const Scalar elapsed = static_cast<Scalar>(r.frontier.stamp - ts.front()) / 1e9;
                    const bool outlier = (!outlier_injected && elapsed > 0.7 && elapsed < 0.9);
                    if (outlier) { target += Vec3(25.0, 0.0, 0.0); outlier_injected = true; }
                    gps.submit_fix(enu_to_fix(lat0, lon0, alt0, target, consume_stamp, Scalar(0.05)));
                }
                next_fix += fix_period;
            }
        }
        return (tail_n > 0) ? tail_err_sum / static_cast<Scalar>(tail_n) : 0.0;
    };

    // ---- baseline (no GPS): the drift ---------------------------------------------------------
    const Scalar baseline_tail = run(/*with_gps=*/false, nullptr, nullptr);
    CHECK(baseline_tail > 0.05);   // the planted scale really does drift

    // ---- GPS run: the corrected trajectory ----------------------------------------------------
    bool outlier_rejected = false;
    int applied = 0;
    const Scalar gps_tail = run(/*with_gps=*/true, &outlier_rejected, &applied);

    CHECK(gps_tail < baseline_tail * 0.5);  // the payoff: GPS halves (or better) the tail drift
    CHECK(applied > 0);                     // fixes were actually applied
    CHECK(outlier_rejected);                // the gate rejected the gross fix (NIS > per-n gate)

    MESSAGE("baseline tail-mean err = " << baseline_tail << " m, GPS tail-mean err = " << gps_tail << " m");
}

// ---------------------------------------------------------------------------------------------
// 7. Innovation-adaptive robust R (GPS_R_NEES_SWEEP.md follow-up). A scalar cov_floor cannot be
//    cross-drive consistent; adaptive_r sets R from the robust (MAD) scale of recent innovations.
//    Default OFF == the cov_enu+cov_floor path (byte-identical); ON tracks the BULK spread, is
//    floored on a tight stream, and IMMUNE to a gross outlier (the multipath case). Identity pose
//    -> innovation == z == the fix ENU, so the fixes control the innovation directly.
// ---------------------------------------------------------------------------------------------
TEST_CASE("GPS adaptive-R: default off byte-identical; on = bulk spread, floored, outlier-immune") {
    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    State x;  // identity pose -> innovation == z == fix ENU

    auto feed = [&](GpsCorrection& g, double dlat, double dlon, double dalt, Vec3& rdiag) {
        g.submit_fix(make_fix(lat0 + dlat, lon0 + dlon, alt0 + dalt));
        Measurement m;
        REQUIRE(g.evaluate(x, m));
        rdiag = m.R.topLeftCorner<3, 3>().diagonal();
    };

    // --- default OFF: R == cov_enu(I) + cov_floor*I (the unchanged path) ---
    {
        GpsConfig c = datum_cfg(lat0, lon0, alt0);
        c.cov_floor_m2 = 4.0;                       // 1 + 4 = 5 on each axis
        GpsCorrection g(c);
        Vec3 r; feed(g, 1e-5, 0, 0, r);
        CHECK(r.x() == doctest::Approx(5.0).epsilon(1e-9));
        CHECK(r.y() == doctest::Approx(5.0).epsilon(1e-9));
        CHECK(r.z() == doctest::Approx(5.0).epsilon(1e-9));
    }

    // --- adaptive ON, TIGHT stream (identical fixes -> zero spread) -> R == floor ---
    {
        GpsConfig c = datum_cfg(lat0, lon0, alt0);
        c.adaptive_r = true; c.adaptive_min_samples = 10; c.adaptive_r_floor_m2 = 0.25;
        GpsCorrection g(c);
        Vec3 r; for (int k = 0; k < 14; ++k) feed(g, 1e-5, 1e-5, 0.0, r);
        CHECK(r.x() == doctest::Approx(0.25).epsilon(1e-6));
        CHECK(r.y() == doctest::Approx(0.25).epsilon(1e-6));
        CHECK(r.z() == doctest::Approx(0.25).epsilon(1e-6));
    }

    // --- adaptive ON, ONE gross outlier among a tight stream -> R STAYS floored (MAD-immune) ---
    {
        GpsConfig c = datum_cfg(lat0, lon0, alt0);
        c.adaptive_r = true; c.adaptive_min_samples = 10; c.adaptive_r_floor_m2 = 0.25;
        GpsCorrection g(c);
        Vec3 r;
        for (int k = 0; k < 13; ++k) feed(g, 1e-5, 1e-5, 0.0, r);  // tight bulk
        feed(g, 5e-3, 5e-3, 0.0, r);                                // gross multipath fix (~hundreds m)
        feed(g, 1e-5, 1e-5, 0.0, r);                                // R now uses the ring incl. the outlier
        CHECK(r.x() == doctest::Approx(0.25).epsilon(1e-6));        // median-immune -> still floored
        CHECK(r.y() == doctest::Approx(0.25).epsilon(1e-6));
    }

    // --- adaptive ON, a real bulk SPREAD -> R rises above the floor ---
    {
        GpsConfig c = datum_cfg(lat0, lon0, alt0);
        c.adaptive_r = true; c.adaptive_min_samples = 10; c.adaptive_r_floor_m2 = 0.01;
        GpsCorrection g(c);
        Vec3 r;
        for (int k = 0; k < 20; ++k) feed(g, 0.0, k * 1e-5, 0.0, r);  // East ramp 0..~14 m (unimodal spread)
        CHECK(r.x() > 0.5);                                          // East spread -> real R
        CHECK(r.y() == doctest::Approx(0.01).epsilon(1e-6));         // North flat -> floor
    }
}
