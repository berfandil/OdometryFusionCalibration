// Adapters Slice 13 (real-data ingestion): the generic-CSV ISource + GT track + replay/eval
// harness. These tests are SELF-CONTAINED — the adapters test exe links only ofc_core +
// ofc_adapters + doctest (NOT ofc_sim), so the equivalence rig uses adapter_test_util's
// TwistSource. No RNG (deterministic).
//
// Coverage:
//   * round-trip — write a known motion to a CSV (each of the 3 forms), read it back, assert
//     query(t0,t1) returns the expected Delta. Quaternion round-trips.
//   * EQUIVALENCE (the key test) — a deterministic constant-twist source's deltas written to a
//     CSV, run through an Estimator, vs the analytic TwistSource through an Estimator with the
//     SAME config: the published frontier poses match (CSV ingestion == in-memory path).
//   * GT eval — synth a GT track + an estimator run; assert the harness's drift + NEES are finite,
//     correctly signed, and a hand-computed NEES matches on a sampled step.
//   * GPS via CSV — GPS fixes fed through GpsCorrection reduce drift in the harness.
//   * Parse robustness — comments/blanks, missing optional cov, a malformed row -> Status+error.
#include <doctest/doctest.h>

#include "ofc_adapters/csv_source.hpp"
#include "ofc_adapters/replay_harness.hpp"

#include "ofc/core/eskf.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include "adapter_test_util.hpp"

#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace ofc;
using namespace ofc::adapters;
using namespace adptest;

namespace {

constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;

// w-first quaternion (qw,qx,qy,qz) from a rotation matrix (test-local; mirrors mat3_to_quat).
void R_to_q(const Mat3& R, Scalar q[4]) {
    Scalar w, x, y, z;
    mat3_to_quat(R, w, x, y, z);
    q[0] = w; q[1] = x; q[2] = y; q[3] = z;
}

// Build an `absolute`/`increment` CSV row.
std::string pose_row(Timestamp t, const SE3& p) {
    Scalar q[4]; R_to_q(p.R, q);
    std::ostringstream os;
    os.precision(15);
    os << t << ',' << p.t.x() << ',' << p.t.y() << ',' << p.t.z() << ','
       << q[0] << ',' << q[1] << ',' << q[2] << ',' << q[3] << '\n';
    return os.str();
}

bool se3_close(const SE3& a, const SE3& b, Scalar tol) {
    return (a.t - b.t).norm() < tol && so3::log(a.R.transpose() * b.R).norm() < tol;
}

} // namespace

// =================================================================================================
// 1. Round-trip: ABSOLUTE form. A constant-twist trajectory sampled as cumulative poses; reading
//    it back and querying [t0,t1] reproduces the SE(3) composition delta.
// =================================================================================================
TEST_CASE("CSV round-trip absolute: query(t0,t1) == pose(t0)^-1 o pose(t1)") {
    Vec6 xi; xi << 1.5, 0.0, 0.0, 0.0, 0.0, 0.4;   // forward + yaw
    const Scalar rate = 100.0;
    const Timestamp step = secs_to_ns(1.0 / rate);

    // Cumulative poses pose(k) = exp(xi * k*dt).
    std::string csv = "# form: absolute\n# t_ns,x,y,z,qw,qx,qy,qz\n";
    std::vector<SE3> poses;
    SE3 cum;  // identity
    const int N = 50;
    for (int k = 0; k < N; ++k) {
        const Timestamp t = static_cast<Timestamp>(k) * step;
        if (k > 0) cum = se3::compose(cum, se3::exp(xi * (1.0 / rate)));
        poses.push_back(cum);
        csv += pose_row(t, cum);
    }

    CsvSource src;
    CsvSourceConfig cfg;
    cfg.id = 1;
    cfg.combine = ConfCombine::NativeOnly;   // no modeled inflation -> exact geometry
    REQUIRE(src.load(csv, cfg) == Status::Ok);
    CHECK(src.form() == OdomForm::AbsolutePose);
    CHECK(src.row_count() == N);

    // Query a window aligned to samples 10..30.
    const Timestamp t0 = 10 * step, t1 = 30 * step;
    Expected<Delta> d = src.query(t0, t1);
    REQUIRE(d.ok());
    const SE3 expect = se3::compose(se3::inverse(poses[10]), poses[30]);
    CHECK(se3_close(d.value().motion, expect, 1e-9));
}

// =================================================================================================
// 2. Round-trip: INCREMENT form. Per-step relative motion; query integrates them back.
// =================================================================================================
TEST_CASE("CSV round-trip increment: query integrates the per-step increments") {
    Vec6 xi; xi << 2.0, 0.0, 0.0, 0.0, 0.1, 0.3;
    const Scalar rate = 100.0;
    const Timestamp step = secs_to_ns(1.0 / rate);
    const SE3 incr = se3::exp(xi * (1.0 / rate));   // constant per-step increment

    std::string csv = "# form: increment\n";
    SE3 cum;
    std::vector<SE3> poses;
    const int N = 40;
    for (int k = 0; k < N; ++k) {
        const Timestamp t = static_cast<Timestamp>(k) * step;
        // First row is the identity seed; subsequent rows carry `incr`.
        const SE3 row_pose = (k == 0) ? SE3{} : incr;
        csv += pose_row(t, row_pose);
        if (k > 0) cum = se3::compose(cum, incr);
        poses.push_back(cum);
    }

    CsvSource src;
    CsvSourceConfig cfg; cfg.id = 2; cfg.combine = ConfCombine::NativeOnly;
    REQUIRE(src.load(csv, cfg) == Status::Ok);
    CHECK(src.form() == OdomForm::Increment);

    const Timestamp t0 = 5 * step, t1 = 25 * step;
    Expected<Delta> d = src.query(t0, t1);
    REQUIRE(d.ok());
    const SE3 expect = se3::compose(se3::inverse(poses[5]), poses[25]);
    CHECK(se3_close(d.value().motion, expect, 1e-9));
}

// =================================================================================================
// 3. Round-trip: TWIST form. Constant body twist; query integrates exp(xi*dt) over the window.
// =================================================================================================
TEST_CASE("CSV round-trip twist: query == exp(xi * (t1-t0))") {
    Vec6 xi; xi << 1.0, 0.2, 0.0, 0.0, 0.0, 0.5;
    const Scalar rate = 200.0;
    const Timestamp step = secs_to_ns(1.0 / rate);

    std::ostringstream os; os.precision(15);
    os << "# form: twist\n";
    const int N = 60;
    for (int k = 0; k < N; ++k) {
        const Timestamp t = static_cast<Timestamp>(k) * step;
        os << t << ',' << xi[0] << ',' << xi[1] << ',' << xi[2] << ','
           << xi[3] << ',' << xi[4] << ',' << xi[5] << '\n';
    }

    CsvSource src;
    CsvSourceConfig cfg; cfg.id = 3; cfg.combine = ConfCombine::NativeOnly;
    REQUIRE(src.load(os.str(), cfg) == Status::Ok);
    CHECK(src.form() == OdomForm::Twist);

    const Timestamp t0 = 10 * step, t1 = 40 * step;
    Expected<Delta> d = src.query(t0, t1);
    REQUIRE(d.ok());
    const Scalar dt = static_cast<Scalar>(t1 - t0) / 1e9;
    const SE3 expect = se3::exp(xi * dt);
    CHECK(se3_close(d.value().motion, expect, 1e-7));
}

// =================================================================================================
// 4. Quaternion round-trips through the writer + reader for a non-trivial rotation.
// =================================================================================================
TEST_CASE("CSV quaternion round-trip: a non-trivial rotation survives write -> read") {
    SE3 p; p.R = so3::exp(Vec3(0.3, -0.7, 1.2)); p.t = Vec3(3.0, -1.0, 0.5);
    std::string csv = "# form: absolute\n";
    csv += pose_row(0, SE3{});           // identity seed
    csv += pose_row(secs_to_ns(0.01), p);

    CsvSource src;
    CsvSourceConfig cfg; cfg.id = 4; cfg.combine = ConfCombine::NativeOnly;
    REQUIRE(src.load(csv, cfg) == Status::Ok);
    Expected<Delta> d = src.query(0, secs_to_ns(0.01));
    REQUIRE(d.ok());
    // Delta over [0, 0.01] = identity^-1 o p == p.
    CHECK(se3_close(d.value().motion, p, 1e-9));
}

// =================================================================================================
// 5. EQUIVALENCE (the key test): CSV ingestion == in-memory path through the Estimator.
//    A constant-twist source's deltas are written to an INCREMENT CSV, then both (a) the analytic
//    TwistSource and (b) the CsvSource are run through Estimators with the SAME Config; the
//    published frontier poses match step-for-step.
// =================================================================================================
TEST_CASE("CSV equivalence: CsvSource through the Estimator matches the analytic in-memory path") {
    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_cfg(sensors);   // 2-source rig, ref 0 + src 1, timesync off

    Vec6 xi; xi << 2.0, 0.0, 0.0, 0.0, 0.0, 0.3;   // forward + yaw (matches make_rig_sources)

    // The timeline: 50 Hz steps over ~2 s.
    std::vector<Timestamp> ts;
    const Timestamp tick = secs_to_ns(0.02);
    for (Timestamp now = secs_to_ns(0.2); now <= secs_to_ns(2.0); now += tick) ts.push_back(now);

    // --- (a) analytic in-memory run ------------------------------------------------------------
    std::map<Timestamp, SE3> analytic;
    {
        std::vector<TwistSource> srcs;
        srcs.emplace_back(0, xi, Scalar(1.0));
        srcs.emplace_back(1, xi, Scalar(1.0));
        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);
        for (Timestamp now : ts) {
            if (est.step(now) != Status::Ok) continue;
            const Result& r = est.latest();
            analytic[r.frontier.stamp] = r.frontier.pose;
        }
    }

    // --- write the SAME motion to INCREMENT CSVs (sampled at 200 Hz so windows are well covered) -
    auto build_csv = [&](SourceId id) {
        const Scalar rate = 200.0;
        const Timestamp step = secs_to_ns(1.0 / rate);
        const SE3 incr = se3::exp(xi * (1.0 / rate));
        std::string csv = "# form: increment\n";
        // Cover the whole timeline (a touch before 0 to a touch after 2 s) so every query window
        // is bracketed.
        const int N = static_cast<int>(2.1 * rate) + 2;
        for (int k = 0; k < N; ++k) {
            const Timestamp t = static_cast<Timestamp>(k) * step;
            csv += pose_row(t, (k == 0) ? SE3{} : incr);
        }
        (void)id;
        return csv;
    };

    // --- (b) CSV-backed run --------------------------------------------------------------------
    std::map<Timestamp, SE3> from_csv;
    {
        CsvSource s0, s1;
        CsvSourceConfig c0; c0.id = 0; c0.combine = ConfCombine::NativeOnly; c0.native_confidence = true;
        CsvSourceConfig c1 = c0; c1.id = 1;
        REQUIRE(s0.load(build_csv(0), c0) == Status::Ok);
        REQUIRE(s1.load(build_csv(1), c1) == Status::Ok);
        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        REQUIRE(est.add_source(&s0) == Status::Ok);
        REQUIRE(est.add_source(&s1) == Status::Ok);
        for (Timestamp now : ts) {
            if (est.step(now) != Status::Ok) continue;
            const Result& r = est.latest();
            from_csv[r.frontier.stamp] = r.frontier.pose;
        }
    }

    // --- compare published frontier poses step-for-step ----------------------------------------
    REQUIRE(!analytic.empty());
    REQUIRE(analytic.size() == from_csv.size());
    int compared = 0;
    for (const auto& kv : analytic) {
        auto it = from_csv.find(kv.first);
        REQUIRE(it != from_csv.end());
        // The CSV path samples at 200 Hz and the SourceBuffer interpolates the window endpoints;
        // the analytic path is exact, so allow a tight (sub-mm / sub-mrad) interpolation tolerance.
        CHECK(se3_close(kv.second, it->second, 1e-6));
        ++compared;
    }
    CHECK(compared > 50);
}

// =================================================================================================
// 6. GT eval: the harness computes finite, correctly-signed drift + a hand-checked NEES.
//    A perfect-tracking run (CSV truth == GT) should read ~0 drift; a known offset reads a known
//    error; the NEES on a sampled step matches the documented se3::log convention by hand.
// =================================================================================================
TEST_CASE("Replay GT eval: drift + 6-DOF NEES are finite, correctly signed, hand-checked") {
    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_cfg(sensors);
    cfg.q_floor[0] = cfg.q_floor[1] = cfg.q_floor[2] = 1e-4;   // a little position uncertainty

    Vec6 xi; xi << 2.0, 0.0, 0.0, 0.0, 0.0, 0.3;
    const Scalar rate = 200.0;
    const Timestamp step = secs_to_ns(1.0 / rate);
    const SE3 incr = se3::exp(xi * (1.0 / rate));

    // Two identical increment sources (no drift) + a GT track that is the SAME constant-twist
    // trajectory expressed as cumulative absolute poses (so the fused frontier should track GT).
    auto inc_csv = [&]() {
        std::string csv = "# form: increment\n";
        const int N = static_cast<int>(2.1 * rate) + 2;
        for (int k = 0; k < N; ++k)
            csv += pose_row(static_cast<Timestamp>(k) * step, (k == 0) ? SE3{} : incr);
        return csv;
    };
    std::string gt_csv = "# GT track t_ns,x,y,z,qw,qx,qy,qz\n";
    {
        SE3 cum;
        const int N = static_cast<int>(2.1 * rate) + 2;
        for (int k = 0; k < N; ++k) {
            const Timestamp t = static_cast<Timestamp>(k) * step;
            if (k > 0) cum = se3::compose(cum, incr);
            gt_csv += pose_row(t, cum);
        }
    }

    CsvSource s0, s1;
    CsvSourceConfig c0; c0.id = 0; c0.combine = ConfCombine::NativeOnly;
    CsvSourceConfig c1 = c0; c1.id = 1;
    REQUIRE(s0.load(inc_csv(), c0) == Status::Ok);
    REQUIRE(s1.load(inc_csv(), c1) == Status::Ok);
    CsvGtTrack gt;
    REQUIRE(gt.load(gt_csv) == Status::Ok);
    CHECK(gt.size() > 100);

    ReplayInputs in;
    in.cfg = &cfg;
    in.sources = {&s0, &s1};
    in.gt = &gt;
    ReplayHarness h;
    REQUIRE(h.run(in) == Status::Ok);

    const RunSummary& s = h.summary();
    CHECK(s.fused_steps > 50);
    CHECK(s.has_gt);
    // Perfect tracking: drift is small and finite (the GT trajectory == the fused odometry).
    CHECK(std::isfinite(s.max_trans_m));
    CHECK(std::isfinite(s.mean_pose_nees));
    CHECK(s.max_trans_m < 0.05);     // sub-5cm: the source IS the GT motion (interp residual only)
    CHECK(s.max_rot_rad < 0.02);
    CHECK(s.nees_count > 0);
    CHECK(s.mean_pose_nees > 0.0);   // a genuine (positive-definite) NEES

    // Hand-check the NEES convention on one fused-with-GT record (the documented se3::log form).
    const std::vector<ReplayRecord>& recs = h.records();
    int checked = 0;
    for (const ReplayRecord& r : recs) {
        if (!r.fused || !r.has_gt) continue;
        const SE3 errT = se3::compose(se3::inverse(r.result.frontier.pose), r.gt_pose);
        const Vec6 e   = se3::log(errT);
        const Mat6 Ppp = r.result.frontier.cov.block<6, 6>(0, 0);
        const Scalar hand = e.dot(Ppp.ldlt().solve(e));
        CHECK(r.pose_nees == doctest::Approx(hand).epsilon(1e-9));
        // trans/rot error signs match the recorded magnitudes.
        CHECK(r.trans_err_m == doctest::Approx((r.result.frontier.pose.t - r.gt_pose.t).norm()));
        if (++checked >= 5) break;
    }
    CHECK(checked > 0);
}

// =================================================================================================
// 6b. Local (GT-anchored fixed-window) relative-pose error: a constant per-step over-scale makes
//     the GLOBAL drift grow without bound with run length, but each fixed-length window re-anchors
//     to GT at its start, so the local metric measures only intra-window drift -> it is bounded and
//     length-independent (the point of the metric). Also asserts it is OFF by default.
// =================================================================================================
TEST_CASE("Replay local windows: GT-anchored relative-pose error is bounded vs the global drift") {
    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_cfg(sensors);

    Vec6 xi; xi << 2.0, 0.0, 0.0, 0.0, 0.0, 0.3;
    const Scalar rate = 200.0;
    const Timestamp step = secs_to_ns(1.0 / rate);
    const Scalar planted_scale = 1.05;                  // 5% over-scale -> linearly growing global drift
    const int N = static_cast<int>(4.1 * rate) + 2;     // long run so global drift >> a single window

    const SE3 incr = se3::exp(xi * (1.0 / rate));
    SE3 incr_scaled = incr; incr_scaled.t *= planted_scale;

    // Two identical over-scale increment sources.
    auto inc_csv = [&]() {
        std::string csv = "# form: increment\n";
        for (int k = 0; k < N; ++k)
            csv += pose_row(static_cast<Timestamp>(k) * step, (k == 0) ? SE3{} : incr_scaled);
        return csv;
    };
    // GT = the TRUE (scale=1) cumulative trajectory.
    std::string gt_csv = "# GT t_ns,x,y,z,qw,qx,qy,qz\n";
    {
        SE3 cum;
        for (int k = 0; k < N; ++k) {
            const Timestamp t = static_cast<Timestamp>(k) * step;
            if (k > 0) cum = se3::compose(cum, incr);
            gt_csv += pose_row(t, cum);
        }
    }

    CsvSource s0, s1;
    CsvSourceConfig c0; c0.id = 0; c0.combine = ConfCombine::NativeOnly;
    CsvSourceConfig c1 = c0; c1.id = 1;
    REQUIRE(s0.load(inc_csv(), c0) == Status::Ok);
    REQUIRE(s1.load(inc_csv(), c1) == Status::Ok);
    CsvGtTrack gt; REQUIRE(gt.load(gt_csv) == Status::Ok);

    // (a) OFF by default (local_batch_len == 0 -> no windows).
    {
        ReplayInputs in; in.cfg = &cfg; in.sources = {&s0, &s1}; in.gt = &gt;
        ReplayHarness h; REQUIRE(h.run(in) == Status::Ok);
        CHECK(h.summary().local_batch_len == 0);
        CHECK(h.summary().local_batch_count == 0);
    }
    // (b) ON: windows measured; errors finite & positive; the local error is far smaller than the
    //     global accumulated drift (length-normalization).
    {
        const int L = 30;   // harness ticks at cfg.tick_rate_hz (50) -> ~180 post-warmup records
        ReplayInputs in; in.cfg = &cfg; in.sources = {&s0, &s1}; in.gt = &gt; in.local_batch_len = L;
        ReplayHarness h; REQUIRE(h.run(in) == Status::Ok);
        const RunSummary& s = h.summary();
        CHECK(s.local_batch_len == L);
        CHECK(s.local_batch_count >= 3);
        CHECK(std::isfinite(s.local_mean_trans_m));
        CHECK(s.local_mean_trans_m > 0.0);          // the over-scale leaks into every window
        CHECK(s.local_med_trans_m  > 0.0);
        // The whole-run drift accumulates far more than any single GT-anchored window.
        CHECK(s.local_max_trans_m < s.max_trans_m);
    }
}

// =================================================================================================
// 7. GPS via CSV-driven harness: GPS fixes reduce the drift of a scale-mis-calibrated rig
//    (mirrors the GPS e2e test, sourced through the harness). The drift-free truth is a scale=1
//    twin run; the GPS targets the twin's frontier position. With GPS the tail drift halves.
// =================================================================================================
TEST_CASE("Replay GPS via CSV: GPS fixes reduce the harness drift vs a no-GPS baseline") {
    std::vector<SensorConfig> sensors;
    Config cfg = make_rig_cfg(sensors);
    cfg.mahalanobis_chi2 = 100.0;                              // test artifact (see GPS e2e test)
    cfg.q_floor[0] = cfg.q_floor[1] = cfg.q_floor[2] = 1e-3;   // position uncertainty -> GPS pulls

    Vec6 xi; xi << 2.0, 0.0, 0.0, 0.0, 0.0, 0.3;
    const Scalar planted_scale = 1.04;   // 4% over-scale -> growing drift
    const Scalar rate = 200.0;
    const Timestamp step = secs_to_ns(1.0 / rate);
    const int Nrows = static_cast<int>(2.1 * rate) + 2;

    // Increment CSV with a translation scale applied (the planted drift).
    auto inc_csv = [&](Scalar scale) {
        SE3 incr = se3::exp(xi * (1.0 / rate));
        incr.t *= scale;
        std::string csv = "# form: increment\n";
        for (int k = 0; k < Nrows; ++k)
            csv += pose_row(static_cast<Timestamp>(k) * step, (k == 0) ? SE3{} : incr);
        return csv;
    };

    const double lat0 = 47.0, lon0 = 8.0, alt0 = 400.0;
    auto datum_cfg = [&]() {
        GpsConfig c; c.has_datum = true;
        c.datum_lat_deg = lat0; c.datum_lon_deg = lon0; c.datum_alt_m = alt0;
        return c;
    };
    // Inverse local-tangent: an ENU target (== odom position, identity alignment) -> a geodetic fix.
    auto enu_to_fix = [&](const Vec3& enu, Timestamp stamp, Scalar sigma) {
        const double a = 6378137.0;
        GpsFix f;
        f.lat_deg = lat0 + (enu.y() / a) / kDeg2Rad;
        f.lon_deg = lon0 + (enu.x() / (a * std::cos(lat0 * kDeg2Rad))) / kDeg2Rad;
        f.alt_m   = alt0 + enu.z();
        f.cov_enu = Mat3::Identity() * (sigma * sigma);
        f.stamp   = stamp;
        f.valid   = true;
        return f;
    };

    // Run the harness on a `scale` rig, optionally GPS-corrected. Returns the tail-mean position
    // error vs a scale=1 twin run (the same odom frame/anchor). The twin's per-frontier position
    // is captured first, then reused as the GPS target.
    auto run = [&](Scalar scale, bool with_gps, std::map<Timestamp, Vec3>* capture_truth,
                   const std::map<Timestamp, Vec3>* truth, bool* applied) -> Scalar {
        CsvSource a, b;
        CsvSourceConfig ca; ca.id = 0; ca.combine = ConfCombine::NativeOnly;
        CsvSourceConfig cb = ca; cb.id = 1;
        REQUIRE(a.load(inc_csv(scale), ca) == Status::Ok);
        REQUIRE(b.load(inc_csv(scale), cb) == Status::Ok);

        GpsCorrection gps(datum_cfg());

        ReplayInputs in;
        in.cfg = &cfg;
        in.sources = {&a, &b};
        if (with_gps) {
            in.gps = &gps;
            // Build fixes against the truth twin at each truth frontier stamp.
            if (truth) {
                for (const auto& kv : *truth)
                    in.gps_fixes.push_back(enu_to_fix(kv.second, kv.first, Scalar(0.05)));
            }
        }
        ReplayHarness h;
        REQUIRE(h.run(in) == Status::Ok);
        if (applied) *applied = h.summary().gps_applied > 0;

        // Capture / compare the frontier positions vs the truth twin.
        Scalar tail_sum = 0; int tail_n = 0;
        const Timestamp last = h.records().empty() ? 0 : h.records().back().now;
        const Timestamp tail_start = last - secs_to_ns(0.5);
        for (const ReplayRecord& r : h.records()) {
            if (!r.fused) continue;
            if (capture_truth) (*capture_truth)[r.result.frontier.stamp] = r.result.frontier.pose.t;
            if (truth) {
                auto it = truth->find(r.result.frontier.stamp);
                if (it != truth->end() && r.now >= tail_start) {
                    tail_sum += (r.result.frontier.pose.t - it->second).norm();
                    ++tail_n;
                }
            }
        }
        return (tail_n > 0) ? tail_sum / tail_n : Scalar(0);
    };

    // 1) capture the drift-free truth (scale=1).
    std::map<Timestamp, Vec3> truth;
    run(Scalar(1.0), /*with_gps=*/false, &truth, nullptr, nullptr);
    REQUIRE(!truth.empty());

    // 2) baseline drift (scale=1.04, no GPS) vs the truth.
    const Scalar baseline = run(planted_scale, /*with_gps=*/false, nullptr, &truth, nullptr);
    CHECK(baseline > 0.05);   // the planted scale really does drift

    // 3) GPS-corrected drift.
    bool applied = false;
    const Scalar with_gps = run(planted_scale, /*with_gps=*/true, nullptr, &truth, &applied);
    CHECK(applied);
    CHECK(with_gps < baseline * 0.5);   // GPS halves (or better) the tail drift
    MESSAGE("harness baseline tail=" << baseline << " m, GPS tail=" << with_gps << " m");
}

// =================================================================================================
// 8. Parse robustness: comments / blank lines / a header row / missing optional cov all parse;
//    a malformed row -> a clean Status + a human error string (no crash).
// =================================================================================================
TEST_CASE("CSV parse robustness: comments/blank/header/optional-cov ok; malformed -> Status+error") {
    // (a) comments, blanks, a header row, mixed delimiters, and the optional 6-cov on some rows.
    std::string good =
        "# form: increment\n"
        "// a slash comment\n"
        "\n"
        "t_ns x y z qw qx qy qz\n"                       // header row (non-numeric first field)
        "0, 0,0,0, 1,0,0,0\n"                            // identity seed (no cov)
        "10000000 0.1 0 0 1 0 0 0 1e-4 1e-4 1e-4 1e-5 1e-5 1e-5\n"   // whitespace + 6 cov
        "20000000,0.1,0,0,1,0,0,0\n";                    // comma, no cov
    CsvSource src;
    CsvSourceConfig cfg; cfg.id = 0;
    CHECK(src.load(good, cfg) == Status::Ok);
    CHECK(src.row_count() == 3);

    // (b) a malformed quaternion (zero-norm) -> InvalidConfig + a non-empty error.
    std::string bad_q =
        "# form: absolute\n"
        "0,0,0,0,1,0,0,0\n"
        "10000000,1,0,0,0,0,0,0\n";   // qw=qx=qy=qz=0 -> degenerate
    CsvSource s2;
    CHECK(s2.load(bad_q, cfg) == Status::InvalidConfig);
    CHECK_FALSE(s2.error().empty());

    // (c) wrong column count -> InvalidConfig.
    std::string bad_cols =
        "# form: absolute\n"
        "0,0,0,0,1,0,0,0\n"
        "10000000,1,0,0,1,0,0\n";     // only 7 columns
    CsvSource s3;
    CHECK(s3.load(bad_cols, cfg) == Status::InvalidConfig);
    CHECK_FALSE(s3.error().empty());

    // (d) non-increasing timestamp -> InvalidConfig.
    std::string bad_ts =
        "# form: increment\n"
        "0,0,0,0,1,0,0,0\n"
        "0,0.1,0,0,1,0,0,0\n";        // duplicate stamp
    CsvSource s4;
    CHECK(s4.load(bad_ts, cfg) == Status::InvalidConfig);
    CHECK_FALSE(s4.error().empty());

    // (e) too few rows -> NoData.
    std::string tiny = "# form: increment\n0,0,0,0,1,0,0,0\n";
    CsvSource s5;
    CHECK(s5.load(tiny, cfg) == Status::NoData);

    // (f) GT track robustness: a malformed row -> InvalidConfig.
    CsvGtTrack gt;
    CHECK(gt.load("0,0,0,0,1,0,0,0\n5,1,0,0\n") == Status::InvalidConfig);   // 2nd row too short
    CHECK_FALSE(gt.error().empty());
    // valid GT.
    CsvGtTrack gt2;
    CHECK(gt2.load("0,0,0,0,1,0,0,0\n5,1,0,0,1,0,0,0\n") == Status::Ok);
    CHECK(gt2.size() == 2);
}

// =================================================================================================
// 9. The forced ctor form overrides a tag; a forced form disagreeing with a tag is an error.
// =================================================================================================
TEST_CASE("CSV form selection: ctor force_form overrides the tag; a disagreement is an error") {
    // Tag says twist, but we force increment -> the rows must be 8-col pose rows.
    std::string csv = "# form: increment\n0,0,0,0,1,0,0,0\n10000000,0.1,0,0,1,0,0,0\n";
    CsvSource s;
    CsvSourceConfig cfg; cfg.id = 0; cfg.form = OdomForm::Increment; cfg.force_form = true;
    CHECK(s.load(csv, cfg) == Status::Ok);
    CHECK(s.form() == OdomForm::Increment);

    // Forced increment but the tag says twist -> disagreement -> InvalidConfig.
    std::string csv2 = "# form: twist\n0,0,0,0,0,0,0\n10000000,0.1,0,0,0,0,0\n";
    CsvSource s2;
    CsvSourceConfig cfg2; cfg2.id = 0; cfg2.form = OdomForm::Increment; cfg2.force_form = true;
    CHECK(s2.load(csv2, cfg2) == Status::InvalidConfig);
    CHECK_FALSE(s2.error().empty());
}
