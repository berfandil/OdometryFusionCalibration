// Slice 5 tests: TimeSync — ‖ω‖ cross-correlation -> per-source offset histogram.
//
// Unit: two synthetic ‖ω‖ signals, one shifted by a known lag (integer AND fractional)
//   -> xcorr recovers it to sub-sample accuracy under each of L1/L2/Ratio/NCC; a flat /
//   low-variance signal is rejected by the excitation gate (no confident estimate).
// Sim e2e (the oracle): omega_varying() + a source with a planted time_offset_s
//   -> TimeSync recovers it within sub-sample tol AND with the CORRECT SIGN (positive
//   planted recovered positive). A straight() (no rotation) scenario -> gated, no
//   confident estimate. Determinism (same input -> identical estimate).
#include <doctest/doctest.h>

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/timesync.hpp"

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

// Absolute-tolerance comparison (doctest 2.4.11's Approx has no .margin()).
bool near_abs(Scalar a, Scalar b, Scalar tol) { return std::abs(a - b) <= tol; }

// A Config with a sensible offset histogram for time-sync (fine bins over the lag range).
Config make_cfg(MatchMetric metric, Scalar tick_hz = 100.0, Scalar max_lag_s = 0.1) {
    Config c;
    c.tick_rate_hz     = tick_hz;
    c.match_metric     = metric;
    c.max_lag_s        = max_lag_s;
    c.timesync_enabled = true;
    // Offset histogram: linear over [-max_lag, +max_lag], many bins so the sub-bin
    // mode resolves a few-ms offset. SlidingK so a few votes already give a sharp peak.
    c.offset_hist.bins        = 256;
    c.offset_hist.range_min   = -max_lag_s;
    c.offset_hist.range_max   =  max_lag_s;
    c.offset_hist.circular    = false;
    c.offset_hist.aging       = Aging::SlidingK;
    c.offset_hist.sliding_k   = 64;
    c.offset_hist.vote_split  = true;
    c.offset_hist.subbin      = true;
    return c;
}

// TimeSync holds large fixed-capacity buffers (the strict-core no-heap-after-init
// contract sizes them from compile-time maxima). The core keeps it inside the
// heap-allocated Estimator::Impl; here in the relaxed test edge we likewise heap-
// allocate it (a stack instance would overflow the default 1 MB test stack).
std::unique_ptr<TimeSync> make_ts() { return std::unique_ptr<TimeSync>(new TimeSync()); }

// Push a sampled signal (one value per grid tick starting at t=0) into a channel id.
void push_signal(TimeSync& ts, SourceId id, const std::vector<Scalar>& sig, Scalar dt_s) {
    for (std::size_t k = 0; k < sig.size(); ++k) {
        ts.push(id, secs(static_cast<double>(k) * dt_s), sig[k]);
    }
}

// A distinctive ‖ω‖-like signal: a couple of Gaussian bumps (smooth, well-excited).
Scalar bump(Scalar t, Scalar c, Scalar w) {
    const Scalar x = (t - c) / w;
    return std::exp(-0.5 * x * x);
}
Scalar shape(Scalar t) {
    return 0.2 + 1.0 * bump(t, 1.0, 0.25) + 0.6 * bump(t, 2.2, 0.18)
               + 0.9 * bump(t, 3.0, 0.30);
}
} // namespace

// ===========================================================================
// configure() validation
// ===========================================================================
TEST_CASE("timesync configure validates and derives the grid") {
    auto tsp = make_ts(); TimeSync& ts = *tsp;
    Config c = make_cfg(MatchMetric::L2, 100.0, 0.1);
    REQUIRE(ts.configure(c, /*reference=*/0) == Status::Ok);
    CHECK(ts.configured());
    CHECK(ts.reference() == 0);
    CHECK(ts.sample_dt() == doctest::Approx(0.01));        // 1/100 Hz
    CHECK(ts.max_lag_samples() == 10);                     // 0.1 s * 100 Hz

    Config bad_tick = c; bad_tick.tick_rate_hz = 0.0;
    CHECK(ts.configure(bad_tick, 0) == Status::OutOfRange);

    Config bad_lag = c; bad_lag.max_lag_s = 0.0;
    CHECK(ts.configure(bad_lag, 0) == Status::OutOfRange);
}

// ===========================================================================
// Unit: known-lag recovery under each metric (sub-sample accuracy)
// ===========================================================================
TEST_CASE("timesync recovers a known integer lag under every metric") {
    const Scalar dt = 0.01;            // 100 Hz
    const int    N  = 400;             // 4 s of signal
    const int    planted_lag = 3;      // grid steps: source leads base by 3 samples

    // base/reference signal omega_ref(t) = shape(t).
    // source signal omega_src(t) = omega_ref(t + off) with off = planted_lag*dt.
    // (source observes the base signal ADVANCED by off — the canonical sign.)
    const Scalar off = planted_lag * dt;
    std::vector<Scalar> ref(N), src(N);
    for (int k = 0; k < N; ++k) {
        const Scalar t = k * dt;
        ref[k] = shape(t);
        src[k] = shape(t + off);
    }

    const MatchMetric metrics[4] = {MatchMetric::L1, MatchMetric::L2,
                                    MatchMetric::Ratio, MatchMetric::NCC};
    for (MatchMetric m : metrics) {
        auto tsp = make_ts(); TimeSync& ts = *tsp;
        REQUIRE(ts.configure(make_cfg(m, 100.0, 0.1), 0) == Status::Ok);
        push_signal(ts, /*ref id*/ 0, ref, dt);
        push_signal(ts, /*src id*/ 1, src, dt);
        ts.update();
        const Scalar est = ts.offset(1);
        INFO("metric=", static_cast<int>(m), " est=", est, " want=", off);
        CHECK(near_abs(est, off, 0.5 * dt));
        CHECK(ts.confidence(1) > 0.0);
        // The reference itself reports no offset / no confidence.
        CHECK(ts.offset(0) == doctest::Approx(0.0));
        CHECK(ts.confidence(0) == doctest::Approx(0.0));
    }
}

TEST_CASE("timesync resolves a FRACTIONAL lag to sub-sample accuracy") {
    const Scalar dt = 0.01;
    const int    N  = 400;
    const Scalar off = 0.027;          // 2.7 grid steps — needs the parabolic refine
    std::vector<Scalar> ref(N), src(N);
    for (int k = 0; k < N; ++k) {
        const Scalar t = k * dt;
        ref[k] = shape(t);
        src[k] = shape(t + off);
    }
    // L2 + NCC are the smooth metrics; both should sub-sample resolve.
    for (MatchMetric m : {MatchMetric::L2, MatchMetric::NCC}) {
        auto tsp = make_ts(); TimeSync& ts = *tsp;
        REQUIRE(ts.configure(make_cfg(m, 100.0, 0.1), 0) == Status::Ok);
        push_signal(ts, 0, ref, dt);
        push_signal(ts, 1, src, dt);
        ts.update();
        const Scalar est = ts.offset(1);
        INFO("metric=", static_cast<int>(m), " est=", est, " want=", off);
        CHECK(near_abs(est, off, 0.4 * dt));
    }
}

TEST_CASE("timesync recovers a NEGATIVE lag (source clock behind base)") {
    const Scalar dt = 0.01;
    const int    N  = 400;
    const Scalar off = -0.04;          // source reads an EARLIER base slice
    std::vector<Scalar> ref(N), src(N);
    for (int k = 0; k < N; ++k) {
        const Scalar t = k * dt;
        ref[k] = shape(t);
        src[k] = shape(t + off);
    }
    auto tsp = make_ts(); TimeSync& ts = *tsp;
    REQUIRE(ts.configure(make_cfg(MatchMetric::L2, 100.0, 0.1), 0) == Status::Ok);
    push_signal(ts, 0, ref, dt);
    push_signal(ts, 1, src, dt);
    ts.update();
    CHECK(near_abs(ts.offset(1), off, 0.5 * dt));
}

// ===========================================================================
// Unit: excitation gate
// ===========================================================================
TEST_CASE("timesync rejects a flat / low-variance signal (excitation gate)") {
    const Scalar dt = 0.01;
    const int    N  = 400;
    std::vector<Scalar> flat(N, 0.0);   // no rotation at all -> zero variance
    std::vector<Scalar> ref = flat;
    std::vector<Scalar> src = flat;

    auto tsp = make_ts(); TimeSync& ts = *tsp;
    Config c = make_cfg(MatchMetric::L2, 100.0, 0.1);
    c.excitation_min_var = 1e-3;        // a real gate
    REQUIRE(ts.configure(c, 0) == Status::Ok);
    push_signal(ts, 0, ref, dt);
    push_signal(ts, 1, src, dt);
    ts.update();
    // No vote was cast -> the histogram is empty -> no confident estimate.
    CHECK(ts.confidence(1) == doctest::Approx(0.0));
    CHECK(ts.offset(1) == doctest::Approx(0.0));
}

TEST_CASE("timesync gate lets an excited signal through") {
    const Scalar dt = 0.01;
    const int    N  = 400;
    const int    planted_lag = 2;
    const Scalar off = planted_lag * dt;
    std::vector<Scalar> ref(N), src(N);
    for (int k = 0; k < N; ++k) {
        ref[k] = shape(k * dt);
        src[k] = shape(k * dt + off);
    }
    auto tsp = make_ts(); TimeSync& ts = *tsp;
    Config c = make_cfg(MatchMetric::L2, 100.0, 0.1);
    c.excitation_min_var = 1e-3;
    REQUIRE(ts.configure(c, 0) == Status::Ok);
    push_signal(ts, 0, ref, dt);
    push_signal(ts, 1, src, dt);
    ts.update();
    CHECK(ts.confidence(1) > 0.0);
    CHECK(near_abs(ts.offset(1), off, 0.5 * dt));
}

// ===========================================================================
// Sim end-to-end (the oracle): planted offset recovered with the correct sign
// ===========================================================================
namespace {
// Sample a source's ‖ω‖ from its reported SE(3) delta over a tiny sub-interval
// ending at the (base-clock) sample time t: ‖ω‖ ≈ ‖log(ΔR)‖ / h.
Scalar omega_norm_at(const SyntheticSource& s, double t_s, double h_s) {
    const Expected<Delta> q = s.query(secs(t_s - h_s), secs(t_s));
    if (!q.ok()) return 0.0;
    const Vec3 w = so3::log(q.value().motion.R);
    return w.norm() / h_s;
}
} // namespace

TEST_CASE("timesync sim e2e: recovers a planted +0.03 s offset with the CORRECT sign") {
    Trajectory tr = Trajectory::omega_varying();   // distinctive ‖ω‖(t) shape

    // Reference source: aligned, no offset.
    SourceParams pr; pr.id = 0; pr.scale = 1.0; pr.time_offset_s = 0.0;
    SyntheticSource ref(tr, pr);

    // Source 1: a planted POSITIVE offset (source clock ahead of base by 0.03 s).
    const Scalar planted = 0.03;
    SourceParams ps; ps.id = 1; ps.scale = 1.0; ps.time_offset_s = planted;
    SyntheticSource src(tr, ps);

    const Scalar tick_hz = 100.0;
    const Scalar dt = 1.0 / tick_hz;
    const Scalar h  = dt;                          // sub-interval for ‖ω‖

    auto tsp = make_ts(); TimeSync& ts = *tsp;
    REQUIRE(ts.configure(make_cfg(MatchMetric::L2, tick_hz, 0.1), 0) == Status::Ok);

    // Sample both sources across the trajectory at the common rate (base clock).
    const double t0 = tr.t0_s() + 2 * h;           // leave room for the sub-interval
    const double t1 = tr.end_s() - 0.05;
    for (double t = t0; t <= t1; t += dt) {
        ts.push(0, secs(t), omega_norm_at(ref, t, h));
        ts.push(1, secs(t), omega_norm_at(src, t, h));
    }
    ts.update();

    const Scalar est = ts.offset(1);
    INFO("planted=", planted, " recovered=", est);
    CHECK(ts.confidence(1) > 0.0);
    CHECK(est > 0.0);                               // SIGN: positive planted -> positive
    CHECK(near_abs(est, planted, 1.5 * dt));
}

TEST_CASE("timesync sim e2e: straight() (no rotation) is gated -> no confident estimate") {
    Trajectory tr = Trajectory::straight();         // ‖ω‖ == 0 everywhere

    SourceParams pr; pr.id = 0;
    SyntheticSource ref(tr, pr);
    SourceParams ps; ps.id = 1; ps.time_offset_s = 0.03;
    SyntheticSource src(tr, ps);

    const Scalar tick_hz = 100.0;
    const Scalar dt = 1.0 / tick_hz;
    const Scalar h  = dt;

    auto tsp = make_ts(); TimeSync& ts = *tsp;
    Config c = make_cfg(MatchMetric::L2, tick_hz, 0.1);
    c.excitation_min_var = 1e-4;                    // a meaningful gate
    REQUIRE(ts.configure(c, 0) == Status::Ok);

    for (double t = tr.t0_s() + 2 * h; t <= tr.end_s() - 0.02; t += dt) {
        ts.push(0, secs(t), omega_norm_at(ref, t, h));
        ts.push(1, secs(t), omega_norm_at(src, t, h));
    }
    ts.update();
    CHECK(ts.confidence(1) == doctest::Approx(0.0));   // gated: no vote
}

TEST_CASE("timesync sim e2e: deterministic (same input -> identical estimate)") {
    Trajectory tr = Trajectory::omega_varying();
    SourceParams pr; pr.id = 0;
    SourceParams ps; ps.id = 1; ps.time_offset_s = 0.025;
    SyntheticSource ref(tr, pr);
    SyntheticSource src(tr, ps);

    const Scalar tick_hz = 100.0, dt = 1.0 / tick_hz, h = dt;

    auto run_once = [&]() {
        auto tsp = make_ts(); TimeSync& ts = *tsp;
        ts.configure(make_cfg(MatchMetric::L2, tick_hz, 0.1), 0);
        for (double t = tr.t0_s() + 2 * h; t <= tr.end_s() - 0.05; t += dt) {
            ts.push(0, secs(t), omega_norm_at(ref, t, h));
            ts.push(1, secs(t), omega_norm_at(src, t, h));
        }
        ts.update();
        return ts.offset(1);
    };
    const Scalar a = run_once();
    const Scalar b = run_once();
    CHECK(a == b);     // bit-identical replay
}

// ===========================================================================
// Estimator wiring (the integration the slice ships): the Estimator feeds ‖ω‖
// into TimeSync each step, recovers a source's RELATIVE clock offset, applies it to
// that source's fusion query, and surfaces it in CalibSnapshot. With time-sync OFF
// the offsets stay at the configured priors (behaves exactly as Slice 2).
// ===========================================================================
namespace {
// Drive the estimator over `tr` at `tick_hz` from `from_s`..`to_s`, returning the
// LAST successful Result. `enable_ts` toggles time-sync; sensors must outlive the call.
bool drive(Estimator& est, const Trajectory& tr,
           const std::vector<const ISource*>& srcs, Scalar from_s, Scalar to_s,
           Scalar tick_hz, Result& last_out) {
    const Timestamp tick = secs(1.0 / tick_hz);
    bool any = false;
    for (Timestamp now = secs(from_s); now <= secs(to_s); now += tick) {
        if (ok(est.step(now))) { last_out = est.latest(); any = true; }
    }
    (void)srcs;
    return any;
}
} // namespace

TEST_CASE("estimator e2e: recovers + applies a RELATIVE planted offset, surfaces it") {
    Trajectory tr = Trajectory::omega_varying();

    // Reference id 0 (offset 0) + source id 1 with a planted +0.04 s RELATIVE offset.
    const Scalar planted = 0.04;
    SourceParams pr; pr.id = 0; pr.time_offset_s = 0.0;
    SourceParams ps; ps.id = 1; ps.time_offset_s = planted;
    SyntheticSource s0(tr, pr);
    SyntheticSource s1(tr, ps);

    // Sensor configs: priors at 0 offset (so any recovered offset is from time-sync).
    SensorConfig sc[2];
    sc[0].id = 0; sc[0].is_reference = true;
    sc[1].id = 1;

    Config cfg = make_cfg(MatchMetric::L2, 100.0, 0.1);
    cfg.max_sources         = 2;
    cfg.reference_sensor_id = 0;
    cfg.sensors             = sc;
    cfg.sensor_count        = 2;
    cfg.commit_concentration = 0.5;     // a few sharp votes clear this

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    REQUIRE(est.add_source(&s0) == Status::Ok);
    REQUIRE(est.add_source(&s1) == Status::Ok);

    std::vector<const ISource*> srcs = {&s0, &s1};
    Result last;
    REQUIRE(drive(est, tr, srcs, 0.3, tr.end_s() - 0.1, 100.0, last));

    // Find source 1's calibration snapshot.
    int idx = -1;
    for (int i = 0; i < last.source_count; ++i) if (last.calib[i].id == 1) idx = i;
    REQUIRE(idx >= 0);
    INFO("recovered offset=", last.calib[idx].time_offset_s,
         " conf=", last.calib[idx].confidence);
    CHECK(last.calib[idx].confidence > 0.0);
    CHECK(last.calib[idx].time_offset_s > 0.0);                 // correct SIGN
    CHECK(near_abs(last.calib[idx].time_offset_s, planted, 0.02));
    CHECK(last.calib[idx].committed);                           // cleared the commit gate
}

TEST_CASE("estimator e2e: timesync OFF leaves offsets at the prior (Slice-2 behavior)") {
    Trajectory tr = Trajectory::omega_varying();
    SourceParams pr; pr.id = 0;
    SourceParams ps; ps.id = 1; ps.time_offset_s = 0.04;
    SyntheticSource s0(tr, pr);
    SyntheticSource s1(tr, ps);

    SensorConfig sc[2];
    sc[0].id = 0; sc[0].is_reference = true;
    sc[1].id = 1; sc[1].prior_time_offset_s = 0.01;   // a non-zero PRIOR

    Config cfg = make_cfg(MatchMetric::L2, 100.0, 0.1);
    cfg.timesync_enabled    = false;                  // OFF
    cfg.max_sources         = 2;
    cfg.reference_sensor_id = 0;
    cfg.sensors             = sc;
    cfg.sensor_count        = 2;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    REQUIRE(est.add_source(&s0) == Status::Ok);
    REQUIRE(est.add_source(&s1) == Status::Ok);

    std::vector<const ISource*> srcs = {&s0, &s1};
    Result last;
    REQUIRE(drive(est, tr, srcs, 0.3, tr.end_s() - 0.1, 100.0, last));

    int idx = -1;
    for (int i = 0; i < last.source_count; ++i) if (last.calib[i].id == 1) idx = i;
    REQUIRE(idx >= 0);
    // Off => the applied offset is exactly the configured prior, confidence 0, uncommitted.
    CHECK(last.calib[idx].time_offset_s == doctest::Approx(0.01));
    CHECK(last.calib[idx].confidence == doctest::Approx(0.0));
    CHECK_FALSE(last.calib[idx].committed);
}
