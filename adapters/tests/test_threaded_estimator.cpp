// Adapters Slice 13: ThreadedEstimator determinism. The wrapper owns a worker thread that
// pumps step() over an injected timestamp queue. Driving a FIXED set of stamps through it and
// joining (drain_and_stop) must leave the final snapshot EQUAL to a single-threaded reference
// run over the same stamps — proving the wrapper does not perturb the core's determinism and
// publishes a tear-free Result. No wall-clock timing is asserted.
#include <doctest/doctest.h>

#include "ofc_adapters/threaded_estimator.hpp"

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/status.hpp"

#include "adapter_test_util.hpp"

#include <vector>

using namespace ofc;
using namespace adptest;

namespace {

// The fixed stamp schedule both runs use (50 Hz over ~1.5 s).
std::vector<Timestamp> stamp_schedule() {
    std::vector<Timestamp> ts;
    const Timestamp tick = secs_to_ns(0.02);
    for (Timestamp now = secs_to_ns(0.2); now <= secs_to_ns(1.5); now += tick) ts.push_back(now);
    return ts;
}

// Single-threaded reference: step() every stamp, return the final latest().
Result reference_run(const Config& cfg, const std::vector<Timestamp>& ts) {
    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);
    for (Timestamp now : ts) (void)est.step(now);
    return est.latest();
}

// Compare the load-bearing fields of two Results (the fusion frontier pose + phase + per-source
// calibration). Bit determinism is expected (same inputs, no noise, no wall-clock), so an exact
// compare is appropriate; we use a tiny tolerance only to be robust to any benign reordering.
void check_same(const Result& a, const Result& b) {
    CHECK(a.phase == b.phase);
    CHECK(a.source_count == b.source_count);
    CHECK(a.frontier.stamp == b.frontier.stamp);
    CHECK((a.frontier.pose.t - b.frontier.pose.t).norm() == doctest::Approx(0.0));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            CHECK(a.frontier.pose.R(i, j) == doctest::Approx(b.frontier.pose.R(i, j)));
    for (int i = 0; i < a.source_count; ++i) {
        CHECK(a.calib[i].id == b.calib[i].id);
        CHECK(a.calib[i].scale == doctest::Approx(b.calib[i].scale));
        CHECK(a.calib[i].committed == b.calib[i].committed);
        CHECK(a.calib[i].extrinsic_committed == b.calib[i].extrinsic_committed);
    }
}

} // namespace

TEST_CASE("ThreadedEstimator determinism: batch-submit + drain matches the single-threaded run") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);
    const std::vector<Timestamp> ts = stamp_schedule();

    const Result ref = reference_run(cfg, ts);

    // Threaded run: start, submit the whole batch, drain_and_stop (joins the worker after every
    // stamp is step()'d). The final snapshot must equal the reference.
    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);

    adapters::ThreadedEstimator te(est);
    CHECK(te.start() == Status::Ok);
    CHECK(te.running());
    CHECK(te.submit_batch(ts.data(), static_cast<int>(ts.size())) == Status::Ok);
    CHECK(te.drain_and_stop() == Status::Ok);
    CHECK_FALSE(te.running());

    CHECK(te.steps_done() > 0);                       // it actually pumped
    check_same(te.snapshot(), ref);
}

TEST_CASE("ThreadedEstimator determinism: incremental submit + drain matches the reference") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);
    const std::vector<Timestamp> ts = stamp_schedule();
    const Result ref = reference_run(cfg, ts);

    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);

    adapters::ThreadedEstimator te(est);
    CHECK(te.start() == Status::Ok);
    for (Timestamp now : ts) CHECK(te.submit(now) == Status::Ok);  // one at a time
    CHECK(te.drain_and_stop() == Status::Ok);
    check_same(te.snapshot(), ref);
}

TEST_CASE("ThreadedEstimator: submit-before-start then drain (inline drain path) matches") {
    // A wrapper that never start()s still honors drain_and_stop by draining inline on the
    // caller thread (still single-threaded into the core). Same final state as the reference.
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);
    const std::vector<Timestamp> ts = stamp_schedule();
    const Result ref = reference_run(cfg, ts);

    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);

    adapters::ThreadedEstimator te(est);
    CHECK(te.submit_batch(ts.data(), static_cast<int>(ts.size())) == Status::Ok);  // before start
    CHECK(te.drain_and_stop() == Status::Ok);        // inline drain (no worker)
    check_same(te.snapshot(), ref);
}

TEST_CASE("ThreadedEstimator: snapshot before any step is the default Result; stop is idempotent") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);
    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);

    adapters::ThreadedEstimator te(est);
    const Result before = te.snapshot();
    CHECK(before.phase == Phase::Init);              // default-constructed, tear-free
    CHECK(before.source_count == 0);

    CHECK(te.start() == Status::Ok);
    CHECK(te.start() == Status::Ok);                 // idempotent
    CHECK(te.stop() == Status::Ok);
    CHECK(te.stop() == Status::Ok);                  // idempotent
    // After stop, further submits are rejected.
    CHECK(te.submit(secs_to_ns(2.0)) == Status::Rejected);
}

TEST_CASE("ThreadedEstimator: drain_and_stop() after stop() honors the drop (no resurrection)") {
    // Contract: stop() DROPS pending stamps; a later drain_and_stop() must NOT replay them via
    // the inline-drain branch. Submit before start, stop() (drops), then drain_and_stop(): the
    // estimator must have stepped ZERO stamps.
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);
    const std::vector<Timestamp> ts = stamp_schedule();

    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);

    adapters::ThreadedEstimator te(est);
    CHECK(te.submit_batch(ts.data(), static_cast<int>(ts.size())) == Status::Ok);  // before start
    CHECK(te.stop() == Status::Ok);                  // drops the pending batch
    CHECK(te.drain_and_stop() == Status::Ok);        // must NOT resurrect the dropped stamps
    CHECK(te.steps_done() == 0);                      // nothing was processed
    CHECK_FALSE(te.running());
}
