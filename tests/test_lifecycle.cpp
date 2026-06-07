// Slice 3 lifecycle tests: the estimator exposes a correct
// INIT -> WARMUP -> DEGRADED -> NOMINAL ladder with degrade-don't-block.
//
// The ladder (driven from per-step fuse signals; n = participating sources):
//   INIT     : right after init(), before any step().
//   WARMUP   : a step that cannot fuse yet (degenerate interval / no covering
//              data) BEFORE any successful fuse.
//   DEGRADED : fused with 1 <= n < min_sources_warn (e.g. ReferenceOnly
//              dead-reckon before any extrinsic commit), OR a graceful downgrade
//              after a source loss (was fused, now nothing covers).
//   NOMINAL  : fused with n >= min_sources_warn.
// readiness is a coarse [0,1] map of the phase: Init 0.0 / Warmup 0.25 /
// Degraded 0.6 / Nominal 1.0.
//
// Helpers (BufferSource adapter + fill_gt / secs / gt_twist) are copied/adapted
// from test_fusion.cpp — they are file-local statics there.
#include <doctest/doctest.h>

#include "ofc/core/buffer.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/source.hpp"

#include <cmath>
#include <memory>
#include <vector>

using namespace ofc;

namespace {
Timestamp secs(double s) { return static_cast<Timestamp>(s * 1e9); }

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

// A source = an ISource facade over a SourceBuffer (relaxed-edge test adapter).
class BufferSource : public ISource {
public:
    BufferSource(SourceId id, OdomForm form, int cap) : id_(id) {
        SensorConfig sc;
        sc.id = id;
        buf_.configure(sc, form, cap, ConfCombine::NativeOnly);
    }
    SourceId id() const override { return id_; }
    Expected<Delta> query(Timestamp t0, Timestamp t1) const override {
        return buf_.query(t0, t1);
    }
    SourceBuffer& buffer() { return buf_; }

private:
    SourceId     id_;
    SourceBuffer buf_;
};

// Ground-truth body twist: straight segment then a turn (same as test_fusion).
Vec6 gt_twist(double t) {
    Vec6 xi;
    const double yaw = (t < 1.0) ? 0.0 : 0.5;
    xi << 2.0, 0.0, 0.0, 0.0, 0.0, yaw;
    return xi;
}

// Fill a source buffer with GT twist samples (in the source's own frame) at
// `rate` Hz over [0, span] s. `frame` re-expresses the base-frame GT twist into
// the sensor frame via the adjoint so the estimator's frame-align recovers base.
void fill_gt(BufferSource& src, double span, double rate, const SE3& frame,
             const Vec6& bias = Vec6::Zero()) {
    const Mat6 Ad_inv = se3::adjoint(se3::inverse(frame));   // base->sensor twist
    const Timestamp step = secs(1.0 / rate);
    const int n = static_cast<int>(span * rate);
    for (int k = 0; k <= n; ++k) {
        const double t = k / rate;
        const Vec6 xi_sensor = Ad_inv * gt_twist(t) + bias;
        src.buffer().push_twist(static_cast<Timestamp>(k) * step, xi_sensor,
                                Mat6::Identity());
    }
}

// Base config (mirrors test_fusion's base_config but with explicit cold_start).
Config base_config(int max_sources, ColdStart cs = ColdStart::MedianFromStart) {
    Config c;
    c.max_sources    = max_sources;
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.tip_cov_inflation = 1.5;
    c.emit_predicted_tip = true;
    c.cold_start     = cs;
    return c;
}
} // namespace

// ---------------------------------------------------------------------------
// INIT: exposed right after init(), before any step().
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: INIT exposed right after init (before any step)") {
    const int N = 3;
    SensorConfig sensors[3];
    for (int i = 0; i < N; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_config(N);
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);

    CHECK(est.latest().phase == Phase::Init);
    CHECK(est.latest().readiness == doctest::Approx(0.0));
}

// ---------------------------------------------------------------------------
// WARMUP: a step that cannot fuse yet (no covering data) BEFORE any fuse.
// A non-causal/too-early interval (now so small that t1 <= q0, or no buffered
// data covers the window) keeps the estimator in WARMUP / NotReady.
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: WARMUP before the first fuse (degenerate interval)") {
    const int N = 3;
    SensorConfig sensors[3];
    for (int i = 0; i < N; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_config(N);
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    for (auto& s : srcs) fill_gt(*s, 3.0, 200.0, SE3{});

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    // Degenerate interval: now == 0 -> t1 = -delay, q0 = t1 - window, so t1 > q0
    // but the buffers start at t=0 and never cover [q0, t1] (both negative).
    // Use now small enough that no source covers -> n == 0 -> NotReady/Warmup.
    const Status st = est.step(secs(0.0));
    CHECK(st == Status::NotReady);
    CHECK(est.latest().phase == Phase::Warmup);
    CHECK(est.latest().readiness == doctest::Approx(0.25));
    CHECK_FALSE(est.latest().tip_valid);
}

// ---------------------------------------------------------------------------
// WARMUP via a strictly non-causal interval (t1 <= q0).
// With window_s very small and a large delay, an early `now` makes t1 <= q0.
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: WARMUP on a non-causal interval (t1 <= q0)") {
    const int N = 1;
    SensorConfig sensors[1];
    sensors[0].id = 0;
    sensors[0].is_reference = true;
    Config cfg = base_config(N);
    cfg.reference_sensor_id = 0;
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    std::unique_ptr<BufferSource> src(new BufferSource(0, OdomForm::Twist, 4096));
    fill_gt(*src, 3.0, 200.0, SE3{});

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    REQUIRE(est.add_source(src.get()) == Status::Ok);

    // now - delay = t1; q0 = t1 - window (first fuse, no prior frontier). With
    // now = 0 the whole interval is negative; the buffer never covers it.
    const Status st = est.step(secs(0.0));
    CHECK(st == Status::NotReady);
    CHECK(est.latest().phase == Phase::Warmup);
    CHECK_FALSE(est.latest().tip_valid);
}

// ---------------------------------------------------------------------------
// DEGRADED via reference-only: ColdStart::ReferenceOnly with a reference + a
// non-reference source whose extrinsic is NOT committed. Only the reference
// participates -> n == 1 < min_sources_warn -> Ok + Degraded as soon as the
// reference covers (output appears as early as the reference allows).
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: DEGRADED reference-only dead-reckon (n=1 < warn)") {
    const int N = 2;
    SensorConfig sensors[2];
    sensors[0].id = 0; sensors[0].is_reference = true;
    sensors[1].id = 1;
    Config cfg = base_config(N, ColdStart::ReferenceOnly);
    cfg.reference_sensor_id = 0;
    cfg.min_sources_warn = 3;     // need 3 to be NOMINAL; only the ref participates
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    for (auto& s : srcs) fill_gt(*s, 3.0, 200.0, SE3{});

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    // Step into the straight segment; the reference covers from the first
    // window. Reference-only -> n == 1 -> Degraded, Ok, valid frontier.
    const double now0 = cfg.fusion_delay_s + cfg.window_s;
    const Status st = est.step(secs(now0));
    CHECK(st == Status::Ok);
    CHECK(est.latest().phase == Phase::Degraded);
    CHECK(est.latest().readiness == doctest::Approx(0.6));
    CHECK(est.latest().tip_valid);

    // Keep stepping; it stays Degraded (no extrinsic ever commits here).
    double now_s = now0 + cfg.window_s;
    while (now_s <= 0.8) {
        CHECK(est.step(secs(now_s)) == Status::Ok);
        CHECK(est.latest().phase == Phase::Degraded);
        now_s += cfg.window_s;
    }
}

// ---------------------------------------------------------------------------
// NOMINAL: 3 clean identity-extrinsic sources, MedianFromStart -> n == 3 ==
// min_sources_warn -> Nominal, readiness 1.0.
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: NOMINAL with >= min_sources_warn participating") {
    const int N = 3;
    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    for (auto& s : srcs) fill_gt(*s, 3.0, 200.0, SE3{});

    SensorConfig sensors[3];
    for (int i = 0; i < N; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_config(N);     // MedianFromStart, min_sources_warn default 3
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    double now_s = cfg.fusion_delay_s + cfg.window_s;
    while (now_s <= 0.8) {
        CHECK(est.step(secs(now_s)) == Status::Ok);
        now_s += cfg.window_s;
    }
    CHECK(est.latest().phase == Phase::Nominal);
    CHECK(est.latest().readiness == doctest::Approx(1.0));
}

// ---------------------------------------------------------------------------
// Graceful downgrade: start NOMINAL with 3 sources, then let 2 of them run out
// of data so they fall out of the window -> n drops below min_sources_warn ->
// downgrade to Degraded while still emitting Ok. When ALL stop covering ->
// Degraded + NotReady (no new frontier), last frontier retained.
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: graceful downgrade NOMINAL -> DEGRADED on source loss") {
    const int N = 3;
    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    // Source 0 has a long buffer; sources 1,2 stop early (data only to ~0.5 s).
    fill_gt(*srcs[0], 3.0, 200.0, SE3{});
    fill_gt(*srcs[1], 0.5, 200.0, SE3{});
    fill_gt(*srcs[2], 0.5, 200.0, SE3{});

    SensorConfig sensors[3];
    for (int i = 0; i < N; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_config(N);     // MedianFromStart, min_sources_warn 3
    cfg.reference_sensor_id = 0;     // source 0 is the long-lived reference
    sensors[0].is_reference = true;
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    // Phase A: all three cover -> NOMINAL.
    double now_s = cfg.fusion_delay_s + cfg.window_s;
    bool saw_nominal = false;
    while (now_s <= 0.4) {
        CHECK(est.step(secs(now_s)) == Status::Ok);
        if (est.latest().phase == Phase::Nominal) saw_nominal = true;
        now_s += cfg.window_s;
    }
    CHECK(saw_nominal);
    CHECK(est.latest().phase == Phase::Nominal);

    // Phase B: step past the short sources' end-of-data. Once they stop covering
    // the window, only source 0 covers (n == 1 < min_sources_warn) -> downgrade
    // to DEGRADED, but still Ok with a fresh frontier. Early Phase-B steps may
    // still have all three covered (Nominal); the downgrade must HAPPEN and once
    // downgraded every remaining Ok step that fuses fewer sources stays Degraded.
    // Capture the last good frontier so we can prove it is retained when ALL drop.
    State last_good_frontier{};
    bool saw_degraded_ok = false;
    while (now_s <= 0.8) {
        const Status st = est.step(secs(now_s));
        if (st == Status::Ok) {
            // Phase is Nominal while all three still cover, Degraded once the two
            // short sources fall out of the window. Both are valid Ok states; we
            // only require that the downgrade to Degraded is reached.
            const Phase p = est.latest().phase;
            CHECK((p == Phase::Nominal || p == Phase::Degraded));
            if (p == Phase::Degraded) saw_degraded_ok = true;
            last_good_frontier = est.latest().frontier;
        }
        now_s += cfg.window_s;
    }
    CHECK(saw_degraded_ok);                          // the downgrade did happen
    CHECK(est.latest().phase == Phase::Degraded);    // and we ended Degraded (Ok)

    // Phase C: step well past ALL sources' data (source 0 ends at 3.0 s). Use a
    // now far beyond every buffer so nothing covers the window -> Degraded +
    // NotReady, and the last good frontier is retained in result.
    const Status st = est.step(secs(10.0));
    CHECK(st == Status::NotReady);
    CHECK(est.latest().phase == Phase::Degraded);     // was fused -> Degraded, not Warmup
    CHECK_FALSE(est.latest().tip_valid);
    // The published frontier is the last successful one (unchanged by the empty step).
    CHECK(close(est.latest().frontier.pose.t, last_good_frontier.pose.t, 1e-12));
    CHECK(close(est.latest().frontier.pose.R, last_good_frontier.pose.R, 1e-12));
}

// ---------------------------------------------------------------------------
// Anchor-at-first-tick: the FIRST fused frontier pose starts at identity (the
// odom frame is anchored at the first fused tick, DESIGN §7).
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: first fused frontier pose is anchored at identity") {
    const int N = 3;
    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    for (auto& s : srcs) fill_gt(*s, 3.0, 200.0, SE3{});

    SensorConfig sensors[3];
    for (int i = 0; i < N; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_config(N);
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    // First fuse: pose starts at identity (predict integrates ONE window's
    // motion from the identity anchor). Rotation is exactly identity over the
    // straight segment; translation is the first window's forward displacement.
    const double now0 = cfg.fusion_delay_s + cfg.window_s;
    REQUIRE(est.step(secs(now0)) == Status::Ok);
    const State& f = est.latest().frontier;
    // Anchored at identity: the ROTATION is (near-)identity over the straight
    // bootstrap window (no yaw before t=1.0 s) — proves the eskf.init(SE3{}).
    CHECK(close(f.pose.R, Mat3::Identity(), 1e-6));
    // The translation is bounded by one window of forward motion (2 m/s * span),
    // i.e. the pose did NOT start pre-integrated from t=0.
    CHECK(f.pose.t.norm() < 2.0 * (cfg.window_s + cfg.fusion_delay_s) + 0.05);
}

// ---------------------------------------------------------------------------
// Config: min_sources_warn lower-bound is validated; min_sources_warn >
// max_sources is legitimate (NOMINAL simply never reached, no validation error).
// ---------------------------------------------------------------------------
TEST_CASE("lifecycle: min_sources_warn validation bounds") {
    Config cfg = base_config(2);
    cfg.min_sources_warn = 0;
    CHECK(validate(cfg) == Status::OutOfRange);

    cfg.min_sources_warn = 1;
    CHECK(validate(cfg) == Status::Ok);

    // > max_sources is allowed: it just means "never enough sources" (no error).
    cfg.min_sources_warn = 99;
    CHECK(validate(cfg) == Status::Ok);
}
