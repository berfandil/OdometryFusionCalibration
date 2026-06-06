// Slice 2 end-to-end tracer-bullet tests: N synthetic ISource's (backed by
// SourceBuffer) fed a known trajectory -> Estimator::step -> fused frontier tracks
// the integrated ground truth; an outlier source is rejected by the median; replay
// is deterministic; the predicted tip leads the frontier; frame-alignment via a
// non-identity prior extrinsic is correct.
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
constexpr Scalar kPi = 3.14159265358979323846;
constexpr Timestamp kSec = 1000000000LL;
Timestamp secs(double s) { return static_cast<Timestamp>(s * 1e9); }

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

// A source = an ISource façade over a SourceBuffer (relaxed-edge test adapter).
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

// Ground-truth body twist as a function of time: straight segment, then a turn.
//   t < 1.0 s : 2 m/s forward, no rotation (straight)
//   t >= 1.0 s: 2 m/s forward + 0.5 rad/s yaw (turn)
Vec6 gt_twist(double t) {
    Vec6 xi;
    const double yaw = (t < 1.0) ? 0.0 : 0.5;
    xi << 2.0, 0.0, 0.0, 0.0, 0.0, yaw;
    return xi;
}

// Integrate the GT twist from 0 to T by fine const-twist composition (the same
// model the buffers integrate with -> a fair reference for the fused pose).
SE3 gt_pose_at(double T, double dt = 0.001) {
    SE3 p;
    for (double t = 0.0; t + dt / 2 < T; t += dt) {
        p = se3::compose(p, se3::exp(gt_twist(t) * dt));
    }
    return p;
}

// Fill a source buffer with GT twist samples (in the source's own frame) at `rate`
// Hz over [0, span] seconds. `frame` re-expresses the base-frame GT twist into the
// sensor frame via the adjoint (so that after the estimator frame-aligns with the
// same extrinsic, the recovered motion is the base GT). `bias` adds a constant
// body twist (used to plant an outlier source).
void fill_gt(BufferSource& src, double span, double rate, const SE3& frame,
             const Vec6& bias = Vec6::Zero()) {
    const Mat6 Ad_inv = se3::adjoint(se3::inverse(frame));   // base->sensor twist map
    const Timestamp step = secs(1.0 / rate);
    const int n = static_cast<int>(span * rate);
    for (int k = 0; k <= n; ++k) {
        const double t = k / rate;
        const Vec6 xi_sensor = Ad_inv * gt_twist(t) + bias;
        src.buffer().push_twist(static_cast<Timestamp>(k) * step, xi_sensor,
                                Mat6::Identity());
    }
}

Config base_config(int max_sources) {
    Config c;
    c.max_sources    = max_sources;
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.tip_cov_inflation = 1.5;
    c.emit_predicted_tip = true;
    return c;
}
} // namespace

// ---------------------------------------------------------------------------
// Fused frontier tracks the integrated ground truth (3 clean sources).
// ---------------------------------------------------------------------------
TEST_CASE("fused frontier tracks integrated GT over a straight+turn trajectory") {
    const int N = 3;
    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    // All identity extrinsics -> source frame == base frame.
    for (auto& s : srcs) fill_gt(*s, 3.0, 200.0, SE3{});

    SensorConfig sensors[3];
    for (int i = 0; i < N; ++i) { sensors[i].id = static_cast<SourceId>(i); }
    Config cfg = base_config(N);
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    // Pump steps. The frontier at `now` is t1 = now - delay; after stepping up to
    // now = T, the integrated fused pose should equal GT integrated to t1.
    const Scalar dt = cfg.window_s;             // one window per step (no overlap)
    double now_s = cfg.fusion_delay_s + cfg.window_s;   // first step that covers a window
    double last_t1 = 0.0;
    while (now_s <= 2.5) {
        const Status st = est.step(secs(now_s));
        CHECK(st == Status::Ok);
        last_t1 = now_s - cfg.fusion_delay_s;
        now_s += dt;
    }

    const SE3 gt = gt_pose_at(last_t1);
    const State& f = est.latest().frontier;
    // Track within a small tolerance (const-twist integration vs the buffer's
    // per-sample composition + endpoint interpolation).
    CHECK(close(f.pose.t, gt.t, 1e-2));
    CHECK(close(f.pose.R, gt.R, 1e-2));
    CHECK(est.latest().phase == Phase::Nominal);
}

// ---------------------------------------------------------------------------
// An outlier source is rejected by the median (3 inliers + 1 gross outlier).
// ---------------------------------------------------------------------------
TEST_CASE("median rejects a biased outlier source and still tracks GT") {
    const int N = 4;
    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    // Three clean sources...
    for (int i = 0; i < 3; ++i) fill_gt(*srcs[i], 3.0, 200.0, SE3{});
    // ...one badly biased source (large spurious forward + yaw rate).
    Vec6 bias; bias << 5.0, 3.0, 0.0, 0.0, 0.0, 2.0;
    fill_gt(*srcs[3], 3.0, 200.0, SE3{}, bias);

    SensorConfig sensors[4];
    for (int i = 0; i < N; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_config(N);
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    const Scalar dt = cfg.window_s;
    double now_s = cfg.fusion_delay_s + cfg.window_s;
    double last_t1 = 0.0;
    while (now_s <= 2.5) {
        CHECK(est.step(secs(now_s)) == Status::Ok);
        last_t1 = now_s - cfg.fusion_delay_s;
        now_s += dt;
    }

    const SE3 gt = gt_pose_at(last_t1);
    const State& f = est.latest().frontier;
    // The outlier (5 m/s side bias) would push the fused x/y far off GT if it were
    // not rejected. With the median + 3 inliers it tracks GT closely.
    CHECK(close(f.pose.t, gt.t, 5e-2));
    CHECK(close(f.pose.R, gt.R, 5e-2));

    // Prove REJECTION (not lucky averaging): the biased source (id 3) must sit far
    // from the consensus — a large residual — while the three inliers sit close. The
    // median's 1/d reweight is what collapses the outlier's pull.
    const Result& r = est.latest();
    REQUIRE(r.source_count == N);
    Scalar outlier_resid = -1.0;
    Scalar max_inlier_resid = 0.0;
    for (int i = 0; i < r.source_count; ++i) {
        const SourceHealth& h = r.health[i];
        REQUIRE(h.in_window);
        if (h.id == 3) outlier_resid = h.residual;
        else           max_inlier_resid = std::max(max_inlier_resid, h.residual);
    }
    // The biased source's residual is large in absolute terms...
    CHECK(outlier_resid > 0.1);
    // ...and dwarfs every inlier's residual (rejection, not averaging).
    CHECK(outlier_resid > 10.0 * max_inlier_resid);
}

// ---------------------------------------------------------------------------
// Determinism: identical inputs -> bit-identical output.
// ---------------------------------------------------------------------------
TEST_CASE("replay is deterministic (bit-identical fused output)") {
    auto run = []() {
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
        est.init(cfg);
        for (auto& s : srcs) est.add_source(s.get());

        double now_s = cfg.fusion_delay_s + cfg.window_s;
        while (now_s <= 2.0) { est.step(secs(now_s)); now_s += cfg.window_s; }
        return est.latest().frontier;
    };

    const State a = run();
    const State b = run();
    CHECK((a.pose.R.array() == b.pose.R.array()).all());   // exactly equal
    CHECK((a.pose.t.array() == b.pose.t.array()).all());
    CHECK((a.twist.xi.array() == b.twist.xi.array()).all());
    CHECK((a.cov.array() == b.cov.array()).all());
}

// ---------------------------------------------------------------------------
// The predicted tip leads the frontier (extrapolated forward in time/space).
// ---------------------------------------------------------------------------
TEST_CASE("predicted tip leads the frontier") {
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

    // Step into the straight segment (moving forward in +x).
    double now_s = cfg.fusion_delay_s + cfg.window_s;
    while (now_s <= 0.6) { CHECK(est.step(secs(now_s)) == Status::Ok); now_s += cfg.window_s; }

    const Result& r = est.latest();
    REQUIRE(r.tip_valid);
    // Forward motion -> tip x ahead of frontier x.
    CHECK(r.tip.pose.t.x() > r.frontier.pose.t.x());
    // Tip stamp is `now`; frontier stamp is the lagged frontier t1 < now.
    CHECK(r.tip.stamp > r.frontier.stamp);
    // Tip covariance is inflated relative to the frontier.
    CHECK(r.tip.cov.trace() > r.frontier.cov.trace());
}

// ---------------------------------------------------------------------------
// Frame-alignment: a source mounted with a non-identity extrinsic is aligned back
// to the base frame (so it agrees with an identity-mounted source on the same GT).
// ---------------------------------------------------------------------------
TEST_CASE("non-identity prior extrinsic is correctly frame-aligned") {
    const int N = 3;
    // Source 0 aligned (identity); sources 1,2 mounted with a yaw+offset extrinsic.
    SE3 X;
    X.R = so3::exp(Vec3(0, 0, kPi / 3));   // 60 deg yaw mount
    X.t = Vec3(0.3, -0.2, 0.1);

    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < N; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i),
                                           OdomForm::Twist, 4096));
    fill_gt(*srcs[0], 3.0, 200.0, SE3{});   // identity mount
    fill_gt(*srcs[1], 3.0, 200.0, X);       // tilted mount (GT re-expressed in sensor frame)
    fill_gt(*srcs[2], 3.0, 200.0, X);

    SensorConfig sensors[3];
    sensors[0].id = 0;  sensors[0].prior_extrinsic = SE3{};
    sensors[1].id = 1;  sensors[1].prior_extrinsic = X;
    sensors[2].id = 2;  sensors[2].prior_extrinsic = X;

    Config cfg = base_config(N);
    cfg.sensors = sensors;
    cfg.sensor_count = N;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    const Scalar dt = cfg.window_s;
    double now_s = cfg.fusion_delay_s + cfg.window_s;
    double last_t1 = 0.0;
    while (now_s <= 2.5) {
        CHECK(est.step(secs(now_s)) == Status::Ok);
        last_t1 = now_s - cfg.fusion_delay_s;
        now_s += dt;
    }

    // After alignment all three sources represent the SAME base motion -> the fused
    // frontier tracks the base-frame GT (not the tilted sensor-frame motion).
    const SE3 gt = gt_pose_at(last_t1);
    const State& f = est.latest().frontier;
    CHECK(close(f.pose.t, gt.t, 1e-2));
    CHECK(close(f.pose.R, gt.R, 1e-2));
}
