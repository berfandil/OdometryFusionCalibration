// Slice 1 unit tests: SourceBuffer normalization, delta query, interpolation,
// covariance accumulation + native-(+)-modeled combine, and ring-buffer edges.
#include <doctest/doctest.h>

#include "ofc/core/buffer.hpp"
#include "ofc/core/lie.hpp"

#include <cmath>

using namespace ofc;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

// 1 second in the Timestamp (nanosecond) unit.
constexpr Timestamp kSec = 1000000000LL;
Timestamp secs(double s) { return static_cast<Timestamp>(s * 1e9); }

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

SensorConfig make_sensor(Scalar per_m = 0.0, Scalar per_rad = 0.0,
                         bool native = true) {
    SensorConfig s;
    s.id = 0;
    s.modeled_noise_per_m   = per_m;
    s.modeled_noise_per_rad = per_rad;
    s.native_confidence     = native;
    return s;
}
} // namespace

// ---------------------------------------------------------------------------
// configure / basic accessors
// ---------------------------------------------------------------------------
TEST_CASE("configure rejects tiny capacity and reports geometry") {
    SourceBuffer b;
    CHECK(b.configure(make_sensor(), OdomForm::Increment, 1, ConfCombine::Sum)
          == Status::OutOfRange);

    CHECK(b.configure(make_sensor(), OdomForm::Increment, 4, ConfCombine::Sum)
          == Status::Ok);
    CHECK(b.capacity() == 4);
    CHECK(b.size() == 0);
    CHECK(b.configured());
    CHECK(b.oldest() == 0);
    CHECK(b.newest() == 0);
}

TEST_CASE("push wrong form for configuration is rejected") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Twist, 8, ConfCombine::Sum);
    SE3 incr;
    CHECK(b.push_increment(kSec, incr, Mat6::Identity()) == Status::InvalidConfig);
    CHECK(b.push_twist(kSec, Vec6::Zero(), Mat6::Zero()) == Status::Ok);
}

TEST_CASE("non-increasing stamps are rejected") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 8, ConfCombine::Sum);
    SE3 incr;
    CHECK(b.push_increment(secs(1.0), incr, Mat6::Zero()) == Status::Ok);
    CHECK(b.push_increment(secs(1.0), incr, Mat6::Zero()) == Status::OutOfRange);
    CHECK(b.push_increment(secs(0.5), incr, Mat6::Zero()) == Status::OutOfRange);
    CHECK(b.push_increment(secs(2.0), incr, Mat6::Zero()) == Status::Ok);
}

// ---------------------------------------------------------------------------
// Twist integration: constant velocity gives the analytic delta on a sub-window
// ---------------------------------------------------------------------------
TEST_CASE("constant-velocity twist integrates to analytic delta over a sub-window") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Twist, 256, ConfCombine::NativeOnly);

    // Constant body twist: 2 m/s forward (x), yaw rate 0.5 rad/s.
    Vec6 xi;
    xi << 2.0, 0.0, 0.0, 0.0, 0.0, 0.5;

    // Sample at 100 Hz from t=0 to t=2s (201 samples; capacity holds them all).
    for (int k = 0; k <= 200; ++k) {
        const Timestamp t = static_cast<Timestamp>(k) * (kSec / 100);
        CHECK(b.push_twist(t, xi, Mat6::Identity()) == Status::Ok);
    }

    // Query a sub-window [0.3s, 1.1s] -> duration 0.8s. Analytic: exp(xi * 0.8).
    // For a constant twist the per-sample composition is EXACT (xi commutes with
    // itself: exp(xi*dt)^n = exp(xi*n*dt)), so the only residual is the endpoint
    // interpolation at the non-sample window bounds -> the small tolerances below.
    const Delta d = b.query(secs(0.3), secs(1.1)).value();
    const SE3 expected = se3::exp(xi * 0.8);
    CHECK(close(d.motion.R, expected.R, 1e-6));
    CHECK(close(d.motion.t, expected.t, 1e-5));
    CHECK(d.t0 == secs(0.3));
    CHECK(d.t1 == secs(1.1));
}

// ---------------------------------------------------------------------------
// Increment composition
// ---------------------------------------------------------------------------
TEST_CASE("increment composition: query spans the composed product") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 16, ConfCombine::NativeOnly);

    SE3 i1; i1.t = Vec3(1.0, 0.0, 0.0);
    SE3 i2; i2.R = so3::exp(Vec3(0.0, 0.0, kPi / 2)); i2.t = Vec3(0.0, 1.0, 0.0);
    SE3 i3; i3.t = Vec3(0.5, 0.0, 0.0);

    SE3 origin;  // identity start
    b.push_increment(secs(1.0), origin, Mat6::Zero());     // seed at identity
    b.push_increment(secs(2.0), i1, Mat6::Identity());
    b.push_increment(secs(3.0), i2, Mat6::Identity());
    b.push_increment(secs(4.0), i3, Mat6::Identity());

    // Whole span [1,4] = i1 o i2 o i3.
    const Delta d = b.query(secs(1.0), secs(4.0)).value();
    const SE3 expected = se3::compose(se3::compose(i1, i2), i3);
    CHECK(close(d.motion.R, expected.R, 1e-9));
    CHECK(close(d.motion.t, expected.t, 1e-9));
}

// ---------------------------------------------------------------------------
// Absolute pose differencing
// ---------------------------------------------------------------------------
TEST_CASE("absolute pose: delta is the difference of two absolute poses") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::AbsolutePose, 16, ConfCombine::NativeOnly);

    SE3 p0; p0.t = Vec3(1.0, 2.0, 0.0);
    SE3 p1; p1.R = so3::exp(Vec3(0, 0, kPi / 2)); p1.t = Vec3(3.0, 2.0, 0.0);
    SE3 p2; p2.R = so3::exp(Vec3(0, 0, kPi));      p2.t = Vec3(3.0, 5.0, 0.0);

    b.push_absolute(secs(1.0), p0, Mat6::Identity());
    b.push_absolute(secs(2.0), p1, Mat6::Identity());
    b.push_absolute(secs(3.0), p2, Mat6::Identity());

    const Delta d = b.query(secs(1.0), secs(3.0)).value();
    const SE3 expected = se3::compose(se3::inverse(p0), p2);
    CHECK(close(d.motion.R, expected.R, 1e-9));
    CHECK(close(d.motion.t, expected.t, 1e-9));
}

// ---------------------------------------------------------------------------
// Interpolation accuracy at non-sample timestamps
// ---------------------------------------------------------------------------
TEST_CASE("query interpolates cum_pose at non-sample timestamps") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Twist, 64, ConfCombine::NativeOnly);

    // Pure straight motion 1 m/s along x, sampled coarsely at 1 Hz.
    Vec6 xi; xi << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    for (int k = 0; k <= 5; ++k) {
        b.push_twist(static_cast<Timestamp>(k) * kSec, xi, Mat6::Zero());
    }

    // [0.5s, 3.25s] -> translation 2.75 m along x (straight line lerp is exact).
    const Delta d = b.query(secs(0.5), secs(3.25)).value();
    CHECK(d.motion.t.x() == doctest::Approx(2.75).epsilon(1e-9));
    CHECK(std::abs(d.motion.t.y()) < 1e-9);
    CHECK(close(d.motion.R, Mat3::Identity(), 1e-9));
}

TEST_CASE("interpolation of a pure rotation is exact at the midpoint") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::AbsolutePose, 16, ConfCombine::NativeOnly);

    SE3 a;                                  // identity
    SE3 c; c.R = so3::exp(Vec3(0, 0, kPi / 2));   // 90 deg about z
    b.push_absolute(secs(0.0), a, Mat6::Zero());
    b.push_absolute(secs(2.0), c, Mat6::Zero());

    // Midpoint t=1.0 should be 45 deg about z; delta [0,1] == 45 deg.
    const Delta d = b.query(secs(0.0), secs(1.0)).value();
    const Mat3 expected = so3::exp(Vec3(0, 0, kPi / 4));
    CHECK(close(d.motion.R, expected, 1e-9));
}

// ---------------------------------------------------------------------------
// Covariance grows with window length (additive native model)
// ---------------------------------------------------------------------------
TEST_CASE("native covariance accumulates with window length") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 64, ConfCombine::NativeOnly);

    SE3 step; step.t = Vec3(1.0, 0.0, 0.0);
    b.push_increment(secs(0.0), SE3{}, Mat6::Zero());     // seed
    for (int k = 1; k <= 10; ++k) {
        b.push_increment(static_cast<Timestamp>(k) * kSec, step, Mat6::Identity());
    }

    const Delta d_short = b.query(secs(0.0), secs(2.0)).value();   // 2 increments
    const Delta d_long  = b.query(secs(0.0), secs(6.0)).value();   // 6 increments
    const Scalar tr_short = d_short.cov.trace();
    const Scalar tr_long  = d_long.cov.trace();
    CHECK(tr_long > tr_short);
    // Each unit increment contributes Identity (trace 6); 2 vs 6 increments.
    CHECK(tr_short == doctest::Approx(12.0));
    CHECK(tr_long  == doctest::Approx(36.0));
}

TEST_CASE("fractional end increments scale the accumulated covariance") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 64, ConfCombine::NativeOnly);

    SE3 step; step.t = Vec3(1.0, 0.0, 0.0);
    b.push_increment(secs(0.0), SE3{}, Mat6::Zero());     // seed at t=0
    b.push_increment(secs(1.0), step, Mat6::Identity());  // increment over (0,1]
    b.push_increment(secs(2.0), step, Mat6::Identity());  // increment over (1,2]

    // Window [0.5, 1.5] overlaps half of each increment -> total ~ Identity.
    const Delta d = b.query(secs(0.5), secs(1.5)).value();
    CHECK(d.cov.trace() == doctest::Approx(6.0));
}

// ---------------------------------------------------------------------------
// native (+) modeled combine for each ConfCombine
// ---------------------------------------------------------------------------
TEST_CASE("native(+)modeled combine honors every ConfCombine rule") {
    // One increment of 1 m along x, no rotation. native = 2*Identity.
    // modeled: per_m = 0.5 over 1 m -> trans var = 0.25; per_rad = 0 -> rot var 0.
    SE3 step; step.t = Vec3(1.0, 0.0, 0.0);
    const Mat6 native = 2.0 * Mat6::Identity();

    auto build = [&](ConfCombine combine, Scalar blend, bool native_conf) {
        SourceBuffer b;
        b.configure(make_sensor(0.5, 0.0, native_conf), OdomForm::Increment, 8,
                    combine, blend);
        b.push_increment(secs(0.0), SE3{}, Mat6::Zero());
        b.push_increment(secs(1.0), step, native);
        return b.query(secs(0.0), secs(1.0)).value().cov;
    };

    // Expected modeled translation variance (per-m)^2 * mag^2 = 0.25; rot = 0.
    // Implementation floors zero-magnitude axes with a tiny epsilon; check only
    // the load-bearing translation block (index 0).
    const Scalar nat_t = 2.0;        // native trans diagonal
    const Scalar mod_t = 0.25;       // modeled trans diagonal

    {
        const Mat6 c = build(ConfCombine::NativeOnly, 0.5, true);
        CHECK(c(0, 0) == doctest::Approx(nat_t));
    }
    {
        const Mat6 c = build(ConfCombine::ModeledOnly, 0.5, true);
        CHECK(c(0, 0) == doctest::Approx(mod_t));
    }
    {
        const Mat6 c = build(ConfCombine::Sum, 0.5, true);
        CHECK(c(0, 0) == doctest::Approx(nat_t + mod_t));
    }
    {
        const Mat6 c = build(ConfCombine::Max, 0.5, true);
        CHECK(c(0, 0) == doctest::Approx(nat_t));   // max(2, 0.25)
    }
    {
        const Mat6 c = build(ConfCombine::Weighted, 0.25, true);
        // blend*native + (1-blend)*modeled = 0.25*2 + 0.75*0.25 = 0.6875
        CHECK(c(0, 0) == doctest::Approx(0.25 * nat_t + 0.75 * mod_t));
    }
    {
        // native_confidence = false -> per D7 the missing native Sigma is the
        // identity, then the configured rule still runs: Sum -> Identity + modeled.
        const Mat6 c = build(ConfCombine::Sum, 0.5, false);
        CHECK(c(0, 0) == doctest::Approx(1.0 + mod_t));
    }
}

// ---------------------------------------------------------------------------
// D7 missing-native: absent native Sigma is the IDENTITY, then the configured
// ConfCombine rule runs (NOT modeled-only). Exercise native-present vs
// native_confidence=false for at least Sum and Max.
// ---------------------------------------------------------------------------
TEST_CASE("missing native covariance combines identity (not modeled-only) per D7") {
    // One increment of 1 m along x, no rotation. native = 2*Identity (when present).
    // modeled: per_m = 0.5 over 1 m -> trans var 0.25; per_rad = 0 -> rot var ~0.
    SE3 step; step.t = Vec3(1.0, 0.0, 0.0);
    const Mat6 native = 2.0 * Mat6::Identity();

    auto build = [&](ConfCombine combine, bool native_conf) {
        SourceBuffer b;
        b.configure(make_sensor(0.5, 0.0, native_conf), OdomForm::Increment, 8,
                    combine);
        b.push_increment(secs(0.0), SE3{}, Mat6::Zero());
        b.push_increment(secs(1.0), step, native);
        return b.query(secs(0.0), secs(1.0)).value().cov;
    };

    const Scalar nat_t = 2.0;     // native trans diagonal (present)
    const Scalar mod_t = 0.25;    // modeled trans diagonal
    const Scalar id_t  = 1.0;     // identity native operand when absent

    // --- Sum --------------------------------------------------------------
    {
        // native present: native + modeled.
        const Mat6 c = build(ConfCombine::Sum, true);
        CHECK(c(0, 0) == doctest::Approx(nat_t + mod_t));
    }
    {
        // native absent (native_confidence=false): Identity + modeled, NOT modeled.
        const Mat6 c = build(ConfCombine::Sum, false);
        CHECK(c(0, 0) == doctest::Approx(id_t + mod_t));
        // The identity native operand also shows through the rotation block, where
        // the modeled rot var is ~0 here -> ~1 from Identity alone.
        CHECK(c(3, 3) == doctest::Approx(id_t));
    }

    // --- Max --------------------------------------------------------------
    {
        // native present: max(2, 0.25) = 2.
        const Mat6 c = build(ConfCombine::Max, true);
        CHECK(c(0, 0) == doctest::Approx(nat_t));
    }
    {
        // native absent: max(Identity, modeled) = max(1, 0.25) = 1, NOT modeled.
        const Mat6 c = build(ConfCombine::Max, false);
        CHECK(c(0, 0) == doctest::Approx(id_t));
    }
}

TEST_CASE("modeled rotation noise scales with rotation magnitude") {
    SourceBuffer b;
    b.configure(make_sensor(0.0, 0.5, false), OdomForm::Increment, 8,
                ConfCombine::ModeledOnly);
    SE3 rot; rot.R = so3::exp(Vec3(0, 0, 1.0));   // 1 rad about z
    b.push_increment(secs(0.0), SE3{}, Mat6::Zero());
    b.push_increment(secs(1.0), rot, Mat6::Zero());

    const Mat6 c = b.query(secs(0.0), secs(1.0)).value().cov;
    // rot variance axis (index 3..5): (per_rad * 1 rad)^2 = 0.25.
    CHECK(c(3, 3) == doctest::Approx(0.25));
}

// ---------------------------------------------------------------------------
// covers() / NoData / OutOfRange
// ---------------------------------------------------------------------------
TEST_CASE("covers and query edge cases") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 16, ConfCombine::NativeOnly);
    b.push_increment(secs(1.0), SE3{}, Mat6::Zero());
    SE3 step; step.t = Vec3(1, 0, 0);
    b.push_increment(secs(2.0), step, Mat6::Identity());
    b.push_increment(secs(3.0), step, Mat6::Identity());

    CHECK(b.covers(secs(1.5), secs(2.5)));
    CHECK(b.covers(secs(1.0), secs(3.0)));   // exact endpoints
    CHECK_FALSE(b.covers(secs(0.5), secs(2.0)));  // before oldest
    CHECK_FALSE(b.covers(secs(2.0), secs(3.5)));  // after newest
    CHECK_FALSE(b.covers(secs(2.0), secs(2.0)));  // empty window

    // t1 <= t0 -> OutOfRange
    CHECK(b.query(secs(2.0), secs(2.0)).status() == Status::OutOfRange);
    CHECK(b.query(secs(2.5), secs(1.5)).status() == Status::OutOfRange);
    // uncovered -> NoData
    CHECK(b.query(secs(0.5), secs(2.0)).status() == Status::NoData);
    CHECK(b.query(secs(2.0), secs(9.0)).status() == Status::NoData);
}

TEST_CASE("query on empty buffer returns NoData") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 8, ConfCombine::NativeOnly);
    CHECK(b.query(secs(0.0), secs(1.0)).status() == Status::NoData);
}

// ---------------------------------------------------------------------------
// Ring-buffer wrap: pushing beyond capacity drops the oldest
// ---------------------------------------------------------------------------
TEST_CASE("ring buffer wraps and drops the oldest entry") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 4, ConfCombine::NativeOnly);

    SE3 step; step.t = Vec3(1, 0, 0);
    for (int k = 0; k < 7; ++k) {
        b.push_increment(static_cast<Timestamp>(k) * kSec, step, Mat6::Identity());
    }
    // Capacity 4 -> only stamps 3,4,5,6 remain.
    CHECK(b.size() == 4);
    CHECK(b.oldest() == secs(3.0));
    CHECK(b.newest() == secs(6.0));
    CHECK_FALSE(b.covers(secs(2.0), secs(5.0)));   // 2.0 was dropped
    CHECK(b.covers(secs(3.0), secs(6.0)));
}

TEST_CASE("reset clears entries but keeps configuration") {
    SourceBuffer b;
    b.configure(make_sensor(), OdomForm::Increment, 8, ConfCombine::NativeOnly);
    b.push_increment(secs(1.0), SE3{}, Mat6::Zero());
    b.push_increment(secs(2.0), SE3{}, Mat6::Identity());
    CHECK(b.size() == 2);

    b.reset();
    CHECK(b.size() == 0);
    CHECK(b.configured());
    CHECK(b.oldest() == 0);
    // can push again after reset
    CHECK(b.push_increment(secs(5.0), SE3{}, Mat6::Zero()) == Status::Ok);
}
