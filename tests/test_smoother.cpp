// Slice 10 unit tests: per-sensor fixed-lag RTS twist smoother (D18).
//
// These pin the TWO headline D18 properties (and are deliberately DISCRIMINATING — each
// must fail for the wrong filter):
//   * VARIANCE REDUCTION (peak sharpening): on a known twist + seeded Gaussian noise the
//     smoothed-output variance about the truth is materially LOWER (< 0.6x) than the raw
//     input variance. (A no-op pass-through would read ~1.0x and fail.)
//   * ZERO PHASE (no peak shift): on a known TIME-VARYING twist (a ramp) the smoothed
//     estimate at the EMITTED (lag-L) timestamp matches the TRUE twist at THAT timestamp —
//     a two-sided RTS does NOT lag. A purely CAUSAL filter (forward-KF / EMA) WOULD lag a
//     ramp by ~L·slope; the test asserts the RTS error is far below that causal-lag floor
//     (so it fails for a causal filter).
//   * DETERMINISM (same input -> same output) + a NOT-READY guard (before the ring fills).
#include <doctest/doctest.h>

#include "ofc/core/lie.hpp"
#include "ofc/core/smoother.hpp"

#include <cmath>
#include <cstdint>
#include <memory>

using namespace ofc;

namespace {
// A tiny deterministic Gaussian PRNG (Box-Muller on a 64-bit LCG) — no <random> ordering
// surprises across platforms; seeded so tests are bit-reproducible.
struct Rng {
    std::uint64_t s;
    explicit Rng(std::uint64_t seed) : s(seed ? seed : 0x9E3779B97F4A7C15ull) {}
    Scalar uniform() {                 // (0, 1)
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        const std::uint64_t x = (s >> 11);
        return (static_cast<Scalar>(x) + 1.0) / (static_cast<Scalar>(1ull << 53) + 2.0);
    }
    Scalar normal() {                  // N(0, 1)
        const Scalar u1 = uniform();
        const Scalar u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979323846 * u2);
    }
};

std::unique_ptr<TwistSmoother> make_smoother(int lag, Scalar q, Scalar r = 1.0) {
    std::unique_ptr<TwistSmoother> sm(new TwistSmoother());
    REQUIRE(ok(sm->configure(1, lag, q, r)));
    return sm;
}
} // namespace

// ---------------------------------------------------------------------------
// Config / degenerate guards
// ---------------------------------------------------------------------------
TEST_CASE("smoother: configure clamps lag to [1, kMaxLag] and rejects bad args") {
    // Heap-allocate: the fixed-capacity rings make TwistSmoother too large for the test
    // stack (same pattern as the calibrator/TimeSync tests). It lives in the estimator's
    // heap-allocated Impl in production.
    std::unique_ptr<TwistSmoother> smp(new TwistSmoother());
    TwistSmoother& sm = *smp;
    CHECK(ok(sm.configure(1, 5, 1.0)));
    CHECK(sm.lag_steps() == 5);

    CHECK(ok(sm.configure(1, 0, 1.0)));      // L=0 clamps UP to 1 (no degenerate pass)
    CHECK(sm.lag_steps() == 1);

    CHECK(ok(sm.configure(1, 10000, 1.0)));
    CHECK(sm.lag_steps() == TwistSmoother::kMaxLag);   // clamped DOWN to the cap

    CHECK(sm.configure(0, 5, 1.0) == Status::OutOfRange);        // max_sources < 1
    CHECK(sm.configure(99, 5, 1.0) == Status::OutOfRange);       // > kMaxSources
    CHECK(sm.configure(1, 5, 0.0) == Status::OutOfRange);        // q <= 0
    CHECK(sm.configure(1, 5, 1.0, 0.0) == Status::OutOfRange);   // r <= 0
}

TEST_CASE("smoother: not ready until lag_steps+1 pushes, then ready") {
    const int L = 4;
    auto sm = make_smoother(L, 1.0);
    Vec6 z = Vec6::Zero(); z(0) = 1.0;     // constant twist
    const Scalar dt = 0.02;

    // The first L pushes leave the ring unfilled -> not ready.
    for (int i = 0; i < L; ++i) {
        CHECK(ok(sm->push(0, z, dt)));
        CHECK_FALSE(sm->ready(0));
    }
    // The (L+1)-th push fills the ring -> ready; the emitted (lag-L) twist is the smoothed
    // estimate of the OLDEST entry, which on a constant signal is the constant.
    CHECK(ok(sm->push(0, z, dt)));
    CHECK(sm->ready(0));
    CHECK((sm->smoothed(0) - z).norm() < 1e-6);
}

TEST_CASE("smoother: out-of-range slot / unconfigured guards") {
    std::unique_ptr<TwistSmoother> smp(new TwistSmoother());
    TwistSmoother& sm = *smp;
    Vec6 z = Vec6::Zero();
    CHECK(sm.push(0, z, 0.02) == Status::NotInitialized);   // before configure
    REQUIRE(ok(sm.configure(2, 3, 1.0)));
    CHECK(sm.push(-1, z, 0.02) == Status::OutOfRange);
    CHECK(sm.push(2, z, 0.02) == Status::OutOfRange);       // slot == max_sources
    CHECK(ok(sm.push(0, z, 0.02)));
    CHECK_FALSE(sm.ready(5));                                // bad slot -> not ready
    CHECK(sm.smoothed(5).norm() == doctest::Approx(0.0));    // bad slot -> zero
}

// ---------------------------------------------------------------------------
// HEADLINE 1 — variance reduction (peak sharpening)
// ---------------------------------------------------------------------------
TEST_CASE("smoother: smoothed-output variance about truth is materially below the raw "
          "input variance (peak sharpening)") {
    const int L = 8;
    const Scalar q = 4.0;          // accel random-walk; tuned so the constant-truth fits
    const Scalar r = 1.0;
    auto sm = make_smoother(L, q, r);

    // TRUTH: a constant body twist (v along x, a steady yaw rate). Inject seeded Gaussian
    // noise per axis; compare the RAW measurement scatter to the SMOOTHED scatter, both
    // measured about the SAME truth at the emitted (lag-L) timestamp.
    Vec6 truth; truth << 2.0, 0.0, 0.0, 0.0, 0.0, 0.5;
    const Scalar sigma = 0.3;      // per-axis measurement noise
    const Scalar dt = 0.02;
    Rng rng(424242u);

    Scalar raw_sse = 0.0, sm_sse = 0.0;
    int raw_n = 0, sm_n = 0;
    const int steps = 4000;
    bool all_ok = true;   // DIGEST: fold the per-push status into one assertion (budget).
    for (int i = 0; i < steps; ++i) {
        Vec6 z = truth;
        for (int c = 0; c < 6; ++c) z(c) += sigma * rng.normal();
        // Score the RAW noise of every measurement about truth.
        raw_sse += (z - truth).squaredNorm();
        ++raw_n;

        all_ok = all_ok && ok(sm->push(0, z, dt));
        if (sm->ready(0)) {
            // Discard the first few ready outputs (filter transient) before scoring (the
            // emitted output corresponds to truth at the lag-L sample = `truth`, constant).
            if (i > L + 200) {
                sm_sse += (sm->smoothed(0) - truth).squaredNorm();
                ++sm_n;
            }
        }
    }
    REQUIRE(all_ok);
    REQUIRE(raw_n > 1000);
    REQUIRE(sm_n > 1000);
    const Scalar raw_var = raw_sse / static_cast<Scalar>(raw_n);
    const Scalar sm_var  = sm_sse  / static_cast<Scalar>(sm_n);
    // The two-sided pass over L=8 steps must cut the variance to well under 0.6x the raw.
    CHECK(sm_var < 0.6 * raw_var);
    // Sanity: it actually reduces (not a no-op), and is positive (not collapsed to truth).
    CHECK(sm_var > 0.0);
    CHECK(sm_var < raw_var);
}

// ---------------------------------------------------------------------------
// HEADLINE 2 — zero phase (no peak shift). DISCRIMINATING vs a causal filter.
// ---------------------------------------------------------------------------
TEST_CASE("smoother: a time-varying (ramp) twist is tracked with ZERO phase lag at the "
          "emitted timestamp (a causal filter would lag)") {
    const int L = 10;
    const Scalar q = 50.0;     // high process noise -> the filter trusts the kinematics,
    const Scalar r = 1.0;      // so a constant-slope ramp is tracked exactly by the CV model
    auto sm = make_smoother(L, q, r);

    // TRUTH: a linear ramp in the yaw-rate channel (index 5): w5(t) = slope * t. The CV
    // model's accel state captures a constant slope, so a NOISELESS ramp is reconstructed
    // exactly by the two-sided pass — at the EMITTED (lag-L) timestamp, with no lag.
    const Scalar slope = 1.0;        // rad/s per second
    const Scalar dt = 0.02;
    const int steps = 2000;

    Scalar max_phase_err = 0.0;
    Scalar causal_lag_ref = 0.0;     // the lag a perfect causal filter would show: L*slope*dt
    bool all_ok = true;              // DIGEST: fold per-push status into one assertion.
    for (int i = 0; i < steps; ++i) {
        const Scalar t_now = static_cast<Scalar>(i) * dt;
        Vec6 z = Vec6::Zero();
        z(5) = slope * t_now;        // noiseless ramp
        all_ok = all_ok && ok(sm->push(0, z, dt));
        if (sm->ready(0) && i > L + 50) {
            // The emitted output corresponds to the sample L steps in the past:
            const Scalar t_emitted = static_cast<Scalar>(i - L) * dt;
            const Scalar truth_emitted = slope * t_emitted;
            const Scalar got = sm->smoothed(0)(5);
            max_phase_err = std::max(max_phase_err, std::abs(got - truth_emitted));
        }
    }
    REQUIRE(all_ok);
    causal_lag_ref = static_cast<Scalar>(L) * slope * dt;   // = 0.2 rad/s
    REQUIRE(causal_lag_ref > 0.1);
    // The two-sided RTS reconstructs the ramp at the emitted timestamp with NO lag: the
    // phase error is FAR below the causal-lag floor (a forward-only filter evaluated at the
    // same emitted time would be off by ~causal_lag_ref). This is the discriminating check.
    CHECK(max_phase_err < 0.1 * causal_lag_ref);
}

TEST_CASE("smoother: DISCRIMINATING — the smoother's OWN causal (forward-only) output DOES "
          "lag the ramp (guards the zero-phase test from being vacuous)") {
    // We read the SMOOTHER'S OWN one-sided causal estimate (its latest FORWARD-filtered twist)
    // and show it lags — rather than comparing an analytic causal reference (the old version's
    // `slope*t_now`, a tautology about the ramp). We obtain the causal output WITHOUT a new
    // accessor: smoothed() is documented to fall back to the latest filtered twist BEFORE the
    // ring fills (the causal pass-through). So we run a SHORT-lag (Lc=2) smoother and read its
    // causal warm-up value on the same ramp: at push i (during its brief warm-up) smoothed()
    // is the latest filtered twist, which tracks ~t_now; compared against the truth Lc steps in
    // the PAST (the sample a two-sided estimate would report) it lags by ~Lc*slope*dt. The
    // genuine lag-L two-sided smoother does NOT lag at its emitted timestamp (re-checked here),
    // so the causal lag is a real property of THIS smoother's forward filter, not an artifact.
    const int L = 10;            // the genuine two-sided lag (matches the zero-phase test)
    const int Lc = 2;            // a short lag whose warm-up exposes the causal fallback
    const Scalar slope = 1.0, dt = 0.02;
    const int steps = 2000;
    auto sm_lag    = make_smoother(L,  50.0, 1.0);
    auto sm_causal = make_smoother(Lc, 50.0, 1.0);

    // Causal-fallback readout: the Lc smoother's latest filtered twist at its LAST warm-up push
    // (push index Lc-1, still not ready), compared to the truth Lc steps in the past.
    Scalar causal_err = 0.0;
    bool causal_scored = false;
    Scalar lag_err_at_emitted = 0.0;
    int scored = 0;
    bool all_ok = true;              // DIGEST: fold per-push status into one assertion.
    for (int i = 0; i < steps; ++i) {
        const Scalar t_now = static_cast<Scalar>(i) * dt;
        Vec6 z = Vec6::Zero(); z(5) = slope * t_now;
        all_ok = all_ok && ok(sm_lag->push(0, z, dt));
        if (i <= Lc) {
            all_ok = all_ok && ok(sm_causal->push(0, z, dt));
            if (i == Lc - 1) {
                // Last warm-up push of the Lc smoother: smoothed() == latest FILTERED twist
                // (causal fallback), tracking ~t_now. The two-sided sample at this point would
                // be the truth Lc steps in the past; the causal readout lags it.
                all_ok = all_ok && !sm_causal->ready(0);
                const Scalar truth_past = slope * static_cast<Scalar>(i - Lc) * dt;
                causal_err = std::abs(sm_causal->smoothed(0)(5) - truth_past);
                causal_scored = true;
            }
        }
        if (sm_lag->ready(0) && i > L + 50) {
            const Scalar t_emitted     = static_cast<Scalar>(i - L) * dt;
            const Scalar truth_emitted = slope * t_emitted;
            lag_err_at_emitted = std::max(lag_err_at_emitted,
                                          std::abs(sm_lag->smoothed(0)(5) - truth_emitted));
            ++scored;
        }
    }
    REQUIRE(all_ok);
    REQUIRE(causal_scored);
    REQUIRE(scored > 100);
    // The smoother's OWN causal readout lags the (Lc-steps-past) truth by ~Lc*slope*dt ...
    CHECK(causal_err > 0.5 * (static_cast<Scalar>(Lc) * slope * dt));
    // ... while its genuine two-sided lag-L output does NOT lag at its emitted timestamp.
    CHECK(lag_err_at_emitted < 0.05);
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------
TEST_CASE("smoother: identical input stream -> bit-identical smoothed output") {
    auto run = [](Vec6* out, int outcap) {
        auto sm = make_smoother(6, 3.0, 1.0);
        Rng rng(13579u);
        int produced = 0;
        for (int i = 0; i < 500 && produced < outcap; ++i) {
            Vec6 z; z << 1.0, 0.0, 0.0, 0.0, 0.0, 0.3;
            for (int c = 0; c < 6; ++c) z(c) += 0.2 * rng.normal();
            sm->push(0, z, 0.02);
            if (sm->ready(0)) out[produced++] = sm->smoothed(0);
        }
        return produced;
    };
    Vec6 a[256], b[256];
    const int na = run(a, 256);
    const int nb = run(b, 256);
    REQUIRE(na == nb);
    REQUIRE(na > 100);
    for (int i = 0; i < na; ++i) {
        CHECK((a[i].array() == b[i].array()).all());
    }
}

TEST_CASE("smoother: multiple independent source slots do not cross-contaminate") {
    auto sm = make_smoother(5, 2.0, 1.0);
    REQUIRE(ok(sm->configure(2, 5, 2.0, 1.0)));
    Vec6 za; za << 2.0, 0, 0, 0, 0, 0.5;
    Vec6 zb; zb << -1.0, 0, 0, 0, 0, -0.3;
    bool all_ok = true;              // DIGEST: fold per-push status into one assertion.
    for (int i = 0; i < 100; ++i) {
        all_ok = all_ok && ok(sm->push(0, za, 0.02));
        all_ok = all_ok && ok(sm->push(1, zb, 0.02));
    }
    REQUIRE(all_ok);
    REQUIRE(sm->ready(0));
    REQUIRE(sm->ready(1));
    CHECK((sm->smoothed(0) - za).norm() < 1e-5);
    CHECK((sm->smoothed(1) - zb).norm() < 1e-5);
}
