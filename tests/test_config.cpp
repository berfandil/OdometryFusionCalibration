// Slice 0 unit tests: Config::validate() bounds checking.
#include <doctest/doctest.h>

#include "ofc/core/config.hpp"

using namespace ofc;

namespace {
Config good() {
    Config c;                 // defaults are valid
    return c;
}
} // namespace

TEST_CASE("validate accepts default config") {
    CHECK(validate(good()) == Status::Ok);
}

TEST_CASE("validate rejects bad source counts") {
    Config c = good();
    c.max_sources = 0;
    CHECK(validate(c) == Status::OutOfRange);

    c = good();
    c.max_sources = 4;
    c.sensor_count = 5;       // exceeds max_sources
    CHECK(validate(c) == Status::CapacityExceeded);
}

TEST_CASE("validate rejects non-positive rates/windows") {
    Config c = good();
    c.tick_rate_hz = 0.0;
    CHECK(validate(c) == Status::OutOfRange);

    c = good();
    c.window_s = -1.0;
    CHECK(validate(c) == Status::OutOfRange);
}

TEST_CASE("validate rejects degenerate median/metric") {
    Config c = good();
    c.weiszfeld_max_iters = 0;
    CHECK(validate(c) == Status::OutOfRange);

    c = good();
    c.metric_lambda = 0.0;
    CHECK(validate(c) == Status::OutOfRange);
}

TEST_CASE("validate rejects inverted hysteresis") {
    Config c = good();
    c.commit_drop = c.commit_concentration;   // must be strictly less
    CHECK(validate(c) == Status::InvalidConfig);
}

TEST_CASE("validate rejects out-of-range reference sensor") {
    Config c = good();
    c.reference_sensor_id = 99;               // >= max_sources
    CHECK(validate(c) == Status::OutOfRange);
}

TEST_CASE("validate range-checks fusion/median knobs") {
    // confidence_blend must be in [0, 1].
    Config c = good();
    c.confidence_blend = 1.5;
    CHECK(validate(c) == Status::OutOfRange);
    c = good();
    c.confidence_blend = -0.1;
    CHECK(validate(c) == Status::OutOfRange);

    // weiszfeld_tol must be in (0, 1).
    c = good();
    c.weiszfeld_tol = 0.0;
    CHECK(validate(c) == Status::OutOfRange);
    c = good();
    c.weiszfeld_tol = 1.0;
    CHECK(validate(c) == Status::OutOfRange);

    // weight_cap must be >= 1.
    c = good();
    c.weight_cap = 0.5;
    CHECK(validate(c) == Status::OutOfRange);

    // tip_cov_inflation must be >= 1.
    c = good();
    c.tip_cov_inflation = 0.9;
    CHECK(validate(c) == Status::OutOfRange);

    // q_scale must be > 0; q_floor entries must be >= 0.
    c = good();
    c.q_scale = 0.0;
    CHECK(validate(c) == Status::OutOfRange);
    c = good();
    c.q_floor[3] = -1e-6;
    CHECK(validate(c) == Status::OutOfRange);

    // fusion_delay_s must be in [0, 2].
    c = good();
    c.fusion_delay_s = 2.5;
    CHECK(validate(c) == Status::OutOfRange);
}
