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
