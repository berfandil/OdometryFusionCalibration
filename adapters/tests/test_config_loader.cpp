// Adapters Slice 13: ConfigLoader. A dependency-free key=value / INI text -> Config, run
// through the core validate(). Covers: a sample config parses to the expected fields + validates
// Ok; the per-sensor SensorConfig storage is owned by the loader (ownership contract); bad input
// (unknown key, malformed value, validate() rejection) is reported as a clear Status.
#include <doctest/doctest.h>

#include "ofc_adapters/config_loader.hpp"

#include "ofc/core/config.hpp"
#include "ofc/core/status.hpp"

#include <cmath>
#include <string>
#include <vector>

using namespace ofc;
using namespace ofc::adapters;

TEST_CASE("ConfigLoader: a sample config parses to the expected fields and validates Ok") {
    const std::string text =
        "# a sample rig config (comments + blank lines ignored)\n"
        "\n"
        "[global]\n"
        "max_sources = 4\n"
        "reference_sensor_id = 0\n"
        "window_s = 0.12\n"
        "fusion_delay_s = 0.04\n"
        "tick_rate_hz = 100\n"
        "calib_lag_s = 0.25\n"
        "cold_start = median_from_start    ; per-DOF cold start\n"
        "vote_weight = one\n"
        "commit_concentration = 0.7\n"
        "commit_min_votes = 150\n"
        "commit_drop = 0.4\n"
        "straight_omega_max = 0.05\n"
        "straight_trans_min = 0.03\n"
        "turn_omega_min = 0.2\n"
        "timesync_enabled = false\n"
        "\n"
        "[sensor.0]\n"
        "id = 0\n"
        "is_reference = true\n"
        "prior_scale = 1.0\n"
        "\n"
        "[sensor.1]\n"
        "id = 1\n"
        "prior_extrinsic = 0.1 0.0 0.0  0.5 0.0 0.2\n"   // yaw pitch roll x y z
        "prior_scale = 1.2\n"
        "prior_time_offset_s = 0.03\n"
        "weight_prior = 2.0\n"
        "per_sensor_kf = on\n"
        "scale_calib = yes\n";

    ConfigLoader loader;
    REQUIRE(loader.parse(text) == Status::Ok);
    const Config& c = loader.config();

    // Top-level knobs.
    CHECK(c.max_sources == 4);
    CHECK(c.reference_sensor_id == 0);
    CHECK(c.window_s == doctest::Approx(0.12));
    CHECK(c.fusion_delay_s == doctest::Approx(0.04));
    CHECK(c.tick_rate_hz == doctest::Approx(100.0));
    CHECK(c.calib_lag_s == doctest::Approx(0.25));
    CHECK(c.cold_start == ColdStart::MedianFromStart);
    CHECK(c.vote_weight == VoteWeight::One);
    CHECK(c.commit_concentration == doctest::Approx(0.7));
    CHECK(c.commit_min_votes == 150);
    CHECK(c.commit_drop == doctest::Approx(0.4));
    CHECK(c.turn_omega_min == doctest::Approx(0.2));
    CHECK_FALSE(c.timesync_enabled);

    // Per-sensor: owned storage, pointed-to by Config::sensors.
    REQUIRE(c.sensor_count == 2);
    REQUIRE(c.sensors != nullptr);
    CHECK(c.sensors == loader.sensors().data());     // ownership: Config points at the loader's vec
    CHECK(c.sensors[0].id == 0);
    CHECK(c.sensors[0].is_reference);
    CHECK(c.sensors[1].id == 1);
    CHECK(c.sensors[1].prior_scale == doctest::Approx(1.2));
    CHECK(c.sensors[1].prior_time_offset_s == doctest::Approx(0.03));
    CHECK(c.sensors[1].weight_prior == doctest::Approx(2.0));
    CHECK(c.sensors[1].per_sensor_kf);
    CHECK(c.sensors[1].scale_calib);
    // prior_extrinsic: yaw 0.1 about Z -> R(1,0) = sin(0.1); t = (0.5, 0, 0.2).
    CHECK(c.sensors[1].prior_extrinsic.R(1, 0) == doctest::Approx(std::sin(0.1)));
    CHECK(c.sensors[1].prior_extrinsic.t.x() == doctest::Approx(0.5));
    CHECK(c.sensors[1].prior_extrinsic.t.z() == doctest::Approx(0.2));
}

TEST_CASE("ConfigLoader ownership: the built Config outlives via the loader's owned vector") {
    const std::string text =
        "[sensor.0]\nid = 0\nis_reference = true\n"
        "[sensor.1]\nid = 1\n";
    ConfigLoader loader;
    REQUIRE(loader.parse(text) == Status::Ok);
    // The Config's sensors pointer is the loader's vector data — stable for the loader's life.
    const SensorConfig* p = loader.config().sensors;
    CHECK(p == loader.sensors().data());
    CHECK(loader.config().sensor_count == 2);
}

TEST_CASE("ConfigLoader: q_scale / q_floor / adaptive_q (the real-data covariance enabler)") {
    // Scalar q_floor -> applied to ALL 6 axes; adaptive_q + q_scale parsed.
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "adaptive_q = false\nq_scale = 0.5\nq_floor = 0.01\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        const Config& c = loader.config();
        CHECK_FALSE(c.adaptive_q);
        CHECK(c.q_scale == doctest::Approx(0.5));
        for (int i = 0; i < 6; ++i) CHECK(c.q_floor[i] == doctest::Approx(0.01));
    }
    // 6-number q_floor -> per-axis [trans; rot].
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "q_floor = 0.1 0.2 0.3 0.4 0.5 0.6\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        const Scalar* qf = loader.config().q_floor;
        CHECK(qf[0] == doctest::Approx(0.1));
        CHECK(qf[3] == doctest::Approx(0.4));
        CHECK(qf[5] == doctest::Approx(0.6));
    }
    // Malformed q_floor (2 numbers, neither 1 nor 6) -> clear error, not Ok.
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "q_floor = 0.1 0.2\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        CHECK(loader.parse(text) != Status::Ok);
        CHECK(loader.error().find("q_floor") != std::string::npos);
    }
}

TEST_CASE("ConfigLoader: rot3d_enabled parses (Slice 17); default off; bad value rejects") {
    // Enabled.
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "rot3d_enabled = true\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        CHECK(loader.config().rot3d_enabled);
    }
    // Key absent -> default OFF (byte-identical legacy behavior).
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        CHECK_FALSE(loader.config().rot3d_enabled);
    }
    // Non-bool value -> loud error naming the line.
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "rot3d_enabled = banana\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        CHECK(loader.parse(text) == Status::InvalidConfig);
        CHECK(loader.error().find("rot3d_enabled") != std::string::npos);
    }
}

TEST_CASE("ConfigLoader: joint_lever_scale parses (Slice 17b); default off; bad value "
          "rejects") {
    // Enabled.
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "joint_lever_scale = true\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        CHECK(loader.config().joint_lever_scale);
    }
    // Key absent -> default OFF (byte-identical 3-unknown lever numerics).
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        CHECK_FALSE(loader.config().joint_lever_scale);
    }
    // Non-bool value -> loud error naming the line.
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "joint_lever_scale = banana\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        CHECK(loader.parse(text) == Status::InvalidConfig);
        CHECK(loader.error().find("joint_lever_scale") != std::string::npos);
    }
}

TEST_CASE("ConfigLoader: subbin_centroid round-trips into all five calibration histograms") {
    // Slice 16: one [global] key switches the centroid sub-bin readout on for
    // every calibration histogram (so3/roll/xyz/scale/offset).
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "subbin_centroid = true\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        const Config& c = loader.config();
        CHECK(c.so3_hist.subbin_centroid);
        CHECK(c.roll_hist.subbin_centroid);
        CHECK(c.xyz_hist.subbin_centroid);
        CHECK(c.scale_hist.subbin_centroid);
        CHECK(c.offset_hist.subbin_centroid);
    }
    // Explicit false parses too (and matches the default).
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "subbin_centroid = false\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        const Config& c = loader.config();
        CHECK_FALSE(c.so3_hist.subbin_centroid);
        CHECK_FALSE(c.offset_hist.subbin_centroid);
    }
    // Key absent -> default OFF everywhere (bit-identical legacy behavior).
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        const Config& c = loader.config();
        CHECK_FALSE(c.so3_hist.subbin_centroid);
        CHECK_FALSE(c.roll_hist.subbin_centroid);
        CHECK_FALSE(c.xyz_hist.subbin_centroid);
        CHECK_FALSE(c.scale_hist.subbin_centroid);
        CHECK_FALSE(c.offset_hist.subbin_centroid);
    }
    // Malformed value -> clear error.
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nsubbin_centroid = maybe\n") == Status::InvalidConfig);
        CHECK_FALSE(loader.error().empty());
    }
}

TEST_CASE("ConfigLoader free function: copies sensor storage to the caller-owned vector") {
    const std::string text = "[sensor.0]\nid = 0\nis_reference = true\n";
    Config cfg;
    std::vector<SensorConfig> sensors;
    std::string err;
    REQUIRE(load_config(text, cfg, sensors, &err) == Status::Ok);
    // Config::sensors points at the CALLER's vector (the documented ownership for the free fn).
    CHECK(cfg.sensors == sensors.data());
    CHECK(cfg.sensor_count == 1);
    CHECK(sensors[0].id == 0);
}

TEST_CASE("ConfigLoader bad input: unknown key / section / malformed value are rejected") {
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nnonsense_key = 3\n") == Status::InvalidConfig);
        CHECK_FALSE(loader.error().empty());
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[unknown_section]\n") == Status::InvalidConfig);
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nwindow_s = not_a_number\n") == Status::InvalidConfig);
    }
    {
        ConfigLoader loader;   // missing '='
        CHECK(loader.parse("[global]\nwindow_s 0.1\n") == Status::InvalidConfig);
    }
    {
        ConfigLoader loader;   // sensor section that skips a slot
        CHECK(loader.parse("[sensor.0]\nid=0\n[sensor.2]\nid=2\n") == Status::InvalidConfig);
    }
    {
        ConfigLoader loader;   // prior_extrinsic with too few numbers
        CHECK(loader.parse("[sensor.0]\nprior_extrinsic = 0.1 0.0\n") == Status::InvalidConfig);
    }
}

TEST_CASE("ConfigLoader: a duplicate sensor id across sections is rejected (InvalidConfig)") {
    // Two distinct [sensor.N] storage slots that bind the SAME id -> silent mis-bind; reject it.
    // (The core validate() does not catch this — its per-sensor checks are a TODO.)
    ConfigLoader loader;
    const Status st = loader.parse("[sensor.0]\nid = 0\nis_reference = true\n[sensor.1]\nid = 0\n");
    CHECK(st == Status::InvalidConfig);
    CHECK(loader.error().find("duplicate sensor id") != std::string::npos);
}

TEST_CASE("ConfigLoader: is_reference disagreeing with reference_sensor_id is rejected") {
    {
        // is_reference flagged on sensor id 1 but the global reference_sensor_id is 0 -> disagree.
        ConfigLoader loader;
        const Status st = loader.parse(
            "[global]\nmax_sources = 2\nreference_sensor_id = 0\n"
            "[sensor.0]\nid = 0\n"
            "[sensor.1]\nid = 1\nis_reference = true\n");
        CHECK(st == Status::InvalidConfig);
        CHECK(loader.error().find("disagrees with reference_sensor_id") != std::string::npos);
    }
    {
        // More than one sensor flagged is_reference -> rejected.
        ConfigLoader loader;
        const Status st = loader.parse(
            "[sensor.0]\nid = 0\nis_reference = true\n"
            "[sensor.1]\nid = 1\nis_reference = true\n");
        CHECK(st == Status::InvalidConfig);
        CHECK(loader.error().find("more than one sensor") != std::string::npos);
    }
    {
        // Agreement: is_reference on the sensor whose id == reference_sensor_id -> Ok.
        ConfigLoader loader;
        const Status st = loader.parse(
            "[global]\nmax_sources = 2\nreference_sensor_id = 1\n"
            "[sensor.0]\nid = 0\n"
            "[sensor.1]\nid = 1\nis_reference = true\n");
        CHECK(st == Status::Ok);
    }
}

TEST_CASE("ConfigLoader: an out-of-range int is reported distinctly from junk") {
    {
        ConfigLoader loader;   // a value far beyond `long` range -> distinct overflow message
        CHECK(loader.parse("[global]\nmax_sources = 99999999999999999999999\n")
              == Status::InvalidConfig);
        CHECK(loader.error().find("out of range") != std::string::npos);
    }
    {
        ConfigLoader loader;   // plain junk -> the generic "expected int"
        CHECK(loader.parse("[global]\nmax_sources = abc\n") == Status::InvalidConfig);
        CHECK(loader.error().find("expected int") != std::string::npos);
    }
}

TEST_CASE("ConfigLoader: a parse that succeeds but fails validate() surfaces the validate status") {
    // max_sources out of range (validate requires 1..32). Parsing is fine; validate rejects it.
    ConfigLoader loader;
    const Status st = loader.parse("[global]\nmax_sources = 99\n");
    CHECK(st == Status::OutOfRange);                 // validate()'s status, not InvalidConfig
    CHECK_FALSE(loader.error().empty());
}

TEST_CASE("ConfigLoader: an empty config validates Ok (all core defaults)") {
    ConfigLoader loader;
    CHECK(loader.parse("# only comments\n\n") == Status::Ok);
    CHECK(loader.config().sensor_count == 0);
    CHECK(loader.config().sensors == nullptr);
}

TEST_CASE("ConfigLoader: multi_bias_enabled + per-sensor bias_states/bias_process_noise parse "
          "(Slice 18); default off; the multi-bias guard surfaces through validate()") {
    // Default off: an empty config leaves the flag false and the per-sensor bias fields at
    // their SensorConfig defaults.
    {
        ConfigLoader loader;
        REQUIRE(loader.parse("# defaults\n") == Status::Ok);
        CHECK_FALSE(loader.config().multi_bias_enabled);
    }

    // The flag + per-sensor keys parse and land in the right fields; TWO bias_states sensors
    // validate Ok only because the flag lifts the Option-A single-bias guard.
    const std::string text =
        "[global]\n"
        "max_sources = 2\n"
        "reference_sensor_id = 0\n"
        "multi_bias_enabled = true\n"
        "[sensor.0]\n"
        "id = 0\n"
        "is_reference = true\n"
        "bias_states = true\n"
        "bias_process_noise = 0.002\n"
        "[sensor.1]\n"
        "id = 1\n"
        "bias_states = on\n";
    {
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        const Config& c = loader.config();
        CHECK(c.multi_bias_enabled);
        REQUIRE(c.sensor_count == 2);
        CHECK(c.sensors[0].bias_states);
        CHECK(c.sensors[0].bias_process_noise == doctest::Approx(0.002));
        CHECK(c.sensors[1].bias_states);
        CHECK(c.sensors[1].bias_process_noise == doctest::Approx(1e-4));  // the core default
    }

    // WITHOUT the flag, the same two bias_states sensors are rejected by the core validate()
    // (the Option-A InvalidConfig guard is intact when the flag is off).
    {
        std::string off_text = text;
        const std::string flag_line = "multi_bias_enabled = true\n";
        off_text.replace(off_text.find(flag_line), flag_line.size(), "");
        ConfigLoader loader;
        CHECK(loader.parse(off_text) == Status::InvalidConfig);
    }

    // Bad values reject with a line diagnostic.
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nmulti_bias_enabled = maybe\n") == Status::InvalidConfig);
        CHECK(loader.error().find("expected bool") != std::string::npos);
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[sensor.0]\nid = 0\nbias_states = maybe\n") == Status::InvalidConfig);
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[sensor.0]\nid = 0\nbias_process_noise = junk\n")
              == Status::InvalidConfig);
    }
    // A negative bias_process_noise parses but fails the core validate() (OutOfRange).
    {
        ConfigLoader loader;
        CHECK(loader.parse("[sensor.0]\nid = 0\nbias_process_noise = -1\n")
              == Status::OutOfRange);
    }
}

TEST_CASE("ConfigLoader: per-DOF bias_process_noise (1 or 6 numbers) + multi_bias_cov0 "
          "(Slice 18 review/B2)") {
    // SIX numbers -> per-DOF rates ([v; omega] order); needs multi_bias_enabled when the
    // vector is non-uniform on a bias source (Option A is uniform-only).
    const std::string per_dof =
        "[global]\n"
        "max_sources = 1\n"
        "reference_sensor_id = 0\n"
        "multi_bias_enabled = true\n"
        "multi_bias_cov0 = 0.09\n"
        "[sensor.0]\n"
        "id = 0\n"
        "is_reference = true\n"
        "bias_states = true\n"
        "bias_process_noise = 0 0 0 0 0 1e-5\n";   // pin all but the yaw-rate DOF
    {
        ConfigLoader loader;
        REQUIRE(loader.parse(per_dof) == Status::Ok);
        const Config& c = loader.config();
        CHECK(c.multi_bias_cov0 == doctest::Approx(0.09));
        REQUIRE(c.sensor_count == 1);
        for (int d = 0; d < 5; ++d) CHECK(c.sensors[0].bias_process_noise[d] == 0.0);
        CHECK(c.sensors[0].bias_process_noise[5] == doctest::Approx(1e-5));
        CHECK_FALSE(c.sensors[0].bias_process_noise.uniform());
    }

    // ONE number still means uniform (the legacy form).
    {
        ConfigLoader loader;
        REQUIRE(loader.parse("[sensor.0]\nid = 0\nbias_process_noise = 2e-5\n")
                == Status::Ok);
        for (int d = 0; d < 6; ++d) {
            CHECK(loader.config().sensors[0].bias_process_noise[d]
                  == doctest::Approx(2e-5));
        }
        CHECK(loader.config().sensors[0].bias_process_noise.uniform());
    }

    // Wrong arity rejects with a line diagnostic.
    {
        ConfigLoader loader;
        CHECK(loader.parse("[sensor.0]\nid = 0\nbias_process_noise = 1e-5 1e-5\n")
              == Status::InvalidConfig);
        CHECK(loader.error().find("expected 1 or 6 numbers") != std::string::npos);
    }

    // A NON-UNIFORM vector on an Option-A (multi_bias_enabled off) bias source parses but
    // fails the core validate() — the per-DOF form (incl. pinning) is an Option-B feature.
    {
        std::string off_text = per_dof;
        const std::string flag_line = "multi_bias_enabled = true\n";
        off_text.replace(off_text.find(flag_line), flag_line.size(), "");
        ConfigLoader loader;
        CHECK(loader.parse(off_text) == Status::InvalidConfig);
    }

    // multi_bias_cov0 must be positive (validate(): OutOfRange) and numeric (parse).
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nmulti_bias_cov0 = 0\n") == Status::OutOfRange);
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nmulti_bias_cov0 = junk\n") == Status::InvalidConfig);
        CHECK(loader.error().find("expected number") != std::string::npos);
    }
}

TEST_CASE("ConfigLoader: split_median / split_veto / rot_weight_prior parse (Slice 19); "
          "defaults preserved; bad values reject") {
    // All three set: the split path enabled, veto disabled, a 10x FOG rotation prior.
    {
        const std::string text =
            "[global]\nmax_sources=2\nreference_sensor_id=0\n"
            "split_median = true\n"
            "split_veto = false\n"
            "[sensor.0]\nid=0\nis_reference=true\n"
            "[sensor.1]\nid=1\nrot_weight_prior = 10.0\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        CHECK(loader.config().split_median);
        CHECK_FALSE(loader.config().split_veto);
        CHECK(loader.sensors()[0].rot_weight_prior == doctest::Approx(1.0));   // default
        CHECK(loader.sensors()[1].rot_weight_prior == doctest::Approx(10.0));
    }
    // Keys absent -> defaults (split off = byte-identical coupled; veto on; prior 1).
    {
        const std::string text =
            "[global]\nmax_sources=1\nreference_sensor_id=0\n"
            "[sensor.0]\nid=0\nis_reference=true\n";
        ConfigLoader loader;
        REQUIRE(loader.parse(text) == Status::Ok);
        CHECK_FALSE(loader.config().split_median);
        CHECK(loader.config().split_veto);
        CHECK(loader.sensors()[0].rot_weight_prior == doctest::Approx(1.0));
    }
    // Bad values: non-bool / non-number reject loudly; a NEGATIVE rot_weight_prior parses
    // but fails the core validate() (OutOfRange) — surfaced as the validate status.
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nsplit_median = banana\n") == Status::InvalidConfig);
        CHECK(loader.error().find("split_median") != std::string::npos);
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[global]\nsplit_veto = banana\n") == Status::InvalidConfig);
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[sensor.0]\nid=0\nrot_weight_prior = junk\n")
              == Status::InvalidConfig);
    }
    {
        ConfigLoader loader;
        CHECK(loader.parse("[sensor.0]\nid=0\nrot_weight_prior = -1.0\n")
              == Status::OutOfRange);
    }
}
