// Slice 10 wiring tests: the per-sensor fixed-lag RTS smoother integrated into the
// estimator's CALIBRATION feed (the deeper, two-sided-smoothed frontier — D18).
//
// Two properties, both end-to-end through the Estimator/Rig (the sim-GT oracle):
//   * NON-REGRESSION (hard guard): with per_sensor_kf OFF for every source (the default),
//     the smoother branch never activates — the calibration feed is byte-identical to the
//     pre-Slice-10 path, AND the FUSION frontier is byte-identical OFF vs ON (fusion is
//     untouched by the smoother). Turning the smoother ON for a source DOES change the
//     calibration snapshot (proving the deeper feed is actually wired, not a no-op).
//   * PEAK SHARPENING (the D18 headline, carried through to calibration): a NOISY source on
//     a straight+turn trajectory calibrates with HIGHER histogram concentration
//     (extrinsic_confidence / scale_confidence) when the smoother is ON, with NO estimate
//     bias (the zero-phase property carries through to no calibration bias).
#include <doctest/doctest.h>

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"

#include "ofc_sim/rig.hpp"
#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;
bool near_abs(Scalar a, Scalar b, Scalar tol) { return std::abs(a - b) <= tol; }

// A mixed bootstrap trajectory (mirrors test_calib_feedback): repeated STRAIGHT stretches
// (yaw/pitch + scale observable) + multi-axis TURNS (roll + xyz) so calibration accrues
// many votes.
Trajectory mixed_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0, 0.35,  0.6;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.35, 0.6;
    for (int rep = 0; rep < 4; ++rep) {
        tr.add_segment(straight, 2.0);
        tr.add_segment(turnA,    1.6);
        tr.add_segment(straight, 1.6);
        tr.add_segment(turnB,    1.6);
    }
    return tr;
}

// Histogram knobs (vote_weight = One so totals are vote COUNTS), matching the feedback test.
void set_hists(Config& c) {
    c.vote_weight = VoteWeight::One;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.turn_omega_min     = 0.20;

    c.so3_hist.bins = 256; c.so3_hist.range_min = -0.8; c.so3_hist.range_max = 0.8;
    c.so3_hist.circular = false; c.so3_hist.aging = Aging::SlidingK; c.so3_hist.sliding_k = 300;
    c.so3_hist.vote_split = true; c.so3_hist.subbin = true;

    c.scale_hist.bins = 256; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.circular = false; c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 300;
    c.scale_hist.vote_split = true; c.scale_hist.subbin = true;

    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 300;
    c.roll_hist.vote_split = true; c.roll_hist.subbin = true;

    c.xyz_hist.bins = 256; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.circular = false; c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 300;
    c.xyz_hist.vote_split = true; c.xyz_hist.subbin = true;
}

const CalibSnapshot* snap(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.calib[i].id == id) return &r.calib[i];
    return nullptr;
}

// The forward axis a source with planted extrinsic X reports, in the BASE frame: X.R*e_x.
Vec3 planted_forward(const SE3& X) { return X.R * Vec3(1, 0, 0); }
} // namespace

// ===========================================================================
// NON-REGRESSION: per_sensor_kf OFF == byte-identical calibration; fusion untouched ON/OFF
// ===========================================================================
TEST_CASE("smoother wiring: per_sensor_kf OFF is byte-identical; ON leaves fusion identical "
          "but changes calibration") {
    const Trajectory tr = mixed_traj();
    SE3 X1; X1.R = so3::exp(Vec3(0.0, 0.0, kPi / 6)); X1.t = Vec3(0.3, -0.2, 0.1);

    // 3 sources; source 1 has a real mount + a little noise so the smoother has something to
    // smooth. Priors == planted (no feedback needed — we are testing the FEED, not commits).
    auto build = [&](bool smoother_on,
                     std::vector<std::unique_ptr<SyntheticSource>>& srcs,
                     std::vector<SensorConfig>& sensors) {
        std::vector<SourceParams> planted(3);
        for (int i = 0; i < 3; ++i) planted[i].id = static_cast<SourceId>(i);
        planted[1].X = X1;
        planted[1].noise_trans_per_m = 0.02; planted[1].noise_rot_per_rad = 0.02;
        planted[1].noise_trans_floor = 0.005; planted[1].noise_rot_floor = 0.005;
        planted[1].seed = 4242u;
        srcs.clear();
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

        sensors.assign(3, SensorConfig{});
        for (int i = 0; i < 3; ++i) { sensors[i].id = static_cast<SourceId>(i); }
        sensors[1].prior_extrinsic = X1;
        sensors[1].per_sensor_kf   = smoother_on;     // the toggle under test
        sensors[1].kf_process_noise = 5.0;
    };

    auto run = [&](bool smoother_on) {
        auto srcs    = std::make_shared<std::vector<std::unique_ptr<SyntheticSource>>>();
        auto sensors = std::make_shared<std::vector<SensorConfig>>();
        build(smoother_on, *srcs, *sensors);
        Config cfg;
        cfg.max_sources    = 3;
        cfg.fusion_delay_s = 0.05;
        cfg.window_s       = 0.10;
        cfg.calib_lag_s    = 0.20;        // L = 0.20 * 50 = 10 steps
        cfg.tick_rate_hz   = 50.0;
        cfg.timesync_enabled = false;
        cfg.cold_start     = ColdStart::MedianFromStart;
        set_hists(cfg);
        // Disable the calibration->fusion FEEDBACK (Slice 8) for this test: an unreachable
        // commit gate means NO calibrated DOF ever swaps into the fusion prior. That ISOLATES
        // the smoother's effect to the calibration FEED — so the FUSION frontier must be
        // byte-identical ON vs OFF (the smoother never touches fusion directly). With feedback
        // ENABLED the deeper-frontier votes legitimately change the commits and hence fusion
        // (the intended D18 latency-tolerant behavior — exercised in the peak-sharpening test).
        cfg.commit_min_votes = 1000000000;
        cfg.sensors        = sensors->data();
        cfg.sensor_count   = 3;

        auto rig = std::make_shared<Rig>();
        rig->set_trajectory(tr);
        REQUIRE(rig->init(cfg) == Status::Ok);
        for (auto& sp : *srcs) REQUIRE(rig->add_source(sp.get()) == Status::Ok);
        const int fuses = rig->run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 200);
        // Keep the storage alive for as long as the rig is queried.
        return std::make_tuple(rig, srcs, sensors);
    };

    auto off_a = run(false);
    auto off_b = run(false);
    auto on    = run(true);
    const std::vector<Record>& roff_a = std::get<0>(off_a)->records();
    const std::vector<Record>& roff_b = std::get<0>(off_b)->records();
    const std::vector<Record>& ron    = std::get<0>(on)->records();
    REQUIRE(roff_a.size() == roff_b.size());
    REQUIRE(roff_a.size() == ron.size());

    // (1) OFF is DETERMINISTIC + byte-identical run-to-run (the deeper path inert).
    bool off_identical = true;
    for (std::size_t i = 0; i < roff_a.size(); ++i) {
        if (roff_a[i].fused != roff_b[i].fused) { off_identical = false; break; }
        if (!roff_a[i].fused) continue;
        const State& fa = roff_a[i].result.frontier;
        const State& fb = roff_b[i].result.frontier;
        if (!(fa.pose.t.array() == fb.pose.t.array()).all()) off_identical = false;
        if (!(fa.cov.array() == fb.cov.array()).all())       off_identical = false;
        // Calibration snapshot byte-identical too.
        const CalibSnapshot* ca = snap(roff_a[i].result, 1);
        const CalibSnapshot* cb = snap(roff_b[i].result, 1);
        if (ca && cb) {
            if (!(ca->extrinsic.R.array() == cb->extrinsic.R.array()).all()) off_identical = false;
            if (ca->extrinsic_confidence != cb->extrinsic_confidence)        off_identical = false;
        }
        if (!off_identical) break;
    }
    CHECK(off_identical);

    // (2) FUSION frontier is byte-identical OFF vs ON — the smoother NEVER touches fusion.
    bool fusion_untouched = true;
    for (std::size_t i = 0; i < ron.size(); ++i) {
        if (roff_a[i].fused != ron[i].fused) { fusion_untouched = false; break; }
        if (!ron[i].fused) continue;
        const State& fo = roff_a[i].result.frontier;
        const State& fn = ron[i].result.frontier;
        if (!(fo.pose.R.array() == fn.pose.R.array()).all()) fusion_untouched = false;
        if (!(fo.pose.t.array() == fn.pose.t.array()).all()) fusion_untouched = false;
        if (!(fo.twist.xi.array() == fn.twist.xi.array()).all()) fusion_untouched = false;
        if (!(fo.cov.array() == fn.cov.array()).all())       fusion_untouched = false;
        if (!fusion_untouched) break;
    }
    CHECK(fusion_untouched);

    // (3) CALIBRATION did change with the smoother ON (proves the deeper feed is wired). The
    // smoothed source's recovered extrinsic differs from the OFF run's at the END of the run.
    const CalibSnapshot* cs_off = snap(std::get<0>(off_a)->estimator().latest(), 1);
    const CalibSnapshot* cs_on  = snap(std::get<0>(on)->estimator().latest(), 1);
    REQUIRE(cs_off != nullptr);
    REQUIRE(cs_on  != nullptr);
    const Scalar dR = (cs_off->extrinsic.R - cs_on->extrinsic.R).cwiseAbs().maxCoeff();
    const Scalar dconf = std::abs(cs_off->extrinsic_confidence - cs_on->extrinsic_confidence);
    CHECK((dR > 1e-9 || dconf > 1e-9));     // the deeper feed materially altered calibration
}

// ===========================================================================
// PEAK SHARPENING: a noisy source calibrates sharper with the smoother ON, no bias (D18)
// ===========================================================================
TEST_CASE("smoother wiring: a noisy source calibrates with SHARPER histogram peaks ON vs "
          "OFF, with no estimate bias (zero-phase carries to no calibration bias)") {
    const Trajectory tr = mixed_traj();
    // Source 1: a real yaw/pitch mount + scale + substantial noise (the regime the smoother
    // helps). Priors == planted (we measure CONCENTRATION + bias of the residual histogram,
    // not bootstrap convergence).
    const Scalar yaw_t = 0.20, pitch_t = -0.12, scale_t = 1.15;
    SE3 X1;
    {
        Mat3 Rz; Rz << std::cos(yaw_t), -std::sin(yaw_t), 0,
                       std::sin(yaw_t),  std::cos(yaw_t), 0, 0, 0, 1;
        Mat3 Ry; Ry << std::cos(pitch_t), 0, std::sin(pitch_t),
                       0, 1, 0, -std::sin(pitch_t), 0, std::cos(pitch_t);
        X1.R = Rz * Ry; X1.t = Vec3::Zero();
    }

    auto run = [&](bool smoother_on,
                   std::vector<std::unique_ptr<SyntheticSource>>& srcs,
                   std::vector<SensorConfig>& sensors) {
        std::vector<SourceParams> planted(3);
        for (int i = 0; i < 3; ++i) planted[i].id = static_cast<SourceId>(i);
        planted[1].X = X1; planted[1].scale = scale_t;
        // Heavy noise on the calibrated source — this is what a two-sided smoother cleans up.
        planted[1].noise_trans_per_m = 0.05; planted[1].noise_rot_per_rad = 0.05;
        planted[1].noise_trans_floor = 0.01; planted[1].noise_rot_floor = 0.01;
        planted[1].seed = 9001u;
        srcs.clear();
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

        sensors.assign(3, SensorConfig{});
        for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
        sensors[1].prior_extrinsic = X1;
        sensors[1].prior_scale     = scale_t;
        sensors[1].per_sensor_kf   = smoother_on;
        sensors[1].kf_process_noise = 5.0;

        Config cfg;
        cfg.max_sources    = 3;
        cfg.fusion_delay_s = 0.05;
        cfg.window_s       = 0.10;
        cfg.calib_lag_s    = 0.20;
        cfg.tick_rate_hz   = 50.0;
        cfg.timesync_enabled = false;
        cfg.cold_start     = ColdStart::MedianFromStart;
        set_hists(cfg);
        cfg.sensors        = sensors.data();
        cfg.sensor_count   = 3;

        auto rig = std::unique_ptr<Rig>(new Rig());
        rig->set_trajectory(tr);
        REQUIRE(rig->init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(rig->add_source(sp.get()) == Status::Ok);
        const int fuses = rig->run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 200);
        return rig;
    };

    std::vector<std::unique_ptr<SyntheticSource>> srcs_off, srcs_on;
    std::vector<SensorConfig> sens_off, sens_on;
    auto rig_off = run(false, srcs_off, sens_off);
    auto rig_on  = run(true,  srcs_on,  sens_on);

    const CalibSnapshot* cs_off = snap(rig_off->estimator().latest(), 1);
    const CalibSnapshot* cs_on  = snap(rig_on->estimator().latest(),  1);
    REQUIRE(cs_off != nullptr);
    REQUIRE(cs_on  != nullptr);

    // --- PEAK SHARPENING: concentration is HIGHER with the smoother ON ----------------
    INFO("extrinsic_confidence OFF=" << cs_off->extrinsic_confidence
         << "  ON=" << cs_on->extrinsic_confidence);
    INFO("scale_confidence     OFF=" << cs_off->scale_confidence
         << "  ON=" << cs_on->scale_confidence);
    CHECK(cs_on->extrinsic_confidence > cs_off->extrinsic_confidence);
    CHECK(cs_on->scale_confidence     > cs_off->scale_confidence);

    // --- NO ESTIMATE BIAS: the recovered forward axis + scale stay on truth (both runs) ---
    // The two-sided RTS is zero-phase, so sharpening the peak must NOT shift it off truth.
    const Vec3 f_true = planted_forward(X1);
    const Vec3 f_off  = cs_off->extrinsic.R * Vec3(1, 0, 0);
    const Vec3 f_on   = cs_on->extrinsic.R  * Vec3(1, 0, 0);
    INFO("fwd-axis err OFF=" << (f_off - f_true).norm()
         << "  ON=" << (f_on - f_true).norm());
    // ON is well within tolerance of truth (no bias) AND no worse than OFF (sharper, not
    // shifted). A small absolute tolerance covers the histogram bin/segment-dynamics floor.
    CHECK((f_on - f_true).norm() < 0.05);
    CHECK((f_on - f_true).norm() <= (f_off - f_true).norm() + 0.02);
    CHECK(near_abs(cs_on->scale, scale_t, 3e-2));     // scale unbiased with the smoother ON
}
