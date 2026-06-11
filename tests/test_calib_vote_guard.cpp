// test_calib_vote_guard.cpp — Slice 18 review/B1: the Slice-17b out-of-range vote SKIP
// guard (scale2 precedent: skip, never edge-clamp, count skips) GENERALIZED to the other
// calibration vote channels — Phase-1 so(3) direction (3 channels, whole-vote), Phase-2
// xyz lever (per channel), Phase-2 rot3d (3 channels, whole-vote).
//
// WHY (the urban12 real-data corruption channel): weakly-observable bias DOFs wander ->
// the de-biased deltas carry a growing fake rotation -> calibration votes fall OUTSIDE
// the configured histogram ranges -> pre-B1 they edge-CLAMPED into the boundary bin ->
// deterministic boundary mass reads concentrated -> COMMITS at the edge-bin center
// (±0.984375 at 64 bins over [-1,1] — the exact committed values measured on urban07/12)
// -> feedback applies absurd extrinsics -> fusion destroyed. The guard makes boundary-bin
// mass UNABLE to commit: an out-of-range vote deposits nothing and increments a per-slot
// per-family skip counter (so3_skipped / xyz_skipped / rot3d_skipped, mirroring
// scale2_skipped) — fail-safe (calib stays at the prior, conf 0) and diagnosable
// ("votes 0, skipped HIGH" = out-of-regime; "votes 0, skipped 0" = unexcited).
//
// Coverage:
//   * Phase-1 so(3): an out-of-range direction vote (57°..90° off +e_x — past the range,
//     short of the π-guard) deposits NOTHING + counter grows; an in-range control deposits
//     with counter 0 (in-range behavior byte-identical).
//   * Phase-2 xyz: an out-of-range lever solve channel deposits nothing + counter grows;
//     in-range channels still vote (per-channel guard); conf gated by the empty channel.
//   * Phase-2 rot3d: an out-of-range δφ residual deposits nothing + counter grows.
//   * ESTIMATOR-LEVEL urban12 signature: a large fake rotation on a source leaves the
//     calibration UNCOMMITTED at confidence 0 — NOT committed at the ±0.984 edge-bin
//     center. Removing any guard (the clamp mutation) re-commits boundary mass and the
//     corresponding case fails loudly.
#include <doctest/doctest.h>

#include "ofc/core/calibration.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"

#include "ofc_sim/rig.hpp"
#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <cmath>
#include <memory>
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

Mat3 yaw_rot(Scalar yaw) {
    Mat3 R;
    R << std::cos(yaw), -std::sin(yaw), 0,
         std::sin(yaw),  std::cos(yaw), 0,
         0,              0,             1;
    return R;
}

std::unique_ptr<Phase1Calibrator> make_calib1() {
    return std::unique_ptr<Phase1Calibrator>(new Phase1Calibrator());
}
std::unique_ptr<Phase2Calibrator> make_calib2() {
    return std::unique_ptr<Phase2Calibrator>(new Phase2Calibrator());
}

// Drive Phase-1 directly: a reference reporting straight-forward plus one source whose
// reported direction sits `theta` rad off the base forward axis (through an identity
// prior, g_obs lands exactly theta off +e_x). Returns after `steps` observed steps.
void drive_phase1(Phase1Calibrator& cal, Scalar theta, int steps) {
    SourceId ids[2] = {0, 1};
    SE3 rep[2];
    rep[0].t = Vec3(0.2, 0, 0);                       // reference: straight forward
    rep[1].t = yaw_rot(theta).transpose() * Vec3(0.2, 0, 0);   // dir_B = R^-1 e_x * 0.2
    const Vec3 omega = Vec3::Zero();
    const Vec3 trans(0.2, 0, 0);
    for (int i = 0; i < steps; ++i) {
        cal.observe(2, ids, rep, omega, trans, nullptr);
    }
}
} // namespace

// ===========================================================================
// Phase-1 so(3): skip not clamp + counter; in-range byte-identical.
// ===========================================================================
TEST_CASE("vote guard / Phase-1 so3: an out-of-range direction vote deposits NOTHING and "
          "is counted; an in-range vote deposits with the counter at 0") {
    Config cfg;                                      // default so3_hist range [-1, 1]
    cfg.straight_omega_max = 0.05;
    cfg.straight_trans_min = 0.02;

    // OUT of range: 1.3 rad (~74°) off +e_x — beyond the [-1, 1] so(3) channel range but
    // short of the 90° π-singularity guard (which would have skipped it silently pre-B1;
    // 57°..90° was the edge-CLAMP corruption window).
    {
        auto calp = make_calib1(); Phase1Calibrator& cal = *calp;
        REQUIRE(cal.configure(cfg, 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        drive_phase1(cal, Scalar(1.3), 300);
        CHECK(cal.vote_count(1) == doctest::Approx(0.0));        // nothing deposited
        CHECK(cal.extrinsic_confidence(1) == doctest::Approx(0.0));
        CHECK(cal.so3_skipped(1) == doctest::Approx(300.0));     // every vote counted
        // The estimate stays at the prior — not a boundary-bin ~0.984 commit candidate.
        const Vec3 f = cal.forward_axis(1);
        CHECK((f - Vec3(1, 0, 0)).norm() < 1e-12);
        // The reference's own (in-range, straight-forward) votes were unaffected.
        CHECK(cal.vote_count(0) > Scalar(0));
        CHECK(cal.so3_skipped(0) == doctest::Approx(0.0));
    }

    // IN range: 0.5 rad — deposits exactly as pre-B1 (in-range behavior unchanged).
    {
        auto calp = make_calib1(); Phase1Calibrator& cal = *calp;
        REQUIRE(cal.configure(cfg, 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        drive_phase1(cal, Scalar(0.5), 300);
        CHECK(cal.vote_count(1) > Scalar(0));
        CHECK(cal.so3_skipped(1) == doctest::Approx(0.0));
        const Vec3 f = cal.forward_axis(1);
        const Vec3 expect = yaw_rot(Scalar(0.5)).transpose() * Vec3(1, 0, 0);
        CHECK((f - expect).norm() < 0.05);
    }
}

// ===========================================================================
// Phase-2 xyz: per-channel skip not clamp + counter (legacy 3-unknown path).
// ===========================================================================
TEST_CASE("vote guard / Phase-2 xyz: an out-of-range lever channel deposits NOTHING and is "
          "counted; in-range channels still vote; conf stays 0 via the empty channel") {
    // A planted lever with x = 1.5 m — beyond the default xyz_hist [-1, 1] — y/z in range.
    const SE3 X = [] {
        SE3 x; x.R = Mat3::Identity(); x.t = Vec3(1.5, 0.2, -0.1); return x;
    }();
    Config cfg;                                      // default xyz_hist range [-1, 1]
    cfg.turn_omega_min = 0.20;

    auto calp = make_calib2(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(cfg, 0) == Status::Ok);
    SE3 prior;                                       // identity prior (zero lever)
    REQUIRE(cal.set_prior(1, prior) == Status::Ok);
    REQUIRE(cal.set_yaw_pitch(1, Mat3::Identity()) == Status::Ok);

    // Multi-axis turning conjugated stream: B = X^-1 A X (the exact hand-eye), enough
    // excitation that all 3 lever axes are observable — only the RANGE withholds x.
    const Scalar dt = 0.1;
    int k = 0;
    for (int i = 0; i < 300; ++i) {
        Vec6 tw;
        const Scalar ph = Scalar(0.7) * static_cast<Scalar>(k++);
        tw << 1.0, 0.2 * std::sin(ph), 0.1 * std::cos(ph),
              0.4 * std::sin(ph), 0.4 * std::cos(ph), 0.5;
        const SE3 A = se3::exp(tw * dt);
        const SE3 B = se3::compose(se3::compose(se3::inverse(X), A), X);
        SourceId ids[1] = {1};
        SE3 rep[1] = {B};
        const Vec3 omega = tw.tail<3>();
        cal.observe(1, ids, rep, A, omega, nullptr);
    }

    // The x channel (true value 1.5, out of range) never voted -> empty -> conf 0; the
    // y/z channels (in range) did. The committed-readout fallback for x is the PRIOR.
    CHECK(cal.xyz_skipped(1) > Scalar(100));
    CHECK(cal.translation_confidence(1) == doctest::Approx(0.0));
    const Vec3 lever = cal.lever_arm(1);
    CHECK(lever.x() == doctest::Approx(0.0));                  // prior, not ~0.984
    CHECK(std::abs(lever.y() - X.t.y()) < 0.05);               // in-range channels intact
    CHECK(std::abs(lever.z() - X.t.z()) < 0.05);
}

// ===========================================================================
// Phase-2 rot3d: whole-vote skip not clamp + counter.
// ===========================================================================
TEST_CASE("vote guard / Phase-2 rot3d: an out-of-range rotation residual deposits NOTHING "
          "and is counted; conf 0; readout stays at the basepoint") {
    // Planted rotation 1.3 rad about z (beyond the default so3_hist [-1, 1] residual
    // range vs the identity basepoint), with multi-axis excitation so the BBw two-axis
    // gate is OPEN — only the RANGE withholds the votes.
    SE3 X; X.R = yaw_rot(Scalar(1.3)); X.t = Vec3(0.1, -0.05, 0.2);
    Config cfg;                                      // default so3_hist range [-1, 1]
    cfg.turn_omega_min = 0.20;
    cfg.rot3d_enabled  = true;

    auto calp = make_calib2(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(cfg, 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);  // identity prior == rot3d basepoint

    const Scalar dt = 0.1;
    int k = 0;
    for (int i = 0; i < 300; ++i) {
        Vec6 tw;
        const Scalar ph = Scalar(0.7) * static_cast<Scalar>(k++);
        tw << 1.0, 0.0, 0.0, 0.5 * std::sin(ph), 0.5 * std::cos(ph), 0.4;
        const SE3 A = se3::exp(tw * dt);
        const SE3 B = se3::compose(se3::compose(se3::inverse(X), A), X);
        SourceId ids[1] = {1};
        SE3 rep[1] = {B};
        const Vec3 omega = tw.tail<3>();
        cal.observe(1, ids, rep, A, omega, nullptr);
    }

    CHECK(cal.rot3d_observable(1));                  // the two-axis gate IS open
    CHECK(cal.rot3d_vote_count(1) == doctest::Approx(0.0));    // skip, NOT clamp
    CHECK(cal.rot3d_confidence(1) == doctest::Approx(0.0));
    CHECK(cal.rot3d_skipped(1) > Scalar(100));
    // Readout falls back to the basepoint — not a boundary-bin partial rotation.
    CHECK((cal.rot3d(1) - Mat3::Identity()).norm() < 1e-12);
}

// ===========================================================================
// ESTIMATOR-LEVEL urban12 signature: a large fake rotation on a source must leave the
// calibration UNCOMMITTED (conf 0) — not committed at the [-1,1] edge-bin center ±0.984.
// ===========================================================================
TEST_CASE("vote guard / estimator urban12 signature: a large fake rotation on a source "
          "leaves calib UNCOMMITTED at conf 0 — never committed-at-0.984") {
    // Straight drive; source 1 carries a LARGE rotation offset vs its (identity) prior —
    // the urban12 shape, where the wandering bias DOFs make a source's de-biased deltas
    // read a big fake rotation. Its Phase-1 direction votes are ~1.31 rad off +e_x:
    // outside the default [-1, 1] so(3) range, inside the <90° π-guard — pre-B1 these
    // edge-clamped, concentrated in the boundary bin, and COMMITTED ~0.984.
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    tr.add_segment(straight, 40.0);

    std::vector<SourceParams> planted(2);
    planted[0].id = 0;                               // reference, identity mount
    planted[1].id = 1;
    planted[1].X.R = yaw_rot(Scalar(0.75) * kPi / Scalar(1.8));   // yaw ~75° (1.31 rad)
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(2);
    sensors[0].id = 0; sensors[0].is_reference = true;
    sensors[1].id = 1;                               // identity PRIOR (the fake-rotation gap)

    Config cfg;                                      // default histograms ([-1,1] so3)
    cfg.max_sources = 2; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    // Straight gate vs the tick cadence: at 50 Hz the per-step displacement is 2.0 m/s *
    // 0.02 s = 0.04 m — below the default straight_trans_min (0.05, cadence-dependent per
    // calibration.hpp), which would close the gate after the first wide bootstrap window
    // and starve Phase-1 entirely. Open it for this cadence.
    cfg.straight_trans_min = 0.02;
    cfg.reference_sensor_id = 0;
    cfg.cold_start = ColdStart::ReferenceOnly;       // consensus = the clean reference
    cfg.timesync_enabled = false;
    cfg.min_sources_warn = 1;
    cfg.commit_concentration = 0.5;
    cfg.commit_drop = 0.3;
    cfg.commit_min_votes = 60;                       // commit EASILY — the guard must win
    cfg.sensors = sensors.data(); cfg.sensor_count = 2;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    const Result& r = rig.estimator().latest();
    REQUIRE(r.source_count == 2);
    const CalibSnapshot* cs = nullptr;
    for (int i = 0; i < r.source_count; ++i) {
        if (r.calib[i].id == 1) { cs = &r.calib[i]; break; }
    }
    REQUIRE(cs != nullptr);
    INFO("ext conf=" << cs->extrinsic_confidence
         << "  committed=" << cs->extrinsic_committed
         << "  so3_skipped=" << rig.estimator().so3_skipped(1));
    // UNCOMMITTED at conf 0 — the corruption channel is dead (the clamp mutation commits
    // a confidently-wrong ~0.984 rotation here and every CHECK fails)...
    CHECK_FALSE(cs->extrinsic_committed);
    CHECK(cs->extrinsic_confidence == doctest::Approx(0.0));
    // ...the published extrinsic stays at the (identity) prior, NOT an edge-bin value...
    CHECK(so3::log(cs->extrinsic.R).norm() < 1e-9);
    // ...and the situation is diagnosable, not silent.
    CHECK(rig.estimator().so3_skipped(1) > Scalar(100));
    // The clean reference keeps voting in range.
    CHECK(rig.estimator().so3_skipped(0) == doctest::Approx(0.0));
}
