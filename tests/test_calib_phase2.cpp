// Slice 7 tests: Phase-2 calibration — turn-regime roll (about forward) + xyz lever-arm
// via hand-eye, for BOTH strategies (VsFusedBase, PairwisePinnedRef).
//
// The oracle is the sim: synthetic sources with a PLANTED extrinsic (yaw/pitch + ROLL +
// lever-arm translation X.t) and scale, driven on turning() (and a multi-axis turning
// trajectory so all 3 t_X components are observable) vs straight() (the observability
// self-test). The hand-eye is A·X = X·B with A = base motion, B = X^{-1} A X (the sim's
// measurement model), so the recovery is exact up to noise.
//
// Coverage:
//   * Convergence (BOTH strategies) — planted roll + t_X recovered under turning;
//     confidence rises.
//   * Observability self-test (LOAD-BEARING) — a straight()-only run does NOT converge:
//     roll/xyz stay near prior, low confidence / LS ill-conditioned + skipped.
//   * Lever-arm needs rotation — a higher turn rate converges t_X tighter, and the
//     near-R_A=I (low-omega) rows are rejected.
//   * Yaw-only observability — pure yaw turning recovers t_X.x/.y but leaves t_X.z
//     unobservable (low translation confidence; the LS z-axis is the null space).
//   * Turn gate (unit) — accepts turning, rejects straight.
//   * Determinism.
//   * End-to-end wiring — the full Estimator/Rig over a turning run populates
//     CalibSnapshot.extrinsic (roll + lever arm) + translation_confidence.
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
Timestamp secs(double s) { return static_cast<Timestamp>(std::llround(s * 1e9)); }
bool near_abs(Scalar a, Scalar b, Scalar tol) { return std::abs(a - b) <= tol; }

Mat3 Rz(Scalar a) { Mat3 R; R << std::cos(a), -std::sin(a), 0,
                                 std::sin(a),  std::cos(a), 0,
                                 0, 0, 1; return R; }
Mat3 Ry(Scalar a) { Mat3 R; R << std::cos(a), 0, std::sin(a),
                                 0, 1, 0,
                                -std::sin(a), 0, std::cos(a); return R; }
Mat3 Rx(Scalar a) { Mat3 R; R << 1, 0, 0,
                                 0, std::cos(a), -std::sin(a),
                                 0, std::sin(a),  std::cos(a); return R; }

// Sensor->base extrinsic = Rz(yaw) Ry(pitch) Rx(roll), with lever-arm translation t.
// The yaw/pitch part (roll=0) is R_yp = Rz(yaw) Ry(pitch) — what Phase 1 recovers.
SE3 make_extrinsic(Scalar yaw, Scalar pitch, Scalar roll, const Vec3& t) {
    SE3 X; X.R = Rz(yaw) * Ry(pitch) * Rx(roll); X.t = t;
    return X;
}
Mat3 yaw_pitch_R(Scalar yaw, Scalar pitch) { return Rz(yaw) * Ry(pitch); }

// Heap-allocate the calibrator (large fixed-capacity histograms overflow the test stack).
std::unique_ptr<Phase2Calibrator> make_calib() {
    return std::unique_ptr<Phase2Calibrator>(new Phase2Calibrator());
}

// A Config with sensible Phase-2 histograms. roll = circular S¹ over (−π, π]; xyz = a
// generous metre range. SlidingK so a few votes give a sharp peak.
Config make_cfg() {
    Config c;
    c.tick_rate_hz   = 50.0;
    c.reference_sensor_id = 0;
    c.turn_omega_min = 0.20;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.phase2_strategy = Phase2Strat::VsFusedBase;

    // Phase-1 histograms (so the estimator's Phase-1 path converges in the wiring test).
    c.so3_hist.bins       = 512;
    c.so3_hist.range_min  = -0.8;
    c.so3_hist.range_max  =  0.8;
    c.so3_hist.circular   = false;
    c.so3_hist.aging      = Aging::SlidingK;
    c.so3_hist.sliding_k  = 256;
    c.so3_hist.vote_split = true;
    c.so3_hist.subbin     = true;

    c.scale_hist.bins       = 512;
    c.scale_hist.range_min  = 0.5;
    c.scale_hist.range_max  = 1.5;
    c.scale_hist.circular   = false;
    c.scale_hist.aging      = Aging::SlidingK;
    c.scale_hist.sliding_k  = 256;
    c.scale_hist.vote_split = true;
    c.scale_hist.subbin     = true;

    c.roll_hist.bins       = 360;
    c.roll_hist.range_min  = -kPi;
    c.roll_hist.range_max  =  kPi;
    c.roll_hist.circular   = true;
    c.roll_hist.aging      = Aging::SlidingK;
    c.roll_hist.sliding_k  = 256;
    c.roll_hist.vote_split = true;
    c.roll_hist.subbin     = true;

    c.xyz_hist.bins       = 512;
    c.xyz_hist.range_min  = -1.5;
    c.xyz_hist.range_max  =  1.5;
    c.xyz_hist.circular   = false;
    c.xyz_hist.aging      = Aging::SlidingK;
    c.xyz_hist.sliding_k  = 256;
    c.xyz_hist.vote_split = true;
    c.xyz_hist.subbin     = true;
    return c;
}

// A multi-axis turning trajectory (yaw + pitch rates) so all 3 lever-arm components are
// observable (pure yaw leaves t_X.z in the null space of R_A − I).
Trajectory multiaxis_turn(Scalar v = 2.0, Scalar wy = 0.3, Scalar wz = 0.6,
                          Scalar dur = 6.0) {
    Trajectory tr;
    Vec6 a; a << v, 0, 0, 0, wy,  wz;     // forward + pitch-rate + yaw-rate
    Vec6 b; b << v, 0, 0, 0, -wy, wz;     // flip pitch-rate (more rotation-axis variety)
    tr.add_segment(a, dur * 0.5);
    tr.add_segment(b, dur * 0.5);
    return tr;
}

// Drive the Phase-2 calibrator over a trajectory by querying a reference + one planted
// source at a fixed cadence, feeding observe() the de-scaled reported deltas + the GT
// fused motion/omega (clean oracle drive — no median). Registers R_yp for the planted
// source (the Phase-1 result) so roll is parameterized on top of yaw/pitch. Returns the
// number of voted (turn-gated, non-skipped) steps.
struct DriveResult { int voted = 0; };
DriveResult drive(Phase2Calibrator& cal, const Trajectory& traj,
                  const SyntheticSource& ref, const SyntheticSource& planted,
                  Scalar prior_scale_planted,
                  Scalar from_s, Scalar to_s, Scalar rate_hz) {
    DriveResult dr;
    const Timestamp step = secs(1.0 / rate_hz);
    SourceId ids[2];
    SE3      rep[2];
    for (Timestamp t = secs(from_s); t + step <= secs(to_s); t += step) {
        const Expected<Delta> qr = ref.query(t, t + step);
        const Expected<Delta> qp = planted.query(t, t + step);
        if (!qr.ok() || !qp.ok()) continue;
        // GT fused motion A + body omega over [t, t+step].
        const SE3 A = se3::compose(se3::inverse(traj.pose(t)), traj.pose(t + step));
        const Scalar dt = static_cast<Scalar>(step) / 1e9;
        const Vec3 fused_omega = so3::log(A.R) / dt;
        // De-scale each source's reported translation by its prior_scale.
        SE3 br = qr.value().motion;
        SE3 bp = qp.value().motion; bp.t = bp.t / prior_scale_planted;
        ids[0] = ref.id();     rep[0] = br;
        ids[1] = planted.id(); rep[1] = bp;
        const Status st = cal.observe(2, ids, rep, A, fused_omega);
        if (ok(st)) ++dr.voted;
    }
    return dr;
}
} // namespace

// ===========================================================================
// configure() validation
// ===========================================================================
TEST_CASE("phase2 configure validates and reports state") {
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), /*reference=*/0) == Status::Ok);
    CHECK(cal.configured());
    CHECK(cal.reference() == 0);

    Config bad = make_cfg();
    bad.roll_hist.bins = 2;            // < 4 -> rejected
    auto cal2 = make_calib();
    CHECK(cal2->configure(bad, 0) != Status::Ok);
}

// ===========================================================================
// Turn gate (unit) — accepts turning, rejects straight
// ===========================================================================
TEST_CASE("phase2 turn gate: accepts turning, rejects straight") {
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}) == Status::Ok);

    SourceId ids[1] = {0};
    SE3      rep[1];
    rep[0].R = so3::exp(Vec3(0, 0, 0.01));
    rep[0].t = Vec3(0.04, 0, 0);
    SE3 A; A.R = so3::exp(Vec3(0, 0, 0.01)); A.t = Vec3(0.04, 0, 0);

    // Turning: omega above the gate -> voted (Ok).
    CHECK(cal.observe(1, ids, rep, A, Vec3(0, 0, 0.5)) == Status::Ok);
    CHECK(cal.gated_turning());

    // Straight: omega below the gate -> gated out.
    CHECK(cal.observe(1, ids, rep, A, Vec3(0, 0, 0.0)) == Status::NotReady);
    CHECK_FALSE(cal.gated_turning());
}

// ===========================================================================
// Convergence — VsFusedBase: planted roll + lever arm recovered under turning
// ===========================================================================
TEST_CASE("phase2 convergence VsFusedBase: recovers planted roll + lever arm") {
    const Trajectory traj = multiaxis_turn();
    const Scalar yaw_p = 0.20, pitch_p = -0.12, roll_p = 0.35;
    const Vec3 t_p(0.40, -0.30, 0.15);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_p);

    SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.0;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    Config c = make_cfg(); c.phase2_strategy = Phase2Strat::VsFusedBase;
    REQUIRE(cal.configure(c, 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
    // Feed Phase 1's recovered yaw/pitch (here the planted truth — Phase 1 converged).
    REQUIRE(cal.set_yaw_pitch(1, yaw_pitch_R(yaw_p, pitch_p)) == Status::Ok);

    const DriveResult dr = drive(cal, traj, ref, planted, 1.0, 0.05, 5.95, 50.0);
    REQUIRE(dr.voted > 100);

    CHECK(near_abs(cal.roll(1), roll_p, 2e-2));
    bool obs = false;
    const Vec3 tx = cal.solve_lever_arm(1, &obs);
    CHECK(obs);
    CHECK(near_abs(tx.x(), t_p.x(), 3e-2));
    CHECK(near_abs(tx.y(), t_p.y(), 3e-2));
    CHECK(near_abs(tx.z(), t_p.z(), 3e-2));
    // Histogram-mode lever arm matches too.
    const Vec3 la = cal.lever_arm(1);
    CHECK(near_abs(la.x(), t_p.x(), 5e-2));
    CHECK(near_abs(la.y(), t_p.y(), 5e-2));
    CHECK(near_abs(la.z(), t_p.z(), 5e-2));
    // Full extrinsic rotation matches the planted X.R.
    const SE3 X = cal.extrinsic(1);
    const Vec3 rerr = so3::log(X_planted.R.transpose() * X.R);
    CHECK(rerr.norm() < 3e-2);

    CHECK(cal.extrinsic_confidence(1) > 0.3);       // roll concentration
    CHECK(cal.translation_confidence(1) > 0.2);     // xyz concentration
}

// ===========================================================================
// Convergence — PairwisePinnedRef: same planted roll + lever arm recovered
// ===========================================================================
TEST_CASE("phase2 convergence PairwisePinnedRef: recovers planted roll + lever arm") {
    const Trajectory traj = multiaxis_turn();
    const Scalar yaw_p = 0.20, pitch_p = -0.12, roll_p = 0.35;
    const Vec3 t_p(0.40, -0.30, 0.15);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_p);

    // Reference: IDENTITY mount (pinned at the prior == identity, the gauge).
    SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.0;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    Config c = make_cfg(); c.phase2_strategy = Phase2Strat::PairwisePinnedRef;
    REQUIRE(cal.configure(c, 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}) == Status::Ok);   // reference pinned at identity
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
    REQUIRE(cal.set_yaw_pitch(1, yaw_pitch_R(yaw_p, pitch_p)) == Status::Ok);

    const DriveResult dr = drive(cal, traj, ref, planted, 1.0, 0.05, 5.95, 50.0);
    REQUIRE(dr.voted > 100);

    CHECK(near_abs(cal.roll(1), roll_p, 3e-2));
    bool obs = false;
    const Vec3 tx = cal.solve_lever_arm(1, &obs);
    CHECK(obs);
    CHECK(near_abs(tx.x(), t_p.x(), 3e-2));
    CHECK(near_abs(tx.y(), t_p.y(), 3e-2));
    CHECK(near_abs(tx.z(), t_p.z(), 3e-2));
    CHECK(cal.extrinsic_confidence(1) > 0.3);
    CHECK(cal.translation_confidence(1) > 0.2);
}

// ===========================================================================
// Observability self-test (LOAD-BEARING): straight-only -> NO convergence
// ===========================================================================
TEST_CASE("phase2 observability: straight-only does NOT converge (frozen at prior)") {
    Trajectory traj = Trajectory::straight(2.0, 6.0);     // ||omega|| == 0
    const Scalar yaw_p = 0.20, pitch_p = -0.12, roll_p = 0.35;
    const Vec3 t_p(0.40, -0.30, 0.15);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_p);

    SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.0;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
    REQUIRE(cal.set_yaw_pitch(1, yaw_pitch_R(yaw_p, pitch_p)) == Status::Ok);

    const DriveResult dr = drive(cal, traj, ref, planted, 1.0, 0.05, 5.95, 50.0);

    // The turn gate rejected every step: no votes, no LS rows.
    CHECK(dr.voted == 0);
    CHECK(cal.roll_vote_count(1) == doctest::Approx(0.0));
    CHECK(cal.xyz_vote_count(1)  == doctest::Approx(0.0));
    // Roll stays at the prior (0); lever arm stays at the prior (identity t == 0).
    CHECK(near_abs(cal.roll(1), 0.0, 1e-9));
    bool obs = true;
    const Vec3 tx = cal.solve_lever_arm(1, &obs);
    CHECK_FALSE(obs);                                  // LS has no rows -> not observable
    CHECK(near_abs(tx.norm(), 0.0, 1e-9));             // prior t == 0
    CHECK(cal.extrinsic_confidence(1)   == doctest::Approx(0.0));
    CHECK(cal.translation_confidence(1) == doctest::Approx(0.0));
}

// ===========================================================================
// Yaw-only turning: t_X.x/.y observable, t_X.z NOT (the singular axis)
// ===========================================================================
TEST_CASE("phase2 yaw-only: recovers t_X.x/.y but z is unobservable (low z confidence)") {
    Trajectory traj = Trajectory::turning(2.0, 0.6, 6.0);     // PURE yaw rate (about z)
    const Scalar yaw_p = 0.0, pitch_p = 0.0, roll_p = 0.25;
    const Vec3 t_p(0.40, -0.30, 0.20);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_p);

    SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.0;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    Config c = make_cfg();
    REQUIRE(cal.configure(c, 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
    REQUIRE(cal.set_yaw_pitch(1, yaw_pitch_R(yaw_p, pitch_p)) == Status::Ok);

    const DriveResult dr = drive(cal, traj, ref, planted, 1.0, 0.05, 5.95, 50.0);
    REQUIRE(dr.voted > 100);

    // Roll still recovered (rotation hand-eye is fine under yaw turning).
    CHECK(near_abs(cal.roll(1), roll_p, 3e-2));

    // The LS is rank-deficient (z is the null space of R_A − I for a z-axis rotation):
    // the conditioning guard rejects the full solve.
    bool obs = true;
    cal.solve_lever_arm(1, &obs);
    CHECK_FALSE(obs);

    // The observable x/y channels still concentrate near the planted values; z does not.
    const Vec3 la = cal.lever_arm(1);
    CHECK(near_abs(la.x(), t_p.x(), 5e-2));
    CHECK(near_abs(la.y(), t_p.y(), 5e-2));
    // z confidence is low (never concentrated) -> the joint translation_confidence is low.
    const Scalar cz = cal.translation_confidence(1);
    CHECK(cz < 0.2);
}

// ===========================================================================
// Lever-arm needs rotation: more omega -> tighter t_X; low-omega rows rejected
// ===========================================================================
TEST_CASE("phase2 lever-arm needs rotation: higher turn rate converges tighter") {
    const Scalar yaw_p = 0.10, pitch_p = 0.08, roll_p = 0.2;
    const Vec3 t_p(0.35, -0.25, 0.18);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_p);

    auto run = [&](Scalar wy, Scalar wz) {
        Trajectory tr;
        Vec6 a; a << 2.0, 0, 0, 0, wy,  wz;
        Vec6 b; b << 2.0, 0, 0, 0, -wy, wz;
        tr.add_segment(a, 3.0);
        tr.add_segment(b, 3.0);
        SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
        SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.0;
        SyntheticSource ref(tr, pr);
        SyntheticSource planted(tr, pp);
        auto calp = make_calib();
        calp->configure(make_cfg(), 0);
        calp->set_prior(0, SE3{});
        calp->set_prior(1, SE3{});
        calp->set_yaw_pitch(1, yaw_pitch_R(yaw_p, pitch_p));
        // Add a small per-window noise so the "tighter with more excitation" effect is
        // genuine (a noiseless oracle recovers t_X exactly at any turn rate).
        SourceParams pp2 = pp; pp2.noise_trans_floor = 0.002; pp2.noise_rot_floor = 0.002;
        pp2.seed = 7;
        SyntheticSource planted_noisy(tr, pp2);
        drive(*calp, tr, ref, planted_noisy, 1.0, 0.05, 5.95, 50.0);
        bool obs = false;
        const Vec3 tx = calp->solve_lever_arm(1, &obs);
        Scalar err = (tx - t_p).norm();
        return std::make_pair(obs, err);
    };

    const auto lo = run(0.15, 0.25);    // gentle turn
    const auto hi = run(0.5, 0.9);      // hard turn (more omega)
    REQUIRE(lo.first);
    REQUIRE(hi.first);
    // More rotation -> a better-conditioned LS -> a tighter lever-arm estimate.
    CHECK(hi.second < lo.second);
}

TEST_CASE("phase2 near-singular rows rejected: a low-omega window adds no LS rows") {
    // Two sources; feed turning windows interleaved with low-omega (below the gate /
    // below the rotation-row floor) windows. The low-omega windows must NOT add rows.
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}) == Status::Ok);

    SourceId ids[1] = {0};
    SE3      rep[1];
    rep[0].t = Vec3(0.04, 0, 0);

    // A turning window: omega above the gate, real rotation in A.
    SE3 A_turn; A_turn.R = so3::exp(Vec3(0, 0, 0.02)); A_turn.t = Vec3(0.04, 0, 0);
    rep[0].R = A_turn.R;
    for (int k = 0; k < 10; ++k)
        REQUIRE(cal.observe(1, ids, rep, A_turn, Vec3(0, 0, 1.0)) == Status::Ok);
    const Scalar rows_after_turn = cal.xyz_vote_count(0);
    CHECK(rows_after_turn == doctest::Approx(10.0));

    // A near-straight window: omega below the gate -> gated out, no rows added.
    SE3 A_straight; A_straight.R = Mat3::Identity(); A_straight.t = Vec3(0.04, 0, 0);
    rep[0].R = Mat3::Identity();
    for (int k = 0; k < 10; ++k)
        CHECK(cal.observe(1, ids, rep, A_straight, Vec3(0, 0, 0.0)) == Status::NotReady);
    CHECK(cal.xyz_vote_count(0) == doctest::Approx(rows_after_turn));   // unchanged
}

// ===========================================================================
// Determinism
// ===========================================================================
TEST_CASE("phase2 determinism: identical input -> identical estimate") {
    const Trajectory traj = multiaxis_turn();
    const SE3 X_planted = make_extrinsic(0.2, -0.1, 0.3, Vec3(0.4, -0.3, 0.15));

    auto run_once = [&]() {
        SourceParams pr; pr.id = 0; pr.X = SE3{};     pr.scale = 1.0;
        SourceParams pp; pp.id = 1; pp.X = X_planted; pp.scale = 1.0;
        SyntheticSource ref(traj, pr);
        SyntheticSource planted(traj, pp);
        auto calp = make_calib();
        calp->configure(make_cfg(), 0);
        calp->set_prior(0, SE3{});
        calp->set_prior(1, SE3{});
        calp->set_yaw_pitch(1, yaw_pitch_R(0.2, -0.1));
        drive(*calp, traj, ref, planted, 1.0, 0.05, 5.95, 50.0);
        const Vec3 la = calp->lever_arm(1);
        return std::vector<Scalar>{calp->roll(1), la.x(), la.y(), la.z(),
                                   calp->extrinsic_confidence(1),
                                   calp->translation_confidence(1)};
    };

    const std::vector<Scalar> a = run_once();
    const std::vector<Scalar> b = run_once();
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(a[i] == b[i]);
}

// ===========================================================================
// Vote weight (D5): the configured mode actually changes the effective vote mass.
// Under `Rotation`, high-omega samples dominate the mode vs `One`.
// ===========================================================================
TEST_CASE("phase2 vote_weight: Rotation lets high-omega samples dominate the roll mode") {
    // Plant a known roll into R_B via R_B = Rx(-roll) R_A Rx(roll) (with R_yp = I,
    // best_roll then recovers exactly `roll`). Feed TWO interleaved roll regimes:
    //   - roll rA at a LOW turn rate, with MANY windows (more raw COUNT), and
    //   - roll rB at a HIGH turn rate, with FEWER windows (more ROTATION mass).
    // The turn-gate / Rotation factor read `fused_omega` (passed independently of R_A's
    // own rotation), so the two regimes differ only in their ω magnitude + repetition.
    //   * One:      every vote weighs 1 -> the more-numerous rA wins the mode.
    //   * Rotation: each vote weighs ‖ω‖ -> the high-ω rB carries more mass and wins,
    //               despite fewer windows. The mode flips — proving the weight is honored.
    const Scalar rA = -0.30, rB = 0.40;     // distinct rolls, distinct bins
    // A real (small) base rotation so the per-window ‖log R_A‖ floor passes.
    const Mat3 R_A = so3::exp(Vec3(0.0, 0.05, 0.10));
    // Return Mat3 (NOT auto) — an Eigen product expression would dangle on the Rx()
    // temporaries. best_roll(I, R_A, R_B=Rx(-roll)·R_A·Rx(roll)) recovers exactly `roll`.
    auto planted_RB = [&](Scalar roll) -> Mat3 {
        return Rx(-roll) * R_A * Rx(roll);
    };

    auto run = [&](VoteWeight mode) {
        auto calp = make_calib();
        Config c = make_cfg();
        c.vote_weight = mode;
        REQUIRE(calp->configure(c, 0) == Status::Ok);
        REQUIRE(calp->set_prior(1, SE3{}) == Status::Ok);
        REQUIRE(calp->set_yaw_pitch(1, Mat3::Identity()) == Status::Ok);  // R_yp = I

        SourceId ids[1] = {1};
        SE3 rep[1];
        rep[0].t = Vec3(0.0, 0, 0);          // translation irrelevant to the roll vote here
        // rA regime: low omega (just above the 0.20 gate), MANY windows.
        rep[0].R = planted_RB(rA);
        SE3 A; A.R = R_A; A.t = Vec3(0.1, 0, 0);
        for (int k = 0; k < 60; ++k)
            REQUIRE(calp->observe(1, ids, rep, A, Vec3(0, 0, 0.25)) == Status::Ok);
        // rB regime: HIGH omega, FEWER windows.
        rep[0].R = planted_RB(rB);
        for (int k = 0; k < 20; ++k)
            REQUIRE(calp->observe(1, ids, rep, A, Vec3(0, 0, 3.0)) == Status::Ok);
        return calp->roll(1);
    };

    const Scalar roll_one = run(VoteWeight::One);
    const Scalar roll_rot = run(VoteWeight::Rotation);

    // One: the more-numerous rA wins (60 vs 20 votes, all weight 1).
    CHECK(near_abs(roll_one, rA, 2e-2));
    // Rotation: rB's high-omega mass dominates (20 * 3.0 = 60 vs 60 * 0.25 = 15) -> mode rB.
    CHECK(near_abs(roll_rot, rB, 2e-2));
    // The two modes are genuinely different -> the weight changed the effective vote mass.
    CHECK(std::abs(roll_rot - roll_one) > 0.5);
}

// ===========================================================================
// Committed lever_arm() ridge bias (MINOR): the COMMITTED histogram-mode lever_arm()
// uses a prior-centred ridge solve (kVoteRidgeRel). Pin it with a NONZERO prior to
// bound the prior-pull bias — a larger ridge regression would be caught here.
// ===========================================================================
TEST_CASE("phase2 committed lever_arm ridge bias is bounded with a nonzero prior") {
    const Trajectory traj = multiaxis_turn();
    const Scalar yaw_p = 0.20, pitch_p = -0.12, roll_p = 0.35;
    const Vec3 t_true(0.40, -0.30, 0.15);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_true);

    SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.0;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    Config c = make_cfg(); c.phase2_strategy = Phase2Strat::VsFusedBase;
    REQUIRE(cal.configure(c, 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}) == Status::Ok);
    // NONZERO prior translation, deliberately OFFSET from the truth by ~0.1 m/axis, so the
    // prior-centred ridge would pull the committed estimate toward this wrong prior. The
    // committed lever_arm() must still land near the TRUE t_X (ridge pull << observable
    // resolution); a materially larger ridge would drag it toward t_prior and fail.
    SE3 prior1; prior1.t = t_true + Vec3(0.10, -0.10, 0.10);
    REQUIRE(cal.set_prior(1, prior1) == Status::Ok);
    REQUIRE(cal.set_yaw_pitch(1, yaw_pitch_R(yaw_p, pitch_p)) == Status::Ok);

    const DriveResult dr = drive(cal, traj, ref, planted, 1.0, 0.05, 5.95, 50.0);
    REQUIRE(dr.voted > 100);

    // The COMMITTED (histogram-mode) lever arm — the path exercised by the estimator — is
    // bounded near the TRUE value despite the offset prior (the ~0.1% ridge pull on a
    // ~0.1 m offset is ~1e-4 m, far inside this tolerance; a large ridge would not be).
    const Vec3 la = cal.lever_arm(1);
    CHECK(near_abs(la.x(), t_true.x(), 3e-2));
    CHECK(near_abs(la.y(), t_true.y(), 3e-2));
    CHECK(near_abs(la.z(), t_true.z(), 3e-2));
    // It is NOT sitting at the (wrong) prior — the observable axes resolved off the data.
    CHECK(std::abs(la.x() - prior1.t.x()) > 5e-2);
    CHECK(std::abs(la.y() - prior1.t.y()) > 5e-2);
    CHECK(std::abs(la.z() - prior1.t.z()) > 5e-2);
}

// ===========================================================================
// End-to-end wiring: the full Estimator/Rig populates the Phase-2 extrinsic.
// ===========================================================================
TEST_CASE("phase2 wiring: rig run fills CalibSnapshot roll + lever arm") {
    // A MIXED trajectory: straight segments (so Phase 1 recovers yaw/pitch) + multi-axis
    // turning segments (so Phase 2 recovers roll + the fully-observable lever arm). This is
    // the faithful deployment scenario — each DOF gets its observability regime.
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0, 0.3,  0.6;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.3, 0.6;
    tr.add_segment(straight, 2.0);
    tr.add_segment(turnA,    2.5);
    tr.add_segment(straight, 1.5);
    tr.add_segment(turnB,    2.5);
    tr.add_segment(straight, 1.5);

    // Sources: an identity reference (id 0) + two clean identity sources (ids 1,2) to give
    // the median a true-base consensus, + a planted source (id 3) with a full mount rotation
    // + a lever arm. Slice 7 is ESTIMATE-ONLY (no feedback / bootstrap loop — that is Slice
    // 8), so to isolate the Phase-2 path end-to-end we give the planted source a PRIOR equal
    // to its true ROTATION (Phase 1 then needs no rotation correction and feeds Phase 2 the
    // correct yaw/pitch basepoint) but a ZERO prior TRANSLATION — so the headline Phase-2
    // DOF, the lever arm, is what the full Estimator path must recover from turning.
    const Scalar yaw_p = 0.12, pitch_p = -0.08, roll_p = 0.30;
    const Vec3 t_p(0.35, -0.25, 0.18);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_p);

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X = X_planted;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    // Config: identity priors for the clean sources; the planted source's prior carries its
    // true ROTATION (zero translation). Phase-2 histograms from make_cfg().
    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) { sensors[i].id = static_cast<SourceId>(i); }  // identity priors
    sensors[3].prior_extrinsic.R = X_planted.R;     // true rotation; t stays zero (recover it)
    Config cfg = make_cfg();
    cfg.max_sources    = 4;
    cfg.fusion_delay_s = 0.05;
    cfg.window_s       = 0.10;
    cfg.timesync_enabled = false;       // isolate the calibration path
    cfg.sensors        = sensors.data();
    cfg.sensor_count   = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 50);

    // Read the planted source's calibration snapshot from the latest result.
    const Result& res = rig.estimator().latest();
    const CalibSnapshot* cs = nullptr;
    for (int i = 0; i < res.source_count; ++i) {
        if (res.calib[i].id == 3) { cs = &res.calib[i]; break; }
    }
    REQUIRE(cs != nullptr);

    // The extrinsic rotation recovered the planted yaw/pitch (Phase 1) ∘ roll (Phase 2).
    const Vec3 rerr = so3::log(X_planted.R.transpose() * cs->extrinsic.R);
    CHECK(rerr.norm() < 6e-2);
    // The lever arm recovered the planted translation.
    CHECK(near_abs(cs->extrinsic.t.x(), t_p.x(), 6e-2));
    CHECK(near_abs(cs->extrinsic.t.y(), t_p.y(), 6e-2));
    CHECK(near_abs(cs->extrinsic.t.z(), t_p.z(), 6e-2));
    // The Phase-2 confidences are exposed and rose.
    CHECK(cs->extrinsic_confidence   > 0.2);
    CHECK(cs->translation_confidence > 0.1);
}

// ===========================================================================
// End-to-end wiring (ROLL-FREE PRIOR): the prior carries the correct yaw/pitch but
// ROLL = 0, while the planted mount has a NONZERO roll. The full Estimator path must
// recover the planted roll through the Phase-1 -> Phase-2 seam — so a broken roll seam
// (e.g. roll left at the prior 0) cannot pass this e2e (the original wiring baked the
// true roll into the prior, making roll recovery vacuous).
// ===========================================================================
TEST_CASE("phase2 wiring roll-free prior: recovers a planted roll through the full seam") {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0, 0.3,  0.6;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.3, 0.6;
    tr.add_segment(straight, 2.0);
    tr.add_segment(turnA,    2.5);
    tr.add_segment(straight, 1.5);
    tr.add_segment(turnB,    2.5);
    tr.add_segment(straight, 1.5);

    const Scalar yaw_p = 0.12, pitch_p = -0.08, roll_p = 0.30;
    const Vec3 t_p(0.35, -0.25, 0.18);
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p, roll_p, t_p);
    const Mat3 R_yp_only = yaw_pitch_R(yaw_p, pitch_p);     // the ROLL-FREE prior rotation

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X = X_planted;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    // PRIOR carries the correct yaw/pitch but ZERO roll (and zero translation). Roll about
    // the forward axis leaves the forward axis fixed, so Phase-1 (forward-axis only) still
    // recovers this yaw/pitch and feeds Phase-2 a roll-free R_yp basepoint — the planted
    // roll must come entirely from the Phase-2 hand-eye residual.
    sensors[3].prior_extrinsic.R = R_yp_only;
    Config cfg = make_cfg();
    cfg.max_sources    = 4;
    cfg.fusion_delay_s = 0.05;
    cfg.window_s       = 0.10;
    cfg.timesync_enabled = false;
    cfg.sensors        = sensors.data();
    cfg.sensor_count   = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 50);

    const Result& res = rig.estimator().latest();
    const CalibSnapshot* cs = nullptr;
    for (int i = 0; i < res.source_count; ++i) {
        if (res.calib[i].id == 3) { cs = &res.calib[i]; break; }
    }
    REQUIRE(cs != nullptr);

    // The recovered FULL rotation matches the planted (yaw/pitch ∘ planted roll).
    const Vec3 rerr = so3::log(X_planted.R.transpose() * cs->extrinsic.R);
    CHECK(rerr.norm() < 6e-2);

    // NON-VACUOUS roll check: extract the recovered roll about the forward axis (the
    // residual rotation R_yp^T · R_recovered is a roll-about-x), and assert it matches the
    // PLANTED roll — not 0 (which a broken seam leaving roll at the prior would give).
    const Vec3 roll_vec = so3::log(R_yp_only.transpose() * cs->extrinsic.R);
    CHECK(near_abs(roll_vec.x(), roll_p, 6e-2));     // recovered roll ≈ planted
    CHECK(std::abs(roll_vec.x()) > 0.15);            // and is genuinely nonzero (not vacuous)
    CHECK(std::abs(roll_vec.y()) < 5e-2);            // residual is a pure roll-about-forward
    CHECK(std::abs(roll_vec.z()) < 5e-2);
}
