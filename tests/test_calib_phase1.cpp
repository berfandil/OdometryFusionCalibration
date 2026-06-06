// Slice 6 tests: Phase-1 calibration — straight-regime yaw/pitch + per-source scale.
//
// The oracle is the sim rig: synthetic sources with PLANTED extrinsic (yaw/pitch) +
// scale, driven on straight() (forward AND reverse) vs turning() trajectories.
//
// Coverage:
//   * Convergence — planted yaw/pitch + scale on a straight fwd+reverse run -> the
//     recovered forward axis / scale converge to the planted values; confidence rises.
//   * Observability self-test (LOAD-BEARING) — a turning()-only run does NOT converge:
//     Phase-1 stays at the prior, low confidence (proves the straight gate guards the
//     whole spine).
//   * Reverse-fold — a reverse straight segment contributes to the SAME forward peak
//     (not a 180°-off second peak): fwd+reverse converges to the same axis as fwd-only.
//   * Reference cross-check on/off behaves.
//   * Determinism — identical input -> identical estimate.
//   * Straight-gate unit test — gate accepts straight, rejects turning.
#include <doctest/doctest.h>

#include "ofc/core/calibration.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"

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

// A Config with sensible Phase-1 histograms: fine bins, SlidingK so a few votes already
// give a sharp peak. so(3) channels span a small tilt range (the mount error is small);
// scale spans a generous ratio range.
Config make_cfg() {
    Config c;
    c.tick_rate_hz   = 50.0;
    c.reference_sensor_id = 0;
    // Straight gate (defaults are fine; make explicit for clarity).
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.reverse_fold       = true;
    c.ref_cross_check    = false;

    c.so3_hist.bins       = 512;
    c.so3_hist.range_min  = -0.8;     // rad — covers a generous yaw/pitch tilt
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
    return c;
}

// Heap-allocate the calibrator (large fixed-capacity histograms would overflow the test
// stack — same pattern as the TimeSync tests).
std::unique_ptr<Phase1Calibrator> make_calib() {
    return std::unique_ptr<Phase1Calibrator>(new Phase1Calibrator());
}

// Yaw-then-pitch sensor->base extrinsic rotation (Rz(yaw) * Ry(pitch)), identity t.
SE3 make_extrinsic(Scalar yaw, Scalar pitch) {
    Mat3 Rz; Rz << std::cos(yaw), -std::sin(yaw), 0,
                   std::sin(yaw),  std::cos(yaw), 0,
                   0,              0,             1;
    Mat3 Ry; Ry << std::cos(pitch), 0, std::sin(pitch),
                   0,               1, 0,
                  -std::sin(pitch), 0, std::cos(pitch);
    SE3 X; X.R = Rz * Ry; X.t = Vec3::Zero();
    return X;
}

// The forward axis a source with planted extrinsic X_true reports, expressed in the BASE
// frame through its PRIOR extrinsic X_prior — the quantity the calibrator recovers:
//   g = R_Xprior * R_Xtrue^{-1} * e_x .
// (With prior == identity this is just R_Xtrue^{-1} * e_x; the oracle pins it.)
Vec3 expected_forward_axis(const SE3& X_true, const SE3& X_prior) {
    const Vec3 ex(1, 0, 0);
    return X_prior.R * X_true.R.transpose() * ex;
}

// Drive the calibrator directly over a straight trajectory by querying two synthetic
// sources (reference + one planted) at a fixed cadence and feeding observe() the
// de-scaled reported deltas + the GT fused twist/translation (clean oracle drive — no
// median). `prior` is the prior extrinsic registered for the planted source. Returns the
// number of voted (straight-gated) steps.
struct DriveResult { int voted = 0; };
DriveResult drive(Phase1Calibrator& cal, const Trajectory& traj,
                  const SyntheticSource& ref, const SyntheticSource& planted,
                  const SE3& prior_ref, const SE3& prior_planted,
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
        // GT fused twist/translation from the trajectory over [t, t+step].
        const SE3 A = se3::compose(se3::inverse(traj.pose(t)), traj.pose(t + step));
        const Scalar dt = static_cast<Scalar>(step) / 1e9;
        const Vec3 fused_omega = so3::log(A.R) / dt;
        const Vec3 fused_trans = A.t;
        // De-scale each source's reported translation by its prior_scale (matches the
        // estimator's pre-median de-scale). Reference prior_scale = 1.
        SE3 br = qr.value().motion;
        SE3 bp = qp.value().motion; bp.t = bp.t / prior_scale_planted;
        ids[0] = ref.id();     rep[0] = br;
        ids[1] = planted.id(); rep[1] = bp;
        const Status st = cal.observe(2, ids, rep, fused_omega, fused_trans);
        if (ok(st)) ++dr.voted;
    }
    return dr;
}
} // namespace

// ===========================================================================
// configure() validation
// ===========================================================================
TEST_CASE("phase1 configure validates and reports state") {
    auto calp = make_calib(); Phase1Calibrator& cal = *calp;
    Config c = make_cfg();
    REQUIRE(cal.configure(c, /*reference=*/0) == Status::Ok);
    CHECK(cal.configured());
    CHECK(cal.reference() == 0);

    // A bad so3 histogram (bins < 4) is rejected.
    Config bad = make_cfg();
    bad.so3_hist.bins = 2;
    auto cal2p = make_calib();
    CHECK(cal2p->configure(bad, 0) != Status::Ok);
}

// ===========================================================================
// Straight gate (unit) — accepts straight, rejects turning
// ===========================================================================
TEST_CASE("phase1 straight gate: accepts straight, rejects turning") {
    auto calp = make_calib(); Phase1Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);

    SourceId ids[1] = {0};
    SE3      rep[1];
    rep[0].R = Mat3::Identity();
    rep[0].t = Vec3(0.1, 0, 0);          // forward translation

    // Straight: small omega, sizable translation -> voted.
    CHECK(cal.observe(1, ids, rep, Vec3(0, 0, 0.0), Vec3(0.1, 0, 0)) == Status::Ok);
    CHECK(cal.gated_straight());

    // Turning: large omega -> gated out (no vote).
    CHECK(cal.observe(1, ids, rep, Vec3(0, 0, 0.5), Vec3(0.1, 0, 0)) == Status::NotReady);
    CHECK_FALSE(cal.gated_straight());

    // Too little translation -> gated out.
    CHECK(cal.observe(1, ids, rep, Vec3(0, 0, 0.0), Vec3(0.001, 0, 0)) == Status::NotReady);
    CHECK_FALSE(cal.gated_straight());
}

// ===========================================================================
// Convergence (the oracle): planted yaw/pitch + scale on straight fwd+reverse
// ===========================================================================
TEST_CASE("phase1 convergence: recovers planted yaw/pitch + scale on straight motion") {
    // Straight forward then reverse (so reverse-fold is exercised in the convergence run).
    Trajectory traj;
    Vec6 fwd; fwd << 2.0, 0, 0, 0, 0, 0;
    Vec6 rev; rev << -2.0, 0, 0, 0, 0, 0;
    traj.add_segment(fwd, 3.0);
    traj.add_segment(rev, 3.0);

    // Reference: identity mount, scale 1. Planted: yaw=0.20, pitch=-0.12, scale=1.25.
    const Scalar yaw_p = 0.20, pitch_p = -0.12, scale_p = 1.25;
    const SE3 X_ref = SE3{};
    const SE3 X_planted = make_extrinsic(yaw_p, pitch_p);

    SourceParams pr; pr.id = 0; pr.X = X_ref;     pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted; pp.scale = scale_p;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase1Calibrator& cal = *calp;
    Config c = make_cfg();
    REQUIRE(cal.configure(c, 0) == Status::Ok);
    // Register priors == IDENTITY (uncalibrated bootstrap): the calibrator recovers the
    // full mount from the prior basepoint.
    REQUIRE(cal.set_prior(0, SE3{}, /*scale_calib=*/true) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}, /*scale_calib=*/true) == Status::Ok);

    const DriveResult dr =
        drive(cal, traj, ref, planted, X_ref, SE3{}, /*prior_scale_planted=*/1.0,
              0.05, 5.95, 50.0);
    REQUIRE(dr.voted > 50);

    // Recovered forward axis matches the oracle's R_Xtrue^{-1} * e_x (prior == identity).
    const Vec3 want = expected_forward_axis(X_planted, SE3{});
    const Vec3 got  = cal.forward_axis(1);
    CHECK(near_abs(got.x(), want.x(), 1e-2));
    CHECK(near_abs(got.y(), want.y(), 1e-2));
    CHECK(near_abs(got.z(), want.z(), 1e-2));

    // Recovered scale matches the planted ratio (ref scale 1), within tol.
    CHECK(near_abs(cal.scale(1), scale_p, 2e-2));

    // Confidence rose well above zero (the histogram concentrated).
    CHECK(cal.extrinsic_confidence(1) > 0.3);
    CHECK(cal.scale_confidence(1)     > 0.3);

    // The reference stays pinned (scale == 1, forward axis == e_x).
    CHECK(near_abs(cal.scale(0), 1.0, 1e-9));
    const Vec3 rf = cal.forward_axis(0);
    CHECK(near_abs(rf.x(), 1.0, 1e-2));
}

// ===========================================================================
// Observability self-test (LOAD-BEARING): turning-only -> NO convergence
// ===========================================================================
TEST_CASE("phase1 observability: turning-only does NOT converge (stays near prior)") {
    Trajectory traj = Trajectory::turning(2.0, 0.6, 6.0);   // ||omega|| well above the gate

    const SE3 X_planted = make_extrinsic(0.20, -0.12);
    SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.25;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase1Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}, true) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}, true) == Status::Ok);

    const DriveResult dr =
        drive(cal, traj, ref, planted, SE3{}, SE3{}, 1.0, 0.05, 5.95, 50.0);

    // The straight gate rejected every step: no votes, estimate stays at the prior.
    CHECK(dr.voted == 0);
    CHECK(cal.vote_count(1) == doctest::Approx(0.0));
    CHECK(cal.extrinsic_confidence(1) == doctest::Approx(0.0));
    // Forward axis is the prior forward (e_x), NOT the planted axis.
    const Vec3 f = cal.forward_axis(1);
    CHECK(near_abs(f.x(), 1.0, 1e-9));
    CHECK(near_abs(f.y(), 0.0, 1e-9));
    CHECK(near_abs(f.z(), 0.0, 1e-9));
    CHECK(near_abs(cal.scale(1), 1.0, 1e-9));   // scale stays at prior
}

// ===========================================================================
// Reverse-fold: forward + reverse converge to the SAME axis as forward-only
// ===========================================================================
TEST_CASE("phase1 reverse-fold: reverse segment hits the same forward peak") {
    const SE3 X_planted = make_extrinsic(0.18, 0.0);

    auto run = [&](const Trajectory& traj, Scalar to_s) {
        SourceParams pr; pr.id = 0; pr.X = SE3{};     pr.scale = 1.0;
        SourceParams pp; pp.id = 1; pp.X = X_planted; pp.scale = 1.0;
        SyntheticSource ref(traj, pr);
        SyntheticSource planted(traj, pp);
        auto calp = make_calib();
        calp->configure(make_cfg(), 0);
        calp->set_prior(0, SE3{}, true);
        calp->set_prior(1, SE3{}, true);
        drive(*calp, traj, ref, planted, SE3{}, SE3{}, 1.0, 0.05, to_s, 50.0);
        return calp->forward_axis(1);
    };

    // Forward-only.
    Trajectory fwd_only = Trajectory::straight(2.0, 4.0);
    const Vec3 axis_fwd = run(fwd_only, 3.95);

    // Forward then reverse: the reverse samples must fold into the SAME forward peak.
    Trajectory fwd_rev;
    Vec6 f; f << 2.0, 0, 0, 0, 0, 0;
    Vec6 r; r << -2.0, 0, 0, 0, 0, 0;
    fwd_rev.add_segment(f, 2.0);
    fwd_rev.add_segment(r, 2.0);
    const Vec3 axis_fr = run(fwd_rev, 3.95);

    // Same recovered axis (no 180°-off second peak dragging it away).
    CHECK(near_abs(axis_fr.x(), axis_fwd.x(), 1e-2));
    CHECK(near_abs(axis_fr.y(), axis_fwd.y(), 1e-2));
    CHECK(near_abs(axis_fr.z(), axis_fwd.z(), 1e-2));
    // And it points forward (+x dominant), confirming the fold (a 180°-off peak would
    // pull x negative or split the histogram).
    CHECK(axis_fr.x() > 0.9);

    // Contrast: with the fold DISABLED the reverse-segment samples point backward
    // (g_obs.x < 0 for this forward-ish mount). They are NOT folded into the forward peak,
    // and — because they are then ≥90° off +e_x — the pi-singularity guard SKIPS them
    // (rather than voting a near-pi so(3) log that would land as edge-clamped boundary-bin
    // mass; see calibration.hpp). So fold-OFF DROPS the reverse half of the votes while
    // fold-ON keeps them: fold-ON has strictly MORE vote mass, yet BOTH recover the correct
    // forward axis (fold-OFF just from the forward votes alone). This replaces the old
    // "180°-off second peak lowers confidence" signature — with the guard there is no second
    // peak, only skipped votes.
    {
        struct FoldRun { Scalar votes; Vec3 axis; Scalar conf; };
        auto run_fold = [&](bool fold) {
            Config c = make_cfg();
            c.reverse_fold = fold;
            Trajectory fr;
            Vec6 f2; f2 << 2.0, 0, 0, 0, 0, 0;
            Vec6 r2; r2 << -2.0, 0, 0, 0, 0, 0;
            fr.add_segment(f2, 2.0);
            fr.add_segment(r2, 2.0);
            SourceParams pr; pr.id = 0; pr.X = SE3{};     pr.scale = 1.0;
            SourceParams pp; pp.id = 1; pp.X = X_planted; pp.scale = 1.0;
            SyntheticSource ref(fr, pr);
            SyntheticSource planted(fr, pp);
            auto calp = make_calib();
            calp->configure(c, 0);
            calp->set_prior(0, SE3{}, true);
            calp->set_prior(1, SE3{}, true);
            drive(*calp, fr, ref, planted, SE3{}, SE3{}, 1.0, 0.05, 3.95, 50.0);
            return FoldRun{calp->vote_count(1), calp->forward_axis(1),
                           calp->extrinsic_confidence(1)};
        };
        const FoldRun on  = run_fold(true);
        const FoldRun off = run_fold(false);
        // Fold ON keeps the reverse votes; fold OFF skips them (guarded) -> fewer votes.
        CHECK(on.votes > off.votes);
        CHECK(off.votes > 0.0);             // forward votes still land
        // Both converge to the SAME forward peak (no wrong 180°-off mass either way).
        CHECK(on.conf  > 0.3);
        CHECK(off.conf > 0.3);
        CHECK(near_abs(on.axis.x(),  off.axis.x(), 1e-2));
        CHECK(near_abs(on.axis.y(),  off.axis.y(), 1e-2));
        CHECK(on.axis.x() > 0.9);
        CHECK(off.axis.x() > 0.9);
    }
}

// ===========================================================================
// Consensus-sign reverse-fold + pi-singularity guard
//  * fold decided by the FUSED translation sign (robust for sideways mounts), and
//  * a folded direction ≥90deg off +e_x is SKIPPED (never logs a pi-rotation -> no NaN).
// ===========================================================================
TEST_CASE("phase1 consensus-sign fold: reverse segment folds by FUSED sign, not g_obs.x") {
    // Directly exercise the consensus-sign fold via observe() with hand-built inputs, so the
    // fold rule is tested in ISOLATION (no oracle geometry coupling fused_trans to g_obs).
    // The planted source mounts near-sideways: forward it reports a direction d with x ~ 0,
    // so its own g_obs.x sign is meaningless. On the REVERSE segment its small x lands on the
    // WRONG (positive) side from measurement noise — exactly the case where the OLD g_obs.x
    // fold breaks. The consensus-sign fold must decide forward-vs-reverse from the FUSED
    // translation sign, not g_obs.x, keeping the recovered axis on the single forward peak.
    auto calp = make_calib(); Phase1Calibrator& cal = *calp;
    Config c = make_cfg();
    c.so3_hist.range_min = -2.0;                 // hold the near-sideways tilt (~pi/2)
    c.so3_hist.range_max =  2.0;
    REQUIRE(cal.configure(c, 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}, true) == Status::Ok);   // reference, identity prior
    REQUIRE(cal.set_prior(1, SE3{}, true) == Status::Ok);   // planted, identity prior

    // Near-sideways direction the mount reports moving FORWARD (just inside the pi-guard's
    // 90deg cone). x is small, so its SIGN is meaningless — the regime where folding by
    // g_obs.x is unreliable and the consensus (fused) sign must decide instead.
    Vec3 d(0.10, -0.995, 0.0); d.normalize();          // forward report -> -e_y hemisphere
    // What the SAME mount reports moving REVERSE, with the small x landing on the WRONG
    // (positive) side due to measurement noise: y flips to +e_y, x reads positive. A real
    // antipode would have x<0, but near-sideways the x-sign is dominated by noise.
    Vec3 e(0.15, 0.988, 0.0); e.normalize();           // reverse report, noisy x>0

    SourceId ids[2] = {0, 1};
    SE3 rep[2];
    rep[0].R = Mat3::Identity();
    rep[1].R = Mat3::Identity();
    const Vec3 fwd_t(0.1, 0, 0);
    const Vec3 rev_t(-0.1, 0, 0);

    // Reverse-weighted run (40 forward, 80 reverse) so the reverse votes DOMINATE the
    // histogram. The OLD `g_obs.x<0` fold sees the reverse reports' x>0 and does NOT fold
    // them, piling 80 votes into the WRONG (+e_y) hemisphere -> the recovered axis flips to
    // +e_y (an antipodal split won by the un-folded reverse mass). The consensus-sign fold
    // folds every reverse vote by the fused -x sign; those that fall outside the 90deg cone
    // are pi-guard-SKIPPED, so only the (correct, -e_y) forward votes survive -> a clean
    // -e_y peak on d. This is the behavior the OLD fold cannot produce.
    for (int k = 0; k < 40; ++k) {
        rep[0].t = fwd_t;  rep[1].t = d * Scalar(0.1);
        REQUIRE(cal.observe(2, ids, rep, Vec3::Zero(), fwd_t) == Status::Ok);
    }
    for (int k = 0; k < 80; ++k) {
        rep[0].t = rev_t;  rep[1].t = e * Scalar(0.1);
        REQUIRE(cal.observe(2, ids, rep, Vec3::Zero(), rev_t) == Status::Ok);
    }

    // Single peak in the -e_y hemisphere (the forward report d), NOT pulled to +e_y by the
    // reverse mass. Finite throughout (pi-guard active -> no NaN).
    const Vec3 got = cal.forward_axis(1);
    CHECK(std::isfinite(got.x()));
    CHECK(std::isfinite(got.y()));
    CHECK(std::isfinite(got.z()));
    CHECK(near_abs(got.y(), d.y(), 5e-2));
    CHECK(got.y() < -0.9);                        // consensus fold kept the peak at -e_y
}

TEST_CASE("phase1 pi-guard: a beyond-90deg-off folded direction is skipped (no NaN vote)") {
    // A mount so far off prior that, even after the consensus fold, the observed direction
    // lands in the BACKWARD hemisphere relative to +e_x (cos(e_x, g) < 0). That is outside
    // Phase-1's small-deviation regime and feeding it to rotation_between would build a
    // near-pi rotation whose so3::log is singular (NaN). The guard SKIPS such votes: no
    // NaN ever enters the histogram, and a fully-skipped source stays at its prior.
    //
    // yaw ~ 170 deg => R_Xtrue^{-1} * e_x ~ (cos 170, -sin 170, 0) ~ (-0.98, -0.17, 0):
    // forward-hemisphere fold leaves g ~ (-0.98, ...), i.e. cos(e_x,g) < 0 -> skip.
    const Scalar yaw_p = kPi - 0.17;            // ~170 deg
    const SE3 X_planted = make_extrinsic(yaw_p, 0.0);

    Trajectory traj = Trajectory::straight(2.0, 4.0);
    SourceParams pr; pr.id = 0; pr.X = SE3{};      pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted;  pp.scale = 1.0;
    SyntheticSource ref(traj, pr);
    SyntheticSource planted(traj, pp);

    auto calp = make_calib(); Phase1Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(0, SE3{}, true) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}, true) == Status::Ok);

    drive(cal, traj, ref, planted, SE3{}, SE3{}, 1.0, 0.05, 3.95, 50.0);

    // Every beyond-90deg sample was skipped: no so(3) votes for the planted source, so its
    // estimate stays pinned at the prior (forward axis e_x) and confidence is zero. Crucially
    // the readout is FINITE (a leaked pi-log would have made the mode/axis NaN).
    CHECK(cal.vote_count(1) == doctest::Approx(0.0));
    const Vec3 f = cal.forward_axis(1);
    CHECK(std::isfinite(f.x()));
    CHECK(std::isfinite(f.y()));
    CHECK(std::isfinite(f.z()));
    CHECK(near_abs(f.x(), 1.0, 1e-9));
    CHECK(near_abs(f.y(), 0.0, 1e-9));
    CHECK(near_abs(f.z(), 0.0, 1e-9));
    CHECK(cal.extrinsic_confidence(1) == doctest::Approx(0.0));
}

// ===========================================================================
// Reference cross-check on/off
// ===========================================================================
TEST_CASE("phase1 reference cross-check: on accepts a straight ref, gates a non-straight ref") {
    Trajectory traj = Trajectory::straight(2.0, 4.0);
    const SE3 X_planted = make_extrinsic(0.15, 0.0);
    SourceParams pr; pr.id = 0; pr.X = SE3{};     pr.scale = 1.0;
    SourceParams pp; pp.id = 1; pp.X = X_planted; pp.scale = 1.0;

    // With the cross-check ON and a clean straight reference, votes still happen.
    {
        Config c = make_cfg();
        c.ref_cross_check = true;
        SyntheticSource ref(traj, pr);
        SyntheticSource planted(traj, pp);
        auto calp = make_calib();
        REQUIRE(calp->configure(c, 0) == Status::Ok);
        calp->set_prior(0, SE3{}, true);
        calp->set_prior(1, SE3{}, true);
        const DriveResult dr =
            drive(*calp, traj, ref, planted, SE3{}, SE3{}, 1.0, 0.05, 3.95, 50.0);
        CHECK(dr.voted > 50);
        CHECK(calp->vote_count(1) > 0.0);
    }

    // With the cross-check ON but the reference badly mounted (so its OWN reported
    // direction is far from base-forward, > 30°), the cross-check gates every step out.
    {
        Config c = make_cfg();
        c.ref_cross_check = true;
        // Reference reads ~50° off forward in its own frame (a sideways-pointing mount,
        // with an identity prior so the prior cannot rotate it back to forward).
        SourceParams pr_bad; pr_bad.id = 0; pr_bad.X = make_extrinsic(0.9, 0.0);
        SyntheticSource ref(traj, pr_bad);
        SyntheticSource planted(traj, pp);
        auto calp = make_calib();
        REQUIRE(calp->configure(c, 0) == Status::Ok);
        calp->set_prior(0, SE3{}, true);     // prior identity != the planted bad mount
        calp->set_prior(1, SE3{}, true);
        const DriveResult dr =
            drive(*calp, traj, ref, planted, SE3{}, SE3{}, 1.0, 0.05, 3.95, 50.0);
        CHECK(dr.voted == 0);                 // cross-check rejected the off-axis reference
    }
}

// ===========================================================================
// Determinism
// ===========================================================================
TEST_CASE("phase1 determinism: identical input -> identical estimate") {
    Trajectory traj;
    Vec6 fwd; fwd << 2.0, 0, 0, 0, 0, 0;
    Vec6 rev; rev << -2.0, 0, 0, 0, 0, 0;
    traj.add_segment(fwd, 2.0);
    traj.add_segment(rev, 2.0);
    const SE3 X_planted = make_extrinsic(0.2, -0.1);

    auto run_once = [&]() {
        SourceParams pr; pr.id = 0; pr.X = SE3{};     pr.scale = 1.0;
        SourceParams pp; pp.id = 1; pp.X = X_planted; pp.scale = 1.2;
        SyntheticSource ref(traj, pr);
        SyntheticSource planted(traj, pp);
        auto calp = make_calib();
        calp->configure(make_cfg(), 0);
        calp->set_prior(0, SE3{}, true);
        calp->set_prior(1, SE3{}, true);
        drive(*calp, traj, ref, planted, SE3{}, SE3{}, 1.0, 0.05, 3.95, 50.0);
        Vec3 f = calp->forward_axis(1);
        return std::vector<Scalar>{f.x(), f.y(), f.z(), calp->scale(1),
                                   calp->extrinsic_confidence(1)};
    };

    const std::vector<Scalar> a = run_once();
    const std::vector<Scalar> b = run_once();
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(a[i] == b[i]);
}

// ===========================================================================
// Vote weight (D5): the configured mode changes the effective vote mass in Phase 1 too.
// Under `Confidence`, fewer high-confidence votes outweigh many low-confidence ones.
// ===========================================================================
TEST_CASE("phase1 vote_weight: Confidence lets high-confidence votes dominate the mode") {
    // Drive observe() directly with two competing forward-direction regimes for one source:
    //   - direction dA (small +y tilt) for MANY steps at LOW confidence, then
    //   - direction dB (small -y tilt) for FEWER steps at HIGH confidence.
    // The straight gate passes (omega 0; sizable forward translation). The per-source
    // confidence array carries the Σ-confidence the D5 Confidence/Combo modes weight by.
    //   * One:        all votes weigh 1 -> the more-numerous dA wins the recovered axis.
    //   * Confidence: votes weigh the confidence -> the high-confidence dB wins, despite
    //                 fewer steps. The recovered forward axis flips -> the weight is honored.
    Vec3 dA(0.999, 0.045, 0.0); dA.normalize();      // mostly +x, slight +y
    Vec3 dB(0.999, -0.045, 0.0); dB.normalize();     // mostly +x, slight -y

    auto run = [&](VoteWeight mode) {
        auto calp = make_calib(); Phase1Calibrator& cal = *calp;
        Config c = make_cfg();
        c.vote_weight = mode;
        REQUIRE(cal.configure(c, 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}, true) == Status::Ok);

        SourceId ids[1] = {1};
        SE3 rep[1];
        rep[0].R = Mat3::Identity();
        const Vec3 fused_trans(0.1, 0, 0);
        const Vec3 fused_omega(0, 0, 0);
        // dA regime: MANY steps, LOW confidence.
        Scalar conf_lo[1] = {0.1};
        for (int k = 0; k < 60; ++k) {
            rep[0].t = dA * Scalar(0.1);
            REQUIRE(cal.observe(1, ids, rep, fused_omega, fused_trans, conf_lo)
                        == Status::Ok);
        }
        // dB regime: FEWER steps, HIGH confidence.
        Scalar conf_hi[1] = {1.0};
        for (int k = 0; k < 20; ++k) {
            rep[0].t = dB * Scalar(0.1);
            REQUIRE(cal.observe(1, ids, rep, fused_omega, fused_trans, conf_hi)
                        == Status::Ok);
        }
        return cal.forward_axis(1).y();      // sign of the recovered y-tilt distinguishes A/B
    };

    const Scalar y_one  = run(VoteWeight::One);
    const Scalar y_conf = run(VoteWeight::Confidence);

    // One: the more-numerous dA (+y) wins -> recovered y > 0.
    CHECK(y_one > 0.01);
    // Confidence: dB's high-confidence mass (20*1.0=20 vs 60*0.1=6) wins -> recovered y < 0.
    CHECK(y_conf < -0.01);
    // The recovered axis genuinely flipped -> the weight changed the effective vote mass.
    CHECK(y_one * y_conf < 0.0);
}
