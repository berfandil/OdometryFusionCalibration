// Slice 17b tests: turn-regime JOINT lever + scale (4-unknown hand-eye LS).
// With the observed (prior-de-scaled) t_B = s_res·t_true the hand-eye translation
// identity is LINEAR in u = [t_X; κ], κ = 1/s_res:
//     (R_A − I)·t_X − κ·(R_X t_B) = −t_A      (row block J = [(R_A − I) | −(R_X t_B)])
// Phase 2 accumulates the 4×4 joint normal equations where today's 3-unknown lever rows
// accumulate; observable lever axes vote the existing xyz channels, an observable κ votes
// s_res = 1/κ̂ into a NEW per-source scale2 histogram; the estimator commits scale2 via
// commit_gate_reanchor and folds EITHER scale path's residual into prior_scale on its
// rising edge (both scale histograms reset — sequential edges, no double-fold).
//
// Coverage (SLICE17B acceptance §3, items 1-6):
//   1. Joint recovery on a clean multi-axis conjugated stream (planted lever + scale):
//      lever < 1e-3 m AND scale < 1e-3 of truth; κ prior-pinned (scale2 conf 0) when
//      translation is negligible.
//   2. Planar (yaw-only) stream: lever xy + scale recovered; lever z stays at prior
//      (per-axis gate) — extends, never weakens, the observability self-tests.
//   3. Scale feedback: scale2 commit folds into prior_scale (two-rate re-anchor), BOTH
//      scale hists reset, no commit thrash / double-fold (mutation guard: dropping the
//      either-path cross-reset lets the second estimator re-fold the same residual and
//      the published scale spikes ~residual² — the bound below fails); Phase-1's own
//      fold path untouched (its tests stay green).
//   4. joint_lever_scale=false (default): the knob defaults false and the scale2 surface
//      is inert (1/0/0) — the 3-unknown path is untouched code (pinned by the existing
//      suite's exact-value tests).
//   5. Slice-17 lever regression with the flag ON: the lever-coupling shape (turn-only,
//      Phase-1 starved) still recovers the lever at unit scale, and scale2 reads ≈ 1.
//   6. Persistence: scale2_committed round-trips (format v3, blob byte-equal after
//      restore); old v2-stamped blobs reject; the config-hash covers the flag (a flip
//      rejects). (The loader key test lives in adapters/tests/test_config_loader.cpp.)
#include <doctest/doctest.h>

#include "ofc/core/calibration.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/persistence.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/source.hpp"

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

Mat3 Rz(Scalar a) { Mat3 R; R << std::cos(a), -std::sin(a), 0,
                                 std::sin(a),  std::cos(a), 0,
                                 0, 0, 1; return R; }
Mat3 Ry(Scalar a) { Mat3 R; R << std::cos(a), 0, std::sin(a),
                                 0, 1, 0,
                                -std::sin(a), 0, std::cos(a); return R; }
Mat3 Rx(Scalar a) { Mat3 R; R << 1, 0, 0,
                                 0, std::cos(a), -std::sin(a),
                                 0, std::sin(a),  std::cos(a); return R; }

// Sensor->base extrinsic = Rz(yaw) Ry(pitch) Rx(roll) — the config_loader convention.
SE3 make_extrinsic(Scalar yaw, Scalar pitch, Scalar roll, const Vec3& t) {
    SE3 X; X.R = Rz(yaw) * Ry(pitch) * Rx(roll); X.t = t;
    return X;
}

// Heap-allocate the calibrator (large fixed-capacity histograms overflow the test stack).
std::unique_ptr<Phase2Calibrator> make_calib() {
    return std::unique_ptr<Phase2Calibrator>(new Phase2Calibrator());
}

// Config with the joint path (+ rot3d, the full-3D-recipe pairing) enabled and the
// Phase-2 histogram shapes the other calib tests use. vote_weight One so vote counts are
// literal counts.
Config make_cfg(bool joint = true, bool rot3d = true, bool centroid = true) {
    Config c;
    c.tick_rate_hz   = 50.0;
    c.reference_sensor_id = 0;
    c.turn_omega_min = 0.20;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.phase2_strategy = Phase2Strat::VsFusedBase;
    c.vote_weight     = VoteWeight::One;
    c.rot3d_enabled   = rot3d;
    c.joint_lever_scale = joint;

    c.so3_hist.bins       = 512;
    c.so3_hist.range_min  = -0.8;
    c.so3_hist.range_max  =  0.8;
    c.so3_hist.aging      = Aging::SlidingK;
    c.so3_hist.sliding_k  = 256;

    c.scale_hist.bins = 512; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 256;

    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 256;

    c.xyz_hist.bins = 512; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 256;

    // The full-3D-recipe readout (design §2: rot3d + joint_lever_scale + subbin_centroid
    // + vote_weight one): the centroid sub-bin estimator is exact for a vote-split point
    // mass, which the <1e-3 acceptance bounds lean on. The Slice-17 regression case
    // turns it OFF to mirror the exact Slice-17 test configuration.
    c.so3_hist.subbin_centroid   = centroid;
    c.roll_hist.subbin_centroid  = centroid;
    c.xyz_hist.subbin_centroid   = centroid;
    c.scale_hist.subbin_centroid = centroid;
    return c;
}

// One conjugated hand-eye window with a PLANTED RESIDUAL SCALE: base motion A (rotation
// `theta` about `axis` + translation `tA`), sensor report B = X^-1 o A o X, then
// B.t *= s_res (the sim's scale model — exactly the unrecovered residual the joint solve
// must absorb into κ = 1/s_res).
Status feed_window(Phase2Calibrator& cal, SourceId id, const SE3& X, Scalar s_res,
                   const Vec3& axis, Scalar theta, const Vec3& tA) {
    SE3 A;
    A.R = so3::exp(axis.normalized() * theta);
    A.t = tA;
    SE3 B = se3::compose(se3::compose(se3::inverse(X), A), X);
    B.t *= s_res;
    const Scalar dt = Scalar(0.5);                      // 0.5 s windows
    const Vec3 omega = so3::log(A.R) / dt;              // > turn gate for theta >= 0.15
    SourceId ids[1] = { id };
    SE3      rep[1] = { B };
    return cal.observe(1, ids, rep, A, omega);
}

// Feed `count` windows. multiaxis cycles 3 genuinely distinct rotation axes; planar
// keeps every window about +z with an xy-only base translation (the ground regime).
void feed_stream(Phase2Calibrator& cal, SourceId id, const SE3& X, Scalar s_res,
                 int count, bool multiaxis, const Vec3& tA = Vec3(0.5, 0.05, -0.02)) {
    const Vec3 axes[3] = { Vec3(0, 0, 1), Vec3(0, 1, 0.2), Vec3(0.3, -0.4, 0.85) };
    for (int k = 0; k < count; ++k) {
        const Vec3 axis = multiaxis ? axes[k % 3] : Vec3(0, 0, 1);
        const Scalar theta = Scalar(0.25) + Scalar(0.10) * std::sin(Scalar(k) * Scalar(0.7));
        feed_window(cal, id, X, s_res, axis, theta, tA);
    }
}

void set_hists(Config& c) {
    c.vote_weight = VoteWeight::One;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.turn_omega_min     = 0.20;
    c.so3_hist.bins = 256; c.so3_hist.range_min = -0.8; c.so3_hist.range_max = 0.8;
    c.so3_hist.aging = Aging::SlidingK; c.so3_hist.sliding_k = 200;
    c.scale_hist.bins = 256; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 200;
    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 200;
    c.xyz_hist.bins = 256; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 200;
    c.offset_hist.bins = 256; c.offset_hist.range_min = -0.1; c.offset_hist.range_max = 0.1;
    c.offset_hist.aging = Aging::SlidingK; c.offset_hist.sliding_k = 200;
}

const CalibSnapshot* snap(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.calib[i].id == id) return &r.calib[i];
    return nullptr;
}

// TURN-ONLY multi-axis trajectory — the EuRoC/drone regime: NO straight segments, so
// Phase-1 (and with it the straight-regime scale path) STARVES; only the turn-regime
// joint path can ever recover scale here. RICHLY EXCITED (4 distinct turn types with
// varying speed/axis): the joint [t_X; κ] system needs the κ column −R_X t_B to vary
// relative to the rotation excitation — two constant-twist window types leave a
// lever/κ trade direction near-null (the ridge then dumps a planted scale into the
// lever; measured: lever err 0.12 m on a two-type trajectory vs 0.004 m here).
Trajectory multiaxis_turnonly_traj() {
    Trajectory tr;
    Vec6 t1; t1 << 2.0,  0,   0,    0,     0.35,  0.6;
    Vec6 t2; t2 << 1.2,  0.2, 0,    0,    -0.35,  0.6;
    Vec6 t3; t3 << 2.6,  0,   0.25, 0,     0.20, -0.5;
    Vec6 t4; t4 << 1.5, -0.2, 0,    0.25, -0.30, -0.45;
    for (int rep = 0; rep < 3; ++rep) {
        tr.add_segment(t1, 0.9);
        tr.add_segment(t2, 0.9);
        tr.add_segment(t3, 0.9);
        tr.add_segment(t4, 0.9);
    }
    return tr;
}
} // namespace

// ===========================================================================
// 1. Joint recovery: clean multi-axis stream, planted lever + scale
// ===========================================================================
TEST_CASE("joint lever+scale: clean multi-axis stream recovers planted lever < 1e-3 m "
          "and scale < 1e-3; kappa prior-pins at negligible translation") {
    // EuRoC-magnitude planted mount + the Slice-17 finding's scale: 1.08.
    const SE3 X = make_extrinsic(8 * kPi / 180, 5 * kPi / 180, 4 * kPi / 180,
                                 Vec3(0.10, -0.05, 0.20));
    const Scalar s_res = 1.08;

    {
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
        // Prior: TRUE rotation (the rot3d/R_yp driver — isolates the translation LS),
        // ZERO lever (wrong), unit scale prior (κ prior = 1; the 1.08 is the residual).
        SE3 Xp; Xp.R = X.R; Xp.t = Vec3::Zero();
        REQUIRE(cal.set_prior(1, Xp) == Status::Ok);

        feed_stream(cal, 1, X, s_res, 300, /*multiaxis=*/true);

        // Direct joint solve (diagnostics readout): essentially exact on clean data (κ is
        // informed, so the solve is the raw joint LS — no ridge bias; the ~1e-6 floor is
        // the FIRST window's pre-rot3d-gate roll-refine resolution in its row R_X).
        bool obs = false;
        const Vec3 t_solve = cal.solve_lever_arm(1, &obs);
        CHECK(obs);
        INFO("solve_lever_arm err = " << (t_solve - X.t).norm() << " m");
        CHECK((t_solve - X.t).norm() < 1e-5);

        // Committed readouts (histogram modes): the acceptance bounds.
        const Vec3 lever = cal.lever_arm(1);
        INFO("lever err = " << (lever - X.t).norm() << " m");
        CHECK((lever - X.t).norm() < 1e-3);
        INFO("scale2 = " << cal.scale2(1));
        CHECK(std::abs(cal.scale2(1) - s_res) < 1e-3);
        CHECK(cal.scale2_confidence(1) > Scalar(0.5));
        CHECK(cal.scale2_vote_count(1) > Scalar(100));
        CHECK(cal.translation_confidence(1) > Scalar(0.5));
    }

    // κ PRIOR-PIN: negligible translation (zero base translation AND zero lever ⇒
    // t_B ≡ 0): the κ diagonal carries no information, so the scale2 channel NEVER
    // votes — scale2 stays at the prior 1 with confidence 0 (the per-DOF spine).
    {
        const SE3 Xrot = make_extrinsic(0.2, -0.12, 0.3, Vec3::Zero());
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, Xrot) == Status::Ok);

        feed_stream(cal, 1, Xrot, /*s_res=*/1.3, 150, /*multiaxis=*/true,
                    /*tA=*/Vec3::Zero());

        CHECK(cal.scale2(1) == doctest::Approx(1.0));
        CHECK(cal.scale2_confidence(1) == doctest::Approx(0.0));
        CHECK(cal.scale2_vote_count(1) == doctest::Approx(0.0));
    }
}

// ===========================================================================
// 2. Planar observability: lever xy + scale recovered; lever z stays at prior
// ===========================================================================
TEST_CASE("joint lever+scale planar: yaw-only recovers lever xy + scale; z pinned at "
          "prior (per-axis gate)") {
    // Ground rig: yaw-only mount, lever with a NON-ZERO z the planar stream cannot see,
    // KAIST-magnitude residual scale 1.10. rot3d stays gated (rank-1), so the row R_X is
    // the R_yp ∘ Rx(roll) composition at the (true-rotation) prior.
    const SE3 X = make_extrinsic(0.30, 0, 0, Vec3(0.20, -0.15, 0.10));
    const Scalar s_res = 1.10;

    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    SE3 Xp; Xp.R = X.R; Xp.t = Vec3::Zero();             // true rotation, zero lever
    REQUIRE(cal.set_prior(1, Xp) == Status::Ok);

    feed_stream(cal, 1, X, s_res, 300, /*multiaxis=*/false,
                /*tA=*/Vec3(0.5, 0.05, 0.0));            // planar (xy) base translation

    CHECK_FALSE(cal.rot3d_observable(1));                 // genuinely rank-1 motion

    // xy + scale recovered; z NEVER voted -> exactly the prior component (0).
    const Vec3 lever = cal.lever_arm(1);
    INFO("lever = [" << lever.x() << ", " << lever.y() << ", " << lever.z() << "]");
    CHECK(std::abs(lever.x() - X.t.x()) < 5e-3);
    CHECK(std::abs(lever.y() - X.t.y()) < 5e-3);
    CHECK(lever.z() == Scalar(0));                        // prior, exactly (empty channel)
    INFO("scale2 = " << cal.scale2(1));
    CHECK(std::abs(cal.scale2(1) - s_res) < 5e-3);
    CHECK(cal.scale2_confidence(1) > Scalar(0.5));
    // The weakest (unobservable z) channel bounds the joint translation confidence — the
    // existing observability self-test shape, extended, never weakened.
    CHECK(cal.translation_confidence(1) == doctest::Approx(0.0));
    // The full-solve guard reports unobservable (the z null direction survives the
    // κ-marginalization) and falls back to the prior.
    bool obs = true;
    const Vec3 t_solve = cal.solve_lever_arm(1, &obs);
    CHECK_FALSE(obs);
    CHECK(t_solve.x() == Scalar(0));
    CHECK(t_solve.z() == Scalar(0));
}

// ===========================================================================
// 1b. Mutation guard on the κ column: WITHOUT the joint κ unknown (the 3-unknown
// status quo) the same scaled stream drags the lever by centimetres — pinning that
// the κ column is load-bearing for the recovery in (1).
// ===========================================================================
TEST_CASE("joint lever+scale: the kappa column is load-bearing (3-unknown control "
          "is dragged by an unrecovered scale; joint is immune)") {
    const SE3 X = make_extrinsic(8 * kPi / 180, 5 * kPi / 180, 4 * kPi / 180,
                                 Vec3(0.10, -0.05, 0.20));
    const Scalar s_res = 1.08;
    SE3 Xp; Xp.R = X.R; Xp.t = Vec3::Zero();

    Scalar err_joint = 0, err_control = 0;
    for (int pass = 0; pass < 2; ++pass) {
        const bool joint = (pass == 0);
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_cfg(joint), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, Xp) == Status::Ok);
        feed_stream(cal, 1, X, s_res, 300, /*multiaxis=*/true);
        const Scalar err = (cal.lever_arm(1) - X.t).norm();
        if (joint) err_joint = err; else err_control = err;
    }
    INFO("lever err: joint = " << err_joint << " m, 3-unknown control = " << err_control);
    CHECK(err_joint < 1e-3);
    CHECK(err_control > 0.01);                 // the Slice-17 3.5 cm class of bias
    CHECK(err_joint < err_control * Scalar(0.1));
}

// ===========================================================================
// 1c. 2-EXCITATION HONESTY (17b review MAJOR-1): exactly TWO constant-twist window
// types leave a MIXED lever/κ trade direction near-null in the joint 4×4 — per-
// diagonal information looks healthy (every raw diagonal is large), but the ridge
// resolves the near-null toward the prior and deterministically dumps the planted
// scale into the lever (measured pre-fix: lever err ~0.12 m, votes CONCENTRATED →
// confidently wrong, worse than flag-off). The joint VOTE path must therefore gate
// on the JOINT conditioning (the per-axis Schur-marginalized information), not the
// raw diagonals: under this trajectory NO joint votes may land — conf 0, vote count
// 0, prior-pinned readouts (honest "cannot observe"), NOT a confident wrong answer.
// MUTATION GUARD: dropping the Schur vote gates in solve_ridge4 re-deposits the
// biased votes and the vote-count/confidence pins below fail.
// ===========================================================================
TEST_CASE("joint lever+scale 2-excitation: constant-twist PAIR yields NO joint votes "
          "(honest prior-pin) instead of a confidently wrong lever/scale") {
    const SE3 X = make_extrinsic(8 * kPi / 180, 5 * kPi / 180, 4 * kPi / 180,
                                 Vec3(0.10, -0.05, 0.20));
    const Scalar s_res = 1.08;
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    SE3 Xp; Xp.R = X.R; Xp.t = Vec3::Zero();            // true rotation, zero lever
    REQUIRE(cal.set_prior(1, Xp) == Status::Ok);

    // EXACTLY two distinct constant-twist window types, alternating: two rotation axes
    // (the rot3d two-axis gate OPENS — the rotation is fine) but only two (A, t_B) row
    // shapes — SAME forward translation, mirrored pitch rate (the measured-0.12 m
    // trajectory pair, see 1d below). The mixed near-null is structural: with axes
    // n± = (0, ±a, b) the relative axis w = axis(R₂ᵀR₁) ≈ e_y and (R₁−I)·w ≈ −θb·e_x
    // is PARALLEL to the x-only translation, so to first order in θ the joint system
    // admits the null family [δt; δκ] = α[w + δκ·s·t_X; δκ] — a lever/κ trade the
    // per-diagonal information cannot see (every raw diagonal stays large).
    for (int k = 0; k < 300; ++k) {
        const bool odd = (k % 2 == 1);
        feed_window(cal, 1, X, s_res,
                    odd ? Vec3(0, -0.35, 0.6) : Vec3(0, 0.35, 0.6), Scalar(0.12),
                    Vec3(0.35, 0, 0));
    }
    CHECK(cal.rot3d_observable(1));                      // rotation IS observable...

    // ...but the null-involved joint axes are NOT (here the null is along e_y + κ —
    // see the geometry note above): κ NEVER votes (scale2 inert, conf 0) and the
    // dragged lever axis y NEVER votes (stays at the prior 0 exactly — pre-fix it read
    // ≈ −0.13, the planted scale dumped into the lever). The honestly-observable x/z
    // marginals (their Schur information survives the null) may keep voting, but their
    // values stay near truth — nothing is confidently wrong, and the empty y channel
    // pins translation_confidence at 0 (the per-DOF spine), so no lever commit either.
    CHECK(cal.scale2_vote_count(1) == doctest::Approx(0.0));
    CHECK(cal.scale2_confidence(1) == doctest::Approx(0.0));
    CHECK(cal.scale2(1) == doctest::Approx(1.0));        // prior (empty histogram)
    CHECK(cal.translation_confidence(1) == doctest::Approx(0.0));
    const Vec3 lever = cal.lever_arm(1);
    INFO("lever = [" << lever.x() << ", " << lever.y() << ", " << lever.z() << "]");
    CHECK(lever.y() == Scalar(0));                       // withheld axis: prior, exactly
    CHECK(std::abs(lever.x() - X.t.x()) < 0.02);         // voted axes: not poisoned
    CHECK(std::abs(lever.z() - X.t.z()) < 0.02);
    // The full-solve guard (kCondMin on the κ-marginalized Schur) agrees: unobservable.
    bool obs = true;
    cal.solve_lever_arm(1, &obs);
    CHECK_FALSE(obs);
}

// ===========================================================================
// 1d. The estimator-level measured failure (17b review MAJOR-1): the SAME
// 2-excitation regime through the full feedback loop. Pre-fix the biased joint
// votes CONCENTRATE, the lever commits ~0.12 m off truth and scale2 commits a
// κ-pinned ≈1 residual — flag-on confidently WORSE than flag-off. Post-fix the
// joint votes are withheld: the mount (prior == truth) stays put and the scale
// stays prior-pinned (honest: this trajectory cannot observe the joint DOFs).
// ===========================================================================
TEST_CASE("joint lever+scale estimator 2-excitation: lever stays at truth and scale "
          "stays prior-pinned (no confident wrong commit)") {
    Trajectory tr;
    Vec6 t1; t1 << 2.0, 0, 0, 0,  0.35, 0.6;             // the measured 0.12 m pair:
    Vec6 t2; t2 << 2.0, 0, 0, 0, -0.35, 0.6;             // same v, mirrored pitch rate
    for (int rep = 0; rep < 6; ++rep) {                  // 2-type constant-twist pair
        tr.add_segment(t1, 0.9);
        tr.add_segment(t2, 0.9);
    }
    const SE3 X3 = make_extrinsic(0.10, -0.07, 0.09, Vec3(0.15, -0.10, 0.05));
    const Scalar planted_scale = 1.08;

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X     = X3;
    planted[3].scale = planted_scale;
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = X3;        // prior == truth: any motion is vote bias

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 120;
    cfg.rot3d_enabled     = true;
    cfg.joint_lever_scale = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    // HONEST: the under-excited joint DOFs never commit a wrong value. The mount
    // holds the (true) prior — pre-fix it walked ~0.12 m off — and the published
    // scale stays at the (unobservable-here) prior 1.
    CHECK_FALSE(rig.estimator().scale2_committed(3));
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    INFO("published scale = " << cs->scale
         << ", lever err = " << (cs->extrinsic.t - X3.t).norm() << " m");
    CHECK((cs->extrinsic.t - X3.t).norm() < 0.02);
    CHECK(std::abs(cs->scale - Scalar(1)) < 0.02);
}

// ===========================================================================
// 5. Slice-17 lever regression with the flag ON: unit scale, turn-only,
// Phase-1 starved — the joint solve must not degrade the unit-scale cases.
// ===========================================================================
TEST_CASE("joint lever+scale: Slice-17 lever-coupling shape at UNIT scale still "
          "recovers (flag ON does not degrade unit-scale levers)") {
    const SE3 X = make_extrinsic(0.25, -0.18, 0.12, Vec3(0.20, -0.15, 0.10));
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    // EXACT Slice-17 test configuration (centroid off) + the joint flag — the regression
    // the acceptance judges (test_calib_rot3d's lever-coupling bound, flag ON).
    REQUIRE(cal.configure(make_cfg(true, true, /*centroid=*/false), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);      // IDENTITY prior — both wrong

    // Turn-only conjugated stream at UNIT scale (the Slice-17 rig, joint flag ON).
    // 400 windows (vs the Slice-17 case's 250; review MINOR — margin banking): with
    // SlidingK 256 the first ~150 pre-R̂-convergence transient votes age fully OUT of
    // the ring, so the mode reads the converged tail instead of mixing the wrong-R_X
    // bootstrap votes in. Same excitation SHAPE as Slice 17 (the regression intent),
    // just enough extra windows that the 0.02 bound carries a real margin instead of
    // passing at ~0.019.
    feed_stream(cal, 1, X, /*s_res=*/1.0, 400, /*multiaxis=*/true);

    CHECK(cal.rot3d_observable(1));
    CHECK(so3::log(cal.rot3d(1) * X.R.transpose()).norm() < 1e-3);
    const Vec3 lever = cal.lever_arm(1);
    INFO("lever err = " << (lever - X.t).norm() << " m");
    // The Slice-17 bound. Measured 0.0148 at 400 windows (0.019 at the original 250 —
    // a 5% margin the review flagged as brittle): the residual is the shape's floor
    // (centroid-off bin-center quantization ~6 mm bins + the R̂-driven row residual),
    // not joint-path bias, and 0.02 now carries ~35% headroom for benign numeric drift.
    CHECK((lever - X.t).norm() < 0.02);
    CHECK(cal.translation_confidence(1) > Scalar(0.5));
    // κ reads "no residual scale": scale2 ≈ 1.
    CHECK(std::abs(cal.scale2(1) - Scalar(1)) < 0.02);
}

// ===========================================================================
// Range guard, calibrator level (17b review MAJOR-2/MAJOR-3): a TRUE residual scale
// outside scale_hist's range is SKIPPED — never deposited as edge-clamped mass. The
// "replace skip with clamp" mutation deposits ~300 deterministic edge votes here:
// vote_count explodes, the mode reads ~1.5 at high confidence (a confidently WRONG
// scale that would commit and irreversibly poison prior_scale at the estimator) and
// the skipped counter stays 0 — every CHECK below fails. Also pins the documented
// LIMITATION: the out-of-range truth is never recovered by this path (votes 0), but
// it is DIAGNOSABLE (skipped high) and costs only the scale vote — the joint solve
// absorbs κ̂ = 0.5 internally, so the lever still recovers (scale-immune).
// ===========================================================================
TEST_CASE("joint lever+scale range guard: out-of-range TRUE residual (s_res = 2.0) "
          "never votes — skip not clamp; skipped counter grows; lever unharmed") {
    const SE3 X = make_extrinsic(8 * kPi / 180, 5 * kPi / 180, 4 * kPi / 180,
                                 Vec3(0.10, -0.05, 0.20));
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    SE3 Xp; Xp.R = X.R; Xp.t = Vec3::Zero();
    REQUIRE(cal.set_prior(1, Xp) == Status::Ok);

    // RICHLY excited stream (κ is fully observable — only the RANGE withholds votes)
    // with the planted residual at 2.0, outside the configured (0.5, 1.5). 400 windows
    // so the early running-solve transient (κ̂ converging toward 0.5 biases the first
    // lever votes) ages out of the SlidingK-256 ring before the mode is read.
    feed_stream(cal, 1, X, /*s_res=*/2.0, 400, /*multiaxis=*/true);

    CHECK(cal.scale2_vote_count(1) == doctest::Approx(0.0));   // skip, NOT clamp
    CHECK(cal.scale2(1) == doctest::Approx(1.0));              // prior (empty histogram)
    CHECK(cal.scale2_confidence(1) == doctest::Approx(0.0));
    CHECK(cal.scale2_skipped(1) > Scalar(350));                // diagnosable, not silent
    // The lever is still recovered (the joint solve is scale-immune): the 2× scale
    // costs only the (out-of-range) scale vote, not centimetres of lever.
    INFO("lever err = " << (cal.lever_arm(1) - X.t).norm() << " m");
    CHECK((cal.lever_arm(1) - X.t).norm() < 0.02);
}

// ===========================================================================
// reset_scale2 / lever-gate decoupling (17b review): a scale fold drops the scale2
// votes + the joint ACCUMULATOR epoch, but must KEEP the xyz histograms AND the
// lever commit gate's vote-mass input xyz_vote_count — pre-fix reset_scale2 zeroed
// rows_ (which doubles as xyz_vote_count), erasing a pending lever commit's entire
// gate progress on every fold even though the votes it measures were kept.
// ===========================================================================
TEST_CASE("joint lever+scale reset_scale2: keeps xyz histograms AND the lever gate's "
          "vote mass; drops only the scale2 votes + the accumulator epoch") {
    const SE3 X = make_extrinsic(8 * kPi / 180, 5 * kPi / 180, 4 * kPi / 180,
                                 Vec3(0.10, -0.05, 0.20));
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    SE3 Xp; Xp.R = X.R; Xp.t = Vec3::Zero();
    REQUIRE(cal.set_prior(1, Xp) == Status::Ok);

    feed_stream(cal, 1, X, /*s_res=*/1.08, 150, /*multiaxis=*/true);
    const Scalar votes_before = cal.xyz_vote_count(1);
    const Vec3   lever_before = cal.lever_arm(1);
    REQUIRE(votes_before > Scalar(100));
    REQUIRE(cal.scale2_vote_count(1) > Scalar(0));

    cal.reset_scale2(1);

    // The lever side is untouched: gate vote mass + histogram readout EXACTLY equal.
    CHECK(cal.xyz_vote_count(1) == votes_before);        // pre-fix: dropped to 0
    CHECK(cal.lever_arm(1).x() == lever_before.x());
    CHECK(cal.lever_arm(1).y() == lever_before.y());
    CHECK(cal.lever_arm(1).z() == lever_before.z());
    // The scale2 votes are gone (fall back to the prior 1)...
    CHECK(cal.scale2_vote_count(1) == doctest::Approx(0.0));
    CHECK(cal.scale2(1) == doctest::Approx(1.0));
    // ...and so is the accumulator EPOCH: the direct solve falls back to the prior
    // until post-fold rows land (the single-κ model cannot mix de-scale epochs).
    bool obs = true;
    const Vec3 t_solve = cal.solve_lever_arm(1, &obs);
    CHECK_FALSE(obs);
    CHECK(t_solve.norm() < 1e-12);                       // prior (zero lever)

    // Post-fold rows (residual now ≈ 1, the re-anchored convention) resume cleanly:
    // scale2 re-votes ≈ 1 with no stale-epoch blending.
    feed_stream(cal, 1, X, /*s_res=*/1.0, 100, /*multiaxis=*/true);
    CHECK(cal.scale2_vote_count(1) > Scalar(50));
    CHECK(std::abs(cal.scale2(1) - Scalar(1)) < 5e-3);
    CHECK(cal.xyz_vote_count(1) == votes_before + Scalar(100));
}

// ===========================================================================
// 4. Default-off: the knob defaults false; the scale2 surface is inert; reset_scale2
// is a no-op on the legacy accumulator. (The 3-unknown path is untouched CODE — its
// byte-identical behavior is pinned by the existing suite's exact-value tests.)
// ===========================================================================
TEST_CASE("joint lever+scale default-off: knob false by default; scale2 surface inert") {
    CHECK_FALSE(Config{}.joint_lever_scale);

    const SE3 X = make_extrinsic(0.2, -0.1, 0.3, Vec3(0.4, -0.3, 0.15));
    auto cal_off = make_calib();
    auto cal_ref = make_calib();
    Config off = make_cfg(/*joint=*/false);
    REQUIRE(cal_off->configure(off, 0) == Status::Ok);
    REQUIRE(cal_ref->configure(off, 0) == Status::Ok);
    REQUIRE(cal_off->set_prior(1, SE3{}) == Status::Ok);
    REQUIRE(cal_ref->set_prior(1, SE3{}) == Status::Ok);

    // A stream that WOULD vote scale2 with the knob on (scaled, translated, turning).
    feed_stream(*cal_off, 1, X, /*s_res=*/1.2, 150, /*multiaxis=*/true);

    CHECK(cal_off->scale2(1) == Scalar(1));               // inert surface, exactly
    CHECK(cal_off->scale2_confidence(1) == Scalar(0));
    CHECK(cal_off->scale2_vote_count(1) == Scalar(0));

    // reset_scale2 with the knob off must NOT disturb the legacy lever state: the
    // readouts before/after are EXACTLY equal (and equal to an untouched twin run).
    const Vec3   lever_before = cal_off->lever_arm(1);
    const Scalar votes_before = cal_off->xyz_vote_count(1);
    cal_off->reset_scale2(1);
    const Vec3 lever_after = cal_off->lever_arm(1);
    CHECK(lever_after.x() == lever_before.x());
    CHECK(lever_after.y() == lever_before.y());
    CHECK(lever_after.z() == lever_before.z());
    CHECK(cal_off->xyz_vote_count(1) == votes_before);

    feed_stream(*cal_ref, 1, X, /*s_res=*/1.2, 150, /*multiaxis=*/true);
    const Vec3 lever_twin = cal_ref->lever_arm(1);
    CHECK(lever_after.x() == lever_twin.x());             // exact-equality twin pin
    CHECK(lever_after.y() == lever_twin.y());
    CHECK(lever_after.z() == lever_twin.z());
    CHECK(cal_off->roll(1) == cal_ref->roll(1));
    CHECK(cal_off->translation_confidence(1) == cal_ref->translation_confidence(1));
}

// ===========================================================================
// 3. Estimator scale feedback: turn-only rig (Phase-1 scale path STARVED) — the
// scale2 commit folds the planted scale into prior_scale; no thrash. Plus (6) the
// persistence v3 round-trip + v2 reject + config-hash flip.
// ===========================================================================
TEST_CASE("joint lever+scale estimator: scale2 commits + folds on a turn-only rig; "
          "persists (v3) and rejects v2/cross-flag blobs") {
    const Trajectory tr = multiaxis_turnonly_traj();
    const SE3 X3 = make_extrinsic(0.10, -0.07, 0.09, Vec3(0.15, -0.10, 0.05));
    const Scalar planted_scale = 1.08;

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X     = X3;                  // sources 0,1,2 identity (consensus anchor)
    planted[3].scale = planted_scale;
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = X3;        // mount prior == truth: isolate the SCALE DOF

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    // commit_min_votes high enough that the scale2 commit sees votes spanning several
    // distinct turn types (the running κ̂'s convergence transient must have left the
    // SlidingK ring's mode before the rising edge folds it — the fold is one-shot).
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 120;
    cfg.rot3d_enabled     = true;
    cfg.joint_lever_scale = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    // COMMITTED + FOLDED: with NO straight regime Phase-1's scale path never fires, so
    // the published absolute scale reaching ~1.08 REQUIRES the scale2 fold.
    CHECK(rig.estimator().scale2_committed(3));
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    INFO("published scale = " << cs->scale);
    CHECK(std::abs(cs->scale - planted_scale) < 0.02);
    CHECK(cs->scale_committed);                     // the scale2 commit is OR-ed in
    // The mount stayed put (prior == truth; the joint solve is scale-immune).
    CHECK((cs->extrinsic.t - X3.t).norm() < 0.02);

    // NO THRASH: once committed the latch holds, and the published scale never
    // overshoots (a double-fold would read ≈ residual² ≈ 1.17).
    bool seen_commit = false;
    for (const auto& rec : rig.records()) {
        if (!rec.fused) continue;
        const CalibSnapshot* c3 = snap(rec.result, 3);
        if (c3 == nullptr) continue;
        if (c3->scale_committed) seen_commit = true;
        else CHECK_FALSE(seen_commit);              // never un-commits after the edge
        CHECK(c3->scale < planted_scale * Scalar(1.06));   // no over-fold, ever
    }
    CHECK(seen_commit);

    // --- (6) Persistence: current-version round-trip; old-stamp reject; config-hash flip ---
    // The format version is now v4 (Slice 20b appended the per-axis lever flags). This test's
    // round-trip + the old-stamp reject below are version-agnostic; we just track the live
    // constant so a future bump that forgot this round-trip stands out.
    REQUIRE(persist::kFormatVersion == 4u);
    unsigned char blob[8192];
    const Expected<int> wr = rig.estimator().serialize(blob, sizeof(blob));
    REQUIRE(wr.ok());

    Estimator fresh;
    REQUIRE(fresh.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(fresh.add_source(sp.get()) == Status::Ok);
    CHECK_FALSE(fresh.scale2_committed(3));          // cold before restore
    REQUIRE(fresh.deserialize(blob, wr.value()) == Status::Ok);
    CHECK(fresh.scale2_committed(3));                // restored + held through refill

    unsigned char blob2[8192];
    const Expected<int> wr2 = fresh.serialize(blob2, sizeof(blob2));
    REQUIRE(wr2.ok());
    REQUIRE(wr2.value() == wr.value());
    bool same = true;
    for (int i = 0; i < wr.value(); ++i) {
        if (blob[i] != blob2[i]) { same = false; break; }
    }
    CHECK(same);                                     // byte-identical re-serialize

    // v2-STAMP REJECT: the previous format version must reject with VersionMismatch
    // (old blobs cold-start; the established no-migration precedent).
    {
        unsigned char old_blob[8192];
        for (int i = 0; i < wr.value(); ++i) old_blob[i] = blob[i];
        old_blob[4] = 2u;                            // version word (LE) -> 2
        old_blob[5] = old_blob[6] = old_blob[7] = 0u;
        CHECK(fresh.deserialize(old_blob, wr.value()) == Status::VersionMismatch);
    }

    // CONFIG-HASH FLIP: the same rig with joint_lever_scale=false must REJECT the blob
    // written with the flag on (InvalidConfig) — the flag is calibration-shaping.
    {
        Config cfg_off = cfg;
        cfg_off.joint_lever_scale = false;
        auto est_off = std::unique_ptr<Estimator>(new Estimator());
        REQUIRE(est_off->init(cfg_off) == Status::Ok);
        CHECK(est_off->deserialize(blob, wr.value()) == Status::InvalidConfig);
    }
}

// ===========================================================================
// 3b. Either-path fold + cross-reset (the no-double-fold mutation guard): BOTH
// scale estimators fill with the SAME ~1.10 residual before either commits; the
// first rising edge folds ONCE and resets BOTH histograms, so the second path can
// only ever re-commit the post-fold ≈1 residual. Dropping the cross-reset (the
// mutation) leaves the other histogram full of stale 1.10 votes — its later edge
// folds the SAME residual AGAIN (prior_scale ≈ 1.21) and the bound below fails.
// ===========================================================================
TEST_CASE("joint lever+scale: either-path rising edge folds once — no double-fold "
          "across the Phase-1 and scale2 paths (both fold orders)") {
    // Alternating straight/turn so BOTH regimes vote from the start; commit_min_votes
    // high enough that both histograms are FULL of the planted residual before either
    // path's first commit. Two passes flip which regime LEADS each pair, so the first
    // fold comes from Phase-1 in one pass and from scale2 in the other — exercising
    // BOTH cross-resets (dropping either one re-folds the same residual: ~1.21 spike).
    for (int pass = 0; pass < 2; ++pass) {
        const bool straight_first = (pass == 0);
        CAPTURE(straight_first);
        Trajectory tr;
        Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
        Vec6 turn;     turn     << 2.0, 0, 0, 0, 0, 0.6;
        for (int rep = 0; rep < 10; ++rep) {
            if (straight_first) {
                tr.add_segment(straight, 0.8);
                tr.add_segment(turn, 0.8);
            } else {
                tr.add_segment(turn, 0.8);
                tr.add_segment(straight, 0.8);
            }
        }

        const Scalar planted_scale = 1.10;
        std::vector<SourceParams> planted(4);
        for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
        planted[3].scale = planted_scale;       // identity mount — pure scale error
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

        std::vector<SensorConfig> sensors(4);
        for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);

        Config cfg;
        cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
        cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
        set_hists(cfg);
        // Both paths accumulate ~40 votes/segment; 150 votes ≈ 4 segment-pairs of
        // BOTH-regime voting before the first commit is possible.
        cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 150;
        cfg.joint_lever_scale = true;           // rot3d off: planar turn, identity mounts
        cfg.sensors = sensors.data(); cfg.sensor_count = 4;

        Rig rig;
        rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
        const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 500);

        // The committed absolute scale lands on the planted value ONCE...
        const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
        REQUIRE(cs != nullptr);
        INFO("final published scale = " << cs->scale);
        CHECK(cs->scale_committed);
        CHECK(std::abs(cs->scale - planted_scale) < 0.03);

        // ...and NEVER overshoots toward residual² ≈ 1.21 at ANY step (the double-fold
        // signature). The bound sits midway between one fold (1.10) and two (1.21).
        Scalar max_scale = 0;
        for (const auto& rec : rig.records()) {
            if (!rec.fused) continue;
            const CalibSnapshot* c3 = snap(rec.result, 3);
            if (c3 != nullptr) max_scale = std::max(max_scale, c3->scale);
        }
        INFO("max published scale over the run = " << max_scale);
        CHECK(max_scale < Scalar(1.155));
    }
}

// ===========================================================================
// Range guard, estimator level (17b review MAJOR-3 i): the POISONING repro. With a
// WRONG prior rotation the pre-rot3d-gate lever rows are built with the wrong R_X
// and the κ sink absorbs the pollution — the running s_res leaves the histogram
// range ON CLEAN UNIT-SCALE DATA. Those residuals must deposit NOTHING: under the
// "replace skip with clamp" mutation the deterministic transient mass concentrates
// in the edge bin, commits, and folds ~1.5 into prior_scale — an IRREVERSIBLE value
// poison no later reset can heal (the per-step scale bound below fails).
// ===========================================================================
TEST_CASE("joint lever+scale estimator: pre-gate wrong-R_X transient deposits no "
          "edge mass and never commits — prior_scale never poisoned") {
    const Trajectory tr = multiaxis_turnonly_traj();
    // Large planted mount ROTATION with an IDENTITY prior (the rot3d bootstrap shape):
    // every pre-gate row's R_X is badly wrong. UNIT planted scale — any scale the κ
    // sink reads during the transient is pure pollution.
    const SE3 X3 = make_extrinsic(0.25, -0.18, 0.12, Vec3(0.20, -0.15, 0.10));

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X = X3;                                   // unit scale, rotated mount
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    // sensors[3].prior_extrinsic left at identity — the WRONG prior.

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 120;
    cfg.rot3d_enabled     = true;
    cfg.joint_lever_scale = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    // The transient DID produce out-of-range residuals (the repro is live, and the
    // guard — not luck — is what kept them out of the histogram)...
    CHECK(rig.estimator().scale2_skipped(3) > Scalar(0));
    // ...and the transient garbage never looked concentrated enough to COMMIT: on
    // every step where the scale is committed the published value sits at unit (a
    // legitimate late ≈1.0 fold at most). Under the clamp mutation the deterministic
    // edge mass concentrates, commits EARLY, and folds ~1.4–1.5 into prior_scale —
    // the committed-step bound fails. (UNcommitted snapshot values may wobble while
    // the in-range slice of the sweeping transient sits in the histogram: that is the
    // pre-existing "turn residual surfaces at Phase-1 conf 0" snapshot multiply, not
    // a fold — prior_scale is only ever moved by a COMMIT edge.)
    Scalar max_committed_dev = 0;
    for (const auto& rec : rig.records()) {
        if (!rec.fused) continue;
        const CalibSnapshot* c3 = snap(rec.result, 3);
        if (c3 == nullptr || !c3->scale_committed) continue;
        max_committed_dev = std::max(max_committed_dev,
                                     std::abs(c3->scale - Scalar(1)));
    }
    INFO("max |published scale - 1| over committed steps = " << max_committed_dev);
    CHECK(max_committed_dev < Scalar(0.05));
    // Final state: unit scale (prior_scale unpoisoned), mount rotation recovered.
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    INFO("final published scale = " << cs->scale);
    CHECK(std::abs(cs->scale - Scalar(1)) < 0.05);
    CHECK(so3::log(cs->extrinsic.R * X3.R.transpose()).norm() < 0.02);
}

// ===========================================================================
// Range guard, estimator level (17b review MAJOR-3 ii): a LEGITIMATE out-of-range
// true scale (2.0) on a TURN-ONLY rig — the documented limitation, pinned: the turn
// path never votes/commits/folds (prior_scale untouched, published scale stays at
// the prior 1) and Phase-1 cannot help (no straight regime). NOT silent: the
// skipped counter reads high, distinguishing out-of-regime from under-excitation.
// ===========================================================================
TEST_CASE("joint lever+scale estimator: out-of-range TRUE scale on a turn-only rig "
          "never commits and prior_scale stays untouched — but reads as skipped") {
    const Trajectory tr = multiaxis_turnonly_traj();
    const SE3 X3 = make_extrinsic(0.10, -0.07, 0.09, Vec3(0.15, -0.10, 0.05));

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X     = X3;
    planted[3].scale = 2.0;                              // outside (0.5, 1.5)
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = X3;                     // isolate the scale DOF

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 120;
    cfg.rot3d_enabled     = true;
    cfg.joint_lever_scale = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    // Never voted, never committed, never folded — at ANY step.
    CHECK_FALSE(rig.estimator().scale2_committed(3));
    for (const auto& rec : rig.records()) {
        if (!rec.fused) continue;
        const CalibSnapshot* c3 = snap(rec.result, 3);
        if (c3 == nullptr) continue;
        CHECK_FALSE(c3->scale_committed);
        CHECK(c3->scale == doctest::Approx(1.0));        // prior_scale untouched
    }
    // The limitation is DIAGNOSABLE: rows accumulated, votes 0, skipped HIGH.
    CHECK(rig.estimator().scale2_skipped(3) > Scalar(50));
}

// ===========================================================================
// Lever commit across the scale2 fold (17b review, reset_scale2/rows_ coupling):
// prior lever != truth (the shape the existing estimator test sidesteps with
// prior == truth) — the lever must commit, HOLD through the scale2 fold's
// reset_scale2 (which keeps the xyz histograms + their gate vote mass), and both
// DOFs land on truth. Pins the previously comment-only "committed lever held
// through the fold" claim.
// ===========================================================================
TEST_CASE("joint lever+scale estimator: lever commit holds across the scale2 fold "
          "(prior lever wrong; no commit thrash)") {
    const Trajectory tr = multiaxis_turnonly_traj();
    const SE3 X3 = make_extrinsic(0.10, -0.07, 0.09, Vec3(0.15, -0.10, 0.05));
    const Scalar planted_scale = 1.08;

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X     = X3;
    planted[3].scale = planted_scale;
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic.R = X3.R;                 // rotation prior = truth,
    sensors[3].prior_extrinsic.t = Vec3::Zero();         // lever prior WRONG (zero)

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 120;
    cfg.rot3d_enabled     = true;
    cfg.joint_lever_scale = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    // Both DOFs committed and landed on truth...
    CHECK(rig.estimator().scale2_committed(3));
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    CHECK(cs->translation_committed);
    CHECK(cs->scale_committed);
    INFO("published scale = " << cs->scale
         << ", lever err = " << (cs->extrinsic.t - X3.t).norm() << " m");
    CHECK(std::abs(cs->scale - planted_scale) < 0.02);
    CHECK((cs->extrinsic.t - X3.t).norm() < 0.02);

    // ...and the lever commit NEVER thrashed: once up, it stayed up — in particular
    // across the scale2 fold's reset_scale2 (which must keep the gate's vote mass).
    bool lever_seen = false, scale_seen = false;
    for (const auto& rec : rig.records()) {
        if (!rec.fused) continue;
        const CalibSnapshot* c3 = snap(rec.result, 3);
        if (c3 == nullptr) continue;
        if (c3->translation_committed) lever_seen = true;
        else CHECK_FALSE(lever_seen);                    // no un-commit, ever
        if (c3->scale_committed) scale_seen = true;
    }
    CHECK(lever_seen);
    CHECK(scale_seen);                                   // the fold DID happen
}

// ===========================================================================
// Premature fold (17b review MINOR): with a LOW commit_min_votes the scale2 commit
// can fire while the running κ̂ convergence transient still owns the SlidingK mode —
// the one-shot rising-edge fold then moves prior_scale by a PARTIAL residual. Pinned
// benign-by-design behavior: the latch HOLDS (no re-edge, no second fold), the
// post-fold votes record the REMAINING residual vs the new prior, and the PUBLISHED
// scale (prior × live residual) still converges to the planted truth.
// ===========================================================================
TEST_CASE("joint lever+scale estimator: premature fold at low commit_min_votes — "
          "latch holds and the published prior x residual still reaches truth") {
    const Trajectory tr = multiaxis_turnonly_traj();
    const SE3 X3 = make_extrinsic(0.10, -0.07, 0.09, Vec3(0.15, -0.10, 0.05));
    const Scalar planted_scale = 1.08;

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X     = X3;
    planted[3].scale = planted_scale;
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = X3;                     // isolate the scale DOF

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    // LOW N_min — the production-config shape the 120-vote main test tunes around.
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 30;
    cfg.rot3d_enabled     = true;
    cfg.joint_lever_scale = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    // Latch: exactly one rising edge, never re-opened (a re-edge would re-fold and
    // overshoot — the over-fold bound below also guards it).
    CHECK(rig.estimator().scale2_committed(3));
    bool seen_commit = false;
    Scalar max_scale = 0;
    for (const auto& rec : rig.records()) {
        if (!rec.fused) continue;
        const CalibSnapshot* c3 = snap(rec.result, 3);
        if (c3 == nullptr) continue;
        if (c3->scale_committed) seen_commit = true;
        else CHECK_FALSE(seen_commit);                   // latch never drops
        max_scale = std::max(max_scale, c3->scale);
    }
    CHECK(seen_commit);
    CHECK(max_scale < planted_scale * Scalar(1.06));     // no over-fold past truth

    // Whether or not the fold was partial, the SNAPSHOT scale is prior × the LIVE
    // residual — the leftover surfaces multiplicatively and lands on the truth.
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    INFO("final published scale = " << cs->scale);
    CHECK(std::abs(cs->scale - planted_scale) < 0.02);
}
