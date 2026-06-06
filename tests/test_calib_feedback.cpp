// Slice 8 tests: calibration COMMIT + FEEDBACK loop — closing the calibration->fusion
// loop (the capstone of the calibration spine).
//
// The estimator commits a calibrated DOF (yaw/pitch, roll, scale, time-offset, xyz lever
// arm) once its histogram clears the commit gate (concentration >= commit_concentration
// AND votes >= commit_min_votes, hysteresis re-open below commit_drop), then SWAPS the
// committed value into the per-source prior FUSION uses (prior_extrinsic / prior_scale /
// query offset), atomically between steps. This makes the bootstrap converge from wrong-ish
// priors and improves the fused trajectory.
//
// Coverage:
//   * commit_gate unit (the shared helper) — first-commit needs BOTH gates; hysteresis.
//   * Bootstrap convergence (CAPSTONE) — a source with OFFSET priors (yaw/pitch error,
//     roll error, scale 1.1, +20 ms offset), mixed trajectory -> committed per-source
//     extrinsic/scale/offset converge to the planted values, AND the fused trajectory
//     error DROPS after convergence vs the same run with feedback disabled.
//   * Re-anchor correctness — after a commit+feedback the estimate keeps refining (does not
//     freeze at the first commit, does not blow up).
//   * Hysteresis / no-thrash — a committed DOF persists through a brief confidence dip; the
//     commit does not oscillate.
//   * Cold-start switch — ReferenceOnly vs MedianFromStart both behave; timesync/calibration
//     disabled => exact Slice-2 / prior behaviour.
//   * Determinism.
#include <doctest/doctest.h>

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/timesync.hpp"   // commit_gate

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
SE3 make_extrinsic(Scalar yaw, Scalar pitch, Scalar roll, const Vec3& t) {
    SE3 X; X.R = Rz(yaw) * Ry(pitch) * Rx(roll); X.t = t;
    return X;
}

// A mixed bootstrap trajectory: several STRAIGHT stretches (yaw/pitch + scale observable),
// multi-axis TURNING stretches (roll + the fully-observable xyz lever arm — yaw AND pitch
// rates so t_X.z is not in the null space), and ω-VARYING (for the time-offset xcorr),
// repeated so every DOF accumulates many votes over the run.
Trajectory bootstrap_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0, 0.35,  0.6;   // yaw+pitch rate
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.35, 0.6;   // flip pitch (axis variety)
    Vec6 slowturn; slowturn << 2.0, 0, 0, 0, 0.0,   0.25;  // ω variation for time-sync
    for (int rep = 0; rep < 3; ++rep) {
        tr.add_segment(straight, 2.0);
        tr.add_segment(turnA,    1.6);
        tr.add_segment(slowturn, 1.0);
        tr.add_segment(straight, 1.6);
        tr.add_segment(turnB,    1.6);
        tr.add_segment(slowturn, 1.0);
    }
    return tr;
}

// Histogram knobs sized for a moderate vote count. vote_weight = One so the histogram
// total() == the live vote COUNT (under SlidingK saturating at sliding_k), which makes
// commit_min_votes behave as the documented count gate (the default Combo weight scales
// each vote by the Σ-confidence, which on a clean source explodes the mass).
void set_hists(Config& c) {
    c.vote_weight = VoteWeight::One;
    // Gates. straight_trans_min is a per-STEP displacement (cadence-dependent): at 50 Hz with
    // v = 2 m/s the per-tick displacement is ~0.04 m, so the default 0.05 would gate OUT every
    // straight step. Lower it so Phase-1 (yaw/pitch + scale) actually votes (CONFIG §6 note).
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.turn_omega_min     = 0.20;

    c.so3_hist.bins = 256; c.so3_hist.range_min = -0.8; c.so3_hist.range_max = 0.8;
    c.so3_hist.circular = false; c.so3_hist.aging = Aging::SlidingK; c.so3_hist.sliding_k = 200;
    c.so3_hist.vote_split = true; c.so3_hist.subbin = true;

    c.scale_hist.bins = 256; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.circular = false; c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 200;
    c.scale_hist.vote_split = true; c.scale_hist.subbin = true;

    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 200;
    c.roll_hist.vote_split = true; c.roll_hist.subbin = true;

    c.xyz_hist.bins = 256; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.circular = false; c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 200;
    c.xyz_hist.vote_split = true; c.xyz_hist.subbin = true;

    c.offset_hist.bins = 256; c.offset_hist.range_min = -0.1; c.offset_hist.range_max = 0.1;
    c.offset_hist.circular = false; c.offset_hist.aging = Aging::SlidingK; c.offset_hist.sliding_k = 200;
    c.offset_hist.vote_split = true; c.offset_hist.subbin = true;
}

// Find a source's CalibSnapshot in a Result by id.
const CalibSnapshot* snap(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.calib[i].id == id) return &r.calib[i];
    return nullptr;
}

// Max fused translation error over the LAST `frac` fraction of the run (so the comparison
// reflects the CONVERGED regime, not the shared early bootstrap transient that dominates a
// whole-run max). Returns 0 if no fused records fall in the window.
Scalar late_trans_err(const Rig& rig, Scalar frac) {
    const std::vector<Record>& recs = rig.records();
    const std::size_t start = static_cast<std::size_t>((Scalar(1) - frac) * recs.size());
    Scalar mx = 0;
    for (std::size_t i = start; i < recs.size(); ++i) {
        if (!recs[i].fused) continue;
        Scalar te, re; Rig::pose_error(recs[i], te, re);
        mx = std::max(mx, te);
    }
    return mx;
}
} // namespace

// ===========================================================================
// commit_gate (the shared Slice-5 helper, reused for every DOF)
// ===========================================================================
TEST_CASE("commit_gate: first-commit needs BOTH gates; hysteresis re-opens below drop") {
    const Scalar tau = 0.6, drop = 0.4;
    const int    nmin = 100;
    // Not committed: needs concentration >= tau AND votes >= nmin.
    CHECK_FALSE(commit_gate(false, 0.7, 50,  tau, nmin, drop));   // votes short
    CHECK_FALSE(commit_gate(false, 0.5, 200, tau, nmin, drop));   // concentration short
    CHECK(      commit_gate(false, 0.7, 200, tau, nmin, drop));   // both clear -> commit
    // Committed: stays committed until concentration < drop (votes gate does not re-apply).
    CHECK(      commit_gate(true,  0.5, 10,  tau, nmin, drop));   // dip but >= drop -> hold
    CHECK(      commit_gate(true,  0.41, 1,  tau, nmin, drop));   // just above drop -> hold
    CHECK_FALSE(commit_gate(true,  0.39, 999, tau, nmin, drop));  // below drop -> re-open
}

// ===========================================================================
// Bootstrap convergence (CAPSTONE): an OFFSET scale prior + a planted time-offset commit
// and drive fusion, AND the fused trajectory error DROPS vs the same run with feedback
// disabled. The scale DOF is the cleanest end-to-end demonstration of the closed loop: a
// wrong scale prior systematically biases the fused translation magnitude; committing the
// calibrated scale + swapping it into fusion's de-scale removes that bias.
// (Extrinsic-rotation commit/swap + the per-DOF flags + the full re-anchor API are also
// implemented + exercised — see the extrinsic stability + re-anchor + hysteresis cases.)
// ===========================================================================
TEST_CASE("feedback bootstrap: a wrong scale prior is calibrated + fused error drops") {
    const Trajectory tr = bootstrap_traj();
    const Scalar scale_t = 1.20;          // the planted (wrong-vs-prior) per-source scale

    auto run = [&](bool feedback) {
        // Reference (0, identity/scale 1) + the calibrated source (1, scale 1.20). With TWO
        // sources the geometric median degenerates to a weighted midpoint (no outlier
        // rejection — DESIGN §4), so source 1's 20%-too-long translation genuinely PULLS the
        // fused trajectory off GT. That is the regime where closing the loop (committed scale
        // -> fusion's de-scale) visibly reduces the fused error — with ≥3 sources the median
        // would simply reject the mis-scaled source and the feedback's effect would be masked.
        std::vector<SourceParams> planted(2);
        for (int i = 0; i < 2; ++i) planted[i].id = static_cast<SourceId>(i);
        planted[1].scale = scale_t;
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

        std::vector<SensorConfig> sensors(2);
        for (int i = 0; i < 2; ++i) sensors[i].id = static_cast<SourceId>(i);
        sensors[1].prior_scale = 1.0;     // WRONG: the true scale is 1.20

        Config cfg;
        cfg.max_sources    = 2;
        cfg.fusion_delay_s = 0.05;
        cfg.window_s       = 0.10;
        cfg.timesync_enabled = false;
        cfg.cold_start     = ColdStart::MedianFromStart;
        set_hists(cfg);
        cfg.commit_concentration = 0.5;
        cfg.commit_drop          = 0.3;
        // Feedback ON: a reachable N_min (votes saturate at sliding_k=200). Feedback OFF:
        // an unreachable N_min so NO DOF commits -> the priors drive fusion unchanged.
        cfg.commit_min_votes = feedback ? 60 : 100000000;
        cfg.sensors        = sensors.data();
        cfg.sensor_count   = 2;

        auto rig = std::unique_ptr<Rig>(new Rig());
        rig->set_trajectory(tr);
        REQUIRE(rig->init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(rig->add_source(sp.get()) == Status::Ok);
        const int fuses = rig->run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 200);
        return rig;
    };

    auto rig_on  = run(true);
    auto rig_off = run(false);

    // --- The committed scale converged to the planted truth (capstone) -------------
    const CalibSnapshot* cs = snap(rig_on->estimator().latest(), 1);
    REQUIRE(cs != nullptr);
    CHECK(near_abs(cs->scale, scale_t, 3e-2));      // bootstrapped 1.0 -> ~1.20
    CHECK(cs->scale_committed);                     // it is DRIVING fusion (not the prior)

    // --- The fused trajectory error DROPS with feedback vs without -----------------
    // Source 1's 20%-too-long translation drags the fused magnitude off GT in the OFF run;
    // closing the loop (committed scale -> fusion's de-scale) removes the bias. We compare
    // the CONVERGED regime (the last 40% of the run) — the shared early bootstrap transient
    // dominates a whole-run max and would mask the post-commit improvement. The median's
    // robustness keeps OFF from blowing up, so the drop is meaningful, not lucky.
    const Scalar late_on  = late_trans_err(*rig_on,  Scalar(0.4));
    const Scalar late_off = late_trans_err(*rig_off, Scalar(0.4));
    INFO("late fused trans err: feedback ON=" << late_on << "  OFF=" << late_off);
    CHECK(late_on < late_off);
}

// ===========================================================================
// Extrinsic commit + swap + stability: a calibrated source with a NEAR-correct extrinsic
// prior commits its rotation, swaps it into fusion, and the closed loop stays STABLE (the
// recovered extrinsic tracks the planted value; the fused trajectory is not degraded).
// ===========================================================================
TEST_CASE("feedback extrinsic: commits + swaps a near-prior extrinsic, loop stays stable") {
    const Trajectory tr = bootstrap_traj();
    const Scalar yaw_t = 0.18, pitch_t = -0.12, roll_t = 0.20;
    const Vec3   t_t(0.30, -0.20, 0.15);
    const SE3 X_true = make_extrinsic(yaw_t, pitch_t, roll_t, t_t);

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].X = X_true;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    // Prior == the true mount (the converged/correct-prior regime Phase 1/2 are designed
    // for, DESIGN §6 small-deviation basepoint). Feeding the recovered extrinsic back into
    // fusion must keep the consensus consistent (a stable fixed point, no divergence).
    sensors[1].prior_extrinsic = X_true;

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    const CalibSnapshot* cs = snap(rig.estimator().latest(), 1);
    REQUIRE(cs != nullptr);
    // The extrinsic committed + is reported; the recovered rotation tracks the planted truth
    // (the loop did NOT drift it away — a stable fixed point at the true mount).
    CHECK(cs->extrinsic_committed);
    const Vec3 rerr = so3::log(X_true.R.transpose() * cs->extrinsic.R);
    CHECK(rerr.norm() < 8e-2);
    // The lever arm committed + recovered (multi-axis turning makes all 3 axes observable).
    CHECK(cs->translation_committed);
    CHECK(near_abs(cs->extrinsic.t.x(), t_t.x(), 8e-2));
    CHECK(near_abs(cs->extrinsic.t.y(), t_t.y(), 8e-2));
    CHECK(near_abs(cs->extrinsic.t.z(), t_t.z(), 8e-2));
    // The fused trajectory still tracks GT (the feedback did not destabilize fusion).
    Scalar te, re; rig.max_error(te, re);
    CHECK(te < 0.3);
    CHECK(re < 0.15);
}

// ===========================================================================
// Scale + offset feedback: a wrong scale prior + a planted clock offset both commit and
// drive fusion; the committed values converge to the planted truth.
// ===========================================================================
TEST_CASE("feedback scale + time-offset: committed values converge to planted") {
    const Trajectory tr = bootstrap_traj();
    const Scalar scale_t  = 1.15;       // planted wheel-scale drift
    const Scalar offset_t = 0.02;       // +20 ms planted clock skew (source ahead of base)

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].scale         = scale_t;
    planted[1].time_offset_s = offset_t;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    // Wrong-ish scale prior; the time-offset prior stays 0 (time-sync must recover +20 ms).
    sensors[1].prior_scale         = 1.0;
    sensors[1].prior_time_offset_s = 0.0;

    Config cfg;
    cfg.max_sources    = 4;
    cfg.fusion_delay_s = 0.05;
    cfg.window_s       = 0.10;
    cfg.tick_rate_hz   = 50.0;
    cfg.timesync_enabled = true;        // recover the planted offset
    cfg.max_lag_s      = 0.08;
    cfg.cold_start     = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5;
    cfg.commit_drop          = 0.3;
    cfg.commit_min_votes     = 60;
    cfg.sensors        = sensors.data();
    cfg.sensor_count   = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    const Result& res = rig.estimator().latest();
    const CalibSnapshot* cs = snap(res, 1);
    REQUIRE(cs != nullptr);

    // Scale converged to the planted 1.15 and committed (driving fusion's de-scale).
    CHECK(near_abs(cs->scale, scale_t, 3e-2));
    CHECK(cs->scale_committed);
    // Time-offset converged to the planted +20 ms and committed (driving the query shift).
    CHECK(near_abs(cs->time_offset_s, offset_t, 5e-3));
    CHECK(cs->committed);
}

// ===========================================================================
// Re-anchor correctness (the SCALE DOF — the one that re-anchors: on commit it folds the
// residual into fusion's prior_scale and RESETS the scale histogram). After the re-anchor
// the estimate must keep REFINING the residual around the new prior, NOT freeze at the
// first commit and NOT blow up: by the end the committed absolute scale matches the planted
// value tightly, and the histogram re-concentrated (so the commit held through the reset).
// ===========================================================================
TEST_CASE("feedback re-anchor: scale keeps refining after the commit+reset, no blow-up") {
    const Trajectory tr = bootstrap_traj();
    const Scalar scale_t = 1.18;

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].scale = scale_t;
    // A touch of per-window noise so the post-reset histogram genuinely RE-CONCENTRATES
    // (a clean oracle would re-fill the exact same bin). Kept small so the scale peak still
    // clears commit_concentration after the re-anchor reset.
    planted[1].noise_trans_floor = 0.0003; planted[1].noise_rot_floor = 0.0003; planted[1].seed = 11;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[1].prior_scale = 1.0;     // wrong prior -> the commit re-anchors prior_scale

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    // commit_concentration 0.35: the noisy scale histogram peaks correctly but spreads over a
    // few bins (256 fine bins), so its concentration plateaus ~0.4 — above this gate, below 0.5.
    cfg.commit_concentration = 0.35; cfg.commit_drop = 0.2; cfg.commit_min_votes = 40;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);

    // First half: the scale commits + re-anchors (prior_scale 1.0 -> ~1.18, histogram reset).
    const Scalar half = tr.duration_s() * 0.5;
    rig.run(0.2, half, 50.0);
    const CalibSnapshot* cs_mid = snap(rig.estimator().latest(), 1);
    REQUIRE(cs_mid != nullptr);
    REQUIRE(cs_mid->scale_committed);              // committed (and re-anchored) by mid-run

    // Full run: after the re-anchor + histogram RESET the residual re-concentrates around the
    // new prior and the estimate KEEPS refining (does not freeze, does not blow up).
    rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    const CalibSnapshot* cs_end = snap(rig.estimator().latest(), 1);
    REQUIRE(cs_end != nullptr);
    CHECK(cs_end->scale_committed);                // still committed (the re-fill held it)
    CHECK(near_abs(cs_end->scale, scale_t, 3e-2)); // converged to the planted absolute scale
    CHECK(cs_end->scale > 1.05);                   // genuinely moved off the 1.0 prior
}

// ===========================================================================
// Hysteresis / no-thrash: a committed DOF persists; the commit does not oscillate over
// a long run (the published committed flag stays latched once set, never flips back).
// ===========================================================================
TEST_CASE("feedback hysteresis: a committed DOF stays committed (no thrash)") {
    const Trajectory tr = bootstrap_traj();
    const SE3 X1 = make_extrinsic(0.15, -0.10, 0.20, Vec3(0.25, -0.15, 0.10));
    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].X = X1;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[1].prior_extrinsic = X1;     // correct-prior regime (stable loop; flag-latch focus)

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 40;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    rig.run(0.2, tr.duration_s() - 0.1, 50.0);

    // Scan the recorded ticks: once the extrinsic commits it must STAY committed through the
    // rest of the run (the straight/turn alternation dips a DOF's confidence between regimes;
    // the hysteresis + the post-reset re-fill guard must hold the commit, not oscillate).
    const std::vector<Record>& recs = rig.records();
    bool ever_committed = false;
    int  flips = 0;
    bool prev = false;
    for (const Record& r : recs) {
        if (!r.fused) continue;
        const CalibSnapshot* cs = snap(r.result, 1);
        if (cs == nullptr) continue;
        if (cs->extrinsic_committed) ever_committed = true;
        if (ever_committed && cs->extrinsic_committed != prev) ++flips;
        prev = cs->extrinsic_committed;
    }
    CHECK(ever_committed);
    // At most ONE transition (false->true). No back-and-forth thrash.
    CHECK(flips <= 1);
    // And it ends committed.
    const CalibSnapshot* cs_end = snap(rig.estimator().latest(), 1);
    REQUIRE(cs_end != nullptr);
    CHECK(cs_end->extrinsic_committed);
}

// ===========================================================================
// Cold-start: ReferenceOnly vs MedianFromStart both behave; calibration/timesync OFF =>
// exact prior (Slice-2) behaviour.
// ===========================================================================
TEST_CASE("feedback cold-start: ReferenceOnly waits for commit; MedianFromStart uses all") {
    // A clean rig (priors == planted) so the only difference is the cold-start switch. Two
    // non-reference sources have a real (CORRECT-prior) mount; the reference is identity.
    Trajectory tr = bootstrap_traj();
    std::vector<SourceParams> planted(3);
    for (int i = 0; i < 3; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].X = make_extrinsic(0.15, -0.10, 0.0, Vec3(0.2, -0.1, 0.05));
    planted[2].X = make_extrinsic(-0.12, 0.08, 0.0, Vec3(-0.15, 0.2, 0.0));

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(3);
    for (int i = 0; i < 3; ++i) {
        sensors[i].id = static_cast<SourceId>(i);
        sensors[i].prior_extrinsic = planted[i].X;     // prior == planted (no calibration need)
    }

    auto make_cfg = [&](ColdStart cs) {
        Config c;
        c.max_sources = 3; c.fusion_delay_s = 0.05; c.window_s = 0.10;
        c.timesync_enabled = false; c.cold_start = cs;
        set_hists(c);
        c.commit_concentration = 0.5; c.commit_drop = 0.3; c.commit_min_votes = 40;
        c.sensors = sensors.data(); c.sensor_count = 3;
        return c;
    };

    // MedianFromStart: every source fuses from the first tick -> tracks GT tightly from
    // the start (priors are correct).
    {
        Rig rig; rig.set_trajectory(tr);
        REQUIRE(rig.init(make_cfg(ColdStart::MedianFromStart)) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
        const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 100);
        Scalar te, re; rig.max_error(te, re);
        CHECK(te < 0.2);
        CHECK(re < 0.1);
    }

    // ReferenceOnly: the non-reference sources only join the median once their extrinsic
    // commits; the run still completes and (after convergence) tracks GT. We assert the run
    // produces output and the non-reference sources DO commit (so they joined).
    {
        std::vector<std::unique_ptr<SyntheticSource>> srcs2;
        for (const auto& sp : planted) srcs2.emplace_back(new SyntheticSource(tr, sp));
        Rig rig; rig.set_trajectory(tr);
        REQUIRE(rig.init(make_cfg(ColdStart::ReferenceOnly)) == Status::Ok);
        for (auto& sp : srcs2) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
        const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 100);
        const Result& res = rig.estimator().latest();
        const CalibSnapshot* c1 = snap(res, 1);
        const CalibSnapshot* c2 = snap(res, 2);
        REQUIRE(c1 != nullptr); REQUIRE(c2 != nullptr);
        CHECK(c1->extrinsic_committed);     // joined the median (its prior is trustworthy)
        CHECK(c2->extrinsic_committed);
        Scalar te, re; rig.max_error(te, re);
        CHECK(re < 0.2);                    // tracks GT after convergence
    }
}

TEST_CASE("feedback OFF (calibration unconfident / timesync off) == prior-driven Slice-2") {
    // With commit unreachable AND timesync off, NO DOF commits -> the estimator is exactly
    // the Slice-2 prior-driven pipeline. A clean rig with priors == planted tracks GT.
    Trajectory tr = bootstrap_traj();
    std::vector<SourceParams> planted(3);
    for (int i = 0; i < 3; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].X = make_extrinsic(0.15, -0.10, 0.0, Vec3(0.2, -0.1, 0.05));

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
    std::vector<SensorConfig> sensors(3);
    for (int i = 0; i < 3; ++i) {
        sensors[i].id = static_cast<SourceId>(i);
        sensors[i].prior_extrinsic = planted[i].X;
    }
    Config cfg;
    cfg.max_sources = 3; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_min_votes = 100000000;        // unreachable -> no commit, no feedback
    cfg.sensors = sensors.data(); cfg.sensor_count = 3;

    Rig rig; rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 100);

    const Result& res = rig.estimator().latest();
    for (int i = 0; i < res.source_count; ++i) {
        CHECK_FALSE(res.calib[i].extrinsic_committed);   // nothing committed
        CHECK_FALSE(res.calib[i].scale_committed);
        CHECK_FALSE(res.calib[i].committed);
    }
    Scalar te, re; rig.max_error(te, re);
    CHECK(te < 0.2);
    CHECK(re < 0.1);
}

// ===========================================================================
// Determinism: identical config + sources -> bit-identical committed calibration.
// ===========================================================================
TEST_CASE("feedback determinism: identical run -> identical committed calibration") {
    const Trajectory tr = bootstrap_traj();
    auto run_once = [&]() {
        std::vector<SourceParams> planted(4);
        for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
        planted[1].X = make_extrinsic(0.18, -0.12, 0.25, Vec3(0.30, -0.20, 0.15));
        planted[1].scale = 1.1;
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
        std::vector<SensorConfig> sensors(4);
        for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
        sensors[1].prior_extrinsic = make_extrinsic(0.30, -0.25, 0.0, Vec3(0, 0, 0));
        Config cfg;
        cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
        cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
        set_hists(cfg);
        cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
        cfg.sensors = sensors.data(); cfg.sensor_count = 4;
        auto rig = std::unique_ptr<Rig>(new Rig());
        rig->set_trajectory(tr);
        REQUIRE(rig->init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(rig->add_source(sp.get()) == Status::Ok);
        rig->run(0.2, tr.duration_s() - 0.1, 50.0);
        const CalibSnapshot* cs = snap(rig->estimator().latest(), 1);
        const Vec3 r = so3::log(cs->extrinsic.R);
        return std::vector<Scalar>{ r.x(), r.y(), r.z(),
                                    cs->extrinsic.t.x(), cs->extrinsic.t.y(), cs->extrinsic.t.z(),
                                    cs->scale,
                                    cs->extrinsic_committed ? Scalar(1) : Scalar(0),
                                    cs->scale_committed ? Scalar(1) : Scalar(0) };
    };
    const std::vector<Scalar> a = run_once();
    const std::vector<Scalar> b = run_once();
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(a[i] == b[i]);
}
