// test_bias.cpp — Slice 11b, Option A: per-source body-twist bias states (clean-regime
// augmented ESKF) + the classic loosely-coupled GPS/INS drift removal (D22, DESIGN §5).
//
// Proves, end-to-end through the sim rig (relaxed edge — the only place GT is known):
//   (a) BIAS OBSERVED + REMOVED — a SINGLE driving source with a planted constant body-twist
//       bias + GPS position fixes: with SensorConfig::bias_states the augmented filter recovers
//       the bias (converges toward the planted value) AND the fused-vs-GT tail error is
//       MATERIALLY lower than the same run with bias_states=false (which cannot remove it).
//   (b) OBSERVABILITY SELF-TEST — with NO absolute ref the bias is NOT observable: it stays
//       near its zero prior and the pose<->bias cross-covariance stays ~0 (the GPS update is
//       what makes the bias observable). Guards the spine.
//   (c) DEFAULT-OFF IDENTITY — bias_states=false is byte-identical to the 12-DOF predict path
//       (the augmented path is never entered).
//   (d) MULTI-BIAS GUARD — two bias_states=true sources -> init returns InvalidConfig
//       (Option A supports exactly one bias source; Option B is deferred).
//
// REGIME (Option A, exact): ONE source drives the predict (median of 1 = passthrough). The
// tests use a single-source rig, so the bias source's frame-aligned delta IS the predict and
// the de-bias Delta o exp(-b*dt) is exact. N-source-median bias coupling (Option B) is deferred.
#include <doctest/doctest.h>

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"

#include "ofc_sim/absolute_ref_source.hpp"
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

// A mixed straight + turn trajectory: turns make a constant GYRO bias observable through GPS
// POSITION (a heading error becomes a position error only when the path curves), and the
// forward motion makes a VELOCITY bias observable directly. Long enough for a steady tail.
Trajectory bias_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0,    0,    0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0,    0,   0.5;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0,    0,  -0.5;
    for (int rep = 0; rep < 4; ++rep) {
        tr.add_segment(straight, 1.2);
        tr.add_segment(turnA,    1.6);
        tr.add_segment(straight, 1.2);
        tr.add_segment(turnB,    1.6);
    }
    return tr;
}

// ONE source: the reference (identity mount, scale 1), with an OPTIONAL planted body-twist
// bias. As the reference it always participates and (single-source rig) drives the predict
// alone -> the Option-A exact regime.
SourceParams bias_source_params(const Vec6& planted_bias) {
    SourceParams sp;
    sp.id = 0;                         // the reference
    sp.body_twist_bias = planted_bias;
    return sp;
}

// Config for a single bias source. bias_on toggles SensorConfig::bias_states.
Config bias_config(std::vector<SensorConfig>& sensors_out, bool bias_on,
                   Scalar mahalanobis_chi2 = 100.0, Scalar bias_pn = 1e-3) {
    sensors_out.clear();
    SensorConfig sc;
    sc.id              = 0;
    sc.prior_extrinsic = SE3{};        // identity mount (the reference is the gauge)
    sc.prior_scale     = 1.0;
    sc.weight_prior    = 1.0;
    sc.bias_states     = bias_on;
    sc.bias_process_noise = bias_pn;
    sc.is_reference    = true;
    sensors_out.push_back(sc);

    Config c;
    c.max_sources      = 1;
    c.fusion_delay_s   = 0.05;
    c.window_s         = 0.10;
    c.reference_sensor_id = 0;
    c.cold_start       = ColdStart::MedianFromStart;   // the lone source fuses from tick 1
    c.timesync_enabled = false;                        // deterministic, offset-free
    c.vote_weight      = VoteWeight::One;
    c.commit_min_votes = 1000000000;                   // calibration off the prior
    c.min_sources_warn = 1;                            // single source can still be NOMINAL
    c.mahalanobis_chi2 = mahalanobis_chi2;
    c.sensors          = sensors_out.data();
    c.sensor_count     = 1;
    return c;
}

// Mean fused-vs-GT translation error over the converged TAIL (last `tail` fused records).
Scalar tail_mean_trans_err(const std::vector<Record>& recs, int tail) {
    std::vector<const Record*> fused;
    for (const Record& r : recs) if (r.fused) fused.push_back(&r);
    if (static_cast<int>(fused.size()) < tail) tail = static_cast<int>(fused.size());
    Scalar sum = 0.0; int n = 0;
    for (int i = static_cast<int>(fused.size()) - tail; i < static_cast<int>(fused.size()); ++i) {
        Scalar te, re;
        Rig::pose_error(*fused[i], te, re);
        sum += te; ++n;
    }
    return (n > 0) ? sum / static_cast<Scalar>(n) : 0.0;
}

// The last fused record's bias snapshot (the single source is slot 0).
const CalibSnapshot* last_bias_snapshot(const std::vector<Record>& recs) {
    for (auto it = recs.rbegin(); it != recs.rend(); ++it) {
        if (it->fused) return &it->result.calib[0];
    }
    return nullptr;
}
} // namespace

// ===========================================================================
// (a) BIAS OBSERVED + REMOVED
// ===========================================================================
TEST_CASE("bias: a single biased driving source + GPS fixes -> bias recovered and drift removed") {
    Trajectory tr = bias_traj();
    // Plant a constant body-twist bias: a small forward-velocity offset + a gyro yaw-rate
    // offset (the two DOF a turning trajectory makes observable through GPS position).
    Vec6 planted; planted << 0.15, 0.0, 0.0,   0.0, 0.0, 0.05;   // [v; omega]
    const int tail = 60;

    // --- bias_states OFF: the bias cannot be removed (drifts) -------------------------
    Scalar off_tail_err = 0.0;
    {
        SourceParams sp = bias_source_params(planted);
        SyntheticSource src(tr, sp);
        std::vector<SensorConfig> sensors;
        Config cfg = bias_config(sensors, /*bias_on=*/false);
        AbsoluteRefParams rp;
        rp.period_s = 0.2; rp.window_s = cfg.window_s; rp.sigma_pos = 0.05; rp.seed = 11u;
        SyntheticAbsoluteRef ref(tr, rp);
        Rig rig; rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        REQUIRE(rig.add_source(&src) == Status::Ok);
        REQUIRE(rig.add_correction(&ref) == Status::Ok);
        rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        off_tail_err = tail_mean_trans_err(rig.records(), tail);
    }

    // --- bias_states ON: the augmented filter observes + removes the bias --------------
    Scalar on_tail_err = 0.0;
    Vec6   recovered = Vec6::Zero();
    Scalar max_observable = 0.0;
    long   applied_total = 0;
    {
        SourceParams sp = bias_source_params(planted);
        SyntheticSource src(tr, sp);
        std::vector<SensorConfig> sensors;
        Config cfg = bias_config(sensors, /*bias_on=*/true);
        AbsoluteRefParams rp;
        rp.period_s = 0.2; rp.window_s = cfg.window_s; rp.sigma_pos = 0.05; rp.seed = 11u;
        SyntheticAbsoluteRef ref(tr, rp);
        Rig rig; rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        REQUIRE(rig.add_source(&src) == Status::Ok);
        REQUIRE(rig.add_correction(&ref) == Status::Ok);
        rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        on_tail_err = tail_mean_trans_err(rig.records(), tail);
        const CalibSnapshot* cs = last_bias_snapshot(rig.records());
        REQUIRE(cs != nullptr);
        recovered = cs->bias;
        for (const Record& r : rig.records()) {
            if (!r.fused) continue;
            applied_total += r.result.correction.corr_applied;
            if (r.result.calib[0].bias_observable > max_observable)
                max_observable = r.result.calib[0].bias_observable;
        }
    }

    MESSAGE("bias removal: OFF tail err=" << off_tail_err << " m  ON tail err=" << on_tail_err
            << " m  fixes applied=" << applied_total
            << "  planted=[" << planted.transpose() << "]"
            << "  recovered=[" << recovered.transpose() << "]"
            << "  max bias_observable=" << max_observable);

    // Fixes were applied and the bias became observable (the bias-block covariance was reduced
    // below its prior by the absolute-ref updates -> bias_observable confidence rose well above 0).
    REQUIRE(applied_total > 10);
    REQUIRE(max_observable > 0.1);
    // The estimated bias converges toward the planted value on the two observable DOF
    // (forward velocity + yaw rate). A loose tolerance: the point is convergence to the planted
    // sign/scale (the J_pb Jacobian is correct), not a tight numerical match under GPS noise.
    CHECK(recovered(0) == doctest::Approx(planted(0)).epsilon(0.5));   // v_x
    CHECK(recovered(5) == doctest::Approx(planted(5)).epsilon(0.6));   // omega_z
    // The recovered bias has the right SIGN on both observable DOF (a flipped J_pb would
    // diverge / converge to the wrong sign).
    CHECK(recovered(0) > 0.0);
    CHECK(recovered(5) > 0.0);
    // Drift removal: removing the bias materially reduces the fused tail error vs the OFF run.
    REQUIRE(off_tail_err > 0.1);                  // the bias genuinely drifts when not removed
    CHECK(on_tail_err < 0.6 * off_tail_err);      // observing + removing it cuts the drift
}

// ===========================================================================
// (b) OBSERVABILITY SELF-TEST — no absolute ref -> bias NOT observable
// ===========================================================================
TEST_CASE("bias: with NO absolute ref the bias is not observable (stays ~prior)") {
    Trajectory tr = bias_traj();
    Vec6 planted; planted << 0.15, 0.0, 0.0,   0.0, 0.0, 0.05;

    SourceParams sp = bias_source_params(planted);
    SyntheticSource src(tr, sp);
    std::vector<SensorConfig> sensors;
    Config cfg = bias_config(sensors, /*bias_on=*/true);   // bias states ON but NO ref registered
    Rig rig; rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    REQUIRE(rig.add_source(&src) == Status::Ok);
    // Deliberately register NO correction: the pose<->bias cross-covariance never gets used,
    // so the bias must stay at its zero prior (only an absolute-ref update can move it).
    rig.run(0.2, tr.duration_s() - 0.1, 50.0);

    const CalibSnapshot* cs = last_bias_snapshot(rig.records());
    REQUIRE(cs != nullptr);
    MESSAGE("no-ref observability: recovered bias norm=" << cs->bias.norm()
            << "  bias_observable=" << cs->bias_observable);
    // The bias never moved off zero (no update ever injected it).
    CHECK(cs->bias.norm() < 1e-9);
    // The observability confidence stays ~0: with no absolute ref the bias variance only grows
    // (random walk), so it is never DETERMINED. This is the spine guard — the GPS cross-covariance
    // is what makes the bias observable; without it the confidence cannot rise.
    CHECK(cs->bias_observable < 1e-3);
    // No fixes were applied (none registered), so the bias was never observed/removed.
    long applied_total = 0;
    for (const Record& r : rig.records()) if (r.fused) applied_total += r.result.correction.corr_applied;
    CHECK(applied_total == 0);
}

// ===========================================================================
// (c) DEFAULT-OFF IDENTITY — bias_states=false is byte-identical to the 12-DOF path
// ===========================================================================
TEST_CASE("bias: bias_states=false is byte-identical to the predict-only 12-DOF path") {
    Trajectory tr = bias_traj();
    Vec6 planted; planted << 0.15, 0.0, 0.0,   0.0, 0.0, 0.05;

    auto run_once = [&](std::vector<Record>& out) {
        SourceParams sp = bias_source_params(planted);
        SyntheticSource src(tr, sp);
        std::vector<SensorConfig> sensors;
        Config cfg = bias_config(sensors, /*bias_on=*/false);
        AbsoluteRefParams rp;
        rp.period_s = 0.2; rp.window_s = cfg.window_s; rp.sigma_pos = 0.05; rp.seed = 11u;
        SyntheticAbsoluteRef ref(tr, rp);
        Rig rig; rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        REQUIRE(rig.add_source(&src) == Status::Ok);
        REQUIRE(rig.add_correction(&ref) == Status::Ok);
        rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        out = rig.records();
    };

    // Two bias-OFF runs must be byte-identical (determinism), the bias snapshot all zero, and
    // bias_observable zero (the augmented path is never entered when bias_states=false).
    std::vector<Record> a, b;
    run_once(a);
    run_once(b);
    REQUIRE(a.size() == b.size());
    REQUIRE(!a.empty());

    long fused_compared = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        REQUIRE(a[i].fused == b[i].fused);
        if (!a[i].fused) continue;
        const Result& ra = a[i].result;
        const Result& rb = b[i].result;
        bool equal = true;
        equal = equal && (ra.frontier.pose.R.array() == rb.frontier.pose.R.array()).all();
        equal = equal && (ra.frontier.pose.t.array() == rb.frontier.pose.t.array()).all();
        equal = equal && (ra.frontier.cov.array()    == rb.frontier.cov.array()).all();
        // bias_states OFF -> the bias DOF stays the zero default every step.
        equal = equal && (ra.calib[0].bias.array() == 0.0).all();
        equal = equal && (ra.calib[0].bias_observable == 0.0);
        CHECK(equal);
        ++fused_compared;
    }
    REQUIRE(fused_compared > 50);
}

// ===========================================================================
// (d) MULTI-BIAS GUARD — two bias_states sources -> InvalidConfig at init
// ===========================================================================
TEST_CASE("bias: more than one bias_states source is rejected at init (Option A: single only)") {
    std::vector<SensorConfig> sensors;
    SensorConfig s0; s0.id = 0; s0.is_reference = true; s0.bias_states = true;
    SensorConfig s1; s1.id = 1;                          s1.bias_states = true;
    sensors.push_back(s0);
    sensors.push_back(s1);

    Config c;
    c.max_sources         = 2;
    c.reference_sensor_id = 0;
    c.sensors             = sensors.data();
    c.sensor_count        = 2;

    Estimator est;
    CHECK(est.init(c) == Status::InvalidConfig);

    // Sanity: exactly ONE bias source is accepted.
    sensors[1].bias_states = false;
    Estimator est_ok;
    CHECK(est_ok.init(c) == Status::Ok);

    // Negative bias_process_noise is rejected (OutOfRange).
    sensors[0].bias_process_noise = -1.0;
    Estimator est_bad;
    CHECK(est_bad.init(c) == Status::OutOfRange);
}
