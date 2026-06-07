// test_cov_calibration.cpp — Slice-14 covariance-knob calibration guard (the q_scale sweep).
//
// The Slice-14 NEES Monte-Carlo case (test_validation.cpp) calibrates the COVARIANCE on ONE
// trajectory (nees_traj). This file is the PERMANENT guard for the multi-trajectory, safety-first
// calibration of the core default `q_scale` (Config::q_scale, now 0.5, was the un-calibrated 1.0
// placeholder). It asserts, at the CALIBRATED default, the two properties the calibration must hold
// across a SET of trajectories — the same selection rule used offline to choose 0.5:
//
//   (a) NEVER OVERCONFIDENT (the hard safety constraint). On EVERY trajectory in the set the
//       ensemble-mean pose NEES stays BELOW the DOF count (6), with a conservative margin
//       (< 5.5). Overconfidence (NEES > DOF) on a safety filter is the one thing we never allow;
//       the sim noise model under-states real-world model mismatch, so we stay mildly pessimistic.
//
//   (b) PESSIMISM MATERIALLY REDUCED (the calibration actually bit). On EVERY trajectory the
//       ensemble-mean NEES is materially ABOVE the old un-calibrated value (~1.0 at q_scale=1):
//       the calibrated 0.5 lifts every trajectory's NEES into ~[2, 3] (worst-case ~2.9), cutting
//       the residual predict-only pessimism roughly in half. A regression of the default back
//       toward 1.0 (less calibrated) trips the lower bound on every trajectory.
//
// HOW 0.5 WAS CHOSEN (offline grid, deleted): swept q_scale in {1.0,0.5,0.3,0.2,0.15,0.1} x
// {nees_traj, mixed, turning, straight} x {1x,2x} noise, M=30. Worst-case ensemble-mean NEES per
// q_scale at 1x: 1.0->1.48, 0.5->2.93, 0.3->4.81, 0.2->7.09(OVERCONFIDENT). At 2x the same crossing
// is sharper (0.3->5.8). 0.5 is the largest pessimism cut whose worst case stays inside [2,4] AND
// never exceeds 6 at either noise level -> the safety-first pick. Lower values overfit the no-ref
// NEES toward 6 on the optimistic sim and risk overconfidence under real model mismatch.
//
// NOTE: like the validation NEES case, this drives the CORE default q_scale through the sim rig
// (sim/, the only place GT is known, D24). It reuses test_validation's noise level + 3-source rig
// and a SUBSET of trajectories kept cheap for CI (M=30, 3 trajectories).
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
#include <string>
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

// Same mixed straight + multi-axis-turn trajectory the validation NEES case uses.
Trajectory cal_nees_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0,    0,    0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0,  0.35,  0.6;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.35,  0.6;
    for (int rep = 0; rep < 3; ++rep) {
        tr.add_segment(straight, 1.6);
        tr.add_segment(turnA,    1.4);
        tr.add_segment(straight, 1.2);
        tr.add_segment(turnB,    1.4);
    }
    return tr;
}

// e = se3::log(T_est^-1 T_gt) (full SE(3) right tangent, [trans;rot]); nees = e^T P_pp^-1 e.
Scalar pose_nees(const Record& r) {
    const SE3& est = r.result.frontier.pose;
    const SE3& gt  = r.gt_frontier;
    const SE3 err_T = se3::compose(se3::inverse(est), gt);
    const Vec6 e    = se3::log(err_T);
    const Mat6 Ppp  = r.result.frontier.cov.block<6, 6>(0, 0);
    const Vec6 Pinv_e = Ppp.ldlt().solve(e);
    return e.dot(Pinv_e);
}

// 3 mounted sources (== test_validation's NEES rig) with q_scale taken from the passed cfg knob.
Config cal_config(const std::vector<SourceParams>& planted,
                  std::vector<SensorConfig>& sensors_out, Scalar q_scale) {
    sensors_out.clear();
    for (const SourceParams& sp : planted) {
        SensorConfig sc;
        sc.id              = sp.id;
        sc.prior_extrinsic = sp.X;
        sc.prior_scale     = sp.scale;
        sc.weight_prior    = 1.0;
        sensors_out.push_back(sc);
    }
    Config c;
    c.max_sources    = static_cast<int>(planted.size());
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.cold_start     = ColdStart::MedianFromStart;
    c.adaptive_q     = true;
    c.q_scale        = q_scale;
    for (int i = 0; i < 6; ++i) c.q_floor[i] = 1e-6;
    c.sensors        = sensors_out.data();
    c.sensor_count   = static_cast<int>(sensors_out.size());
    return c;
}

std::vector<SourceParams> cal_base() {
    std::vector<SourceParams> base(3);
    base[0].id = 0;
    base[1].id = 1;
    base[1].X.R = so3::exp(Vec3(0, 0, kPi / 6));
    base[1].X.t = Vec3(0.3, -0.2, 0.1);
    base[2].id = 2;
    base[2].X.R = so3::exp(Vec3(0.05, 0.1, -kPi / 7));
    base[2].X.t = Vec3(-0.25, 0.15, 0.05);
    for (auto& sp : base) {
        sp.noise_trans_per_m = 0.02;
        sp.noise_rot_per_rad = 0.02;
        sp.noise_trans_floor = 0.005;
        sp.noise_rot_floor   = 0.005;
    }
    return base;
}

// Ensemble-mean steady-state pose NEES over M seeds on one trajectory at the given q_scale.
Scalar ensemble_nees(const Trajectory& tr, Scalar q_scale, int M) {
    std::vector<SourceParams> base = cal_base();
    Scalar sum = 0.0; long N = 0;
    for (int run = 0; run < M; ++run) {
        std::vector<SourceParams> planted = base;
        for (auto& sp : planted) sp.seed = 5000u + run * 11u + sp.id;
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
        std::vector<SensorConfig> sensors;
        Config cfg = cal_config(planted, sensors, q_scale);
        Rig rig;
        rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
        rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        int fused_seen = 0; const int warmup = 20;
        for (const Record& r : rig.records()) {
            if (!r.fused) continue;
            ++fused_seen;
            if (fused_seen <= warmup) continue;
            sum += pose_nees(r); ++N;
        }
    }
    REQUIRE(N > 1000);
    return sum / static_cast<Scalar>(N);
}
} // namespace

TEST_CASE("cov calibration: at the calibrated default q_scale every trajectory's ensemble-mean "
          "pose NEES is never overconfident (<6) AND materially less pessimistic than q_scale=1") {
    struct Entry { const char* name; Trajectory tr; };
    std::vector<Entry> trajs;
    trajs.push_back({"nees_traj", cal_nees_traj()});
    trajs.push_back({"turning",   Trajectory::turning(2.0, 0.5, 16.0)});   // worst case in the sweep
    trajs.push_back({"straight",  Trajectory::straight(2.0, 16.0)});

    const Scalar q_scale_default = Config{}.q_scale;   // the CALIBRATED core default (0.5)
    const int    M = 30;

    // The DOF and the safety margin: NEES must stay below DOF with headroom (never overconfident).
    const Scalar dof          = 6.0;
    const Scalar overconf_cap = 5.5;     // hard upper guard (conservative margin below DOF=6)
    // The calibration must MATERIALLY beat the old un-calibrated value (~1.0 at q_scale=1): every
    // trajectory's NEES must clear this. Observed worst (nees_traj) ~2.07 -> 1.6 leaves headroom.
    const Scalar pessimism_lo = 1.6;

    for (const Entry& e : trajs) {
        const Scalar nees = ensemble_nees(e.tr, q_scale_default, M);
        const std::string name = e.name;
        MESSAGE("cov-cal NEES @ q_scale=" << q_scale_default << " traj=" << name
                << " ensemble-mean=" << nees << " (DOF=6; safe band ~[2,4], must be <"
                << overconf_cap << " and >" << pessimism_lo << ")");
        // (a) NEVER OVERCONFIDENT — the hard safety constraint, on every trajectory.
        CHECK(nees < overconf_cap);
        CHECK(nees < dof);
        // (b) PESSIMISM REDUCED — the calibration bit (NEES well above the old ~1.0), every traj.
        CHECK(nees > pessimism_lo);
    }

    // PROOF THE KNOB BIT: at the OLD un-calibrated default (q_scale=1) the worst-case trajectory's
    // NEES is ~1.0-1.5 (pessimistic). The calibrated 0.5 is materially larger on the SAME rig, so a
    // regression of the default back toward 1.0 trips the (b) lower bound above. We pin the
    // contrast directly on nees_traj: 1.0 -> ~1.04, 0.5 -> ~2.07 (a ~2x pessimism cut).
    const Scalar nees_old = ensemble_nees(cal_nees_traj(), /*q_scale=*/1.0, M);
    const Scalar nees_new = ensemble_nees(cal_nees_traj(), q_scale_default, M);
    MESSAGE("cov-cal nees_traj: old q_scale=1 NEES=" << nees_old
            << "  calibrated q_scale=" << q_scale_default << " NEES=" << nees_new);
    CHECK(nees_old < 1.5);                 // the old default was pessimistic (~1.0)
    CHECK(nees_new > nees_old * 1.6);      // the calibration is a material, real improvement
    CHECK(nees_new < overconf_cap);        // still never overconfident
}
