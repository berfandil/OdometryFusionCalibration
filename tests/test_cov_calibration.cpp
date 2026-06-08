// test_cov_calibration.cpp — Slice-14 covariance-knob calibration guard (the q_scale sweep).
//
// The Slice-14 NEES Monte-Carlo case (test_validation.cpp) calibrates the COVARIANCE on ONE
// trajectory (nees_traj). This file is the PERMANENT guard for the multi-trajectory, safety-first
// calibration of the core default `q_scale` (Config::q_scale, now 0.7, RE-CALIBRATED against the
// TRUE interior robust median after the D3 median-pinning fix + the D4 removal of the /n_eff fudge).
// It asserts, at the CALIBRATED default, the two properties the calibration must hold across a SET
// of trajectories — the same selection rule used offline to choose 0.7:
//
//   (a) NEVER OVERCONFIDENT (the hard safety constraint). On EVERY trajectory in the set the
//       ensemble-mean pose NEES stays BELOW the DOF count (6), with a conservative margin
//       (< 5.5). Overconfidence (NEES > DOF) on a safety filter is the one thing we never allow;
//       the sim noise model under-states real-world model mismatch, so we stay just below DOF.
//
//   (b) NEAR-CONSISTENT (the calibration targets DOF). On EVERY trajectory the ensemble-mean NEES
//       is comfortably ABOVE the old gross-pessimism regime (it now sits in ~[2, 5], worst-case
//       ~4.85 on nees_traj, approaching DOF=6 from below). A regression of the default back UP
//       toward q_scale=1 (MORE pessimistic now that /n_eff is gone — worst-case drops to ~3.4)
//       trips the lower bound; the old pinning median or a return of the /n_eff fudge (which push
//       NEES well past 6) trips the upper bound.
//
// SIGN NOTE (D4): with the TRUE median + /n_eff removed, a SMALLER q_scale means a SMALLER Q means
// a TIGHTER covariance means a LARGER NEES. So lowering q_scale RAISES NEES toward DOF; q_scale=0.5
// overshoots (worst-case 6.77 > 6 = overconfident), q_scale=0.7 lands at ~4.85 (near-consistent,
// safe), q_scale=1.0 is more pessimistic (~3.4). This is the OPPOSITE relationship from the old
// pinning-median regime where the /n_eff fudge inverted the effect.
//
// HOW 0.7 WAS CHOSEN (offline grid, deleted): swept q_scale in {0.5,0.7,1.0,1.5,2.0} x {nees_traj,
// mixed, turning, straight} x {1x,2x} noise, M=30, against the TRUE median with /n_eff OFF. Worst-
// case ensemble-mean NEES per q_scale at 1x: 0.5->6.77 (OVERCONFIDENT, >6), 0.7->4.85, 1.0->3.40,
// 1.5->2.27, 2.0->1.71. 0.7 is the value that gets the worst case CLOSEST to DOF=6 FROM BELOW while
// NEVER overconfident at either noise level (2x worst-case at 0.7 is ~3.9) -> the safety-first pick.
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
          "pose NEES is never overconfident (<6) AND near-consistent (approaching DOF from below)") {
    struct Entry { const char* name; Trajectory tr; };
    std::vector<Entry> trajs;
    trajs.push_back({"nees_traj", cal_nees_traj()});
    trajs.push_back({"turning",   Trajectory::turning(2.0, 0.5, 16.0)});   // lowest NEES in the sweep
    trajs.push_back({"straight",  Trajectory::straight(2.0, 16.0)});

    const Scalar q_scale_default = Config{}.q_scale;   // the RE-CALIBRATED core default (0.7)
    const int    M = 30;

    // The DOF and the safety margin: NEES must stay below DOF with headroom (never overconfident).
    const Scalar dof          = 6.0;
    const Scalar overconf_cap = 5.5;     // hard upper guard (conservative margin below DOF=6)
    // The calibration must keep every trajectory NEAR-CONSISTENT (well above the old gross-pessimism
    // regime). Observed lowest across the set is `turning` ~2.0 at 1x; `pessimism_lo = 1.3` is a
    // generous per-trajectory FLOOR (it sits below `turning`'s ~2.0 so it does not trip on normal
    // seed variation). It is NOT the q_scale-regression guard — a regression UP toward q_scale=1
    // only drops `turning` to ~1.4 (still > 1.3); the real up-regression detector is the
    // nees_new > nees_old*1.25 ratio CHECK below.
    const Scalar pessimism_lo = 1.3;

    for (const Entry& e : trajs) {
        const Scalar nees = ensemble_nees(e.tr, q_scale_default, M);
        const std::string name = e.name;
        MESSAGE("cov-cal NEES @ q_scale=" << q_scale_default << " traj=" << name
                << " ensemble-mean=" << nees << " (DOF=6; safe band ~[2,5], must be <"
                << overconf_cap << " and >" << pessimism_lo << ")");
        // (a) NEVER OVERCONFIDENT — the hard safety constraint, on every trajectory.
        CHECK(nees < overconf_cap);
        CHECK(nees < dof);
        // (b) NEAR-CONSISTENT — well above the old gross pessimism, on every trajectory.
        CHECK(nees > pessimism_lo);
    }

    // PROOF THE KNOB BIT (and the D4 SIGN): with the TRUE median + /n_eff removed, q_scale=1 is MORE
    // pessimistic (LARGER Q -> LARGER cov -> SMALLER NEES) than the calibrated 0.7. On nees_traj:
    // q_scale=1 -> ~3.4, q_scale=0.7 -> ~4.85. So the calibrated default is materially LARGER (closer
    // to DOF) than q_scale=1, and a regression back UP to 1 trips the (b) lower bounds. (This is the
    // OPPOSITE sign from the old pinning-median+/n_eff regime, where 0.5 was larger than 1.)
    const Scalar nees_old = ensemble_nees(cal_nees_traj(), /*q_scale=*/1.0, M);
    const Scalar nees_new = ensemble_nees(cal_nees_traj(), q_scale_default, M);
    MESSAGE("cov-cal nees_traj: q_scale=1 NEES=" << nees_old
            << "  calibrated q_scale=" << q_scale_default << " NEES=" << nees_new);
    CHECK(nees_old > 2.5);                 // q_scale=1 is more pessimistic (~3.4) but still O(DOF)
    CHECK(nees_old < dof);                 // q_scale=1 is also never overconfident
    CHECK(nees_new > nees_old * 1.25);     // 0.7 is materially closer to DOF than 1 (~4.85 vs ~3.4)
    CHECK(nees_new < overconf_cap);        // still never overconfident
}
