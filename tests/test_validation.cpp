// test_validation.cpp — Slice 14 validation harness: NEES consistency (Monte-Carlo) +
// deterministic golden regression (DESIGN §10, DECISIONS D24).
//
// This slice adds the REMAINING pieces of the trust apparatus on top of the sim rig and
// the per-slice observability self-tests that already exist (test_sim.cpp +
// test_calib_*.cpp): a Monte-Carlo NEES consistency check and a golden-regression drift
// detector. The sim rig (sim/, relaxed edge) is the only place ground truth is known, so
// it is the only place filter-covariance calibration is checkable (D24).
//
// =============================================================================
// NEES (Normalized Estimation Error Squared) — covariance consistency
// =============================================================================
// A filter that PUBLISHES a covariance Σ must prove that Σ is neither over- nor under-
// confident. The 6-DOF pose NEES at a fused record is
//       nees = e^T (P_pp)^{-1} e
// where P_pp is the 6x6 POSE block (rows/cols 0..5) of result.frontier.cov and e is the
// 6-vector pose error between the fused frontier pose and the GT frontier pose, in the
// SAME tangent convention Σ is propagated in.
//
// ERROR CONVENTION (load-bearing — must match the ESKF's covariance, D21 + eskf.cpp):
//   The ESKF propagates a RIGHT-invariant error eta defined by
//         T_true = T_est o exp(eta)          (eskf.cpp:7-13)
//   i.e. exp(eta) = T_est^{-1} o T_true, with eta in [trans; rot] order (pose 0..5). The
//   propagation Jacobian is the SE(3) adjoint Ad(delta^-1), so the covariance lives in the
//   full SE(3) tangent, NOT a decoupled SO(3)xR^3 split. We therefore form the error with
//   the SAME full SE(3) log used for propagation:
//         e = se3::log( T_est^{-1} o T_gt )      (Vec6 = [v; omega], trans-first)
//   se3::log returns [trans; rot] (lie.hpp), so e lines up index-for-index with P_pp.
//   Rotation rows 3..5 are so3::log(R_est^T R_gt) embedded in the SE(3) log; translation
//   rows 0..2 are the V^{-1}-mapped residual consistent with that same right tangent (NOT
//   a raw t_gt - t_est, which would mix conventions). The convention is validated
//   ALGEBRAICALLY (the reconstruction identity T_est o exp(e) == T_gt — see the first NEES
//   case) rather than by an "ensemble-mean NEES ~ DOF" band: because this filter's
//   covariance is grossly pessimistic (init P = I_12, ~100x the actual error), EVERY
//   sub-block's NEES is ~0.01..0.1 regardless of convention, so a NEES magnitude cannot
//   distinguish right from wrong here. For a CORRECTLY-CALIBRATED linear-Gaussian filter the
//   ensemble-mean NEES would sit near the DOF count (6); this one sits at ~0.13 (pessimistic).
//
// MONTE-CARLO: the SAME scenario is run across an ensemble of M independent seeds (each
// SyntheticSource's `seed` varies per run) with genuine injected noise so the covariance
// is actually exercised. We average the per-record pose NEES over steady-state ticks (drop
// warmup) and over the ensemble, then test the ensemble-mean NEES against a chi-square
// consistency interval for 6 DOF.
//
// HONESTY CLAUSE (DESIGN §10 / the brief): the integrator is PREDICT-ONLY with a HEURISTIC
// adaptive Q (q = q_scale*spread^2 + q_floor), NOT a calibrated noise model. So the
// published Σ can be genuinely mis-calibrated and the strict chi-square test may fail.
// That is a VALID scientific outcome. We first try to achieve consistency by choosing the
// TEST's sim Config Q knobs (q_scale, q_floor) and noise levels (never the core defaults);
// if the filter is still inconsistent we pin a DOCUMENTED ACHIEVABLE bound on the value the
// filter actually produces (a non-vacuous regression guard) and the report surfaces the
// mis-calibration as an open question. See the case body for the empirical outcome.
//
// =============================================================================
// NIS deferral
// =============================================================================
// NIS (Normalized INNOVATION Squared) is defined only at a MEASUREMENT UPDATE (innovation
// = measurement - predicted measurement, gated by its innovation covariance). The core is
// PREDICT-ONLY: there are NO absolute-reference corrections yet (those land in Slice 11,
// ICorrection / Mahalanobis-gated updates, D22). There is therefore no innovation to
// normalize, so NIS is NOT computable now and is DEFERRED to Slice 11. We deliberately do
// NOT fake it. (Reported to the orchestrator as a documented deferral.)
//
// =============================================================================
// Golden regression — deterministic drift detector
// =============================================================================
// Style chosen: (A) SELF-CONTAINED determinism + committed numeric values. A fully-seeded,
// fixed scenario is run through the Rig twice and asserted BYTE-IDENTICAL across the two
// in-test runs (fused pose + covariance + the per-source CALIB snapshot — extending the
// existing test_sim determinism case, which pins only the fused pose/cov, to the calib
// snapshot). Then a small set of committed expected NUMERIC values is compared within a
// TIGHT tolerance. The committed-value scenario is NOISE-FREE so the mt19937 /
// normal_distribution cross-stdlib non-determinism (HANDOFF §6 caveat: std::normal_-
// distribution's transform of the engine stream is unspecified across libstdc++/libc++/
// MSVC) cannot make the committed values flaky — the noise-free path is pure SE(3)
// arithmetic and portable. Style (B) (a committed tests/golden/ fixture file) was NOT
// chosen: a fixture buys nothing over in-source committed constants for a snapshot this
// small, and would add a file to keep in sync. The golden pins the fused frontier pose
// (t + R) at sampled ticks + the final-tick calib snapshot (extrinsic/scale/offset +
// confidences) so a deliberate change to fusion OR calibration output breaks it.
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
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

// A mixed straight + multi-axis turn trajectory: straight stretches (so the integrator
// sees translation), turns about y and z (so rotation covariance is exercised on all
// axes), repeated to give a long steady-state run after warmup. Mirrors the shape used by
// the feedback/weights tests so the regimes are representative.
Trajectory nees_traj() {
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

// The 6-DOF pose NEES at a fused record (see the file header for the convention):
//   e = se3::log( T_est^{-1} o T_gt )   (full SE(3) right tangent, [trans; rot])
//   nees = e^T (P_pp)^{-1} e            P_pp = pose block (0..5) of frontier.cov
Scalar pose_nees(const Record& r) {
    const SE3& est = r.result.frontier.pose;
    const SE3& gt  = r.gt_frontier;
    const SE3 err_T = se3::compose(se3::inverse(est), gt);   // T_est^{-1} o T_gt
    const Vec6 e    = se3::log(err_T);                        // [trans; rot]
    const Mat6 Ppp  = r.result.frontier.cov.block<6, 6>(0, 0);
    // P_pp is small (12x12 dense, symmetric PSD) — a direct ldlt solve is fine and stable.
    const Vec6 Pinv_e = Ppp.ldlt().solve(e);
    return e.dot(Pinv_e);
}

// Build a 3-source rig config whose priors EQUAL the planted params (no calibration needed
// — this isolates the COVARIANCE-vs-error question from calibration transients). q_scale /
// q_floor are TEST knobs we are allowed to tune to chase consistency (never core defaults).
Config nees_config(const std::vector<SourceParams>& planted,
                   std::vector<SensorConfig>& sensors_out,
                   Scalar q_scale, const Scalar (&q_floor)[6]) {
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
    c.cold_start     = ColdStart::MedianFromStart;   // priors == planted -> no calib needed
    c.adaptive_q     = true;
    c.q_scale        = q_scale;
    for (int i = 0; i < 6; ++i) c.q_floor[i] = q_floor[i];
    c.sensors        = sensors_out.data();
    c.sensor_count   = static_cast<int>(sensors_out.size());
    return c;
}
} // namespace

// ===========================================================================
// DELIVERABLE 1 — NEES Monte-Carlo consistency
// ===========================================================================

TEST_CASE("validation NEES: error convention is the ESKF right-error (reconstruction "
          "identity + it is the minimal tangent)") {
    // The error convention is validated ALGEBRAICALLY, not by an O(DOF) NEES band. We
    // originally tried "rotation sub-NEES ~ DOF=3" as the convention probe, but the filter's
    // covariance is so PESSIMISTIC (init P = I_12, ~100x the actual error — see the next
    // case) that EVERY sub-block's NEES is ~0.01..0.1 regardless of convention. So a NEES
    // magnitude cannot distinguish a right convention from a wrong one here. Instead we pin
    // the convention by the two properties that actually define it:
    //
    //  (1) RECONSTRUCTION IDENTITY. The ESKF defines its error by T_true = T_est o exp(eta)
    //      (eskf.cpp:7-13, right-invariant, eta in [trans;rot]). Our e = se3::log(T_est^{-1}
    //      o T_gt) must therefore satisfy  T_est o exp(e) == T_gt  EXACTLY. That is the
    //      definition of "e is the right-error tangent the covariance is expressed in". A
    //      left-error (e' = log(T_gt o T_est^{-1})) does NOT reconstruct via T_est o exp(e').
    //
    //  (2) MINIMAL-TANGENT / ORDER. se3::log returns [trans; rot], matching the cov's pose
    //      block order (lie.hpp + eskf.cpp). We confirm the rotation rows are so3::log(R_est^T
    //      R_gt) (NOT the transpose R_gt^T R_est, which negates the tangent) and the trans
    //      rows are the SE(3)-log (V^{-1}-mapped) residual, not a raw t_gt - t_est.
    //
    // These are deterministic algebraic checks on a representative fused record — no noise,
    // no Monte-Carlo needed (the convention is a property of the math, not the statistics).
    Trajectory tr = nees_traj();

    // A lightly-noised, seeded rig run; priors == planted. Noise-FREE here would track GT to
    // machine precision (the const-twist GT is exactly what the estimator integrates), giving
    // a zero error that proves nothing — so we inject a small seeded noise to create a genuine
    // residual to round-trip. The convention checks below are algebraic and hold for ANY e.
    std::vector<SourceParams> planted(3);
    planted[0].id = 0;
    planted[1].id = 1;
    planted[1].X.R = so3::exp(Vec3(0, 0, kPi / 6));
    planted[1].X.t = Vec3(0.3, -0.2, 0.1);
    planted[2].id = 2;
    planted[2].X.R = so3::exp(Vec3(0.05, 0.1, -kPi / 7));
    planted[2].X.t = Vec3(-0.25, 0.15, 0.05);
    for (auto& sp : planted) {
        sp.noise_trans_per_m = 0.02;
        sp.noise_rot_per_rad = 0.02;
        sp.noise_trans_floor = 0.005;
        sp.noise_rot_floor   = 0.005;
        sp.seed = 314u + sp.id;
    }

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
    std::vector<SensorConfig> sensors;
    const Scalar qf[6] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
    Config cfg = nees_config(planted, sensors, 1.0, qf);

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
    rig.run(0.2, tr.duration_s() - 0.1, 50.0);

    // Pick the fused record with the LARGEST pose error, so the round-trip and the
    // left-vs-right / raw-vs-SE(3)log distinctions are numerically meaningful (a near-zero-
    // error record would make every convention reconstruct to within round-off and prove
    // nothing). A genuine non-zero error exists because the windowed composition does not
    // exactly equal the const-twist GT, especially through the turns.
    const Record* probe = nullptr;
    Scalar best_err = 0.0;
    for (const Record& r : rig.records()) {
        if (!r.fused) continue;
        const Vec6 ee = se3::log(se3::compose(se3::inverse(r.result.frontier.pose),
                                              r.gt_frontier));
        if (ee.norm() > best_err) { best_err = ee.norm(); probe = &r; }
    }
    REQUIRE(probe != nullptr);

    const SE3& est = probe->result.frontier.pose;
    const SE3& gt  = probe->gt_frontier;

    // (1) RECONSTRUCTION IDENTITY: e = log(T_est^{-1} T_gt) must reconstruct T_gt via the
    //     RIGHT composition T_est o exp(e) (the ESKF's error definition).
    const Vec6 e = se3::log(se3::compose(se3::inverse(est), gt));
    REQUIRE(e.norm() > 1e-6);                          // a genuine non-zero error to test
    const SE3 reconstructed = se3::compose(est, se3::exp(e));
    CHECK((reconstructed.R - gt.R).cwiseAbs().maxCoeff() < 1e-12);
    CHECK((reconstructed.t - gt.t).cwiseAbs().maxCoeff() < 1e-12);

    // The LEFT-error tangent does NOT reconstruct via the right composition (proves the
    // convention is genuinely right-, not left-, error — they differ on a non-commuting SE3).
    const Vec6 e_left = se3::log(se3::compose(gt, se3::inverse(est)));
    const SE3 wrong   = se3::compose(est, se3::exp(e_left));
    CHECK((wrong.t - gt.t).cwiseAbs().maxCoeff() > 1e-6);

    // (2a) ORDER: the rotation rows (3..5) equal so3::log(R_est^T R_gt) — NOT the transpose.
    const Vec3 rot_err = so3::log(est.R.transpose() * gt.R);
    CHECK((e.tail<3>() - rot_err).cwiseAbs().maxCoeff() < 1e-12);
    CHECK((e.tail<3>() + so3::log(gt.R.transpose() * est.R)).cwiseAbs().maxCoeff() < 1e-12);

    // (2b) The trans rows (0..2) are the SE(3)-log residual, which for a non-zero rotation
    //      error differs from the naive (t_gt - t_est) expressed in the est frame — confirming
    //      we use the full SE(3) tangent the covariance is propagated in, not a raw residual.
    const Vec3 raw_trans = est.R.transpose() * (gt.t - est.t);
    CHECK((e.head<3>() - raw_trans).cwiseAbs().maxCoeff() > 1e-9);

    MESSAGE("convention OK: T_est o exp(log(T_est^-1 T_gt)) == T_gt; |e|=" << e.norm());
}

TEST_CASE("validation NEES: Monte-Carlo covariance consistency (chi-square interval, 6 "
          "DOF) — documents a PESSIMISTIC covariance") {
    // The non-vacuous consistency test. Ensemble-mean pose NEES over N steady-state samples
    // is compared to a chi-square consistency interval for 6 DOF. For N averaged samples the
    // sum N*nees ~ chi2(N*6) under the consistency hypothesis, so the bound on the MEAN is
    //     [ chi2_inv(alpha/2, N*6)/N , chi2_inv(1-alpha/2, N*6)/N ].
    // With N in the thousands the interval is TIGHT (a Wilson-Hilferty normal approximation
    // to chi2 gives, for large k=N*6, mean ~ 6*(1 +/- z*sqrt(2/k))) — so a consistent filter
    // would be pinned hard, and an inconsistent one is flagged clearly.
    //
    // EMPIRICAL OUTCOME (the honesty clause, DESIGN §10): this filter is NOT consistent. The
    // ensemble-mean pose NEES is ~0.13, far BELOW the ~6 the chi-square interval demands —
    // the published covariance is GROSSLY PESSIMISTIC (too large) relative to the actual
    // fused-vs-GT error. Root cause (diagnosed, see report):
    //   (1) The ESKF initializes P = I_12 (eskf.cpp / estimator.cpp:910), so the pose block
    //       starts at 1.0 rad^2 / 1.0 m^2 — already ~100x the steady-state squared error
    //       (|e| ~ 0.09, e^T e ~ 0.009). A predict-only filter has NO correction to shrink
    //       P, so it can only grow from there.
    //   (2) The right-error Ad(delta^-1) propagation INFLATES the translation block with
    //       distance travelled (P_tt grows to O(3..20) m^2 over the run), compounding (1).
    // This CANNOT be fixed by tuning the test's Q knobs (q_scale/q_floor only ADD to P,
    // making it MORE pessimistic) — it needs a SMALLER initial P and/or a correction step,
    // both of which are STRICT-CORE changes out of scope for this validation slice. So per
    // the brief's honesty clause we: (a) compute + report the actual NEES and the chi-square
    // bound it violates, and (b) assert a DOCUMENTED ACHIEVABLE band on the value the filter
    // actually produces, making this a non-vacuous regression guard (a 3x change in either P
    // or the tracking error breaks it) while explicitly recording the mis-calibration.
    // Surfaced to the orchestrator as an open question (covariance calibration → likely a
    // core init-P fix and feeds the CONFIG Q sweep).
    Trajectory tr = nees_traj();
    const int M = 30;

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

    // TEST-TUNED Q (allowed — never core defaults). Chosen to best match the published Σ to
    // the injected noise; see the report for the chase.
    const Scalar q_scale = 1.0;
    const Scalar qf[6]   = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};

    Scalar sum_nees = 0.0;
    long   N        = 0;
    for (int run = 0; run < M; ++run) {
        std::vector<SourceParams> planted = base;
        for (auto& sp : planted) sp.seed = 5000u + run * 11u + sp.id;

        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
        std::vector<SensorConfig> sensors;
        Config cfg = nees_config(planted, sensors, q_scale, qf);

        Rig rig;
        rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
        rig.run(0.2, tr.duration_s() - 0.1, 50.0);

        const std::vector<Record>& recs = rig.records();
        int fused_seen = 0;
        const int warmup = 20;
        for (const Record& r : recs) {
            if (!r.fused) continue;
            ++fused_seen;
            if (fused_seen <= warmup) continue;
            sum_nees += pose_nees(r);
            ++N;
        }
    }
    REQUIRE(N > 1000);
    const Scalar mean_nees = sum_nees / static_cast<Scalar>(N);

    // Two-sided 99% chi-square interval on the mean for k = N*6 DOF, via Wilson-Hilferty:
    //   chi2_inv(p, k) ~ k * (1 - 2/(9k) + z_p*sqrt(2/(9k)))^3 ,  z = +/-2.576 for 99%.
    const Scalar dof = 6.0;
    const Scalar k   = static_cast<Scalar>(N) * dof;
    const Scalar z   = 2.5758;     // 99% two-sided
    auto wh = [&](Scalar zp) {
        const Scalar a = Scalar(2) / (Scalar(9) * k);
        const Scalar f = Scalar(1) - a + zp * std::sqrt(a);
        return k * f * f * f;
    };
    const Scalar lo_mean = wh(-z) / static_cast<Scalar>(N);
    const Scalar hi_mean = wh(+z) / static_cast<Scalar>(N);

    const bool truly_consistent = (mean_nees >= lo_mean && mean_nees <= hi_mean);
    MESSAGE("NEES consistency: ensemble-mean=" << mean_nees << " DOF=6 N=" << N
            << " chi2 99% interval on mean=[" << lo_mean << ", " << hi_mean << "]"
            << "  truly_consistent=" << (truly_consistent ? "YES" : "NO")
            << "  (mean << 6 => PESSIMISTIC / over-conservative covariance)");

    // TRUE-CONSISTENCY verdict, recorded as a (currently-failing-by-design) expectation so
    // the day the core gains a calibrated init-P / correction step this flips to PASS and the
    // ACHIEVABLE band below should be tightened toward the chi-square interval. We assert the
    // verdict is FALSE here to PIN the documented inconsistency (a silent flip to consistent —
    // e.g. someone fixes init-P — should prompt revisiting this test's bounds).
    CHECK_FALSE(truly_consistent);     // documents: filter is NOT chi-square-consistent today

    // DOCUMENTED ACHIEVABLE BOUND (non-vacuous regression guard on the value the filter
    // ACTUALLY produces). Observed ensemble-mean ~0.13..0.16 across seed bases; the band
    // [0.04, 0.40] absorbs Monte-Carlo seed variation yet a ~3x drift in EITHER the published
    // covariance OR the fused-vs-GT tracking error breaks it. This is NOT the consistency
    // interval — it pins the (pessimistic) covariance the filter emits today.
    CHECK(mean_nees > 0.04);
    CHECK(mean_nees < 0.40);
}

// ===========================================================================
// DELIVERABLE 2 — Golden regression (deterministic drift detector)
// ===========================================================================
// Style (A): self-contained determinism + committed numeric values (see file header for the
// rationale vs a committed-fixture file). Two parts:
//   PART 1 — byte-identical replay of the fused pose + covariance AND the per-source calib
//            snapshot across two in-test runs (extends test_sim's determinism case, which
//            pins only pose/cov, to the calibration output).
//   PART 2 — a NOISE-FREE scenario's committed numeric golden: the fused frontier pose at
//            sampled ticks + the final-tick calib snapshot, compared within a TIGHT tolerance.
//            Noise-free => no std::normal_distribution in the path => portable across stdlib
//            (HANDOFF §6 caveat avoided); the values are pure SE(3) arithmetic.
//
// WHAT THE GOLDEN PINS (its drift-detection contract): the exact fused-frontier trajectory
// the integrator emits for this fixed (trajectory + planted sources + config + seeds)
// scenario, plus the calibrator's committed extrinsic/scale/offset estimates and per-DOF
// confidences. A deliberate change to fusion (median, frame-align, ESKF predict, Q) OR to
// calibration (histogram votes, commit gates, feedback) shifts these numbers and fails the
// committed-value checks. REGENERATION: if an intended algorithm change moves the golden,
// re-run this case once with the MESSAGE() lines below uncommented (or read the failing
// expected-vs-actual) and paste the new constants into kGolden* / the snapshot expectations.

namespace {
// The fixed golden scenario: a deterministic trajectory + 3 planted sources with real mounts
// + a fully-specified config. Shared by the determinism and committed-value parts.
Trajectory golden_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0,    0,   0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0,  0.3,  0.5;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.3,  0.5;
    tr.add_segment(straight, 2.0);
    tr.add_segment(turnA,    1.5);
    tr.add_segment(straight, 1.5);
    tr.add_segment(turnB,    1.5);
    return tr;
}

void golden_planted(std::vector<SourceParams>& planted, bool with_noise) {
    planted.assign(3, SourceParams{});
    planted[0].id = 0;                                   // reference, identity mount
    planted[1].id = 1;
    planted[1].X.R = so3::exp(Vec3(0, 0, kPi / 6));
    planted[1].X.t = Vec3(0.3, -0.2, 0.1);
    planted[1].scale = 1.1;
    planted[2].id = 2;
    planted[2].X.R = so3::exp(Vec3(0.05, 0.1, -kPi / 7));
    planted[2].X.t = Vec3(-0.25, 0.15, 0.05);
    planted[2].scale = 0.95;
    if (with_noise) {
        for (auto& sp : planted) {
            sp.noise_trans_per_m = 0.01;
            sp.noise_rot_per_rad = 0.005;
            sp.noise_trans_floor = 0.002;
            sp.seed = 700u + sp.id;
        }
    }
}

// Config with priors == planted (so fusion tracks GT; the calibrator still OBSERVES and
// publishes a snapshot). Calibration histograms configured so the snapshot is populated.
Config golden_config(const std::vector<SourceParams>& planted,
                     std::vector<SensorConfig>& sensors_out) {
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
    c.max_sources    = 3;
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.cold_start     = ColdStart::MedianFromStart;
    c.timesync_enabled = false;        // deterministic, offset-free golden
    c.vote_weight    = VoteWeight::One;
    c.sensors        = sensors_out.data();
    c.sensor_count   = 3;
    return c;
}

const CalibSnapshot* find_calib(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.calib[i].id == id) return &r.calib[i];
    return nullptr;
}
} // namespace

TEST_CASE("validation golden: byte-identical replay of fused pose, covariance, AND calib "
          "snapshot") {
    // PART 1 — determinism over the calibration output too (test_sim pins pose+cov only).
    Trajectory tr = golden_traj();
    auto run_once = [&](std::vector<Record>& out) {
        std::vector<SourceParams> planted;
        golden_planted(planted, /*with_noise=*/true);     // determinism must hold WITH noise
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
        std::vector<SensorConfig> sensors;
        Config cfg = golden_config(planted, sensors);
        Rig rig;
        rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
        rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        out = rig.records();
    };
    std::vector<Record> a, b;
    run_once(a);
    run_once(b);
    REQUIRE(a.size() == b.size());
    REQUIRE(!a.empty());
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].fused == b[i].fused);
        if (!a[i].fused) continue;
        const Result& ra = a[i].result;
        const Result& rb = b[i].result;
        // Fused pose + covariance (as test_sim does)...
        CHECK((ra.frontier.pose.R.array() == rb.frontier.pose.R.array()).all());
        CHECK((ra.frontier.pose.t.array() == rb.frontier.pose.t.array()).all());
        CHECK((ra.frontier.cov.array()    == rb.frontier.cov.array()).all());
        // ...AND the per-source calib snapshot (the extension this golden adds).
        REQUIRE(ra.source_count == rb.source_count);
        for (int s = 0; s < ra.source_count; ++s) {
            const CalibSnapshot& ca = ra.calib[s];
            const CalibSnapshot& cb = rb.calib[s];
            CHECK(ca.id == cb.id);
            CHECK((ca.extrinsic.R.array() == cb.extrinsic.R.array()).all());
            CHECK((ca.extrinsic.t.array() == cb.extrinsic.t.array()).all());
            CHECK(ca.scale         == cb.scale);
            CHECK(ca.time_offset_s == cb.time_offset_s);
            CHECK(ca.extrinsic_confidence   == cb.extrinsic_confidence);
            CHECK(ca.scale_confidence       == cb.scale_confidence);
            CHECK(ca.translation_confidence == cb.translation_confidence);
            CHECK(ca.extrinsic_committed   == cb.extrinsic_committed);
            CHECK(ca.scale_committed       == cb.scale_committed);
            CHECK(ca.translation_committed == cb.translation_committed);
        }
    }
}

TEST_CASE("validation golden: committed numeric values (noise-free, portable) — fused pose "
          "+ final calib snapshot") {
    // PART 2 — the committed-value drift detector. Noise-free, so the values are exact SE(3)
    // arithmetic with no std::normal_distribution (portable across stdlib; HANDOFF §6).
    Trajectory tr = golden_traj();
    std::vector<SourceParams> planted;
    golden_planted(planted, /*with_noise=*/false);
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
    std::vector<SensorConfig> sensors;
    Config cfg = golden_config(planted, sensors);

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 50);

    const std::vector<Record>& recs = rig.records();

    // --- Sampled fused-frontier poses at fixed fused-record indices (golden) -------------
    // Collect the fused records in order; sample a few indices spanning the run.
    std::vector<const Record*> fused;
    for (const Record& r : recs) if (r.fused) fused.push_back(&r);
    REQUIRE(fused.size() > 60);

    // GOLDEN CONSTANTS (captured from a clean run; regenerate per the header if an intended
    // algorithm change moves them). Pinned at FIXED fused-record indices (not size-relative,
    // so the golden is stable). idx 10/40 are in the opening straight (pose is pure +x, all
    // else ~0); idx 300 is deep in the run after both turns (a full 6-DOF pose).
    struct PoseSample { int idx; Scalar tx, ty, tz; Scalar rx, ry, rz; };
    const PoseSample kGoldenPose[] = {
        { 10, 0.6, 0.0,     0.0,     0.0,      0.0,       0.0      },
        { 40, 1.8, 0.0,     0.0,     0.0,      0.0,       0.0      },
        {300, 9.554356, 4.960752, -2.212240, 0.264728, 0.112838, 1.288650},
    };
    const Scalar kAbsTol = 1e-5;     // absolute tolerance (machine-precision tracking; 1e-5
                                     // absorbs the ~1e-14 round-off + 6-digit literal rounding)
    for (const PoseSample& g : kGoldenPose) {
        REQUIRE(g.idx < static_cast<int>(fused.size()));
        const State& f = fused[g.idx]->result.frontier;
        const Vec3 lr = so3::log(f.pose.R);
        CHECK(std::abs(f.pose.t.x() - g.tx) < kAbsTol);
        CHECK(std::abs(f.pose.t.y() - g.ty) < kAbsTol);
        CHECK(std::abs(f.pose.t.z() - g.tz) < kAbsTol);
        CHECK(std::abs(lr.x() - g.rx) < kAbsTol);
        CHECK(std::abs(lr.y() - g.ry) < kAbsTol);
        CHECK(std::abs(lr.z() - g.rz) < kAbsTol);
    }

    // Overall tracking (noise-free, priors == planted => GT to machine precision).
    Scalar max_te, max_re;
    rig.max_error(max_te, max_re);
    CHECK(max_te < 1e-9);
    CHECK(max_re < 1e-9);

    // Final-tick calib snapshot for source 1 (pins the calibrator's running output: the
    // scale estimate mid-convergence + the per-DOF confidences, all deterministic here).
    const Result& last = fused.back()->result;
    const CalibSnapshot* c1 = find_calib(last, 1);
    REQUIRE(c1 != nullptr);
    CHECK(std::abs(c1->scale - 1.08281) < 1e-4);
    CHECK(std::abs(c1->extrinsic_confidence   - 1.0) < 1e-4);
    CHECK(std::abs(c1->scale_confidence       - 1.0) < 1e-4);
    CHECK(std::abs(c1->translation_confidence - 1.0) < 1e-4);
}
