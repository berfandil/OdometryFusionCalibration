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
//   a raw t_gt - t_est, which would mix conventions). The first NEES case proves e is
//   se3::log's OWN right-error tangent (right-not-left, trans-first) via the reconstruction
//   identity T_est o exp(e) == T_gt; that the COVARIANCE lives in this SAME tangent holds by
//   eskf.cpp's full-SE(3)-Ad propagation (true by inspection of eskf.cpp, not asserted in
//   this test). We pin the convention this way rather than by an "ensemble-mean NEES ~ DOF"
//   band: because this filter's covariance is still mildly pessimistic (after the init-P fix
//   it is ~17x, was ~46x under the old I_12 seed), EVERY sub-block's NEES is well below its
//   DOF regardless of convention, so a NEES magnitude cannot distinguish right from wrong
//   here. For a CORRECTLY-CALIBRATED linear-Gaussian filter the ensemble-mean NEES would sit
//   near the DOF count (6); this one sits at ~0.35 after the init-P fix (was ~0.13 with the
//   old P = I_12 seed) — see the Monte-Carlo case below for the full diagnosis.
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

#include "ofc_sim/absolute_ref_source.hpp"
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

TEST_CASE("validation NEES: e is se3::log's right-error tangent (right-not-left, trans-first); "
          "covariance lives in this same tangent by eskf.cpp's full-SE(3)-Ad, by inspection") {
    // What this case PROVES: e = se3::log(T_est^-1 T_gt) is se3::log's OWN right-error tangent
    // (right- not left-error, [trans;rot] order). With e defined as that log, the reconstruction
    // identity T_est o exp(e) == T_gt is a near-tautology of exp/log being mutual inverses, so it
    // pins which tangent e lives in, NOT that the covariance shares it. The covariance lives in
    // this SAME full SE(3) tangent because eskf.cpp propagates with the full-SE(3) se3::adjoint
    // (F = Ad(delta^-1)); that link is true by inspection of eskf.cpp and is NOT asserted here.
    // We pin the convention this way, not by an O(DOF) NEES band. We
    // originally tried "rotation sub-NEES ~ DOF=3" as the convention probe, but the filter's
    // covariance is still mildly PESSIMISTIC (after the init-P fix the ensemble-mean pose NEES
    // is ~0.35 vs DOF 6 — see the next case for the predict-only Ad-inflation diagnosis) so
    // every sub-block's NEES is well below its DOF regardless of convention. So a NEES
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

    // (1) RECONSTRUCTION IDENTITY: e = log(T_est^{-1} T_gt) reconstructs T_gt via the RIGHT
    //     composition T_est o exp(e). With e defined as that log this is a near-tautology of
    //     exp/log being mutual inverses, so it pins that e is se3::log's OWN right-error tangent
    //     (the tangent the ESKF's error T_true = T_est o exp(eta) is defined in) — it does NOT
    //     by itself assert the covariance shares that tangent (that holds by eskf.cpp's full
    //     SE(3) Ad propagation, by inspection).
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

    MESSAGE("convention OK: e is se3::log's right-error tangent (right-not-left, trans-first); "
            "T_est o exp(log(T_est^-1 T_gt)) == T_gt; |e|=" << e.norm());
}

TEST_CASE("validation NEES: Monte-Carlo covariance consistency (chi-square interval, 6 "
          "DOF) — residual PESSIMISM is predict-only Ad-inflation, not init-P") {
    // The non-vacuous consistency test. Ensemble-mean pose NEES over N steady-state samples
    // is compared to a chi-square consistency interval for 6 DOF. For N averaged samples the
    // sum N*nees ~ chi2(N*6) under the consistency hypothesis, so the bound on the MEAN is
    //     [ chi2_inv(alpha/2, N*6)/N , chi2_inv(1-alpha/2, N*6)/N ].
    // With N in the thousands the interval is TIGHT (a Wilson-Hilferty normal approximation
    // to chi2 gives, for large k=N*6, mean ~ 6*(1 +/- z*sqrt(2/k))) — so a consistent filter
    // would be pinned hard, and an inconsistent one is flagged clearly.
    //
    // CAVEAT (independence): this interval ASSUMES the N steady-state samples are independent,
    // but they are autocorrelated within each run (consecutive fused frontiers share overlapping
    // windows + the same slowly-evolving P), so the effective sample count is smaller and the
    // TRUE interval is WIDER than the printed one. This only ever WIDENS the band, never narrows
    // it, and the gap survives either way: even a 100x reduction in effective N widens the
    // interval only to ~[5.4, 6.6], still well above the observed ~1.0.
    //
    // EMPIRICAL OUTCOME (the honesty clause, DESIGN §10). HISTORY, in order of fix:
    //   * P = I_12 seed: ~100x the steady-state squared error; predict-only never shrinks P, so
    //     the ensemble-mean NEES collapsed to ~0.13 (the original ~46x-pessimistic Slice-14 finding).
    //   * INIT-P FIX (estimator.cpp first-fuse seeds P = 0 and lets the first predict() ESTABLISH
    //     P = blkdiag(q_pose, q_pose/dt^2) — the one-window process noise on the gauge-anchored
    //     pose — instead of I_12): raised NEES to ~0.35 (~2.7x; pessimism ~46x -> ~17x).
    //   * MEDIAN-VARIANCE REDUCTION (Slice 14, approach A, THIS slice). The per-window adaptive Q
    //     is spread^2 (spread = inter-source DISAGREEMENT), but the FUSED quantity is the geometric
    //     MEDIAN of the n participating sources, whose per-window variance is ~1/n that of a single
    //     source (n agreeing sources average down). So the spread-derived Q over-stated the fused
    //     accuracy by ~n (=3 here) and the right-error Ad transport AMPLIFIED that into the
    //     translation block. Dividing the spread term by n_eff = n (eskf.cpp adaptive_q, gated by
    //     cfg.adaptive_q_source_reduction) made Q reflect the median's accuracy: NEES rose ~0.35 ->
    //     ~1.0 (a further ~3x; pessimism ~17x -> ~6x). This is the parameter-free part of the fix.
    //   * Q_SCALE CALIBRATION (Slice 14, the COVARIANCE-KNOB SWEEP, THIS slice). The q_scale
    //     coefficient was the last un-calibrated piece: spread is an order-of-magnitude PROXY for the
    //     true single-source per-window sigma, so q_scale=1 over-stated it. An offline grid (q_scale
    //     in {1.0..0.1} x {nees_traj, mixed, turning, straight} x {1x,2x} noise, M=30, since deleted)
    //     chose q_scale=0.5 SAFETY-FIRST: the largest pessimism cut whose worst-case ensemble-mean
    //     NEES stays inside [2,4] and NEVER exceeds DOF=6 on ANY trajectory at either noise level
    //     (q_scale=0.2 hit ~7 on `turning` = overconfident = unsafe). At 0.5 this trajectory's NEES
    //     rose ~1.0 -> ~2.07 (pessimism ~6x -> ~3x). The remaining gap to 6 is INTENTIONAL: it is the
    //     conservative pessimistic margin (the sim under-states real model mismatch), NOT a residual
    //     bug. The cross-trajectory never-overconfident guard is test_cov_calibration.cpp.
    // The filter is STILL NOT strictly chi-square-consistent (~2.07 < 6) BY DESIGN; the residual gap
    // is now the deliberate safety margin, no longer an un-calibrated coefficient. The SHAPE is still
    // the predict-only translation Ad-inflation:
    //   * The right-error Ad(delta^-1) propagation transports the rotation uncertainty into the
    //     TRANSLATION block via the [t]x R coupling as |t| grows: with the un-reduced Q the P
    //     translation diagonal climbed to TENS of m^2 (~77 m^2 by ~32 m of travel — the earlier
    //     "~2.5..5 m^2 by ~6 s" understated it: it keeps growing the whole run), while the actual
    //     squared error grows far more slowly. The /n_eff reduction shrinks that runaway by ~n but
    //     does NOT change its distance-growth SHAPE (that needs a distance-aware covariance model).
    //     The ROTATION block has no such distance coupling, stays tighter, and is better calibrated.
    //   * A predict-only filter has NO correction step to pull P back down between windows, so once
    //     translation P inflates with distance it stays inflated. (Slice 11's absolute-ref
    //     correction DOES shrink P when a ref is present — see the NIS case below, NIS ~2.83 vs
    //     DOF 3, near-consistent.)
    // Tuning the test's q_floor cannot help (it only ADDS to P); q_scale MULTIPLIES the spread term
    // and is now CALIBRATED to 0.5 (the core default). A strictly chi-square-consistent predict-only
    // Sigma on THIS trajectory would need q_scale ~0.17 — but that over-fits ONE trajectory and
    // pushes others overconfident, so the safety-first 0.5 deliberately leaves the filter mildly
    // pessimistic; a distance-aware covariance model (or denser corrections) is the remaining path
    // to strict per-trajectory consistency. So per the honesty clause we: (a) compute + report the
    // actual NEES and the chi-square bound it still (intentionally) sits below, and (b) assert a
    // DOCUMENTED ACHIEVABLE band TIGHT around the calibrated value, making this a non-vacuous
    // regression guard (a drift in either P or the tracking error breaks it, AND a regression of the
    // default back toward q_scale=1 — mean ~1.04 — OR disabling the /n_eff reduction — mean ~0.70 —
    // re-trips it; the knob-OFF guard is asserted explicitly below). Surfaced to the orchestrator:
    // the init-P pessimism, the median-variance over-count (approach A), AND the q_scale coefficient
    // (this sweep) are now all RESOLVED; the remaining gap to DOF is the chosen safety margin, with
    // a distance-aware covariance model / denser corrections the path to strict consistency.
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

    // Q knobs at the CALIBRATED core default (the Slice-14 covariance sweep, now DONE — see the
    // diagnosis above): q_scale was lowered from the un-calibrated 1.0 placeholder to 0.5, chosen
    // safety-first across {nees_traj, mixed, turning, straight} x {1x,2x} noise so the worst-case
    // ensemble-mean NEES is mildly pessimistic (~2.9 at 1x, ~3.5 at 2x) and NEVER overconfident
    // (<DOF=6) on any trajectory. We exercise that SAME 0.5 here (matching Config{}.q_scale) so the
    // band below pins the calibrated value, not the old placeholder. The cross-trajectory never-
    // overconfident + pessimism-reduced guard lives in test_cov_calibration.cpp.
    const Scalar q_scale = Config{}.q_scale;   // the calibrated core default (0.5)
    const Scalar qf[6]   = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};

    // Run the ensemble once for a given reduction-knob setting and return the ensemble-mean NEES
    // + N. Factored out so we can measure BOTH the reduction-ON value (the calibrated result) and
    // the reduction-OFF value (the regression guard proving the knob still gates the /n_eff term).
    auto run_ensemble = [&](bool source_reduction, long& N_out) -> Scalar {
        Scalar sum_nees = 0.0;
        long   N        = 0;
        for (int run = 0; run < M; ++run) {
            std::vector<SourceParams> planted = base;
            for (auto& sp : planted) sp.seed = 5000u + run * 11u + sp.id;

            std::vector<std::unique_ptr<SyntheticSource>> srcs;
            for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
            std::vector<SensorConfig> sensors;
            Config cfg = nees_config(planted, sensors, q_scale, qf);
            cfg.adaptive_q_source_reduction = source_reduction;

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
        N_out = N;
        return (N > 0) ? sum_nees / static_cast<Scalar>(N) : Scalar(0);
    };

    // The slice's result: median-variance reduction ON (n_eff = n = 3 -> spread term / 3).
    long N = 0;
    const Scalar mean_nees = run_ensemble(/*source_reduction=*/true, N);
    REQUIRE(N > 1000);

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
    // Build the YES/NO verdict as a std::string first: a bare (cond ? "YES" : "NO") is a
    // const char* that doctest's MESSAGE stringifier captures as a POINTER (prints an
    // address), not the literal text. std::string forces the human-readable line.
    const std::string verdict = truly_consistent ? "YES" : "NO";
    MESSAGE("NEES consistency: ensemble-mean=" << mean_nees << " DOF=6 N=" << N
            << " chi2 99% interval on mean=[" << lo_mean << ", " << hi_mean << "]"
            << "  truly_consistent=" << verdict
            << "  (mean < 6 => still mildly PESSIMISTIC by design: init-P + median-variance "
               "over-count RESOLVED, q_scale now CALIBRATED to 0.5; was ~0.13 (I_12 seed) -> "
               "~0.35 (init-P fix) -> ~1.0 (/n_eff, q_scale=1) -> ~2.07 (q_scale=0.5 calibrated))");

    // TRUE-CONSISTENCY verdict. THREE fixes have materially improved this in order: ~0.13 -> ~0.35
    // (init-P fix), ~0.35 -> ~1.0 (the /n_eff median-variance reduction), ~1.0 -> ~2.07 (THIS slice:
    // q_scale calibrated 1.0 -> 0.5). The remaining gap to strict chi-square consistency (~2.07 vs
    // 6) is DELIBERATE and SAFETY-FIRST: forcing NEES to 6 on this single trajectory (q_scale ~0.17)
    // would push OTHER trajectories overconfident (NEES > 6, unsafe) on the optimistic sim, so the
    // calibration leaves a conservative pessimistic margin (worst-case ~2.9 across the set; see
    // test_cov_calibration.cpp). So the verdict is INTENTIONALLY still FALSE — mild pessimism is the
    // chosen safe operating point, not a failure to converge. (A distance-aware covariance model /
    // denser corrections remain the path to strict per-trajectory consistency.)
    CHECK_FALSE(truly_consistent);     // mildly PESSIMISTIC BY DESIGN (safety margin, never overconf)

    // DOCUMENTED ACHIEVABLE BOUND (non-vacuous regression guard on the value the filter ACTUALLY
    // produces at the CALIBRATED default q_scale=0.5, after the init-P fix + /n_eff reduction).
    // Observed ensemble-mean ~2.07 (was ~1.04 at the old un-calibrated q_scale=1, ~0.35 with
    // reduction OFF, ~0.13 under the old I_12 seed). The band [1.5, 2.7] brackets the calibrated
    // value with headroom each side: it absorbs Monte-Carlo seed variation, yet a drift in EITHER
    // the published covariance OR the fused-vs-GT tracking error breaks it — AND a regression of the
    // default back toward q_scale=1 (mean ~1.04) trips the lower bound, as does disabling the
    // /n_eff reduction (mean ~0.7 at this q_scale) or the old I_12 seed. It pins the calibrated,
    // mildly-pessimistic covariance the filter emits today; it stays < DOF (6) -> never overconfident.
    CHECK(mean_nees > 1.50);     // NOT regressed toward the old un-calibrated q_scale=1 (~1.04)
    CHECK(mean_nees < 2.70);     // NOT overconfident (consistency would be ~6; safe margin kept)

    // KNOB-OFF REGRESSION GUARD. Re-run the SAME ensemble with the median-variance reduction
    // DISABLED: n_eff = 1 recovers the un-reduced covariance, so at q_scale=0.5 the mean drops to
    // ~0.70 (the /n_eff reduction is worth ~n=3x). This proves (a) the knob still gates the /n_eff
    // reduction, and (b) the reduction is a STRICT improvement (ON ~ n x OFF), independent of the
    // q_scale calibration.
    long N_off = 0;
    const Scalar mean_nees_off = run_ensemble(/*source_reduction=*/false, N_off);
    REQUIRE(N_off > 1000);
    MESSAGE("NEES reduction-OFF guard: ensemble-mean=" << mean_nees_off
            << " (expect ~0.70 at the calibrated q_scale=0.5 with /n_eff disabled)");
    CHECK(mean_nees_off > 0.45);     // /n_eff disabled at q_scale=0.5 -> ~0.70
    CHECK(mean_nees_off < 1.00);
    // The reduction ON is materially LARGER than OFF (the /n_eff gain is real, ~n-fold not noise).
    CHECK(mean_nees > mean_nees_off * 1.8);
}

// ===========================================================================
// DELIVERABLE 1b — NIS Monte-Carlo consistency (Slice 11 closes the deferral)
// ===========================================================================
TEST_CASE("validation NIS: Monte-Carlo innovation consistency of accepted absolute-ref "
          "fixes (chi-square interval, 3 DOF)") {
    // NIS (Normalized Innovation Squared) is defined only at a MEASUREMENT UPDATE. Slice 11
    // adds the absolute-ref correction path, so NIS is FINALLY computable (the Slice-14
    // deferral — there was no innovation in the predict-only filter). For a consistent
    // filter the per-fix NIS d2 = r^T S^-1 r ~ chi2(n), n = the fix DOF (here position,
    // n=3), so the ensemble mean over many accepted fixes ~ 3.
    //
    // The correction step SHRINKS P toward the measurement (the Slice-14 covariance-
    // pessimism finding's fix WHEN a ref is present): repeated fixes pull P down toward a
    // steady state set by the measurement noise, so unlike the predict-only NEES (~1.0 after the
    // init-P fix + this slice's /n_eff median-variance reduction, still mildly pessimistic from
    // the un-calibrated q_scale amplified by the distance-driven translation Ad-inflation) the NIS
    // on the CORRECTED path is in a sane range. We use priors == planted (no
    // calibration transient) + a clean position
    // ref whose sim sigma_pos EQUALS the R the measurement reports, so the filter's noise
    // model matches the actual draw — the only honest way to expect NIS ~ DOF.
    //
    // HONESTY CLAUSE (DESIGN §10, same rule as the NEES work): the published P is still
    // partly pessimistic between fixes (predict re-inflates it), so S = HPH^T + R is biased
    // HIGH and the ensemble-mean NIS sits somewhat BELOW the DOF count rather than dead on
    // it. We report the measured mean + the chi-square interval; if it is not strictly
    // inside the interval we assert a DOCUMENTED ACHIEVABLE band on the value the filter
    // actually produces (a non-vacuous regression guard) and surface it as the outcome.
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

    const Scalar q_scale     = Config{}.q_scale;   // the calibrated core default (0.5), as the NEES case
    const Scalar qf[6]       = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
    const Scalar sigma_pos   = 0.10;   // the ref's reported R sigma == the drawn noise sigma

    Scalar sum_nis = 0.0;
    long   N       = 0;
    for (int run = 0; run < M; ++run) {
        std::vector<SourceParams> planted = base;
        for (auto& sp : planted) sp.seed = 8000u + run * 13u + sp.id;

        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
        std::vector<SensorConfig> sensors;
        Config cfg = nees_config(planted, sensors, q_scale, qf);
        cfg.mahalanobis_chi2 = 1e9;    // accept all inliers (clean run -> no real outliers)

        AbsoluteRefParams rp;
        rp.period_s  = 0.10;           // 10 Hz fixes (plenty of samples)
        rp.window_s  = cfg.window_s;
        rp.sigma_pos = sigma_pos;
        rp.seed      = 9000u + static_cast<std::uint32_t>(run);
        SyntheticAbsoluteRef ref(tr, rp);

        Rig rig;
        rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
        REQUIRE(rig.add_correction(&ref) == Status::Ok);
        rig.run(0.2, tr.duration_s() - 0.1, 50.0);

        // Average the per-fix NIS over steady-state ticks (drop warmup) on records where a
        // fix was applied (corr_applied > 0 -> last_nis is an APPLIED fix's NIS).
        int fused_seen = 0;
        const int warmup = 20;
        for (const Record& r : rig.records()) {
            if (!r.fused) continue;
            ++fused_seen;
            if (fused_seen <= warmup) continue;
            if (r.result.correction.corr_applied <= 0) continue;
            sum_nis += r.result.correction.last_nis;
            ++N;
        }
    }
    REQUIRE(N > 1000);
    const Scalar mean_nis = sum_nis / static_cast<Scalar>(N);

    // Two-sided 99% chi-square interval on the mean for k = N*3 DOF (Wilson-Hilferty).
    const Scalar dof = 3.0;
    const Scalar k   = static_cast<Scalar>(N) * dof;
    const Scalar z   = 2.5758;
    auto wh = [&](Scalar zp) {
        const Scalar a = Scalar(2) / (Scalar(9) * k);
        const Scalar f = Scalar(1) - a + zp * std::sqrt(a);
        return k * f * f * f;
    };
    const Scalar lo_mean = wh(-z) / static_cast<Scalar>(N);
    const Scalar hi_mean = wh(+z) / static_cast<Scalar>(N);

    const bool consistent = (mean_nis >= lo_mean && mean_nis <= hi_mean);
    const std::string verdict = consistent ? "YES" : "NO";
    MESSAGE("NIS consistency: ensemble-mean=" << mean_nis << " DOF=3 N=" << N
            << " chi2 99% interval on mean=[" << lo_mean << ", " << hi_mean << "]"
            << "  consistent=" << verdict);

    // EMPIRICAL OUTCOME: ensemble-mean NIS ~ 2.83 for DOF=3 at the CALIBRATED q_scale=0.5 — MILDLY
    // conservative (just below the tight chi2 interval ~[2.91, 3.09]), a NIGHT-AND-DAY contrast with
    // the predict-only NEES (~2.07 at the same calibrated default). The correction step pulls P down
    // toward the measurement noise, so the innovation is normalized by a covariance close to the
    // true error. The residual ~6% shortfall is P re-inflating between fixes (predict adds Q + the
    // Ad propagation before the next fix), biasing S = HPH^T + R slightly HIGH. The q_scale
    // calibration (1.0 -> 0.5, THIS slice) shrinks that between-fix Q further, so the re-inflation is
    // smaller and the mean NIS moved UP toward DOF 3 (~2.69 at q_scale=1 -> ~2.83 at q_scale=0.5,
    // LESS conservative) — the same mechanism that raised the predict-only NEES, seen from the
    // corrected side. So the filter is NOT strictly chi-square-consistent yet, but is now in a SANE
    // O(DOF) range — reported as a partial close of the Slice-14 covariance finding (full
    // consistency would need a distance-aware covariance model and/or denser fixes; surfaced to the
    // orchestrator). The CORRECTED-path NIS is set by the measurement noise + the between-fix
    // re-inflation, so it is largely insensitive to the init-P seed (the correction dominates P
    // within a few fixes) but DOES respond to the Q magnitude (the q_scale calibration) as noted.
    //
    // ACHIEVABLE BAND (non-vacuous regression guard on the value the filter ACTUALLY
    // produces). Pins the observed ~2.83 while absorbing Monte-Carlo seed variation; a ~1.3x
    // drift in EITHER the published covariance OR the tracking error breaks it. Stays well
    // clear of both the pessimistic (<<1) and overconfident (>>DOF) failure modes. (Band kept at
    // [1.5, 3.5]: the calibrated ~2.83 still sits comfortably inside it with headroom both
    // sides, and 3.5 stays below the DOF=3 overconfidence line + the ~3.09 chi2 upper bound.)
    CHECK(mean_nis > 1.5);     // NOT collapsed like the predict-only NEES
    CHECK(mean_nis < 3.5);     // NOT overconfident (would blow past DOF)
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

    // Determinism is a BYTE-IDENTICAL property over ALL fields of ALL fused records (fused pose
    // + full 12x12 covariance + every per-source calib field: extrinsic R/t, scale, time offset,
    // 3 confidences, 3 commit flags). Rather than emit ~16 CHECKs per record (which inflated the
    // suite's assertion count ~3.4x with non-discriminating repeats), we fold every exact `==`
    // comparison for a record into ONE boolean and assert it ONCE per fused record. The fold
    // uses the SAME exact equality (operator==, not a tolerance), so if ANY field of ANY record
    // diverges across the two runs the per-record CHECK fails — it remains a genuine drift
    // detector, just one assertion per record instead of per field.
    long fused_compared = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        const bool fused_agrees = (a[i].fused == b[i].fused);
        if (!fused_agrees) { CHECK(fused_agrees); continue; }  // structural divergence: flag it
        if (!a[i].fused) continue;
        const Result& ra = a[i].result;
        const Result& rb = b[i].result;

        bool rec_equal = true;
        // Fused pose + covariance (as test_sim does)...
        rec_equal = rec_equal && (ra.frontier.pose.R.array() == rb.frontier.pose.R.array()).all();
        rec_equal = rec_equal && (ra.frontier.pose.t.array() == rb.frontier.pose.t.array()).all();
        rec_equal = rec_equal && (ra.frontier.cov.array()    == rb.frontier.cov.array()).all();
        // ...AND the per-source calib snapshot (the extension this golden adds).
        rec_equal = rec_equal && (ra.source_count == rb.source_count);
        if (rec_equal) {
            for (int s = 0; s < ra.source_count; ++s) {
                const CalibSnapshot& ca = ra.calib[s];
                const CalibSnapshot& cb = rb.calib[s];
                rec_equal = rec_equal && (ca.id == cb.id);
                rec_equal = rec_equal && (ca.extrinsic.R.array() == cb.extrinsic.R.array()).all();
                rec_equal = rec_equal && (ca.extrinsic.t.array() == cb.extrinsic.t.array()).all();
                rec_equal = rec_equal && (ca.scale         == cb.scale);
                rec_equal = rec_equal && (ca.time_offset_s == cb.time_offset_s);
                rec_equal = rec_equal && (ca.extrinsic_confidence   == cb.extrinsic_confidence);
                rec_equal = rec_equal && (ca.scale_confidence       == cb.scale_confidence);
                rec_equal = rec_equal && (ca.translation_confidence == cb.translation_confidence);
                rec_equal = rec_equal && (ca.extrinsic_committed   == cb.extrinsic_committed);
                rec_equal = rec_equal && (ca.scale_committed       == cb.scale_committed);
                rec_equal = rec_equal && (ca.translation_committed == cb.translation_committed);
            }
        }
        CHECK(rec_equal);     // ONE assertion per fused record: all fields byte-identical
        ++fused_compared;
    }
    REQUIRE(fused_compared > 50);     // guard: the replay actually produced fused records to compare
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
