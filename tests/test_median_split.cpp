// Slice 19 unit tests: per-channel split median (D3 amendment) — solve_split().
//
// Coverage (SLICE19 §2 items 1-2, solver level):
//   * HEADLINE — a rotation-channel outlier is rejected by the rotation median while that
//     SAME source's translation fully participates (and vice versa) — the new capability
//     the coupled solver cannot express (pinned by direct contrast against solve()).
//   * n == 0 / 1 / 2 closed forms (identity / passthrough / per-channel weighted
//     interpolation with DIFFERENT channel weights).
//   * D3 safeguards mirrored per channel: high-weight-outlier guard (the old vertex-init
//     pinning blind spot), interior-not-pinned, VZ coincident-vertex guard, WCET caps.
//   * Per-channel spreads (rad / m, no lambda unit-mixing).
//   * Cross-channel veto: a gross both-channel fault is rejected as completely as the
//     coupled solver rejects it (parity bounds); veto OFF leaves the clean channel's
//     contribution; the final-weight outputs expose the veto scaling.
//   * Determinism.
#include <doctest/doctest.h>

#include "ofc/core/lie.hpp"
#include "ofc/core/median.hpp"

#include <cmath>

using namespace ofc;

namespace {

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

median::Params params() {
    median::Params p;
    p.max_iters = 100;
    p.tol       = 1e-6;     // the production weiszfeld_tol (see test_median.cpp on why
                            // tighter step demands are unreachable near a vertex)
    p.eps       = 1e-9;
    p.lambda    = 1.0;
    return p;
}

SE3 make_se3(const Vec3& w, const Vec3& t) {
    SE3 T;
    T.R = so3::exp(w);
    T.t = t;
    return T;
}

Scalar geo_dist(const Mat3& a, const Mat3& b) { return so3::log(a.transpose() * b).norm(); }

} // namespace

// ---------------------------------------------------------------------------
// Degenerate counts / closed forms
// ---------------------------------------------------------------------------
TEST_CASE("split median: zero inputs -> identity, zero spreads") {
    const median::SplitResult r = median::solve_split(nullptr, nullptr, nullptr, 0, params());
    CHECK(close(r.value.R, Mat3::Identity(), 1e-12));
    CHECK(close(r.value.t, Vec3::Zero(), 1e-12));
    CHECK(r.spread_rot == doctest::Approx(0.0));
    CHECK(r.spread_trans == doctest::Approx(0.0));
    CHECK(r.converged_rot);
    CHECK(r.converged_trans);
}

TEST_CASE("split median: one input is a passthrough in both channels") {
    const SE3 x = make_se3(Vec3(0.1, -0.2, 0.3), Vec3(1.0, 2.0, -1.0));
    const Scalar w = 1.0;
    const median::SplitResult r = median::solve_split(&x, &w, &w, 1, params());
    CHECK(close(r.value.R, x.R, 1e-12));
    CHECK(close(r.value.t, x.t, 1e-12));
    CHECK(r.spread_rot == doctest::Approx(0.0));
    CHECK(r.spread_trans == doctest::Approx(0.0));
}

TEST_CASE("split median: n==2 interpolates PER CHANNEL with that channel's weights") {
    // The split solver's n==2 closed form honors DIFFERENT weights per channel — the
    // coupled solver can only interpolate both channels at one shared fraction.
    constexpr Scalar kPi = 3.14159265358979323846;
    const SE3 a = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    const SE3 b = make_se3(Vec3(0, 0, kPi / 2), Vec3(4.0, 0.0, 0.0));
    const SE3 xs[2] = {a, b};
    const Scalar w_rot[2]   = {1.0, 3.0};   // rotation pulled 3/4 toward b
    const Scalar w_trans[2] = {3.0, 1.0};   // translation pulled only 1/4 toward b

    const median::SplitResult r = median::solve_split(xs, w_rot, w_trans, 2, params());
    const Mat3 expect_R = so3::exp(Vec3(0, 0, kPi * 3.0 / 8.0));   // u_rot = 3/4 of pi/2
    CHECK(close(r.value.R, expect_R, 1e-9));
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(1e-9));    // u_trans = 1/4 of 4.0
    CHECK(r.converged_rot);
    CHECK(r.converged_trans);
}

TEST_CASE("split median: all-<=0 weights fall back to uniform PER CHANNEL") {
    SE3 xs[3];
    xs[0] = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.0), Vec3(1.0, 0.0, 0.0));
    xs[2] = make_se3(Vec3(0, 0, 0.0), Vec3(2.0, 0.0, 0.0));
    const Scalar w_zero[3] = {0.0, 0.0, 0.0};
    const Scalar w_unit[3] = {1.0, 1.0, 1.0};

    // Rotation weights all-zero, translation weights valid: the rotation channel uses the
    // uniform fallback (well-defined median), the translation channel the given weights.
    const median::SplitResult r_zero = median::solve_split(xs, w_zero, w_unit, 3, params());
    const median::SplitResult r_unit = median::solve_split(xs, w_unit, w_unit, 3, params());
    CHECK(close(r_zero.value.R, r_unit.value.R, 1e-9));
    CHECK(close(r_zero.value.t, r_unit.value.t, 1e-9));
    // Spread is reported with the effective (uniform) weights — not understated to 0.
    CHECK(r_zero.spread_trans == doctest::Approx(r_unit.spread_trans).epsilon(1e-9));
    CHECK(r_zero.spread_trans > 0.1);
}

// ---------------------------------------------------------------------------
// HEADLINE (SLICE19 §2 item 1): a rotation-channel outlier is rejected while the
// SAME source's translation FULLY participates — impossible for the coupled solver.
// ---------------------------------------------------------------------------
TEST_CASE("split median HEADLINE: rotation outlier rejected, same source's translation "
          "fully participates (coupled solver pinned as unable to do both)") {
    // Source 3: GROSS rotation fault (~1.6 rad off the inlier cluster) but the CORRECT
    // translation (1.0, 0, 0), carried at translation weight 4 (> the inlier mass 3).
    // Sources 0-2: clean rotations ~ yaw 0.30, translations biased to ~1.4 m.
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.00, 0.00, 0.30), Vec3(1.40, 0.00, 0.0));
    xs[1] = make_se3(Vec3(0.01, 0.00, 0.31), Vec3(1.42, 0.01, 0.0));
    xs[2] = make_se3(Vec3(0.00, 0.01, 0.29), Vec3(1.38, -0.01, 0.0));
    xs[3] = make_se3(Vec3(1.20, -0.80, 0.60), Vec3(1.00, 0.00, 0.0));  // rot FAULT, best t
    const Scalar w_rot[4]   = {1, 1, 1, 1};
    const Scalar w_trans[4] = {1, 1, 1, 4};

    median::Params p = params();

    // --- The split solver achieves BOTH (veto OFF: the per-channel-grace policy fork;
    //     the veto-ON companion below documents the whole-source-rejection default).
    Scalar wr_fin[4], wt_fin[4];
    const median::SplitResult sr =
        median::solve_split(xs, w_rot, w_trans, 4, p, /*veto=*/false, wr_fin, wt_fin);
    // Rotation channel: the gross rotation fault is REJECTED (consensus on the inliers).
    const Vec3 log_r = so3::log(sr.value.R);
    CHECK(std::abs(log_r.z() - 0.30) < 0.05);
    CHECK(std::abs(log_r.x()) < 0.05);
    CHECK(std::abs(log_r.y()) < 0.05);
    // Translation channel: the SAME source's translation FULLY participates — its weight
    // is untouched and its (dominant-mass) value carries the consensus to (1.0, 0, 0).
    CHECK(wt_fin[3] == doctest::Approx(4.0));
    CHECK(sr.value.t.x() == doctest::Approx(1.0).epsilon(0.05));
    CHECK(std::abs(sr.value.t.y()) < 0.05);

    // --- The coupled solver CANNOT do both (the "impossible today" pin). One scalar
    //     weight per source forces a choice:
    //     (a) weight the source high -> its rotation mass (4 > 3) drags/pins the coupled
    //         rotation onto the FAULT;
    {
        const Scalar w_high[4] = {1, 1, 1, 4};
        const median::Result cr = median::solve(xs, w_high, 4, p);
        CHECK(geo_dist(cr.value.R, so3::exp(Vec3(0, 0, 0.30))) > 0.5);   // rotation lost
    }
    //     (b) weight it equal -> the gross rotation inflates its split-distance, its
    //         translation is discarded with the rest of the source, and the fused
    //         translation stays on the biased inlier cluster (~1.4), NOT 1.0.
    {
        const Scalar w_eq[4] = {1, 1, 1, 1};
        const median::Result cr = median::solve(xs, w_eq, 4, p);
        CHECK(std::abs(cr.value.t.x() - 1.0) > 0.2);                     // translation lost
    }

    // --- Veto-ON companion (the policy fork made visible): the gross rotation fault
    //     trips the cross-channel veto, which scales the source's translation weight by
    //     kVetoWeightScale — whole-source rejection (the safe default for hard faults),
    //     deliberately giving up this scenario's cross-channel grace.
    {
        Scalar wr2[4], wt2[4];
        const median::SplitResult sv =
            median::solve_split(xs, w_rot, w_trans, 4, p, /*veto=*/true, wr2, wt2);
        CHECK(wt2[3] == doctest::Approx(4.0 * median::kVetoWeightScale));
        // The fused translation falls back toward the inlier cluster (no longer 1.0).
        CHECK(sv.value.t.x() > 1.2);
    }
}

TEST_CASE("split median HEADLINE (vice versa): translation outlier rejected, same source's "
          "rotation fully participates") {
    // Source 3: GROSS translation fault but the CORRECT rotation (yaw 0.30), carried at
    // rotation weight 4. Sources 0-2: clean translations ~ (1, 0, 0), rotations biased to
    // yaw ~0.45.
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.00, 0.00, 0.45), Vec3(1.00, 0.00, 0.0));
    xs[1] = make_se3(Vec3(0.01, 0.00, 0.46), Vec3(1.02, 0.01, 0.0));
    xs[2] = make_se3(Vec3(0.00, 0.01, 0.44), Vec3(0.98, -0.01, 0.0));
    xs[3] = make_se3(Vec3(0.00, 0.00, 0.30), Vec3(50.0, -40.0, 7.0));  // trans FAULT, best R
    const Scalar w_rot[4]   = {1, 1, 1, 4};
    const Scalar w_trans[4] = {1, 1, 1, 1};

    median::Params p = params();
    Scalar wr_fin[4], wt_fin[4];
    const median::SplitResult sr =
        median::solve_split(xs, w_rot, w_trans, 4, p, /*veto=*/false, wr_fin, wt_fin);
    // Translation channel: the gross translation fault is rejected.
    CHECK(sr.value.t.x() == doctest::Approx(1.0).epsilon(0.05));
    CHECK(std::abs(sr.value.t.y()) < 0.1);
    CHECK(std::abs(sr.value.t.z()) < 0.1);
    // Rotation channel: the SAME source's rotation fully participates (weight untouched,
    // dominant mass carries the consensus rotation to its yaw 0.30).
    CHECK(wr_fin[3] == doctest::Approx(4.0));
    const Vec3 log_r = so3::log(sr.value.R);
    CHECK(std::abs(log_r.z() - 0.30) < 0.05);
}

// ---------------------------------------------------------------------------
// D3 safeguards mirrored per channel (port of the coupled solver's pinning tests)
// ---------------------------------------------------------------------------
TEST_CASE("split median: translation channel rejects a HIGHEST-WEIGHT outlier (D3 blind "
          "spot, per channel)") {
    // Four translation inliers ~ (1,0,0) + one gross outlier carrying the LARGEST single
    // weight (2.0 < the inlier mass 4.0). The OLD vertex-init pinning bug would return the
    // outlier verbatim; the off-vertex interior init + 1/d reweight must reject it.
    SE3 xs[5];
    xs[0] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.00,  0.00, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.01),  Vec3(1.02,  0.01, 0.0));
    xs[2] = make_se3(Vec3(0, 0, -0.01), Vec3(0.98, -0.01, 0.0));
    xs[3] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.01,  0.00, 0.0));
    xs[4] = make_se3(Vec3(0, 0, 0.02),  Vec3(50.0, -40.0, 7.0));   // trans OUTLIER, w 2.0
    const Scalar w[5] = {1, 1, 1, 1, 2.0};

    const median::SplitResult r = median::solve_split(xs, w, w, 5, params(), /*veto=*/false);
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(0.1));
    CHECK(std::abs(r.value.t.y()) < 0.2);
    CHECK(std::abs(r.value.t.z()) < 0.2);
    // INTERIOR-NOT-PINNED: the channel genuinely iterated and the consensus is far off the
    // high-weight outlier vertex.
    CHECK(r.iters_trans > 1);
    CHECK((r.value.t - xs[4].t).norm() > 1.0);
}

TEST_CASE("split median: rotation channel rejects a HIGHEST-WEIGHT outlier (D3 blind spot, "
          "per channel)") {
    SE3 xs[5];
    xs[0] = make_se3(Vec3(0.00, 0.00, 0.30), Vec3(1.0, 0, 0));
    xs[1] = make_se3(Vec3(0.01, 0.00, 0.31), Vec3(1.0, 0, 0));
    xs[2] = make_se3(Vec3(0.00, 0.01, 0.29), Vec3(1.0, 0, 0));
    xs[3] = make_se3(Vec3(-0.01, 0.00, 0.30), Vec3(1.0, 0, 0));
    xs[4] = make_se3(Vec3(1.20, -0.80, 0.60), Vec3(1.0, 0, 0));   // rot OUTLIER, w 2.0
    const Scalar w[5] = {1, 1, 1, 1, 2.0};

    const median::SplitResult r = median::solve_split(xs, w, w, 5, params(), /*veto=*/false);
    const Vec3 log_r = so3::log(r.value.R);
    CHECK(std::abs(log_r.z() - 0.30) < 0.05);
    CHECK(std::abs(log_r.x()) < 0.05);
    CHECK(std::abs(log_r.y()) < 0.05);
    CHECK(r.iters_rot > 1);
    CHECK(geo_dist(r.value.R, xs[4].R) > 1.0);
}

TEST_CASE("split median: VZ coincident-vertex guard per channel (no NaN, weight-dominant "
          "coincident pair wins)") {
    // Two EXACTLY coincident inputs jointly weight-dominant in both channels + one off
    // input: the channel medians sit on the coincident pair; the VZ d<=eps skip keeps the
    // iteration finite (no w/eps self-weight blow-up).
    SE3 xs[3];
    xs[0] = make_se3(Vec3(0, 0, 0.20), Vec3(1.0, 0.0, 0.0));
    xs[1] = xs[0];                                          // exact coincidence
    xs[2] = make_se3(Vec3(0, 0, -0.40), Vec3(-2.0, 1.0, 0.0));
    const Scalar w[3] = {1.0, 1.0, 0.5};

    const median::SplitResult r = median::solve_split(xs, w, w, 3, params(), /*veto=*/false);
    CHECK(r.value.t.allFinite());
    CHECK(r.value.R.allFinite());
    CHECK(close(r.value.t, xs[0].t, 1e-3));
    CHECK(geo_dist(r.value.R, xs[0].R) < 1e-3);
    CHECK(std::isfinite(r.spread_rot));
    CHECK(std::isfinite(r.spread_trans));
}

TEST_CASE("split median WCET: production caps (max_iters=10), veto ON — bounded iterations, "
          "well-formed value") {
    // At the production cap the split path runs at most 2 base solves + 2 veto re-solves,
    // each <= max_iters — iters_* <= 2*max_iters by construction. The value is already
    // decisively on the inlier cluster inside the cap (the D3 WCET pin, per channel).
    SE3 xs[5];
    xs[0] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.00,  0.00, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.01),  Vec3(1.02,  0.01, 0.0));
    xs[2] = make_se3(Vec3(0, 0, -0.01), Vec3(0.98, -0.01, 0.0));
    xs[3] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.01,  0.00, 0.0));
    xs[4] = make_se3(Vec3(1.20, -0.80, 0.60), Vec3(50.0, -40.0, 7.0));  // both-channel fault
    const Scalar w[5] = {1, 1, 1, 1, 2.0};

    median::Params p = params();
    p.max_iters = 10;       // production weiszfeld_max_iters (the WCET cap)
    const median::SplitResult r = median::solve_split(xs, w, w, 5, p, /*veto=*/true);
    CHECK(r.iters_rot <= 2 * p.max_iters);
    CHECK(r.iters_trans <= 2 * p.max_iters);
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(0.1));
    CHECK(std::abs(r.value.t.y()) < 0.2);
    CHECK(so3::log(r.value.R).norm() < 0.2);
}

// ---------------------------------------------------------------------------
// Per-channel spreads (rad / m — no lambda unit-mixing)
// ---------------------------------------------------------------------------
TEST_CASE("split median: spreads are per-channel — rotation disagreement does not leak "
          "into the translation spread (and vice versa)") {
    // Rotation-only disagreement: spread_rot > 0, spread_trans == 0.
    SE3 rot_only[3];
    rot_only[0] = make_se3(Vec3(0, 0, 0.00), Vec3(1.0, 0, 0));
    rot_only[1] = make_se3(Vec3(0, 0, 0.30), Vec3(1.0, 0, 0));
    rot_only[2] = make_se3(Vec3(0, 0, -0.30), Vec3(1.0, 0, 0));
    const Scalar w[3] = {1, 1, 1};
    const median::SplitResult rr = median::solve_split(rot_only, w, w, 3, params(),
                                                       /*veto=*/false);
    CHECK(rr.spread_rot > 0.1);
    CHECK(rr.spread_trans == doctest::Approx(0.0).epsilon(1e-9));

    // Translation-only disagreement: spread_trans > 0 (meters), spread_rot == 0.
    SE3 trans_only[3];
    trans_only[0] = make_se3(Vec3(0, 0, 0.2), Vec3(0.0, 0, 0));
    trans_only[1] = make_se3(Vec3(0, 0, 0.2), Vec3(1.0, 0, 0));
    trans_only[2] = make_se3(Vec3(0, 0, 0.2), Vec3(-1.0, 0, 0));
    const median::SplitResult rt = median::solve_split(trans_only, w, w, 3, params(),
                                                       /*veto=*/false);
    CHECK(rt.spread_trans > 0.1);
    CHECK(rt.spread_rot == doctest::Approx(0.0).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// Cross-channel veto (SLICE19 §2 item 2)
// ---------------------------------------------------------------------------
namespace {
// The veto parity rig: 3 inliers with REAL translation scatter (~0.1 m) + one HIGH-WEIGHT
// source faulty in BOTH channels — rotation grossly (1.9 rad), translation well outside
// the scatter (0.9 m). The per-channel 1/d reweight alone leaves a measurable translation
// drag (the fault's weight 2.0 vs the scattered inlier mass 3.0); the coupled solver's
// shared distance — dominated by the gross rotation — rejects the WHOLE source. The veto
// must recover that whole-source rejection.
void veto_parity_rig(SE3 xs[4], Scalar w[4]) {
    xs[0] = make_se3(Vec3(0.00, 0.00, 0.30), Vec3(1.00, 0.00, 0.0));
    xs[1] = make_se3(Vec3(0.01, 0.00, 0.31), Vec3(1.10, 0.08, 0.0));
    xs[2] = make_se3(Vec3(0.00, 0.01, 0.29), Vec3(0.90, -0.08, 0.0));
    xs[3] = make_se3(Vec3(1.20, -0.80, 0.60), Vec3(1.80, 0.50, 0.0));  // BOTH-channel fault
    w[0] = 1.0; w[1] = 1.0; w[2] = 1.0; w[3] = 2.0;
}
} // namespace

TEST_CASE("split median veto: a gross BOTH-channel fault is rejected as completely as the "
          "coupled solver rejects it (parity); veto OFF leaves the clean-channel drag") {
    SE3 xs[4];
    Scalar w[4];
    veto_parity_rig(xs, w);
    median::Params p = params();

    // The truth the inliers agree on (their cluster center / rotation).
    const Vec3 t_true(1.0, 0.0, 0.0);
    const Mat3 R_true = so3::exp(Vec3(0, 0, 0.30));

    // The coupled solver's whole-source rejection (the parity REFERENCE): the fault's
    // gross rotation inflates its split-distance, so its translation pull collapses too —
    // a bounded residual drag remains (the high weight against real inlier scatter).
    const median::Result coupled = median::solve(xs, w, 4, p);
    const Scalar coupled_terr = (coupled.value.t - t_true).norm();
    const Scalar coupled_rerr = geo_dist(coupled.value.R, R_true);
    CHECK(coupled_terr < 0.1);                           // the coupled outlier-test bound
    CHECK(coupled_rerr < 0.1);

    // Veto ON (the split-path default): PARITY — the fault is rejected AT LEAST as
    // completely as the coupled solver rejects it (the veto's whole-source rejection).
    Scalar wr_on[4], wt_on[4];
    const median::SplitResult v_on =
        median::solve_split(xs, w, w, 4, p, /*veto=*/true, wr_on, wt_on);
    CHECK(wt_on[3] == doctest::Approx(2.0 * median::kVetoWeightScale));  // veto fired
    const Scalar von_terr = (v_on.value.t - t_true).norm();
    CHECK(von_terr <= coupled_terr + 0.01);                              // parity bound
    CHECK(geo_dist(v_on.value.R, R_true) <= coupled_rerr + 0.01);
    CHECK(von_terr < 0.1);                               // and within the coupled bound

    // Veto OFF: the fault's translation weight is untouched, and its (per-channel-only
    // down-weighted) contribution drags the translation median MEASURABLY further off
    // truth than the coupled rejection — this is exactly what the veto-removal mutation
    // reproduces (the parity CHECK above then fails).
    Scalar wr_off[4], wt_off[4];
    const median::SplitResult v_off =
        median::solve_split(xs, w, w, 4, p, /*veto=*/false, wr_off, wt_off);
    CHECK(wt_off[3] == doctest::Approx(2.0));                            // untouched
    const Scalar voff_terr = (v_off.value.t - t_true).norm();
    MESSAGE("veto parity: coupled terr=" << coupled_terr << " veto-on terr=" << von_terr
            << " veto-off terr=" << voff_terr);
    CHECK(voff_terr > coupled_terr + 0.02);              // the drag the veto removes
}

TEST_CASE("split median veto: does NOT fire on clean quality differences (graceful "
          "per-channel weighting preserved)") {
    // Sources that disagree only within their normal scatter: no channel distance exceeds
    // 3x the leave-one-out spread, so the veto must leave every weight untouched and the
    // result must equal the veto-OFF solve bit-exactly.
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.00, 0.00, 0.30), Vec3(1.00, 0.00, 0.0));
    xs[1] = make_se3(Vec3(0.02, 0.00, 0.32), Vec3(1.06, 0.04, 0.0));
    xs[2] = make_se3(Vec3(0.00, 0.02, 0.28), Vec3(0.94, -0.05, 0.0));
    xs[3] = make_se3(Vec3(-0.02, 0.01, 0.31), Vec3(1.03, 0.02, 0.0));
    const Scalar w[4] = {1.0, 0.8, 1.2, 1.0};

    Scalar wr_on[4], wt_on[4];
    const median::SplitResult v_on =
        median::solve_split(xs, w, w, 4, params(), /*veto=*/true, wr_on, wt_on);
    const median::SplitResult v_off =
        median::solve_split(xs, w, w, 4, params(), /*veto=*/false);
    for (int i = 0; i < 4; ++i) {
        CHECK(wr_on[i] == doctest::Approx(w[i]));
        CHECK(wt_on[i] == doctest::Approx(w[i]));
    }
    CHECK((v_on.value.t.array() == v_off.value.t.array()).all());     // bit-identical
    CHECK((v_on.value.R.array() == v_off.value.R.array()).all());
    CHECK(v_on.spread_rot == v_off.spread_rot);
    CHECK(v_on.spread_trans == v_off.spread_trans);
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------
TEST_CASE("split median is deterministic: identical inputs give identical output") {
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.0, 0.1, 0.3), Vec3(1.0, 0.2, -0.3));
    xs[1] = make_se3(Vec3(0.1, 0.0, 0.2), Vec3(1.1, 0.1, -0.2));
    xs[2] = make_se3(Vec3(-0.1, 0.1, 0.4), Vec3(0.9, 0.3, -0.4));
    xs[3] = make_se3(Vec3(0.0, -0.1, 0.3), Vec3(1.0, 0.2, -0.3));
    const Scalar wr[4] = {1.0, 0.7, 1.3, 0.9};
    const Scalar wt[4] = {0.9, 1.1, 0.6, 1.4};

    const median::SplitResult a = median::solve_split(xs, wr, wt, 4, params());
    const median::SplitResult b = median::solve_split(xs, wr, wt, 4, params());
    CHECK((a.value.R.array() == b.value.R.array()).all());   // bit-identical
    CHECK((a.value.t.array() == b.value.t.array()).all());
    CHECK(a.spread_rot == b.spread_rot);
    CHECK(a.spread_trans == b.spread_trans);
}
