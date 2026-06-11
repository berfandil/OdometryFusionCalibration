// test_multi_bias.cpp — Slice 18 FD PIN (acceptance §3 item 1) and its OUTCOME RECORD.
//
// Slice 18 (median-coupled multi-source bias, 11b Option B) mandates: pin the coupling
// Jacobian with a finite-difference test against the ACTUAL de-bias + frame-align + median
// pipeline BEFORE wiring F, and STOP if FD disagrees with the design's analytic block
//
//     J_spec_i = -dt * omega_i * Ad(X_i),
//     omega_i  = (w_i/max(d_i,eps)) / sum_j (w_j/max(d_j,eps))     (SLICE18 §2, spike 5a)
//
// OUTCOME (2026-06-11, this file): the scalar-omega block is FALSIFIED for the TRUE interior
// Weiszfeld median — FD relative error 0.38..0.98 across randomized states (28 blocks).
// Spike 5a's FD verification predates the D3 median-pinning fix: against the OLD pinning
// median the weights were degenerate (1,0,...,0) and the median was an input PASSTHROUGH,
// for which omega*I is exact — that regime no longer exists. The interior median's input
// sensitivity carries a RADIAL-PROJECTOR structure the scalar omega cannot express (a
// perturbation of an input ALONG its offset from the median does not move the median).
//
// The EXACT first-order sensitivity was derived here by the implicit function theorem at the
// converged Weiszfeld fixed point and is FD-VERIFIED to ~1% (worst 1.5%) over the same
// randomized states — using ONLY ingredients the estimator already has (per-source
// split_distance d_i to med.value + the fusion weights; no median API change):
//
//     fixed point:  F(m, x) = sum_i u_i xi_i = 0
//         xi_i = [t_i - t_m ; log(R_m^T R_i)]            ([trans; rot], tangent at the median)
//         u_i  = w_i / d_i ,  d_i^2 = xi_i^T W xi_i ,  W = diag(lambda*I3, I3)
//     dF/dx_i =  u_i P_i ,  P_i = I6 - xi_i (W xi_i)^T / d_i^2     (W-radial projector)
//     dF/dm   = -sum_j u_j P_j = -M
//  => the MATRIX-VALUED median influence   Omega_i = M^-1 u_i P_i   with  sum_i Omega_i = I
//     (the generalization of the scalar omega_i, sum omega_i = 1; omega_i*I is its isotropic
//     approximation), and the pose<->bias coupling block
//
//     J_exact_i = blkdiag(R_m^T, I) * Omega_i * blkdiag(R_{A_i}, I) * (-dt * Ad(X_i))
//
//     (the blkdiag factors are the right-tangent translation transports between the input
//     A_i, the ambient median equations, and the output median tangent; the -dt*Ad(X_i)
//     factor — the sign/frame convention — is UNCHANGED from the design).
//
// Regime pins (all FD-verified below):
//   * n == 1 (sole driver):    J = -dt * Ad(X)              -> Omega = I (Option A exactly).
//   * n == 2 (geodesic midpt): Omega_i = (w_i / sum w) * I   (the solver interpolates with
//     FIXED weights at n=2 — NOT the distance-based omega formula, which gives
//     w_i^2/(w_0^2+w_1^2) and mismatches FD for unequal weights).
//   * n >= 3 (Weiszfeld):      Omega_i = M^-1 u_i P_i; absent source -> u_i = 0 -> Omega_i = 0
//     (subsumes the frozen path); a vertex-dominant source -> Omega_i -> I.
//
// Per the slice brief, implementation STOPPED at this gate (no F wiring on a falsified
// Jacobian). This file is the permanent record: it asserts the verified exact block AND the
// falsification of the scalar block, so a design amendment inherits a ready-made acceptance
// item 1 (swap the candidate under test) and a regression of either finding is loud.
#include <doctest/doctest.h>

#include "ofc/core/eskf.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/median.hpp"
#include "ofc/core/types.hpp"

#include <cmath>
#include <random>
#include <vector>

using namespace ofc;

namespace {

// The pipeline under test: per-source de-bias (SOURCE-frame B_i' = B_i o exp(-b_i*dt)),
// frame-align A_i = X_i o B_i' o X_i^-1, weighted split-metric median.
SE3 debias_align_median(const std::vector<SE3>& B, const std::vector<SE3>& X,
                        const std::vector<Vec6>& bias, const std::vector<Scalar>& w,
                        Scalar dt, const median::Params& mp, std::vector<SE3>* A_out) {
    const int n = static_cast<int>(B.size());
    std::vector<SE3> A(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const SE3 Bp = se3::compose(B[static_cast<size_t>(i)],
                                    se3::exp(-bias[static_cast<size_t>(i)] * dt));
        A[static_cast<size_t>(i)] =
            se3::compose(se3::compose(X[static_cast<size_t>(i)], Bp),
                         se3::inverse(X[static_cast<size_t>(i)]));
    }
    if (A_out != nullptr) *A_out = A;
    return median::solve(A.data(), w.data(), n, mp).value;
}

// FD Jacobian d(right-tangent of the median)/d(bias_i) at epsilon e (central differences).
Mat6 fd_jacobian(const std::vector<SE3>& B, const std::vector<SE3>& X,
                 const std::vector<Vec6>& bias, const std::vector<Scalar>& w,
                 Scalar dt, const median::Params& mp, const SE3& med0, int i, Scalar e) {
    Mat6 J;
    for (int j = 0; j < 6; ++j) {
        std::vector<Vec6> bp = bias, bm = bias;
        bp[static_cast<size_t>(i)](j) += e;
        bm[static_cast<size_t>(i)](j) -= e;
        const SE3 gp = debias_align_median(B, X, bp, w, dt, mp, nullptr);
        const SE3 gm = debias_align_median(B, X, bm, w, dt, mp, nullptr);
        const Vec6 xp = se3::log(se3::compose(se3::inverse(med0), gp));
        const Vec6 xm = se3::log(se3::compose(se3::inverse(med0), gm));
        J.col(j) = (xp - xm) / (2.0 * e);
    }
    return J;
}

// Scalar Weiszfeld fixed-point weights omega_i (the SLICE18 §2 formula, under test).
void omega_weights(const SE3& med, const std::vector<SE3>& A, const std::vector<Scalar>& w,
                   Scalar lambda, Scalar eps, std::vector<Scalar>& omega) {
    const int n = static_cast<int>(A.size());
    omega.assign(static_cast<size_t>(n), 0.0);
    Scalar usum = 0.0;
    for (int i = 0; i < n; ++i) {
        const Scalar d = se3::split_distance(med, A[static_cast<size_t>(i)], lambda);
        const Scalar u = w[static_cast<size_t>(i)] / std::max(d, eps);
        omega[static_cast<size_t>(i)] = u;
        usum += u;
    }
    if (usum > 0.0) for (auto& o : omega) o /= usum;
}

// EXACT first-order median influence Omega_i = M^-1 u_i P_i (implicit function theorem at
// the converged fixed point; see the file header for the derivation).
void exact_influence(const SE3& med, const std::vector<SE3>& A, const std::vector<Scalar>& w,
                     Scalar lambda, std::vector<Mat6>& Omega) {
    const int n = static_cast<int>(A.size());
    Mat6 W = Mat6::Identity();
    W.topLeftCorner(3, 3) *= lambda;
    Mat6 M = Mat6::Zero();
    std::vector<Mat6> uP(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const SE3& Ai = A[static_cast<size_t>(i)];
        Vec6 xi;
        xi.head<3>() = Ai.t - med.t;
        xi.tail<3>() = so3::log(med.R.transpose() * Ai.R);
        const Scalar d2 = std::max(xi.dot(W * xi), 1e-30);
        const Scalar u  = w[static_cast<size_t>(i)] / std::sqrt(d2);
        const Mat6 P = Mat6::Identity() - (xi * (W * xi).transpose()) / d2;
        uP[static_cast<size_t>(i)] = u * P;
        M += u * P;
    }
    Omega.assign(static_cast<size_t>(n), Mat6::Zero());
    for (int i = 0; i < n; ++i) {
        // M is W-self-adjoint, not symmetric — full-pivot solve for robustness.
        Omega[static_cast<size_t>(i)] = M.fullPivLu().solve(uP[static_cast<size_t>(i)]);
    }
}

// The verified coupling block: tangent transports * exact influence * the design's
// (unchanged) sign/frame factor -dt*Ad(X_i).
Mat6 exact_block(const SE3& med, const SE3& Ai, const Mat6& Omega_i, const SE3& Xi, Scalar dt) {
    Mat6 Tg = Mat6::Identity();  Tg.topLeftCorner(3, 3) = med.R.transpose();
    Mat6 Ta = Mat6::Identity();  Ta.topLeftCorner(3, 3) = Ai.R;
    return Tg * Omega_i * Ta * (-dt * se3::adjoint(Xi));
}

std::mt19937_64& rng() { static std::mt19937_64 e(0xC0FFEEu); return e; }
Scalar urand(Scalar lo, Scalar hi) {
    std::uniform_real_distribution<double> d(lo, hi);
    return static_cast<Scalar>(d(rng()));
}
Vec3 vrand(Scalar s) { return Vec3(urand(-s, s), urand(-s, s), urand(-s, s)); }
Vec6 wrand(Scalar st, Scalar sr) {
    Vec6 v; v.head<3>() = vrand(st); v.tail<3>() = vrand(sr); return v;
}

// One randomized rig: shared base motion + per-source spread (sources roughly agree, the
// median sits interior — the realistic fusion-window shape).
void random_rig(int n, Scalar dt, std::vector<SE3>& B, std::vector<SE3>& X,
                std::vector<Vec6>& bias, std::vector<Scalar>& w) {
    B.resize(static_cast<size_t>(n));
    X.resize(static_cast<size_t>(n));
    bias.assign(static_cast<size_t>(n), Vec6::Zero());
    w.resize(static_cast<size_t>(n));
    const Vec6 base_twist = wrand(2.0, 0.6);
    for (int i = 0; i < n; ++i) {
        X[static_cast<size_t>(i)].R = so3::exp(vrand(0.6));
        X[static_cast<size_t>(i)].t = vrand(0.4);
        const Vec6 tw = base_twist + wrand(0.3, 0.1);
        const SE3 Abase = se3::exp(tw * dt);
        // Source-frame report B = X^-1 A X (the estimator's frame-align inverts this).
        B[static_cast<size_t>(i)] = se3::compose(
            se3::compose(se3::inverse(X[static_cast<size_t>(i)]), Abase),
            X[static_cast<size_t>(i)]);
        bias[static_cast<size_t>(i)] = wrand(0.2, 0.08);
        w[static_cast<size_t>(i)] = urand(0.5, 2.0);
    }
}

median::Params tight_params(Scalar lambda) {
    median::Params mp;
    mp.max_iters = 50000;     // run to the numerical fixed point (FD needs a converged solver)
    mp.tol       = 1e-14;
    mp.eps       = 1e-12;
    mp.lambda    = lambda;
    return mp;
}

// FD epsilon: signal ~dt*eps = 1e-5 sits far above the solver's slow-tail position noise
// (~1e-8 on hard configurations) and far below the nonlinearity scale (min d_i ~ 1e-2).
constexpr Scalar kFdEps = 1e-4;

} // namespace

// ===========================================================================
// n >= 3 (Weiszfeld): the exact influence matches FD; the scalar omega does NOT.
// ===========================================================================
TEST_CASE("slice18 FD pin: n>=3 median coupling — exact implicit-function block verified, "
          "scalar -dt*omega_i*Ad(X_i) falsified") {
    const Scalar dt = 0.1, lambda = 1.0;
    const median::Params mp = tight_params(lambda);

    Scalar worst_exact = 0.0, worst_spec_best_case = 1e30, worst_fd_consistency = 0.0;
    Scalar worst_sum_rule = 0.0;
    int blocks = 0;

    for (int trial = 0; trial < 8; ++trial) {
        const int n = 3 + (trial % 2);
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(n, dt, B, X, bias, w);

        std::vector<SE3> A;
        const SE3 med0 = debias_align_median(B, X, bias, w, dt, mp, &A);
        std::vector<Scalar> omega;
        omega_weights(med0, A, w, lambda, mp.eps, omega);
        std::vector<Mat6> Omg;
        exact_influence(med0, A, w, lambda, Omg);

        // Influence decomposition sanity: sum_i Omega_i = I (the matrix sum rule that
        // generalizes sum omega_i = 1).
        Mat6 sumO = Mat6::Zero();
        for (int i = 0; i < n; ++i) sumO += Omg[static_cast<size_t>(i)];
        worst_sum_rule = std::max(worst_sum_rule, (sumO - Mat6::Identity()).norm());

        for (int i = 0; i < n; ++i) {
            const Mat6 J_fd  = fd_jacobian(B, X, bias, w, dt, mp, med0, i, kFdEps);
            const Mat6 J_fd2 = fd_jacobian(B, X, bias, w, dt, mp, med0, i, kFdEps * 10.0);
            const Scalar scale = std::max(J_fd.norm(), 1e-12);
            const Scalar fd_consistency = (J_fd - J_fd2).norm() / scale;   // FD trustworthy?

            const Mat6 J_spec = -dt * omega[static_cast<size_t>(i)] *
                                se3::adjoint(X[static_cast<size_t>(i)]);
            const Mat6 J_exact = exact_block(med0, A[static_cast<size_t>(i)],
                                             Omg[static_cast<size_t>(i)],
                                             X[static_cast<size_t>(i)], dt);

            const Scalar e_spec  = (J_fd - J_spec).norm() / scale;
            const Scalar e_exact = (J_fd - J_exact).norm() / scale;
            worst_exact          = std::max(worst_exact, e_exact);
            worst_spec_best_case = std::min(worst_spec_best_case, e_spec);
            worst_fd_consistency = std::max(worst_fd_consistency, fd_consistency);
            ++blocks;

            MESSAGE("trial " << trial << " src " << i
                    << " omega=" << omega[static_cast<size_t>(i)]
                    << "  rel-err spec=" << e_spec << "  rel-err exact=" << e_exact
                    << "  fd-consistency=" << fd_consistency);
        }
    }

    MESSAGE("FD summary over " << blocks << " blocks: worst exact-err=" << worst_exact
            << "  BEST-case spec-err=" << worst_spec_best_case
            << "  worst fd-consistency=" << worst_fd_consistency
            << "  worst |sum Omega - I|=" << worst_sum_rule);

    // FD itself is trustworthy at this epsilon (two epsilons agree to ~1%).
    REQUIRE(worst_fd_consistency < 0.02);
    // The influence decomposition is exact: sum_i Omega_i = I.
    CHECK(worst_sum_rule < 1e-9);
    // VERIFIED: the exact implicit-function coupling block matches FD (residual = the
    // small-angle tangent-transport terms, ~1% at window-sized rotations).
    CHECK(worst_exact < 0.05);
    // FALSIFIED (the slice-18 STOP finding): the design's scalar block does not match FD
    // anywhere — even its BEST block across all randomized states is >30% off. If this ever
    // starts passing tightly, the median solver's character changed (e.g. pinning returned)
    // — investigate before trusting either result.
    CHECK(worst_spec_best_case > 0.3);
}

// ===========================================================================
// n == 1 (sole driver / Option-A regime): passthrough, J = -dt*Ad(X).
// ===========================================================================
TEST_CASE("slice18 FD pin: n==1 sole driver — J = -dt*Ad(X) (Option A generalizes exactly)") {
    const Scalar dt = 0.1, lambda = 1.0;
    const median::Params mp = tight_params(lambda);

    Scalar worst = 0.0;
    for (int trial = 0; trial < 4; ++trial) {
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(1, dt, B, X, bias, w);
        // Keep the bias small so the O(|b|*dt) right-Jacobian-of-exp refinement (present in
        // Option A's -dt*I too; accepted) stays below the pin tolerance.
        bias[0] = wrand(0.02, 0.02);

        const SE3 med0 = debias_align_median(B, X, bias, w, dt, mp, nullptr);
        const Mat6 J_fd   = fd_jacobian(B, X, bias, w, dt, mp, med0, 0, kFdEps);
        const Mat6 J_spec = -dt * se3::adjoint(X[0]);   // omega = 1 (passthrough)
        const Scalar err = (J_fd - J_spec).norm() / std::max(J_fd.norm(), 1e-12);
        worst = std::max(worst, err);
        MESSAGE("n==1 trial " << trial << " rel-err=" << err);
    }
    // Tight: validates the sign + frame + Ad-order conventions of the whole harness.
    CHECK(worst < 0.02);
}

// ===========================================================================
// n == 2 (geodesic midpoint): influence = the FIXED interpolation weight w_i/sum(w),
// NOT the distance-based omega formula.
// ===========================================================================
TEST_CASE("slice18 FD pin: n==2 midpoint — influence is w_i/sum(w), not the d-based omega") {
    const Scalar dt = 0.1, lambda = 1.0;
    const median::Params mp = tight_params(lambda);

    Scalar worst_interp = 0.0;
    Scalar best_dbased  = 1e30;
    for (int trial = 0; trial < 4; ++trial) {
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(2, dt, B, X, bias, w);
        w[0] = 2.0; w[1] = 0.5;        // deliberately unequal: the two formulas separate
                                        // (interp 0.8/0.2 vs d-based ~0.94/0.06)

        std::vector<SE3> A;
        const SE3 med0 = debias_align_median(B, X, bias, w, dt, mp, &A);
        std::vector<Scalar> om_d;
        omega_weights(med0, A, w, lambda, mp.eps, om_d);

        for (int i = 0; i < 2; ++i) {
            const Mat6 J_fd = fd_jacobian(B, X, bias, w, dt, mp, med0, i, kFdEps);
            const Scalar scale = std::max(J_fd.norm(), 1e-12);
            const Scalar om_interp = w[static_cast<size_t>(i)] / (w[0] + w[1]);
            const Mat6 J_interp = -dt * om_interp * se3::adjoint(X[static_cast<size_t>(i)]);
            const Mat6 J_dbased = -dt * om_d[static_cast<size_t>(i)] *
                                  se3::adjoint(X[static_cast<size_t>(i)]);
            const Scalar e_int = (J_fd - J_interp).norm() / scale;
            const Scalar e_db  = (J_fd - J_dbased).norm() / scale;
            worst_interp = std::max(worst_interp, e_int);
            best_dbased  = std::min(best_dbased, e_db);
            MESSAGE("n==2 trial " << trial << " src " << i << " interp-err=" << e_int
                    << "  d-based-err=" << e_db);
        }
    }
    // The solver interpolates with FIXED weights at n==2, so the interpolation weight is the
    // correct influence (small-angle transport residual only)...
    CHECK(worst_interp < 0.06);
    // ...and the d-based omega formula systematically mismatches for unequal weights.
    CHECK(best_dbased > 0.10);
}

// ===========================================================================
// PRODUCTION == PIN (Slice 18 implementation guard). The estimator-side production blocks
// (Eskf::median_influence + Eskf::bias_coupling) must equal THIS FILE's FD-verified reference
// construction (exact_influence + exact_block) on randomized states — guards drift between
// the permanent FD pin above and the production implementation the filter actually runs.
// ===========================================================================
TEST_CASE("slice18 production: Eskf::median_influence/bias_coupling == the FD-pinned reference "
          "construction on randomized states") {
    const Scalar dt = 0.1, lambda = 1.0;
    const median::Params mp = tight_params(lambda);

    Scalar worst_omega = 0.0, worst_J = 0.0;
    for (int trial = 0; trial < 8; ++trial) {
        const int n = 3 + (trial % 2);
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(n, dt, B, X, bias, w);

        std::vector<SE3> A;
        const SE3 med0 = debias_align_median(B, X, bias, w, dt, mp, &A);

        // Reference (the FD-verified construction in this file).
        std::vector<Mat6> Omg_ref;
        exact_influence(med0, A, w, lambda, Omg_ref);

        // Production (what the estimator computes inside step()).
        std::vector<Mat6> Omg_prod(static_cast<size_t>(n));
        Eskf::median_influence(med0, A.data(), w.data(), n, lambda, mp.eps,
                               Omg_prod.data());

        for (int i = 0; i < n; ++i) {
            const Scalar eo = (Omg_prod[static_cast<size_t>(i)] -
                               Omg_ref[static_cast<size_t>(i)]).norm();
            const Mat6 J_ref  = exact_block(med0, A[static_cast<size_t>(i)],
                                            Omg_ref[static_cast<size_t>(i)],
                                            X[static_cast<size_t>(i)], dt);
            const Mat6 J_prod = Eskf::bias_coupling(med0, A[static_cast<size_t>(i)],
                                                    Omg_prod[static_cast<size_t>(i)],
                                                    X[static_cast<size_t>(i)], dt);
            const Scalar ej = (J_prod - J_ref).norm();
            worst_omega = std::max(worst_omega, eo);
            worst_J     = std::max(worst_J, ej);
        }
    }
    MESSAGE("production-vs-pin: worst |Omega diff|=" << worst_omega
            << "  worst |J diff|=" << worst_J);
    // Same math, same conventions -> agreement to numerical noise (NOT a percent tolerance:
    // any re-derivation drift, sign flip, or transport reorder fails loudly).
    CHECK(worst_omega < 1e-9);
    CHECK(worst_J < 1e-9);
}

// ===========================================================================
// Influence edge cases (acceptance §3 item 5, unit level): absent source (w = 0 -> Omega = 0),
// sole participant (Omega = I -> J = -dt*Ad(X), Option A exact), n == 2 fixed interpolation,
// VZ coincident vertex, weight-dominant vertex (Omega -> I, no blow-up).
// ===========================================================================
TEST_CASE("slice18 production: influence edge cases — absent, sole, n==2, coincident vertex, "
          "weight-dominant vertex") {
    const Scalar dt = 0.1, lambda = 1.0, eps = 1e-12;
    const median::Params mp = tight_params(lambda);

    // --- sole participant (n == 1): Omega = I, J = -dt*Ad(X) (Option A generalizes exactly).
    {
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(1, dt, B, X, bias, w);
        std::vector<SE3> A;
        const SE3 med0 = debias_align_median(B, X, bias, w, dt, mp, &A);
        Mat6 Om;
        Eskf::median_influence(med0, A.data(), w.data(), 1, lambda, eps, &Om);
        CHECK((Om - Mat6::Identity()).norm() < 1e-12);
        // med == A[0] at n == 1, so the two transports cancel: J = -dt*Ad(X).
        const Mat6 J = Eskf::bias_coupling(med0, A[0], Om, X[0], dt);
        CHECK((J - (-dt * se3::adjoint(X[0]))).norm() < 1e-9);
    }

    // --- n == 2: fixed interpolation weights w_i/sum(w) (NOT the d-based form).
    {
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(2, dt, B, X, bias, w);
        w[0] = 2.0; w[1] = 0.5;
        std::vector<SE3> A;
        const SE3 med0 = debias_align_median(B, X, bias, w, dt, mp, &A);
        Mat6 Om[2];
        Eskf::median_influence(med0, A.data(), w.data(), 2, lambda, eps, Om);
        CHECK((Om[0] - (2.0 / 2.5) * Mat6::Identity()).norm() < 1e-12);
        CHECK((Om[1] - (0.5 / 2.5) * Mat6::Identity()).norm() < 1e-12);
    }

    // --- absent source: a zero-weight participant gets Omega = 0 exactly (u_i = 0); the
    //     remaining influences still sum to I.
    {
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(4, dt, B, X, bias, w);
        w[2] = 0.0;                          // "absent" (the estimator passes it no weight)
        std::vector<SE3> A;
        const SE3 med0 = debias_align_median(B, X, bias, w, dt, mp, &A);
        Mat6 Om[4];
        Eskf::median_influence(med0, A.data(), w.data(), 4, lambda, eps, Om);
        CHECK(Om[2].norm() < 1e-12);
        Mat6 sum = Mat6::Zero();
        for (int i = 0; i < 4; ++i) sum += Om[i];
        CHECK((sum - Mat6::Identity()).norm() < 1e-9);
    }

    // --- VZ coincident vertex: two inputs EXACTLY coincident, jointly weight-dominant -> the
    //     median sits ON them (the solver's Vardi-Zhang regime; we pass the exact vertex pose,
    //     the limit point, so the test does not depend on the solver's vertex asymptotics).
    //     The influence must not blow up: the coincident set splits I by weight, the
    //     off-vertex source gets 0.
    {
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(3, dt, B, X, bias, w);
        std::vector<SE3> A;
        (void)debias_align_median(B, X, bias, w, dt, mp, &A);
        A[1] = A[0];                          // exact coincidence
        w[0] = 1.0; w[1] = 1.0; w[2] = 0.5;   // w0 + w1 > w2 -> median == the coincident pair
        Mat6 Om[3];
        Eskf::median_influence(A[0], A.data(), w.data(), 3, lambda, eps, Om);
        bool finite = true;
        for (int i = 0; i < 3; ++i) finite = finite && Om[i].allFinite();
        CHECK(finite);
        CHECK((Om[0] - 0.5 * Mat6::Identity()).norm() < 1e-9);
        CHECK((Om[1] - 0.5 * Mat6::Identity()).norm() < 1e-9);
        CHECK(Om[2].norm() < 1e-12);
    }

    // --- weight-dominant vertex: w_0 > sum(others) -> the median IS vertex 0 (d_0 = 0). The
    //     documented limit: Omega_0 -> I, others -> 0, no 1/d blow-up. ALSO pin the NEAR-vertex
    //     numerics (d_0 ~ 1e-7, just above the eps guard): finite + bounded (graceful
    //     degradation toward I, the honest-caveat regime of the design).
    {
        std::vector<SE3> B, X;
        std::vector<Vec6> bias;
        std::vector<Scalar> w;
        random_rig(3, dt, B, X, bias, w);
        std::vector<SE3> A;
        (void)debias_align_median(B, X, bias, w, dt, mp, &A);
        w[0] = 10.0; w[1] = 0.3; w[2] = 0.3;
        Mat6 Om[3];
        // Exactly ON the vertex (d_0 = 0 <= eps): the limit branch.
        Eskf::median_influence(A[0], A.data(), w.data(), 3, lambda, eps, Om);
        bool finite = true;
        for (int i = 0; i < 3; ++i) finite = finite && Om[i].allFinite();
        CHECK(finite);
        CHECK((Om[0] - Mat6::Identity()).norm() < 1e-9);
        CHECK(Om[1].norm() < 1e-12);
        CHECK(Om[2].norm() < 1e-12);
        // JUST OFF the vertex (d_0 ~ 1e-7 > eps): the smooth formula must stay finite and
        // bounded (no 1/d blow-up into the covariance), with Omega_0 near I.
        SE3 med_near = A[0];
        med_near.t += Vec3(1e-7, 0, 0);
        Eskf::median_influence(med_near, A.data(), w.data(), 3, lambda, eps, Om);
        finite = true;
        for (int i = 0; i < 3; ++i) finite = finite && Om[i].allFinite();
        CHECK(finite);
        // Bounded — no 1/d blow-up (|I| = sqrt(6) ~ 2.45; the near-vertex limit is I minus a
        // bounded W-radial rank-1 deficiency, the design's documented non-smoothness caveat).
        CHECK(Om[0].norm() < 10.0);
        Mat6 sum = Mat6::Zero();
        for (int i = 0; i < 3; ++i) sum += Om[i];
        CHECK((sum - Mat6::Identity()).norm() < 1e-6);    // the sum rule still holds
    }
}
