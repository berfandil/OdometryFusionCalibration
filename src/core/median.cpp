// ofc/core/median.cpp — weighted split-metric geometric median (Slice 2, D3).
//
// Weiszfeld IRLS. At the current estimate m, for each input x_i:
//   d_i = se3::split_distance(m, x_i, lambda)         (split SO(3)/R^3 metric)
//   u_i = w_i / max(d_i, eps)                          (epsilon-regularized 1/d)
// then update the two parts CONSISTENTLY at the same reweight:
//   translation: t_m <- sum(u_i t_i) / sum(u_i)        (weighted average in R^3)
//   rotation:    R_m <- R_m * exp( sum(u_i log(R_m^T R_i)) / sum(u_i) )
//                (tangent average about the current R_m — a Weiszfeld step on SO(3))
// Iterate to max_iters or until the tangent step norm < tol. The 1/d reweight is
// what makes this the geometric MEDIAN (robust) rather than the Karcher MEAN: a
// gross outlier sits at large d_i, so its u_i collapses and it stops dragging the
// estimate (DECISIONS D3).
//
// STRICT CORE: no heap (all accumulation is in fixed-size Eigen locals); the loop
// is hard-capped by max_iters; double math, no exceptions.
#include "ofc/core/median.hpp"

#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {
namespace median {

namespace {

// Weighted RMS split-distance of the inputs to a pose m: the spread that drives
// the adaptive Q. sqrt( sum(w_i d_i^2) / sum(w_i) ).
Scalar weighted_rms_spread(const SE3& m, const SE3* xs, const Scalar* ws, int n,
                           Scalar lambda) {
    Scalar acc_w  = Scalar(0);
    Scalar acc_d2 = Scalar(0);
    for (int i = 0; i < n; ++i) {
        const Scalar w = ws[i];
        const Scalar d = se3::split_distance(m, xs[i], lambda);
        acc_w  += w;
        acc_d2 += w * d * d;
    }
    if (acc_w <= Scalar(0)) return Scalar(0);
    return std::sqrt(acc_d2 / acc_w);
}

} // namespace

Result solve(const SE3* deltas, const Scalar* weights, int n, const Params& p) {
    Result out;

    if (n <= 0) {
        out.value = SE3{};      // identity
        out.spread = Scalar(0);
        out.converged = true;
        return out;
    }

    // Normalize weights to non-negative; if the set is degenerate (all <= 0),
    // fall back to uniform so the median is still well defined.
    // (No heap: we read weights[] directly and guard each access.)
    auto wsafe = [&](int i) -> Scalar {
        const Scalar w = weights[i];
        return (w > Scalar(0)) ? w : Scalar(0);
    };
    Scalar wsum_all = Scalar(0);
    for (int i = 0; i < n; ++i) wsum_all += wsafe(i);
    const bool uniform = (wsum_all <= Scalar(0));
    auto w_of = [&](int i) -> Scalar { return uniform ? Scalar(1) : wsafe(i); };

    // --- n == 1: passthrough -------------------------------------------------
    if (n == 1) {
        out.value = deltas[0];
        out.spread = Scalar(0);
        out.iters = 0;
        out.converged = true;
        return out;
    }

    // --- n == 2: weighted geodesic midpoint (median degenerate; no rejection) -
    if (n == 2) {
        const Scalar w0 = w_of(0);
        const Scalar w1 = w_of(1);
        const Scalar denom = w0 + w1;
        const Scalar u = (denom > Scalar(0)) ? (w1 / denom) : Scalar(0.5);
        out.value = se3::interpolate(deltas[0], deltas[1], u);
        out.spread = weighted_rms_spread(out.value, deltas, weights, n, p.lambda);
        out.iters = 1;
        out.converged = true;
        return out;
    }

    // --- n >= 3: Weiszfeld IRLS ---------------------------------------------
    // Initialize at the highest-weight input (a finite, data-anchored start that
    // keeps the rotation tangent average near its basepoint).
    int start = 0;
    {
        Scalar best = w_of(0);
        for (int i = 1; i < n; ++i) {
            if (w_of(i) > best) { best = w_of(i); start = i; }
        }
    }
    SE3 m = deltas[start];

    const Scalar eps = std::max(p.eps, Scalar(0));
    const int max_iters = std::max(1, p.max_iters);

    out.converged = false;
    for (int it = 0; it < max_iters; ++it) {
        out.iters = it + 1;

        Scalar usum = Scalar(0);
        Vec3   t_acc = Vec3::Zero();      // sum(u_i t_i)
        Vec3   r_acc = Vec3::Zero();      // sum(u_i log(R_m^T R_i))
        const Mat3 RmT = m.R.transpose();

        for (int i = 0; i < n; ++i) {
            const Scalar d = se3::split_distance(m, deltas[i], p.lambda);
            const Scalar u = w_of(i) / std::max(d, eps);
            usum  += u;
            t_acc += u * deltas[i].t;
            r_acc += u * so3::log(RmT * deltas[i].R);
        }

        if (usum <= Scalar(0)) { out.converged = true; break; }

        const Vec3 t_new   = t_acc / usum;
        const Vec3 r_step  = r_acc / usum;          // tangent update about R_m

        // Tangent step norm (split metric, same lambda) measures convergence.
        const Vec3 t_step = t_new - m.t;
        const Scalar step_norm =
            std::sqrt(r_step.squaredNorm() + p.lambda * t_step.squaredNorm());

        m.t = t_new;
        m.R = m.R * so3::exp(r_step);

        if (step_norm < p.tol) { out.converged = true; break; }
    }

    out.value  = m;
    out.spread = weighted_rms_spread(m, deltas, weights, n, p.lambda);
    return out;
}

} // namespace median
} // namespace ofc
