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
// INTERIOR INIT (D3 fix). The iteration starts at the weighted ARITHMETIC MEAN (an interior
// point with every d_i>0), NOT at a data vertex. A vertex start made iter-0 see d_start=0 ->
// self-weight w/eps (~1e9), which PINNED the estimate on that vertex and "converged" in one
// step -> fusion degenerated to the highest-weight INPUT verbatim and a HIGH-WEIGHT outlier
// (its own d=0 self-weight is immune to the 1/d reweight) was never rejected. The Vardi-Zhang
// coincident-vertex guard (skip a d<=eps self-term) is the safety net should an iterate land
// exactly on a vertex.
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
//
// `ws` are the raw caller weights; the spread MUST use the SAME effective weights
// the solver used: each entry clamped to >= 0, and (when the raw set is all <= 0)
// uniform 1.0 (signalled by `uniform`). Using the raw weights would (a) let a
// negative weight subtract from the accumulators and yield NaN/garbage, and
// (b) report 0 for the all-<=0 uniform-fallback case even though the median is the
// well-defined uniform solution — understating Q. We clamp/substitute here so the
// effective weights match w_of() in solve() exactly.
Scalar weighted_rms_spread(const SE3& m, const SE3* xs, const Scalar* ws, int n,
                           Scalar lambda, bool uniform) {
    Scalar acc_w  = Scalar(0);
    Scalar acc_d2 = Scalar(0);
    for (int i = 0; i < n; ++i) {
        const Scalar w = uniform ? Scalar(1)
                                 : ((ws[i] > Scalar(0)) ? ws[i] : Scalar(0));
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
        out.spread = weighted_rms_spread(out.value, deltas, weights, n, p.lambda,
                                         uniform);
        out.iters = 1;
        out.converged = true;
        return out;
    }

    // --- n >= 3: Weiszfeld IRLS (robust interior geometric median) ----------
    // Highest-weight source: used only as the rotation tangent BASEPOINT below.
    int start = 0;
    {
        Scalar best = w_of(0);
        for (int i = 1; i < n; ++i) {
            if (w_of(i) > best) { best = w_of(i); start = i; }
        }
    }

    const Scalar eps = std::max(p.eps, Scalar(0));
    const int max_iters = std::max(1, p.max_iters);

    // INITIALIZE OFF-VERTEX at the weighted mean (NOT at a data vertex). Vertex init made
    // iter-0 see d_start=0 -> self-weight w/eps (~1e9) -> pinned on that vertex and "converged"
    // in one iteration, returning the highest-weight INPUT verbatim (and never rejecting a
    // high-weight outlier). The weighted mean is interior, so every d_i>0 and the 1/d reweight
    // engages. translation = weighted arithmetic mean; rotation = one weighted tangent (Karcher)
    // step about the highest-weight source's R.
    SE3 m;
    {
        Scalar wsum   = Scalar(0);
        Vec3   t_mean = Vec3::Zero();
        const Mat3 RbT = deltas[start].R.transpose();
        Vec3   r_mean = Vec3::Zero();
        for (int i = 0; i < n; ++i) {
            const Scalar w = w_of(i);
            wsum   += w;
            t_mean += w * deltas[i].t;
            r_mean += w * so3::log(RbT * deltas[i].R);
        }
        if (wsum > Scalar(0)) { t_mean /= wsum; r_mean /= wsum; }
        m.t = t_mean;
        m.R = deltas[start].R * so3::exp(r_mean);
    }

    out.converged = false;
    for (int it = 0; it < max_iters; ++it) {
        out.iters = it + 1;
        Scalar usum = Scalar(0);
        Vec3   t_acc = Vec3::Zero();
        Vec3   r_acc = Vec3::Zero();
        const Mat3 RmT = m.R.transpose();
        for (int i = 0; i < n; ++i) {
            const Scalar d = se3::split_distance(m, deltas[i], p.lambda);
            // Vardi-Zhang coincident-vertex guard: exclude a d<=eps self-term (else its w/eps
            // self-weight re-pins). Off-vertex init means this rarely fires; it is the safety net.
            if (d <= eps) continue;
            const Scalar u = w_of(i) / d;
            usum  += u;
            t_acc += u * deltas[i].t;
            r_acc += u * so3::log(RmT * deltas[i].R);
        }
        if (usum <= Scalar(0)) { out.converged = true; break; }
        const Vec3 t_new  = t_acc / usum;
        const Vec3 r_step = r_acc / usum;
        const Vec3 t_step = t_new - m.t;
        const Scalar step_norm = std::sqrt(r_step.squaredNorm() + p.lambda * t_step.squaredNorm());
        m.t = t_new;
        m.R = m.R * so3::exp(r_step);
        if (step_norm < p.tol) { out.converged = true; break; }
    }

    out.value  = m;
    out.spread = weighted_rms_spread(m, deltas, weights, n, p.lambda, uniform);
    return out;
}

} // namespace median
} // namespace ofc
