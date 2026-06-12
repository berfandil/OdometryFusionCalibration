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

// ===========================================================================================
// Per-channel split median (Slice 19, D3 amendment). Two independent Weiszfeld solvers —
// SO(3) geodesic for rotation, R^3 Euclidean for translation — each carrying the SAME D3
// safeguards as the coupled solver above: OFF-VERTEX init (the vertex-init pinning bug must
// not be reintroduced per channel), the Vardi-Zhang d<=eps guard, eps-regularized 1/d, the
// max_iters hard cap, and the n<=2 closed forms. See median.hpp for the full contract.
// ===========================================================================================
namespace {

// Per-channel effective-weight view: clamp negatives to 0; if the channel's whole set is
// <= 0, fall back to uniform (the per-channel mirror of solve()'s w_of contract).
struct ChannelWeights {
    const Scalar* w = nullptr;
    bool uniform    = false;
    void bind(const Scalar* ws, int n) {
        w = ws;
        Scalar sum = Scalar(0);
        for (int i = 0; i < n; ++i) {
            if (ws != nullptr && ws[i] > Scalar(0)) sum += ws[i];
        }
        uniform = (sum <= Scalar(0));
    }
    Scalar of(int i) const {
        if (uniform || w == nullptr) return Scalar(1);
        return (w[i] > Scalar(0)) ? w[i] : Scalar(0);
    }
};

// Rotation-channel solve: weighted geometric median on SO(3) under the geodesic distance.
struct RotChannel {
    Mat3   R = Mat3::Identity();
    Scalar spread    = Scalar(0);
    int    iters     = 0;
    bool   converged = true;
};

RotChannel solve_rot_channel(const SE3* deltas, const ChannelWeights& cw, int n,
                             const Params& p) {
    RotChannel out;
    if (n <= 0) return out;                          // identity (degenerate; caller avoids)
    if (n == 1) { out.R = deltas[0].R; return out; } // passthrough, spread 0

    const Scalar eps = std::max(p.eps, Scalar(0));

    if (n == 2) {
        // Weighted geodesic interpolation (the SO(3) part of se3::interpolate): the median
        // of two points is degenerate along the connecting geodesic — no rejection.
        const Scalar w0 = cw.of(0), w1 = cw.of(1);
        const Scalar denom = w0 + w1;
        const Scalar u = (denom > Scalar(0)) ? (w1 / denom) : Scalar(0.5);
        out.R = deltas[0].R * so3::exp(u * so3::log(deltas[0].R.transpose() * deltas[1].R));
        out.iters = 1;
    } else {
        // n >= 3: Weiszfeld IRLS on SO(3). OFF-VERTEX INIT (D3): one weighted Karcher
        // (tangent-mean) step about the highest-weight R — an interior point, so every
        // d_i > 0 at iter 0 and the 1/d reweight engages (a vertex start would re-create
        // the d=0 self-weight pin per channel).
        int start = 0;
        {
            Scalar best = cw.of(0);
            for (int i = 1; i < n; ++i) {
                if (cw.of(i) > best) { best = cw.of(i); start = i; }
            }
        }
        Mat3 Rm;
        {
            Scalar wsum   = Scalar(0);
            Vec3   r_mean = Vec3::Zero();
            const Mat3 RbT = deltas[start].R.transpose();
            for (int i = 0; i < n; ++i) {
                const Scalar w = cw.of(i);
                wsum   += w;
                r_mean += w * so3::log(RbT * deltas[i].R);
            }
            if (wsum > Scalar(0)) r_mean /= wsum;
            Rm = deltas[start].R * so3::exp(r_mean);
        }
        const int max_iters = std::max(1, p.max_iters);
        out.converged = false;
        for (int it = 0; it < max_iters; ++it) {
            out.iters = it + 1;
            Scalar usum  = Scalar(0);
            Vec3   r_acc = Vec3::Zero();
            const Mat3 RmT = Rm.transpose();
            for (int i = 0; i < n; ++i) {
                const Vec3 r = so3::log(RmT * deltas[i].R);
                const Scalar d = r.norm();
                if (d <= eps) continue;              // Vardi-Zhang coincident-vertex guard
                const Scalar u = cw.of(i) / d;
                usum  += u;
                r_acc += u * r;
            }
            if (usum <= Scalar(0)) { out.converged = true; break; }
            const Vec3 r_step = r_acc / usum;
            Rm = Rm * so3::exp(r_step);
            if (r_step.norm() < p.tol) { out.converged = true; break; }
        }
        out.R = Rm;
    }

    // Weighted RMS geodesic spread (rad) at the channel median, with the SAME effective
    // weights the solve used (clamped / uniform-resolved).
    Scalar acc_w = Scalar(0), acc_d2 = Scalar(0);
    const Mat3 RmT = out.R.transpose();
    for (int i = 0; i < n; ++i) {
        const Scalar w = cw.of(i);
        const Scalar d = so3::log(RmT * deltas[i].R).norm();
        acc_w  += w;
        acc_d2 += w * d * d;
    }
    out.spread = (acc_w > Scalar(0)) ? std::sqrt(acc_d2 / acc_w) : Scalar(0);
    return out;
}

// Translation-channel solve: weighted geometric median on R^3 (Euclidean).
struct TransChannel {
    Vec3   t = Vec3::Zero();
    Scalar spread    = Scalar(0);
    int    iters     = 0;
    bool   converged = true;
};

TransChannel solve_trans_channel(const SE3* deltas, const ChannelWeights& cw, int n,
                                 const Params& p) {
    TransChannel out;
    if (n <= 0) return out;
    if (n == 1) { out.t = deltas[0].t; return out; }

    const Scalar eps = std::max(p.eps, Scalar(0));

    if (n == 2) {
        const Scalar w0 = cw.of(0), w1 = cw.of(1);
        const Scalar denom = w0 + w1;
        const Scalar u = (denom > Scalar(0)) ? (w1 / denom) : Scalar(0.5);
        out.t = (Scalar(1) - u) * deltas[0].t + u * deltas[1].t;
        out.iters = 1;
    } else {
        // n >= 3: Weiszfeld IRLS on R^3. OFF-VERTEX INIT (D3): the weighted arithmetic
        // mean — interior, every d_i > 0 at iter 0 (same rationale as the rotation channel).
        Vec3 m = Vec3::Zero();
        {
            Scalar wsum = Scalar(0);
            for (int i = 0; i < n; ++i) {
                const Scalar w = cw.of(i);
                wsum += w;
                m    += w * deltas[i].t;
            }
            if (wsum > Scalar(0)) m /= wsum;
        }
        const int max_iters = std::max(1, p.max_iters);
        out.converged = false;
        for (int it = 0; it < max_iters; ++it) {
            out.iters = it + 1;
            Scalar usum  = Scalar(0);
            Vec3   t_acc = Vec3::Zero();
            for (int i = 0; i < n; ++i) {
                const Scalar d = (deltas[i].t - m).norm();
                if (d <= eps) continue;              // Vardi-Zhang coincident-vertex guard
                const Scalar u = cw.of(i) / d;
                usum  += u;
                t_acc += u * deltas[i].t;
            }
            if (usum <= Scalar(0)) { out.converged = true; break; }
            const Vec3 t_new = t_acc / usum;
            const Scalar step = (t_new - m).norm();
            m = t_new;
            if (step < p.tol) { out.converged = true; break; }
        }
        out.t = m;
    }

    Scalar acc_w = Scalar(0), acc_d2 = Scalar(0);
    for (int i = 0; i < n; ++i) {
        const Scalar w = cw.of(i);
        const Scalar d = (deltas[i].t - out.t).norm();
        acc_w  += w;
        acc_d2 += w * d * d;
    }
    out.spread = (acc_w > Scalar(0)) ? std::sqrt(acc_d2 / acc_w) : Scalar(0);
    return out;
}

} // namespace

SplitResult solve_split(const SE3* deltas, const Scalar* w_rot, const Scalar* w_trans,
                        int n, const Params& p, bool veto,
                        Scalar* w_rot_final, Scalar* w_trans_final) {
    SplitResult out;
    if (n <= 0 || deltas == nullptr) {
        out.value = SE3{};   // identity, spreads 0, converged (degenerate; caller avoids)
        return out;
    }

    ChannelWeights cw_rot, cw_trans;
    cw_rot.bind(w_rot, n);
    cw_trans.bind(w_trans, n);

    // ---- base solves (one per channel) ----------------------------------------------------
    RotChannel   rc = solve_rot_channel(deltas, cw_rot, n, p);
    TransChannel tc = solve_trans_channel(deltas, cw_trans, n, p);
    out.iters_rot   = rc.iters;
    out.iters_trans = tc.iters;

    // Final effective weights (what the FINAL channel solves consumed): start at the
    // clamped / uniform-resolved base weights; the veto below may scale entries. Fixed
    // scratch (strict core); n is bounded by the estimator's source cap == kMaxSplitInputs.
    Scalar we_rot[kMaxSplitInputs];
    Scalar we_trans[kMaxSplitInputs];
    const int nc = (n < kMaxSplitInputs) ? n : kMaxSplitInputs;
    for (int i = 0; i < nc; ++i) {
        we_rot[i]   = cw_rot.of(i);
        we_trans[i] = cw_trans.of(i);
    }

    // ---- cross-channel outlier veto (Slice 19; >= 3 inputs, fixed constants) --------------
    // A source gross in EITHER channel (d_i > kVetoNormDist x the LEAVE-ONE-OUT channel
    // spread — see median.hpp for why leave-one-out) gets its OTHER-channel weight scaled by
    // kVetoWeightScale and that other channel is re-solved ONCE. Flags are computed from the
    // BASE solves for both channels before either re-solve (no cascade); WCET is bounded at
    // 2 base solves + 2 re-solves, all iteration-capped. No heap.
    if (veto && n >= 3 && n <= kMaxSplitInputs) {
        // Channel distances at the base medians + the weighted sums for leave-one-out.
        Scalar d_rot[kMaxSplitInputs], d_trans[kMaxSplitInputs];
        Scalar Sw_r = Scalar(0), Sd2_r = Scalar(0);
        Scalar Sw_t = Scalar(0), Sd2_t = Scalar(0);
        const Mat3 RmT = rc.R.transpose();
        for (int i = 0; i < n; ++i) {
            d_rot[i]   = so3::log(RmT * deltas[i].R).norm();
            d_trans[i] = (deltas[i].t - tc.t).norm();
            Sw_r  += we_rot[i];
            Sd2_r += we_rot[i] * d_rot[i] * d_rot[i];
            Sw_t  += we_trans[i];
            Sd2_t += we_trans[i] * d_trans[i] * d_trans[i];
        }
        bool any_rot_outlier = false, any_trans_outlier = false;
        bool rot_flag[kMaxSplitInputs], trans_flag[kMaxSplitInputs];
        for (int i = 0; i < n; ++i) {
            // Leave-one-out spread of the OTHER inputs in each channel (0 if degenerate),
            // FLOORED at the absolute per-channel constant (review MAJOR-1): coincident
            // other inputs (loo -> 0) must not turn the flag into a d_i > ~0 hair trigger
            // that fires on honest noise — see kVetoSpreadFloor* in median.hpp.
            const Scalar wr  = we_rot[i];
            const Scalar den_r = Sw_r - wr;
            const Scalar num_r = Sd2_r - wr * d_rot[i] * d_rot[i];
            const Scalar loo_r = (den_r > Scalar(0) && num_r > Scalar(0))
                                     ? std::sqrt(num_r / den_r) : Scalar(0);
            const Scalar wt  = we_trans[i];
            const Scalar den_t = Sw_t - wt;
            const Scalar num_t = Sd2_t - wt * d_trans[i] * d_trans[i];
            const Scalar loo_t = (den_t > Scalar(0) && num_t > Scalar(0))
                                     ? std::sqrt(num_t / den_t) : Scalar(0);
            const Scalar loo_r_eff = (loo_r > kVetoSpreadFloorRot) ? loo_r
                                                                   : kVetoSpreadFloorRot;
            const Scalar loo_t_eff = (loo_t > kVetoSpreadFloorTrans) ? loo_t
                                                                     : kVetoSpreadFloorTrans;
            rot_flag[i]   = d_rot[i]   > kVetoNormDist * loo_r_eff;
            trans_flag[i] = d_trans[i] > kVetoNormDist * loo_t_eff;
            any_rot_outlier   = any_rot_outlier || rot_flag[i];
            any_trans_outlier = any_trans_outlier || trans_flag[i];
        }
        // A rotation-channel outlier suppresses that source's TRANSLATION weight (and vice
        // versa); the affected channel is re-solved once with the scaled weights.
        //
        // ALL-SOURCES-FLAGGED edge (review MINOR): the veto SCALES by kVetoWeightScale —
        // it never zeroes, so every source keeps voting whatever the flag pattern; and
        // Weiszfeld is scale-invariant, so if every source in one channel were flagged the
        // re-solve would equal the base solve (wasted-but-capped WCET, identical result).
        // That pattern is in fact UNREACHABLE under the LOO normalization: d_i > 3*loo_i
        // for all i implies sum_i w_i^2 d_i^2 > (sum w) * (sum w d^2), which contradicts
        // w_i <= sum w — at most n-1 sources can be flagged per channel. Both facts (never
        // zero, scale-invariant re-solve) are pinned by tests, not just this argument.
        if (any_rot_outlier) {
            for (int i = 0; i < n; ++i) {
                if (rot_flag[i]) we_trans[i] *= kVetoWeightScale;
            }
            ChannelWeights cw2; cw2.bind(we_trans, n);
            tc = solve_trans_channel(deltas, cw2, n, p);
            out.iters_trans += tc.iters;
        }
        if (any_trans_outlier) {
            for (int i = 0; i < n; ++i) {
                if (trans_flag[i]) we_rot[i] *= kVetoWeightScale;
            }
            ChannelWeights cw2; cw2.bind(we_rot, n);
            rc = solve_rot_channel(deltas, cw2, n, p);
            out.iters_rot += rc.iters;
        }
    }

    if (w_rot_final != nullptr)   { for (int i = 0; i < nc; ++i) w_rot_final[i]   = we_rot[i]; }
    if (w_trans_final != nullptr) { for (int i = 0; i < nc; ++i) w_trans_final[i] = we_trans[i]; }

    out.value.R       = rc.R;
    out.value.t       = tc.t;
    out.spread_rot    = rc.spread;        // recomputed by the re-solve when the veto fired
    out.spread_trans  = tc.spread;
    out.converged_rot   = rc.converged;
    out.converged_trans = tc.converged;
    return out;
}

} // namespace median
} // namespace ofc
