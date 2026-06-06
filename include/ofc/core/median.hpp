// ofc/core/median.hpp — weighted geometric median of SE(3) deltas (Slice 2, D3).
//
// Robust consensus motion via Weiszfeld IRLS under the SPLIT metric
//   d(a,b)^2 = ||log(Ra^T Rb)||^2 + lambda * ||ta - tb||^2
// (no bi-invariant metric on SE(3) -> rotation and translation are medianed
// separately, tied by the weight lambda; see DESIGN §4, DECISIONS D3).
//
// STRICT CORE: caller-supplied fixed arrays (no heap); the iteration count is
// hard-capped (bounded WCET); double math, no exceptions.
#ifndef OFC_CORE_MEDIAN_HPP
#define OFC_CORE_MEDIAN_HPP

#include "ofc/core/types.hpp"

namespace ofc {
namespace median {

// Solver knobs (subset of MedianConfig in CONFIG §3, passed explicitly so the
// median module has no dependency on the full Config tree).
struct Params {
    int    max_iters = 10;     // hard iteration cap (bounded WCET)
    Scalar tol       = 1e-6;   // convergence: tangent step norm below this stops
    Scalar eps       = 1e-9;   // epsilon-regularized 1/d weight (vertex guard)
    Scalar lambda    = 1.0;    // rotation-vs-translation split-metric weight
};

// Output of the median solve.
//   value     — the consensus SE(3) delta.
//   spread    — weighted RMS split-distance of the inputs to the median; drives
//               the ESKF's adaptive Q (tight agreement -> small, disagreement -> large).
//   iters     — Weiszfeld iterations actually run (<= Params::max_iters).
//   converged — true iff the iteration met `tol` before the cap (always true for
//               the n<=2 closed forms).
struct Result {
    SE3    value;
    Scalar spread    = Scalar(0);
    int    iters     = 0;
    bool   converged = true;
};

// Weighted geometric median of `n` SE(3) deltas with per-input weights.
//   n == 0 -> identity, spread 0 (degenerate; caller should avoid).
//   n == 1 -> passthrough (the median of one point is that point).
//   n == 2 -> WEIGHTED GEODESIC MIDPOINT. The geometric median of two points is
//             not unique along the connecting geodesic, so it degenerates: we
//             return the weight-interpolated midpoint (NO outlier rejection
//             possible with two inputs — see DESIGN §4, needs >= 3).
//   n >= 3 -> Weiszfeld IRLS; rejects a single gross outlier via the 1/d reweight.
// Weights must be non-negative; a zero/empty weight set falls back to uniform.
// `deltas` and `weights` point to caller-owned arrays of length >= n.
Result solve(const SE3* deltas, const Scalar* weights, int n, const Params& p);

} // namespace median
} // namespace ofc
#endif // OFC_CORE_MEDIAN_HPP
