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

// ---- Per-channel split median (Slice 19, D3 amendment) -----------------------------------
//
// solve_split() runs TWO independent Weiszfeld medians — rotation on SO(3) under the
// geodesic distance ||log(R_m^T R_i)|| with weights `w_rot`, translation on R^3 under the
// Euclidean distance with weights `w_trans` — so fusion can express PER-CHANNEL source
// quality (a heading-grade source can dominate rotation without distorting translation).
// The coupled solve() above stays bit-untouched (the default path); the split solver is
// mathematically different even at equal weights (the coupled IRLS couples the channels
// through the shared 1/d reweight).
//
// Both channel solvers carry the SAME D3 safeguards as the fixed coupled solver:
//   * OFF-VERTEX INIT — rotation: one weighted Karcher (tangent-mean) step about the
//     highest-weight R; translation: the weighted arithmetic mean. A vertex start made
//     iter-0 see d=0 -> w/eps self-weight -> pinned on that vertex (the D3 pinning bug);
//     the interior init keeps every d_i > 0 so the 1/d reweight engages.
//   * VARDI-ZHANG coincident-vertex guard — a d <= eps self-term is skipped per iteration.
//   * eps-regularized 1/d, max_iters hard cap (bounded WCET), n <= 2 closed forms
//     (n == 2 = weighted geodesic/linear interpolation; no rejection possible).
// Per-channel weights are clamped to >= 0 with a PER-CHANNEL uniform fallback when a
// channel's whole set is <= 0 (mirrors solve()'s w_of contract).
//
// Per-channel spreads: weighted RMS geodesic distance (rad) and weighted RMS Euclidean
// distance (m) of the inputs to the channel median — NO lambda unit-mixing (resolves the
// D21 scalar-spread wart for the split path). They drive the per-channel adaptive Q.
//
// CROSS-CHANNEL OUTLIER VETO (`veto`, the Config::split_veto knob — default ON for the
// split path). Hard sensor faults usually corrupt BOTH channels; the veto recovers the
// coupled solver's whole-source rejection for gross outliers while leaving graceful
// per-channel weighting for quality differences. After the two base solves, a source whose
// channel distance d_i exceeds kVetoNormDist x the channel spread in EITHER channel gets
// its weight in the OTHER channel scaled by kVetoWeightScale, and that other channel is
// re-solved ONCE (bounded WCET: at most one extra IRLS per channel per step; at most
// 2 base solves + 2 veto re-solves total).
//   NORMALIZATION NOTE (implementation choice, documented): the veto threshold normalizes
//   d_i by the LEAVE-ONE-OUT spread (the weighted RMS over the OTHER inputs). The naive
//   full-set RMS includes the outlier's own d_i^2, which bounds d_i/spread by
//   sqrt(sum_w/w_i) (< 3 for any rig under ~10 equal-weight sources) — i.e. a single gross
//   outlier could NEVER trip the 3.0 threshold. Leave-one-out is the minimal normalization
//   that makes the guard live at small n while keeping the fixed constants.
//   ABSOLUTE FLOOR (Slice-19 review MAJOR-1): the LOO spread is floored at a per-channel
//   absolute constant (kVetoSpreadFloorRot / kVetoSpreadFloorTrans) BEFORE the ratio.
//   Without it, (near-)coincident OTHER inputs drive loo -> 0 and the flag degenerates to
//   d_i > ~0 — ANY honest noise then flags the source EVERY step. That is the project's
//   own target-rig shape (KAIST: all sources share wheel translation, so the translation
//   channel is near-coincident by construction — the FOG heading source's ROTATION weight
//   would be permanently vetoed by mm-level translation deviations). With the floor, a
//   source must deviate by more than kVetoNormDist x the floor in absolute terms even when
//   the other inputs agree perfectly — only GROSS faults trip it, which is the veto's
//   contract (graceful per-channel weighting handles everything below).
// The veto needs >= 3 inputs (rejection is undefined below that, as in the coupled solver)
// and at most kMaxSplitInputs (the fixed scratch capacity, == the estimator's source cap).
//
// `w_rot_final` / `w_trans_final` (optional, length >= n, n <= kMaxSplitInputs): receive
// the EFFECTIVE per-channel weights the FINAL channel solves consumed — clamped/uniform-
// resolved and veto-scaled. The estimator feeds these to Eskf::median_influence_split so
// the Slice-18 coupling blocks describe the median ACTUALLY solved (veto included).
constexpr int kMaxSplitInputs = 32;

// Fixed veto policy constants (Slice 19 — deliberately NOT config knobs; only the
// split_veto bool is). See the veto contract above.
constexpr Scalar kVetoNormDist    = 3.0;   // flag: d_i > 3 x leave-one-out channel spread
constexpr Scalar kVetoWeightScale = 0.1;   // flagged source's OTHER-channel weight scale
// Absolute per-channel floors on the leave-one-out spread (review MAJOR-1) — in the
// channel's NATIVE distance units on a WINDOWED delta: rad (geodesic) / m (Euclidean).
// Sized for the project's window regime (0.05-0.2 s windows at vehicle dynamics): honest
// per-window sensor noise sits BELOW the floor (wheel/visual odometry per-window scatter
// ~mm-cm and ~mrad; the sim noise floors are 0.005 m / 0.005 rad), while a hard fault is
// orders of magnitude ABOVE the kVetoNormDist x floor trip point (0.06 m / 0.03 rad per
// window ~ 0.6 m/s / 0.3 rad/s at a 0.1 s window). Consequence (deliberate): when the
// other inputs are coincident, a fault must exceed the ABSOLUTE trip point to be vetoed —
// moderate degradation is left to the graceful per-channel weighting, exactly the policy
// split that motivated the veto. Window-length CAVEAT: the floors are window-agnostic
// constants; rigs with much longer windows (>~1 s) accumulate more honest per-window
// noise and would need these revisited (report, do not tune per rig — they are fixed
// policy like kVetoNormDist).
constexpr Scalar kVetoSpreadFloorRot   = 0.01;   // rad, per windowed delta
constexpr Scalar kVetoSpreadFloorTrans = 0.02;   // m,   per windowed delta

// Output of the split solve. `value` = { R = rotation-channel median, t = translation-
// channel median }; the per-channel spreads/iters/convergence mirror Result's fields.
struct SplitResult {
    SE3    value;
    Scalar spread_rot      = Scalar(0);   // weighted RMS geodesic distance (rad)
    Scalar spread_trans    = Scalar(0);   // weighted RMS Euclidean distance (m)
    int    iters_rot       = 0;           // total rotation-channel iterations (base + veto re-solve)
    int    iters_trans     = 0;           // total translation-channel iterations
    bool   converged_rot   = true;
    bool   converged_trans = true;
};

SplitResult solve_split(const SE3* deltas, const Scalar* w_rot, const Scalar* w_trans,
                        int n, const Params& p, bool veto = true,
                        Scalar* w_rot_final = nullptr, Scalar* w_trans_final = nullptr);

} // namespace median
} // namespace ofc
#endif // OFC_CORE_MEDIAN_HPP
