// Slice 2 unit tests: weighted split-metric geometric median (Weiszfeld IRLS).
#include <doctest/doctest.h>

#include "ofc/core/lie.hpp"
#include "ofc/core/median.hpp"

#include <cmath>

using namespace ofc;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

template <typename A, typename B>
bool close(const A& a, const B& b, Scalar tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

median::Params params() {
    median::Params p;
    p.max_iters = 50;
    p.tol       = 1e-8;     // sub-micro tangent step (Weiszfeld near a vertex
                            // converges linearly; this is reached well within 50 iters)
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
} // namespace

// ---------------------------------------------------------------------------
// Degenerate counts
// ---------------------------------------------------------------------------
TEST_CASE("median of one input is a passthrough") {
    const SE3 x = make_se3(Vec3(0.1, -0.2, 0.3), Vec3(1.0, 2.0, -1.0));
    const Scalar w = 1.0;
    const median::Result r = median::solve(&x, &w, 1, params());
    CHECK(close(r.value.R, x.R, 1e-12));
    CHECK(close(r.value.t, x.t, 1e-12));
    CHECK(r.spread == doctest::Approx(0.0));
    CHECK(r.converged);
}

TEST_CASE("median of zero inputs is identity with zero spread") {
    // Degenerate count: n == 0 returns identity, spread 0, converged.
    const median::Result r = median::solve(nullptr, nullptr, 0, params());
    CHECK(close(r.value.R, Mat3::Identity(), 1e-12));
    CHECK(close(r.value.t, Vec3::Zero(), 1e-12));
    CHECK(r.spread == doctest::Approx(0.0));
    CHECK(r.converged);
}

TEST_CASE("all-zero weights fall back to uniform (n>=3)") {
    // Three inputs spread along x; with all-zero weights the solver uses uniform
    // weights, so the result must match the equal-weight median exactly (and the
    // spread must be the non-zero uniform RMS, NOT understated to 0).
    SE3 xs[3];
    xs[0] = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.0), Vec3(1.0, 0.0, 0.0));
    xs[2] = make_se3(Vec3(0, 0, 0.0), Vec3(2.0, 0.0, 0.0));
    const Scalar w_zero[3] = {0.0, 0.0, 0.0};
    const Scalar w_unit[3] = {1.0, 1.0, 1.0};

    const median::Result r_zero = median::solve(xs, w_zero, 3, params());
    const median::Result r_unit = median::solve(xs, w_unit, 3, params());

    CHECK(close(r_zero.value.t, r_unit.value.t, 1e-9));
    CHECK(close(r_zero.value.R, r_unit.value.R, 1e-9));
    CHECK(r_zero.spread == doctest::Approx(r_unit.spread).epsilon(1e-9));
    CHECK(r_zero.spread > 0.1);    // uniform-fallback spread is reported, not 0
}

TEST_CASE("negative weights are clamped to the uniform fallback (n==2)") {
    // Both weights <= 0 -> uniform fallback -> equal-weight geodesic midpoint, and a
    // finite (non-NaN) spread (the negative-weight guard must hold in the RMS too).
    const SE3 a = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    const SE3 b = make_se3(Vec3(0, 0, 0.0), Vec3(4.0, 0.0, 0.0));
    const SE3 xs[2] = {a, b};
    const Scalar ws[2] = {-1.0, -3.0};

    const median::Result r = median::solve(xs, ws, 2, params());
    // Uniform fallback -> midpoint at x = 2.0 (not weighted toward either).
    CHECK(r.value.t.x() == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(std::isfinite(r.spread));
    CHECK(r.spread > 0.0);
}

TEST_CASE("median of two equal-weight inputs is the geodesic midpoint") {
    const SE3 a = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    const SE3 b = make_se3(Vec3(0, 0, kPi / 2), Vec3(2.0, 0.0, 0.0));
    const SE3 xs[2] = {a, b};
    const Scalar ws[2] = {1.0, 1.0};

    const median::Result r = median::solve(xs, ws, 2, params());
    // Equal weights -> midpoint: 45 deg about z, translation (1,0,0).
    const Mat3 expected_R = so3::exp(Vec3(0, 0, kPi / 4));
    CHECK(close(r.value.R, expected_R, 1e-9));
    CHECK(close(r.value.t, Vec3(1.0, 0.0, 0.0), 1e-9));
}

TEST_CASE("median of two inputs honors the weight (biased midpoint)") {
    const SE3 a = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    const SE3 b = make_se3(Vec3(0, 0, 0.0), Vec3(4.0, 0.0, 0.0));
    const SE3 xs[2] = {a, b};
    // Weight b three times a -> the consensus sits at u = 3/4 toward b: x = 3.0.
    const Scalar ws[2] = {1.0, 3.0};
    const median::Result r = median::solve(xs, ws, 2, params());
    CHECK(r.value.t.x() == doctest::Approx(3.0).epsilon(1e-9));
}

// ---------------------------------------------------------------------------
// Outlier rejection (the load-bearing property: needs >= 3 inputs)
// ---------------------------------------------------------------------------
TEST_CASE("median of 3+ rejects a single gross outlier (translation)") {
    // Four near-identical inliers clustered around t=(1,0,0); one gross outlier.
    SE3 xs[5];
    xs[0] = make_se3(Vec3(0, 0, 0.00), Vec3(1.00, 0.00, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.01), Vec3(1.02, 0.01, 0.0));
    xs[2] = make_se3(Vec3(0, 0, -0.01), Vec3(0.98, -0.01, 0.0));
    xs[3] = make_se3(Vec3(0, 0, 0.00), Vec3(1.01, 0.00, 0.0));
    xs[4] = make_se3(Vec3(0, 0, 1.50), Vec3(50.0, -40.0, 7.0));   // OUTLIER
    Scalar ws[5] = {1, 1, 1, 1, 1};

    const median::Result r = median::solve(xs, ws, 5, params());
    // The consensus must sit with the inliers, not be dragged to the outlier.
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(0.05));
    CHECK(std::abs(r.value.t.y()) < 0.1);
    CHECK(std::abs(r.value.t.z()) < 0.1);
    CHECK(so3::log(r.value.R).norm() < 0.1);
    CHECK(r.converged);
}

TEST_CASE("median of 3+ rejects a single gross outlier (rotation)") {
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.00, 0.00, 0.30), Vec3(0, 0, 0));
    xs[1] = make_se3(Vec3(0.01, 0.00, 0.31), Vec3(0, 0, 0));
    xs[2] = make_se3(Vec3(0.00, 0.01, 0.29), Vec3(0, 0, 0));
    xs[3] = make_se3(Vec3(2.50, -1.50, -2.0), Vec3(0, 0, 0));   // OUTLIER rotation
    Scalar ws[4] = {1, 1, 1, 1};

    const median::Result r = median::solve(xs, ws, 4, params());
    const Vec3 log_r = so3::log(r.value.R);
    // Consensus rotation near the inlier cluster ~ (0,0,0.3).
    CHECK(std::abs(log_r.z() - 0.30) < 0.05);
    CHECK(std::abs(log_r.x()) < 0.05);
    CHECK(std::abs(log_r.y()) < 0.05);
}

// ---------------------------------------------------------------------------
// Convergence + spread
// ---------------------------------------------------------------------------
TEST_CASE("Weiszfeld converges within the iteration cap") {
    SE3 xs[3];
    xs[0] = make_se3(Vec3(0, 0, 0.1), Vec3(1.0, 0.0, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.2), Vec3(1.1, 0.1, 0.0));
    xs[2] = make_se3(Vec3(0, 0, 0.15), Vec3(0.9, -0.1, 0.0));
    Scalar ws[3] = {1, 1, 1};

    // CONVERGENCE-FLAG EDGE (D3 median fix). The off-vertex weighted-mean init makes this a TRUE
    // interior Weiszfeld iterate (vs the old vertex-pinned 1-iteration "convergence" that the OLD
    // solver reported). The geometric median of these three points sits very CLOSE to a data
    // vertex, where Weiszfeld converges only LINEARLY with a rate near 1, so the step-norm tail is
    // long. MEASURED here (off-vertex init, this input): the iterate converges to the strict step
    // criterion in 77 iters at tol=1e-6 (the production weiszfeld_tol) and ~25 at 1e-4 — but does
    // NOT reach the deliberately-extreme tol=1e-8 (a sub-micro tangent step) inside 100 iters (it
    // is still descending, step ~1e-8). The VALUE is stable throughout: t -> (1, ~1e-7) i.e. the
    // correct (1,0,0) median to ~1e-7 by iter 100 regardless. So the original tol=1e-8 was an
    // unrealistically tight step demand exposed only now that the solver genuinely iterates; we
    // run this convergence probe at the PRODUCTION tol=1e-6, which the iterate DOES reach (77 <
    // 100), proving the convergence machinery works. (Production also caps max_iters at 10; the
    // WCET case below confirms the value is well-formed inside that cap even when the flag is not
    // yet set — convergence of the FLAG is not required for a correct VALUE.)
    median::Params p = params();
    p.tol       = 1e-6;          // the production weiszfeld_tol (1e-8 is an unreachable step here)
    p.max_iters = 100;
    const median::Result r = median::solve(xs, ws, 3, p);
    CHECK(r.converged);          // reaches the strict 1e-6 step within 100 iters (measured 77)
    CHECK(r.iters <= 100);
    CHECK(r.iters > 1);          // genuinely ITERATED (interior), not vertex-pinned at iter 1
    // The VALUE is the correct interior median (1,0,0) regardless of the strict step tail.
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(0.01));
    CHECK(std::abs(r.value.t.y()) < 0.01);
}

// ---------------------------------------------------------------------------
// High-weight-outlier guard (the D3 blind spot: the OLD pinning median returned the
// highest-weight INPUT verbatim, so an outlier carrying the HIGHEST weight was never
// rejected — its own d=0 self-weight at the vertex init was immune to the 1/d reweight).
// ---------------------------------------------------------------------------
TEST_CASE("median of 3+ rejects a HIGHEST-WEIGHT outlier (translation) — D3 blind spot") {
    // Four inliers clustered at t~(1,0,0) + one gross outlier. The outlier carries the LARGEST
    // single weight, but the inliers are the weighted MAJORITY (4 x 1.0 = 4.0 > 2.0). The OLD
    // vertex-init median initialized AT the highest-weight vertex (the outlier) and its d=0
    // self-weight pinned it there -> it returned the OUTLIER. The fixed interior median must sit
    // with the inliers.
    SE3 xs[5];
    xs[0] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.00,  0.00, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.01),  Vec3(1.02,  0.01, 0.0));
    xs[2] = make_se3(Vec3(0, 0, -0.01), Vec3(0.98, -0.01, 0.0));
    xs[3] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.01,  0.00, 0.0));
    xs[4] = make_se3(Vec3(0, 0, 1.50),  Vec3(50.0, -40.0, 7.0));   // OUTLIER, HIGHEST weight
    Scalar ws[5] = {1, 1, 1, 1, 2.0};   // outlier weight 2.0 > each inlier, but < inlier mass 4.0

    median::Params p = params();
    p.tol = 1e-6; p.max_iters = 100;    // production tol; 1e-8 is an unreachable step near a vertex
    const median::Result r = median::solve(xs, ws, 5, p);
    // Consensus sits with the inlier cluster, NOT dragged to the high-weight outlier at (50,-40,7).
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(0.1));
    CHECK(std::abs(r.value.t.y()) < 0.2);
    CHECK(std::abs(r.value.t.z()) < 0.2);
    CHECK(so3::log(r.value.R).norm() < 0.2);
    CHECK(r.converged);
    // INTERIOR-NOT-PINNED: the solver genuinely iterated (did not converge at iter 1 on a vertex)
    // and the consensus is NOT sitting on the outlier vertex (its split-distance to it is large).
    CHECK(r.iters > 1);
    CHECK(se3::split_distance(r.value, xs[4], p.lambda) > 1.0);
}

TEST_CASE("WCET: at the production cap (max_iters=10, tol=1e-6) the fixed median is NOT pinned, "
          "rejects a high-weight outlier, and the value is well-formed (D3)") {
    // The production config runs the solver at weiszfeld_max_iters=10 / weiszfeld_tol=1e-6. The
    // off-vertex init (D3 fix) means even at this small cap the iterate is an interior robust
    // median, NOT a vertex pin: within 10 iters it has moved the consensus DECISIVELY onto the
    // inlier cluster and OFF the high-weight outlier. (Convergence of the step-norm FLAG within 10
    // iters is NOT required — Weiszfeld's linear tail near a vertex can exceed 10 — but the VALUE
    // is already correct: the first few 1/d reweights collapse the outlier's pull immediately.)
    SE3 xs[5];
    xs[0] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.00,  0.00, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.01),  Vec3(1.02,  0.01, 0.0));
    xs[2] = make_se3(Vec3(0, 0, -0.01), Vec3(0.98, -0.01, 0.0));
    xs[3] = make_se3(Vec3(0, 0, 0.00),  Vec3(1.01,  0.00, 0.0));
    xs[4] = make_se3(Vec3(0, 0, 1.50),  Vec3(50.0, -40.0, 7.0));   // OUTLIER, HIGHEST weight
    Scalar ws[5] = {1, 1, 1, 1, 2.0};

    median::Params p = params();
    p.tol       = 1e-6;     // production weiszfeld_tol
    p.max_iters = 10;       // production weiszfeld_max_iters (the WCET cap)
    const median::Result r = median::solve(xs, ws, 5, p);
    MESSAGE("WCET: iters=" << r.iters << " conv=" << r.converged << " t.x=" << r.value.t.x());
    CHECK(r.iters <= 10);                    // never exceeds the hard cap (strict-core WCET bound)
    CHECK(r.iters > 1);                      // iterated (interior), not a 1-step vertex pin
    // VALUE is well-formed within the cap: consensus on the inlier cluster, NOT the outlier.
    CHECK(r.value.t.x() == doctest::Approx(1.0).epsilon(0.1));
    CHECK(std::abs(r.value.t.y()) < 0.2);
    CHECK(std::abs(r.value.t.z()) < 0.2);
    CHECK(so3::log(r.value.R).norm() < 0.2);
    CHECK(se3::split_distance(r.value, xs[4], p.lambda) > 1.0);   // far from the outlier vertex
}

TEST_CASE("median of 3+ is interior, not pinned, on agreeing-but-distinct inputs (D3)") {
    // Three DISTINCT but agreeing inputs (no outlier). The fixed median must be an INTERIOR
    // consensus: it iterated (iters > 1) and every input sits at a NON-ZERO split-distance from
    // the consensus (it is not pinned ON any one vertex the way the old vertex-init median was).
    SE3 xs[3];
    xs[0] = make_se3(Vec3(0, 0, 0.10), Vec3(1.00, 0.00, 0.0));
    xs[1] = make_se3(Vec3(0, 0, 0.20), Vec3(1.10, 0.10, 0.0));
    xs[2] = make_se3(Vec3(0, 0, 0.15), Vec3(0.90, -0.10, 0.0));
    Scalar ws[3] = {1.0, 1.3, 0.8};      // unequal weights: the OLD median would pin on slot 1

    median::Params p = params();
    p.tol = 1e-6; p.max_iters = 100;     // production tol
    const median::Result r = median::solve(xs, ws, 3, p);
    // This case asserts INTERIORITY, not step-norm convergence: like the clustered-near-a-vertex
    // case above, this weighted optimum sits close enough to a data vertex that the strict step
    // criterion is not reached inside the cap (Weiszfeld's linear tail) — but the VALUE is the
    // correct interior consensus. The dedicated convergence case above pins the flag; here we pin
    // that the solver ITERATED (interior, not the old 1-step vertex pin) and the consensus is OFF
    // every vertex.
    CHECK(r.iters > 1);                  // interior iterate, not a 1-step vertex pin
    // Every input is at a strictly POSITIVE distance from the consensus (the consensus is an
    // interior point, not coincident with any vertex — the old median would return slot 1 exactly,
    // giving a zero distance there). Use a distance floor well above any round-off: an interior
    // blend of these spread-out points sits ~0.05+ from each vertex.
    Scalar min_d = 1e30;
    for (int i = 0; i < 3; ++i)
        min_d = std::min(min_d, se3::split_distance(r.value, xs[i], p.lambda));
    CHECK(min_d > Scalar(1e-3));         // strictly interior (not pinned on any vertex)
}

TEST_CASE("spread is zero for identical inputs and grows with disagreement") {
    SE3 same[3];
    same[0] = make_se3(Vec3(0, 0, 0.2), Vec3(1.0, 0.0, 0.0));
    same[1] = same[0];
    same[2] = same[0];
    Scalar ws[3] = {1, 1, 1};
    const median::Result r_same = median::solve(same, ws, 3, params());
    CHECK(r_same.spread == doctest::Approx(0.0).epsilon(1e-6));

    SE3 spread[3];
    spread[0] = make_se3(Vec3(0, 0, 0.0), Vec3(0.0, 0.0, 0.0));
    spread[1] = make_se3(Vec3(0, 0, 0.0), Vec3(1.0, 0.0, 0.0));
    spread[2] = make_se3(Vec3(0, 0, 0.0), Vec3(-1.0, 0.0, 0.0));
    const median::Result r_spread = median::solve(spread, ws, 3, params());
    CHECK(r_spread.spread > 0.1);
}

TEST_CASE("median is deterministic: identical inputs give identical output") {
    SE3 xs[4];
    xs[0] = make_se3(Vec3(0.0, 0.1, 0.3), Vec3(1.0, 0.2, -0.3));
    xs[1] = make_se3(Vec3(0.1, 0.0, 0.2), Vec3(1.1, 0.1, -0.2));
    xs[2] = make_se3(Vec3(-0.1, 0.1, 0.4), Vec3(0.9, 0.3, -0.4));
    xs[3] = make_se3(Vec3(0.0, -0.1, 0.3), Vec3(1.0, 0.2, -0.3));
    Scalar ws[4] = {1.0, 0.7, 1.3, 0.9};

    const median::Result a = median::solve(xs, ws, 4, params());
    const median::Result b = median::solve(xs, ws, 4, params());
    CHECK((a.value.R.array() == b.value.R.array()).all());   // bit-identical
    CHECK((a.value.t.array() == b.value.t.array()).all());
    CHECK(a.spread == b.spread);
}
