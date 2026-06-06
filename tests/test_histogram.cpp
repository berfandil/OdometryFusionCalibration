// Slice 4 unit tests: Histogram1D — binning, weighted/split voting, decay &
// sliding-K aging, parabolic sub-bin peak, peak-concentration confidence, and the
// circular (wrap-seam) variant. Analytic expected values where possible.
#include <doctest/doctest.h>

#include "ofc/core/histogram.hpp"

#include <cmath>

using namespace ofc;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

// Default linear config over [lo, hi] with `n` bins.
HistogramConfig linear(int n, Scalar lo, Scalar hi, bool split = true,
                       bool sub = true) {
    HistogramConfig c;
    c.bins       = n;
    c.range_min  = lo;
    c.range_max  = hi;
    c.circular   = false;
    c.aging      = Aging::Decay;
    c.decay_gamma = 0.999;
    c.sliding_k  = 1000;
    c.vote_split = split;
    c.subbin     = sub;
    return c;
}

HistogramConfig circular(int n, Scalar lo, Scalar hi, bool split = true,
                         bool sub = true) {
    HistogramConfig c = linear(n, lo, hi, split, sub);
    c.circular = true;
    return c;
}
} // namespace

// ---------------------------------------------------------------------------
// configure() validation
// ---------------------------------------------------------------------------
TEST_CASE("configure validates bins, range, and aging knobs") {
    Histogram1D h;

    HistogramConfig too_few = linear(3, -1.0, 1.0);
    CHECK(h.configure(too_few) == Status::OutOfRange);

    HistogramConfig too_many = linear(Histogram1D::kMaxBins + 1, -1.0, 1.0);
    CHECK(h.configure(too_many) == Status::OutOfRange);

    HistogramConfig bad_range = linear(8, 1.0, 1.0);      // empty span
    CHECK(h.configure(bad_range) == Status::InvalidConfig);

    HistogramConfig bad_range2 = linear(8, 1.0, -1.0);    // inverted
    CHECK(h.configure(bad_range2) == Status::InvalidConfig);

    HistogramConfig bad_gamma = linear(8, -1.0, 1.0);
    bad_gamma.aging = Aging::Decay;
    bad_gamma.decay_gamma = 1.0;                          // must be < 1
    CHECK(h.configure(bad_gamma) == Status::OutOfRange);

    HistogramConfig bad_gamma2 = bad_gamma;
    bad_gamma2.decay_gamma = 0.0;                         // must be > 0
    CHECK(h.configure(bad_gamma2) == Status::OutOfRange);

    HistogramConfig bad_k = linear(8, -1.0, 1.0);
    bad_k.aging = Aging::SlidingK;
    bad_k.sliding_k = Histogram1D::kMaxSlidingK + 1;      // too large for the ring
    CHECK(h.configure(bad_k) == Status::OutOfRange);

    HistogramConfig bad_k2 = bad_k;
    bad_k2.sliding_k = 0;                                 // must be >= 1
    CHECK(h.configure(bad_k2) == Status::OutOfRange);

    // A valid config configures and starts empty.
    CHECK(h.configure(linear(64, -1.0, 1.0)) == Status::Ok);
    CHECK(h.configured());
    CHECK(h.bins() == 64);
    CHECK(h.empty());
    CHECK(h.peak_bin() == -1);
    CHECK(h.confidence() == doctest::Approx(0.0));
}

// ---------------------------------------------------------------------------
// Empty behavior: mode is the range midpoint; confidence 0; peak_bin -1.
// ---------------------------------------------------------------------------
TEST_CASE("empty histogram reports midpoint mode, zero confidence, no peak") {
    Histogram1D h;
    h.configure(linear(64, 2.0, 6.0));
    CHECK(h.empty());
    CHECK(h.total() == doctest::Approx(0.0));
    CHECK(h.peak_bin() == -1);
    CHECK(h.confidence() == doctest::Approx(0.0));
    CHECK(h.mode() == doctest::Approx(4.0));   // (2 + 6) / 2
}

// ---------------------------------------------------------------------------
// A single concentrated cluster -> mode ~ planted value, high confidence.
// ---------------------------------------------------------------------------
TEST_CASE("concentrated cluster recovers planted value with high confidence") {
    Histogram1D h;
    h.configure(linear(101, -5.0, 5.0));     // width ~0.099, bin center at value
    const Scalar planted = 1.3;

    // Tight Gaussian-ish cluster around the planted value.
    for (int k = -3; k <= 3; ++k) {
        h.add(planted + static_cast<Scalar>(k) * 0.01);
    }
    CHECK_FALSE(h.empty());
    CHECK(h.mode() == doctest::Approx(planted).epsilon(0.02));
    CHECK(h.confidence() > 0.9);              // almost all mass at the peak+neighbors
    CHECK(h.confidence() <= 1.0);
}

// ---------------------------------------------------------------------------
// Uniform votes -> low confidence; mode still well-defined and in range.
// ---------------------------------------------------------------------------
TEST_CASE("uniform votes give low confidence and an in-range mode") {
    Histogram1D h;
    // SlidingK with a window larger than the vote count -> exact accumulation
    // (no decay), so the analytic total/flatness holds.
    HistogramConfig c = linear(100, 0.0, 10.0, /*split=*/false);
    c.aging = Aging::SlidingK;
    c.sliding_k = 200;
    h.configure(c);

    // One vote into the center of each bin -> perfectly flat.
    for (int i = 0; i < 100; ++i) {
        h.add(0.05 + static_cast<Scalar>(i) * 0.1);
    }
    CHECK(h.total() == doctest::Approx(100.0));
    // Peak bin + 2 neighbors out of 100 flat bins ~ 3/100.
    CHECK(h.confidence() < 0.05);
    CHECK(h.confidence() >= 0.0);
    // Mode lands inside the range.
    CHECK(h.mode() >= 0.0);
    CHECK(h.mode() <= 10.0);
}

// ---------------------------------------------------------------------------
// vote_split deposits the correct linear fractions between two bins.
// ---------------------------------------------------------------------------
TEST_CASE("vote_split distributes weight linearly between the two nearest bins") {
    Histogram1D h;
    // 10 bins over [0, 10) -> width 1, centers at 0.5, 1.5, ... A value of 2.0
    // sits exactly between bin 1 (center 1.5) and bin 2 (center 2.5): 50/50.
    h.configure(linear(10, 0.0, 10.0, /*split=*/true, /*sub=*/false));
    h.add(2.0, 4.0);
    // Total mass conserved.
    CHECK(h.total() == doctest::Approx(4.0));
    // Equal split -> tie broken to the lower bin (peak_bin == 1), confidence:
    // peak(2) + both neighbors (bins 0 and 2). bin0=0, bin1=2, bin2=2 -> 4/4 = 1.
    CHECK(h.confidence() == doctest::Approx(1.0));

    // A value 3/4 of the way from center 1.5 to center 2.5 (= 2.25) -> 0.25 to
    // bin 1, 0.75 to bin 2.
    Histogram1D h2;
    h2.configure(linear(10, 0.0, 10.0, /*split=*/true, /*sub=*/false));
    h2.add(2.25, 1.0);
    CHECK(h2.total() == doctest::Approx(1.0));
    CHECK(h2.peak_bin() == 2);   // the larger fraction
}

// ---------------------------------------------------------------------------
// Parabolic sub-bin accuracy on a synthetic symmetric/triangular pattern.
//   Around the argmax i with neighbors h[i-1], h[i+1]:
//     offset = 0.5 * (h[i-1] - h[i+1]) / (h[i-1] - 2 h[i] + h[i+1])
//   Plant masses directly (split off, sub on) so the formula is exact.
// ---------------------------------------------------------------------------
TEST_CASE("parabolic sub-bin interpolation matches the analytic offset") {
    Histogram1D h;
    // 10 bins over [0, 10): width 1, centers k+0.5. Use split off so each add
    // lands wholly in one bin; SlidingK (large K) so masses accumulate exactly
    // (no decay between deposits would skew the planted triple).
    HistogramConfig c = linear(10, 0.0, 10.0, /*split=*/false, /*sub=*/true);
    c.aging = Aging::SlidingK;
    c.sliding_k = 100;
    h.configure(c);

    // Plant: h[3] = 2, h[4] = 8, h[5] = 6 (deposit at each bin center).
    for (int k = 0; k < 2; ++k) h.add(3.5);
    for (int k = 0; k < 8; ++k) h.add(4.5);
    for (int k = 0; k < 6; ++k) h.add(5.5);

    CHECK(h.peak_bin() == 4);
    // offset = 0.5*(2 - 6)/(2 - 16 + 6) = 0.5*(-4)/(-8) = 0.25
    const Scalar expected_offset = 0.5 * (2.0 - 6.0) / (2.0 - 2.0 * 8.0 + 6.0);
    const Scalar expected_mode = 4.5 + expected_offset * 1.0;   // center + offset*width
    CHECK(h.mode() == doctest::Approx(expected_mode));
    CHECK(expected_offset == doctest::Approx(0.25));
}

TEST_CASE("subbin=false returns the bin center") {
    Histogram1D h;
    h.configure(linear(10, 0.0, 10.0, /*split=*/false, /*sub=*/false));
    for (int k = 0; k < 5; ++k) h.add(4.5);
    for (int k = 0; k < 2; ++k) h.add(5.5);
    CHECK(h.peak_bin() == 4);
    CHECK(h.mode() == doctest::Approx(4.5));   // bin center, no parabola
}

// ---------------------------------------------------------------------------
// Decay aging: votes at A then many at B -> mode tracks B; the A peak washes out.
// ---------------------------------------------------------------------------
TEST_CASE("decay aging lets the mode migrate from A to B as A decays away") {
    HistogramConfig c = linear(101, -5.0, 5.0);
    c.aging = Aging::Decay;
    c.decay_gamma = 0.97;                      // brisk decay so A washes out fast
    Histogram1D h;
    h.configure(c);

    const Scalar A = -2.0, B = 3.0;
    for (int k = 0; k < 50; ++k) h.add(A);
    // Initially the mode is at A.
    CHECK(h.mode() == doctest::Approx(A).epsilon(0.05));

    for (int k = 0; k < 200; ++k) h.add(B);
    // The old A votes have decayed by gamma^200 (~negligible); mode is at B.
    CHECK(h.mode() == doctest::Approx(B).epsilon(0.05));
    CHECK(h.confidence() > 0.9);

    // The A bin's residual mass is tiny vs the B peak.
    CHECK(h.peak_bin() != -1);
}

// ---------------------------------------------------------------------------
// SlidingK aging: after > K votes the oldest no longer influence mode/total;
// total saturates at K (for unit weights).
// ---------------------------------------------------------------------------
TEST_CASE("sliding-K window saturates total and forgets the oldest votes") {
    HistogramConfig c = linear(101, -5.0, 5.0, /*split=*/false);
    c.aging = Aging::SlidingK;
    c.sliding_k = 20;
    Histogram1D h;
    h.configure(c);

    const Scalar A = -2.0, B = 3.0;
    // 20 votes at A: total fills the window.
    for (int k = 0; k < 20; ++k) h.add(A);
    CHECK(h.total() == doctest::Approx(20.0));
    CHECK(h.mode() == doctest::Approx(A).epsilon(0.05));

    // 5 more at A keeps total saturated (oldest A evicted as newest A enters).
    for (int k = 0; k < 5; ++k) h.add(A);
    CHECK(h.total() == doctest::Approx(20.0));

    // 20 votes at B fully replace the window -> not a single A vote remains.
    for (int k = 0; k < 20; ++k) h.add(B);
    CHECK(h.total() == doctest::Approx(20.0));
    CHECK(h.mode() == doctest::Approx(B).epsilon(0.05));
    CHECK(h.confidence() > 0.9);   // window is pure B now
}

TEST_CASE("sliding-K with weighted votes saturates at K entries (not K mass)") {
    HistogramConfig c = linear(50, 0.0, 5.0, /*split=*/false);
    c.aging = Aging::SlidingK;
    c.sliding_k = 3;
    Histogram1D h;
    h.configure(c);

    h.add(1.0, 2.0);
    h.add(1.0, 2.0);
    h.add(1.0, 2.0);
    CHECK(h.total() == doctest::Approx(6.0));   // 3 entries * weight 2
    h.add(1.0, 2.0);                            // evicts one, adds one
    CHECK(h.total() == doctest::Approx(6.0));   // still 3 live entries
}

// ---------------------------------------------------------------------------
// Circular: votes clustered near the wrap seam (~±pi) peak AT the seam, not 0,
// and neighbor indexing + sub-bin wrap correctly.
// ---------------------------------------------------------------------------
TEST_CASE("circular histogram peaks at the wrap seam, not the middle") {
    Histogram1D h;
    // Roll over [-pi, pi), circular. 72 bins -> 5 deg each.
    h.configure(circular(72, -kPi, kPi));

    // Cluster straddling the seam: a few just under +pi, a few just over -pi.
    h.add(kPi - 0.02);
    h.add(kPi - 0.04);
    h.add(-kPi + 0.02);
    h.add(-kPi + 0.04);
    h.add(kPi - 0.01);
    h.add(-kPi + 0.01);

    // The recovered mode should be near the seam (|mode| ~ pi), NOT near 0.
    const Scalar m = h.mode();
    const Scalar dist_to_seam = std::abs(std::abs(m) - kPi);
    CHECK(dist_to_seam < 0.1);
    CHECK(std::abs(m) > 3.0);          // definitely not near 0
    CHECK(h.confidence() > 0.9);       // tight cluster -> concentrated
}

TEST_CASE("circular folds out-of-range values into the range") {
    Histogram1D h;
    h.configure(circular(72, -kPi, kPi, /*split=*/false));
    // 3*pi wraps to pi -> folds to -pi (the seam); -3*pi likewise.
    h.add(3.0 * kPi);
    h.add(-3.0 * kPi);
    h.add(kPi - 0.001);
    CHECK_FALSE(h.empty());
    // All votes collapse near the seam; mode magnitude ~ pi.
    CHECK(std::abs(h.mode()) > 3.0);
}

TEST_CASE("circular neighbor wrap: peak in bin 0 counts the last bin as a neighbor") {
    Histogram1D h;
    h.configure(circular(8, 0.0, 8.0, /*split=*/false));   // width 1, centers k+0.5
    // Put the peak in bin 0 and mass in the wrap-around neighbor (bin 7).
    h.add(0.5, 5.0);   // bin 0
    h.add(7.5, 3.0);   // bin 7 (wrap neighbor of bin 0)
    h.add(1.5, 1.0);   // bin 1
    CHECK(h.peak_bin() == 0);
    // Confidence = (bin0 + bin7 + bin1) / total = (5 + 3 + 1)/9 = 1.0 (all mass).
    CHECK(h.confidence() == doctest::Approx(1.0));
}

// ---------------------------------------------------------------------------
// Non-circular boundary bin: clamps out-of-range, peak at the edge, no wrap.
// ---------------------------------------------------------------------------
TEST_CASE("non-circular clamps out-of-range values to the boundary bin") {
    Histogram1D h;
    h.configure(linear(10, 0.0, 10.0, /*split=*/false));
    h.add(-100.0, 3.0);   // clamps into bin 0
    h.add(100.0, 2.0);    // clamps into bin 9
    h.add(0.5, 1.0);      // bin 0
    CHECK(h.peak_bin() == 0);            // bin0 = 4 > bin9 = 2
    CHECK(h.mode() >= 0.0);
    CHECK(h.mode() <= 10.0);
}

TEST_CASE("non-circular peak at a boundary bin: missing neighbor contributes 0") {
    Histogram1D h;
    h.configure(linear(10, 0.0, 10.0, /*split=*/false, /*sub=*/true));
    // Peak hard against the lower boundary (bin 0); only the right neighbor exists.
    h.add(0.5, 10.0);   // bin 0
    h.add(1.5, 4.0);    // bin 1
    CHECK(h.peak_bin() == 0);
    // No parabola for a boundary bin (missing left neighbor) -> bin center.
    CHECK(h.mode() == doctest::Approx(0.5));
    // Confidence = (bin0 + bin1 + [missing left = 0]) / total = (10+4)/14 = 1.0.
    CHECK(h.confidence() == doctest::Approx(1.0));
}

// ---------------------------------------------------------------------------
// Confidence is always within [0, 1] across a mixed sequence; reset clears.
// ---------------------------------------------------------------------------
TEST_CASE("confidence stays within [0,1] and reset clears state") {
    Histogram1D h;
    h.configure(linear(64, -3.0, 3.0));
    for (int k = 0; k < 500; ++k) {
        const Scalar v = -3.0 + 6.0 * (static_cast<Scalar>((k * 37) % 100) / 100.0);
        h.add(v, 0.5 + 0.5 * static_cast<Scalar>(k % 3));
        const Scalar c = h.confidence();
        CHECK(c >= 0.0);
        CHECK(c <= 1.0 + 1e-12);
    }
    CHECK_FALSE(h.empty());

    h.reset();
    CHECK(h.empty());
    CHECK(h.total() == doctest::Approx(0.0));
    CHECK(h.peak_bin() == -1);
    CHECK(h.confidence() == doctest::Approx(0.0));
    CHECK(h.mode() == doctest::Approx(0.0));   // midpoint of [-3, 3]
}

// ---------------------------------------------------------------------------
// add() on an unconfigured histogram is a safe no-op.
// ---------------------------------------------------------------------------
TEST_CASE("add before configure is a no-op") {
    Histogram1D h;
    h.add(1.0, 1.0);
    CHECK(h.empty());
    CHECK(h.peak_bin() == -1);
}

// ---------------------------------------------------------------------------
// Non-positive weights are ignored.
// ---------------------------------------------------------------------------
TEST_CASE("non-positive weights are ignored") {
    Histogram1D h;
    h.configure(linear(10, 0.0, 10.0));
    h.add(5.0, 0.0);
    h.add(5.0, -1.0);
    CHECK(h.empty());
    h.add(5.0, 2.0);
    CHECK_FALSE(h.empty());
    CHECK(h.total() == doctest::Approx(2.0));
}
