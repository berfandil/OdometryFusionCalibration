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

// ---------------------------------------------------------------------------
// SlidingK float drift: over many evictions, subtracting stored deposits must
// not let any bin settle negative, and total() must stay equal to the sum of
// bin masses. Runs >> 10*K mixed-value adds to exercise the +/- residual path.
// (Fails before the per-evict clamp/rebuild fix: a bin drifts slightly negative
// and total_ desyncs from sum(bins_).)
// ---------------------------------------------------------------------------
TEST_CASE("sliding-K stays non-negative and total tracks sum(bins) over many evictions") {
    const int K = 17;
    HistogramConfig c = linear(64, -3.0, 3.0, /*split=*/true);
    c.aging = Aging::SlidingK;
    c.sliding_k = K;
    Histogram1D h;
    h.configure(c);

    // >> 10*K mixed-value adds with non-representable split fractions and large
    // weights, so the same bin sees a long, constantly-changing interleave of
    // +deposit / -deposit of *different* fractional masses. Non-associative float
    // addition then accumulates a residual in the running total_: it desyncs from
    // sum(bins_) (here by ~1e-9 over ~850k evictions) unless every eviction
    // rebuilds total_ = sum(bins_) and clamps tiny-negative bin residuals to 0.
    const int N = 50000 * K;   // 850000 adds, ~849983 evictions
    for (int k = 0; k < N; ++k) {
        // Irrational-ish walk over the range so values never repeat exactly and
        // the linear split lands on non-binary fractions. Large weight magnifies
        // the absolute float residual well above the round-off floor.
        const Scalar t = std::fmod(static_cast<Scalar>(k) * 0.61803398875, 1.0);
        const Scalar v = -3.0 + 6.0 * t;
        const Scalar w = (1.0 / 3.0 + (static_cast<Scalar>(k % 7) + 1.0) * 0.07) * 1e6;
        h.add(v, w);
    }

    // Invariant 1: no bin reports negative. The clamp drives any tiny-negative
    // residual to exactly 0, so every bin mass is >= 0.
    Scalar sum = 0.0;
    for (int i = 0; i < h.bins(); ++i) {
        const Scalar m = h.bin_mass(i);
        CHECK(m >= 0.0);
        sum += m;
    }

    // Invariant 2 (the load-bearing one): total() is EXACTLY the sum of bin
    // masses. The per-evict rebuild recomputes total_ = sum(bins_) the same way,
    // so the delta is 0. Without the rebuild, the running total_ has drifted
    // ~1e-9 from sum(bins_): this absolute-tolerance check fails before the fix.
    CHECK(std::abs(h.total() - sum) < 1e-12);

    // Invariant 3: total() still equals the exact analytic mass of the last K
    // live votes (window holds exactly K). The rebuild keeps total_ drift-free.
    Scalar exact_total = 0.0;
    for (int k = N - K; k < N; ++k) {
        exact_total += (1.0 / 3.0 + (static_cast<Scalar>(k % 7) + 1.0) * 0.07) * 1e6;
    }
    CHECK(h.total() == doctest::Approx(exact_total).epsilon(1e-9));

    // Confidence remains a valid probability throughout.
    const Scalar conf = h.confidence();
    CHECK(conf >= 0.0);
    CHECK(conf <= 1.0);
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
// Circular wrap-seam, ANALYTIC: pin the exact wrap-aware split fractions of a
// vote straddling the seam, and the exact parabolic sub-bin offset of a peak
// whose left neighbour is the wrapped last bin. (Strengthens the "near the
// seam" cluster test with planted masses and closed-form expected values.)
// ---------------------------------------------------------------------------
TEST_CASE("circular wrap-seam: exact split fractions and exact parabolic offset") {
    // 8 bins over [0, 8), circular. width 1, centers k+0.5. Seam between bin 7
    // (center 7.5) and bin 0 (center 0.5, i.e. wrapped center 8.5).

    // --- Part A: wrap-aware linear-split fractions straddling the seam. --------
    // A value v = 7.75 lies 0.25 of the way from center 7.5 (bin 7) toward the
    // wrapped center 8.5 (bin 0): f = v - 0.5 = 7.25, floor = 7, frac = 0.25.
    //   -> w_lo = w*(1-frac) = 0.75*w into bin 7
    //   -> w_hi = w*frac     = 0.25*w into bin 0 (hi = 8 wraps to 0)
    {
        Histogram1D h;
        h.configure(circular(8, 0.0, 8.0, /*split=*/true, /*sub=*/false));
        const Scalar w = 4.0;
        h.add(7.75, w);                       // straddles the seam
        CHECK(h.total() == doctest::Approx(w));          // mass conserved across seam
        CHECK(h.bin_mass(7) == doctest::Approx(0.75 * w)); // 3.0
        CHECK(h.bin_mass(0) == doctest::Approx(0.25 * w)); // 1.0
        CHECK(h.peak_bin() == 7);             // larger fraction stays in bin 7
        // Confidence = (peak 7 + neighbours 6 and 0)/total = (3 + 0 + 1)/4 = 1.0.
        CHECK(h.confidence() == doctest::Approx(1.0));
    }

    // --- Part B: parabolic sub-bin offset using the WRAPPED neighbour. ---------
    // Peak in bin 0; its left neighbour is bin 7 (wrapped). Plant masses directly
    // (split off, sub on) so the parabola is exact:
    //   hm = mass(bin 7) = 8, h0 = mass(bin 0) = 10, hp = mass(bin 1) = 2.
    //   offset = 0.5*(hm - hp)/(hm - 2*h0 + hp)
    //          = 0.5*(8 - 2)/(8 - 20 + 2) = 0.5*6/(-10) = -0.3
    //   mode = center(0) + offset*width = 0.5 + (-0.3)*1 = 0.2  (toward the seam)
    {
        Histogram1D h;
        HistogramConfig c = circular(8, 0.0, 8.0, /*split=*/false, /*sub=*/true);
        c.aging = Aging::SlidingK;
        c.sliding_k = 100;                    // exact accumulation, no decay
        h.configure(c);
        for (int k = 0; k < 8;  ++k) h.add(7.5);   // bin 7 (wrapped left neighbour)
        for (int k = 0; k < 10; ++k) h.add(0.5);   // bin 0 (the peak)
        for (int k = 0; k < 2;  ++k) h.add(1.5);   // bin 1 (right neighbour)

        CHECK(h.peak_bin() == 0);
        CHECK(h.bin_mass(7) == doctest::Approx(8.0));
        CHECK(h.bin_mass(0) == doctest::Approx(10.0));
        CHECK(h.bin_mass(1) == doctest::Approx(2.0));

        const Scalar hm = 8.0, h0 = 10.0, hp = 2.0;
        const Scalar expected_offset = 0.5 * (hm - hp) / (hm - 2.0 * h0 + hp);
        CHECK(expected_offset == doctest::Approx(-0.3));
        const Scalar expected_mode = 0.5 + expected_offset * 1.0;   // 0.2
        CHECK(h.mode() == doctest::Approx(expected_mode));
        CHECK(h.mode() == doctest::Approx(0.2));
    }

    // --- Part C: mirror case, peak in the LAST bin with a wrapped RIGHT
    // neighbour (bin 0). Pins the seam neighbour read from the other side.
    //   hm = mass(bin 6) = 2, h0 = mass(bin 7) = 10, hp = mass(bin 0) = 8.
    //   offset = 0.5*(2 - 8)/(2 - 20 + 8) = 0.5*(-6)/(-10) = 0.3
    //   mode = center(7) + 0.3 = 7.5 + 0.3 = 7.8  (toward the seam, stays in range)
    {
        Histogram1D h;
        HistogramConfig c = circular(8, 0.0, 8.0, /*split=*/false, /*sub=*/true);
        c.aging = Aging::SlidingK;
        c.sliding_k = 100;
        h.configure(c);
        for (int k = 0; k < 2;  ++k) h.add(6.5);   // bin 6 (left neighbour)
        for (int k = 0; k < 10; ++k) h.add(7.5);   // bin 7 (the peak)
        for (int k = 0; k < 8;  ++k) h.add(0.5);   // bin 0 (wrapped right neighbour)
        CHECK(h.peak_bin() == 7);
        CHECK(h.bin_mass(6) == doctest::Approx(2.0));
        CHECK(h.bin_mass(7) == doctest::Approx(10.0));
        CHECK(h.bin_mass(0) == doctest::Approx(8.0));   // wrapped right neighbour
        const Scalar hm = 2.0, h0 = 10.0, hp = 8.0;
        const Scalar expected_offset = 0.5 * (hm - hp) / (hm - 2.0 * h0 + hp);
        CHECK(expected_offset == doctest::Approx(0.3));
        CHECK(h.mode() == doctest::Approx(7.5 + 0.3 * 1.0));   // 7.8
    }
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

    // Same mixed stress under a SlidingK config, run well past the window size so
    // confidence is exercised over many evictions (the drift path), not just the
    // default Decay aging.
    Histogram1D hk;
    HistogramConfig ck = linear(64, -3.0, 3.0);
    ck.aging = Aging::SlidingK;
    ck.sliding_k = 25;
    hk.configure(ck);
    for (int k = 0; k < 500; ++k) {           // 20x the window -> many evictions
        const Scalar v = -3.0 + 6.0 * (static_cast<Scalar>((k * 37) % 100) / 100.0);
        hk.add(v, 0.5 + 0.5 * static_cast<Scalar>(k % 3));
        const Scalar c = hk.confidence();
        CHECK(c >= 0.0);
        CHECK(c <= 1.0 + 1e-12);
    }
    CHECK_FALSE(hk.empty());

    h.reset();
    CHECK(h.empty());
    CHECK(h.total() == doctest::Approx(0.0));
    CHECK(h.peak_bin() == -1);
    CHECK(h.confidence() == doctest::Approx(0.0));
    CHECK(h.mode() == doctest::Approx(0.0));   // midpoint of [-3, 3]
}

// ===========================================================================
// Slice 16 — opt-in centroid sub-bin readout (subbin_centroid).
// mode() returns the mass-weighted centroid over the peak bin and its two
// immediate neighbors instead of the parabolic interpolation. Default OFF
// keeps every existing path byte-identical.
// ===========================================================================

namespace {
// Same linear/circular helpers with the centroid readout switched on.
HistogramConfig linear_centroid(int n, Scalar lo, Scalar hi, bool split = true,
                                bool sub = true) {
    HistogramConfig c = linear(n, lo, hi, split, sub);
    c.subbin_centroid = true;
    return c;
}

HistogramConfig circular_centroid(int n, Scalar lo, Scalar hi, bool split = true,
                                  bool sub = true) {
    HistogramConfig c = circular(n, lo, hi, split, sub);
    c.subbin_centroid = true;
    return c;
}
} // namespace

// ---------------------------------------------------------------------------
// Acceptance 1: a SINGLE vote at value v with vote_split lands as a linear
// split across two bins; the 3-bin centroid reconstructs v EXACTLY (<= 1e-12
// of width) at ANY sub-bin position. The parabola does not — pin its >= 0.1-bin
// error at the quarter-bin position (the documented Slice-16 rationale).
// ---------------------------------------------------------------------------
TEST_CASE("centroid readout reconstructs a single split vote exactly at any sub-bin position") {
    // 10 bins over [0, 10): width 1, bin 4 center at 4.5. Sub-bin fractions f
    // measured from the bin-4 center: POSITIVE toward bin 5 (split mass in the
    // RIGHT neighbor, the hp-path) and NEGATIVE toward bin 3 (split mass in the
    // LEFT neighbor, the hm-path), so a sign error in the (hp - hm) numerator
    // is caught by this headline acceptance sweep itself.
    const Scalar fractions[] = {-0.4, -0.25, -0.1, 0.0, 0.1, 0.25, 0.4, 0.5};
    for (const Scalar f : fractions) {
        Histogram1D h;
        REQUIRE(h.configure(linear_centroid(10, 0.0, 10.0)) == Status::Ok);
        const Scalar v = 4.5 + f;          // f bins above/below the bin-4 center
        h.add(v);
        CHECK(std::abs(h.mode() - v) <= 1e-12);   // exact (width = 1)
    }

    // The parabola at the SAME quarter-bin vote reads ~0.1 bins instead of
    // 0.25 — a systematic pull toward the peak-bin center of >= 0.1 bins.
    {
        Histogram1D h;
        REQUIRE(h.configure(linear(10, 0.0, 10.0)) == Status::Ok);   // parabolic (default)
        const Scalar v = 4.75;             // quarter-bin position
        h.add(v);
        const Scalar err_bins = std::abs(h.mode() - v) / 1.0;        // width = 1
        CHECK(err_bins >= 0.1);            // pinned bias (analytic: reads 0.1, true 0.25)
    }
}

// ---------------------------------------------------------------------------
// Acceptance 2: a symmetric multi-vote spread -> centroid == the true mean
// (votes split linearly, the centroid is linear in mass, all mass within
// peak±1); the parabola misreads the same spread by >> the centroid error.
// ---------------------------------------------------------------------------
TEST_CASE("centroid reads a symmetric vote spread at its true mean; parabola does not") {
    const Scalar mu = 4.75;                // quarter-bin position inside bin 4
    const Scalar d  = 0.3;                 // spread stays within bins 3..5

    // SlidingK (window > vote count) -> exact accumulation, no decay between
    // the votes, so the analytic mean recovery holds to round-off.
    HistogramConfig cc = linear_centroid(10, 0.0, 10.0);
    cc.aging = Aging::SlidingK;
    cc.sliding_k = 100;
    Histogram1D hc;
    REQUIRE(hc.configure(cc) == Status::Ok);
    hc.add(mu - d);
    hc.add(mu + d);
    hc.add(mu);
    CHECK(std::abs(hc.mode() - mu) <= 1e-12);   // exact mean recovery

    HistogramConfig cp = linear(10, 0.0, 10.0);
    cp.aging = Aging::SlidingK;
    cp.sliding_k = 100;
    Histogram1D hp;
    REQUIRE(hp.configure(cp) == Status::Ok);    // parabolic
    hp.add(mu - d);
    hp.add(mu + d);
    hp.add(mu);
    CHECK(std::abs(hp.mode() - mu) > 0.1);      // the parabola's pull-to-center bias
}

// ---------------------------------------------------------------------------
// Acceptance 3: circular seam — a peak straddling the wrap reads continuously.
// Neighbor MASS wraps; neighbor CENTERS unwrap (center(p) ± width); the
// centroid folds back into range.
// ---------------------------------------------------------------------------
TEST_CASE("centroid reads continuously across the circular wrap seam") {
    // 8 bins over [0, 8), width 1, seam between bin 7 (center 7.5) and bin 0
    // (center 0.5 == unwrapped 8.5). Single split votes at sub-bin positions on
    // BOTH sides of the seam must read back exactly (modulo the fold).
    const Scalar votes[] = {7.6, 7.75, 7.95, 8.05, 0.1, 0.4};
    for (const Scalar v : votes) {
        Histogram1D h;
        REQUIRE(h.configure(circular_centroid(8, 0.0, 8.0)) == Status::Ok);
        h.add(v);
        const Scalar m = h.mode();
        // In range...
        CHECK(m >= 0.0);
        CHECK(m < 8.0);
        // ...and circularly exact: distance to the folded vote ~ 0.
        Scalar dist = std::fmod(std::abs(m - v), 8.0);
        if (dist > 4.0) dist = 8.0 - dist;
        CHECK(dist <= 1e-12);
    }

    // A 50/50 seam-straddling pair (vote exactly AT the seam value 0 == 8):
    // mass splits bin 7 / bin 0 equally; the centroid must read the seam (0.0),
    // not a mid-range value.
    {
        Histogram1D h;
        REQUIRE(h.configure(circular_centroid(8, 0.0, 8.0)) == Status::Ok);
        h.add(0.0);                        // folds: 0.5 to bin 7, 0.5 to bin 0
        const Scalar m = h.mode();
        Scalar dist = std::fmod(std::abs(m - 0.0), 8.0);
        if (dist > 4.0) dist = 8.0 - dist;
        CHECK(dist <= 1e-12);
    }
}

// ---------------------------------------------------------------------------
// Acceptance 4: non-circular boundary peak — the missing neighbor contributes
// 0 mass (and its unwrapped center to nothing); the one-sided centroid is
// still exact for an in-range split pair, with NO fall-back-to-center cliff
// (the parabola falls back to the bin center there). No NaN ever.
// ---------------------------------------------------------------------------
TEST_CASE("centroid at a non-circular boundary peak degrades gracefully (no cliff, no NaN)") {
    // (a) Split vote near the LOWER boundary, pair fully in range: exact.
    {
        Histogram1D h;
        REQUIRE(h.configure(linear_centroid(10, 0.0, 10.0)) == Status::Ok);
        h.add(0.7);                        // splits bins 0/1; peak = boundary bin 0
        CHECK(h.peak_bin() == 0);
        CHECK(std::abs(h.mode() - 0.7) <= 1e-12);   // graceful one-sided centroid

        // The parabolic path at the same vote falls back to the bin center 0.5
        // (missing left neighbor) — the cliff the centroid removes.
        Histogram1D hp;
        REQUIRE(hp.configure(linear(10, 0.0, 10.0)) == Status::Ok);
        hp.add(0.7);
        CHECK(hp.mode() == doctest::Approx(0.5));
    }
    // (b) Split vote near the UPPER boundary: exact from the other side.
    {
        Histogram1D h;
        REQUIRE(h.configure(linear_centroid(10, 0.0, 10.0)) == Status::Ok);
        h.add(9.3);                        // splits bins 8/9; peak = boundary bin 9
        CHECK(h.peak_bin() == 9);
        CHECK(std::abs(h.mode() - 9.3) <= 1e-12);
    }
    // (c) ALL mass in the boundary bin (a clamped out-of-range vote): the
    // centroid is the bin center; finite, no NaN.
    {
        Histogram1D h;
        REQUIRE(h.configure(linear_centroid(10, 0.0, 10.0)) == Status::Ok);
        h.add(-5.0);                       // clamps wholly into bin 0
        const Scalar m = h.mode();
        CHECK(std::isfinite(m));
        CHECK(m == doctest::Approx(0.5));  // one-sided centroid of a point mass
    }
}

// ---------------------------------------------------------------------------
// Acceptance 5 (default-off): subbin_centroid == false leaves the parabolic
// path bit-identical; subbin == false still returns the bin center even with
// the centroid flag set; an empty histogram still reads the midpoint.
// ---------------------------------------------------------------------------
TEST_CASE("subbin_centroid honors the existing switches: default-off, subbin=false, empty") {
    // Default-off: the planted (2, 8, 6) triple reads the EXACT analytic
    // parabola value — bitwise the same expression as before the flag existed.
    {
        Histogram1D h;
        HistogramConfig c = linear(10, 0.0, 10.0, /*split=*/false, /*sub=*/true);
        c.aging = Aging::SlidingK;
        c.sliding_k = 100;
        REQUIRE(h.configure(c) == Status::Ok);
        for (int k = 0; k < 2; ++k) h.add(3.5);
        for (int k = 0; k < 8; ++k) h.add(4.5);
        for (int k = 0; k < 6; ++k) h.add(5.5);
        const Scalar expected_offset = 0.5 * (2.0 - 6.0) / (2.0 - 2.0 * 8.0 + 6.0);
        CHECK(h.mode() == 4.5 + expected_offset * 1.0);   // EXACT equality, not Approx
    }
    // subbin == false dominates: the centroid flag alone must not enable a
    // sub-bin readout.
    {
        Histogram1D h;
        HistogramConfig c = linear_centroid(10, 0.0, 10.0, /*split=*/true, /*sub=*/false);
        REQUIRE(h.configure(c) == Status::Ok);
        h.add(4.75);
        CHECK(h.mode() == doctest::Approx(4.5));   // bin center
    }
    // Empty: unchanged midpoint.
    {
        Histogram1D h;
        REQUIRE(h.configure(linear_centroid(64, 2.0, 6.0)) == Status::Ok);
        CHECK(h.mode() == doctest::Approx(4.0));
    }
}

// ---------------------------------------------------------------------------
// Config-space interactions of subbin_centroid (Slice 16 review NITs):
//   (a) vote_split == false + centroid: a single UNSPLIT vote carries no
//       sub-bin encoding in the masses (the whole weight lands in the nearest
//       bin), so the centroid must read the BIN CENTER — not invent an offset.
//   (b) Decay aging: gamma scales every bin UNIFORMLY, so it cancels in the
//       centroid's mass ratio — two identical decayed votes still read the
//       vote value exactly (decay-invariance), and a heterogeneous two-vote
//       history reads the analytic DECAY-WEIGHTED mean. All five production
//       histograms may run Decay, so this invariance is pinned here.
// ---------------------------------------------------------------------------
TEST_CASE("centroid with vote_split=false reads the bin center; centroid under Decay is gamma-invariant") {
    // (a) split=false + centroid: vote at 4.75 lands wholly in bin 4 (nearest);
    // both neighbors carry 0 mass -> centroid == bin center 4.5.
    {
        Histogram1D h;
        REQUIRE(h.configure(linear_centroid(10, 0.0, 10.0, /*split=*/false,
                                            /*sub=*/true)) == Status::Ok);
        h.add(4.75);
        CHECK(h.mode() == doctest::Approx(4.5));   // bin center, no sub-bin offset
    }
    // (b1) Decay invariance: two IDENTICAL split votes — the older one is
    // gamma-scaled before the second deposit, but gamma scales BOTH of its bin
    // deposits, so the (hp - hm)/sum ratio (and hence the centroid) is
    // unchanged: still exactly the vote value, independent of gamma.
    {
        HistogramConfig c = linear_centroid(10, 0.0, 10.0);
        c.aging = Aging::Decay;
        c.decay_gamma = 0.5;               // aggressive gamma: invariance, not luck
        Histogram1D h;
        REQUIRE(h.configure(c) == Status::Ok);
        h.add(4.75);
        h.add(4.75);                       // decays the first vote by 0.5 first
        CHECK(std::abs(h.mode() - 4.75) <= 1e-12);
    }
    // (b2) Heterogeneous decayed history: votes v1 = 4.6 then v2 = 4.9 under
    // gamma = 0.5. Masses: bin4 = 0.9*g + 0.6, bin5 = 0.1*g + 0.4 (both votes
    // split bins 4/5, all mass within peak+-1), so the centroid is the analytic
    // decay-weighted mean (g*v1 + v2)/(g + 1) = (0.5*4.6 + 4.9)/1.5 = 4.8.
    {
        HistogramConfig c = linear_centroid(10, 0.0, 10.0);
        c.aging = Aging::Decay;
        c.decay_gamma = 0.5;
        Histogram1D h;
        REQUIRE(h.configure(c) == Status::Ok);
        h.add(4.6);
        h.add(4.9);
        const Scalar expected = (0.5 * 4.6 + 4.9) / 1.5;   // 4.8
        CHECK(std::abs(h.mode() - expected) <= 1e-12);
    }
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
