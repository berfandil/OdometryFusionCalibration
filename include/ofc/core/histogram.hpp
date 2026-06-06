// ofc/core/histogram.hpp — the shared 1-D robustness primitive (Slice 4).
//
// A reusable 1-D histogram that votes values with weights, ages old votes, and
// extracts a sub-bin peak (parabolic interpolation) plus a peak-concentration
// confidence. This is the single robustness estimator reused throughout
// calibration (DESIGN §6, DECISIONS D9/D10/D11): Phase-1 so(3) channels, Phase-2
// roll (circular) + xyz (linear), per-source scale, and time-offset. Slice 4 is
// the 1-D primitive only — multi-channel composition belongs to the calibration
// slices.
//
// STRICT CORE: storage is fixed-capacity with compile-time maxima; configure()
// binds the active sizes and clears state; add()/mode()/confidence() allocate
// nothing and run bounded loops (every scan is bounded by the active bin count
// or the sliding-window size). No exceptions; status-code returns; double math.
#ifndef OFC_CORE_HISTOGRAM_HPP
#define OFC_CORE_HISTOGRAM_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

namespace ofc {

// 1-D histogram primitive. Configure once from a HistogramConfig, then add()
// weighted votes. The peak is read as a continuous value via mode() and its
// reliability via confidence().
//
// Confidence definition (documented, per the brief): the PEAK CONCENTRATION =
//   (mass in the peak bin + its two immediate neighbours) / total mass,
// clamped to [0, 1]; 0 when the histogram is empty. Neighbours wrap when the
// range is circular; for a non-circular range an out-of-range neighbour
// contributes 0 (so a peak hard against a boundary can read slightly lower).
//
// Empty (total mass ~ 0): mode() returns the range MIDPOINT (range_min+range_max)/2,
// confidence() == 0, peak_bin() == -1, empty() == true.
class Histogram1D {
public:
    // Compile-time maxima (strict core, no heap). bins <= kMaxBins matches the
    // upper bound enforced by Config::validate() (CONFIG §8). The sliding-window
    // ring is capped at kMaxSlidingK; configure() rejects sliding_k > kMaxSlidingK
    // with Status::OutOfRange.
    static constexpr int kMaxBins     = 4096;
    static constexpr int kMaxSlidingK = 4096;

    Histogram1D() = default;

    // Preallocate the active bin layout and bind the aging mode. No heap occurs
    // here (storage is fixed-size members) nor in any later call. Validates:
    //   bins in [4, kMaxBins]                  -> else OutOfRange
    //   range_max > range_min (finite span)    -> else InvalidConfig
    //   decay_gamma in (0, 1)  (Decay mode)    -> else OutOfRange
    //   sliding_k in [1, kMaxSlidingK] (SlidingK mode) -> else OutOfRange
    // Clears all state on success (and on a re-configure).
    Status configure(const HistogramConfig& cfg);

    // Drop all votes (keeps the configuration). Total mass -> 0.
    void reset();

    // Vote `value` with `weight` (default 1). Applies aging first (Decay scales
    // all bins by decay_gamma before the deposit; SlidingK evicts the oldest vote
    // once more than sliding_k votes are live), then deposits the (optionally
    // split) vote. Out-of-range values are wrapped into range when circular, else
    // clamped. A non-positive weight is ignored (no-op). No-op if not configured.
    void add(Scalar value, Scalar weight = Scalar(1));

    bool   empty() const;       // total mass ~ 0 (within a tiny epsilon)
    Scalar total() const;       // sum of all bin masses (>= 0)
    int    peak_bin() const;    // argmax bin index, or -1 if empty
    Scalar mode() const;        // sub-bin peak as a continuous value (see header)
    Scalar confidence() const;  // peak concentration in [0, 1]; 0 if empty

    int    bins() const { return nbins_; }
    bool   configured() const { return configured_; }

    // Raw mass of bin `i` (0 for an out-of-range index). Read-only inspection of
    // the internal histogram, used by tests to assert invariants (no bin ever
    // negative; total() == sum of bin masses). No heap; bounds-checked.
    Scalar bin_mass(int i) const;

private:
    // One recorded vote, kept only in SlidingK mode so its exact bin contributions
    // can be subtracted when it is evicted. Stores up to two bin deposits (the
    // split halves; lo == hi with w_hi == 0 for a single-bin vote).
    struct Vote {
        int    bin_lo = 0;
        int    bin_hi = 0;
        Scalar w_lo   = Scalar(0);
        Scalar w_hi   = Scalar(0);
    };

    // Geometry helpers.
    Scalar width()    const { return (range_max_ - range_min_) / static_cast<Scalar>(nbins_); }
    Scalar bin_center(int i) const {
        return range_min_ + (static_cast<Scalar>(i) + Scalar(0.5)) * width();
    }
    // Wrap a value into the half-open [range_min, range_max) (circular) or clamp
    // to that half-open range (non-circular: below -> range_min; at or above
    // range_max -> just below range_max, i.e. into the last bin). Returns the
    // in-range value.
    Scalar fold(Scalar value) const;
    // Wrap-or-clamp a bin index to the active range.
    int    wrap_bin(int i) const;
    // Bin mass for the parabola / confidence: wraps when circular, returns 0 for
    // an out-of-range index when non-circular.
    Scalar mass_at(int i) const;

    // Deposit a vote of `weight` at `value` into `out` (one or two bins); fills
    // the deposit fractions. Used by both add() and SlidingK eviction bookkeeping.
    void   plan_vote(Scalar value, Scalar weight, Vote& out) const;
    void   apply_vote(const Vote& v, Scalar sign);   // sign = +1 deposit, -1 evict
    // SlidingK numerical hygiene after an eviction-subtract: clamp each tiny
    // negative bin residual to 0 and rebuild total_ = sum(bins_) so the running
    // total never drifts from the bin sum. O(nbins_), heap-free.
    void   scrub_after_evict();

    // Active configuration.
    int    nbins_      = 0;
    Scalar range_min_  = Scalar(0);
    Scalar range_max_  = Scalar(0);
    bool   circular_   = false;
    Aging  aging_      = Aging::Decay;
    Scalar decay_gamma_ = Scalar(0.999);
    int    sliding_k_  = 1000;
    bool   vote_split_ = true;
    bool   subbin_     = true;
    bool   configured_ = false;

    // Fixed-capacity storage (no heap).
    Scalar bins_[kMaxBins] = {};
    Scalar total_          = Scalar(0);

    // SlidingK ring of live votes (only used in SlidingK mode).
    Vote ring_[kMaxSlidingK];
    int  ring_head_  = 0;    // physical index of the oldest live vote
    int  ring_count_ = 0;    // number of live votes (<= sliding_k_)
};

} // namespace ofc
#endif // OFC_CORE_HISTOGRAM_HPP
