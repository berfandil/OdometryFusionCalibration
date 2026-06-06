// ofc/core/histogram.cpp — the 1-D robustness primitive (Slice 4).
//
// STRICT CORE: all storage is fixed-capacity member arrays sized in configure();
// add()/mode()/confidence() allocate nothing and run loops bounded by the active
// bin count or the sliding-window size. No exceptions; status-code returns; double.
//
// Layout & math (mirrors the header + the PLAN):
//   * Bins: `nbins_` uniform bins over [range_min, range_max]. width = span/nbins.
//     bin_center(i) = range_min + (i + 0.5) * width. The fractional bin coordinate
//     of a value is f = (value - range_min)/width - 0.5, so f == i at bin i's
//     center; round(f) is the nearest bin, floor(f)..floor(f)+1 the split pair.
//   * Circular: values fold into [range_min, range_max) by modulo; bin neighbors
//     wrap mod nbins. Non-circular: values clamp to the range; neighbors clamp
//     (out-of-range neighbor mass = 0).
//   * Voting: vote_split distributes the weight linearly between the two nearest
//     bins (wrap-aware when circular, edge-clamped when not); else the single
//     nearest bin gets the full weight.
//   * Aging (D9): Decay scales every bin (and the running total) by decay_gamma
//     before each deposit. SlidingK keeps a ring of the last `sliding_k` votes and
//     subtracts the evicted oldest vote's exact bin contributions when the
//     (K+1)-th arrives — fixed memory, exact forgetting.
//   * Sub-bin peak: parabolic interpolation around the argmax i using its two
//     neighbors (wrap-aware if circular):
//       offset = 0.5*(h[i-1]-h[i+1]) / (h[i-1]-2 h[i]+h[i+1])  (guarded denom),
//       mode = bin_center(i) + offset*width, wrapped into range if circular. For a
//     non-circular boundary bin (a neighbor missing) we fall back to the bin
//     center. subbin==false always returns the bin center.
//   * Confidence = peak concentration = (peak bin + its two immediate neighbors) /
//     total, clamped to [0,1]; 0 when empty.
#include "ofc/core/histogram.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
// Mass below which the histogram is considered empty (and totals are treated as
// zero for division guards). Small relative to any meaningful vote weight.
constexpr Scalar kEmptyEps = Scalar(1e-12);
// Denominator guard for the parabolic peak (a ~flat or non-concave triple).
constexpr Scalar kParabEps = Scalar(1e-12);
} // namespace

Status Histogram1D::configure(const HistogramConfig& cfg) {
    if (cfg.bins < 4 || cfg.bins > kMaxBins) return Status::OutOfRange;
    if (!(cfg.range_max > cfg.range_min))    return Status::InvalidConfig;

    if (cfg.aging == Aging::Decay) {
        if (!(cfg.decay_gamma > Scalar(0)) || !(cfg.decay_gamma < Scalar(1)))
            return Status::OutOfRange;
    } else { // SlidingK
        if (cfg.sliding_k < 1 || cfg.sliding_k > kMaxSlidingK)
            return Status::OutOfRange;
    }

    nbins_       = cfg.bins;
    range_min_   = cfg.range_min;
    range_max_   = cfg.range_max;
    circular_    = cfg.circular;
    aging_       = cfg.aging;
    decay_gamma_ = cfg.decay_gamma;
    sliding_k_   = cfg.sliding_k;
    vote_split_  = cfg.vote_split;
    subbin_      = cfg.subbin;
    configured_  = true;

    reset();
    return Status::Ok;
}

void Histogram1D::reset() {
    for (int i = 0; i < nbins_; ++i) bins_[i] = Scalar(0);
    total_      = Scalar(0);
    ring_head_  = 0;
    ring_count_ = 0;
}

Scalar Histogram1D::fold(Scalar value) const {
    const Scalar span = range_max_ - range_min_;
    if (circular_) {
        // Wrap into [range_min, range_max).
        Scalar x = std::fmod(value - range_min_, span);
        if (x < Scalar(0)) x += span;     // fmod can be negative
        return range_min_ + x;
    }
    // Clamp into [range_min, range_max].
    if (value < range_min_) return range_min_;
    if (value > range_max_) return range_max_;
    return value;
}

int Histogram1D::wrap_bin(int i) const {
    if (circular_) {
        int m = i % nbins_;
        if (m < 0) m += nbins_;
        return m;
    }
    if (i < 0)        return 0;
    if (i >= nbins_)  return nbins_ - 1;
    return i;
}

Scalar Histogram1D::mass_at(int i) const {
    if (circular_) {
        int m = i % nbins_;
        if (m < 0) m += nbins_;
        return bins_[m];
    }
    if (i < 0 || i >= nbins_) return Scalar(0);   // missing neighbor
    return bins_[i];
}

void Histogram1D::plan_vote(Scalar value, Scalar weight, Vote& out) const {
    const Scalar v = fold(value);
    const Scalar w = width();
    // Fractional bin coordinate: f == i at the center of bin i.
    const Scalar f = (v - range_min_) / w - Scalar(0.5);

    if (!vote_split_) {
        // Single nearest bin.
        int i = static_cast<int>(std::lround(f));
        i = wrap_bin(i);
        out.bin_lo = i;
        out.bin_hi = i;
        out.w_lo   = weight;
        out.w_hi   = Scalar(0);
        return;
    }

    // Linear split between the two nearest bin centers.
    const Scalar fl   = std::floor(f);
    int   lo  = static_cast<int>(fl);
    Scalar frac = f - fl;                 // in [0, 1): share going to the higher bin
    int   hi  = lo + 1;

    lo = wrap_bin(lo);
    hi = wrap_bin(hi);

    // When non-circular and the pair straddles a boundary, wrap_bin clamps both to
    // the same edge bin; the two fractions then simply re-merge into that bin.
    out.bin_lo = lo;
    out.bin_hi = hi;
    out.w_lo   = weight * (Scalar(1) - frac);
    out.w_hi   = weight * frac;
}

void Histogram1D::apply_vote(const Vote& v, Scalar sign) {
    if (v.bin_lo == v.bin_hi) {
        const Scalar dw = sign * (v.w_lo + v.w_hi);
        bins_[v.bin_lo] += dw;
        total_          += dw;
    } else {
        bins_[v.bin_lo] += sign * v.w_lo;
        bins_[v.bin_hi] += sign * v.w_hi;
        total_          += sign * (v.w_lo + v.w_hi);
    }
}

void Histogram1D::add(Scalar value, Scalar weight) {
    if (!configured_)          return;
    if (!(weight > Scalar(0))) return;   // ignore zero/negative weight votes

    if (aging_ == Aging::Decay) {
        // Scale all live mass by gamma before depositing the fresh vote.
        for (int i = 0; i < nbins_; ++i) bins_[i] *= decay_gamma_;
        total_ *= decay_gamma_;

        Vote v;
        plan_vote(value, weight, v);
        apply_vote(v, Scalar(1));
        return;
    }

    // SlidingK: deposit the new vote and evict the oldest if the window is full.
    Vote v;
    plan_vote(value, weight, v);

    if (ring_count_ == sliding_k_) {
        // Subtract the oldest vote's exact contributions, then overwrite its slot.
        apply_vote(ring_[ring_head_], Scalar(-1));
        ring_[ring_head_] = v;
        ring_head_ = (ring_head_ + 1) % sliding_k_;
    } else {
        const int slot = (ring_head_ + ring_count_) % sliding_k_;
        ring_[slot] = v;
        ++ring_count_;
    }
    apply_vote(v, Scalar(1));

    // Numerical hygiene: repeated +/- of equal masses can leave a tiny negative
    // residual in an emptied bin. Clamp the running total's float noise away when
    // the window is empty so empty() stays exact.
    if (ring_count_ == 0) total_ = Scalar(0);
}

bool Histogram1D::empty() const {
    return total_ <= kEmptyEps;
}

Scalar Histogram1D::total() const {
    return total_ > Scalar(0) ? total_ : Scalar(0);
}

int Histogram1D::peak_bin() const {
    if (empty()) return -1;
    int best = 0;
    Scalar best_mass = bins_[0];
    for (int i = 1; i < nbins_; ++i) {
        if (bins_[i] > best_mass) {
            best_mass = bins_[i];
            best = i;
        }
    }
    return best;
}

Scalar Histogram1D::mode() const {
    if (empty()) {
        return Scalar(0.5) * (range_min_ + range_max_);   // documented midpoint
    }
    const int i = peak_bin();
    const Scalar center = bin_center(i);

    if (!subbin_) return center;

    // Parabolic interpolation needs both neighbors. Non-circular boundary bins
    // miss one -> fall back to the bin center.
    const bool have_left  = circular_ || (i - 1 >= 0);
    const bool have_right = circular_ || (i + 1 < nbins_);
    if (!have_left || !have_right) return center;

    const Scalar hm = mass_at(i - 1);
    const Scalar h0 = mass_at(i);
    const Scalar hp = mass_at(i + 1);
    const Scalar denom = hm - Scalar(2) * h0 + hp;
    if (std::abs(denom) < kParabEps) return center;   // ~flat / non-concave triple

    Scalar offset = Scalar(0.5) * (hm - hp) / denom;
    // A well-formed peak gives |offset| <= 0.5; clamp against numerical overshoot.
    if (offset >  Scalar(0.5)) offset =  Scalar(0.5);
    if (offset < -Scalar(0.5)) offset = -Scalar(0.5);

    Scalar m = center + offset * width();
    if (circular_) m = fold(m);   // keep the continuous mode inside the range
    return m;
}

Scalar Histogram1D::confidence() const {
    if (empty()) return Scalar(0);
    const int i = peak_bin();
    const Scalar peak_mass = mass_at(i - 1) + mass_at(i) + mass_at(i + 1);
    Scalar c = peak_mass / total_;
    if (c < Scalar(0)) c = Scalar(0);
    if (c > Scalar(1)) c = Scalar(1);
    return c;
}

} // namespace ofc
