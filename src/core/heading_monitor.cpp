// ofc/core/heading_monitor.cpp — GPS-course heading-drift monitor (Slice 19c).
// See heading_monitor.hpp for the full pipeline + design contract.
#include "ofc/core/heading_monitor.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

using namespace heading_monitor_const;

namespace {

constexpr Scalar kPi  = Scalar(3.14159265358979323846);
constexpr Scalar k2Pi = Scalar(6.28318530717958647692);

// Wrap an angle to (-pi, pi]. Branch-free fmod form (matches the prototype's wrap_pi).
Scalar wrap_pi(Scalar a) {
    Scalar w = std::fmod(a + kPi, k2Pi);
    if (w < Scalar(0)) w += k2Pi;
    return w - kPi;
}

// Weighted median of `m` (value, weight) samples (prototype wmedian, fixed-capacity). Sorts
// indices by value, returns the value at the first index where the cumulative weight reaches
// half the total. `idx` is caller-owned scratch of length >= m. Precondition: m >= 1, every
// weight >= 0 with a positive total. Returns vals[idx[0]] for m == 1.
// Tie behavior: termination at the FIRST index whose cumulative weight reaches `half` picks the
// LOWER median when the half-weight mark falls exactly on a boundary (the first-half element
// wins) -- matches the prototype's np.searchsorted(cw, 0.5*cw[-1]) and is correct for slope.
Scalar weighted_median(const Scalar* vals, const Scalar* wts, int m, int* idx) {
    for (int k = 0; k < m; ++k) idx[k] = k;
    std::sort(idx, idx + m, [&](int a, int b) { return vals[a] < vals[b]; });
    Scalar total = Scalar(0);
    for (int k = 0; k < m; ++k) total += wts[idx[k]];
    const Scalar half = Scalar(0.5) * total;
    Scalar acc = Scalar(0);
    for (int k = 0; k < m; ++k) {
        acc += wts[idx[k]];
        if (acc >= half) return vals[idx[k]];
    }
    return vals[idx[m - 1]];
}

// Plain (unit-weight) median of `m` values, via the weighted form (scratch idx of length m).
Scalar median_of(const Scalar* vals, int m, int* idx) {
    for (int k = 0; k < m; ++k) idx[k] = k;
    std::sort(idx, idx + m, [&](int a, int b) { return vals[a] < vals[b]; });
    const int mid = m / 2;            // upper-median for even m (matches numpy on the odd
    return vals[idx[mid]];            // n>=3 consensus path the monitor actually uses)
}

// Block MEDIAN over `m` values held in caller-owned scratch `buf` (length >= m; CLOBBERED in
// place). Strict core: an in-place partial selection (std::nth_element -- no heap, no allocation;
// it partitions the fixed stack array) gives O(m) WCET with no per-block sort/alloc. Matches the
// prototype's np.median: the average of the two middle order statistics for even m, the middle
// one for odd m. Precondition: m >= 1.
Scalar block_median(Scalar* buf, int m) {
    const int mid = m / 2;
    std::nth_element(buf, buf + mid, buf + m);
    const Scalar hi = buf[mid];
    if (m % 2 == 1) return hi;
    // Even m: pair the upper-middle with the max of the lower half (the lower-middle statistic).
    Scalar lo = buf[0];
    for (int k = 1; k < mid; ++k) if (buf[k] > lo) lo = buf[k];
    return Scalar(0.5) * (lo + hi);
}

bool finite(Scalar x) { return std::isfinite(x); }

}  // namespace

void HeadingMonitor::reset() {
    have_fix_    = false;
    have_anchor_ = false;
    fix_t_       = Scalar(0);
    fix_pos_     = Vec3::Zero();
    anc_t_       = Scalar(0);
    anc_course_  = Scalar(0);
    t_valid_     = Scalar(0);
    anchors_     = 0;
    pairs_       = 0;
    for (int i = 0; i < kMaxSources; ++i) {
        fix_yaw_[i] = Scalar(0);
        fix_fwd_[i] = Scalar(0);
        anc_yaw_[i] = Scalar(0);
        track_[i]   = Track{};
    }
}

void HeadingMonitor::configure(int n) {
    n_ = (n < 0) ? 0 : (n > kMaxSources ? kMaxSources : n);
    reset();
}

void HeadingMonitor::finalize_block(Track& tr) {
    if (tr.cur_n == 0) return;
    // Reduce the open block to its MEDIAN (t, cum) (spec 1.3 / prototype). cur_n counts every
    // submitted sample but the buffer holds at most kBlockSamples of them (surplus dropped on
    // submit); reduce over what we retained. block_median clobbers its scratch, so median the
    // time buffer first (copying into scratch_med_), then the cumulative buffer in place.
    const int m = (tr.cur_n < kBlockSamples) ? tr.cur_n : kBlockSamples;
    Block blk;
    for (int k = 0; k < m; ++k) scratch_med_[k] = tr.cur_ts[k];
    blk.t = block_median(scratch_med_, m);
    blk.c = block_median(tr.cur_cs, m);   // cur_cs is reset below, safe to clobber in place
    // Pair this block with every earlier block of the segment whose baseline >= min_base.
    for (int q = 0; q < tr.block_count; ++q) {
        const int ix = (tr.block_head + q) % kSlopeReservoir;
        const Scalar base = blk.t - tr.blocks[ix].t;
        if (base >= kMinBase) {
            const Scalar s = (blk.c - tr.blocks[ix].c) / base;
            // Push (s, base) into the bounded slope reservoir (evict oldest when full). FIFO
            // eviction of the OLDEST short-baseline slope is safe: the baseline-weighted median
            // (score()) down-weights small-baseline samples, so an aged-out 120-125 s slope (e.g.
            // the dense pool a post-denial restart produces) carries little weight anyway, while
            // the most recent long-baseline slopes -- the current drift regime -- dominate.
            if (tr.slope_count < kSlopeReservoir) {
                const int w = (tr.slope_head + tr.slope_count) % kSlopeReservoir;
                tr.slope_v[w] = s;
                tr.slope_w[w] = base;
                ++tr.slope_count;
            } else {
                tr.slope_v[tr.slope_head] = s;
                tr.slope_w[tr.slope_head] = base;
                tr.slope_head = (tr.slope_head + 1) % kSlopeReservoir;
            }
        }
    }
    // Append the new block to the segment's block ring (evict oldest when full). block_count is
    // the water level: it increments until it saturates at kSlopeReservoir, after which the ring
    // wraps and the OLDEST block is evicted FIFO (block_head advances). The pairing loop above
    // always iterates exactly the live block_count entries, so saturation never double-counts.
    if (tr.block_count < kSlopeReservoir) {
        const int w = (tr.block_head + tr.block_count) % kSlopeReservoir;
        tr.blocks[w] = blk;
        ++tr.block_count;
    } else {
        tr.blocks[tr.block_head] = blk;
        tr.block_head = (tr.block_head + 1) % kSlopeReservoir;
    }
    tr.cur_t0 = Scalar(0);
    tr.cur_n  = 0;
}

void HeadingMonitor::segment_break() {
    // Close the open block of each track, then clear the segment's block history so no
    // cross-segment slope (over the GPS-denied / multipath staircase step) can ever form.
    for (int i = 0; i < n_; ++i) {
        finalize_block(track_[i]);
        track_[i].block_count = 0;
        track_[i].block_head  = 0;
    }
}

void HeadingMonitor::push_pair(Scalar t_pair, Scalar dt_pair, const Scalar* r) {
    // Cross-source consensus: m = median_i(r_i) removes the common-mode course error.
    if (n_ < 3) return;                          // inert below the median floor
    int idx[kMaxSources];
    const Scalar m = median_of(r, n_, idx);
    t_valid_ += dt_pair;
    for (int i = 0; i < n_; ++i) {
        Track& tr = track_[i];
        const Scalar dev    = r[i] - m;
        const Scalar excess = std::max(std::abs(dev) - kEvent, Scalar(0));
        if (excess > Scalar(0)) tr.event_sum += excess;   // event channel
        const Scalar rc = m + std::max(-kEvent, std::min(kEvent, dev));   // rate channel
        tr.cum += rc;
        // Block bookkeeping: a new sample closes the open block once it spans block_s.
        if (tr.cur_n > 0 && (t_pair - tr.cur_t0) >= kBlockS) finalize_block(tr);
        if (tr.cur_n == 0) tr.cur_t0 = t_pair;
        // Buffer the (t, cum) sample for the per-block MEDIAN; keep the first kBlockSamples and
        // drop any surplus (reservoir stance -- the median of the retained head stays robust).
        if (tr.cur_n < kBlockSamples) {
            tr.cur_ts[tr.cur_n] = t_pair;
            tr.cur_cs[tr.cur_n] = tr.cum;
        }
        ++tr.cur_n;
    }
    ++pairs_;
}

void HeadingMonitor::submit_fix(Scalar t_s, const Vec3& pos,
                                const Scalar* yaw, const Scalar* fwd) {
    if (n_ <= 0) return;
    if (!finite(t_s) || !finite(pos.x()) || !finite(pos.y())) return;
    for (int i = 0; i < n_; ++i) {
        if (!finite(yaw[i]) || !finite(fwd[i])) return;
    }

    if (!have_fix_) {                            // first fix: just seed the previous sample
        have_fix_ = true;
        fix_t_    = t_s;
        fix_pos_  = pos;
        for (int i = 0; i < n_; ++i) { fix_yaw_[i] = yaw[i]; fix_fwd_[i] = fwd[i]; }
        return;
    }

    const Scalar dt = t_s - fix_t_;
    // Course-sample (anchor) gating — straight-at-speed only (1.1). A failed gate still
    // advances the previous-fix sample (the chord telescopes over the rejected interval; the
    // course delta is taken anchor-to-anchor, so the skipped span never injects an offset).
    bool ok = (dt > kDtMin) && (dt <= kDtMax);
    const Scalar dE = pos.x() - fix_pos_.x();
    const Scalar dN = pos.y() - fix_pos_.y();
    const Scalar v_gps = std::hypot(dE, dN) / std::max(dt, Scalar(1e-9));
    if (ok) ok = (v_gps >= kVMin) && (v_gps <= kVMax);
    // Per-source forward + yaw increments over the fix interval -> consensus v_odo / yaw rate.
    // These are unsigned magnitudes feeding two SEPARATE gates: the cross-val gate (below)
    // compares GPS speed v_gps vs the odometry SPEED v_odo (|magnitude| difference, kills
    // multipath jumps); the forward-only gate is separate -- it tests the SIGN of the median
    // forward increment dx_med (a reverse flips the GPS course 180 deg).
    Scalar v_odo = Scalar(0), yaw_rate = Scalar(0), dx_med = Scalar(0);
    if (ok) {
        for (int i = 0; i < n_; ++i) {
            scratch_a_[i] = fwd[i] - fix_fwd_[i];                       // forward increment
            scratch_b_[i] = std::abs(yaw[i] - fix_yaw_[i]) / std::max(dt, Scalar(1e-9));
        }
        int idx[kMaxSources];
        dx_med   = median_of(scratch_a_, n_, idx);
        v_odo    = std::abs(dx_med) / std::max(dt, Scalar(1e-9));
        yaw_rate = median_of(scratch_b_, n_, idx);
        const Scalar xtol = std::max(kXvalAbs, kXvalRel * v_gps);
        if (std::abs(v_gps - v_odo) >= xtol) ok = false;               // multipath cross-val
        if (dx_med <= Scalar(0))             ok = false;               // net forward motion
        if (yaw_rate >= kOmegaMax)           ok = false;               // |yaw rate| gate
    }

    if (ok) {
        const Scalar tm     = Scalar(0.5) * (fix_t_ + t_s);            // midpoint stamp
        const Scalar course = std::atan2(dN, dE);
        Scalar ym[kMaxSources];
        for (int i = 0; i < n_; ++i) ym[i] = Scalar(0.5) * (fix_yaw_[i] + yaw[i]);
        ++anchors_;
        if (have_anchor_) {
            const Scalar gap = tm - anc_t_;
            if (gap > kPairGap) {
                segment_break();                                       // GPS-contiguity break
            } else {
                const Scalar dc = course - anc_course_;
                Scalar r[kMaxSources];
                for (int i = 0; i < n_; ++i)
                    r[i] = wrap_pi(dc - (ym[i] - anc_yaw_[i]));
                push_pair(tm, gap, r);
            }
        }
        have_anchor_ = true;
        anc_t_       = tm;
        anc_course_  = course;
        for (int i = 0; i < n_; ++i) anc_yaw_[i] = ym[i];
    }

    // Advance the previous-fix sample (always — even on a rejected anchor).
    fix_t_   = t_s;
    fix_pos_ = pos;
    for (int i = 0; i < n_; ++i) { fix_yaw_[i] = yaw[i]; fix_fwd_[i] = fwd[i]; }
}

Scalar HeadingMonitor::score(int i) const {
    if (i < 0 || i >= n_) return Scalar(-1);
    const Track& tr = track_[i];
    if (tr.slope_count < 1) return Scalar(-1);                         // no baseline yet
    // Copy the slope ring into scratch (the reservoir may be a wrapped ring).
    Scalar v[kMaxSources > kSlopeReservoir ? kMaxSources : kSlopeReservoir];
    Scalar w[kMaxSources > kSlopeReservoir ? kMaxSources : kSlopeReservoir];
    for (int k = 0; k < tr.slope_count; ++k) {
        const int ix = (tr.slope_head + k) % kSlopeReservoir;
        v[k] = tr.slope_v[ix];
        w[k] = tr.slope_w[ix];
    }
    int idx[kMaxSources > kSlopeReservoir ? kMaxSources : kSlopeReservoir];
    const Scalar slope = weighted_median(v, w, tr.slope_count, idx);
    const Scalar event_rate = (t_valid_ > Scalar(0)) ? tr.event_sum / t_valid_ : Scalar(0);
    return std::abs(slope) + event_rate;
}

bool HeadingMonitor::scored(int i) const {
    if (i < 0 || i >= n_) return false;
    return track_[i].slope_count >= 1;
}

int HeadingMonitor::scored_count() const {
    int c = 0;
    for (int i = 0; i < n_; ++i) if (track_[i].slope_count >= 1) ++c;
    return c;
}

Scalar HeadingMonitor::boost(int i, Scalar boost_max) const {
    const Scalar cap = std::max(Scalar(1), boost_max);
    if (i < 0 || i >= n_) return Scalar(1);
    // ABSTAIN until >= 2 sources are scored (insufficient data -> all boosts 1.0).
    int sc = 0;
    Scalar best = Scalar(0);
    bool   have_best = false;
    for (int j = 0; j < n_; ++j) {
        if (!scored(j)) continue;
        ++sc;
        // Floor each score at kScoreFloor (the GPS-course noise floor): sources at or below it
        // are mutually indistinguishable, so they all rank equal and the ratio stays continuous
        // (no zero-drift discontinuity). Floor BEFORE taking the min so `best` is also floored.
        const Scalar s = std::max(score(j), kScoreFloor);
        if (!have_best || s < best) { best = s; have_best = true; }
    }
    if (sc < 2 || !scored(i)) return Scalar(1);                        // abstain / unscored
    // Floor the denominator too: an exactly-zero (or sub-floor) drift no longer jumps to the
    // full cap -- it floors to kScoreFloor like every other below-resolvability source, so the
    // boost is continuous (boost = clip(cap * best / max(score_i, floor), 1, cap)).
    const Scalar s_i = std::max(score(i), kScoreFloor);
    const Scalar b   = cap * (best / s_i);
    return std::min(cap, std::max(Scalar(1), b));
}

}  // namespace ofc
