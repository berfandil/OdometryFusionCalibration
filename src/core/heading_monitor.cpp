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
    Block blk;
    blk.t = tr.cur_tsum / static_cast<Scalar>(tr.cur_n);
    blk.c = tr.cur_csum / static_cast<Scalar>(tr.cur_n);
    // Pair this block with every earlier block of the segment whose baseline >= min_base.
    for (int q = 0; q < tr.block_count; ++q) {
        const int ix = (tr.block_head + q) % kSlopeReservoir;
        const Scalar base = blk.t - tr.blocks[ix].t;
        if (base >= kMinBase) {
            const Scalar s = (blk.c - tr.blocks[ix].c) / base;
            // Push (s, base) into the bounded slope reservoir (evict oldest when full).
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
    // Append the new block to the segment's block ring (evict oldest when full).
    if (tr.block_count < kSlopeReservoir) {
        const int w = (tr.block_head + tr.block_count) % kSlopeReservoir;
        tr.blocks[w] = blk;
        ++tr.block_count;
    } else {
        tr.blocks[tr.block_head] = blk;
        tr.block_head = (tr.block_head + 1) % kSlopeReservoir;
    }
    tr.cur_t0   = Scalar(0);
    tr.cur_tsum = Scalar(0);
    tr.cur_csum = Scalar(0);
    tr.cur_n    = 0;
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
        tr.cur_tsum += t_pair;
        tr.cur_csum += tr.cum;
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
        const Scalar s = score(j);
        if (!have_best || s < best) { best = s; have_best = true; }
    }
    if (sc < 2 || !scored(i)) return Scalar(1);                        // abstain / unscored
    const Scalar s_i = score(i);
    if (!(s_i > Scalar(0))) return cap;                               // zero-drift -> full cap
    const Scalar b = cap * (best / s_i);
    return std::min(cap, std::max(Scalar(1), b));
}

}  // namespace ofc
