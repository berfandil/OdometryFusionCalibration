// ofc/core/timesync.cpp — per-source clock-offset estimation (Slice 5, D16).
//
// STRICT CORE: all storage is fixed-capacity member arrays bound in configure();
// push()/update() allocate nothing and run loops bounded by the buffered sample
// count, the lag-scan radius, or kMaxSources. No exceptions; status-code returns; double.
//
// ALGORITHM (mirrors the header + the PLAN):
//   * Common grid. Samples are keyed to integer grid indices g = round(stamp / dt_ns)
//     (dt = 1/tick_rate_hz). Each source keeps a monotone ring of the most recent
//     kMaxSamples grid slots; a skipped grid slot is filled by carrying the previous
//     value forward (place()), so a missed push never punches a hole in the signal.
//   * Cross-correlation + SIGN. The source observes the base signal ADVANCED by the
//     clock offset: src(t) = ref(t + off)  (canonical sign — positive off => source
//     clock ahead of base; DECISIONS D21, CONFIG §9). Equivalently ref[g] = src[g - p]
//     with p = off/dt. We therefore score the alignment of ref[i] against src[i - lag]
//     over the overlap and MAXIMIZE over integer lag ∈ [-max_lag, +max_lag]; the peak
//     sits at lag = p, so the recovered offset = lag * dt = off DIRECTLY (positive
//     planted -> positive recovered). offset(id) returns the histogram mode of that.
//   * Pluggable metric (all normalized to "maximize a score"):
//       L1    : score = -Σ |ref - src_shifted|
//       L2    : score = -Σ (ref - src_shifted)^2
//       Ratio : score = -Σ |ref - src_shifted| / (|ref| + |src_shifted| + eps)
//                       (a bounded, scale-robust normalized abs-difference cost)
//       NCC   : score =  Σ (ref-ref̄)(src_shifted-src̄) / (n * σ_ref * σ_src)
//                       (zero-mean normalized cross-correlation in [-1, 1])
//   * Parabolic sub-sample refine. Around the best integer lag L*, fit a parabola to
//     the three scores s(L*-1), s(L*), s(L*+1):
//       δ = 0.5 (s_- - s_+) / (s_- - 2 s_0 + s_+)   (guarded denom, clamped to ±0.5),
//     refined lag = L* + δ, offset = (L* + δ) * dt. (Same formula as the histogram
//     peak refine — the score is concave near its max for all four metrics.)
//   * Excitation gate. Skip the window (no vote) unless BOTH the reference and the
//     source overlap windows have ‖ω‖ variance > excitation_min_var. A straight (no-
//     rotation) trajectory has ~zero ‖ω‖ variance and is rejected -> no confident
//     estimate (the time-offset DOF is observable only when ‖ω‖ VARIES, DESIGN §6).
//   * Vote. The accepted offset (seconds) is voted into the source's Histogram1D
//     (config offset_hist); offset() = mode, confidence() = concentration. Constant
//     offset is the histogram peak; slow drift is tracked by the histogram aging.
#include "ofc/core/timesync.hpp"

#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
constexpr Scalar    kNanosPerSec = Scalar(1e9);
constexpr Scalar    kRatioEps    = Scalar(1e-9);
constexpr Scalar    kParabEps    = Scalar(1e-15);
// A score sentinel for an empty overlap (smaller than any real score). Metrics return
// negated non-negative costs (>= very-negative) or an NCC in [-1,1], so this is safely
// below every valid score.
constexpr Scalar    kNoScore     = -Scalar(1e300);
} // namespace

Status TimeSync::configure(const Config& cfg, SourceId reference_id) {
    if (!(cfg.tick_rate_hz > Scalar(0)))          return Status::OutOfRange;
    if (!(cfg.max_lag_s > Scalar(0)) || cfg.max_lag_s > Scalar(2))
        return Status::OutOfRange;
    if (static_cast<int>(reference_id) >= kMaxSources) return Status::OutOfRange;

    sample_dt_ = Scalar(1) / cfg.tick_rate_hz;
    dt_ns_     = static_cast<long long>(std::llround(sample_dt_ * kNanosPerSec));
    if (dt_ns_ < 1) dt_ns_ = 1;                  // guard a degenerate tiny period

    // Bounded lag scan radius (grid steps), capped so a meaningful overlap remains.
    long long lag = std::llround(cfg.max_lag_s * cfg.tick_rate_hz);
    if (lag < 1)               lag = 1;
    if (lag > kMaxSamples / 4) lag = kMaxSamples / 4;
    max_lag_samp_ = static_cast<int>(lag);

    metric_        = cfg.match_metric;
    excite_min_var_ = std::max(Scalar(0), cfg.excitation_min_var);
    reference_id_   = reference_id;
    hist_cfg_       = cfg.offset_hist;
    min_overlap_    = 4;

    // Configure all per-source histograms up front (one per channel slot). Propagate
    // the first failure so a bad offset_hist is surfaced at configure().
    for (int i = 0; i < kMaxSources; ++i) {
        const Status hs = hists_[i].configure(hist_cfg_);
        if (!ok(hs)) return hs;
    }

    configured_   = true;
    reset();

    // Pre-create the reference channel so a source-only update() still has a ref slot.
    ref_slot_ = ensure_slot(reference_id_);
    return Status::Ok;
}

void TimeSync::reset() {
    for (int i = 0; i < kMaxSources; ++i) {
        chans_[i] = Channel{};
        ids_[i]   = 0;
        hists_[i].reset();
    }
    source_count_ = 0;
    ref_slot_     = -1;
}

int TimeSync::slot_for(SourceId id) const {
    for (int i = 0; i < source_count_; ++i) {
        if (ids_[i] == id) return i;
    }
    return -1;
}

int TimeSync::ensure_slot(SourceId id) {
    const int s = slot_for(id);
    if (s >= 0) return s;
    if (source_count_ >= kMaxSources) return -1;
    const int slot = source_count_++;
    ids_[slot]          = id;
    chans_[slot]        = Channel{};
    chans_[slot].active = true;
    if (id == reference_id_) ref_slot_ = slot;
    return slot;
}

void TimeSync::place(Channel& c, long long g, Scalar value) {
    if (!c.has_base) {
        c.base_grid = g;
        c.has_base  = true;
        c.count     = 1;
        c.grid[0]   = value;
        return;
    }

    // Within the current window: overwrite (latest wins for a duplicate grid slot).
    if (g >= c.base_grid && g < c.base_grid + c.count) {
        c.grid[static_cast<int>(g - c.base_grid)] = value;
        return;
    }

    // Behind the window start: ignore (older than the retained history).
    if (g < c.base_grid) return;

    // Bounded-WCET guard: a gap larger than the whole ring would slide every old sample
    // out anyway. Restart the channel fresh at `g` instead of iterating the gap (keeps
    // the gap-fill loop below bounded by kMaxSamples regardless of how big the jump is).
    const long long window_end = c.base_grid + c.count;   // first empty grid index
    if (g - window_end >= static_cast<long long>(kMaxSamples)) {
        c.base_grid = g;
        c.count     = 1;
        c.grid[0]   = value;
        return;
    }

    // Ahead of the window: extend, carrying the previous value forward over any gap.
    long long next = window_end;                     // next absolute grid index
    while (next < g) {
        const Scalar carry = c.grid[c.count - 1];    // last known value
        if (c.count < kMaxSamples) {
            c.grid[c.count++] = carry;
        } else {
            // Slide the window forward by one (drop the oldest), carry into the tail.
            for (int i = 1; i < c.count; ++i) c.grid[i - 1] = c.grid[i];
            c.grid[c.count - 1] = carry;
            ++c.base_grid;
        }
        ++next;
    }
    // Deposit the real sample at g.
    if (c.count < kMaxSamples) {
        c.grid[c.count++] = value;
    } else {
        for (int i = 1; i < c.count; ++i) c.grid[i - 1] = c.grid[i];
        c.grid[c.count - 1] = value;
        ++c.base_grid;
    }
}

void TimeSync::push(SourceId id, Timestamp stamp, Scalar omega_norm) {
    if (!configured_) return;
    if (!std::isfinite(omega_norm) || omega_norm < Scalar(0)) return;
    const int slot = ensure_slot(id);
    if (slot < 0) return;                            // at capacity
    const long long g = static_cast<long long>(std::llround(
        static_cast<Scalar>(stamp) / static_cast<Scalar>(dt_ns_)));
    place(chans_[slot], g, omega_norm);
}

Scalar TimeSync::sample_at(const Channel& c, long long g) const {
    // Clamp out-of-range indices to the nearest in-range sample (carry the edges).
    long long idx = g - c.base_grid;
    if (idx < 0)            idx = 0;
    if (idx >= c.count)     idx = c.count - 1;
    return c.grid[static_cast<int>(idx)];
}

Scalar TimeSync::window_variance(const Channel& c, long long g0, long long g1) const {
    if (g1 < g0) return Scalar(0);
    const long long n = g1 - g0 + 1;
    if (n <= 1) return Scalar(0);
    Scalar mean = Scalar(0);
    for (long long g = g0; g <= g1; ++g) mean += sample_at(c, g);
    mean /= static_cast<Scalar>(n);
    Scalar var = Scalar(0);
    for (long long g = g0; g <= g1; ++g) {
        const Scalar d = sample_at(c, g) - mean;
        var += d * d;
    }
    return var / static_cast<Scalar>(n);
}

Scalar TimeSync::match_score(const Channel& src, const Channel& ref, int lag) const {
    // Overlap of absolute grid indices i where BOTH ref[i] and src[i - lag] are real
    // samples (not edge-carried): i in [ref.base, ref.end) AND i-lag in [src.base, src.end).
    const long long ref_lo = ref.base_grid;
    const long long ref_hi = ref.base_grid + ref.count - 1;
    const long long src_lo = src.base_grid + lag;          // i such that i-lag == src.base
    const long long src_hi = src.base_grid + src.count - 1 + lag;
    const long long i0 = std::max(ref_lo, src_lo);
    const long long i1 = std::min(ref_hi, src_hi);
    const long long n  = i1 - i0 + 1;
    if (n < static_cast<long long>(min_overlap_)) return kNoScore;

    if (metric_ == MatchMetric::NCC) {
        // Zero-mean normalized cross-correlation over the overlap.
        Scalar ma = Scalar(0), mb = Scalar(0);
        for (long long i = i0; i <= i1; ++i) {
            ma += sample_at(ref, i);
            mb += sample_at(src, i - lag);
        }
        const Scalar inv_n = Scalar(1) / static_cast<Scalar>(n);
        ma *= inv_n; mb *= inv_n;
        Scalar num = Scalar(0), va = Scalar(0), vb = Scalar(0);
        for (long long i = i0; i <= i1; ++i) {
            const Scalar a = sample_at(ref, i) - ma;
            const Scalar b = sample_at(src, i - lag) - mb;
            num += a * b;
            va  += a * a;
            vb  += b * b;
        }
        const Scalar denom = std::sqrt(va * vb);
        if (denom < kRatioEps) return kNoScore;     // a flat window correlates with nothing
        return num / denom;                         // in [-1, 1]
    }

    Scalar cost = Scalar(0);
    for (long long i = i0; i <= i1; ++i) {
        const Scalar a = sample_at(ref, i);
        const Scalar b = sample_at(src, i - lag);
        const Scalar d = a - b;
        if (metric_ == MatchMetric::L1) {
            cost += std::abs(d);
        } else if (metric_ == MatchMetric::L2) {
            cost += d * d;
        } else { // Ratio — normalized abs difference, bounded per-sample in [0, 1].
            cost += std::abs(d) / (std::abs(a) + std::abs(b) + kRatioEps);
        }
    }
    return -cost;   // maximize the negated cost
}

bool TimeSync::estimate_offset(const Channel& src, const Channel& ref,
                               Scalar& out_off) const {
    if (src.count < min_overlap_ || ref.count < min_overlap_) return false;

    // Scan integer lags; keep the best (maximizing) score and its neighbours for the
    // parabolic refine. Bounded by [-max_lag, +max_lag].
    int    best_lag   = 0;
    Scalar best_score = kNoScore;
    for (int lag = -max_lag_samp_; lag <= max_lag_samp_; ++lag) {
        const Scalar s = match_score(src, ref, lag);
        if (s > best_score) {
            best_score = s;
            best_lag   = lag;
        }
    }
    if (best_score <= kNoScore) return false;       // no usable overlap at any lag

    // Excitation gate: require BOTH overlap windows (at the best lag) to vary in ‖ω‖.
    // Use the reference overlap span [i0, i1] and the matching source span [i0-lag, i1-lag].
    const long long ref_lo = ref.base_grid;
    const long long ref_hi = ref.base_grid + ref.count - 1;
    const long long src_lo = src.base_grid + best_lag;
    const long long src_hi = src.base_grid + src.count - 1 + best_lag;
    const long long i0 = std::max(ref_lo, src_lo);
    const long long i1 = std::min(ref_hi, src_hi);
    const Scalar var_ref = window_variance(ref, i0, i1);
    const Scalar var_src = window_variance(src, i0 - best_lag, i1 - best_lag);
    if (var_ref <= excite_min_var_ || var_src <= excite_min_var_) return false;

    // Parabolic sub-sample refine around best_lag (needs both integer neighbours inside
    // the scan range; at the scan boundary fall back to the integer lag).
    Scalar delta = Scalar(0);
    if (best_lag > -max_lag_samp_ && best_lag < max_lag_samp_) {
        const Scalar sm = match_score(src, ref, best_lag - 1);
        const Scalar s0 = best_score;
        const Scalar sp = match_score(src, ref, best_lag + 1);
        if (sm > kNoScore && sp > kNoScore) {
            const Scalar denom = sm - Scalar(2) * s0 + sp;
            if (std::abs(denom) > kParabEps) {
                delta = Scalar(0.5) * (sm - sp) / denom;
                if (delta >  Scalar(0.5)) delta =  Scalar(0.5);
                if (delta < -Scalar(0.5)) delta = -Scalar(0.5);
            }
        }
    }

    out_off = (static_cast<Scalar>(best_lag) + delta) * sample_dt_;
    return true;
}

void TimeSync::update() {
    if (!configured_ || ref_slot_ < 0) return;
    const Channel& ref = chans_[ref_slot_];
    if (ref.count < min_overlap_) return;

    for (int i = 0; i < source_count_; ++i) {
        if (i == ref_slot_) continue;
        const Channel& src = chans_[i];
        if (!src.active) continue;
        Scalar off = Scalar(0);
        if (estimate_offset(src, ref, off)) {
            hists_[i].add(off, Scalar(1));
        }
    }
}

Scalar TimeSync::offset(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0 || s == ref_slot_) return Scalar(0);
    if (hists_[s].empty()) return Scalar(0);        // unvoted -> "no offset"
    return hists_[s].mode();
}

Scalar TimeSync::confidence(SourceId id) const {
    if (!configured_) return Scalar(0);
    const int s = slot_for(id);
    if (s < 0 || s == ref_slot_) return Scalar(0);
    return hists_[s].confidence();
}

int TimeSync::sample_count(SourceId id) const {
    if (!configured_) return 0;
    const int s = slot_for(id);
    if (s < 0) return 0;
    return chans_[s].count;
}

} // namespace ofc
