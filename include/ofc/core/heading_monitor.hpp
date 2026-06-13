// ofc/core/heading_monitor.hpp — GPS-course heading-drift monitor (Slice 19c, split
// policy layer c).
//
// Ranks the per-source heading quality against GPS course-over-ground ONLINE and emits a
// per-source ROTATION-channel weight BOOST that auto-discovers the heading-grade source
// (e.g. the FOG), replacing the hand-tuned rot_weight_prior=10 urban recipe. The weight
// rule needs only a RANKING (relative drift rates), NOT an absolute per-source yaw bias —
// which sidesteps the Slice-18/D28 observability boundary (position-only GPS cannot observe
// per-source yaw bias absolutely; source-vs-course drift RATES are directly comparable).
// Design + acceptance: SLICE19C_HEADING_MONITOR.md; validated reference math:
// tools/proto_heading_monitor.py.
//
// PIPELINE (all online: running accumulators, one previous-fix sample, one previous-anchor
// sample, one open block + a bounded slope reservoir per source):
//   1. Each APPLIED position fix is submitted with (t, odom-frame position, per-source
//      cumulative yaw + forward distance). The chord between CONSECUTIVE fixes yields a
//      course sample; ANCHOR gates (1.1) trust it only when the chord is a clean heading
//      proxy (straight-at-speed). Course is compared via DELTAS so any fixed mounting/frame
//      yaw offset cancels (no odom_from_enu_yaw needed).
//   2. Consecutive ANCHORS within `pairgap` form a residual pair (anchors deliberately
//      BRIDGE turns — wheel-yaw drift is turn-correlated; nothing is hard-dropped after
//      anchoring, or the cumulative residual's telescoping breaks and a permanent offset is
//      injected). Per source: r_i = wrap(dcourse - dyaw_i).
//   3. Cross-source split (1.2): m = median_i(r_i) removes common-mode course error; the
//      per-source deviation dev_i = r_i - m feeds a RATE channel (m + clamp(dev_i, +-event),
//      accumulated cumulatively) and an EVENT channel (the clamp excess — catches slip steps
//      a magnitude gate would reject). Requires n >= 3 sources (median degenerates below;
//      the monitor stays inert otherwise).
//   4. Drift estimator (1.3): within GPS-CONTIGUOUS segments only (a fix gap > pairgap closes
//      the segment), block medians of (t, cumulative residual) over `block_s`; the drift
//      slope is the baseline-weighted median of pairwise block-median slopes with baseline
//      >= `min_base`. score_i = |slope_i| + event_rate_i (rad/s; reported deg/h-agnostic — a
//      pure RANKING quantity).
//   5. Weight rule (1.4): boost_i = clip(boost_max * min_j score_j / score_i, 1, boost_max).
//      ABSTAIN (every boost exactly 1.0) until >= 2 sources are scored.
//
// STRICT CORE: fixed-capacity per-source state, a fixed-size slope reservoir, bounded loops,
// no heap after construction, no exceptions, status-free value readouts (double math).
#ifndef OFC_CORE_HEADING_MONITOR_HPP
#define OFC_CORE_HEADING_MONITOR_HPP

#include "ofc/core/types.hpp"

namespace ofc {

// Fixed gate/estimator constants (Slice-19 veto-constant style: swept in the prototype,
// pinned here — only the boost cap is a Config knob). Units SI (seconds, m/s, rad).
namespace heading_monitor_const {
constexpr int    kMaxSources = 32;                 // matches Result arrays / Config cap
// A bounded per-source slope reservoir (1.3). Each finalized block pairs with EVERY earlier
// block of the current segment whose baseline >= min_base, so the slope count is O(blocks^2)
// per long segment; a fixed reservoir of this depth bounds the WCET. When full, the OLDEST
// slope sample is evicted (a sliding Theil-Sen pool): the most recent blocks — the current
// drift regime — always dominate the baseline-weighted median, and an old multipath block is
// aged out rather than pinned forever. 256 covers a multi-hour drive's contiguous segment at
// the 60-s block cadence with margin (the prototype's longest KAIST segment pooled < 200).
constexpr int    kSlopeReservoir = 256;

constexpr Scalar kDtMin    = Scalar(0.05);         // fix spacing > 0.05 s
constexpr Scalar kDtMax    = Scalar(3.0);          // fix spacing <= 3 s
constexpr Scalar kVMin     = Scalar(3.0);          // horizontal GPS speed >= 3 m/s
constexpr Scalar kVMax     = Scalar(30.0);         // <= 30 m/s
constexpr Scalar kXvalAbs  = Scalar(3.0);          // |v_gps - v_odo| < max(3, 0.5 v)
constexpr Scalar kXvalRel  = Scalar(0.5);
constexpr Scalar kOmegaMax = Scalar(0.05235987755982988);  // 3 deg/s, |yaw rate| gate
constexpr Scalar kPairGap  = Scalar(60.0);         // anchor-pair + segment-close gap (s)
constexpr Scalar kEvent    = Scalar(0.08726646259971647);  // 5 deg, rate/event clamp split
constexpr Scalar kBlockS   = Scalar(60.0);         // block length (s)
constexpr Scalar kMinBase  = Scalar(120.0);        // min slope baseline (s)
}  // namespace heading_monitor_const

// Online per-source GPS-course heading-drift monitor. One instance owned by the estimator
// Impl (allocated once); fed one sample per APPLIED position fix; read for the per-source
// boost each step.
class HeadingMonitor {
public:
    static constexpr int kMaxSources = heading_monitor_const::kMaxSources;

    HeadingMonitor() { reset(); }

    // Bind the active source COUNT (the registered participants whose yaw the monitor ranks)
    // and clear all state. n is clamped to [0, kMaxSources]. The monitor is INERT (every boost
    // 1.0) while n < 3 (the median consensus degenerates below three sources).
    void configure(int n);

    // Drop all accumulated state (keeps the source count). Called from configure() and the
    // estimator's init()/reset paths.
    void reset();

    // Submit ONE applied position fix. `t_s` is the frontier time in seconds; `pos` the fix's
    // odom-frame position (x,y the horizontal plane); `yaw` / `fwd` the per-source CUMULATIVE
    // base-frame yaw (rad, unwrapped) and forward distance (m) at this fix, indexed by the
    // SAME source order configure() counted. The chord to the PREVIOUS fix is gated into a
    // course anchor (1.1); consecutive anchors within pairgap feed the residual/slope
    // machinery. Non-finite inputs are ignored (no NaN ever enters the accumulators). Bounded
    // loops, no heap.
    void submit_fix(Scalar t_s, const Vec3& pos, const Scalar* yaw, const Scalar* fwd);

    // Per-source drift score (rad/s; |slope| + event-rate). Returns < 0 (sentinel) for a
    // source not yet scored (no slope baseline formed). A LOWER score is a BETTER heading
    // source.
    Scalar score(int i) const;

    // Whether source `i` has a finalized drift score (a slope baseline has formed in some
    // GPS-contiguous segment). False for an out-of-range index.
    bool scored(int i) const;

    // Number of sources currently scored.
    int scored_count() const;

    // Per-source ROTATION-channel boost for the given cap (>= 1). The weight rule
    // boost_i = clip(boost_max * min_j score_j / score_i, 1, boost_max) over the SCORED
    // sources; ABSTAIN (returns exactly 1.0 for every source) until >= 2 are scored. An
    // unscored source always returns 1.0 (it cannot be ranked yet). boost_max < 1 is clamped
    // to 1 (a degenerate cap disables boosting). No heap; bounded.
    Scalar boost(int i, Scalar boost_max) const;

    // Anchor / pair counters (diagnostics + tests).
    int anchor_count() const { return anchors_; }
    int pair_count()   const { return pairs_; }

private:
    // A finalized block of the CURRENT segment: median time + median cumulative residual.
    struct Block { Scalar t = 0; Scalar c = 0; };

    // Per-source drift-estimator state (mirrors the prototype's SlopeTracker, fixed-capacity).
    struct Track {
        // Cumulative rate-channel residual (the staircase whose SLOPE is the drift rate).
        Scalar cum = 0;
        // Open block accumulator: running sums for the median-free block reduction. The
        // prototype takes the per-block MEDIAN of (t, cum); for the strict core we use the
        // block MEAN (sums / count) — the block is a short (block_s) near-linear run of the
        // staircase, so mean and median coincide to estimator tolerance, and the mean needs no
        // per-block heap/sort. cur_t0 marks the open block's first sample time.
        Scalar cur_t0   = 0;
        Scalar cur_tsum = 0;
        Scalar cur_csum = 0;
        int    cur_n    = 0;
        // Finalized blocks of the current segment (bounded; pairwise-slope source).
        Block  blocks[heading_monitor_const::kSlopeReservoir];
        int    block_count = 0;          // saturates at the reservoir depth
        int    block_head  = 0;          // ring start (oldest) when saturated
        // Slope reservoir (baseline-weighted-median pool): value + baseline weight.
        Scalar slope_v[heading_monitor_const::kSlopeReservoir];
        Scalar slope_w[heading_monitor_const::kSlopeReservoir];
        int    slope_count = 0;          // filled slots (<= reservoir)
        int    slope_head  = 0;          // ring start (oldest) when saturated
        // Event channel (clamp-excess accumulator) + valid pair time (event-rate denominator).
        Scalar event_sum = 0;            // sum of |dev| - event excess (rad)
    };

    void finalize_block(Track& tr);      // close the open block, pool its slopes
    void segment_break();                // a pairgap exceeded: finalize + clear blocks
    void push_pair(Scalar t_pair, Scalar dt_pair, const Scalar* r);  // one residual pair

    int     n_ = 0;                      // active source count
    bool    have_fix_ = false;           // a previous fix sample is buffered
    bool    have_anchor_ = false;        // a previous anchor is buffered

    // Previous fix sample.
    Scalar  fix_t_ = 0;
    Vec3    fix_pos_ = Vec3::Zero();
    Scalar  fix_yaw_[kMaxSources] = {};
    Scalar  fix_fwd_[kMaxSources] = {};

    // Previous anchor sample (midpoint-stamped).
    Scalar  anc_t_ = 0;
    Scalar  anc_course_ = 0;
    Scalar  anc_yaw_[kMaxSources] = {};

    Track   track_[kMaxSources];
    Scalar  t_valid_ = 0;                // total valid pair time (event-rate denominator)

    int     anchors_ = 0;
    int     pairs_   = 0;

    // Scratch for the per-pair cross-source median + the boost/slope weighted medians (no
    // heap; sized at the cap). Reused sequentially.
    mutable Scalar scratch_a_[kMaxSources];
    mutable Scalar scratch_b_[kMaxSources];
};

}  // namespace ofc
#endif  // OFC_CORE_HEADING_MONITOR_HPP
