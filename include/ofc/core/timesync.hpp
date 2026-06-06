// ofc/core/timesync.hpp — per-source clock-offset estimation (Slice 5, D16).
//
// Estimates each non-reference source's clock offset relative to the pinned
// reference by cross-correlating the EXTRINSIC-INVARIANT scalar signal ‖ω‖(t)
// (the body angular-rate magnitude). ‖ω‖ does not depend on the source's mount
// (a rotation magnitude is invariant under the sensor->base conjugation), so
// time-sync decouples from spatial calibration and runs FIRST (DESIGN §2/§3).
//
// SIGNAL ‖ω‖ EXTRACTION (documented contract). From a source's SE(3) delta over
// a small sub-interval [t-h, t] the instantaneous body angular-rate magnitude is
//      ‖ω‖(t) ≈ ‖log(ΔR)‖ / dt = ‖so3::log(delta.R)‖ / dt ,   dt = h seconds.
// The estimator (or a test) samples this per source at a COMMON rate and feeds it
// to push(). TimeSync itself never queries a source — it only consumes the scalar.
//
// SIGN CONVENTION (the load-bearing piece — DECISIONS D21, CONFIG §9). The planted
// convention is: a POSITIVE offset `off` means the source clock is AHEAD of base,
// so the source's reported window [t0,t1] maps to TRUE base motion over
// [t0+off, t1+off]. Indexed by base/query time t, the source therefore observes the
// base signal ADVANCED by off:
//      omega_src(t) = omega_base(t + off)   ( ≈ omega_ref(t + off) ).
// We search for the lag L that best aligns omega_src(t) to omega_ref(t + L); the
// match peaks at L = off. So the recovered estimate IS the offset with the canonical
// sign (a positive planted offset is recovered positive). offset(id) returns that
// value directly. The estimator REMOVES the offset by shifting a source's query
// interval EARLIER by offset (query [q0-off, t1-off]) so the source's internal +off
// re-shift lands back on true base time.
//
// STRICT CORE: per-source fixed-capacity ‖ω‖ ring buffers + a bounded lag scan +
// bounded loops; no heap after configure(); no exceptions; status-code returns;
// double math. Reuses the Slice-4 Histogram1D for the per-source offset histogram.
#ifndef OFC_CORE_TIMESYNC_HPP
#define OFC_CORE_TIMESYNC_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/histogram.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

namespace ofc {

// Commit-gate state machine for a calibrated DOF (DESIGN §2/§6). Returns the NEW
// committed state given the previous one and the current histogram reading:
//   * Not yet committed: commit only when the peak concentration clears
//     `commit_concentration` (τ_commit) AND the vote count clears `commit_min_votes`
//     (N_min) — so a burst of a few sharp votes cannot commit prematurely.
//   * Already committed: STAY committed (hysteresis) until the concentration falls
//     strictly below `commit_drop` (τ_drop), which re-opens the estimate. The N_min
//     gate does not re-apply once committed.
// Pure function — keeps no state of its own; the caller owns the per-source bool. The
// estimator holds the previous-committed flags in its fixed Impl arrays (strict core).
// Thresholds are validated by Config::validate (commit_drop < commit_concentration).
inline bool commit_gate(bool prev_committed, Scalar confidence, Scalar votes,
                        Scalar commit_concentration, int commit_min_votes,
                        Scalar commit_drop) {
    if (prev_committed) {
        return confidence >= commit_drop;                  // hysteresis: re-open below τ_drop
    }
    return confidence >= commit_concentration &&
           votes >= static_cast<Scalar>(commit_min_votes); // BOTH gates to first commit
}

// Cross-correlation time-sync over ‖ω‖(t). Configure once from a Config + the
// reference source id, push resampled ‖ω‖ samples per source, periodically call
// update() to run the xcorr + vote, then read offset()/confidence() per source.
class TimeSync {
public:
    // Compile-time maxima (strict core, no heap).
    static constexpr int kMaxSources = 32;     // matches Result arrays / Config cap
    static constexpr int kMaxSamples = 2048;   // per-source ‖ω‖ history depth

    TimeSync() = default;

    // Bind the active configuration and the pinned reference source id. Builds the
    // common-rate sample grid (sample period = 1 / tick_rate_hz), the bounded lag
    // scan radius (round(max_lag_s * tick_rate_hz), capped so a full reference window
    // remains), the excitation threshold (cfg.offset_hist drives the per-source
    // histogram), and clears all buffers. Validates:
    //   tick_rate_hz > 0                         -> else OutOfRange
    //   max_lag_s in (0, 2]                       -> else OutOfRange
    //   reference_id < kMaxSources                -> else OutOfRange
    //   offset_hist passes Histogram1D::configure -> propagated
    // No heap occurs here nor in any later call.
    Status configure(const Config& cfg, SourceId reference_id);

    // Drop all buffered samples + histogram votes (keeps the configuration).
    void reset();

    // Push one ‖ω‖ sample for source `id` at base/query time `stamp` (ns). Samples
    // are placed onto the common grid by their stamp; out-of-order or duplicate-grid
    // samples overwrite the matching slot (latest wins) so the per-source signal stays
    // a function of grid time, not of call order. A sample for an unknown id (not the
    // reference and never seen, beyond kMaxSources) is ignored. No-op if not configured
    // or omega_norm is non-finite/negative. No heap.
    void push(SourceId id, Timestamp stamp, Scalar omega_norm);

    // Run the xcorr + vote for every non-reference source that has enough overlapping,
    // EXCITED buffered signal against the reference. For each such source: cross-
    // correlate over lag ∈ [−max_lag, +max_lag] under the configured metric, parabolic
    // sub-sample refine the best lag, excitation-gate (skip unless BOTH the reference
    // and source overlap windows have ‖ω‖ variance > excitation_min_var), and vote the
    // accepted offset (seconds) into the source's histogram. Bounded loops only.
    void update();

    // Histogram MODE for source `id` (seconds), the committed offset estimate. Returns
    // 0 for the reference, an unknown id, or a source with no votes yet (the empty
    // histogram mode is its range midpoint; we special-case empty -> 0 so an unobserved
    // source reports "no offset" rather than the bin midpoint).
    Scalar offset(SourceId id) const;

    // Histogram CONCENTRATION for source `id` in [0,1] (0 for the reference / unknown /
    // unvoted). Drives the estimator's commit decision + diagnostics confidence.
    Scalar confidence(SourceId id) const;

    bool     configured() const { return configured_; }
    SourceId reference()  const { return reference_id_; }
    Scalar   sample_dt()  const { return sample_dt_; }   // common-grid period (s)
    int      max_lag_samples() const { return max_lag_samp_; }

    // Number of grid samples currently buffered for source `id` (0 if unknown).
    int sample_count(SourceId id) const;

    // Total vote mass currently in source `id`'s offset histogram (0 for the
    // reference / unknown / unvoted). Drives the estimator's N_min commit gate
    // (`commit_min_votes`): with unit-weight votes this is the live vote count
    // (under SlidingK it saturates at the window size). No heap; bounded.
    Scalar vote_count(SourceId id) const;

private:
    // Per-source resampled ‖ω‖ ring on the common grid. We key samples by an integer
    // grid index g = round((stamp - epoch) / dt_ns); the ring holds the most recent
    // kMaxSamples contiguous grid slots. `have_` marks whether the source is active.
    struct Channel {
        bool   active = false;
        // Ring storage: grid_[k] is the ‖ω‖ at absolute grid index base_grid_ + k for
        // k in [0, count_). A monotone-advancing window; gaps are filled by carrying
        // the previous sample forward (so a missed push does not punch a hole in the
        // correlation). count_ <= kMaxSamples.
        Scalar  grid[kMaxSamples] = {};
        long long base_grid = 0;     // absolute grid index of grid[0]
        int     count       = 0;     // number of filled slots
        bool    has_base    = false; // base_grid initialized
    };

    // Resolve the channel slot for an id (linear scan over the active set, bounded by
    // source_count_). Returns -1 if not present.
    int  slot_for(SourceId id) const;
    // Get-or-create a channel slot for an id (creates on first push, up to kMaxSources).
    // Returns -1 if at capacity.
    int  ensure_slot(SourceId id);

    // Append/overwrite the grid value at absolute grid index `g` in channel `c`.
    void place(Channel& c, long long g, Scalar value);

    // Cross-correlate channel `c` (source) against the reference channel and, if the
    // excitation gate passes, write the sub-sample offset (seconds) into `out_off` and
    // return true. Returns false (no vote) when there is insufficient overlap or the
    // window fails the excitation gate.
    bool estimate_offset(const Channel& src, const Channel& ref, Scalar& out_off) const;

    // Match SCORE (higher = better alignment) between the reference window ref[] and the
    // source window src[] when the source is shifted by integer `lag` grid steps, over
    // the inclusive overlap [i0, i1] of absolute grid indices. The metric is selected by
    // metric_. All four metrics are normalized to "maximize a score": cost metrics return
    // the NEGATED cost. Returns -infinity-ish (very negative) when the overlap is empty.
    Scalar match_score(const Channel& src, const Channel& ref, int lag) const;

    // Sample value of channel `c` at absolute grid index `g`, carrying the nearest
    // in-range sample for out-of-range indices (clamped). Precondition: c.count > 0.
    Scalar sample_at(const Channel& c, long long g) const;

    // Variance of channel `c` over absolute grid indices [g0, g1] (inclusive). Used by
    // the excitation gate. <= 0 when the span is empty.
    Scalar window_variance(const Channel& c, long long g0, long long g1) const;

    // Active configuration. configure() is the SOLE initializer of these fields and is
    // a precondition for every operation that reads them (push/update/offset/... all
    // short-circuit while !configured_), so the member initializers below are just safe
    // zero/neutral values for a default-constructed-but-unconfigured instance — they are
    // never the values used in a real computation (configure() overwrites them all).
    bool        configured_   = false;
    SourceId    reference_id_  = 0;
    MatchMetric metric_        = MatchMetric::L2;   // overwritten by cfg.match_metric
    Scalar      sample_dt_     = Scalar(0);         // set to 1/tick_rate_hz in configure()
    long long   dt_ns_         = 0;                 // set to round(sample_dt_*1e9)
    int         max_lag_samp_  = 0;                 // set from max_lag_s * tick_rate_hz
    Scalar      excite_min_var_ = Scalar(0);        // set from cfg.excitation_min_var
    int         min_overlap_   = 0;                 // set to 4 in configure()

    HistogramConfig hist_cfg_{};

    // Per-source channels + histograms. ref_slot_ indexes the reference channel.
    Channel      chans_[kMaxSources];
    Histogram1D  hists_[kMaxSources];
    SourceId     ids_[kMaxSources] = {};
    int          source_count_ = 0;
    int          ref_slot_     = -1;
};

} // namespace ofc
#endif // OFC_CORE_TIMESYNC_HPP
