// ofc/core/config.hpp — the single validated config struct (D19).
// File parsing lives in adapters; the core receives this struct and preallocates.
// Full knob documentation: CONFIG.md.
#ifndef OFC_CORE_CONFIG_HPP
#define OFC_CORE_CONFIG_HPP

#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

namespace ofc {

enum class ColdStart   { ReferenceOnly, MedianFromStart };
enum class LateSample  { Drop, Reintegrate };
enum class ConfCombine { NativeOnly, ModeledOnly, Sum, Max, Weighted };
enum class MatchMetric { L1, L2, Ratio, NCC };
enum class Phase2Strat { VsFusedBase, PairwisePinnedRef };
enum class Aging       { Decay, SlidingK };
enum class VoteWeight  { One, Confidence, Rotation, Combo };

struct HistogramConfig {
    int    bins       = 64;
    Scalar range_min  = -1.0;
    Scalar range_max  =  1.0;
    bool   circular   = false;
    Aging  aging      = Aging::Decay;
    Scalar decay_gamma = 0.999;
    int    sliding_k  = 1000;
    bool   vote_split = true;
    bool   subbin     = true;
};

struct SensorConfig {
    SourceId id              = 0;
    SE3      prior_extrinsic;            // also the Phase-1 histogram basepoint
    Scalar   prior_scale     = 1.0;
    Scalar   prior_time_offset_s = 0.0;
    Scalar   weight_prior    = 1.0;
    bool     native_confidence = true;
    Scalar   modeled_noise_per_m   = 0.01;
    Scalar   modeled_noise_per_rad = 0.01;
    bool     per_sensor_kf   = false;
    Scalar   kf_process_noise = 1.0;
    bool     scale_calib     = true;
    bool     bias_states     = false;   // GPS/INS-style drift removal (D22)
    bool     is_reference    = false;   // gauge anchor
};

struct Config {
    // Runtime / timing
    Scalar     tick_rate_hz   = 50.0;
    Scalar     fusion_delay_s = 0.05;
    Scalar     window_s       = 0.10;
    Scalar     calib_lag_s    = 0.20;
    LateSample late_policy    = LateSample::Drop;
    int        max_sources    = 8;       // [strict] sizes preallocation
    int        buffer_capacity = 512;    // [strict]

    // Gauge
    SourceId   reference_sensor_id = 0;
    ColdStart  cold_start     = ColdStart::ReferenceOnly;

    // Median fusion
    int        weiszfeld_max_iters = 10;
    Scalar     weiszfeld_tol  = 1e-6;
    Scalar     weiszfeld_eps  = 1e-9;
    Scalar     metric_lambda  = 1.0;
    bool       adaptive_q     = true;
    Scalar     q_scale        = 1.0;   // multiplier on the spread-derived Q (CONFIG §3)
    Scalar     q_floor[6]     = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};  // per-axis min Q ([trans;rot])

    // Weights
    Scalar      reliability_ema_alpha = 0.02;
    Scalar      weight_floor  = 0.05;
    Scalar      weight_cap    = 10.0;
    ConfCombine conf_combine  = ConfCombine::Sum;
    Scalar      confidence_blend = 0.5;   // blend factor when conf_combine == Weighted (CONFIG §4)

    // Time-sync
    bool         timesync_enabled = true;
    MatchMetric  match_metric = MatchMetric::L2;
    Scalar       max_lag_s    = 0.10;
    Scalar       excitation_min_var = 1e-4;   // min ‖ω‖ variance to accept a window (CONFIG §5)

    // Gates
    Scalar straight_omega_max = 0.05;
    Scalar straight_trans_min = 0.05;
    Scalar turn_omega_min     = 0.20;
    bool   reverse_fold       = true;
    bool   ref_cross_check    = false;

    // Calibration
    Phase2Strat phase2_strategy   = Phase2Strat::VsFusedBase;
    Scalar      commit_concentration = 0.6;
    // commit_min_votes is compared against the histogram TOTAL (vote_count), which is a vote
    // MASS, not a count, under any non-`One` vote_weight: with Confidence/Rotation/Combo each
    // vote's mass is scaled by the per-source Σ-confidence and/or ‖ω‖, so the same int gates
    // a wildly different EFFECTIVE vote count across sources. It is a true COUNT only when
    // vote_weight == One. RE-TUNE this threshold whenever vote_weight changes (the default
    // Combo can explode the mass on a high-confidence source). [Slice-8 review MINOR.]
    int         commit_min_votes  = 200;
    Scalar      commit_drop       = 0.4;
    VoteWeight  vote_weight       = VoteWeight::Combo;

    // Histograms (per quantity)
    HistogramConfig so3_hist;
    HistogramConfig roll_hist;
    HistogramConfig xyz_hist;
    HistogramConfig scale_hist;
    HistogramConfig offset_hist;

    // Absolute refs
    Scalar mahalanobis_chi2 = 9.0;

    // Output
    bool   emit_predicted_tip = true;
    bool   emit_diagnostics   = true;
    Scalar tip_cov_inflation  = 1.5;

    // Per-sensor (size <= max_sources). Storage owned by the caller/adapter.
    const SensorConfig* sensors      = nullptr;
    int                 sensor_count = 0;
};

// Bounds/capacity/consistency check. No throw; returns first violation.
Status validate(const Config& cfg);

} // namespace ofc
#endif // OFC_CORE_CONFIG_HPP
