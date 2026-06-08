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
    bool     bias_states     = false;   // GPS/INS-style drift removal (D22; Slice 11b)
    // Random-walk process-noise rate for this source's body-twist bias (Slice 11b, Option A).
    // Only used when bias_states == true AND this is the single driving source: the augmented
    // ESKF adds Q_bias = bias_process_noise * dt * I6 to the bias block each predict, letting an
    // absolute-ref update slowly re-estimate the bias. Larger => the bias tracks faster but is
    // noisier; >= 0 (validate enforces). Default small (slow constant-bias assumption).
    Scalar   bias_process_noise = 1e-4;
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
    // Lifecycle (Slice 3): below this many PARTICIPATING sources fused in a step,
    // outlier rejection degrades and the phase stays DEGRADED (a health warning);
    // at/above it the phase is NOMINAL (CONFIG §3). Lower-bound only — NOT tied to
    // max_sources: if min_sources_warn > max_sources, NOMINAL is simply never
    // reached (a legitimate "never enough sources" state).
    int        min_sources_warn = 3;

    // Median fusion
    int        weiszfeld_max_iters = 10;
    Scalar     weiszfeld_tol  = 1e-6;
    Scalar     weiszfeld_eps  = 1e-9;
    Scalar     metric_lambda  = 1.0;
    bool       adaptive_q     = true;
    // Multiplier on the spread-derived adaptive Q (CONFIG §3). RE-CALIBRATED against the TRUE
    // interior robust median (D3 median fix) with the /n_eff reduction REMOVED (D4). The fixed
    // median's spread is correctly sized (distance-to-centroid of the participating sources), so
    // the calibration now targets genuine CONSISTENCY (NEES ~ DOF=6 from below) rather than the
    // old pinning-median's pessimism. Chosen SAFETY-FIRST against the sim NEES harness across 4
    // trajectories (nees_traj, mixed, turning, long straight) x 2 noise levels, M=30 seeds:
    //   q_scale  worst-case ensemble-mean pose NEES (1x / 2x noise)
    //     0.5      6.77 / 5.44   <- OVERCONFIDENT (>6) at 1x: REJECTED
    //     0.7      4.85 / 3.89   <- chosen: closest to DOF=6 from below, never overconfident
    //     1.0      3.40 / 2.73
    //     1.5      2.27 / 1.82
    //     2.0      1.71 / 1.36
    // At 0.7 the worst-case is ~4.85 (1x) — near-consistent (NEES approaching DOF=6) yet NEVER
    // overconfident on any trajectory at either noise level, keeping a conservative margin below
    // 6 (the sim under-states real-world model mismatch, so we sit a touch below DOF, not on it).
    Scalar     q_scale        = 0.7;
    // Per-axis minimum (additive) Q on the pose increment ([trans;rot]). Small by default: the
    // no-ref consistency objective wants the spread term to dominate (Slice-14 sweep showed a
    // larger floor only ADDS no-ref pessimism without helping multi-source rigs). DEPLOYMENT
    // GUIDANCE (CONFIG §3): for absolute-ref / GPS rigs whose sources have ~zero spread (identical
    // odometry), RAISE the translation floor (e.g. q_floor[0..2] ~ 1e-3) so the position
    // covariance grows between fixes and the absolute-fix Kalman gain can actually pull — the
    // adapters/GPS drift tests do exactly this test-locally. Left small here as the no-ref default.
    Scalar     q_floor[6]     = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
    // Motion-proportional (distance-aware) dead-reckoning process noise (CONFIG §3). The spread-
    // derived adaptive Q cannot see COMMON-MODE drift (sources that agree yet drift together), so
    // on real data it is overconfident; these grow the per-window pose process noise with the
    // distance/rotation moved. Eskf::add_distance_q adds (q_dist_trans*||delta.t||)^2 to each
    // translation diagonal and (q_dist_rot*||log(delta.R)||)^2 to each rotation diagonal. Default
    // 0 = OFF (byte-identical to the spread-only behavior); a real deployment sets them.
    Scalar     q_dist_trans   = 0.0;   // std growth per metre of per-window translation (>= 0)
    Scalar     q_dist_rot     = 0.0;   // std growth per radian of per-window rotation (>= 0)

    // Weights
    Scalar      reliability_ema_alpha = 0.02;
    // Variance-EMA reliability clamp (Slice 9, D17). The reliability multiplier is
    // ref_var / max(resid_var, eps) clamped to [reliability_floor, reliability_cap]:
    //   * reliability_floor in (0, 1] — a NOISY source can never collapse to 0 (always
    //     recoverable); a multiplier of 1 disables down-weighting.
    //   * reliability_cap   in [1, ∞) — a quiet source never DOMINATES the median.
    Scalar      reliability_floor = 0.2;   // min reliability multiplier (0, 1]
    Scalar      reliability_cap   = 5.0;   // max reliability multiplier [1, ∞)
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
    // SCALE is a strictly-positive RATIO (bn/ref_mag) whose calibrated value clusters around a
    // UNIT ratio (1.0) for a well-mounted source. The generic HistogramConfig default range of
    // [-1, 1] is WRONG for a scale: 1.0 lands exactly on the (half-open) upper boundary, clamps
    // into the last bin, and — being a boundary bin with no right neighbor — the parabolic
    // sub-bin refine in mode() falls back to that bin's CENTER (~0.984 at 64 bins), so a true
    // unit residual COMMITS ~0.984 instead of 1.0 and the rising-edge feedback folds the error
    // into prior_scale. Default this field to a positive-ratio range with 1.0 STRICTLY INTERIOR
    // ([0.5, 1.5] — the convention every dedicated calibration test already uses), keeping the
    // generic bin-count / aging defaults. The other (signed) histograms keep [-1, 1].
    HistogramConfig scale_hist{ /*bins=*/64, /*range_min=*/0.5, /*range_max=*/1.5 };
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
