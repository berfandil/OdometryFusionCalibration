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
    // Opt-in centroid sub-bin readout (Slice 16; honored only when subbin == true). The
    // parabolic interpolation systematically pulls toward the peak-bin center (a quarter-bin
    // offset reads ~0.1 bins instead of 0.25); the mass-weighted centroid over peak+-1 is
    // exact for a vote-split point mass at any sub-bin position. Default OFF = byte-identical
    // legacy behavior.
    bool   subbin_centroid = false;
};

// Per-DOF random-walk process-noise rates for a source's body-twist bias, in the codebase
// [v; omega] = [trans; rot] order (Slice 18 review/B2). The SCALAR form — one rate for all
// 6 DOFs — is accepted implicitly (the pre-B2 knob shape; every existing caller/config is
// unchanged). PINNING CONTRACT (multi_bias_enabled ONLY): a DOF whose rate is EXACTLY 0 is
// PINNED at zero — zero prior variance, zero walk, coupling column zero — excluded from
// estimation entirely. This lets a deployment enable e.g. ONLY the yaw-rate bias
// ("0 0 0 0 0 1e-5", the urban12 wheel −75°/h case) instead of opening 18 weakly-observable
// junk-sink DOFs. Option A (multi_bias_enabled off) keeps its original semantics — a 0 rate
// means "constant bias, still estimated under the prior" — and supports only the UNIFORM
// form (validate() rejects a non-uniform vector on an Option-A bias source). Each rate must
// be >= 0 (validate enforces per DOF).
struct BiasProcessNoise {
    Scalar dof[6];
    // Implicit scalar -> uniform (back-compat: `bias_process_noise = 1e-5` still compiles
    // and means the same thing it always did).
    BiasProcessNoise(Scalar uniform = Scalar(1e-4)) {
        for (int i = 0; i < 6; ++i) dof[i] = uniform;
    }
    Scalar& operator[](int i)       { return dof[i]; }
    Scalar  operator[](int i) const { return dof[i]; }
    bool uniform() const {
        for (int i = 1; i < 6; ++i) { if (dof[i] != dof[0]) return false; }
        return true;
    }
    // Back-compat scalar READ for legacy scalar contexts (comparisons in pre-B2 tests/
    // adapters): the DOF-0 rate, which IS the uniform value for every scalar-form config.
    // Per-DOF consumers must use operator[] — do not add scalar arithmetic on this type.
    operator Scalar() const { return dof[0]; }
};

struct SensorConfig {
    SourceId id              = 0;
    SE3      prior_extrinsic;            // also the Phase-1 histogram basepoint
    Scalar   prior_scale     = 1.0;
    Scalar   prior_time_offset_s = 0.0;
    Scalar   weight_prior    = 1.0;
    // ROTATION-channel weight multiplier (Slice 19, split-median policy layer (a) — config
    // priors). Only consumed when Config::split_median is on; multiplies this source's
    // effective weight in the ROTATION channel only:
    //   w_rot_i   = clamp(weight_prior_i * reliability_i * sigma_conf_i) * rot_weight_prior_i
    //   w_trans_i = clamp(weight_prior_i * reliability_i * sigma_conf_i)
    // (applied OUTSIDE the [weight_floor, weight_cap] clamp: the clamp bounds the noisy
    // Sigma-confidence placeholder, while this is a deliberate datasheet declaration — e.g.
    // a FOG-grade heading source carries rot_weight_prior ~ 10 so fusion's rotation channel
    // tracks it without distorting translation). >= 0 (validate; 0 = "never trust this
    // source's rotation"). INERT (byte-identical) when split_median is off. Hashed.
    Scalar   rot_weight_prior = 1.0;
    bool     native_confidence = true;
    Scalar   modeled_noise_per_m   = 0.01;
    Scalar   modeled_noise_per_rad = 0.01;
    bool     per_sensor_kf   = false;
    Scalar   kf_process_noise = 1.0;
    bool     scale_calib     = true;
    bool     bias_states     = false;   // GPS/INS-style drift removal (D22; Slice 11b)
    // Random-walk process-noise rate(s) for this source's body-twist bias (Slice 11b Option
    // A / Slice 18 Option B). Only used when bias_states == true. Option A adds Q_bias =
    // rate * dt * I6 to the bias block each predict (uniform form only); Option B
    // (multi_bias_enabled) applies the rate PER DOF — Q_bias = diag(rate) * dt — and a DOF
    // whose rate is exactly 0 is PINNED at zero (see BiasProcessNoise above). Larger => the
    // bias tracks faster but is noisier; >= 0 per DOF (validate enforces). Default small
    // uniform (slow constant-bias assumption).
    BiasProcessNoise bias_process_noise = BiasProcessNoise(Scalar(1e-4));
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

    // Per-channel SPLIT median (Slice 19, D3 amendment). When enabled, the fusion solve
    // routes to median::solve_split — TWO independent Weiszfeld medians (rotation on SO(3)
    // geodesic, translation on R^3) with per-channel weights (w_rot additionally carries
    // SensorConfig::rot_weight_prior), per-channel spreads driving a per-channel adaptive Q
    // (q_trans from spread_trans, q_rot from spread_rot — shared q_scale, no lambda
    // unit-mixing), and a block-diagonal Slice-18 median influence. Lets fusion express
    // per-channel source quality: the fused heading can track the best heading source
    // without distorting translation (the KAIST urban12 FOG case). Default OFF = the
    // coupled solver runs BYTE-IDENTICALLY (the split solver is mathematically different
    // even at equal weights — the coupled IRLS couples the channels — hence opt-in,
    // flip-default-later after real-data soak). Hashed (a flip rejects stale restores).
    bool       split_median   = false;
    // Cross-channel outlier veto for the split path (Slice 19; only read when split_median
    // is on). A source whose channel distance exceeds kVetoNormDist x the leave-one-out
    // channel spread in EITHER channel gets its OTHER-channel weight scaled by
    // kVetoWeightScale and that channel is re-solved once (median.hpp) — recovers the
    // coupled solver's whole-source rejection for hard (usually both-channel) faults while
    // keeping graceful per-channel weighting for quality differences. Default ON for the
    // split path (the safe policy); turn OFF only when a deliberately rotation-faulty-but-
    // translation-good source must keep full translation participation. Constants are
    // fixed (median.hpp) — only this bool is a knob. Hashed.
    bool       split_veto     = true;

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
    // Turn-regime FULL rotation-extrinsic via axis-correspondence hand-eye (Slice 17).
    // When enabled, Phase 2 additionally accumulates a per-source Wahba problem over the
    // turn-gated windowed rotation axes (a = log R_A, b = log R_B; a = R_X·b) and, once TWO
    // distinct rotation axes have been excited (the BBw rank gate), votes the full recovered
    // R_X into 3 so(3) channels; a committed rot3d SUPERSEDES the Phase-1 yaw/pitch ∘ roll
    // composition in the published prior_extrinsic rotation. On planar (yaw-only) motion the
    // rank gate never opens, so behavior is unchanged. Default OFF = byte-identical (sim +
    // every existing config); in the persistence config-hash (a flip rejects stale restores).
    bool        rot3d_enabled     = false;
    // Turn-regime JOINT lever + scale via the 4-unknown hand-eye LS (Slice 17b). When
    // enabled, Phase 2's lever rows become J = [(R_A − I) | −(R_X t_B)] with rhs −t_A and
    // unknowns [t_X; κ] (κ = 1/s_res, the residual-scale inverse — t_B arrives already
    // prior-de-scaled, so the κ prior is 1). Observable lever axes vote into the existing
    // xyz histograms; the κ axis votes s_res = 1/κ̂ into a per-source Phase-2 scale
    // histogram (configured from scale_hist), giving a SECOND, turn-regime scale estimator
    // (scale2) that commits/folds into prior_scale alongside Phase-1's straight-regime
    // path. Recovers scale on rigs with no straight regime (drones) and immunizes the
    // lever solve against an unrecovered scale. RANGE SEMANTICS (review MAJOR-2):
    // scale_hist's (range_min, range_max) BOUNDS the residual this path can recover per
    // fold — an out-of-range residual is SKIPPED at the vote site (never edge-clamped:
    // deterministic out-of-regime mass must not concentrate, commit, and irreversibly
    // poison prior_scale), so on a turn-only rig a true scale error beyond the range
    // (default: >50%) stays uncorrected at confidence 0 (diagnosable via
    // Phase2Calibrator::scale2_skipped / Estimator::scale2_skipped). Phase-1's
    // straight-regime scale vote has NO such guard — its edge-clamped fold iterates
    // toward an out-of-range truth on rigs WITH a straight regime; the two scale paths
    // deliberately differ here. Default OFF = byte-identical 3-unknown lever numerics
    // (this changes lever numbers even without scale activity, so it cannot ride an
    // existing knob); in the persistence config-hash (a flip rejects stale restores).
    bool        joint_lever_scale = false;

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
    // With joint_lever_scale, this range also BOUNDS the turn-regime (scale2) path's
    // recoverable residual per fold — see the joint_lever_scale comment above.
    HistogramConfig scale_hist{ /*bins=*/64, /*range_min=*/0.5, /*range_max=*/1.5 };
    HistogramConfig offset_hist;

    // Absolute refs
    Scalar mahalanobis_chi2 = 9.0;
    // Robust correction update (Slice 15): Huber down-weighting of the absolute-ref Kalman gain
    // for large innovations. For a per-DOF RMS Mahalanobis distance dbar = sqrt(NIS/n) above
    // `correction_robust_kappa`, the active measurement-noise block is inflated by dbar/kappa (gain
    // shrinks ~kappa/dbar), smoothly bounding the injected correction (esp. the heading kick fed by
    // the pose trans-rot cross-cov) instead of a binary gate. 0 = DISABLED -> the update is
    // bit-identical to the non-robust path. Typical enabled value ~3 (3 sigma per DOF).
    Scalar correction_robust_kappa = 0.0;
    // Bounded heading injection from position-only corrections (Slice 15b, lever C4). A correction
    // that does NOT observe rotation (H rotation-error columns ~ 0, e.g. a GPS position fix) can
    // still rotate the heading through the pose trans-rot cross-covariance in K. Under a LARGE
    // residual that coupling is outside linear validity and over-rotates (the urban12 t=1929 s 68
    // deg kick). When dbar = sqrt(NIS/n) > this kappa, the ROTATION rows of the gain are scaled by
    // kappa/dbar (translation correction untouched), bounding the injected heading. 0 = DISABLED
    // (bit-identical). Orthogonal to correction_robust_kappa (which scales the WHOLE gain). A
    // typical enabled value is below the gate's per-DOF dbar ceiling sqrt(mahalanobis_chi2/3)
    // (~1.73 at the default gate), e.g. ~0.8, so an in-gate-but-large position residual is protected.
    Scalar correction_rot_suppress_kappa = 0.0;
    // Median-coupled multi-source bias states (Slice 18, 11b Option B). When enabled, EVERY
    // source with SensorConfig::bias_states (up to Eskf::kMaxBiasSources = 4) carries its own
    // body-twist bias state in a generalized augmented ESKF [pose; twist; bias_1..bias_k]:
    // each biased source's SOURCE-FRAME delta is de-biased B_i' = B_i ∘ exp(−b_i·dt) BEFORE
    // the frame-align + median (the calibrators consume the SAME de-biased deltas), and the
    // pose<->bias coupling uses the exact FD-verified median-influence block Ω_i (see
    // eskf.hpp). Biases are observable only with an absolute ref. Default OFF = byte-identical
    // legacy behavior: the Option-A single-bias path is unchanged and >1 bias_states source
    // still rejects with InvalidConfig at validate(). In the persistence config-hash (a flip
    // rejects stale restores).
    //
    // FRAME NOTE (Slice-18 review MINOR): the published CalibSnapshot::bias changes FRAME
    // with this flag — Option A (off) learns a bias of the frame-ALIGNED delta (de-bias
    // after the align), Option B (on) learns a SOURCE-FRAME bias (de-bias before the
    // align; the coupling carries Ad(X_i)). Consumers comparing the vector across modes
    // must account for the extrinsic rotation.
    //
    // RELIABILITY INTERACTION (Slice-18 review MAJOR-4): the Slice-9 variance-EMA
    // reliability and the bias learning transient FEED BACK on each other — a co-source
    // whose bias is transiently mis-attributed scatters vs the consensus, gets
    // down-weighted, loses median influence, and its bias error becomes LESS observable
    // (it persists longer). The combination is covered by a dedicated acceptance case
    // (no divergence; the coast still beats the unmodeled baseline) but converges slower
    // than with reliability disabled; deployments wanting the fastest bias separation can
    // pin reliability_floor == reliability_cap == 1 during a calibration phase.
    bool   multi_bias_enabled = false;
    // Bias-prior variance seed per (unpinned) bias DOF for the multi-bias filter — Slice-18
    // review MAJOR-3 promoted this from a compile-time constant to a knob: it is a
    // RIG-DEPENDENT stability parameter, not a universal constant. sigma ~ 0.2 m/s / rad/s
    // by default, deliberately between Option A's uninformative 1.0 and a tight informative
    // prior; both extremes fail concretely (measured, tests/test_multi_bias_sim.cpp prior
    // sweep):
    //   * 1.0 (Option-A parity): the per-window coupling linearizes a heavily NONLINEAR
    //     median response, so the large prior gain injects spurious corrections along the
    //     weakly-observable bias DIFFERENCES — co-source biases blow up toward ~prior-scale
    //     garbage while the planted one overshoots.
    //   * 0.01: informative against any real bias — the estimate is shrunk by exactly the
    //     unresolved-variance fraction (recovered/planted tracks bias_observable), so the
    //     de-bias plateaus well short of the planted value.
    // The safe scale is set by the coupling's linearization radius: the per-window bias
    // displacement sigma*dt must stay below the per-window noise scale d_i for the rig
    // (sensor noise floors, window length, source count) — RE-TUNE for rigs whose noise
    // floor/cadence differs materially from the sim defaults. > 0 (validate). Hashed
    // (calibration-shaping: it seeds the learned-bias state a restore carries).
    Scalar multi_bias_cov0 = 0.04;

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
