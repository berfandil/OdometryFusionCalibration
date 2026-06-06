// ofc/core/estimator.cpp — the caller-pumped fusion facade (Slice 2).
//
// STRICT CORE: all working memory lives in Impl and is sized from Config at init();
// step() allocates nothing; no exceptions; bounded loops.
//
// Slice 2 pipeline (no calibration, no absolute corrections yet — DESIGN §§4-5):
//   per step(now): frontier t1 = now - fusion_delay
//     integration interval [q0, t1]:
//       q0 = last_t1        (steady state — picks up where the last fuse ended)
//       q0 = t1 - window_s  (bootstrap of the FIRST fuse — no prior frontier yet)
//     dt = (t1 - q0) seconds                   (ACTUAL elapsed motion, not a fixed step)
//     for each registered source:
//       B = query(q0, t1)                     (source-frame delta over the interval)
//       A = X o B o X^-1                       (frame-align to base; X = sensor->base
//                                               prior extrinsic — see convention below)
//       w = clamp(weight_prior * sigma_confidence)
//     median of {A_i} weighted by {w_i}        (split-metric Weiszfeld)
//     ESKF predict on the median delta (dt = t1 - q0), adaptive Q from the spread
//     on success: last_t1 <- t1               (gap/overlap-free across ticks)
//     populate Result.frontier and Result.tip
//   Returns NotReady until >= 1 source covers the interval (full lifecycle is Slice 3).
//   NOTE: `window_s` is the bootstrap/lookback interval for the first fuse only; the
//   integrator no longer assumes one window per tick, so tick cadence may differ from
//   window_s without opening integration gaps/overlaps (review fix — DESIGN §7).
//
// Frame-align convention (documented): each SensorConfig::prior_extrinsic X is the
// sensor->base transform (a point/twist expressed in the sensor frame maps to base
// via X). A motion B measured in the sensor frame becomes, in the base frame, the
// conjugation A = X o B o X^-1. With X = identity (aligned mount) A == B.
#include "ofc/core/estimator.hpp"

#include "ofc/core/buffer.hpp"
#include "ofc/core/eskf.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/median.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
constexpr Scalar kNanosPerSec = Scalar(1e9);
constexpr int    kMaxSourcesCap = 32;   // matches Result::calib/health array sizes

Timestamp secs_to_ns(Scalar s) {
    return static_cast<Timestamp>(std::llround(s * kNanosPerSec));
}

// Sigma-confidence of a windowed delta: an inverse-covariance scalar. We use
// 1 / (mean diagonal + eps) so a tighter (smaller-Sigma) source weighs more. This
// is the prior x Sigma-confidence weight of DESIGN §4.
// PLACEHOLDER (unit-mixing): the mean over all 6 diagonal entries blends
// translation (m^2) and rotation (rad^2) variances into one scalar, so a source
// with tiny rad var but huge m var gets a misleading confidence. This is an
// accepted Slice-2 placeholder until the Slice 9 reliability EMA replaces it with a
// properly unit-separated weight.
Scalar sigma_confidence(const Mat6& cov) {
    const Scalar mean_diag = cov.diagonal().mean();
    constexpr Scalar kEps = Scalar(1e-9);
    return Scalar(1) / (std::max(Scalar(0), mean_diag) + kEps);
}
} // namespace

// All working memory lives here and is sized from Config at init().
struct Estimator::Impl {
    Config cfg{};
    Result result{};
    bool   inited = false;

    // Registered sources (pointers owned by caller). Fixed array sized by the cap;
    // the live count is bounded by cfg.max_sources.
    const ISource* sources[kMaxSourcesCap] = {};
    int            source_count = 0;

    // Per-source prior extrinsic (sensor->base) and prior weight, looked up from the
    // SensorConfig matching the source id (identity / 1.0 if no config provided).
    SE3    prior_extrinsic[kMaxSourcesCap];
    Scalar weight_prior[kMaxSourcesCap];

    // Scratch for the per-step fuse (no heap in step()).
    SE3    aligned[kMaxSourcesCap];
    Scalar weights[kMaxSourcesCap];

    Eskf   eskf;
    bool   eskf_started = false;

    // End of the last successfully-fused window. The next predict integrates the
    // actual motion over [last_t1, t1] (dt = t1 - last_t1), so integration is
    // gap/overlap-free regardless of how the caller spaces step() ticks. Only
    // advanced on a successful fuse; a step with no covering source leaves it put.
    Timestamp last_t1     = 0;
    bool      has_last_t1 = false;

    // Resolve the SensorConfig for a given source id (linear scan, bounded by
    // sensor_count). Returns nullptr if none.
    const SensorConfig* sensor_for(SourceId id) const {
        if (cfg.sensors == nullptr) return nullptr;
        for (int i = 0; i < cfg.sensor_count; ++i) {
            if (cfg.sensors[i].id == id) return &cfg.sensors[i];
        }
        return nullptr;
    }
};

Estimator::Estimator() : impl_(nullptr) {}

Estimator::~Estimator() { delete impl_; }   // sole heap free; nothing after init

Status Estimator::init(const Config& cfg) {
    const Status s = validate(cfg);
    if (!ok(s)) return s;
    if (impl_ == nullptr) impl_ = new Impl();   // allocate-once at init
    impl_->cfg          = cfg;
    impl_->result       = Result{};
    impl_->result.phase = Phase::Init;
    impl_->source_count = 0;
    impl_->eskf_started = false;
    impl_->last_t1      = 0;
    impl_->has_last_t1  = false;
    for (int i = 0; i < kMaxSourcesCap; ++i) {
        impl_->sources[i]         = nullptr;
        impl_->prior_extrinsic[i] = SE3{};
        impl_->weight_prior[i]    = Scalar(1);
    }
    impl_->inited = true;
    return Status::Ok;
}

Status Estimator::add_source(const ISource* src) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    if (src == nullptr)                     return Status::InvalidConfig;
    if (impl_->source_count >= impl_->cfg.max_sources ||
        impl_->source_count >= kMaxSourcesCap) {
        return Status::CapacityExceeded;
    }
    const int slot = impl_->source_count;
    impl_->sources[slot] = src;
    // Bind the prior extrinsic + weight for this source from its SensorConfig.
    const SensorConfig* sc = impl_->sensor_for(src->id());
    impl_->prior_extrinsic[slot] = (sc != nullptr) ? sc->prior_extrinsic : SE3{};
    impl_->weight_prior[slot]    = (sc != nullptr) ? sc->weight_prior   : Scalar(1);
    ++impl_->source_count;
    return Status::Ok;
}

Status Estimator::add_correction(const ICorrection* corr) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)corr;
    return Status::Ok;   // TODO Slice 11
}

Status Estimator::step(Timestamp now) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    Impl& s = *impl_;
    const Config& cfg = s.cfg;

    // Causal frontier at t1 = now - delay. We integrate the ACTUAL motion since the
    // last fused frontier: the query/predict interval is [q0, t1] with
    //   q0 = last_t1          (steady state — gap/overlap-free between ticks)
    //   q0 = t1 - window_s    (bootstrap of the first fused step — no prior frontier)
    // and dt = (t1 - q0) in seconds. So `window_s` is now the BOOTSTRAP / LOOKBACK
    // interval for the very first fuse, not a fixed per-tick step: the integrator is
    // correct regardless of tick cadence (DESIGN §7). last_t1 only advances on a
    // successful fuse, so a skipped (uncovered) step never opens an integration gap.
    const Timestamp t1 = now - secs_to_ns(cfg.fusion_delay_s);
    const Timestamp q0 = s.has_last_t1 ? s.last_t1 : (t1 - secs_to_ns(cfg.window_s));

    // Degenerate / non-causal interval guard: nothing to integrate.
    if (t1 <= q0) {
        s.result.phase     = Phase::Warmup;
        s.result.tip_valid = false;
        return Status::NotReady;
    }
    const Scalar dt = static_cast<Scalar>(t1 - q0) / kNanosPerSec;

    // Collect each covered source's frame-aligned delta + weight.
    int n = 0;
    for (int i = 0; i < s.source_count; ++i) {
        const ISource* src = s.sources[i];
        if (src == nullptr) continue;
        const Expected<Delta> q = src->query(q0, t1);

        // Diagnostics (best-effort; coverage flag set regardless of fuse outcome).
        SourceHealth& h = s.result.health[i];
        h.id        = src->id();
        h.in_window = q.ok();
        h.weight    = Scalar(0);
        h.residual  = Scalar(0);

        if (!q.ok()) continue;

        const Delta& d = q.value();
        // Frame-align: A = X o B o X^-1  (X = sensor->base prior extrinsic).
        const SE3& X = s.prior_extrinsic[i];
        const SE3 A  = se3::compose(se3::compose(X, d.motion), se3::inverse(X));

        // Weight = prior x Sigma-confidence, clamped to [floor, cap].
        Scalar w = s.weight_prior[i] * sigma_confidence(d.cov);
        w = std::max(cfg.weight_floor, std::min(cfg.weight_cap, w));

        s.aligned[n] = A;
        s.weights[n] = w;
        h.weight     = w;
        ++n;
    }

    s.result.source_count = s.source_count;

    // Phase: need the window fully covered by at least one source. (Full lifecycle —
    // INIT/WARMUP/DEGRADED/NOMINAL — is Slice 3; keep this minimal but correct.)
    if (n == 0) {
        s.result.phase     = Phase::Warmup;
        s.result.tip_valid = false;
        return Status::NotReady;
    }

    // Robust consensus delta over the window.
    median::Params mp;
    mp.max_iters = cfg.weiszfeld_max_iters;
    mp.tol       = cfg.weiszfeld_tol;
    mp.eps       = cfg.weiszfeld_eps;
    mp.lambda    = cfg.metric_lambda;
    const median::Result med = median::solve(s.aligned, s.weights, n, mp);

    // Per-source residual to the consensus (diagnostics).
    for (int i = 0, k = 0; i < s.source_count; ++i) {
        if (!s.result.health[i].in_window) continue;
        s.result.health[i].residual =
            se3::split_distance(med.value, s.aligned[k], cfg.metric_lambda);
        ++k;
    }

    // Adaptive process noise from the inter-source spread (DESIGN §4, D4).
    // q_scale / q_floor come from Config (CONFIG §3). Adaptive: floor + q_scale*spread^2.
    // Non-adaptive: just the per-axis floor (no spread term).
    Mat6 q_pose;
    if (cfg.adaptive_q) {
        q_pose = Eskf::adaptive_q(med.spread, cfg.q_scale, cfg.q_floor);
    } else {
        q_pose = Eskf::adaptive_q(Scalar(0), cfg.q_scale, cfg.q_floor);
    }

    // ESKF predict on the median delta over the integrated interval [q0, t1] (dt is
    // the actual elapsed time, NOT a fixed window). Anchor the odom frame at the
    // first fused tick (pose starts at identity — gauge anchored at init, DESIGN §7).
    if (!s.eskf_started) {
        s.eskf.init(SE3{}, Mat12::Identity());
        s.eskf_started = true;
    }
    s.eskf.predict(med.value, dt, q_pose);

    // This fuse succeeded: advance the integration frontier so the next step picks
    // up exactly where this one ended (no gap, no overlap).
    s.last_t1     = t1;
    s.has_last_t1 = true;

    // Frontier state. Drive the published stamp from the real frontier t1.
    s.result.frontier        = s.eskf.state();
    s.result.frontier.stamp  = t1;
    s.result.phase           = Phase::Nominal;

    // Predicted tip: const-velocity extrapolation from the frontier to `now`.
    if (cfg.emit_predicted_tip) {
        const Scalar dt_ahead = cfg.fusion_delay_s;   // frontier (t1) -> now
        s.eskf.predict_tip(dt_ahead, cfg.tip_cov_inflation, s.result.tip);
        s.result.tip.stamp = now;
        s.result.tip_valid = true;
    } else {
        s.result.tip_valid = false;
    }

    return Status::Ok;
}

const Result& Estimator::latest() const { return impl_->result; }

Expected<int> Estimator::serialize(unsigned char* buf, int cap) const {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)buf; (void)cap;
    return Status::NoData;   // TODO Slice 12
}

Status Estimator::deserialize(const unsigned char* buf, int len) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)buf; (void)len;
    return Status::NoData;   // TODO Slice 12 (incl. config-hash guard)
}

// --- validate() lives with the config it checks --------------------------
Status validate(const Config& cfg) {
    if (cfg.max_sources < 1 || cfg.max_sources > 32) return Status::OutOfRange;
    if (cfg.buffer_capacity < 16)                    return Status::OutOfRange;
    if (cfg.sensor_count > cfg.max_sources)          return Status::CapacityExceeded;
    if (cfg.tick_rate_hz <= 0.0)                     return Status::OutOfRange;
    if (cfg.window_s <= 0.0)                         return Status::OutOfRange;
    if (cfg.weiszfeld_max_iters < 1)                 return Status::OutOfRange;
    if (cfg.metric_lambda <= 0.0)                    return Status::OutOfRange;
    if (cfg.commit_drop >= cfg.commit_concentration) return Status::InvalidConfig;
    if (cfg.reference_sensor_id >= cfg.max_sources)  return Status::OutOfRange;

    // Median / fusion knob ranges (CONFIG §§1-4).
    if (cfg.confidence_blend < 0.0 || cfg.confidence_blend > 1.0)
        return Status::OutOfRange;                       // blend factor in [0, 1]
    if (cfg.weiszfeld_tol <= 0.0 || cfg.weiszfeld_tol >= 1.0)
        return Status::OutOfRange;                       // tol in (0, 1)
    if (cfg.weiszfeld_eps <= 0.0 || cfg.weiszfeld_eps >= 1.0)
        return Status::OutOfRange;                       // eps in (0, 1)
    if (cfg.fusion_delay_s < 0.0 || cfg.fusion_delay_s > 2.0)
        return Status::OutOfRange;                       // delay in [0, 2] s
    if (cfg.weight_floor <= 0.0 || cfg.weight_floor >= 1.0)
        return Status::OutOfRange;                       // floor in (0, 1)
    if (cfg.weight_cap < 1.0)
        return Status::OutOfRange;                       // cap >= 1
    if (cfg.tip_cov_inflation < 1.0)
        return Status::OutOfRange;                       // inflation >= 1
    if (cfg.q_scale <= 0.0)
        return Status::OutOfRange;                       // q_scale > 0
    for (int i = 0; i < 6; ++i) {
        if (cfg.q_floor[i] < 0.0) return Status::OutOfRange;   // each q_floor[i] >= 0
    }

    // TODO: per-sensor + histogram range checks.
    return Status::Ok;
}

} // namespace ofc
