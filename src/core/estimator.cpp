// ofc/core/estimator.cpp — the caller-pumped fusion facade (Slice 2).
//
// STRICT CORE: all working memory lives in Impl and is sized from Config at init();
// step() allocates nothing; no exceptions; bounded loops.
//
// Slice 2 pipeline (no calibration, no absolute corrections yet — DESIGN §§4-5):
//   per step(now): frontier t1 = now - fusion_delay; t0 = t1 - window
//     for each registered source:
//       B = query(t0, t1)                     (source-frame windowed delta)
//       A = X o B o X^-1                       (frame-align to base; X = sensor->base
//                                               prior extrinsic — see convention below)
//       w = clamp(weight_prior * sigma_confidence)
//     median of {A_i} weighted by {w_i}        (split-metric Weiszfeld)
//     ESKF predict on the median delta (dt = window), adaptive Q from the spread
//     populate Result.frontier and Result.tip
//   Returns NotReady until >= 1 source covers the window (full lifecycle is Slice 3).
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
// is the prior x Sigma-confidence weight of DESIGN §4 (reliability EMA is Slice 9).
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

    // Causal frontier window [t0, t1].
    const Timestamp t1 = now - secs_to_ns(cfg.fusion_delay_s);
    const Timestamp t0 = t1 - secs_to_ns(cfg.window_s);

    // Collect each covered source's frame-aligned delta + weight.
    int n = 0;
    for (int i = 0; i < s.source_count; ++i) {
        const ISource* src = s.sources[i];
        if (src == nullptr) continue;
        const Expected<Delta> q = src->query(t0, t1);

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
    Mat6 q_pose;
    if (cfg.adaptive_q) {
        const Scalar floor[6] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
        q_pose = Eskf::adaptive_q(med.spread, Scalar(1), floor);
    } else {
        q_pose = Scalar(1e-4) * Mat6::Identity();
    }

    // ESKF predict on the median delta over the window. Anchor the odom frame at the
    // first fused tick (pose starts at identity — gauge anchored at init, DESIGN §7).
    if (!s.eskf_started) {
        s.eskf.init(SE3{}, Mat12::Identity());
        s.eskf_started = true;
    }
    s.eskf.predict(med.value, cfg.window_s, q_pose);

    // Frontier state.
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
    // TODO: per-sensor + histogram range checks.
    return Status::Ok;
}

} // namespace ofc
