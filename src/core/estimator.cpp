// ofc/core/estimator.cpp — the caller-pumped fusion facade (Slice 2 + Slice 5 time-sync).
//
// STRICT CORE: all working memory lives in Impl and is sized from Config at init();
// step() allocates nothing; no exceptions; bounded loops.
//
// TIME-SYNC (Slice 5, D16). When cfg.timesync_enabled, each step() samples every
// source's ‖ω‖ over a small sub-interval ending at the frontier, feeds it (keyed to
// base time) into TimeSync, and runs a bounded update() (‖ω‖ xcorr -> per-source offset
// histogram). A source's committed offset (histogram concentration >= commit_concentration)
// then SHIFTS its fusion query interval earlier by `off` so its internal +off re-shift
// lands back on true base time (sign per D21 — REMOVES a planted offset); below the gate,
// the configured prior is used. The offset + confidence are surfaced in CalibSnapshot.
// With time-sync OFF the per-source offset is exactly the prior (default 0) — Slice-2
// behavior unchanged.
//
// Slice 2 pipeline (no extrinsic/scale calibration, no absolute corrections — DESIGN §§4-5):
//   per step(now): frontier t1 = now - fusion_delay
//     integration interval [q0, t1]:
//       q0 = last_t1        (steady state — picks up where the last fuse ended)
//       q0 = t1 - window_s  (bootstrap of the FIRST fuse — no prior frontier yet)
//     dt = (t1 - q0) seconds                   (ACTUAL elapsed motion, not a fixed step)
//     for each registered source:
//       B = query(q0 - off, t1 - off)         (source-frame delta; off = committed/prior
//                                               clock offset, removed per Slice 5)
//       B_corr = { B.R, B.t / prior_scale }   (de-scale the reported translation, D20)
//       A = X o B_corr o X^-1                  (frame-align to base; X = sensor->base
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
#include "ofc/core/timesync.hpp"

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

// ‖ω‖ extraction for time-sync (DESIGN §6 time-offset row, D16): the body angular-rate
// MAGNITUDE over a small sub-interval [t-h, t]. From the source's reported SE(3) delta
// over that window, ‖ω‖ ≈ ‖log(ΔR)‖ / h. This is extrinsic-invariant (a rotation
// magnitude survives the sensor->base conjugation), so it can be sampled per source
// BEFORE any spatial calibration and cross-correlated against the reference. dt is the
// sub-interval length in seconds (> 0).
Scalar omega_norm_from_delta(const SE3& motion, Scalar dt) {
    if (!(dt > Scalar(0))) return Scalar(0);
    return so3::log(motion.R).norm() / dt;
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

    // Per-source prior extrinsic (sensor->base), prior translation scale, and prior
    // weight, looked up from the SensorConfig matching the source id (identity / 1.0 /
    // 1.0 if no config provided). prior_scale de-scales the source's reported
    // translation before the frame-align (D20 — inverts the sim's per-source scale).
    SE3    prior_extrinsic[kMaxSourcesCap];
    Scalar prior_scale[kMaxSourcesCap];
    Scalar weight_prior[kMaxSourcesCap];
    // Per-source prior clock offset (seconds). Sign per D21 / CONFIG §9: positive =>
    // source clock ahead of base; its reported [t0,t1] reads true [t0+off, t1+off].
    Scalar prior_time_offset[kMaxSourcesCap];

    // Time-sync (Slice 5, D16). When cfg.timesync_enabled, each step() feeds every
    // source's current ‖ω‖ sample into `timesync` and periodically runs update(); the
    // committed per-source offset shifts that source's fusion query interval so the
    // planted clock skew is removed before fusion. Heap-resident inside Impl (the
    // strict-core no-post-init-heap contract: Impl is allocated once in init()).
    TimeSync timesync;
    bool     timesync_active = false;
    int      ticks_since_sync = 0;

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

    // Effective per-source clock offset (seconds) applied to its fusion query. When
    // time-sync is active AND the source's offset estimate is confident enough
    // (histogram concentration >= cfg.commit_concentration), use the estimated offset;
    // otherwise fall back to the configured prior (default 0). Slot `i` indexes the
    // registered-source arrays.
    Scalar effective_offset(int i) const {
        if (timesync_active) {
            const SourceId id = sources[i]->id();
            if (id != cfg.reference_sensor_id &&
                timesync.confidence(id) >= cfg.commit_concentration) {
                return timesync.offset(id);
            }
        }
        return prior_time_offset[i];
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
    impl_->ticks_since_sync = 0;
    for (int i = 0; i < kMaxSourcesCap; ++i) {
        impl_->sources[i]           = nullptr;
        impl_->prior_extrinsic[i]   = SE3{};
        impl_->prior_scale[i]       = Scalar(1);
        impl_->weight_prior[i]      = Scalar(1);
        impl_->prior_time_offset[i] = Scalar(0);
    }

    // Time-sync: configure once at init (preallocates its buffers/histograms). Disabled
    // configs skip it entirely (behave exactly as before — offsets stay at the priors).
    impl_->timesync_active = false;
    if (cfg.timesync_enabled) {
        const Status ts = impl_->timesync.configure(cfg, cfg.reference_sensor_id);
        if (!ok(ts)) return ts;
        impl_->timesync_active = true;
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
    // Bind the prior extrinsic + translation scale + weight for this source from its
    // SensorConfig (defaults: identity / 1.0 / 1.0 when no config matches).
    const SensorConfig* sc = impl_->sensor_for(src->id());
    impl_->prior_extrinsic[slot] = (sc != nullptr) ? sc->prior_extrinsic : SE3{};
    impl_->weight_prior[slot]    = (sc != nullptr) ? sc->weight_prior    : Scalar(1);
    // prior_scale must be positive (validate is expected to enforce this, but be safe:
    // a non-positive scale would blow up the de-scale division in step()).
    const Scalar ps = (sc != nullptr) ? sc->prior_scale : Scalar(1);
    impl_->prior_scale[slot] = (ps > Scalar(0)) ? ps : Scalar(1);
    impl_->prior_time_offset[slot] =
        (sc != nullptr) ? sc->prior_time_offset_s : Scalar(0);
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

    // ---- Time-sync sampling (Slice 5, D16) --------------------------------------
    // Feed each source's current ‖ω‖ sample (over a small sub-interval ending at the
    // frontier t1) into the cross-correlator, keyed to BASE/query time, then run a
    // bounded update() once per tick. The ‖ω‖ sub-interval is one common-grid period
    // (sample_dt) so samples land on the time-sync grid. No-op when time-sync is off.
    if (s.timesync_active) {
        const Scalar h_s  = s.timesync.sample_dt();
        const Timestamp h = secs_to_ns(h_s);
        if (h > 0 && t1 - h > 0) {
            for (int i = 0; i < s.source_count; ++i) {
                const ISource* src = s.sources[i];
                if (src == nullptr) continue;
                const Expected<Delta> wq = src->query(t1 - h, t1);
                if (!wq.ok()) continue;
                // NOTE: ‖ω‖ is sampled at the source's RAW timeline (no offset applied)
                // — the planted skew between two sources' ‖ω‖ streams is exactly what
                // the cross-correlation measures. Applying the offset here would cancel
                // the very signal we are trying to recover.
                const Scalar wn = omega_norm_from_delta(wq.value().motion, h_s);
                s.timesync.push(src->id(), t1, wn);
            }
        }
        s.timesync.update();
        ++s.ticks_since_sync;
    }

    // Collect each covered source's frame-aligned delta + weight.
    int n = 0;
    for (int i = 0; i < s.source_count; ++i) {
        const ISource* src = s.sources[i];
        if (src == nullptr) continue;
        // Apply the source's committed clock offset to its fusion query: shift the
        // interval EARLIER by `off` so the source's internal +off re-shift lands back
        // on true base time [q0, t1] (sign per D21 — removes a planted offset). When
        // time-sync is off / not yet confident this is the configured prior (default 0).
        const Timestamp off = secs_to_ns(s.effective_offset(i));
        const Expected<Delta> q = src->query(q0 - off, t1 - off);

        // Diagnostics (best-effort; coverage flag set regardless of fuse outcome).
        SourceHealth& h = s.result.health[i];
        h.id        = src->id();
        h.in_window = q.ok();
        h.weight    = Scalar(0);
        h.residual  = Scalar(0);

        if (!q.ok()) continue;

        const Delta& d = q.value();
        // De-scale the reported translation by the per-source prior_scale BEFORE the
        // frame-align (D20). The source reports B with B.t already multiplied by `scale`
        // (the sim applies scale to B.t AFTER the X-conjugation), so dividing B.t by the
        // matching prior_scale and then conjugating inverts it exactly:
        //   B_corr = { B.R, B.t / prior_scale } ;  A = X o B_corr o X^-1.
        // prior_scale is guaranteed > 0 by add_source (guarded), so this never /0.
        SE3 B_corr = d.motion;
        B_corr.t   = d.motion.t / s.prior_scale[i];
        // Frame-align: A = X o B_corr o X^-1  (X = sensor->base prior extrinsic).
        const SE3& X = s.prior_extrinsic[i];
        const SE3 A  = se3::compose(se3::compose(X, B_corr), se3::inverse(X));

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

    // Calibration snapshot (Slice 5 fills the time-offset DOF; extrinsic/scale stay at
    // their priors until Slices 6-8). Per source: the offset estimate currently APPLIED
    // to its query (committed time-sync value if confident, else the prior) plus the
    // time-sync histogram confidence. `committed` marks that the estimate cleared the
    // commit-concentration gate (so it is actually driving fusion rather than the prior).
    for (int i = 0; i < s.source_count; ++i) {
        const SourceId id = s.sources[i]->id();
        CalibSnapshot& cs = s.result.calib[i];
        cs.id            = id;
        cs.extrinsic     = s.prior_extrinsic[i];
        cs.scale         = s.prior_scale[i];
        cs.time_offset_s = s.effective_offset(i);
        const Scalar conf = s.timesync_active ? s.timesync.confidence(id) : Scalar(0);
        cs.confidence    = conf;
        cs.committed     = s.timesync_active &&
                           (id != cfg.reference_sensor_id) &&
                           (conf >= cfg.commit_concentration);
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
