// ofc/core/estimator.cpp — facade skeleton.
// STRICT CORE: no allocation after init(); no exceptions; bounded loops.
// Implementation is filled in per ISSUES.md (Slices 0-12). This is scaffold.
#include "ofc/core/estimator.hpp"

namespace ofc {

// All working memory lives here and is sized from Config at init().
struct Estimator::Impl {
    Config cfg{};
    Result result{};
    bool   inited = false;
    // TODO Slice 1: per-source ring buffers (fixed capacity)
    // TODO Slice 2: median solver + ESKF (pose+twist, 12x12)
    // TODO Slice 4: histogram bank
    // TODO Slice 5: time-sync
    // TODO Slices 6-8: calibration phases + commit/feedback
    // TODO Slice 3: lifecycle state machine
};

Estimator::Estimator() : impl_(nullptr) {}

Estimator::~Estimator() { delete impl_; }   // sole heap free; nothing after init

Status Estimator::init(const Config& cfg) {
    const Status s = validate(cfg);
    if (!ok(s)) return s;
    if (impl_ == nullptr) impl_ = new Impl();   // allocate-once at init
    impl_->cfg    = cfg;
    impl_->result = Result{};
    impl_->result.phase = Phase::Init;
    impl_->inited = true;
    return Status::Ok;
}

Status Estimator::add_source(const ISource* src) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)src;
    return Status::Ok;   // TODO Slice 1
}

Status Estimator::add_correction(const ICorrection* corr) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)corr;
    return Status::Ok;   // TODO Slice 11
}

Status Estimator::step(Timestamp now) {
    if (impl_ == nullptr || !impl_->inited) return Status::NotInitialized;
    (void)now;
    // TODO Slice 2: integrate window -> deltas -> median -> ESKF predict
    // TODO Slice 11: absolute-ref updates (Mahalanobis-gated)
    // TODO Slices 5-8: bounded calibration slice
    // TODO Slice 3: lifecycle transition + populate result
    return Status::NotReady;
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
