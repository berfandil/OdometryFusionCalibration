// adapter_test_util.hpp — shared, SELF-CONTAINED test fixtures for the adapters doctest
// target (relaxed edge). Deliberately does NOT depend on ofc_sim: the adapters test executable
// links only ofc_core + ofc_adapters + doctest, so we hand-roll a tiny constant-twist source
// here. This keeps the adapter tests robust against sim/tests subdir configuration order.
#ifndef OFC_ADAPTERS_TEST_UTIL_HPP
#define OFC_ADAPTERS_TEST_UTIL_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/source.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

namespace adptest {

using namespace ofc;

constexpr Scalar kNanosPerSec = Scalar(1e9);

inline Timestamp secs_to_ns(Scalar s) {
    return static_cast<Timestamp>(std::llround(s * kNanosPerSec));
}

// A deterministic constant-twist ISource: the body moves with a fixed twist xi=[v;omega], so
// the motion over [t0,t1] is exp(xi * (t1-t0)). An optional translation `scale` is applied
// (mirrors a per-source scale error) and an identity extrinsic is assumed. No noise → fully
// deterministic, so a threaded run and a single-threaded run produce identical Results.
class TwistSource : public ISource {
public:
    TwistSource(SourceId id, const Vec6& xi, Scalar scale = Scalar(1))
        : id_(id), xi_(xi), scale_(scale) {}

    SourceId id() const override { return id_; }

    Expected<Delta> query(Timestamp t0, Timestamp t1) const override {
        if (t1 <= t0) return Status::NoData;
        const Scalar dt = static_cast<Scalar>(t1 - t0) / kNanosPerSec;
        Vec6 step = xi_ * dt;
        SE3 m = se3::exp(step);
        m.t *= scale_;                 // planted translation scale
        Delta d;
        d.motion = m;
        d.cov = Mat6::Identity() * Scalar(1e-3);
        d.t0 = t0;
        d.t1 = t1;
        return d;
    }

private:
    SourceId id_;
    Vec6     xi_;
    Scalar   scale_;
};

// A two-source rig config (reference 0 + a scaled source 1). MedianFromStart so both sources
// fuse immediately; time-sync off (the threaded determinism + persistence tests do not need
// it). Sensor storage is caller-owned (the contract); the caller holds `sensors`.
inline Config make_rig_cfg(std::vector<SensorConfig>& sensors) {
    sensors.assign(2, SensorConfig{});
    sensors[0].id = 0;
    sensors[0].is_reference = true;
    sensors[1].id = 1;

    Config cfg;
    cfg.max_sources = 2;
    cfg.tick_rate_hz = 50.0;
    cfg.fusion_delay_s = 0.05;
    cfg.window_s = 0.10;
    cfg.timesync_enabled = false;
    cfg.cold_start = ColdStart::MedianFromStart;
    cfg.reference_sensor_id = 0;
    cfg.sensors = sensors.data();
    cfg.sensor_count = 2;
    return cfg;
}

// A pair of constant-twist sources (a gentle curve, so both translation + rotation move).
inline void make_rig_sources(std::vector<TwistSource>& srcs) {
    Vec6 xi;
    xi << 2.0, 0, 0, 0, 0, 0.3;     // forward + yaw
    srcs.clear();
    srcs.emplace_back(0, xi, Scalar(1.0));
    srcs.emplace_back(1, xi, Scalar(1.0));
}

} // namespace adptest
#endif // OFC_ADAPTERS_TEST_UTIL_HPP
