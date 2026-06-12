// Slice 19 estimator-level tests: per-channel split median wired into the fusion spine.
//
// Coverage (SLICE19 §2 items 1, 3, 4, 6 + the multi-bias smoke of item 7):
//   * ESTIMATOR HEADLINE — a source with a gross ROTATION fault but the (uniquely
//     informative) correct translation: under split_median the fused rotation rejects the
//     fault while the fused translation tracks that same source — and the coupled solver
//     is pinned as unable to do both with any single scalar weight (run twice coupled,
//     once per weighting). Removing the split routing fails this test (the mutation).
//   * DEFAULT-OFF BYTE-IDENTITY — split_median=false is bit-identical to the pre-slice
//     path, and the new knobs (split_veto, rot_weight_prior) are INERT while off.
//   * rot_weight_prior — a 10x rotation prior pulls the fused rotation toward that source
//     under disagreement while the same step's fused translation is bit-unchanged.
//   * SPLIT-ON OBSERVABILITY — with split_median on (veto default ON), every calibration
//     DOF still converges in its regime and NOT in others (the calibrators consume the
//     split consensus).
//   * multi_bias + split smoke — the block-diagonal coupling path runs sane (finite, bias
//     unobservable without an absolute ref). The FD pin itself is tests/test_multi_bias.cpp.
#include <doctest/doctest.h>

#include "ofc/core/buffer.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/source.hpp"

#include "ofc_sim/rig.hpp"
#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <cmath>
#include <memory>
#include <vector>

using namespace ofc;

namespace {
Timestamp secs(double s) { return static_cast<Timestamp>(s * 1e9); }

// A source = an ISource façade over a SourceBuffer (the test_fusion.cpp adapter shape).
class BufferSource : public ISource {
public:
    BufferSource(SourceId id, int cap) : id_(id) {
        SensorConfig sc;
        sc.id = id;
        buf_.configure(sc, OdomForm::Twist, cap, ConfCombine::NativeOnly);
    }
    SourceId id() const override { return id_; }
    Expected<Delta> query(Timestamp t0, Timestamp t1) const override {
        return buf_.query(t0, t1);
    }
    SourceBuffer& buffer() { return buf_; }

private:
    SourceId     id_;
    SourceBuffer buf_;
};

// Fill a buffer with a CONSTANT body twist at `rate` Hz over [0, span] s.
void fill_twist(BufferSource& src, double span, double rate, const Vec6& xi) {
    const Timestamp step = secs(1.0 / rate);
    const int n = static_cast<int>(span * rate);
    for (int k = 0; k <= n; ++k) {
        src.buffer().push_twist(static_cast<Timestamp>(k) * step, xi, Mat6::Identity());
    }
}

Config split_base_config(int max_sources) {
    Config c;
    c.max_sources    = max_sources;
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.cold_start     = ColdStart::MedianFromStart;
    // Widen the weight clamp so weight_prior RATIOS survive it (the default [0.05, 10]
    // clamp saturates the sigma-confidence-scaled weights this noise-free rig produces,
    // flattening any per-source prior to equality).
    c.weight_floor   = 0.001;
    c.weight_cap     = 1000.0;
    return c;
}
} // namespace

// ---------------------------------------------------------------------------
// ESTIMATOR HEADLINE (SLICE19 §2 item 1, wired): rotation fault rejected, the SAME
// source's translation drives the fused translation — coupled pinned as unable.
// ---------------------------------------------------------------------------
TEST_CASE("split estimator HEADLINE: rotation-faulty source's translation still drives the "
          "fused translation; coupled solver cannot do both (mutation: split routing)") {
    // GT body motion: 2.0 m/s forward, no rotation.
    //   Sources 0, 1: clean rotation, translation biased +30% (2.6 m/s) — they agree.
    //   Source  2:    gross yaw-rate FAULT (+1.5 rad/s) but the CORRECT 2.0 m/s forward.
    // Per-channel weighting (the split capability): source 2 carries weight_prior 2.5 (its
    // translation is trusted — mass 2.5 > the 2.0 inlier mass) with rot_weight_prior 0.4
    // (equalizing its rotation weight, where the fault is rejected 2-vs-1 by the median).
    Vec6 xi_biased; xi_biased << 2.6, 0, 0, 0, 0, 0;
    Vec6 xi_fault;  xi_fault  << 2.0, 0, 0, 0, 0, 1.5;
    const double kSpan = 2.5, kRate = 200.0, kRunTo = 1.9;   // < reliability warmup

    struct RunOut { Scalar yaw_err; Scalar x_err; };
    // GT at the last frontier: x = 2.0 * t, yaw = 0.
    auto run = [&](bool split, bool veto, Scalar w2, Scalar rw2) -> RunOut {
        std::vector<std::unique_ptr<BufferSource>> srcs;
        for (int i = 0; i < 3; ++i)
            srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 4096));
        fill_twist(*srcs[0], kSpan, kRate, xi_biased);
        fill_twist(*srcs[1], kSpan, kRate, xi_biased);
        fill_twist(*srcs[2], kSpan, kRate, xi_fault);

        SensorConfig sensors[3];
        for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
        sensors[2].weight_prior     = w2;
        sensors[2].rot_weight_prior = rw2;
        Config cfg = split_base_config(3);
        cfg.split_median = split;
        cfg.split_veto   = veto;
        cfg.sensors      = sensors;
        cfg.sensor_count = 3;

        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

        double now_s = cfg.fusion_delay_s + cfg.window_s;
        double last_t1 = 0.0;
        while (now_s <= kRunTo) {
            CHECK(est.step(secs(now_s)) == Status::Ok);
            last_t1 = now_s - cfg.fusion_delay_s;
            now_s += cfg.window_s;
        }
        const State& f = est.latest().frontier;
        RunOut out;
        out.yaw_err = so3::log(f.pose.R).norm();                       // GT yaw = 0
        out.x_err   = std::abs(f.pose.t.x() - 2.0 * last_t1);          // GT x = 2 t
        return out;
    };

    // --- SPLIT (veto OFF — the per-channel-grace policy this scenario needs): BOTH right.
    const RunOut sp = run(/*split=*/true, /*veto=*/false, 2.5, 0.4);
    MESSAGE("split: yaw_err=" << sp.yaw_err << " x_err=" << sp.x_err);
    CHECK(sp.yaw_err < 0.1);     // rotation fault rejected (would be ~2.7 rad if followed)
    CHECK(sp.x_err < 0.15);      // the faulty source's CORRECT translation won (not +30%)

    // --- COUPLED, same weights: the fault's 2.5 mass drags the WHOLE source in — the
    //     fused rotation follows the fault (one scalar weight cannot split the channels).
    const RunOut c_hi = run(/*split=*/false, /*veto=*/false, 2.5, 1.0);
    MESSAGE("coupled w=2.5: yaw_err=" << c_hi.yaw_err << " x_err=" << c_hi.x_err);
    CHECK(c_hi.yaw_err > 0.5);

    // --- COUPLED, equal weights: the fault source is rejected whole, so its correct
    //     translation is LOST and the fused translation keeps the +30% inlier bias.
    const RunOut c_eq = run(/*split=*/false, /*veto=*/false, 1.0, 1.0);
    MESSAGE("coupled w=1: yaw_err=" << c_eq.yaw_err << " x_err=" << c_eq.x_err);
    CHECK(c_eq.yaw_err < 0.1);
    CHECK(c_eq.x_err > 0.6);

    // --- SPLIT with veto ON (the default policy): the gross rotation fault vetoes the
    //     source's translation too — whole-source rejection, the translation falls back
    //     to the biased inliers. Documents the policy fork the split_veto knob selects.
    const RunOut sv = run(/*split=*/true, /*veto=*/true, 2.5, 0.4);
    MESSAGE("split veto-on: yaw_err=" << sv.yaw_err << " x_err=" << sv.x_err);
    CHECK(sv.yaw_err < 0.1);
    CHECK(sv.x_err > 0.6);
}

// ---------------------------------------------------------------------------
// DEFAULT-OFF byte-identity (SLICE19 §2 item 4): split_median=false is bit-identical,
// and the new knobs are inert while off.
// ---------------------------------------------------------------------------
TEST_CASE("split default-off: byte-identical fused output; split_veto/rot_weight_prior "
          "inert while split_median is false") {
    Vec6 xi0; xi0 << 2.0, 0, 0, 0, 0, 0.3;
    Vec6 xi1; xi1 << 2.1, 0.05, 0, 0, 0, 0.32;
    Vec6 xi2; xi2 << 1.9, -0.05, 0, 0, 0, 0.28;

    struct Capture {
        std::vector<Scalar> vals;
    };
    auto run = [&](bool set_new_knobs) -> Capture {
        std::vector<std::unique_ptr<BufferSource>> srcs;
        for (int i = 0; i < 3; ++i)
            srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 4096));
        fill_twist(*srcs[0], 3.0, 200.0, xi0);
        fill_twist(*srcs[1], 3.0, 200.0, xi1);
        fill_twist(*srcs[2], 3.0, 200.0, xi2);

        SensorConfig sensors[3];
        for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
        Config cfg = split_base_config(3);
        cfg.sensors      = sensors;
        cfg.sensor_count = 3;
        if (set_new_knobs) {
            // split_median stays FALSE — every new knob must be inert.
            cfg.split_median = false;
            cfg.split_veto   = false;             // non-default
            sensors[1].rot_weight_prior = 10.0;   // non-default
        }

        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

        Capture cap;
        double now_s = cfg.fusion_delay_s + cfg.window_s;
        while (now_s <= 2.5) {
            CHECK(est.step(secs(now_s)) == Status::Ok);
            const State& f = est.latest().frontier;
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) cap.vals.push_back(f.pose.R(r, c));
            for (int k = 0; k < 3; ++k) cap.vals.push_back(f.pose.t(k));
            for (int r = 0; r < 12; ++r) cap.vals.push_back(f.cov(r, r));
            now_s += cfg.window_s;
        }
        return cap;
    };

    const Capture a = run(/*set_new_knobs=*/false);   // pristine default config
    const Capture b = run(/*set_new_knobs=*/true);    // new knobs set, split OFF
    REQUIRE(a.vals.size() == b.vals.size());
    REQUIRE(a.vals.size() > 100);
    for (std::size_t i = 0; i < a.vals.size(); ++i) {
        REQUIRE(a.vals[i] == b.vals[i]);              // BIT-identical (exact equality)
    }
}

// ---------------------------------------------------------------------------
// rot_weight_prior (SLICE19 §2 item 3): a 10x rotation prior pulls the fused rotation
// toward that source under disagreement; the same step's fused translation is unchanged.
// ---------------------------------------------------------------------------
TEST_CASE("split rot_weight_prior: 10x rotation prior pulls the rotation median; the "
          "translation median is unaffected (single-step bit-equality)") {
    // Three sources, SAME translation (2.0 m/s forward), DISAGREEING yaw rates
    // {0, +0.5, -0.5}. Rotation median at equal weights = the middle (yaw 0); with a 10x
    // rotation prior on source 1 (+0.5) the rotation channel mass (10 > 2) pulls the
    // consensus onto source 1. The translation channel weights are IDENTICAL in both runs,
    // so the single-step fused translation must be bit-equal.
    Vec6 xi0; xi0 << 2.0, 0, 0, 0, 0, 0.0;
    Vec6 xi1; xi1 << 2.0, 0, 0, 0, 0, 0.5;
    Vec6 xi2; xi2 << 2.0, 0, 0, 0, 0, -0.5;

    struct StepOut { Scalar yaw; Vec3 t; };
    auto one_step = [&](Scalar rw1) -> StepOut {
        std::vector<std::unique_ptr<BufferSource>> srcs;
        for (int i = 0; i < 3; ++i)
            srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 1024));
        fill_twist(*srcs[0], 0.5, 200.0, xi0);
        fill_twist(*srcs[1], 0.5, 200.0, xi1);
        fill_twist(*srcs[2], 0.5, 200.0, xi2);

        SensorConfig sensors[3];
        for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
        sensors[1].rot_weight_prior = rw1;
        Config cfg = split_base_config(3);
        cfg.split_median = true;                   // the knob is split-path-only
        // Veto OFF for the bit-equality claim: the cross-channel veto DELIBERATELY couples
        // the channels (a changed rotation median changes the rotation flags, which scale
        // TRANSLATION weights) — item 3's "translation unaffected" is about the per-channel
        // weighting itself, so it is pinned with the veto coupling disabled.
        cfg.split_veto   = false;
        cfg.sensors      = sensors;
        cfg.sensor_count = 3;

        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);
        REQUIRE(est.step(secs(cfg.fusion_delay_s + cfg.window_s)) == Status::Ok);
        const State& f = est.latest().frontier;
        StepOut out;
        out.yaw = so3::log(f.pose.R).z();
        out.t   = f.pose.t;
        return out;
    };

    const StepOut neutral = one_step(1.0);
    const StepOut pulled  = one_step(10.0);
    MESSAGE("rot prior: neutral yaw=" << neutral.yaw << " pulled yaw=" << pulled.yaw);
    // Neutral: the median rotation is the middle source (yaw rate 0) -> ~0 over the window.
    CHECK(std::abs(neutral.yaw) < 0.005);
    // 10x prior: pulled onto source 1 (~0.5 rad/s * 0.1 s window).
    CHECK(pulled.yaw > 0.03);
    // The translation median is UNAFFECTED: bit-equal fused translation.
    CHECK((neutral.t.array() == pulled.t.array()).all());
}

// ---------------------------------------------------------------------------
// SPLIT-ON observability self-tests (SLICE19 §2 item 6): the calibrators consume the
// split consensus — every DOF still converges in its regime and NOT in others.
// ---------------------------------------------------------------------------
namespace {
// Calibration-friendly config (the test_calib_phase1 histogram shape) + split ON.
Config split_calib_config(std::vector<SensorConfig>& sensors_out,
                          const std::vector<sim::SourceParams>& planted) {
    sensors_out.clear();
    for (const sim::SourceParams& sp : planted) {
        SensorConfig sc;
        sc.id = sp.id;                 // priors stay IDENTITY/1.0: the calibrators must
        sensors_out.push_back(sc);     // RECOVER the planted mount from scratch
    }
    Config c;
    c.max_sources    = static_cast<int>(planted.size());
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.tick_rate_hz   = 50.0;
    c.cold_start     = ColdStart::MedianFromStart;
    c.timesync_enabled = false;        // offset DOF: time-sync samples RAW source ‖ω‖
                                       // (never the median), untouched by the split — out
                                       // of scope here.
    c.split_median   = true;           // the path under test (veto stays default ON)
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.so3_hist.bins      = 512;
    c.so3_hist.range_min = -0.8;
    c.so3_hist.range_max =  0.8;
    c.so3_hist.aging     = Aging::SlidingK;
    c.so3_hist.sliding_k = 256;
    c.scale_hist.bins      = 512;
    c.scale_hist.aging     = Aging::SlidingK;
    c.scale_hist.sliding_k = 256;
    c.roll_hist.bins      = 360;
    c.roll_hist.circular  = true;
    c.roll_hist.range_min = -3.14159265358979323846;
    c.roll_hist.range_max =  3.14159265358979323846;
    c.xyz_hist.bins      = 512;
    c.xyz_hist.range_min = -1.5;
    c.xyz_hist.range_max =  1.5;
    c.xyz_hist.aging     = Aging::SlidingK;
    c.xyz_hist.sliding_k = 256;
    c.sensors      = sensors_out.data();
    c.sensor_count = static_cast<int>(sensors_out.size());
    return c;
}

SE3 yaw_pitch_extrinsic(Scalar yaw, Scalar pitch) {
    Mat3 Rz; Rz << std::cos(yaw), -std::sin(yaw), 0,
                   std::sin(yaw),  std::cos(yaw), 0,
                   0, 0, 1;
    Mat3 Ry; Ry << std::cos(pitch), 0, std::sin(pitch),
                   0, 1, 0,
                  -std::sin(pitch), 0, std::cos(pitch);
    SE3 X; X.R = Rz * Ry; X.t = Vec3::Zero();
    return X;
}

struct CalibOut {
    Scalar ext_conf;     // yaw/pitch (+roll once voted) rotation confidence
    Scalar scale_conf;
    Scalar scale;
    Scalar trans_conf;   // lever (xyz) confidence
    Vec3   fwd_err;      // recovered-vs-planted forward axis error
};
CalibOut run_split_calib(const sim::Trajectory& tr, const SE3& X_planted, Scalar scale_p) {
    std::vector<sim::SourceParams> planted(3);
    planted[0].id = 0;
    planted[1].id = 1;
    planted[1].X = X_planted;
    planted[1].scale = scale_p;
    planted[2].id = 2;

    std::vector<std::unique_ptr<sim::SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new sim::SyntheticSource(tr, sp));
    std::vector<SensorConfig> sensors;
    Config cfg = split_calib_config(sensors, planted);

    sim::Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(rig.add_source(s.get()) == Status::Ok);
    rig.run(0.2, tr.duration_s() - 0.1, 50.0);

    const Result& r = rig.estimator().latest();
    CalibOut out;
    out.ext_conf   = r.calib[1].extrinsic_confidence;
    out.scale_conf = r.calib[1].scale_confidence;
    out.scale      = r.calib[1].scale;
    out.trans_conf = r.calib[1].translation_confidence;
    const Vec3 want = X_planted.R.transpose() * Vec3(1, 0, 0);
    out.fwd_err = r.calib[1].extrinsic.R.transpose() * Vec3(1, 0, 0) - want;
    return out;
}
} // namespace

TEST_CASE("split-ON observability: straight regime converges yaw/pitch + scale; turn-only "
          "does NOT (the calibrators consume the split consensus)") {
    const SE3 X_planted = yaw_pitch_extrinsic(0.20, -0.12);
    const Scalar scale_p = 1.25;

    // POSITIVE regime: straight fwd + reverse — Phase-1 (yaw/pitch + scale) converges.
    {
        sim::Trajectory tr;
        Vec6 fwd; fwd << 2.0, 0, 0, 0, 0, 0;
        Vec6 rev; rev << -2.0, 0, 0, 0, 0, 0;
        tr.add_segment(fwd, 3.0);
        tr.add_segment(rev, 3.0);
        const CalibOut o = run_split_calib(tr, X_planted, scale_p);
        MESSAGE("split straight: ext_conf=" << o.ext_conf << " scale=" << o.scale
                << " scale_conf=" << o.scale_conf << " fwd_err=" << o.fwd_err.norm());
        CHECK(o.ext_conf > 0.3);
        CHECK(o.scale_conf > 0.3);
        CHECK(std::abs(o.scale - scale_p) < 0.03);
        CHECK(o.fwd_err.norm() < 0.03);     // recovered forward axis on the planted mount
    }

    // NEGATIVE regime: turning-only — the straight gate starves Phase 1 (stays at prior,
    // zero confidence) even though fusion now runs the split consensus.
    {
        sim::Trajectory tr = sim::Trajectory::turning(2.0, 0.6, 6.0);
        const CalibOut o = run_split_calib(tr, X_planted, scale_p);
        MESSAGE("split turning: ext_conf=" << o.ext_conf << " scale_conf=" << o.scale_conf);
        CHECK(o.scale_conf == doctest::Approx(0.0));
        CHECK(o.scale == doctest::Approx(1.0));      // prior (1.0 * residual 1.0)
    }
}

TEST_CASE("split-ON observability: turn regime converges the lever arm; straight-only does "
          "NOT (Phase 2 under the split consensus)") {
    // Planted LEVER (translation extrinsic) — Phase-2's xyz LS needs the turn regime.
    SE3 X_lever;                       // identity rotation, planted lever arm
    X_lever.t = Vec3(0.30, -0.20, 0.10);

    // POSITIVE: a MULTI-AXIS turning trajectory (yaw + pitch, the canonical Phase-2
    // convergence regime — a planar yaw-only turn leaves the z lever axis unobservable) —
    // translation (lever) confidence rises.
    {
        sim::Trajectory tr;
        Vec6 t1; t1 << 2.0, 0, 0, 0,  0.3,  0.6;
        Vec6 t2; t2 << 2.0, 0, 0, 0, -0.3,  0.6;
        for (int rep = 0; rep < 2; ++rep) {
            tr.add_segment(t1, 2.0);
            tr.add_segment(t2, 2.0);
        }
        const CalibOut o = run_split_calib(tr, X_lever, 1.0);
        MESSAGE("split turn lever: trans_conf=" << o.trans_conf);
        CHECK(o.trans_conf > 0.1);
    }
    // NEGATIVE: straight-only — the turn gate starves Phase 2 (no lever votes).
    {
        sim::Trajectory tr = sim::Trajectory::straight(2.0, 6.0);
        const CalibOut o = run_split_calib(tr, X_lever, 1.0);
        MESSAGE("split straight lever: trans_conf=" << o.trans_conf);
        CHECK(o.trans_conf == doctest::Approx(0.0));
    }
}

// ---------------------------------------------------------------------------
// multi_bias + split smoke (item 7's estimator side; the FD pin is test_multi_bias.cpp):
// the block-diagonal coupling path runs end-to-end — finite state, bias unobservable
// (~0 confidence) without an absolute ref, bias estimate stays ~0.
// ---------------------------------------------------------------------------
TEST_CASE("split + multi_bias: block-diagonal coupling path runs sane (finite, bias "
          "unobservable without an absolute ref)") {
    Vec6 xi0; xi0 << 2.0, 0, 0, 0, 0, 0.3;
    Vec6 xi1; xi1 << 2.05, 0.02, 0, 0, 0, 0.31;
    Vec6 xi2; xi2 << 1.95, -0.02, 0, 0, 0, 0.29;

    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < 3; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 4096));
    fill_twist(*srcs[0], 3.0, 200.0, xi0);
    fill_twist(*srcs[1], 3.0, 200.0, xi1);
    fill_twist(*srcs[2], 3.0, 200.0, xi2);

    SensorConfig sensors[3];
    for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[1].bias_states = true;     // one multi-bias source through the split coupling
    Config cfg = split_base_config(3);
    cfg.split_median       = true;
    cfg.multi_bias_enabled = true;
    cfg.sensors      = sensors;
    cfg.sensor_count = 3;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

    double now_s = cfg.fusion_delay_s + cfg.window_s;
    while (now_s <= 2.5) {
        CHECK(est.step(secs(now_s)) == Status::Ok);
        now_s += cfg.window_s;
    }
    const Result& r = est.latest();
    CHECK(r.frontier.pose.t.allFinite());
    CHECK(r.frontier.pose.R.allFinite());
    CHECK(r.frontier.cov.allFinite());
    // No absolute ref: the bias never becomes observable and its estimate stays ~0 (the
    // multi-bias observability contract, unchanged under the split coupling).
    CHECK(r.calib[1].bias_observable < 0.05);
    CHECK(r.calib[1].bias.norm() < 1e-6);
}
