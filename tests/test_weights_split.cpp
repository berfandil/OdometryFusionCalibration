// Slice 19b tests: PER-CHANNEL scatter reliability under the split median (policy layer b).
//
// Under Config::split_median the Slice-9 residual feed becomes PER-CHANNEL: two independent
// variance-EMA tracks per source — d_rot = ||log(R_consensus^T R_i)|| (rad) and
// d_trans = ||t_i - t_consensus|| (m) against the split consensus — so a source noisy in
// ROTATION but clean in TRANSLATION is downweighted ONLY in the rotation channel (and vice
// versa). The COUPLED path (split_median = false) keeps the mixed split_distance EMA + the
// scalar reliability BYTE-IDENTICALLY (pinned below against HEAD-captured literals).
//
// Coverage (the slice-19b brief's test list):
//   a. HEADLINE (fields): a rot-noisy/trans-clean source reads reliability_rot << 1 with
//      reliability_trans ~ 1, and the symmetric trans-noisy case — via SourceHealth.
//   a'. HEADLINE (effect): the rot-noisy source's TRANSLATION still drives the fused
//      translation after warmup (the mixed residual would have floored BOTH channels and
//      lost it — the mutation this slice exists to kill), while the fused rotation tracks
//      the clean pair.
//   b. D17 PER CHANNEL: a systematic rotation BIAS with low scatter is NOT rot-downweighted
//      (bias -> calibrator, never the weight — now per channel).
//   c. Warmup: per-channel reliabilities stay exactly 1.0 before kRelWarmup samples.
//   d. COUPLED-PATH PIN: split_median=false results EXACTLY equal the pre-change values
//      (HEAD-captured literals on a deterministic seeded rig) and the new SourceHealth
//      fields stay at their defaults (1, 1, 0, 0) on the coupled path.
//   e. rot_weight_prior composes MULTIPLICATIVELY with reliability_rot (a large enough
//      datasheet prior overrides the learned floor — and a neutral one does not).
#include <doctest/doctest.h>

#include "ofc/core/buffer.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/source.hpp"

#include "ofc_sim/rig.hpp"
#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;
bool near_abs(Scalar a, Scalar b, Scalar tol) { return std::abs(a - b) <= tol; }
Timestamp secs(double s) { return static_cast<Timestamp>(s * 1e9); }

// The test_weights.cpp mixed straight+turn trajectory (long enough to pass kRelWarmup).
Trajectory mixed_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0, 0.35,  0.6;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.35, 0.6;
    Vec6 slowturn; slowturn << 2.0, 0, 0, 0, 0.0,   0.25;
    for (int rep = 0; rep < 3; ++rep) {
        tr.add_segment(straight, 2.0);
        tr.add_segment(turnA,    1.6);
        tr.add_segment(slowturn, 1.0);
        tr.add_segment(straight, 1.6);
        tr.add_segment(turnB,    1.6);
        tr.add_segment(slowturn, 1.0);
    }
    return tr;
}

// Histogram knobs (the test_weights.cpp shape; calibrator side only).
void set_hists(Config& c) {
    c.vote_weight = VoteWeight::One;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.turn_omega_min     = 0.20;

    c.so3_hist.bins = 256; c.so3_hist.range_min = -0.8; c.so3_hist.range_max = 0.8;
    c.so3_hist.circular = false; c.so3_hist.aging = Aging::SlidingK; c.so3_hist.sliding_k = 200;
    c.so3_hist.vote_split = true; c.so3_hist.subbin = true;

    c.scale_hist.bins = 256; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.circular = false; c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 200;
    c.scale_hist.vote_split = true; c.scale_hist.subbin = true;

    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 200;
    c.roll_hist.vote_split = true; c.roll_hist.subbin = true;

    c.xyz_hist.bins = 256; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.circular = false; c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 200;
    c.xyz_hist.vote_split = true; c.xyz_hist.subbin = true;

    c.offset_hist.bins = 256; c.offset_hist.range_min = -0.1; c.offset_hist.range_max = 0.1;
    c.offset_hist.circular = false; c.offset_hist.aging = Aging::SlidingK; c.offset_hist.sliding_k = 200;
    c.offset_hist.vote_split = true; c.offset_hist.subbin = true;
}

const SourceHealth* health(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.health[i].id == id) return &r.health[i];
    return nullptr;
}

Config base_cfg(std::vector<SensorConfig>& sensors, int n) {
    Config c;
    c.max_sources    = n;
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.tick_rate_hz   = 50.0;
    c.timesync_enabled = false;
    c.cold_start     = ColdStart::MedianFromStart;
    set_hists(c);
    c.sensors      = sensors.data();
    c.sensor_count = static_cast<int>(sensors.size());
    return c;
}

// A source = an ISource façade over a SourceBuffer (the test_split_median.cpp shape).
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

// Fill with a body twist whose YAW RATE varies per 0.1 s window slot: wz_k = base_wz +
// amp * kNoiseTable[(k + phase) % N] (a fixed zero-mean VARYING-MAGNITUDE table — pure +/-
// alternation would read as constant DISTANCE, i.e. bias, not scatter, to the unsigned
// residual track; varying magnitudes are what the variance EMA keys on).
constexpr double kNoiseTable[8] = {0.9, -0.4, 0.7, -1.0, 0.3, -0.8, 0.5, -0.2};
void fill_twist_noisy_yaw(BufferSource& src, double span, double rate, const Vec6& xi_base,
                          double amp, double base_wz = 0.0, int phase = 0) {
    const Timestamp step = secs(1.0 / rate);
    const int n = static_cast<int>(span * rate);
    const int per_window = static_cast<int>(rate * 0.1);   // pushes per 0.1 s window slot
    for (int k = 0; k <= n; ++k) {
        const int slot = (per_window > 0) ? (k / per_window) : 0;
        Vec6 xi = xi_base;
        xi(5) += base_wz + amp * kNoiseTable[(slot + phase) % 8];
        src.buffer().push_twist(static_cast<Timestamp>(k) * step, xi, Mat6::Identity());
    }
}
void fill_twist(BufferSource& src, double span, double rate, const Vec6& xi) {
    fill_twist_noisy_yaw(src, span, rate, xi, 0.0);
}

Config split_base_config(int max_sources) {
    Config c;
    c.max_sources    = max_sources;
    c.fusion_delay_s = 0.05;
    c.window_s       = 0.10;
    c.cold_start     = ColdStart::MedianFromStart;
    // Widen the weight clamp so prior ratios survive it (the test_split_median.cpp shape).
    c.weight_floor   = 0.001;
    c.weight_cap     = 1000.0;
    // Block calibration COMMITS for these long runs: a mid-run prior re-anchor (scale /
    // extrinsic swap) step-changes every residual-to-consensus track and would conflate the
    // commit transient with the reliability claim under test (the rigs deliberately plant
    // biased/noisy channels that the calibrators would otherwise chase).
    c.commit_min_votes = 1000000;
    return c;
}
} // namespace

// ===========================================================================
// d. COUPLED-PATH PIN: split_median = false must produce EXACTLY the pre-19b values.
// The literals below were captured at HEAD (commit acb94be, BEFORE the per-channel
// change) on the deterministic seeded rig; exact double equality (17 significant
// digits round-trips IEEE-754 bit-exactly). If this test fails after a 19b edit, the
// change LEAKED into the coupled path — fix the leak, never the literals.
// ===========================================================================
namespace {
struct CoupledCapture {
    Scalar rel[4];     // per-source reliability (slots 0..3)
    Scalar bias[4];    // per-source EMA bias
    Scalar weight[4];  // per-source published weight
    Scalar tx, ty, tz; // final frontier translation
    Scalar r01, r12;   // two off-diagonal rotation entries (orientation pin)
    Scalar c0, c7;     // two covariance diagonals (trans + rot blocks)
};
CoupledCapture run_coupled_pin_rig() {
    Trajectory tr = mixed_traj();
    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[2].noise_trans_per_m = 0.08;
    planted[2].noise_rot_per_rad = 0.08;
    planted[2].noise_trans_floor = 0.01;
    planted[2].noise_rot_floor   = 0.01;
    planted[2].seed              = 42u;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    // Non-default split-side knobs that must stay INERT while split_median is false.
    sensors[1].rot_weight_prior = 10.0;
    Config cfg = base_cfg(sensors, 4);
    cfg.split_median = false;            // THE PIN: the coupled path
    cfg.weight_cap   = 1e9;
    cfg.weight_floor = 1e-6;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    const Result& res = rig.estimator().latest();
    CoupledCapture cap{};
    for (int i = 0; i < 4; ++i) {
        const SourceHealth* h = health(res, static_cast<SourceId>(i));
        REQUIRE(h != nullptr);
        cap.rel[i]    = h->reliability;
        cap.bias[i]   = h->bias;
        cap.weight[i] = h->weight;
    }
    cap.tx  = res.frontier.pose.t.x();
    cap.ty  = res.frontier.pose.t.y();
    cap.tz  = res.frontier.pose.t.z();
    cap.r01 = res.frontier.pose.R(0, 1);
    cap.r12 = res.frontier.pose.R(1, 2);
    cap.c0  = res.frontier.cov(0, 0);
    cap.c7  = res.frontier.cov(7, 7);
    return cap;
}
} // namespace

TEST_CASE("19b coupled pin: split_median=false is EXACTLY the pre-change behavior "
          "(HEAD-captured literals; a failure here = the 19b change leaked)") {
    const CoupledCapture cap = run_coupled_pin_rig();
    {
        std::ostringstream os;
        os << std::setprecision(17);
        os << "CAPTURE rel={" << cap.rel[0] << ", " << cap.rel[1] << ", " << cap.rel[2]
           << ", " << cap.rel[3] << "}\n bias={" << cap.bias[0] << ", " << cap.bias[1]
           << ", " << cap.bias[2] << ", " << cap.bias[3] << "}\n w={" << cap.weight[0]
           << ", " << cap.weight[1] << ", " << cap.weight[2] << ", " << cap.weight[3]
           << "}\n t={" << cap.tx << ", " << cap.ty << ", " << cap.tz << "}\n R01="
           << cap.r01 << " R12=" << cap.r12 << " c0=" << cap.c0 << " c7=" << cap.c7;
        MESSAGE(os.str());
    }
    // HEAD-captured expected values (EXACT equality — see the header comment).
    const Scalar exp_rel[4]  = {1.0, 1.0, 0.20000000000000001, 1.0};
    const Scalar exp_bias[4] = {1.0850258221605425e-09, 5.4230971732095468e-10,
                                0.024267219281089773,   5.4230971732095468e-10};
    const Scalar exp_w[4]    = {999000999.00099885, 999000999.00099885,
                                1999.9800001999979, 999000999.00099885};
    for (int i = 0; i < 4; ++i) {
        CHECK(cap.rel[i]    == exp_rel[i]);
        CHECK(cap.bias[i]   == exp_bias[i]);
        CHECK(cap.weight[i] == exp_w[i]);
    }
    CHECK(cap.tx  == Scalar(4.9755902487059549));
    CHECK(cap.ty  == Scalar(0.84798055055101496));
    CHECK(cap.tz  == Scalar(-0.96278390432460448));
    CHECK(cap.r01 == Scalar(-0.67006230862515093));
    CHECK(cap.r12 == Scalar(-0.2446423125561257));
    CHECK(cap.c0  == Scalar(0.087509427637821685));
    CHECK(cap.c7  == Scalar(0.0025003141274100348));
}

// ===========================================================================
// d (continued). Coupled path: the NEW per-channel SourceHealth fields stay at their
// DEFAULTS (1, 1, 0, 0) when split_median is false — the documented coupled-path choice
// (the per-channel tracks only exist on the split path; the scalar track is authoritative).
// ===========================================================================
TEST_CASE("19b coupled defaults: per-channel SourceHealth fields stay (1,1,0,0) on the "
          "coupled path") {
    Trajectory tr = mixed_traj();
    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[2].noise_rot_per_rad = 0.08;
    planted[2].noise_rot_floor   = 0.01;
    planted[2].seed              = 42u;
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_cfg(sensors, 4);
    cfg.split_median = false;                       // COUPLED
    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);                           // well past warmup: scalar track active
    const Result& res = rig.estimator().latest();
    for (int i = 0; i < 4; ++i) {
        const SourceHealth* h = health(res, static_cast<SourceId>(i));
        REQUIRE(h != nullptr);
        CHECK(h->reliability_rot   == Scalar(1));   // defaults — never written when coupled
        CHECK(h->reliability_trans == Scalar(1));
        CHECK(h->bias_rot          == Scalar(0));
        CHECK(h->bias_trans        == Scalar(0));
    }
}

// ===========================================================================
// a. HEADLINE (fields): under split_median a source noisy ONLY in rotation reads
// reliability_rot floored with reliability_trans ~ 1 — and the symmetric trans-noisy
// source reads the mirror. The clean sources stay ~1 in BOTH channels. The legacy
// scalar fields mirror the TRANSLATION channel (back-compat contract).
// ===========================================================================
TEST_CASE("19b split headline (fields): rot-noisy source downweighted ONLY in rotation; "
          "trans-noisy source ONLY in translation") {
    Trajectory tr = mixed_traj();
    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    // Source 2: ROTATION noise only (translation exact). Source 3: TRANSLATION noise only.
    planted[2].noise_rot_per_rad = 0.08;
    planted[2].noise_rot_floor   = 0.01;
    planted[2].seed              = 42u;
    planted[3].noise_trans_per_m = 0.08;
    planted[3].noise_trans_floor = 0.01;
    planted[3].seed              = 7u;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_cfg(sensors, 4);
    cfg.split_median = true;       // the path under test
    cfg.split_veto   = false;      // isolate layer (b): the veto is its own pinned policy
                                   // (it would scale the noisy sources' OTHER-channel
                                   // weights, which is not the per-channel-residual claim)
    cfg.commit_min_votes = 1000000;   // no mid-run prior re-anchors (see split_base_config)

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    const Result& res = rig.estimator().latest();
    const SourceHealth* h0 = health(res, 0);
    const SourceHealth* h1 = health(res, 1);
    const SourceHealth* h2 = health(res, 2);   // rot-noisy
    const SourceHealth* h3 = health(res, 3);   // trans-noisy
    REQUIRE(h0 != nullptr); REQUIRE(h1 != nullptr);
    REQUIRE(h2 != nullptr); REQUIRE(h3 != nullptr);
    MESSAGE("h2 rel_rot=" << h2->reliability_rot << " rel_trans=" << h2->reliability_trans
            << " | h3 rel_rot=" << h3->reliability_rot
            << " rel_trans=" << h3->reliability_trans);

    // The rot-noisy source: downweighted in ROTATION ONLY.
    const Scalar clean_rr = std::min(h0->reliability_rot, h1->reliability_rot);
    const Scalar clean_rt = std::min(h0->reliability_trans, h1->reliability_trans);
    CHECK(h2->reliability_rot < Scalar(0.8) * clean_rr);
    CHECK(near_abs(h2->reliability_trans, Scalar(1.0), Scalar(0.3)));   // trans KEPT
    // The trans-noisy source: downweighted in TRANSLATION ONLY (the symmetric case).
    CHECK(h3->reliability_trans < Scalar(0.8) * clean_rt);
    CHECK(near_abs(h3->reliability_rot, Scalar(1.0), Scalar(0.3)));     // rot KEPT
    // Clean sources ~1 in BOTH channels; nothing collapses below the floor.
    for (const SourceHealth* h : {h0, h1}) {
        CHECK(near_abs(h->reliability_rot,   Scalar(1.0), Scalar(0.3)));
        CHECK(near_abs(h->reliability_trans, Scalar(1.0), Scalar(0.3)));
    }
    for (const SourceHealth* h : {h0, h1, h2, h3}) {
        CHECK(h->reliability_rot   >= cfg.reliability_floor - Scalar(1e-9));
        CHECK(h->reliability_trans >= cfg.reliability_floor - Scalar(1e-9));
        // Legacy scalar fields mirror the TRANSLATION channel under split (back-compat).
        CHECK(h->reliability == h->reliability_trans);
        CHECK(h->bias        == h->bias_trans);
    }
}

// ===========================================================================
// a'. HEADLINE (effect) + e. rot_weight_prior composition. Sources 0, 1: clean rotation
// (yaw 0) but translation biased +30% (2.6 m/s). Source 2: the UNIQUELY CORRECT 2.0 m/s
// translation (weight_prior 2.5 — its translation mass carries the channel) with a
// zero-mean VARYING-magnitude yaw-noise channel.
//   * rot_weight_prior 0.7 (rot mass 1.75 < the clean pair's 2): the rotation channel
//     rejects the noise from step one, the per-channel track FLOORS reliability_rot —
//     and reliability_trans stays ~1, so the fused translation keeps tracking source 2's
//     CORRECT translation through the whole run. THE MUTATION THIS SLICE KILLS: the old
//     MIXED residual would have floored the scalar reliability (its scatter is the rot
//     noise) and dropped source 2's translation mass to 0.625 < 2 — losing the only
//     correct translation (x error ~2 m instead of ~0).
//   * rot_weight_prior 1.2 (rot mass 3.0 > 2.5 = half): the datasheet prior OVERRIDES the
//     channel — the rotation median rides source 2 (and its self-consistent residual
//     CAPS reliability_rot) — pinning that the prior MULTIPLIES the learned per-channel
//     reliability into the same channel weight rather than replacing or bypassing it.
// ===========================================================================
TEST_CASE("19b split headline (effect): rot-noisy source's translation still drives the "
          "fused translation; rot_weight_prior composes with reliability_rot") {
    Vec6 xi_biased; xi_biased << 2.6, 0, 0, 0, 0, 0;   // sources 0, 1 (translation +30%)
    Vec6 xi_good;   xi_good   << 2.0, 0, 0, 0, 0, 0;   // source 2 base (correct trans)
    const double kSpan = 6.0, kRate = 200.0, kRunTo = 5.5, kAmp = 0.8;

    struct RunOut {
        Scalar early_yaw;   // mean |window yaw| over fused windows 5..15
        Scalar late_yaw;    // mean |window yaw| over the last 15 windows (post-warmup)
        Scalar x_err;       // final |x - 2.0 t| (the correct-translation tracking error)
        Scalar rel_rot, rel_trans, rel_legacy;          // source 2's published fields
    };
    auto run = [&](Scalar rw2) -> RunOut {
        std::vector<std::unique_ptr<BufferSource>> srcs;
        for (int i = 0; i < 3; ++i)
            srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 4096));
        fill_twist(*srcs[0], kSpan, kRate, xi_biased);
        fill_twist(*srcs[1], kSpan, kRate, xi_biased);
        fill_twist_noisy_yaw(*srcs[2], kSpan, kRate, xi_good, kAmp);

        SensorConfig sensors[3];
        for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
        sensors[2].weight_prior     = 2.5;
        sensors[2].rot_weight_prior = rw2;
        Config cfg = split_base_config(3);
        cfg.split_median = true;
        cfg.split_veto   = false;   // per-channel grace, not whole-source veto (the
                                    // Slice-19 headline's policy fork, same rationale)
        cfg.sensors      = sensors;
        cfg.sensor_count = 3;

        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);

        std::vector<Scalar> wyaw;   // |window yaw| per fused step
        Mat3 R_prev = Mat3::Identity();
        double now_s = cfg.fusion_delay_s + cfg.window_s;
        double last_t1 = 0.0;
        while (now_s <= kRunTo) {
            REQUIRE(est.step(secs(now_s)) == Status::Ok);
            const Mat3 R_now = est.latest().frontier.pose.R;
            wyaw.push_back(so3::log(R_prev.transpose() * R_now).norm());
            R_prev  = R_now;
            last_t1 = now_s - cfg.fusion_delay_s;
            now_s += cfg.window_s;
        }
        REQUIRE(wyaw.size() > 40);   // past warmup (kRelWarmup = 20) with margin

        RunOut out{};
        out.early_yaw = 0; out.late_yaw = 0;
        for (std::size_t k = 5; k < 15; ++k) out.early_yaw += wyaw[k];
        out.early_yaw /= Scalar(10);
        for (std::size_t k = wyaw.size() - 15; k < wyaw.size(); ++k) out.late_yaw += wyaw[k];
        out.late_yaw /= Scalar(15);
        out.x_err = std::abs(est.latest().frontier.pose.t.x() - 2.0 * last_t1);
        const SourceHealth* h2 = health(est.latest(), 2);
        REQUIRE(h2 != nullptr);
        out.rel_rot    = h2->reliability_rot;
        out.rel_trans  = h2->reliability_trans;
        out.rel_legacy = h2->reliability;
        return out;
    };

    // --- rot prior 0.7: rotation rejected, translation KEPT (the layer-b headline).
    const RunOut a = run(0.7);
    MESSAGE("rw=0.7: early_yaw=" << a.early_yaw << " late_yaw=" << a.late_yaw
            << " x_err=" << a.x_err << " rel_rot=" << a.rel_rot
            << " rel_trans=" << a.rel_trans);
    CHECK(a.late_yaw < 0.005);                  // rotation tracks the clean pair
    // The noisy-rot source's translation WON. The bound is NOT ~0: the eps-regularized
    // Weiszfeld approaches the majority vertex but keeps a small blend remainder toward
    // the coincident minority pair (~0.005 m/window, ~0.26 m over the run). The MUTATION
    // case (mixed residual flooring the source's translation weight too) loses the only
    // correct translation entirely: x_err ~ 0.6 m/s x ~3.5 s ~ 2 m. 0.4 separates cleanly.
    CHECK(a.x_err   < 0.4);
    CHECK(a.rel_rot   <= Scalar(0.25));         // rotation channel floored (0.2)
    CHECK(a.rel_trans >= Scalar(0.75));         // translation channel KEPT (~1)
    CHECK(a.rel_legacy == a.rel_trans);         // legacy mirrors the translation channel

    // --- rot prior 1.2: the prior MULTIPLIES through — rotation mass 2.5 x 1.2 = 3.0
    //     holds the channel (and the self-consistent residual caps reliability_rot), while
    //     the translation channel is unchanged (still correct).
    const RunOut b = run(1.2);
    MESSAGE("rw=1.2: late_yaw=" << b.late_yaw << " x_err=" << b.x_err
            << " rel_rot=" << b.rel_rot);
    CHECK(b.late_yaw > 0.02);                   // rotation rides the prior-boosted source
    CHECK(b.rel_rot  > Scalar(2.0));            // its own-consensus residual reads quiet
    CHECK(b.x_err    < 0.4);                    // translation channel unaffected either way
}

// ===========================================================================
// b. D17 PER CHANNEL: a source with a SYSTEMATIC rotation bias but LOW rotation scatter
// is NOT rot-downweighted (its offset lands in bias_rot, not the weight — the calibrator's
// job); the genuinely rot-NOISY source in the same rig IS. All four sources share the
// translation exactly, so every reliability_trans stays ~1.
// ===========================================================================
TEST_CASE("19b D17 per channel: systematic rotation bias is NOT rot-downweighted; "
          "rotation scatter IS") {
    // ZERO linear velocity: with v != 0 a yaw-rate difference curves the integrated window
    // translation (t_y ~ v*w*dt^2/2), so the planted ROTATION noise would leak into the
    // translation channel and cloud the "reliability_trans stays ~1" claim. v = 0 keeps
    // every source's window translation IDENTICALLY zero — the rotation channel is the
    // only one perturbed, by construction.
    Vec6 xi_base; xi_base << 0, 0, 0, 0, 0, 0;
    const double kSpan = 5.5, kRate = 200.0, kRunTo = 5.0;

    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < 4; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 4096));
    // A small COMMON yaw scatter on every source (different table phases) makes the
    // per-channel variance baseline well-posed (the test_weights.cpp biased-test pattern);
    // source 2 = 10x scatter (the NOISY one); source 3 = +0.4 rad/s CONSTANT yaw offset
    // on top of the common scatter (the BIASED one: large mean, comparable variance).
    fill_twist_noisy_yaw(*srcs[0], kSpan, kRate, xi_base, 0.1, 0.0, 0);
    fill_twist_noisy_yaw(*srcs[1], kSpan, kRate, xi_base, 0.1, 0.0, 3);
    fill_twist_noisy_yaw(*srcs[2], kSpan, kRate, xi_base, 1.0, 0.0, 5);
    fill_twist_noisy_yaw(*srcs[3], kSpan, kRate, xi_base, 0.1, 0.4, 6);

    SensorConfig sensors[4];
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = split_base_config(4);
    cfg.split_median = true;
    cfg.split_veto   = false;   // the biased source's constant 0.04 rad deviation would
                                // trip the cross-channel veto (its own policy, pinned in
                                // test_split_median.cpp) and cloud the reliability claim
    cfg.sensors      = sensors;
    cfg.sensor_count = 4;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);
    double now_s = cfg.fusion_delay_s + cfg.window_s;
    int steps = 0;
    while (now_s <= kRunTo) {
        REQUIRE(est.step(secs(now_s)) == Status::Ok);
        ++steps;
        now_s += cfg.window_s;
    }
    REQUIRE(steps > 40);   // past warmup with margin

    const Result& res = est.latest();
    const SourceHealth* h0 = health(res, 0);
    const SourceHealth* h1 = health(res, 1);
    const SourceHealth* h2 = health(res, 2);   // rot-NOISY
    const SourceHealth* h3 = health(res, 3);   // rot-BIASED
    REQUIRE(h0 != nullptr); REQUIRE(h1 != nullptr);
    REQUIRE(h2 != nullptr); REQUIRE(h3 != nullptr);
    MESSAGE("rel_rot: clean={" << h0->reliability_rot << ", " << h1->reliability_rot
            << "} noisy=" << h2->reliability_rot << " biased=" << h3->reliability_rot
            << " | bias_rot: clean={" << h0->bias_rot << ", " << h1->bias_rot
            << "} biased=" << h3->bias_rot);

    // The BIASED source is KEPT in the rotation channel (D17, per channel): its constant
    // offset lives in bias_rot, not the variance. NOTE the clean sources read ABOVE 1
    // (up to the cap): with q = 4 warmed participants the reference variance is the
    // lower-MIDDLE (3rd-smallest), so the quieter-than-reference sources read > 1 — the
    // assertions are therefore DIRECTIONAL (kept vs floored), not "== 1".
    CHECK(h3->reliability_rot >= Scalar(0.75));                          // kept (~1)
    CHECK(h3->reliability_rot >= cfg.reliability_floor + Scalar(0.3));   // far from collapse
    // ...and its bias_rot diagnostic carries the systematic excess.
    CHECK(h3->bias_rot > Scalar(2.0) * std::max(h0->bias_rot, h1->bias_rot));
    // The NOISY source IS rot-downweighted (the contrast in the same rig): well below the
    // kept biased source and headed for the floor.
    CHECK(h2->reliability_rot <= Scalar(0.45));
    CHECK(h2->reliability_rot < Scalar(0.6) * h3->reliability_rot);
    // The clean sources are never floored either.
    CHECK(h0->reliability_rot >= Scalar(0.75));
    CHECK(h1->reliability_rot >= Scalar(0.75));
    // Identical (zero) translations: every reliability_trans stays ~1.
    for (const SourceHealth* h : {h0, h1, h2, h3})
        CHECK(near_abs(h->reliability_trans, Scalar(1.0), Scalar(0.25)));
}

// ===========================================================================
// c. Warmup: per-channel reliabilities hold EXACTLY 1.0 before kRelWarmup (= 20) samples —
// a short split run never perturbs the Slice-2 weights.
// ===========================================================================
TEST_CASE("19b warmup: per-channel reliabilities stay exactly 1.0 before kRelWarmup") {
    Trajectory tr = mixed_traj();
    std::vector<SourceParams> planted(3);
    for (int i = 0; i < 3; ++i) planted[i].id = static_cast<SourceId>(i);
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(3);
    for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
    Config cfg = base_cfg(sensors, 3);
    cfg.split_median = true;                    // split path (veto default ON is fine here)

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, 0.4, 50.0);  // ~10 fused ticks, well below kRelWarmup
    REQUIRE(fuses > 2);
    REQUIRE(fuses < 20);

    const Result& res = rig.estimator().latest();
    for (int i = 0; i < 3; ++i) {
        const SourceHealth* h = health(res, static_cast<SourceId>(i));
        REQUIRE(h != nullptr);
        CHECK(h->reliability_rot   == Scalar(1));   // exactly neutral pre-warmup
        CHECK(h->reliability_trans == Scalar(1));
        CHECK(h->reliability       == Scalar(1));   // legacy mirror included
    }
}
