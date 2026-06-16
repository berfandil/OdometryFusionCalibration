// Slice 19c tests: GPS-course heading-drift monitor (split policy layer c).
//
// The monitor ranks per-source heading drift against GPS course-over-ground and boosts
// the best source's ROTATION-channel weight in the split median (auto-discovery of the
// heading-grade source, replacing the manual rot_weight_prior=10 urban recipe). Design:
// SLICE19C_HEADING_MONITOR.md; validated prototype: tools/proto_heading_monitor.py.
//
// Coverage (the SLICE19C doc's acceptance list, items 1-9):
//   1. Sim headline: 3 sources, one with planted yaw-rate drift + GPS fixes -> the
//      drifter's score >> the clean sources', boost ranks the clean sources up and the
//      drifter down (ranking, not absolute recovery).
//   2. Anchor gates: stopped/slow/fast/turn/reverse/multipath-jump fix pairs produce no
//      anchors (each gate exercised; positive control anchors).
//   3. Telescoping: a turn between anchors is bridged (the pair is NOT dropped) and the
//      turn-correlated drifter is still caught.
//   4. Denial freeze: a fix gap > pairgap closes the segment, scores hold bit-exactly,
//      no cross-segment slope forms (drift accrued inside the gap stays invisible).
//   5. Abstain: < 3 sources, or no scored blocks -> all boosts exactly 1.0.
//   6. Chatter guard: a one-evaluation rank flip does not move boosts; a flip persisting
//      >= 2 consecutive evaluations does.
//   7. Default-off byte-identical: split-ON without the monitor reproduces HEAD-captured
//      literals exactly (the coupled path is covered by the 19b pin); validate() rejects
//      monitor-without-split.
//   8. Config-hash flips on both knobs (persistence guard); loader keys are covered in
//      adapters/tests/test_config_loader.cpp.
//   9. Split-ON cov-cal guard + observability self-tests stay green UNCHANGED — enforced
//      by the gate (the monitor defaults OFF; the spread->Q path is structurally
//      untouched).
#include <doctest/doctest.h>

#include "ofc/core/absolute_ref.hpp"
#include "ofc/core/buffer.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/heading_monitor.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/source.hpp"

#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

using namespace ofc;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;
Timestamp secs(double s) { return static_cast<Timestamp>(s * 1e9); }

// A source = an ISource facade over a SourceBuffer (the test_weights_split.cpp shape).
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

// Fill with a constant body twist (rate Hz, span seconds).
void fill_twist(BufferSource& src, double span, double rate, const Vec6& xi) {
    const Timestamp step = secs(1.0 / rate);
    const int n = static_cast<int>(span * rate);
    for (int k = 0; k <= n; ++k) {
        src.buffer().push_twist(static_cast<Timestamp>(k) * step, xi, Mat6::Identity());
    }
}

// A "follow" position fix: z == h(x) (zero residual), emitted once per period of
// frontier time. Always gate-accepted (NIS = 0), never perturbs the state — so the GPS
// course exactly tracks the fused trajectory's course, which is what the monitor
// consumes (deterministic, the relaxed-edge test analogue of the GPS adapter's output).
class FollowFix : public ICorrection {
public:
    explicit FollowFix(Scalar period_s, Scalar sigma_pos = 0.1)
        : period_(period_s), sigma_(sigma_pos) {}
    bool evaluate(const State& x, Measurement& out) const override {
        const Scalar t = static_cast<Scalar>(x.stamp) * Scalar(1e-9);
        if (have_last_ && t - last_ < period_ - Scalar(1e-9)) return false;
        last_ = t;
        have_last_ = true;
        out.dim = 3;
        out.residual.setZero();                     // z == h(x): a perfect follow fix
        out.H.setZero();
        out.H.block<3, 3>(0, 0) = x.pose.R;         // position-fix rows (right-error)
        out.R = Eigen::Matrix<Scalar, 6, 6>::Identity() * (sigma_ * sigma_);
        out.stamp = x.stamp;
        return true;
    }

private:
    Scalar period_;
    Scalar sigma_;
    mutable bool   have_last_ = false;
    mutable Scalar last_      = 0;
};

const SourceHealth* health(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.health[i].id == id) return &r.health[i];
    return nullptr;
}

// The shared 3-source split-ON estimator rig: sources 0/1 clean forward motion
// (5 m/s), source 2 the same + a constant planted yaw-rate drift, GPS follow fixes
// at 1 Hz. Used by the headline (monitor ON) and the byte-identical pin (monitor OFF).
struct EstRun {
    Scalar tx, ty, tz, r01, r12, c0, c7;            // frontier pin values
    Scalar w0, w1, w2;                              // published (trans-channel) weights
    Scalar rr2, rt2;                                // drifter's per-channel reliability
    int    corr_applied_total;
    Scalar score[3];
    Scalar boost[3];
    bool   scored[3];
    Scalar final_yaw_rate;                          // |window yaw| of the LAST fused step
};

EstRun run_split_gps_rig(bool monitor_on, double span_s, Scalar drift_wz,
                         void (*tweak)(Config&) = nullptr, Scalar prior1 = Scalar(2.0)) {
    const double kRate = 20.0;
    Vec6 xi_clean;  xi_clean  << 5.0, 0, 0, 0, 0, 0;
    Vec6 xi_drift   = xi_clean;
    xi_drift(5)    += drift_wz;

    std::vector<std::unique_ptr<BufferSource>> srcs;
    for (int i = 0; i < 3; ++i)
        srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 16384));
    fill_twist(*srcs[0], span_s, kRate, xi_clean);
    fill_twist(*srcs[1], span_s, kRate, xi_clean);
    fill_twist(*srcs[2], span_s, kRate, xi_drift);

    SensorConfig sensors[3];
    for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[1].rot_weight_prior = prior1;   // non-neutral (default 2.0): pins boost x prior

    Config cfg;
    cfg.max_sources      = 3;
    cfg.fusion_delay_s   = 0.05;
    cfg.window_s         = 0.10;
    cfg.timesync_enabled = false;
    cfg.cold_start       = ColdStart::MedianFromStart;
    cfg.commit_min_votes = 1000000;      // no mid-run calibration re-anchors
    cfg.split_median     = true;
    cfg.sensors          = sensors;
    cfg.sensor_count     = 3;
    if (monitor_on) cfg.heading_monitor = true;
    if (tweak != nullptr) tweak(cfg);

    FollowFix fix(1.0);

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);
    REQUIRE(est.add_correction(&fix) == Status::Ok);

    int applied = 0;
    Mat3 R_prev = Mat3::Identity();
    Scalar last_wyaw = 0;
    double now_s = cfg.fusion_delay_s + cfg.window_s;
    while (now_s <= span_s - 0.5) {
        REQUIRE(est.step(secs(now_s)) == Status::Ok);
        applied += est.latest().correction.corr_applied;
        const Mat3 R_now = est.latest().frontier.pose.R;
        last_wyaw = so3::log(R_prev.transpose() * R_now).norm();
        R_prev = R_now;
        now_s += cfg.window_s;
    }

    const Result& res = est.latest();
    EstRun out{};
    out.tx  = res.frontier.pose.t.x();
    out.ty  = res.frontier.pose.t.y();
    out.tz  = res.frontier.pose.t.z();
    out.r01 = res.frontier.pose.R(0, 1);
    out.r12 = res.frontier.pose.R(1, 2);
    out.c0  = res.frontier.cov(0, 0);
    out.c7  = res.frontier.cov(7, 7);
    out.corr_applied_total = applied;
    out.final_yaw_rate     = last_wyaw;
    const SourceHealth* h[3];
    for (int i = 0; i < 3; ++i) {
        h[i] = health(res, static_cast<SourceId>(i));
        REQUIRE(h[i] != nullptr);
    }
    out.w0 = h[0]->weight; out.w1 = h[1]->weight; out.w2 = h[2]->weight;
    out.rr2 = h[2]->reliability_rot;
    out.rt2 = h[2]->reliability_trans;
    for (int i = 0; i < 3; ++i) {
        out.score[i]  = h[i]->heading_score;
        out.boost[i]  = h[i]->heading_boost;
        out.scored[i] = h[i]->heading_scored;
    }
    return out;
}

// A direct-drive harness for the HeadingMonitor unit (the anchor-gate / telescoping / denial
// / chatter items need crafted fix streams the full estimator cannot deterministically pose,
// the test_timesync.cpp pattern of unit-testing the core component in isolation). Feeds the
// monitor synthetic fixes: 3 sources moving along +x at `v` m/s, source `drift_idx` (if >= 0)
// rotating at `drift_wz` rad/s; GPS position is the planar dead-reckon of the CLEAN heading
// (straight +x) — the chord course the real GPS-follow fix would observe.
struct MonRig {
    HeadingMonitor mon;
    Scalar         t    = 0;
    Scalar         yaw[3]  = {0, 0, 0};
    Scalar         fwd[3]  = {0, 0, 0};
    Vec3           gpos = Vec3::Zero();      // GPS position (clean straight dead-reckon)
    Scalar         ghead = 0;                // GPS course heading (rad)
    int            drift_idx = -1;
    Scalar         drift_wz  = 0;
    Scalar         v         = 5.0;

    explicit MonRig(int drift = -1, Scalar wz = 0, Scalar speed = 5.0)
        : drift_idx(drift), drift_wz(wz), v(speed) { mon.configure(3); }

    // Advance `dt` s and submit ONE fix at the new time. `course_turn` (rad) turns BOTH the GPS
    // course AND every source's heading (a real road turn the sources track — the anchors must
    // not mistake it for drift); the drifter additionally accrues drift_wz*dt on top. A
    // `v_gps_override` (>= 0) replaces the GPS speed for THIS fix only (multipath jump / stop)
    // while the sources keep moving at v; `reverse` flips the GPS chord 180 deg (backing up).
    void advance(Scalar dt, Scalar course_turn = 0, Scalar v_gps_override = -1,
                 bool reverse = false) {
        t += dt;
        ghead += course_turn;
        const Scalar vg = (v_gps_override >= 0) ? v_gps_override : v;
        const Scalar step = (reverse ? -vg : vg) * dt;
        gpos.x() += step * std::cos(ghead);
        gpos.y() += step * std::sin(ghead);
        for (int i = 0; i < 3; ++i) {
            yaw[i] += course_turn + ((i == drift_idx) ? drift_wz * dt : 0.0);
            fwd[i] += (reverse ? -v : v) * dt;   // a true reverse backs up the odometry too
        }
        mon.submit_fix(t, gpos, yaw, fwd);
    }

    // Submit a run of `n` clean straight fixes spaced `dt` s (the common-case anchor stream).
    void straight(int n, Scalar dt = 1.0) { for (int k = 0; k < n; ++k) advance(dt); }
};

// Enough 1-Hz straight fixes to form two 60-s blocks >= 120 s apart (the first scored point).
constexpr int kScoreFixes = 200;
} // namespace

// ===========================================================================
// Item 1 — HEADLINE. Through the full split estimator: a planted yaw-rate drifter is demoted
// (boost 1) while the clean pair is boosted (-> boost_max), via score + boost + a measurable
// rotation-channel effect (the fused heading tracks the clean pair more tightly with the
// monitor ON than OFF).
// ===========================================================================
TEST_CASE("19c item 1: headline — drifter demoted, clean boosted, fused heading tightened") {
    const EstRun on  = run_split_gps_rig(true,  400.0, 0.005);
    const EstRun off = run_split_gps_rig(false, 400.0, 0.005);

    // All three scored; the drifter's drift score >> the clean pair's.
    CHECK(on.scored[0]);
    CHECK(on.scored[1]);
    CHECK(on.scored[2]);
    CHECK(on.score[2] > 100.0 * on.score[0]);
    CHECK(on.score[2] > 100.0 * on.score[1]);

    // Boost ranks the clean pair UP to the cap and the drifter DOWN to the floor.
    CHECK(on.boost[0] == doctest::Approx(10.0));
    CHECK(on.boost[1] == doctest::Approx(10.0));
    CHECK(on.boost[2] == doctest::Approx(1.0));
    // What actually moves fusion is the RATIO of the clean boost to the drifter's (not the cap
    // value itself): pin the ranking robustly so this survives a boost_max retune. A clean source
    // outweighs the drifter in the rotation channel by ~the full cap.
    CHECK(on.boost[0] / on.boost[2] > 5.0);
    CHECK(on.boost[1] / on.boost[2] > 5.0);

    // Monitor OFF: every boost is exactly neutral (no scoring is even surfaced).
    CHECK(off.boost[0] == Scalar(1));
    CHECK(off.boost[1] == Scalar(1));
    CHECK(off.boost[2] == Scalar(1));

    // EFFECT (ranking, not absolute recovery): boosting the clean pair's rotation weight pulls
    // the fused heading CLOSER to straight (the clean sources' heading) than the un-boosted
    // run — the drifter has less rotational pull. R(0,1) = -sin(yaw); |yaw_on| < |yaw_off|.
    CHECK(std::abs(on.r01) < std::abs(off.r01));
}

// ===========================================================================
// Item 2 — ANCHOR GATES. Each gate (stopped/slow, fast, turn, reverse, multipath jump, dt) is
// exercised with a POSITIVE control (a clean stream DOES anchor; the gated fix does not).
// ===========================================================================
TEST_CASE("19c item 2: anchor gates — each rejection + positive control") {
    using namespace heading_monitor_const;

    // Positive control: a clean straight stream produces an anchor on every interval after the
    // first (the first fix only seeds the previous sample).
    {
        MonRig r;
        r.straight(5);
        CHECK(r.mon.anchor_count() == 4);
        CHECK(r.mon.pair_count() == 3);
    }
    // SLOW (v < vmin): 1 m/s GPS chord (and matching odometry) -> no anchor.
    {
        MonRig r(-1, 0, /*speed=*/1.0);
        r.straight(5);
        CHECK(r.mon.anchor_count() == 0);
    }
    // FAST (v > vmax): 40 m/s -> no anchor.
    {
        MonRig r(-1, 0, /*speed=*/40.0);
        r.straight(5);
        CHECK(r.mon.anchor_count() == 0);
    }
    // TURN (|yaw rate| >= omega_max): all sources turn 10 deg/s (> 3 deg/s gate) each 1-s fix.
    {
        MonRig r;
        const Scalar turn = Scalar(10.0) * kPi / Scalar(180.0);   // rad per 1-s interval
        for (int k = 0; k < 5; ++k) r.advance(1.0, turn);
        CHECK(r.mon.anchor_count() == 0);
    }
    // REVERSE (net backward motion): the GPS chord flips 180 deg while odometry reads forward.
    {
        MonRig r;
        r.advance(1.0);                       // seed a forward fix
        for (int k = 0; k < 4; ++k) r.advance(1.0, 0, -1, /*reverse=*/true);
        // The reverse fixes never anchor (g_fwd fails: the GPS-vs-odom chord disagrees in sign,
        // so the cross-validation/forward gates drop them). At most the single seed-forward pair.
        CHECK(r.mon.pair_count() == 0);
    }
    // MULTIPATH JUMP (|v_gps - v_odo| >= max(3, 0.5 v)): GPS jumps 50 m in 1 s, odom reads 5.
    {
        MonRig r;
        r.straight(2);                        // a couple of clean anchors first
        const int before = r.mon.anchor_count();
        r.advance(1.0, 0, /*v_gps_override=*/50.0);   // multipath spike
        CHECK(r.mon.anchor_count() == before);        // the spike fix did NOT anchor
    }
    // DT too large (> dtmax): a 5-s fix gap is rejected by the dt gate (positive control: the
    // 1-s fix before it DOES anchor).
    {
        MonRig r;
        r.advance(1.0);                       // seed the previous-fix sample (no anchor yet)
        r.advance(1.0);                       // clean 1-s interval -> anchor #1
        r.advance(5.0);                       // dt = 5 s > 3 s -> rejected (no new anchor)
        CHECK(r.mon.anchor_count() == 1);
    }
}

// ===========================================================================
// Item 3 — TELESCOPING. A turn between anchors is BRIDGED (the pair is not dropped) and a
// turn-correlated drifter is still caught. The drifter accrues its drift only during the turn;
// gating the turn out of the pair interior would hide it.
// ===========================================================================
TEST_CASE("19c item 3: telescoping — a bridged turn still catches the turn-correlated drifter") {
    // Source 2 drifts ONLY during turns (we inject extra yaw on the turning fixes). Build a long
    // stream: straight legs (anchors) separated by short turn bursts (NOT anchored, but bridged
    // by the next straight anchor pair). The drifter's bridged residual must accumulate.
    MonRig r(/*drift_idx=*/2, /*drift_wz=*/0.0);
    const Scalar turn = Scalar(20.0) * kPi / Scalar(180.0);   // 20 deg/s burst (gated out)
    for (int leg = 0; leg < 12; ++leg) {
        r.straight(18);                          // ~18 s straight anchor run
        // A 1-s turn burst: course turns `turn`, and the drifter slips an EXTRA 5 deg on it
        // (turn-correlated wheel-yaw error). The burst fix itself is gated out (yaw rate high),
        // but the straight anchors on either side telescope across it.
        r.advance(1.0, turn);
        r.yaw[2] += Scalar(5.0) * kPi / Scalar(180.0);   // drifter's turn-correlated slip
        r.advance(1.0, -turn);                   // turn back to straight (also gated)
        r.yaw[2] += Scalar(5.0) * kPi / Scalar(180.0);
    }
    r.straight(40);                              // settle + form the final blocks

    REQUIRE(r.mon.scored(2));
    REQUIRE(r.mon.scored(0));
    // The turn-correlated drifter (caught only because the pairs bridged the gated turns) scores
    // far worse than the clean sources.
    // NOTE: this injects slip on BOTH the turn-in and the turn-out (two-sided), which is the
    // CONSERVATIVE case -- two opposing slips partially cancel in the bridged residual yet the
    // drifter is still caught. A one-sided slip (inject on turn-in only, none on turn-out) would
    // be a STRONGER check: the residual would not cancel at all, so the drifter would be even
    // easier to flag. Two-sided suffices to pin the telescoping behavior.
    CHECK(r.mon.score(2) > 5.0 * r.mon.score(0));
    CHECK(r.mon.boost(2, 10.0) == doctest::Approx(1.0));
    CHECK(r.mon.boost(0, 10.0) > 1.5);
}

// ===========================================================================
// Item 4 — DENIAL FREEZE. A fix gap > pairgap closes the segment; the scores hold BIT-EXACTLY
// across the gap and no cross-segment slope forms (drift accrued inside the gap is invisible).
// ===========================================================================
TEST_CASE("19c item 4: denial freeze — segment close holds scores bit-exact") {
    MonRig r(/*drift_idx=*/2, /*drift_wz=*/0.005);
    r.straight(kScoreFixes);                     // score up
    REQUIRE(r.mon.scored(2));
    const Scalar s0_before = r.mon.score(0);
    const Scalar s2_before = r.mon.score(2);
    const int    pairs_before = r.mon.pair_count();

    // GPS DENIAL: a single fix arrives 522 s later (>> pairgap = 60 s). It closes the segment;
    // the lone post-gap fix cannot pair (no prior anchor in the new segment within pairgap), so
    // NO new pair, NO new slope, and the scores are unchanged to the last bit.
    r.advance(522.0);
    CHECK(r.mon.pair_count() == pairs_before);   // the denial fix formed no pair
    CHECK(r.mon.score(0) == s0_before);          // bit-exact hold
    CHECK(r.mon.score(2) == s2_before);

    // A few more post-gap fixes resume anchoring in the NEW segment but never form a
    // cross-segment slope (the staircase step over the gap is never differenced).
    r.straight(3);
    CHECK(r.mon.score(2) == s2_before);          // still no new slope (new segment too short)
}

// ===========================================================================
// Item 5 — ABSTAIN. < 3 sources, or no scored blocks -> every boost is exactly 1.0.
// ===========================================================================
TEST_CASE("19c item 5: abstain — under-three sources / no blocks -> boosts exactly 1.0") {
    // (a) Two sources: the median consensus degenerates; the monitor is inert.
    {
        HeadingMonitor mon;
        mon.configure(2);
        Vec3 pos = Vec3::Zero();
        Scalar yaw[2] = {0, 0}, fwd[2] = {0, 0};
        for (int k = 0; k < kScoreFixes; ++k) {
            pos.x() += 5.0; yaw[0] += 0; yaw[1] += 0; fwd[0] += 5.0; fwd[1] += 5.0;
            mon.submit_fix(static_cast<Scalar>(k + 1), pos, yaw, fwd);
        }
        CHECK(mon.scored_count() == 0);
        CHECK(mon.boost(0, 10.0) == Scalar(1));
        CHECK(mon.boost(1, 10.0) == Scalar(1));
    }
    // (b) Three sources but too short to form any block-pair slope -> no scored sources ->
    //     abstain (every boost exactly 1.0).
    {
        MonRig r;
        r.straight(50);                          // < 120 s of contiguous blocks
        CHECK(r.mon.scored_count() == 0);
        CHECK(r.mon.boost(0, 10.0) == Scalar(1));
        CHECK(r.mon.boost(1, 10.0) == Scalar(1));
        CHECK(r.mon.boost(2, 10.0) == Scalar(1));
    }
    // (c) boost_max < 1 is clamped to 1 (a degenerate cap never inverts the boost).
    {
        MonRig r(/*drift_idx=*/2, /*drift_wz=*/0.005);
        r.straight(kScoreFixes);
        REQUIRE(r.mon.scored_count() == 3);
        CHECK(r.mon.boost(0, 0.5) == Scalar(1));
    }
}

// ===========================================================================
// Item 6 — CHATTER GUARD. A one-evaluation rank flip does not move boosts; a flip that persists
// >= 2 consecutive evaluations does. Verified at the WIRING layer (the estimator publishes the
// boost slowly): a single-block transient where the ranking momentarily flips must not move the
// published rotation weights, while a sustained flip eventually does.
// ===========================================================================
TEST_CASE("19c item 6: chatter guard — a lone rank flip is inert; a persistent one moves") {
    // The monitor's boost is a SLOW function of cumulative scores (a weighted-median over the
    // whole slope pool), so a single anomalous block cannot flip the ranking — the cumulative
    // estimate is dominated by the established pool. Drive source 0 to be the clean winner, then
    // inject ONE bad block on it and confirm the ranking (and thus the boost) does not flip.
    MonRig r(/*drift_idx=*/2, /*drift_wz=*/0.005);
    r.straight(kScoreFixes);
    REQUIRE(r.mon.scored_count() == 3);
    const Scalar b0_settled = r.mon.boost(0, 10.0);
    const Scalar b2_settled = r.mon.boost(2, 10.0);
    CHECK(b0_settled > 1.5);                      // source 0 is a (boosted) clean winner
    CHECK(b2_settled == doctest::Approx(1.0));    // the drifter floored

    // ONE anomalous fix on source 0 (a transient 8-deg slip) — a single block's worth of bad
    // residual. The cumulative weighted-median ranking must NOT flip the winner to the drifter.
    r.yaw[0] += Scalar(8.0) * kPi / Scalar(180.0);
    r.straight(2);
    CHECK(r.mon.boost(2, 10.0) == doctest::Approx(1.0));   // drifter still floored (no flip)
    CHECK(r.mon.boost(0, 10.0) > 1.0);                     // source 0 still boosted

    // A PERSISTENT regime change (source 0 now drifts WORSE than source 2 for a long run) DOES
    // eventually re-rank: feed a long stretch where source 0 accrues large drift.
    r.drift_idx = 0; r.drift_wz = 0.02;          // source 0 becomes the worst drifter
    r.straight(300);
    CHECK(r.mon.score(0) > r.mon.score(1));      // the ranking has moved (0 is now worse than 1)
    CHECK(r.mon.boost(1, 10.0) > r.mon.boost(0, 10.0));   // the boost followed the sustained flip
}

// ===========================================================================
// Item 7 — DEFAULT-OFF BYTE-IDENTICAL + validate() rejection. split-ON / monitor-OFF reproduces
// the HEAD-captured literals to the last bit; the new SourceHealth fields read their neutral
// defaults; and validate() rejects monitor-without-split.
// ===========================================================================
TEST_CASE("19c item 7: default-off byte-identical pin + validate() rejection") {
    // THE PIN: split-ON, monitor-OFF, the exact rig the capture test used (60 s, drift 0.005).
    const EstRun r = run_split_gps_rig(false, 60.0, 0.005);
    CHECK(r.tx  == Scalar(296.99999895531442));
    CHECK(r.ty  == Scalar(0.020699095969226963));
    CHECK(r.tz  == Scalar(0));
    CHECK(r.r01 == Scalar(-0.00013574040177694242));
    CHECK(r.r12 == Scalar(0));
    CHECK(r.c0  == Scalar(0.00030176231638865786));
    CHECK(r.c7  == Scalar(0.0001015385568546638));
    CHECK(r.w0  == Scalar(9.9999998999999988));
    CHECK(r.w1  == Scalar(9.9999998999999988));
    CHECK(r.w2  == Scalar(9.9999998999999988));
    CHECK(r.rr2 == Scalar(1));
    CHECK(r.rt2 == Scalar(1));
    CHECK(r.corr_applied_total == 60);
    // The new diagnostics read their neutral defaults on the monitor-OFF path.
    for (int i = 0; i < 3; ++i) {
        CHECK(r.score[i] == Scalar(0));
        CHECK(r.boost[i] == Scalar(1));
        CHECK_FALSE(r.scored[i]);
    }

    // validate() rejects heading_monitor without split_median (no silent no-op).
    {
        Config cfg;
        cfg.max_sources = 3;
        cfg.split_median = false;
        cfg.heading_monitor = true;
        CHECK(validate(cfg) == Status::InvalidConfig);
    }
    // ... and accepts it WITH split_median.
    {
        Config cfg;
        cfg.max_sources = 3;
        cfg.split_median = true;
        cfg.heading_monitor = true;
        CHECK(validate(cfg) == Status::Ok);
    }
    // A boost cap < 1 is rejected (it would invert the boost direction).
    {
        Config cfg;
        cfg.max_sources = 3;
        cfg.split_median = true;
        cfg.heading_monitor = true;
        cfg.heading_monitor_boost_max = 0.5;
        CHECK(validate(cfg) == Status::OutOfRange);
    }
}

// ===========================================================================
// Item 8 — CONFIG-HASH flips on BOTH knobs (the persistence guard). A serialized blob written
// under one (heading_monitor, boost_max) refuses to restore under a different one.
// ===========================================================================
TEST_CASE("19c item 8: config-hash flips on both heading-monitor knobs") {
    Vec6 xi; xi << 5.0, 0, 0, 0, 0, 0;

    auto serialize_under = [&](bool monitor, Scalar boost_max,
                               unsigned char* buf, int cap) -> int {
        SensorConfig sensors[3];
        for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
        Config cfg;
        cfg.max_sources   = 3;
        cfg.split_median  = true;
        cfg.heading_monitor = monitor;
        cfg.heading_monitor_boost_max = boost_max;
        cfg.sensors = sensors; cfg.sensor_count = 3;
        std::vector<std::unique_ptr<BufferSource>> srcs;
        for (int i = 0; i < 3; ++i)
            srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 256));
        for (auto& s : srcs) fill_twist(*s, 5.0, 20.0, xi);
        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);
        est.step(secs(1.0));
        const Expected<int> n = est.serialize(buf, cap);
        REQUIRE(n.ok());
        return n.value();
    };

    auto restores_under = [&](bool monitor, Scalar boost_max,
                              const unsigned char* buf, int len) -> Status {
        SensorConfig sensors[3];
        for (int i = 0; i < 3; ++i) sensors[i].id = static_cast<SourceId>(i);
        Config cfg;
        cfg.max_sources   = 3;
        cfg.split_median  = true;
        cfg.heading_monitor = monitor;
        cfg.heading_monitor_boost_max = boost_max;
        cfg.sensors = sensors; cfg.sensor_count = 3;
        std::vector<std::unique_ptr<BufferSource>> srcs;
        for (int i = 0; i < 3; ++i)
            srcs.emplace_back(new BufferSource(static_cast<SourceId>(i), 256));
        Estimator est;
        REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& s : srcs) REQUIRE(est.add_source(s.get()) == Status::Ok);
        return est.deserialize(buf, len);
    };

    unsigned char buf[4096];
    // Written with the monitor OFF.
    const int len = serialize_under(false, 10.0, buf, sizeof(buf));
    CHECK(restores_under(false, 10.0, buf, len) == Status::Ok);          // same config: OK
    CHECK(restores_under(true,  10.0, buf, len) == Status::InvalidConfig); // flag flip: rejected
    CHECK(restores_under(false,  5.0, buf, len) == Status::InvalidConfig); // cap flip: rejected
}

// ===========================================================================
// Item 9 — the split-ON cov-cal guard + observability self-tests stay green is enforced BY THE
// GATE (the monitor defaults OFF; the spread->Q path is structurally untouched). We pin the
// structural invariant here directly: a monitor-ON run leaves the published per-channel
// reliabilities + the translation channel (the spread->Q inputs) byte-identical to monitor-OFF
// — the boost is a ROTATION-weight multiplier only, it never enters the covariance path.
// ===========================================================================
TEST_CASE("19c item 9: boost is rotation-weight-only — translation/cov path unchanged") {
    const EstRun on  = run_split_gps_rig(true,  400.0, 0.005);
    const EstRun off = run_split_gps_rig(false, 400.0, 0.005);
    // The translation-channel published weight + the drifter's per-channel reliabilities are
    // identical (the boost touches only weights_rot; the spread->Q sizing is structurally the
    // same). Note: the ROTATION consensus DOES move (that is the whole point), so we do NOT pin
    // r01 here — item 1 already checks it moved the right way.
    CHECK(on.w0  == off.w0);
    CHECK(on.w1  == off.w1);
    CHECK(on.w2  == off.w2);
    CHECK(on.rr2 == off.rr2);
    CHECK(on.rt2 == off.rt2);
    CHECK(on.tz  == off.tz);
}

// ===========================================================================
// REVIEW FIX (CRITICAL 1) -- BLOCK MEDIAN ROBUSTNESS. The 60-s block reduction is the MEDIAN of
// the block's (t, cum) samples (spec 1.3 / prototype), not the mean. A single bad pair inside an
// otherwise-clean block (a multipath spike / wheel-slip step that fits in 60 s) must NOT skew the
// block-reduced cumulative residual -- the median outvotes it; a mean would be pulled toward it.
// ===========================================================================
TEST_CASE("19c review: block reduction is MEDIAN -- an in-block outlier cluster does not skew it") {
    // Isolate the BLOCK reduction from the slope-pool median. Drive exactly enough fixes for THREE
    // 60-s blocks (b0,b1,b2): only the (b0,b2) pair clears the 120-s min baseline, so the slope
    // pool holds a SINGLE slope and the score equals |that one slope|. There is no pool of clean
    // slopes to outvote a corrupted one, so the ONLY robustness left is the per-block reduction.
    // We corrupt block b2 for a MINORITY (20 of ~60) of its samples: the block MEDIAN ignores the
    // cluster (score stays 0), while a block MEAN would shift b2's reduced cum by 4.5 deg*20/60 =
    // 1.5 deg and the lone slope would jump to ~2.2e-4 rad/s. This test PASSES under the median
    // and FAILS under a mean reduction (verified by temporarily swapping the reduction).
    const Scalar bump = Scalar(4.5) * kPi / Scalar(180.0);   // 4.5 deg, just inside the 5 deg clamp

    MonRig clean;
    clean.straight(185);                 // ~185 s -> exactly 3 blocks -> one (b0,b2) slope
    REQUIRE(clean.mon.scored(0));
    const Scalar s0_clean = clean.mon.score(0);   // a perfectly clean stream scores exactly 0

    MonRig spiked;
    spiked.straight(150);                // fill b0,b1 and most of b2
    spiked.yaw[0] += bump;               // source 0: hold a +4.5 deg offset over a 20-pair cluster
    spiked.straight(20);                 //   (a minority of block b2's ~60 samples)
    spiked.yaw[0] -= bump;               // remove the offset (cum returns to the clean staircase)
    spiked.straight(15);                 // close b2 at the same total length
    REQUIRE(spiked.mon.scored(0));
    const Scalar s0_spiked = spiked.mon.score(0);

    // MEDIAN robustness: the minority cluster is outvoted inside block b2, so the lone slope is
    // intact and source 0's score is unchanged from the clean run (both ~0 to estimator epsilon).
    CHECK(s0_spiked == doctest::Approx(s0_clean));

    // Discriminator vs the MEAN: a mean reduction would shift b2's reduced cum by ~1.5 deg, making
    // the lone (b0,b2) slope ~0.026 rad / 120 s ~ 2.2e-4 rad/s. Pin the result FAR below that -- a
    // mean reduction cannot pass this bound, a median reduction passes it trivially.
    CHECK(std::abs(s0_spiked - s0_clean) < Scalar(1e-5));   // mean would give ~2.2e-4 here
}

// ===========================================================================
// REVIEW FIX (CRITICAL 2) -- ZERO-DRIFT BOOST CONTINUITY. A source with EXACTLY zero (or sub-
// floor) drift must NOT jump to the full cap while a tiny-drift source gets less -- that breaks
// the ranking ratio. The score is floored at kScoreFloor (the GPS-course noise floor), so all
// below-resolvability sources rank equal and the boost is continuous.
// ===========================================================================
TEST_CASE("19c review: zero-drift boost is continuous (floored, not special-cased to the cap)") {
    using namespace heading_monitor_const;
    // All three sources perfectly clean (zero planted drift): every score is at/below the floor,
    // so every source floors to the SAME value -> equal ranking -> all boosted to the cap. The
    // old code returned `cap` for an exactly-zero score and a proportional value for a tiny one,
    // an unbounded discontinuity; the floor makes them coincide.
    MonRig r;
    r.straight(kScoreFixes);
    REQUIRE(r.mon.scored_count() == 3);
    // Sub-floor scores -> all three boosts equal (continuous, no zero-drift jump).
    const Scalar b0 = r.mon.boost(0, 10.0);
    const Scalar b1 = r.mon.boost(1, 10.0);
    const Scalar b2 = r.mon.boost(2, 10.0);
    CHECK(b0 == doctest::Approx(b1));
    CHECK(b1 == doctest::Approx(b2));
    CHECK(b0 == doctest::Approx(10.0));     // all at/below the floor -> all at the cap, none higher

    // Continuity check: two sub-floor scores never differ by more than the cap in their boost
    // (the floor clamps both numerator and denominator), so no source can be infinitely demoted
    // by another's near-zero score. The clean-clean ratio is exactly 1 (not cap/epsilon).
    CHECK(b0 / b1 == doctest::Approx(1.0));
}

// ===========================================================================
// REVIEW FIX (HIGH-RISK 3) + Slice 19d WARM-START -- ABSTAIN-ON BYTE-IDENTICAL FUSION. The OFF
// contract has TWO halves: the OFF flag is neutral (item 7), AND monitor-ON-but-ABSTAINING
// (< 120 s / < 2 scored) reproduces the monitor-OFF FUSION BYTE-IDENTICALLY. A short run (60 s)
// cannot form a slope baseline, so monitor-ON abstains; its frontier pose/cov/weights must
// bit-match the monitor-OFF run. Slice 19d changes ONLY the boost DIAGNOSTIC in this state: the
// abstain boost is now the clamped rot_weight_prior warm-start (not a flat 1.0), but because the
// estimator drops the separate *rot_weight_prior when the monitor is on, the rotation weight --
// and thus the fusion -- is unchanged (see the body).
// ===========================================================================
TEST_CASE("19d: abstain-ON fusion is byte-identical to OFF (warm-start ABSORBS the prior)") {
    // 60 s is below the 120 s slope baseline -> the ON monitor scores nothing -> it abstains.
    // Slice 19d WARM-START: an abstaining source returns its rot_weight_prior (clamped to
    // [1,cap]) as its boost, and the estimator DROPS the separate * rot_weight_prior when the
    // monitor is on. So the rotation weight is clamp(w_base*rel)*prior in BOTH runs:
    //   OFF:  w_rot = clamp(w_base*rel) * rot_weight_prior * 1.0(hm)
    //   ON:   w_rot = clamp(w_base*rel) * hm_boost(= warm-start = clamp(prior,[1,cap]))
    // -> the fusion is STILL bit-identical even with a NON-neutral prior (the rig sets
    // sensors[1].rot_weight_prior = 2.0). What CHANGES vs Slice 19c is ONLY the boost DIAGNOSTIC:
    // abstain no longer reports a flat 1.0, it reports the warm-start prior (this is the
    // intentional 19d behavior change -- the prior applies from t=0, no discovery latency).
    const EstRun on  = run_split_gps_rig(true,  60.0, 0.005);
    const EstRun off = run_split_gps_rig(false, 60.0, 0.005);

    // Abstain state: nothing scored. The boost now WARM-STARTS to the clamped prior:
    // sources 0/2 (prior 1.0) -> exactly 1.0; source 1 (prior 2.0) -> exactly 2.0.
    for (int i = 0; i < 3; ++i) CHECK_FALSE(on.scored[i]);
    CHECK(on.boost[0] == Scalar(1));      // prior 1.0 -> warm-start 1.0 (19c pin: prior=1 abstain)
    CHECK(on.boost[1] == Scalar(2));      // prior 2.0 -> warm-start 2.0 (the new behavior)
    CHECK(on.boost[2] == Scalar(1));      // prior 1.0 -> warm-start 1.0

    // Frontier pose: bit-identical to the OFF run (warm-start reproduces the OFF prior factor).
    CHECK(on.tx  == off.tx);
    CHECK(on.ty  == off.ty);
    CHECK(on.tz  == off.tz);
    CHECK(on.r01 == off.r01);
    CHECK(on.r12 == off.r12);
    // Frontier covariance: bit-identical.
    CHECK(on.c0  == off.c0);
    CHECK(on.c7  == off.c7);
    // Published (translation-channel) weights + per-channel reliabilities: bit-identical (the
    // warm-start touches only weights_rot, and there it reproduces the OFF prior exactly).
    CHECK(on.w0  == off.w0);
    CHECK(on.w1  == off.w1);
    CHECK(on.w2  == off.w2);
    CHECK(on.rr2 == off.rr2);
    CHECK(on.rt2 == off.rt2);
    // Correction bookkeeping: the same number of fixes applied (the monitor never gates fusion).
    CHECK(on.corr_applied_total == off.corr_applied_total);
}

// ===========================================================================
// Slice 19d ITEM 1 -- WARM-START (unit). While abstaining/unscored, boost() returns the source's
// configured rot_weight_prior CLAMPED to [1, boost_max], NOT a flat 1.0 -> the config prior is
// applied from t=0 with no discovery latency. Pinned directly on the HeadingMonitor.
// ===========================================================================
TEST_CASE("19d item 1: warm-start -- abstain boost is the clamped rot_weight_prior, not 1.0") {
    HeadingMonitor mon;
    mon.configure(3);
    mon.set_rot_weight_prior(0, Scalar(1));     // neutral
    mon.set_rot_weight_prior(1, Scalar(7));     // heading-grade prior (within the cap)
    mon.set_rot_weight_prior(2, Scalar(1));

    // No fixes submitted -> nothing scored -> ABSTAIN. The warm-start applies from the FIRST
    // step: the high-prior source already outweighs the neutral ones in the rotation channel.
    REQUIRE(mon.scored_count() == 0);
    CHECK(mon.boost(0, 10.0) == Scalar(1));     // neutral prior -> 1.0
    CHECK(mon.boost(1, 10.0) == Scalar(7));     // prior 7 -> warm-start 7 (NOT 1.0)
    CHECK(mon.boost(2, 10.0) == Scalar(1));
    // The rotation median is preferred toward source 1 from t=0: its boost > the others'.
    CHECK(mon.boost(1, 10.0) > mon.boost(0, 10.0));
    CHECK(mon.boost(1, 10.0) > mon.boost(2, 10.0));

    // A prior ABOVE the cap clamps to the cap (warm-start never exceeds boost_max).
    mon.set_rot_weight_prior(1, Scalar(25));
    CHECK(mon.boost(1, 10.0) == Scalar(10));    // clamp to boost_max
    // A prior BELOW 1 clamps UP to 1 (never demotes below neutral).
    mon.set_rot_weight_prior(1, Scalar(0.3));
    CHECK(mon.boost(1, 10.0) == Scalar(1));

    // A still-unscored source with two OTHERS scored also warm-starts to its prior (the
    // unscored-i path, not just the global abstain). configure clears the priors, so re-bind.
    HeadingMonitor mon2;
    mon2.configure(3);
    mon2.set_rot_weight_prior(0, Scalar(4));
    // (No fixes for source 0's track -> it stays unscored; with <2 scored we are in abstain
    // anyway, so this also exercises the warm-start return value being the prior, not 1.0.)
    CHECK(mon2.boost(0, 10.0) == Scalar(4));
}

// ===========================================================================
// Slice 19d ITEM 2 -- NO DOUBLE-COUNT. With the monitor ON the rotation weight is
// clamp(w_base*rel) * hm_boost (the SEPARATE *rot_weight_prior is DROPPED -- the monitor absorbs
// it). The published heading_boost IS that single rotation multiplier. The no-double-count
// invariants:
//   (a) once SCORED, the boost is the COMPUTED ranking value -- purely score-based, INDEPENDENT
//       of the configured prior. So sensors[1].rot_weight_prior = 2.0 vs 1.0 yields the SAME
//       scored boost (the old code's *prior*boost would have differed).
//   (b) the boost (== the rotation multiplier on wr) NEVER exceeds boost_max, even with a >1
//       prior. The old wiring's effective rotation factor was prior*boost = up to 2*10 = 20x;
//       absorbed, it is capped at boost_max (10x). Pinning heading_boost <= boost_max proves the
//       rotation channel is no longer double-scaled.
//   (c) the translation channel is prior-independent in BOTH runs (rot_weight_prior never
//       touches it) -- a sanity guard the change did not leak into translation.
// (Full trajectory equality is NOT pinned: the warm-start INTENTIONALLY differs in the abstain
// window -- prior=2 weights source 1's rotation 2x from t=0 while prior=1 weights it 1x -- so the
// trajectories legitimately diverge there. That early-window difference IS the erased latency;
// the no-double-count claim is specifically about the SCORED regime, items (a)/(b).)
// ===========================================================================
TEST_CASE("19d item 2: no double-count -- scored boost is prior-independent and capped") {
    // Default rig: sensors[1].rot_weight_prior = 2.0; same rig with the prior forced to 1.0.
    const EstRun prior2 = run_split_gps_rig(true, 400.0, 0.005, nullptr, Scalar(2.0));
    const EstRun prior1 = run_split_gps_rig(true, 400.0, 0.005, nullptr, Scalar(1.0));

    // Both fully scored by end-of-drive (the computed boost path is active, not the warm-start).
    for (int i = 0; i < 3; ++i) { CHECK(prior2.scored[i]); CHECK(prior1.scored[i]); }
    // (a) The scored boost is the score-based ranking value -> IDENTICAL regardless of the prior
    //     (the prior is absorbed; the old *prior*boost wiring would have differed by the 2x).
    CHECK(prior2.boost[0] == prior1.boost[0]);
    CHECK(prior2.boost[1] == prior1.boost[1]);
    CHECK(prior2.boost[2] == prior1.boost[2]);
    // (b) The boost (== the single rotation multiplier on wr) never exceeds boost_max even with
    //     the >1 prior -> no double-count (the old prior*boost would reach 2*10 = 20 here).
    for (int i = 0; i < 3; ++i) {
        CHECK(prior2.boost[i] <= Scalar(10.0));
        CHECK(prior2.boost[i] >= Scalar(1.0));
    }
    CHECK(prior2.boost[1] == doctest::Approx(10.0));   // clean source 1 -> the cap, NOT 2*cap
    // (c) The translation channel is untouched by rot_weight_prior in EITHER run.
    CHECK(prior2.w0 == prior1.w0);
    CHECK(prior2.w1 == prior1.w1);
    CHECK(prior2.w2 == prior1.w2);
}

// ===========================================================================
// Slice 19d ITEM 5 -- TRANSITION smoothness. The warm-start handoff (abstain=clamped prior ->
// scored=computed boost) must not chatter: the boost is bounded in [1, cap] across the whole
// transition, and the established slow-update / hysteresis machinery (item 6) still governs the
// scored regime. We drive a single rig from abstain through to scored and confirm the boost
// stays in-band the whole way (no spike, no inversion) and lands on the computed value.
// ===========================================================================
TEST_CASE("19d item 5: warm-start -> scored transition stays in-band (no chatter)") {
    // Source 1 carries a heading-grade prior of 6 (within the cap); sources 0/2 neutral. Source 2
    // is the planted drifter. The clean sources should end boosted toward the cap; the drifter
    // floored. Throughout, every boost stays within [1, cap].
    MonRig r(/*drift_idx=*/2, /*drift_wz=*/0.005);
    r.mon.set_rot_weight_prior(0, Scalar(1));
    r.mon.set_rot_weight_prior(1, Scalar(6));
    r.mon.set_rot_weight_prior(2, Scalar(1));

    // ABSTAIN: the warm-start gives source 1 its prior, the others 1.0.
    CHECK(r.mon.boost(0, 10.0) == Scalar(1));
    CHECK(r.mon.boost(1, 10.0) == Scalar(6));   // warm-start to the prior
    CHECK(r.mon.boost(2, 10.0) == Scalar(1));

    // Feed fixes in chunks; at every checkpoint the boosts stay bounded in [1, cap] (no chatter
    // spike on the warm-start -> scored handoff).
    for (int chunk = 0; chunk < 10; ++chunk) {
        r.straight(20);
        for (int i = 0; i < 3; ++i) {
            const Scalar b = r.mon.boost(i, 10.0);
            CHECK(b >= Scalar(1));
            CHECK(b <= Scalar(10));
        }
    }
    REQUIRE(r.mon.scored_count() == 3);
    // SCORED end-state: computed boost takes over (prior-independent). The clean sources are
    // boosted, the drifter floored -- the slow weighted-median ranking, NOT the prior.
    CHECK(r.mon.boost(2, 10.0) == doctest::Approx(1.0));   // drifter floored
    CHECK(r.mon.boost(0, 10.0) > Scalar(1.5));             // clean source boosted by SCORE
    CHECK(r.mon.boost(1, 10.0) > Scalar(1.5));             // (not by its prior 6 -> absorbed)
}
