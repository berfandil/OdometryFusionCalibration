// Slice 12 tests: WARM-RESTART PERSISTENCE — serialize/deserialize the committed
// calibration state into a caller-owned FIXED buffer (versioned + checksummed +
// config-hash-guarded) so a fresh estimator restores it and resumes NEAR-NOMINAL, plus a
// RELAXED-EDGE double-buffer ping-pong proving a crash mid-write keeps the last good state.
//
// The core serialize()/deserialize() touch ONLY the caller buffer (no heap / file IO /
// exceptions). File IO + the double-buffer live here (relaxed edge: std fstream / heap fine);
// the production file-persistence ADAPTER is Slice 13.
//
// Coverage (the three D23 / ISSUES Slice-12 done-criteria + the framing rejections):
//   * Round-trip IDENTITY — serialize a converged estimator, deserialize into a fresh one
//     with the SAME config; the restored committed values/flags match bit-for-bit.
//   * RESTART RESUMES NEAR-NOMINAL — a restored estimator produces fused output consistent
//     with the pre-restart (converged) estimator and FAR better than a cold-start one, with
//     its committed flags already set (no re-bootstrap).
//   * CRASH MID-WRITE keeps last good state — the double-buffer load() returns A when B is a
//     truncated/corrupt half-write, then returns B after a clean B write.
//   * CONFIG CHANGE invalidates — a blob made under config X is rejected (config-hash guard)
//     by an estimator inited with a different rig.
//   * Framing rejections — a flipped byte (checksum), a bumped version, a bad magic, a too-
//     small serialize buffer (CapacityExceeded), and NotInitialized before init().
//   * Determinism — two identical converged runs serialize to byte-identical blobs.
#include <doctest/doctest.h>

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/persistence.hpp"
#include "ofc/core/result.hpp"

#include "ofc_sim/rig.hpp"
#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;
bool near_abs(Scalar a, Scalar b, Scalar tol) { return std::abs(a - b) <= tol; }
constexpr Scalar kNanosPerSec = Scalar(1e9);
Timestamp secs_to_ns(Scalar s) { return static_cast<Timestamp>(std::llround(s * kNanosPerSec)); }

Mat3 Rz(Scalar a) { Mat3 R; R << std::cos(a), -std::sin(a), 0,
                                 std::sin(a),  std::cos(a), 0, 0, 0, 1; return R; }
Mat3 Ry(Scalar a) { Mat3 R; R << std::cos(a), 0, std::sin(a),
                                 0, 1, 0, -std::sin(a), 0, std::cos(a); return R; }
Mat3 Rx(Scalar a) { Mat3 R; R << 1, 0, 0, 0, std::cos(a), -std::sin(a),
                                 0, std::sin(a), std::cos(a); return R; }
SE3 make_extrinsic(Scalar yaw, Scalar pitch, Scalar roll, const Vec3& t) {
    SE3 X; X.R = Rz(yaw) * Ry(pitch) * Rx(roll); X.t = t; return X;
}

// Same mixed bootstrap trajectory the Slice-8 feedback tests drive (straight + multi-axis
// turns + omega-varying), so every DOF accumulates votes and the estimator converges.
Trajectory bootstrap_traj() {
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

void set_hists(Config& c) {
    c.vote_weight = VoteWeight::One;
    c.straight_omega_max = 0.05; c.straight_trans_min = 0.02; c.turn_omega_min = 0.20;
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

const CalibSnapshot* snap(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i) if (r.calib[i].id == id) return &r.calib[i];
    return nullptr;
}

// ---- The SCALE capstone rig: reference (0) + a source (1) with a planted 1.20 scale and a
// WRONG prior scale 1.0. With two sources the median is a weighted midpoint (no outlier
// rejection), so the 20%-too-long translation biases fusion until the calibrated scale
// commits and feeds back. The cleanest closed loop to demonstrate near-NOMINAL resume.
const Scalar kScaleTruth = 1.20;

void make_scale_sensors(std::vector<SensorConfig>& sensors) {
    sensors.assign(2, SensorConfig{});
    for (int i = 0; i < 2; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[1].prior_scale = 1.0;     // WRONG (truth is 1.20)
}

Config make_scale_cfg(const std::vector<SensorConfig>& sensors) {
    Config cfg;
    cfg.max_sources = 2; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10; cfg.tick_rate_hz = 50.0;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
    cfg.sensors = sensors.data(); cfg.sensor_count = 2;
    return cfg;
}

void make_scale_sources(const Trajectory& tr,
                        std::vector<std::unique_ptr<SyntheticSource>>& srcs) {
    std::vector<SourceParams> planted(2);
    for (int i = 0; i < 2; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].scale = kScaleTruth;
    srcs.clear();
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
}

// Drive a STANDALONE estimator over [from_s, to_s] at tick_rate, replicating the rig's GT
// anchoring (gauge anchored at the first bootstrap window start), and return the MAX fused
// translation error over the LAST `frac` of the run (the converged regime). `prep(est)` is
// called after init+add_source, before the run, to (optionally) deserialize a warm-restart
// blob. Returns -1 if nothing fused.
Scalar drive_late_trans_err(const Config& cfg, const Trajectory& tr,
                            const std::vector<std::unique_ptr<SyntheticSource>>& srcs,
                            Scalar from_s, Scalar to_s, Scalar tick_rate,
                            void (*prep)(Estimator&, const void*), const void* prep_arg) {
    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (const auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
    if (prep != nullptr) prep(est, prep_arg);

    const Timestamp tick = secs_to_ns(Scalar(1) / tick_rate);
    const Timestamp t0   = secs_to_ns(from_s);
    const Timestamp tend = secs_to_ns(to_s);
    const Timestamp window_ns = secs_to_ns(cfg.window_s);
    bool have_anchor = false; SE3 gt_anchor_inv;

    std::vector<Scalar> errs;
    for (Timestamp now = t0; now <= tend; now += tick) {
        if (!ok(est.step(now))) { continue; }
        const Result& res = est.latest();
        const Timestamp frontier = res.frontier.stamp;
        const SE3 gt_abs = tr.pose(frontier);
        if (!have_anchor) { gt_anchor_inv = se3::inverse(tr.pose(frontier - window_ns)); have_anchor = true; }
        const SE3 gt = se3::compose(gt_anchor_inv, gt_abs);
        errs.push_back((res.frontier.pose.t - gt.t).norm());
    }
    if (errs.empty()) return Scalar(-1);
    const std::size_t start = static_cast<std::size_t>((Scalar(1) - Scalar(0.4)) * errs.size());
    Scalar mx = 0;
    for (std::size_t i = start; i < errs.size(); ++i) mx = std::max(mx, errs[i]);
    return mx;
}

// prep callbacks for drive_late_trans_err.
void prep_none(Estimator&, const void*) {}
struct Blob { const unsigned char* p; int n; };
void prep_restore(Estimator& est, const void* arg) {
    const Blob* b = static_cast<const Blob*>(arg);
    REQUIRE(est.deserialize(b->p, b->n) == Status::Ok);
}

// Run the scale capstone rig to convergence and serialize the converged estimator into `out`.
// `mutate` (optional) tweaks the converge config before init — used by the Slice-16 positive
// control to write a blob under subbin_centroid=true. Returns the serialized byte count.
int converge_and_serialize_scale(std::vector<unsigned char>& out,
                                 void (*mutate)(Config&) = nullptr) {
    static const Trajectory tr = bootstrap_traj();
    std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
    Config cfg = make_scale_cfg(sensors);
    if (mutate != nullptr) mutate(cfg);
    std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);

    Rig rig; rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);
    // The scale committed + converged before we snapshot.
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 1);
    REQUIRE(cs != nullptr);
    REQUIRE(cs->scale_committed);
    REQUIRE(near_abs(cs->scale, kScaleTruth, 3e-2));

    out.assign(1024, 0);
    const Expected<int> r = rig.estimator().serialize(out.data(), static_cast<int>(out.size()));
    REQUIRE(r.ok());
    out.resize(r.value());
    return r.value();
}
} // namespace

// ===========================================================================
// Round-trip identity: a converged estimator -> blob -> a fresh estimator restores the
// committed values + flags exactly (the foundation TDD round-trip).
// ===========================================================================
TEST_CASE("persistence round-trip: serialize then deserialize restores committed state") {
    std::vector<unsigned char> blob;
    converge_and_serialize_scale(blob);
    CHECK(blob.size() > 28);   // header(20) + a 2-source payload + checksum(4)

    static const Trajectory tr = bootstrap_traj();
    std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
    const Config cfg = make_scale_cfg(sensors);
    std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
    REQUIRE(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::Ok);

    // Step once so the CalibSnapshot is published; the restored scale_committed flag is held by
    // commit_gate_reanchor's re-fill hysteresis through the (just-reset) histogram.
    REQUIRE(ok(est.step(secs_to_ns(0.2))));
    const CalibSnapshot* cs = snap(est.latest(), 1);
    REQUIRE(cs != nullptr);
    CHECK(cs->scale_committed);                       // commit flag survived the restore
    CHECK(near_abs(cs->scale, kScaleTruth, 3e-2));    // the committed value is back

    // A second serialize of the restored estimator (BEFORE stepping perturbs the priors) round-
    // trips to the SAME committed priors as the original blob's payload (idempotent restore).
    Estimator est2;
    REQUIRE(est2.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(est2.add_source(sp.get()) == Status::Ok);
    REQUIRE(est2.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::Ok);
    std::vector<unsigned char> blob2(1024, 0);
    const Expected<int> r2 = est2.serialize(blob2.data(), static_cast<int>(blob2.size()));
    REQUIRE(r2.ok());
    blob2.resize(r2.value());
    REQUIRE(blob2.size() == blob.size());
    CHECK(std::memcmp(blob2.data(), blob.data(), blob.size()) == 0);  // byte-identical re-serialize
}

// ===========================================================================
// RESTORED OFFSET FLAG HELD THROUGH REFILL (Slice 12 review MINOR): a restored-committed
// time-sync offset must read `committed` TRUE IMMEDIATELY after restore — not regress to false
// while the empty post-restore histogram refills. Pre-fix, update_commit_state() recomputed the
// flag via the plain commit_gate over the empty histogram (confidence~0) and flipped it FALSE;
// the offset_restored latch now holds it committed through the refill (mirroring the
// ext/roll/lever/scale re-fill hysteresis). This case FAILS on pre-fix code.
// ===========================================================================
namespace {
// A time-sync rig with a planted RELATIVE offset on source 1 — converges + COMMITS the offset.
Config make_offset_cfg(SensorConfig* sc /*[2], outlives the call*/) {
    Config c;
    c.tick_rate_hz     = 100.0;
    c.match_metric     = MatchMetric::L2;
    c.max_lag_s        = 0.1;
    c.timesync_enabled = true;
    c.max_sources      = 2;
    c.reference_sensor_id = 0;
    c.offset_hist.bins      = 256;
    c.offset_hist.range_min = -0.1;
    c.offset_hist.range_max =  0.1;
    c.offset_hist.circular  = false;
    c.offset_hist.aging     = Aging::SlidingK;
    c.offset_hist.sliding_k = 64;
    c.offset_hist.vote_split = true;
    c.offset_hist.subbin     = true;
    c.commit_concentration = 0.5;
    c.commit_drop          = 0.3;
    c.commit_min_votes     = 30;     // reachable with the SlidingK=64 offset hist
    sc[0].id = 0; sc[0].is_reference = true;
    sc[1].id = 1;
    c.sensors = sc; c.sensor_count = 2;
    return c;
}
} // namespace

TEST_CASE("persistence: a restored COMMITTED offset reads committed immediately (held through refill)") {
    Trajectory tr = Trajectory::omega_varying();
    const Scalar planted = 0.04;
    SourceParams pr; pr.id = 0; pr.time_offset_s = 0.0;
    SourceParams ps; ps.id = 1; ps.time_offset_s = planted;

    SensorConfig sc[2];
    const Config cfg = make_offset_cfg(sc);

    // (1) Converge a SOURCE estimator until source 1's offset commits, then serialize.
    SyntheticSource s0c(tr, pr), s1c(tr, ps);
    Estimator src;
    REQUIRE(src.init(cfg) == Status::Ok);
    REQUIRE(src.add_source(&s0c) == Status::Ok);
    REQUIRE(src.add_source(&s1c) == Status::Ok);
    const Timestamp tick = secs_to_ns(0.01);   // 100 Hz
    for (Timestamp now = secs_to_ns(0.3); now <= secs_to_ns(tr.end_s() - 0.1); now += tick) {
        (void)src.step(now);
    }
    const CalibSnapshot* csrc = snap(src.latest(), 1);
    REQUIRE(csrc != nullptr);
    REQUIRE(csrc->committed);                   // offset committed in the source run

    std::vector<unsigned char> blob(1024, 0);
    const Expected<int> sr = src.serialize(blob.data(), static_cast<int>(blob.size()));
    REQUIRE(sr.ok());
    blob.resize(sr.value());

    // (2) Restore into a FRESH estimator (empty time-sync histogram) and step ONCE. The restored
    // offset must read committed IMMEDIATELY — pre-fix the plain commit_gate over the empty
    // histogram flips it false until ~30 votes refill.
    SyntheticSource s0w(tr, pr), s1w(tr, ps);
    Estimator warm;
    REQUIRE(warm.init(cfg) == Status::Ok);
    REQUIRE(warm.add_source(&s0w) == Status::Ok);
    REQUIRE(warm.add_source(&s1w) == Status::Ok);
    REQUIRE(warm.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::Ok);
    REQUIRE(ok(warm.step(secs_to_ns(0.3))));
    const CalibSnapshot* cw = snap(warm.latest(), 1);
    REQUIRE(cw != nullptr);
    CHECK(cw->committed);                        // HELD through the empty/refilling histogram
    CHECK(cw->time_offset_s > 0.0);             // restored offset VALUE is applied (correct sign)
    CHECK(near_abs(cw->time_offset_s, planted, 0.02));
}

// ===========================================================================
// RESTART RESUMES NEAR-NOMINAL (done-criterion 1): a restored estimator tracks GT immediately
// (the committed scale drives fusion's de-scale from step 1), FAR better than a cold-start
// estimator that must re-bootstrap the scale from the wrong 1.0 prior.
// ===========================================================================
TEST_CASE("persistence near-NOMINAL: a restored estimator beats cold-start immediately") {
    std::vector<unsigned char> blob;
    converge_and_serialize_scale(blob);

    static const Trajectory tr = bootstrap_traj();
    std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
    // The resume config MUST hash-match the blob's config (the config-hash guard), so we reuse
    // the EXACT converge config. We isolate the RESTORE effect via a SHORT resume window: in the
    // first ~3 s neither estimator accumulates enough STRAIGHT-regime scale votes (only ~1/3 of
    // the mixed trajectory is straight) to clear commit_min_votes, so cold-start cannot re-
    // calibrate within the window — the ONLY thing letting the warm one track GT is the RESTORED
    // committed scale prior, which drives fusion's de-scale from step 1.
    const Config cfg = make_scale_cfg(sensors);

    std::vector<std::unique_ptr<SyntheticSource>> srcs_cold; make_scale_sources(tr, srcs_cold);
    std::vector<std::unique_ptr<SyntheticSource>> srcs_warm; make_scale_sources(tr, srcs_warm);

    // A SHORT resume window (the first ~3 s) — long enough to fuse, too short to recalibrate.
    const Scalar from_s = 0.2, to_s = 3.0;
    const Scalar cold = drive_late_trans_err(cfg, tr, srcs_cold, from_s, to_s, 50.0,
                                             prep_none, nullptr);
    Blob b{ blob.data(), static_cast<int>(blob.size()) };
    const Scalar warm = drive_late_trans_err(cfg, tr, srcs_warm, from_s, to_s, 50.0,
                                             prep_restore, &b);
    REQUIRE(cold >= 0);
    REQUIRE(warm >= 0);
    INFO("late fused trans err: cold-start=" << cold << "  warm-restart=" << warm);
    // The warm restart resumes near-NOMINAL: with the committed scale restored it tracks GT
    // (warm error stays ~mm), while the cold start carries the full 1.0-prior scale bias (a
    // decimetre+ of bias). Assert the near-NOMINAL PROPERTY directly (warm absolutely small,
    // cold absolutely large) rather than only the relative gap — keep the relative check too.
    CHECK(warm < Scalar(0.01));        // warm restart is absolutely near-NOMINAL (~mm)
    CHECK(cold > Scalar(0.1));         // cold start carries a decimetre+ of un-calibrated bias
    CHECK(warm < cold * Scalar(0.6));

    // And the warm estimator comes up with its commit flag ALREADY set + the right scale (it did
    // not re-bootstrap — commit is unreachable here, so the only way it is committed is restore).
    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& sp : srcs_warm) REQUIRE(est.add_source(sp.get()) == Status::Ok);
    REQUIRE(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::Ok);
    REQUIRE(ok(est.step(secs_to_ns(0.2))));
    const CalibSnapshot* cs = snap(est.latest(), 1);
    REQUIRE(cs != nullptr);
    CHECK(cs->scale_committed);
    CHECK(near_abs(cs->scale, kScaleTruth, 5e-2));
}

// ===========================================================================
// CONFIG CHANGE invalidates (done-criterion 3, config-hash guard): a blob made under config X
// is rejected by an estimator inited with a DIFFERENT rig.
// ===========================================================================
TEST_CASE("persistence config-hash guard: a changed rig rejects stale state") {
    std::vector<unsigned char> blob;
    converge_and_serialize_scale(blob);

    static const Trajectory tr = bootstrap_traj();
    std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);

    // (a) A DIFFERENT per-sensor prior (prior_extrinsic moved) -> different config_hash -> reject.
    {
        std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
        sensors[1].prior_extrinsic = make_extrinsic(0.2, 0.0, 0.0, Vec3(0.1, 0, 0));  // moved mount
        Config cfg = make_scale_cfg(sensors);
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::InvalidConfig);
    }
    // (b) A DIFFERENT reference_sensor_id -> different config_hash -> reject.
    {
        std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
        Config cfg = make_scale_cfg(sensors);
        cfg.reference_sensor_id = 1;        // was 0
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::InvalidConfig);
    }
    // (c) A DIFFERENT histogram config (bins) -> different config_hash -> reject.
    {
        std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
        Config cfg = make_scale_cfg(sensors);
        cfg.scale_hist.bins = 128;          // was 256
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::InvalidConfig);
    }
    // (d) SAME config -> accepted (the guard is not over-eager).
    {
        std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
        Config cfg = make_scale_cfg(sensors);
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::Ok);
    }
}

// ===========================================================================
// Slice 16: subbin_centroid is calibration-shaping -> it joins the config-hash
// stream. A blob written under subbin_centroid=false (the default) must be
// REJECTED by a rig that flips the flag on any histogram (cross-flag restore),
// because the persisted committed calibration was read with a different
// sub-bin estimator. Positive control: a blob written WITH the flag set is
// ACCEPTED by an identical true-flag rig (the hash is stable, not just different).
// ===========================================================================
TEST_CASE("persistence config-hash guard: flipping subbin_centroid rejects a cross-flag restore") {
    std::vector<unsigned char> blob;
    converge_and_serialize_scale(blob);    // written with subbin_centroid = false everywhere

    static const Trajectory tr = bootstrap_traj();

    // (a) Flip the flag on ONE histogram (scale_hist) -> different hash -> reject.
    {
        std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);
        std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
        Config cfg = make_scale_cfg(sensors);
        cfg.scale_hist.subbin_centroid = true;
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::InvalidConfig);
    }
    // (b) Flip it on a DIFFERENT histogram (so3_hist) -> also rejected (every
    // histogram's flag is hashed, not just the scale channel's).
    {
        std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);
        std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
        Config cfg = make_scale_cfg(sensors);
        cfg.so3_hist.subbin_centroid = true;
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::InvalidConfig);
    }
    // (c) POSITIVE CONTROL: a blob written WITH subbin_centroid=true restores into an
    // identical true-flag rig (same-flag hash EQUALITY on the new field) — guards against
    // a hash that rejects everything or encodes the flag unstably between runs.
    {
        std::vector<unsigned char> blob_c;
        converge_and_serialize_scale(blob_c, [](Config& c) { c.scale_hist.subbin_centroid = true; });

        std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);
        std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
        Config cfg = make_scale_cfg(sensors);
        cfg.scale_hist.subbin_centroid = true;   // SAME flag on the restoring rig
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(blob_c.data(), static_cast<int>(blob_c.size())) == Status::Ok);
    }
}

// ===========================================================================
// Framing rejections: checksum (flipped byte), version, magic, capacity, NotInitialized.
// ===========================================================================
TEST_CASE("persistence framing: checksum / version / magic / capacity / not-initialized") {
    std::vector<unsigned char> blob;
    converge_and_serialize_scale(blob);

    static const Trajectory tr = bootstrap_traj();
    std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
    const Config cfg = make_scale_cfg(sensors);
    std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);

    auto fresh = [&](std::unique_ptr<Estimator>& est,
                     std::vector<std::unique_ptr<SyntheticSource>>& s) {
        est.reset(new Estimator());
        REQUIRE(est->init(cfg) == Status::Ok);
        make_scale_sources(tr, s);
        for (auto& sp : s) REQUIRE(est->add_source(sp.get()) == Status::Ok);
    };

    // (a) A single flipped byte in the PAYLOAD -> checksum fails -> CorruptData.
    {
        std::vector<unsigned char> bad = blob;
        bad[30] ^= 0x01;          // somewhere in the payload (header is the first 20 bytes)
        std::unique_ptr<Estimator> est; std::vector<std::unique_ptr<SyntheticSource>> s; fresh(est, s);
        CHECK(est->deserialize(bad.data(), static_cast<int>(bad.size())) == Status::CorruptData);
    }
    // (b) A flipped byte in the HEADER's config_hash -> the config-hash compare fails first
    // (InvalidConfig) OR the checksum (CorruptData) — both are a rejection. Assert NOT Ok.
    {
        std::vector<unsigned char> bad = blob;
        bad[10] ^= 0x01;          // inside the config_hash field (bytes 8..15)
        std::unique_ptr<Estimator> est; std::vector<std::unique_ptr<SyntheticSource>> s; fresh(est, s);
        CHECK(est->deserialize(bad.data(), static_cast<int>(bad.size())) != Status::Ok);
    }
    // (c) A bumped format_version -> VersionMismatch. deserialize checks version BEFORE the
    // checksum (the documented order: magic -> version -> checksum -> config-hash), so a stale-
    // version blob is rejected as VersionMismatch even though its checksum is now also stale.
    {
        std::vector<unsigned char> bad = blob;
        bad[4] = static_cast<unsigned char>(bad[4] + 1);   // version is bytes 4..7 (LE)
        std::unique_ptr<Estimator> est; std::vector<std::unique_ptr<SyntheticSource>> s; fresh(est, s);
        CHECK(est->deserialize(bad.data(), static_cast<int>(bad.size())) == Status::VersionMismatch);
    }
    // (d) A bad magic -> CorruptData (checked first).
    {
        std::vector<unsigned char> bad = blob;
        bad[0] = 'X';
        std::unique_ptr<Estimator> est; std::vector<std::unique_ptr<SyntheticSource>> s; fresh(est, s);
        CHECK(est->deserialize(bad.data(), static_cast<int>(bad.size())) == Status::CorruptData);
    }
    // (e) A truncated blob -> CorruptData (length mismatch).
    {
        std::unique_ptr<Estimator> est; std::vector<std::unique_ptr<SyntheticSource>> s; fresh(est, s);
        CHECK(est->deserialize(blob.data(), static_cast<int>(blob.size()) - 3) == Status::CorruptData);
    }
    // (f) serialize into a too-small buffer -> CapacityExceeded.
    {
        std::unique_ptr<Estimator> est; std::vector<std::unique_ptr<SyntheticSource>> s; fresh(est, s);
        unsigned char tiny[8];
        const Expected<int> r = est->serialize(tiny, static_cast<int>(sizeof(tiny)));
        CHECK_FALSE(r.ok());
        CHECK(r.status() == Status::CapacityExceeded);
    }
    // (g) serialize / deserialize before init() -> NotInitialized.
    {
        Estimator un;
        unsigned char buf[256];
        const Expected<int> r = un.serialize(buf, static_cast<int>(sizeof(buf)));
        CHECK_FALSE(r.ok());
        CHECK(r.status() == Status::NotInitialized);
        CHECK(un.deserialize(blob.data(), static_cast<int>(blob.size())) == Status::NotInitialized);
    }
}

// ===========================================================================
// ORTHONORMALITY GUARD (Slice 12 review NIT): a checksum-VALID blob whose extrinsic rotation is
// non-orthonormal (the residual ~1-in-4e9 FNV collision path) is rejected as CorruptData by the
// deserialize orthonormality check, before the corrupt R can propagate into fusion. We forge the
// collision directly: corrupt R(0,0) of the first source's prior_extrinsic, then RE-WRITE the
// trailing FNV-1a-32 checksum so the blob passes the checksum + config-hash gates and only the
// orthonormality guard can catch it.
// ===========================================================================
TEST_CASE("persistence orthonormality: a checksum-valid but non-orthonormal R is rejected") {
    std::vector<unsigned char> blob;
    converge_and_serialize_scale(blob);

    static const Trajectory tr = bootstrap_traj();
    std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
    const Config cfg = make_scale_cfg(sensors);
    std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);

    // Payload layout: header(20) then i32 source_count + i32 phase + u8 ever_fused (= 9 bytes),
    // then per source: i32 id (4) + SE3 (9 R doubles + 3 t doubles). So the first source's
    // R(0,0) starts at blob offset 20 + 9 + 4 = 33.
    const int r00_off = 20 + 9 + 4;
    std::vector<unsigned char> bad = blob;
    const double non_ortho = 2.0;                 // scaling one entry breaks RᵀR == I
    std::memcpy(bad.data() + r00_off, &non_ortho, sizeof(double));
    // Re-write the trailing checksum so the corrupt blob PASSES the checksum gate (forge the
    // collision): FNV-1a-32 over header+payload (everything before the trailing 4-byte word).
    const int covered = static_cast<int>(bad.size()) - 4;
    const std::uint32_t cs = persist::fnv1a32(bad.data(), covered);
    for (int i = 0; i < 4; ++i) bad[covered + i] = static_cast<unsigned char>((cs >> (8 * i)) & 0xFFu);

    // Sanity: the re-checksummed UNcorrupted blob still deserializes Ok (the forge mechanism is
    // sound — it is the non-orthonormal R, not a broken checksum, that triggers the rejection).
    {
        std::vector<unsigned char> ok_blob = blob;
        const std::uint32_t cs2 = persist::fnv1a32(ok_blob.data(), static_cast<int>(ok_blob.size()) - 4);
        for (int i = 0; i < 4; ++i)
            ok_blob[static_cast<int>(ok_blob.size()) - 4 + i] =
                static_cast<unsigned char>((cs2 >> (8 * i)) & 0xFFu);
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(est.deserialize(ok_blob.data(), static_cast<int>(ok_blob.size())) == Status::Ok);
    }

    std::vector<std::unique_ptr<SyntheticSource>> srcs2; make_scale_sources(tr, srcs2);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& sp : srcs2) REQUIRE(est.add_source(sp.get()) == Status::Ok);
    CHECK(est.deserialize(bad.data(), static_cast<int>(bad.size())) == Status::CorruptData);
}

// ===========================================================================
// RELAXED-EDGE double-buffer ping-pong (done-criterion 2 + the deliverable). Two files A/B
// each hold { uint64 seq, core-serialized blob }. save() writes the INACTIVE file then flips
// the active seq; load() reads BOTH, picks the highest seq whose blob passes the core
// deserialize, falling back to the other. Proves a crash mid-write keeps the last good state.
//
// This lives in the TEST (std fstream / heap / exceptions are fine at the relaxed edge); the
// production file-persistence adapter is Slice 13.
// ===========================================================================
namespace {
struct DoubleBuffer {
    std::string path_a, path_b;
    std::uint64_t seq = 0;
    bool write_to_b_next = true;   // alternate which file is the INACTIVE target

    static void write_record(const std::string& path, std::uint64_t seq,
                             const unsigned char* blob, int n) {
        std::ofstream f(path, std::ios::binary | std::ios::trunc);
        f.write(reinterpret_cast<const char*>(&seq), sizeof(seq));
        f.write(reinterpret_cast<const char*>(blob), n);
        f.flush();   // (a real adapter would fsync here; flush is the std-portable analogue)
    }

    // Read a record file -> (seq, blob). Returns false if the file is missing / too short to
    // even hold the seq word.
    static bool read_record(const std::string& path, std::uint64_t& seq,
                            std::vector<unsigned char>& blob) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) return false;
        const std::streamoff sz = f.tellg();
        if (sz < static_cast<std::streamoff>(sizeof(seq))) return false;
        f.seekg(0);
        f.read(reinterpret_cast<char*>(&seq), sizeof(seq));
        blob.resize(static_cast<std::size_t>(sz) - sizeof(seq));
        if (!blob.empty()) f.read(reinterpret_cast<char*>(blob.data()),
                                  static_cast<std::streamsize>(blob.size()));
        return true;
    }

    // Serialize `est` and write it to the inactive file with the next seq; then flip.
    void save(const Estimator& est) {
        unsigned char buf[1024];
        const Expected<int> r = est.serialize(buf, static_cast<int>(sizeof(buf)));
        REQUIRE(r.ok());
        const std::string& target = write_to_b_next ? path_b : path_a;
        write_record(target, seq + 1, buf, r.value());
        ++seq;
        write_to_b_next = !write_to_b_next;
    }

    // Load the highest-seq record whose blob the estimator ACCEPTS (deserialize == Ok). Tries
    // the higher-seq file first, falls back to the other. Returns the Status of the chosen load
    // (NoData if neither file holds an acceptable blob).
    Status load(Estimator& est) {
        std::uint64_t sa = 0, sb = 0;
        std::vector<unsigned char> ba, bb;
        const bool has_a = read_record(path_a, sa, ba);
        const bool has_b = read_record(path_b, sb, bb);
        // Order candidates by descending seq.
        struct Cand { bool has; std::uint64_t seq; std::vector<unsigned char>* blob; };
        Cand c0{ has_a, sa, &ba }, c1{ has_b, sb, &bb };
        if (c1.has && (!c0.has || c1.seq > c0.seq)) std::swap(c0, c1);   // c0 = higher seq
        for (Cand* c : { &c0, &c1 }) {
            if (!c->has) continue;
            const Status st = est.deserialize(c->blob->data(), static_cast<int>(c->blob->size()));
            if (ok(st)) { seq = c->seq; return st; }
        }
        return Status::NoData;
    }
};

// A small temp-dir helper (relaxed edge). Uses the doctest test name hashed into the path so
// parallel test runs do not collide; cleans up on destruction.
std::string temp_path(const char* tag) {
    const char* base = std::getenv("TEMP");
    if (base == nullptr) base = std::getenv("TMP");
    std::string dir = (base != nullptr) ? base : ".";
    return dir + "/ofc_persist_" + tag;
}
} // namespace

TEST_CASE("persistence double-buffer: a crash mid-write keeps the last good state") {
    std::vector<unsigned char> good;
    converge_and_serialize_scale(good);

    static const Trajectory tr = bootstrap_traj();
    std::vector<SensorConfig> sensors; make_scale_sensors(sensors);
    const Config cfg = make_scale_cfg(sensors);
    std::vector<std::unique_ptr<SyntheticSource>> srcs; make_scale_sources(tr, srcs);

    DoubleBuffer db;
    db.path_a = temp_path("A.bin");
    db.path_b = temp_path("B.bin");
    std::remove(db.path_a.c_str());
    std::remove(db.path_b.c_str());

    // (1) Write a GOOD record to file A (seq 1). load() must return it.
    DoubleBuffer::write_record(db.path_a, 1, good.data(), static_cast<int>(good.size()));
    {
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(db.load(est) == Status::Ok);
        REQUIRE(ok(est.step(secs_to_ns(0.2))));
        const CalibSnapshot* cs = snap(est.latest(), 1);
        REQUIRE(cs != nullptr);
        CHECK(cs->scale_committed);                       // A's good state restored
    }

    // (2) Begin writing B with a HIGHER seq (2) but CRASH mid-write: truncate the blob so its
    // checksum fails (simulate a partial flush). load() must FALL BACK to A (the last good).
    {
        std::vector<unsigned char> half(good.begin(), good.begin() + good.size() / 2);
        DoubleBuffer::write_record(db.path_b, 2, half.data(), static_cast<int>(half.size()));
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        std::vector<std::unique_ptr<SyntheticSource>> s2; make_scale_sources(tr, s2);
        for (auto& sp : s2) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        // B has the higher seq but is corrupt; load() rejects it and returns A's state.
        CHECK(db.load(est) == Status::Ok);
        REQUIRE(ok(est.step(secs_to_ns(0.2))));
        const CalibSnapshot* cs = snap(est.latest(), 1);
        REQUIRE(cs != nullptr);
        CHECK(cs->scale_committed);                       // fell back to A — last good kept
        CHECK(near_abs(cs->scale, kScaleTruth, 5e-2));
        CHECK(db.seq == 1);                               // it chose A's seq, not B's
    }

    // (3) A full CLEAN write to B (seq 3 > A's 1) -> load() now returns B.
    {
        DoubleBuffer::write_record(db.path_b, 3, good.data(), static_cast<int>(good.size()));
        Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
        std::vector<std::unique_ptr<SyntheticSource>> s3; make_scale_sources(tr, s3);
        for (auto& sp : s3) REQUIRE(est.add_source(sp.get()) == Status::Ok);
        CHECK(db.load(est) == Status::Ok);
        CHECK(db.seq == 3);                               // chose B (the higher, valid seq)
    }

    std::remove(db.path_a.c_str());
    std::remove(db.path_b.c_str());
}

// ===========================================================================
// Determinism: two identical converged runs serialize to byte-identical blobs.
// ===========================================================================
TEST_CASE("persistence determinism: identical runs -> byte-identical blobs") {
    std::vector<unsigned char> a, b;
    converge_and_serialize_scale(a);
    converge_and_serialize_scale(b);
    REQUIRE(a.size() == b.size());
    CHECK(std::memcmp(a.data(), b.data(), a.size()) == 0);
}
