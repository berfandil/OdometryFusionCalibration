// Slice 20 tests: TRANSLATION-ONLY source lever (velocity/Doppler sensors).
//
// A translation-only source (Doppler radar, optical-flow) measures its own VELOCITY but
// not its rotation (R_B ~ I, heading-blind). Today such a source gets lever conf 0 (never
// commits) AND commits a confidently-WRONG rotation extrinsic (the nuScenes front radar
// committed rx=-0.984 rad at conf 0.95). The translation_only flag PINS the rotation
// extrinsic to the prior and runs the existing hand-eye lever LS with that trusted R_X,
// recovering the observable lever axes.
//
// THE MATH (the prototype tools/proto_radar_lever.py pins it). The hand-eye translation
// identity (R_A - I) t_X = R_X t_B - t_A never uses the source's own rotation. For a
// heading-blind sensor reporting R_B = I and the sensor-frame displacement t_B that
// physically carries the lever (v_sensor = v_ego + omega x lever), the IDEAL t_B is
//   t_B = R_X^T [ t_A + (R_A - I) t_X ]
// (solve the identity for t_B). Accumulate over turn windows -> 3x3 normal eqs -> t_X.
// On planar (yaw-only) driving (R_A - I) has the z-axis in its null space -> z UNOBSERVABLE.
//
// Coverage (SLICE20 acceptance section 4, items 1-8):
//   1. HEADLINE: flag OFF -> garbage rotation / lever conf 0; flag ON -> observable lever
//      axes recover + rotation extrinsic stays = prior (not garbage). Both pinned.
//   2. Footgun: flag ON -> the rotation extrinsic never commits off the prior.
//   3. Along-track / noisy axis: a source clean on one axis, noisy on another -> commits
//      only the clean axis (the existing concentration gate withholds the noisy one).
//   4. z / planar null: yaw-only turns -> z withheld, x,y observable.
//   5. Default-off byte-identical (the untouched existing suite + an exact-equality pin).
//   6. Scale preserved for a translation_only source.
//   7. Loader key + config-hash flip; validate() rejects nonsense. (Loader key lives in
//      adapters/tests/test_config_loader.cpp; here: validate() + the hash flip.)
//   8. Observability self-test: a translation_only case converges in its regime.
#include <doctest/doctest.h>

#include "ofc/core/calibration.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/persistence.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/source.hpp"

#include "ofc_sim/rig.hpp"
#include "ofc_sim/synthetic_source.hpp"
#include "ofc_sim/trajectory.hpp"

#include <cmath>
#include <memory>
#include <random>
#include <vector>

using namespace ofc;
using namespace ofc::sim;

namespace {
constexpr Scalar kPi = 3.14159265358979323846;

Mat3 Rz(Scalar a) { Mat3 R; R << std::cos(a), -std::sin(a), 0,
                                 std::sin(a),  std::cos(a), 0,
                                 0, 0, 1; return R; }
Mat3 Ry(Scalar a) { Mat3 R; R << std::cos(a), 0, std::sin(a),
                                 0, 1, 0,
                                -std::sin(a), 0, std::cos(a); return R; }
Mat3 Rx(Scalar a) { Mat3 R; R << 1, 0, 0,
                                 0, std::cos(a), -std::sin(a),
                                 0, std::sin(a),  std::cos(a); return R; }

// Sensor->base extrinsic = Rz(yaw) Ry(pitch) Rx(roll) (config_loader convention).
SE3 make_extrinsic(Scalar yaw, Scalar pitch, Scalar roll, const Vec3& t) {
    SE3 X; X.R = Rz(yaw) * Ry(pitch) * Rx(roll); X.t = t;
    return X;
}

std::unique_ptr<Phase2Calibrator> make_p2() {
    return std::unique_ptr<Phase2Calibrator>(new Phase2Calibrator());
}
std::unique_ptr<Phase1Calibrator> make_p1() {
    return std::unique_ptr<Phase1Calibrator>(new Phase1Calibrator());
}

// Phase-2 config with the histogram shapes the other calib tests use. vote_weight One so
// vote_count is a literal count.
Config make_p2_cfg(bool rot3d = false) {
    Config c;
    c.tick_rate_hz   = 50.0;
    c.reference_sensor_id = 0;
    c.turn_omega_min = 0.20;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.phase2_strategy = Phase2Strat::VsFusedBase;
    c.vote_weight     = VoteWeight::One;
    c.rot3d_enabled   = rot3d;

    c.so3_hist.bins = 512; c.so3_hist.range_min = -0.8; c.so3_hist.range_max = 0.8;
    c.so3_hist.aging = Aging::SlidingK; c.so3_hist.sliding_k = 256; c.so3_hist.vote_split = true;
    c.so3_hist.subbin = true;
    c.scale_hist.bins = 512; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 256;
    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 256;
    c.xyz_hist.bins = 1024; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 256; c.xyz_hist.subbin = true;
    return c;
}

// One IDEAL translation-only hand-eye window. The base motion A turns about `axis` by
// `theta` with translation `t_A`; the heading-blind sensor reports R_B = I and the
// lever-bearing displacement t_B = R_X^T [ t_A + (R_A - I) t_X ]. Optional per-axis
// SENSOR noise on t_B in the SENSOR frame (the radar Doppler drift) -> sigma per axis.
SE3 trans_only_window(const SE3& X, const Vec3& t_X, const Vec3& axis, Scalar theta,
                      const Vec3& t_A, const Vec3& sensor_noise = Vec3::Zero(),
                      std::mt19937_64* rng = nullptr) {
    SE3 A;
    A.R = so3::exp(axis.normalized() * theta);
    A.t = t_A;
    Vec3 t_B = X.R.transpose() * (t_A + (A.R - Mat3::Identity()) * t_X);
    if (rng != nullptr && sensor_noise.squaredNorm() > Scalar(0)) {
        std::normal_distribution<Scalar> nx(0.0, sensor_noise.x());
        std::normal_distribution<Scalar> ny(0.0, sensor_noise.y());
        std::normal_distribution<Scalar> nz(0.0, sensor_noise.z());
        const Scalar e0 = nx(*rng), e1 = ny(*rng), e2 = nz(*rng);
        t_B += Vec3(e0, e1, e2);
    }
    SE3 B; B.R = Mat3::Identity(); B.t = t_B;
    return B;
}

// Feed `count` ideal translation-only windows of base motion A (yaw-only when !multiaxis,
// else cycling 3 distinct axes for a non-planar z-observable regime), returning A per
// window so the estimator-path mirror can reuse it. `sensor_noise` is the per-axis t_B
// sensor noise.
void feed_trans_only(Phase2Calibrator& cal, SourceId id, const SE3& X, const Vec3& t_X,
                     int count, bool multiaxis, const Vec3& sensor_noise = Vec3::Zero(),
                     unsigned seed = 7) {
    std::mt19937_64 rng(seed);
    const Vec3 axes[3] = { Vec3(0, 0, 1), Vec3(0, 1, 0.2), Vec3(0.3, -0.4, 0.85) };
    const Scalar dt = Scalar(0.5);
    for (int k = 0; k < count; ++k) {
        const Vec3 axis = multiaxis ? axes[k % 3] : Vec3(0, 0, 1);
        const Scalar theta = Scalar(0.30) + Scalar(0.12) * std::sin(Scalar(k) * Scalar(0.7));
        // Vary the base translation a little so the rows are not perfectly collinear.
        const Vec3 t_A(Scalar(0.8) + Scalar(0.1) * std::cos(Scalar(k) * 0.5),
                       Scalar(0.03) * std::sin(Scalar(k) * 0.3), Scalar(-0.01));
        SE3 A; A.R = so3::exp(axis.normalized() * theta); A.t = t_A;
        const Vec3 omega = so3::log(A.R) / dt;
        SE3 B = trans_only_window(X, t_X, axis, theta, t_A, sensor_noise,
                                  sensor_noise.squaredNorm() > 0 ? &rng : nullptr);
        SourceId ids[1] = { id };
        SE3      rep[1] = { B };
        cal.observe(1, ids, rep, A, omega);
    }
}

// A HEADING-BLIND radar ISource for the ESTIMATOR-level path. Over a query window it reads
// the TRUE base delta A from the trajectory and reports R_B = I (no rotation) + the
// lever-bearing sensor-frame displacement t_B = R_X^T [ t_A + (R_A - I) t_X ] (the physical
// Doppler model: v_sensor = v_ego + omega x lever). The mount (R_X = X.R) is the TRUSTED
// prior; t_X = X.t the lever. Optional per-axis SENSOR-frame noise on t_B. This is the
// real measurement a translation_only source produces — the estimator pins R_X = prior.
class RadarSource : public ISource {
public:
    RadarSource(SourceId id, const Trajectory& traj, const SE3& X,
                const Vec3& sensor_noise = Vec3::Zero(), std::uint32_t seed = 0)
        : id_(id), traj_(&traj), X_(X), noise_(sensor_noise), seed_(seed) {}
    SourceId id() const override { return id_; }
    Expected<Delta> query(Timestamp t0, Timestamp t1) const override {
        if (t1 <= t0) return Status::OutOfRange;
        const SE3 p0 = traj_->pose(t0), p1 = traj_->pose(t1);
        SE3 A; A.R = p0.R.transpose() * p1.R; A.t = p0.R.transpose() * (p1.t - p0.t);
        Vec3 t_B = X_.R.transpose() * (A.t + (A.R - Mat3::Identity()) * X_.t);
        if (noise_.squaredNorm() > Scalar(0)) {
            std::mt19937_64 g(seed_ ^ (static_cast<std::uint64_t>(t0) * 0x9E3779B97F4A7C15ull)
                              ^ static_cast<std::uint64_t>(t1));
            std::normal_distribution<Scalar> nx(0, noise_.x()), ny(0, noise_.y()), nz(0, noise_.z());
            const Scalar e0 = nx(g), e1 = ny(g), e2 = nz(g);
            t_B += Vec3(e0, e1, e2);
        }
        Delta d;
        d.t0 = t0; d.t1 = t1;
        d.motion.R = Mat3::Identity();      // HEADING-BLIND: no rotation reported
        d.motion.t = t_B;
        d.cov = Mat6::Identity() * Scalar(1e-2);
        return d;
    }
private:
    SourceId          id_;
    const Trajectory* traj_;
    SE3               X_;
    Vec3              noise_;
    std::uint32_t     seed_;
};

const CalibSnapshot* snap(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.calib[i].id == id) return &r.calib[i];
    return nullptr;
}

// A MULTI-AXIS turning trajectory (the lever needs rotation excitation; multi-axis so all
// 3 lever axes — incl. z — are observable) with a few straight segments (a well-formed
// reference-driven consensus). yaw + pitch turns give two distinct rotation axes.
Trajectory turn_rich_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0,     0;
    Vec6 yawL;     yawL     << 2.0, 0, 0, 0, 0,     0.6;
    Vec6 yawR;     yawR     << 2.0, 0, 0, 0, 0,    -0.6;
    Vec6 pitch;    pitch    << 2.0, 0, 0, 0, 0.5,   0.0;
    for (int rep = 0; rep < 12; ++rep) {
        tr.add_segment(straight, 0.6);
        tr.add_segment(yawL, 1.4);
        tr.add_segment(pitch, 1.2);
        tr.add_segment(yawR, 1.4);
    }
    return tr;
}

void set_hists(Config& c) {
    c.vote_weight = VoteWeight::One;
    c.straight_omega_max = 0.05; c.straight_trans_min = 0.02; c.turn_omega_min = 0.20;
    c.so3_hist.bins = 256; c.so3_hist.range_min = -0.8; c.so3_hist.range_max = 0.8;
    c.so3_hist.aging = Aging::SlidingK; c.so3_hist.sliding_k = 200;
    c.scale_hist.bins = 256; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 200;
    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 200;
    c.xyz_hist.bins = 256; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 200;
    c.offset_hist.bins = 256; c.offset_hist.range_min = -0.1; c.offset_hist.range_max = 0.1;
    c.offset_hist.aging = Aging::SlidingK; c.offset_hist.sliding_k = 200;
}
} // namespace

// ===========================================================================
// 1. HEADLINE — flag OFF garbage rotation; flag ON lever recovers + rotation = prior
// ===========================================================================
TEST_CASE("trans-only headline: flag OFF garbage rotation, flag ON lever recovers + rotation pinned") {
    // A heading-blind radar mount: yaw 25 deg, pitch 0, roll 0; lever [0.6, -0.3, 0.0].
    // The z lever is 0 (planar null on the multi-axis set z IS observable, so plant a true
    // value too to prove z recovers when observable).
    const SE3 X = make_extrinsic(25 * kPi / 180, 0, 0, Vec3(0.6, -0.3, 0.10));
    const Vec3 t_X = X.t;

    // --- flag OFF (the existing path): R_B = I -> rot3d never opens; Phase-2 roll fits the
    //     absent rotation. The lever rows use R_X = R_yp o Rx(roll) with R_yp = prior (no
    //     Phase-1 here), and a wrong roll pollutes the rows. We check the OFF path does NOT
    //     recover the lever (the failure this slice fixes). ---
    {
        auto calp = make_p2(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_p2_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, X, /*scale_calib=*/true, /*translation_only=*/false) == Status::Ok);
        // No set_yaw_pitch fed (bootstrap path) -> R_yp falls back to yaw_pitch_of(prior).
        feed_trans_only(cal, 1, X, t_X, 240, /*multiaxis=*/true);
        const Vec3 lev = cal.lever_arm(1);
        const Scalar err = (lev - t_X).norm();
        INFO("flag OFF lever err = " << err << " m  (recovered " << lev.transpose() << ")");
        // The OFF path does NOT cleanly recover the lever (roll fit of an absent rotation
        // corrupts the rows). Pin that it is materially worse than the ON path below.
        CHECK(err > Scalar(0.02));
    }

    // --- flag ON: rotation pinned to prior, lever recovers ---
    {
        auto calp = make_p2(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_p2_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, X, /*scale_calib=*/true, /*translation_only=*/true) == Status::Ok);
        feed_trans_only(cal, 1, X, t_X, 240, /*multiaxis=*/true);

        // Lever recovers (all 3 axes observable on the multi-axis set).
        bool obs = false;
        const Vec3 lev_ls = cal.solve_lever_arm(1, &obs);
        CHECK(obs);
        const Vec3 lev = cal.lever_arm(1);
        const Scalar err = (lev - t_X).norm();
        INFO("flag ON lever err = " << err << " m  (recovered " << lev.transpose() << ")");
        CHECK(err < Scalar(2e-3));
        CHECK(cal.translation_confidence(1) > Scalar(0.0));

        // Rotation extrinsic stays = prior (NOT garbage). roll never voted; rot3d empty.
        const SE3 Xr = cal.extrinsic(1);
        const Scalar rot_off = so3::log(Xr.R * X.R.transpose()).norm();
        INFO("flag ON rotation deviation from prior = " << rot_off << " rad");
        CHECK(rot_off < Scalar(1e-12));
        CHECK(cal.roll_vote_count(1) == Scalar(0));        // roll never voted
        CHECK(cal.rot3d_vote_count(1) == Scalar(0));       // rot3d never voted
        CHECK(cal.extrinsic_confidence(1) == Scalar(0));   // roll concentration 0
    }
}

// ===========================================================================
// 2. FOOTGUN — even with rot3d ENABLED, a translation_only source never votes a rotation
// ===========================================================================
TEST_CASE("trans-only footgun: rot3d enabled but translation_only never votes a rotation") {
    const SE3 X = make_extrinsic(15 * kPi / 180, -8 * kPi / 180, 5 * kPi / 180,
                                 Vec3(0.4, 0.2, -0.1));
    auto calp = make_p2(); Phase2Calibrator& cal = *calp;
    // rot3d ON — the footgun path (a heading-blind source could try to vote a full rotation).
    REQUIRE(cal.configure(make_p2_cfg(/*rot3d=*/true), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, X, true, /*translation_only=*/true) == Status::Ok);
    feed_trans_only(cal, 1, X, X.t, 240, /*multiaxis=*/true);

    // The rot3d gate must NEVER open and NO rotation vote may land.
    CHECK_FALSE(cal.rot3d_observable(1));
    CHECK(cal.rot3d_vote_count(1) == Scalar(0));
    CHECK(cal.rot3d_confidence(1) == Scalar(0));
    CHECK(cal.roll_vote_count(1) == Scalar(0));
    // rot3d() falls back to the basepoint (= prior rotation): no garbage rotation.
    const Scalar dev = so3::log(cal.rot3d(1) * X.R.transpose()).norm();
    CHECK(dev < Scalar(1e-12));
    // Lever still recovers with the clean R_X = prior.
    CHECK((cal.lever_arm(1) - X.t).norm() < Scalar(2e-3));
}

// ===========================================================================
// 3. ALONG-TRACK / NOISY AXIS — the existing concentration gate measurement
// ===========================================================================
// MEASURE (SLICE20 section "THE ALONG-TRACK GUARD"): a source clean on the LATERAL (y)
// axis but carrying along-track (x) sensor drift. Does the noisy x axis COMMIT a wrong
// value, or does the existing histogram-concentration gate already withhold it?
TEST_CASE("trans-only along-track: concentration gate already withholds the noisy axis") {
    // Yaw-only base motion (planar) so x and y are the two observable lever axes and z is
    // the null. Mount yaw 20 deg; lever [0.7 (along-track), 0.25 (lateral), 0].
    const SE3 X = make_extrinsic(20 * kPi / 180, 0, 0, Vec3(0.7, 0.25, 0.0));
    const Vec3 t_X = X.t;
    const Config cfg = make_p2_cfg();

    // CLEAN control: NO sensor noise -> both observable axes concentrate, the joint
    // translation confidence rises ABOVE the commit gate (a wrong value would commit — so
    // the gate is NOT a blanket "never commit", it is selective).
    {
        auto calp = make_p2(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(cfg, 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, X, true, true) == Status::Ok);
        feed_trans_only(cal, 1, X, t_X, 400, /*multiaxis=*/false);
        // x,y recover; z is the planar null (withheld by the conditioning info gate, so it
        // falls back to the prior = truth here). Slice-20 min-over-OBSERVED translation
        // confidence EXCLUDES the empty z, so the OBSERVABLE axes' clean concentration is
        // reported -> the planar lever COMMITS (the strict min-over-3 would have been pinned
        // at 0 by the empty z and never committed -- the fix's whole point).
        CHECK(std::abs(cal.lever_arm(1).x() - t_X.x()) < Scalar(2e-3));
        CHECK(std::abs(cal.lever_arm(1).y() - t_X.y()) < Scalar(2e-3));
        INFO("clean trans_conf = " << cal.translation_confidence(1));
        CHECK(cal.translation_confidence(1) >= cfg.commit_concentration);
    }

    // NOISY along-track: heavy SENSOR-frame drift on the along-track (sensor x) axis (the
    // ~0.5 m/window Doppler drift the prototype found), ~clean on lateral. The MEASUREMENT
    // the slice's along-track guard hinges on: does the noisy along-track axis COMMIT a
    // wrong value, or does the existing histogram-CONCENTRATION gate already withhold it?
    {
        auto calp = make_p2(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(cfg, 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, X, true, true) == Status::Ok);
        const Vec3 sensor_noise(Scalar(0.5), Scalar(0.01), Scalar(0.0));
        feed_trans_only(cal, 1, X, t_X, 400, /*multiaxis=*/false, sensor_noise);

        const Vec3 lev = cal.lever_arm(1);
        const Scalar ex = std::abs(lev.x() - t_X.x());
        const Scalar ey = std::abs(lev.y() - t_X.y());
        const Scalar tconf = cal.translation_confidence(1);
        INFO("noisy along-track x err = " << ex << " m ; lateral y err = " << ey
             << " m ; trans_conf = " << tconf);
        // THE FINDING (Slice-20 min-over-OBSERVED, the HONEST version): the noisy along-track
        // x IS an observed axis (well-conditioned, gets votes) but its scatter destroys its
        // per-channel concentration cx. min-over-OBSERVED = min(cx, cy) is therefore dragged
        // down by the noisy cx (NOT by the empty z — z is now correctly excluded), so the
        // JOINT translation_confidence collapses below the commit gate and the lever is
        // WITHHELD. (Confirms: a noisy-but-OBSERVED along-track axis still blocks the joint
        // commit -> committing only the clean lateral axis would need PER-AXIS commit, a
        // separate change; the min-over-observed fix only unpins the empty unobservable axis.)
        CHECK(tconf < cfg.commit_concentration);
    }
}

// ===========================================================================
// 4. z / PLANAR NULL — yaw-only turns leave z unobservable, x,y observable
// ===========================================================================
TEST_CASE("trans-only planar null: yaw-only -> z withheld, x,y recover") {
    const SE3 X = make_extrinsic(10 * kPi / 180, 0, 0, Vec3(0.5, -0.2, 0.33));
    const Vec3 t_X = X.t;
    auto calp = make_p2(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_p2_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, X, true, /*translation_only=*/true) == Status::Ok);
    feed_trans_only(cal, 1, X, t_X, 300, /*multiaxis=*/false);   // yaw-only

    // x,y observable -> recovered tightly; z stays at the prior (the conditioning guard
    // pins the null direction). solve_lever_arm reports the conditioning observability.
    bool obs = false;
    cal.solve_lever_arm(1, &obs);
    CHECK_FALSE(obs);                  // the 3x3 is rank-deficient (z null) -> not observable
    const Vec3 lev = cal.lever_arm(1);
    CHECK(std::abs(lev.x() - t_X.x()) < Scalar(2e-3));
    CHECK(std::abs(lev.y() - t_X.y()) < Scalar(2e-3));
    // z stays at the prior (= truth here, but the point is it is NOT driven off by a blown
    // solve): pin that the z channel never concentrated.
    CHECK(std::abs(lev.z() - t_X.z()) < Scalar(2e-3));   // == prior (unvoted) by construction
}

// ===========================================================================
// 6. SCALE preserved — a translation_only source still recovers its straight-regime scale
// ===========================================================================
// Phase-1 keeps the scale vote for a translation_only source (the straight-regime magnitude
// ratio is valid for a velocity sensor) while skipping the direction (yaw/pitch) vote.
TEST_CASE("trans-only scale preserved: Phase-1 scale votes, direction skipped") {
    Config c = make_p2_cfg();
    c.straight_omega_max = 0.05; c.straight_trans_min = 0.02;
    auto p1p = make_p1(); Phase1Calibrator& p1 = *p1p;
    REQUIRE(p1.configure(c, 0) == Status::Ok);
    // Reference (id 0) + a translation_only source (id 1) with a true residual scale 1.2.
    REQUIRE(p1.set_prior(0, SE3{}, /*scale_calib=*/false, /*translation_only=*/false) == Status::Ok);
    const SE3 X1 = make_extrinsic(30 * kPi / 180, 0, 0, Vec3(0.5, 0.1, 0.0));
    REQUIRE(p1.set_prior(1, X1, /*scale_calib=*/true, /*translation_only=*/true) == Status::Ok);

    const Scalar true_scale = Scalar(1.2);
    // Straight forward windows: A.R = I, A.t = d * e_x. Reference reports A directly; the
    // translation_only source reports B = X1^-1 o A o X1 then *scale (R_B carries the mount
    // rotation but is IGNORED by the direction skip; only |t_B| feeds the scale ratio).
    for (int k = 0; k < 200; ++k) {
        const Scalar d = Scalar(1.0) + Scalar(0.05) * std::sin(Scalar(k) * 0.3);
        SE3 A; A.R = Mat3::Identity(); A.t = Vec3(d, 0, 0);
        SE3 Bref = A;                                   // reference: identity mount
        SE3 B1 = se3::compose(se3::compose(se3::inverse(X1), A), X1);
        B1.t *= true_scale;
        SourceId ids[2] = { 0, 1 };
        SE3      rep[2] = { Bref, B1 };
        const Vec3 omega = Vec3::Zero();                // straight
        const Vec3 trans = A.t;
        p1.observe(2, ids, rep, omega, trans);
    }

    // Scale recovers; the DIRECTION (yaw/pitch) is SKIPPED -> so(3) empty -> extrinsic = prior.
    const Scalar s = p1.scale(1);
    INFO("trans-only recovered scale = " << s << " (true " << true_scale << ")");
    CHECK(std::abs(s - true_scale) < Scalar(0.02));
    CHECK(p1.scale_confidence(1) > Scalar(0.5));
    CHECK(p1.vote_count(1) == Scalar(0));              // so(3) direction never voted
    CHECK(p1.extrinsic_confidence(1) == Scalar(0));
    const Scalar rot_off = so3::log(p1.extrinsic(1).R * X1.R.transpose()).norm();
    CHECK(rot_off < Scalar(1e-12));                   // rotation = prior, not recovered
}

// ===========================================================================
// 1b/2b. ESTIMATOR-LEVEL HEADLINE + FOOTGUN — the full feedback loop
// ===========================================================================
// A heading-blind RadarSource in a multi-source rig: the estimator pins its rotation
// extrinsic to the prior and recovers its lever via the feedback loop. Pins (1) the lever
// commits/recovers, AND (2) the published rotation extrinsic NEVER moves off the prior and
// extrinsic_committed stays false (the footgun: today it would commit a garbage rotation).
TEST_CASE("trans-only estimator: radar lever recovers, rotation extrinsic never commits off prior") {
    const Trajectory tr = turn_rich_traj();

    // Rig: reference (id 0, identity mount) + 2 normal odom sources (ids 1,2, small mounts)
    // forming a good turning consensus + the radar (id 3, heading-blind, yaw-30 mount,
    // lever [0.8, -0.35, 0]). z lever 0 (the planar null on yaw-only turns).
    const SE3 X1 = make_extrinsic(3 * kPi / 180, 0, 0, Vec3(0.1, 0.05, 0.0));
    const SE3 X2 = make_extrinsic(-4 * kPi / 180, 2 * kPi / 180, 0, Vec3(-0.1, 0.0, 0.05));
    // Radar mount: yaw-30 mount, lever [0.8, -0.35, 0.12]. The multi-axis (yaw+pitch)
    // trajectory makes all 3 lever axes (incl. z) observable.
    const SE3 Xr = make_extrinsic(30 * kPi / 180, 0, 0, Vec3(0.8, -0.35, 0.12));

    SourceParams p0; p0.id = 0;                       // reference, identity
    SourceParams p1; p1.id = 1; p1.X = X1;
    SourceParams p2; p2.id = 2; p2.X = X2;
    SyntheticSource s0(tr, p0), s1(tr, p1), s2(tr, p2);
    RadarSource radar(3, tr, Xr);

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[0].is_reference = true;
    sensors[1].prior_extrinsic = X1;
    sensors[2].prior_extrinsic = X2;
    sensors[3].prior_extrinsic = Xr;                  // the TRUSTED radar mount (prior)
    sensors[3].translation_only = true;               // <-- the slice flag
    sensors[3].scale_calib = false;                   // radar: scale pinned (focus on lever)

    Config cfg; set_hists(cfg);
    cfg.reference_sensor_id = 0;
    cfg.cold_start = ColdStart::MedianFromStart;
    cfg.commit_concentration = 0.5; cfg.commit_min_votes = 30; cfg.commit_drop = 0.2;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;
    REQUIRE(validate(cfg) == Status::Ok);

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    REQUIRE(est.add_source(&s0) == Status::Ok);
    REQUIRE(est.add_source(&s1) == Status::Ok);
    REQUIRE(est.add_source(&s2) == Status::Ok);
    REQUIRE(est.add_source(&radar) == Status::Ok);

    const Timestamp step_ns = static_cast<Timestamp>(std::llround(1e9 / 50.0));
    const Timestamp t_end = static_cast<Timestamp>(std::llround((tr.duration_s() - 0.1) * 1e9));
    for (Timestamp now = static_cast<Timestamp>(std::llround(0.2 * 1e9)); now < t_end; now += step_ns) {
        est.step(now);
    }

    const CalibSnapshot* cs = snap(est.latest(), 3);
    REQUIRE(cs != nullptr);

    // (1) The radar LEVER recovered: all 3 axes near truth (multi-axis turns make z
    //     observable too). The lever DOF committed.
    INFO("radar lever recovered = " << cs->extrinsic.t.transpose()
         << "  truth = " << Xr.t.transpose());
    CHECK((cs->extrinsic.t - Xr.t).norm() < Scalar(0.04));
    CHECK(cs->translation_committed);

    // (2) THE FOOTGUN FIX: the published rotation extrinsic NEVER moved off the prior, and
    //     extrinsic_committed is FALSE (trusted, not recovered). The nuScenes rx=-0.984
    //     case must not occur.
    const Scalar rot_dev = so3::log(cs->extrinsic.R * Xr.R.transpose()).norm();
    INFO("radar rotation deviation from prior = " << rot_dev << " rad");
    CHECK(rot_dev < Scalar(1e-9));
    CHECK_FALSE(cs->extrinsic_committed);
}

// ===========================================================================
// 7. validate() + config-hash flip
// ===========================================================================
TEST_CASE("trans-only validate: requires an orthonormal prior rotation; hash flips") {
    std::vector<SensorConfig> sensors(2);
    sensors[0].id = 0; sensors[0].is_reference = true;
    sensors[1].id = 1;
    sensors[1].translation_only = true;
    sensors[1].prior_extrinsic = make_extrinsic(0.3, -0.1, 0.05, Vec3(0.5, 0, 0));

    Config c;
    c.reference_sensor_id = 0;
    c.sensors = sensors.data();
    c.sensor_count = 2;

    // A valid orthonormal prior rotation validates Ok.
    CHECK(validate(c) == Status::Ok);

    // A GARBAGE (non-orthonormal) prior rotation on a translation_only source is rejected
    // (its rotation is trusted, not recovered — a bad prior would corrupt every lever row).
    SUBCASE("garbage prior rotation rejected") {
        sensors[1].prior_extrinsic.R = Mat3::Identity() * Scalar(2.0);   // not SO(3)
        CHECK(validate(c) == Status::InvalidConfig);
    }
    // A NON-translation_only source with the same garbage rotation still validates (the
    // check is gated on the flag — the recovered-rotation path tolerates a rough prior).
    SUBCASE("garbage prior tolerated when flag off") {
        sensors[1].translation_only = false;
        sensors[1].prior_extrinsic.R = Mat3::Identity() * Scalar(2.0);
        CHECK(validate(c) == Status::Ok);
    }
    // The config-hash includes the flag: flipping it on one source changes the hash, so a
    // blob written under translation_only=true REJECTS when restored into the otherwise-
    // identical translation_only=false rig (pre-flag blobs cold-start by design).
    SUBCASE("config-hash flips with the flag") {
        Estimator est_on;
        REQUIRE(est_on.init(c) == Status::Ok);               // flag ON

        std::vector<SensorConfig> sensors_off = sensors;
        sensors_off[1].translation_only = false;
        Config c_off = c;
        c_off.sensors = sensors_off.data();
        Estimator est_off;
        REQUIRE(est_off.init(c_off) == Status::Ok);          // flag OFF, same rig otherwise

        std::vector<unsigned char> blob(8192);
        const Expected<int> w = est_on.serialize(blob.data(), static_cast<int>(blob.size()));
        REQUIRE(w.ok());
        blob.resize(w.value());
        // Restoring the flag-ON blob into the flag-OFF rig must reject on the config-hash.
        CHECK(est_off.deserialize(blob.data(), static_cast<int>(blob.size()))
              == Status::InvalidConfig);
    }
}

// ===========================================================================
// 8. OBSERVABILITY SELF-TEST — a translation_only case converges in its regime
// ===========================================================================
// Each DOF converges only where it is observable: the lever needs rotation excitation
// (turn windows) and a non-degenerate axis set; the rotation extrinsic is pinned (never
// recovered). Mirrors the existing per-calibrator observability self-tests.
TEST_CASE("trans-only observability: lever needs turns; pure-straight starves; multi-axis full") {
    const SE3 X = make_extrinsic(12 * kPi / 180, 4 * kPi / 180, 0, Vec3(0.55, -0.15, 0.22));
    const Vec3 t_X = X.t;

    // (a) NO turn windows (theta below the turn gate) -> the lever never accumulates.
    {
        auto calp = make_p2(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_p2_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, X, true, true) == Status::Ok);
        const Scalar dt = Scalar(0.5);
        for (int k = 0; k < 100; ++k) {
            SE3 A; A.R = so3::exp(Vec3(0, 0, 1) * Scalar(0.02)); A.t = Vec3(1.0, 0, 0);
            const Vec3 omega = so3::log(A.R) / dt;   // below turn_omega_min
            SE3 B = trans_only_window(X, t_X, Vec3(0, 0, 1), Scalar(0.02), A.t);
            SourceId ids[1] = { 1 }; SE3 rep[1] = { B };
            cal.observe(1, ids, rep, A, omega);
        }
        CHECK(cal.xyz_vote_count(1) == Scalar(0));     // no turn -> no lever rows
        CHECK((cal.lever_arm(1) - X.t).norm() < Scalar(1e-12));  // == prior (unvoted)
    }

    // (b) multi-axis turns -> all 3 lever axes observable, converge.
    {
        auto calp = make_p2(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_p2_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, X, true, true) == Status::Ok);
        feed_trans_only(cal, 1, X, t_X, 300, /*multiaxis=*/true);
        bool obs = false;
        cal.solve_lever_arm(1, &obs);
        CHECK(obs);
        CHECK((cal.lever_arm(1) - t_X).norm() < Scalar(2e-3));
    }
}

// ===========================================================================
// Slice 20b. PER-AXIS commit machinery
// ===========================================================================

// A PURE yaw-only (planar) turning trajectory: x = along-track, y = lateral, z = the planar
// null. Straight segments give a clean reference-driven consensus; the yaw turns excite the
// lever's x,y axes. No pitch (so z stays unobservable, the null axis stays at the prior).
Trajectory yaw_only_traj() {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0,  0;
    Vec6 yawL;     yawL     << 2.0, 0, 0, 0, 0,  0.6;
    Vec6 yawR;     yawR     << 2.0, 0, 0, 0, 0, -0.6;
    for (int rep = 0; rep < 16; ++rep) {
        tr.add_segment(straight, 0.6);
        tr.add_segment(yawL, 1.6);
        tr.add_segment(straight, 0.6);
        tr.add_segment(yawR, 1.6);
    }
    return tr;
}

// 9. PER-AXIS HEADLINE (sim): a translation_only source on yaw-only planar motion, one lever
//    axis CLEAN and one with heavy sensor drift (z = the planar null). The CLEAN axis commits
//    while the NOISY axis + null z are each WITHHELD -- the capability the whole-lever min
//    could not deliver. CONTRAST: the whole-lever min-over-observed confidence is BELOW the
//    commit gate (the noisy axis drags the min down -> the whole-lever gate would NOT have
//    committed), yet the per-axis clean-axis commit fires.
//
// AXIS NOTE: the lever LS couples the recovered xyz channels through the yaw rotation + the
// mount, so a sensor-frame noise direction does NOT map 1:1 to a recovered channel. With this
// mount (yaw 20) + this geometry the sensor-x Doppler drift leaves the RECOVERED x channel
// concentrated (clean) and scatters the RECOVERED y channel (measured below: cx ~ 0.61,
// cy ~ 0.17). So the CLEAN committed axis here is x, the WITHHELD noisy axis is y. The slice
// capability (one axis commits, the other is held) is what is pinned -- not a fixed label.
TEST_CASE("trans-only 20b per-axis headline: clean axis commits, noisy axis + null z withheld") {
    const SE3 X = make_extrinsic(20 * kPi / 180, 0, 0, Vec3(0.7, 0.25, 0.0));
    const Vec3 t_X = X.t;
    const Config cfg = make_p2_cfg();
    const Scalar gate = cfg.commit_concentration;

    // --- (A) CALIBRATOR-LEVEL per-axis separation (deterministic; pins the accessors) ----
    // The same noisy along-track rig as the existing case-3 measurement. Sensor-x drift 0.5 m.
    auto calxp = make_p2(); Phase2Calibrator& calx = *calxp;
    REQUIRE(calx.configure(cfg, 0) == Status::Ok);
    REQUIRE(calx.set_prior(1, X, /*scale_calib=*/true, /*translation_only=*/true) == Status::Ok);
    feed_trans_only(calx, 1, X, t_X, 400, /*multiaxis=*/false, Vec3(Scalar(0.5), Scalar(0.01), 0));

    const Scalar cx = calx.translation_confidence_axis(1, 0);
    const Scalar cy = calx.translation_confidence_axis(1, 1);
    const Scalar cz = calx.translation_confidence_axis(1, 2);
    INFO("per-axis conf: cx=" << cx << " cy=" << cy << " cz=" << cz
         << "  whole=" << calx.translation_confidence(1) << "  gate=" << gate);
    // The CLEAN axis (x) clears the gate; the NOISY axis (y) does not; z is the empty null.
    CHECK(cx >= gate);                          // clean axis -> per-axis commit eligible
    CHECK(cy <  gate);                          // noisy axis -> per-axis WITHHELD
    CHECK(cz == Scalar(0));                      // planar null z -> empty
    // The accessors expose REAL vote mass on the observed axes, none on the null z.
    CHECK(calx.xyz_axis_mass(1, 0) > Scalar(0));
    CHECK(calx.xyz_axis_mass(1, 1) > Scalar(0));
    CHECK(calx.xyz_axis_mass(1, 2) == Scalar(0));
    // Out-of-range axis -> 0 (defensive).
    CHECK(calx.translation_confidence_axis(1, 3) == Scalar(0));
    CHECK(calx.xyz_axis_mass(1, -1) == Scalar(0));
    // The recovered CLEAN-axis lever is near truth (x err < 2 cm); the recovered value lives
    // in lever_arm()[0]. (The noisy y is scattered -- it is exactly why it is withheld.)
    CHECK(std::abs(calx.lever_arm(1).x() - t_X.x()) < Scalar(0.02));
    // THE CONTRAST: the whole-lever min-over-observed confidence is BELOW the gate (cy drags
    // it down) -> the whole-lever gate would NOT have committed, yet the per-axis x clears it.
    CHECK(calx.translation_confidence(1) < gate);
    CHECK(cx >= gate);

    // --- (B) ESTIMATOR-LEVEL per-axis commit + publish (the feedback loop) ----------------
    // Drive the full estimator with the same kind of source; assert the per-axis COMMIT flags
    // split (one axis committed, the other withheld) and the whole flag is the OR. The per-axis
    // PUBLISH into the fusion prior (committed axis folded, withheld axis held at prior) is
    // pinned via the persistence blob in the 20b persistence test (the blob carries the real
    // prior_extrinsic; the snapshot's extrinsic.t is the calibrator lever readout, not the
    // fusion prior). The estimator path is noisier than the calibrator, so use a drift level
    // that splits the commit cleanly (sensor-x 0.3 m).
    const Trajectory tr = yaw_only_traj();
    const SE3 X1 = make_extrinsic(3 * kPi / 180, 0, 0, Vec3(0.1, 0.05, 0.0));
    const SE3 X2 = make_extrinsic(-4 * kPi / 180, 0, 0, Vec3(-0.1, 0.0, 0.05));
    SourceParams p0; p0.id = 0;
    SourceParams p1; p1.id = 1; p1.X = X1;
    SourceParams p2; p2.id = 2; p2.X = X2;
    SyntheticSource s0(tr, p0), s1(tr, p1), s2(tr, p2);
    RadarSource radar(3, tr, X, Vec3(Scalar(0.3), Scalar(0.002), Scalar(0.0)), /*seed=*/4242);

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[0].is_reference = true;
    sensors[1].prior_extrinsic = X1;
    sensors[2].prior_extrinsic = X2;
    sensors[3].prior_extrinsic = X;
    sensors[3].translation_only = true;
    sensors[3].scale_calib = false;

    Config ecfg; set_hists(ecfg);
    ecfg.reference_sensor_id = 0;
    ecfg.cold_start = ColdStart::MedianFromStart;
    ecfg.commit_concentration = 0.5; ecfg.commit_min_votes = 30; ecfg.commit_drop = 0.2;
    ecfg.sensors = sensors.data(); ecfg.sensor_count = 4;
    REQUIRE(validate(ecfg) == Status::Ok);

    Estimator est;
    REQUIRE(est.init(ecfg) == Status::Ok);
    REQUIRE(est.add_source(&s0) == Status::Ok);
    REQUIRE(est.add_source(&s1) == Status::Ok);
    REQUIRE(est.add_source(&s2) == Status::Ok);
    REQUIRE(est.add_source(&radar) == Status::Ok);

    const Timestamp step_ns = static_cast<Timestamp>(std::llround(1e9 / 50.0));
    const Timestamp t_end = static_cast<Timestamp>(std::llround((tr.duration_s() - 0.1) * 1e9));
    for (Timestamp now = static_cast<Timestamp>(std::llround(0.2 * 1e9)); now < t_end; now += step_ns) {
        est.step(now);
    }

    const CalibSnapshot* cs = snap(est.latest(), 3);
    REQUIRE(cs != nullptr);
    INFO("estimator per-axis flags = {" << cs->translation_committed_xyz[0] << ","
         << cs->translation_committed_xyz[1] << "," << cs->translation_committed_xyz[2] << "}");
    // The per-axis commit SPLIT: exactly the clean axis (x) committed; the noisy (y) + null z
    // are withheld. The whole flag is the OR.
    CHECK(cs->translation_committed_xyz[0]);          // clean axis committed
    CHECK_FALSE(cs->translation_committed_xyz[1]);    // noisy axis withheld
    CHECK_FALSE(cs->translation_committed_xyz[2]);    // null z withheld
    CHECK(cs->translation_committed);                  // whole flag = OR = true
    // The whole-lever min could NOT have committed (the noisy axis keeps the min under gate),
    // yet a per-axis commit fired -> the 20b unlock.
    CHECK(cs->translation_confidence < ecfg.commit_concentration);
    // The footgun stays fixed: rotation extrinsic never moved off the prior, not committed.
    CHECK(so3::log(cs->extrinsic.R * X.R.transpose()).norm() < Scalar(1e-9));
    CHECK_FALSE(cs->extrinsic_committed);
}

// 10. NON-FLAGGED BYTE-IDENTICAL: a NORMAL (non-translation_only) source's lever commit +
//     published prior_extrinsic.t + per-axis snapshot flags follow the WHOLE-lever path
//     exactly. The per-axis array mirrors the whole flag (all 3 == translation_committed),
//     and the whole-lever publish is unchanged. (The wider byte-identical guarantee is the
//     untouched existing calib + golden suites; here we pin the per-axis mirror invariant.)
TEST_CASE("trans-only 20b non-flagged: per-axis array mirrors the whole-lever flag") {
    const Trajectory tr = turn_rich_traj();

    // A normal rig: reference + 3 NORMAL odom sources (NO translation_only anywhere).
    const SE3 X1 = make_extrinsic(3 * kPi / 180, 0, 0, Vec3(0.1, 0.05, 0.0));
    const SE3 X2 = make_extrinsic(-4 * kPi / 180, 2 * kPi / 180, 0, Vec3(-0.1, 0.0, 0.05));
    const SE3 X3 = make_extrinsic(6 * kPi / 180, -3 * kPi / 180, 0, Vec3(0.2, -0.15, 0.08));

    SourceParams p0; p0.id = 0;
    SourceParams p1; p1.id = 1; p1.X = X1;
    SourceParams p2; p2.id = 2; p2.X = X2;
    SourceParams p3; p3.id = 3; p3.X = X3;
    SyntheticSource s0(tr, p0), s1(tr, p1), s2(tr, p2), s3(tr, p3);

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[0].is_reference = true;
    sensors[1].prior_extrinsic = X1;
    sensors[2].prior_extrinsic = X2;
    sensors[3].prior_extrinsic = X3;
    sensors[3].scale_calib = false;
    // NO translation_only flag set on any source.

    Config cfg; set_hists(cfg);
    cfg.reference_sensor_id = 0;
    cfg.cold_start = ColdStart::MedianFromStart;
    cfg.commit_concentration = 0.5; cfg.commit_min_votes = 30; cfg.commit_drop = 0.2;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;
    REQUIRE(validate(cfg) == Status::Ok);

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    REQUIRE(est.add_source(&s0) == Status::Ok);
    REQUIRE(est.add_source(&s1) == Status::Ok);
    REQUIRE(est.add_source(&s2) == Status::Ok);
    REQUIRE(est.add_source(&s3) == Status::Ok);

    const Timestamp step_ns = static_cast<Timestamp>(std::llround(1e9 / 50.0));
    const Timestamp t_end = static_cast<Timestamp>(std::llround((tr.duration_s() - 0.1) * 1e9));
    // Walk EVERY published snapshot: for a non-flagged source the per-axis array must, at
    // every step, equal {whole, whole, whole} -- the whole-lever path drives all three.
    for (Timestamp now = static_cast<Timestamp>(std::llround(0.2 * 1e9)); now < t_end; now += step_ns) {
        est.step(now);
        const Result& r = est.latest();
        for (int k = 0; k < r.source_count; ++k) {
            const CalibSnapshot& c = r.calib[k];
            const bool whole = c.translation_committed;
            CHECK(c.translation_committed_xyz[0] == whole);
            CHECK(c.translation_committed_xyz[1] == whole);
            CHECK(c.translation_committed_xyz[2] == whole);
        }
    }

    // Source 3's lever committed via the WHOLE-lever path (the multi-axis turns observe all 3).
    const CalibSnapshot* cs = snap(est.latest(), 3);
    REQUIRE(cs != nullptr);
    CHECK(cs->translation_committed);
    CHECK(cs->translation_committed_xyz[0]);
    CHECK(cs->translation_committed_xyz[1]);
    CHECK(cs->translation_committed_xyz[2]);
    CHECK((cs->extrinsic.t - X3.t).norm() < Scalar(0.04));
}

// 11. PERSISTENCE: serialize a rig with a translation_only source MID per-axis commit (lateral
//     committed, along-track not) -> deserialize into a fresh estimator -> the per-axis flags +
//     the per-axis folded prior round-trip; a pre-v4 (stamped v3) blob cold-starts (rejected by
//     the version guard).
TEST_CASE("trans-only 20b persistence: per-axis flags + folded prior round-trip; pre-v4 rejects") {
    const Trajectory tr = yaw_only_traj();

    const SE3 X1 = make_extrinsic(3 * kPi / 180, 0, 0, Vec3(0.1, 0.05, 0.0));
    const SE3 X2 = make_extrinsic(-4 * kPi / 180, 0, 0, Vec3(-0.1, 0.0, 0.05));
    const SE3 Xr = make_extrinsic(20 * kPi / 180, 0, 0, Vec3(0.7, 0.25, 0.0));

    SourceParams p0; p0.id = 0;
    SourceParams p1; p1.id = 1; p1.X = X1;
    SourceParams p2; p2.id = 2; p2.X = X2;
    SyntheticSource s0(tr, p0), s1(tr, p1), s2(tr, p2);
    // Sensor-x 0.3 m drift -> the recovered x channel commits, the scattered y is withheld
    // (same split as the 20b per-axis headline), so the rig is genuinely MID per-axis commit.
    RadarSource radar(3, tr, Xr, Vec3(Scalar(0.3), Scalar(0.002), Scalar(0.0)), /*seed=*/4242);

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[0].is_reference = true;
    sensors[1].prior_extrinsic = X1;
    sensors[2].prior_extrinsic = X2;
    sensors[3].prior_extrinsic = Xr;
    sensors[3].translation_only = true;
    sensors[3].scale_calib = false;

    Config cfg; set_hists(cfg);
    cfg.reference_sensor_id = 0;
    cfg.cold_start = ColdStart::MedianFromStart;
    cfg.commit_concentration = 0.5; cfg.commit_min_votes = 30; cfg.commit_drop = 0.2;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;
    REQUIRE(validate(cfg) == Status::Ok);

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    REQUIRE(est.add_source(&s0) == Status::Ok);
    REQUIRE(est.add_source(&s1) == Status::Ok);
    REQUIRE(est.add_source(&s2) == Status::Ok);
    REQUIRE(est.add_source(&radar) == Status::Ok);

    const Timestamp step_ns = static_cast<Timestamp>(std::llround(1e9 / 50.0));
    const Timestamp t_end = static_cast<Timestamp>(std::llround((tr.duration_s() - 0.1) * 1e9));
    for (Timestamp now = static_cast<Timestamp>(std::llround(0.2 * 1e9)); now < t_end; now += step_ns) {
        est.step(now);
    }

    const CalibSnapshot* cs = snap(est.latest(), 3);
    REQUIRE(cs != nullptr);
    // Confirm the pre-condition the round-trip must preserve: the clean x axis committed, the
    // noisy y NOT, the null z NOT (the rig is genuinely MID per-axis commit).
    REQUIRE(cs->translation_committed_xyz[0]);
    REQUIRE_FALSE(cs->translation_committed_xyz[1]);
    REQUIRE_FALSE(cs->translation_committed_xyz[2]);

    // Format must be v4 (the Slice-20b bump).
    REQUIRE(persist::kFormatVersion == 4u);

    unsigned char blob[8192];
    const Expected<int> wr = est.serialize(blob, sizeof(blob));
    REQUIRE(wr.ok());

    // PER-AXIS FOLD PIN (parse the blob's fusion prior_extrinsic directly): the persisted
    // prior_extrinsic.t must have the COMMITTED x axis FOLDED off the prior (the recovered
    // lever) while the withheld y + null z stay EXACTLY at the prior (Xr.t). This pins the
    // per-axis publish (the snapshot's extrinsic.t is the calibrator readout, not the fusion
    // prior -- so parse the blob, which carries the real prior). Schema: estimator.cpp serialize().
    Vec3 prior_t_persisted;
    {
        persist::Reader pr(blob + 20, wr.value() - 24);   // skip header; drop checksum
        const int n_rec = pr.get_i32();
        (void)pr.get_i32();          // phase
        (void)pr.get_bool();         // ever_fused
        bool found = false;
        for (int i = 0; i < n_rec; ++i) {
            const SourceId rid = static_cast<SourceId>(pr.get_i32());
            const SE3 Xp = pr.get_se3();
            (void)pr.get_f64();      // prior_scale
            (void)pr.get_f64();      // time_offset
            // 5 base + rot3d + scale2 + 3 per-axis lever (v4) = 10 commit flags.
            bool flags[10];
            for (int b = 0; b < 10; ++b) flags[b] = pr.get_bool();
            (void)pr.get_f64(); (void)pr.get_f64(); (void)pr.get_f64();   // EMA triple
            (void)pr.get_i32();      // resid_n
            if (rid == 3) {
                found = true;
                prior_t_persisted = Xp.t;
                // The persisted per-axis flags (positions 7,8,9) match the snapshot split.
                CHECK(flags[7]);            // x committed
                CHECK_FALSE(flags[8]);      // y withheld
                CHECK_FALSE(flags[9]);      // z withheld
            }
        }
        REQUIRE(found);
        REQUIRE_FALSE(pr.underflow);
    }
    // x folded OFF the prior (the recovered lever); y, z held EXACTLY at the prior.
    CHECK(std::abs(prior_t_persisted.x() - Xr.t.x()) > Scalar(0.02));   // x moved off prior
    CHECK(std::abs(prior_t_persisted.y() - Xr.t.y()) < Scalar(1e-12));  // y == prior
    CHECK(std::abs(prior_t_persisted.z() - Xr.t.z()) < Scalar(1e-12));  // z == prior

    // Restore into a fresh estimator with the SAME rig.
    Estimator fresh;
    REQUIRE(fresh.init(cfg) == Status::Ok);
    REQUIRE(fresh.add_source(&s0) == Status::Ok);
    REQUIRE(fresh.add_source(&s1) == Status::Ok);
    REQUIRE(fresh.add_source(&s2) == Status::Ok);
    REQUIRE(fresh.add_source(&radar) == Status::Ok);
    REQUIRE(fresh.deserialize(blob, wr.value()) == Status::Ok);

    // Drive ONE step so the restored snapshot republishes from the restored state (the commit
    // flags are held through the empty refilling histogram by the per-axis re-fill hysteresis).
    fresh.step(static_cast<Timestamp>(std::llround(0.2 * 1e9)));
    const CalibSnapshot* cf = snap(fresh.latest(), 3);
    REQUIRE(cf != nullptr);
    INFO("restored per-axis flags = {" << cf->translation_committed_xyz[0] << ","
         << cf->translation_committed_xyz[1] << "," << cf->translation_committed_xyz[2] << "}");
    // The per-axis flags round-tripped (held through the empty refilling histogram).
    CHECK(cf->translation_committed_xyz[0]);
    CHECK_FALSE(cf->translation_committed_xyz[1]);
    CHECK_FALSE(cf->translation_committed_xyz[2]);
    CHECK(cf->translation_committed);                 // whole = OR

    // RE-SERIALIZE round-trip: the restored estimator serializes the per-axis flags + folded
    // prior BYTE-IDENTICAL (before any new votes shift the histogram), proving the v4 flags +
    // the per-axis prior survive a full save -> restore -> save cycle.
    unsigned char blob2[8192];
    {
        Estimator fresh2;
        REQUIRE(fresh2.init(cfg) == Status::Ok);
        REQUIRE(fresh2.add_source(&s0) == Status::Ok);
        REQUIRE(fresh2.add_source(&s1) == Status::Ok);
        REQUIRE(fresh2.add_source(&s2) == Status::Ok);
        REQUIRE(fresh2.add_source(&radar) == Status::Ok);
        REQUIRE(fresh2.deserialize(blob, wr.value()) == Status::Ok);
        const Expected<int> wr2 = fresh2.serialize(blob2, sizeof(blob2));
        REQUIRE(wr2.ok());
        REQUIRE(wr2.value() == wr.value());
        bool same = true;
        for (int i = 0; i < wr.value(); ++i) if (blob[i] != blob2[i]) { same = false; break; }
        CHECK(same);
    }

    // A pre-v4 (stamped v3) blob cold-starts (rejected by the version guard -- old blobs never
    // restore the new per-axis flags; the established no-migration precedent).
    {
        unsigned char old_blob[8192];
        for (int i = 0; i < wr.value(); ++i) old_blob[i] = blob[i];
        old_blob[4] = 3u;                             // version word (LE) -> 3 (pre-20b)
        old_blob[5] = old_blob[6] = old_blob[7] = 0u;
        Estimator v3est;
        REQUIRE(v3est.init(cfg) == Status::Ok);
        REQUIRE(v3est.add_source(&s0) == Status::Ok);
        REQUIRE(v3est.add_source(&s1) == Status::Ok);
        REQUIRE(v3est.add_source(&s2) == Status::Ok);
        REQUIRE(v3est.add_source(&radar) == Status::Ok);
        CHECK(v3est.deserialize(old_blob, wr.value()) == Status::VersionMismatch);
    }
}
