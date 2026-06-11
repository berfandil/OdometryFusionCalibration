// Slice 17 tests: turn-regime FULL rotation extrinsic via axis-correspondence hand-eye
// (rot3d). Per turn-gated window the hand-eye rotation identity R_A R_X = R_X R_B is the
// axis correspondence a = R_X·b (a = log R_A, b = log R_B); Phase 2 accumulates the Wahba
// problem Mw = Σ w·a bᵀ + the axis Gram BBw = Σ w·b bᵀ and, once TWO distinct rotation
// axes have been excited (the BBw λ_mid gate — the observability spine), votes the
// basepoint-relative Kabsch residual δφ = log(R̂·R_bpᵀ) into 3 so(3) channels.
//
// Coverage (SLICE17 acceptance §3, items 1-7):
//   1. Clean conjugated multi-axis stream recovers a planted full R_X to < 1e-3 rad;
//      votes ≈ 0 at prior == truth; the a = R_X·b CONVENTION + the +yaw SIGN are pinned.
//   2. Observability self-test (LOAD-BEARING — never weaken): yaw-only -> BBw rank 1 ->
//      channels empty / conf 0 / never commits; multi-axis -> converges + commits.
//   3. Commit + contractive re-anchor: from a large wrong prior the basepoint walks to
//      truth (calibrator-level re-anchor + the estimator-level bootstrap mirror of the
//      Slice-8 ext tests).
//   4. Publish precedence: committed -> the full R drives prior/snapshot; uncommitted
//      (planar) -> the existing path, bit-identical to rot3d_enabled=false.
//   5. rot3d_enabled=false default: pinned byte-identical by (4)'s exact-equality run +
//      the untouched existing suite (the flag is default-off in Config).
//   6. Noise sanity: 2 mrad/step synthetic rotation noise -> recovery within a documented
//      loose bound; the rank gate stays honest (noise does not fake a second axis).
//   7. Persistence: committed flag/value round-trips (blob byte-equal after restore);
//      the config-hash covers the flag (a flip rejects); old format v1 blobs reject.
// (The loader key test lives in adapters/tests/test_config_loader.cpp.)
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
bool near_abs(Scalar a, Scalar b, Scalar tol) { return std::abs(a - b) <= tol; }

Mat3 Rz(Scalar a) { Mat3 R; R << std::cos(a), -std::sin(a), 0,
                                 std::sin(a),  std::cos(a), 0,
                                 0, 0, 1; return R; }
Mat3 Ry(Scalar a) { Mat3 R; R << std::cos(a), 0, std::sin(a),
                                 0, 1, 0,
                                -std::sin(a), 0, std::cos(a); return R; }
Mat3 Rx(Scalar a) { Mat3 R; R << 1, 0, 0,
                                 0, std::cos(a), -std::sin(a),
                                 0, std::sin(a),  std::cos(a); return R; }

// Sensor->base extrinsic = Rz(yaw) Ry(pitch) Rx(roll) — the config_loader
// parse_extrinsic / tools/inject_calib.py convention the sign pin is judged against.
SE3 make_extrinsic(Scalar yaw, Scalar pitch, Scalar roll, const Vec3& t) {
    SE3 X; X.R = Rz(yaw) * Ry(pitch) * Rx(roll); X.t = t;
    return X;
}

// Heap-allocate the calibrator (large fixed-capacity histograms overflow the test stack).
std::unique_ptr<Phase2Calibrator> make_calib() {
    return std::unique_ptr<Phase2Calibrator>(new Phase2Calibrator());
}

// Config with rot3d enabled + the Phase-2 histogram shapes the other calib tests use.
// vote_weight One so vote_count is a literal count. The rot3d channels reuse so3_hist
// (basepoint-relative residual, +-0.8 rad — holds even the large-wrong-prior votes).
Config make_cfg() {
    Config c;
    c.tick_rate_hz   = 50.0;
    c.reference_sensor_id = 0;
    c.turn_omega_min = 0.20;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.phase2_strategy = Phase2Strat::VsFusedBase;
    c.vote_weight     = VoteWeight::One;
    c.rot3d_enabled   = true;

    c.so3_hist.bins       = 512;
    c.so3_hist.range_min  = -0.8;
    c.so3_hist.range_max  =  0.8;
    c.so3_hist.circular   = false;
    c.so3_hist.aging      = Aging::SlidingK;
    c.so3_hist.sliding_k  = 256;
    c.so3_hist.vote_split = true;
    c.so3_hist.subbin     = true;

    c.scale_hist.bins = 512; c.scale_hist.range_min = 0.5; c.scale_hist.range_max = 1.5;
    c.scale_hist.aging = Aging::SlidingK; c.scale_hist.sliding_k = 256;

    c.roll_hist.bins = 360; c.roll_hist.range_min = -kPi; c.roll_hist.range_max = kPi;
    c.roll_hist.circular = true; c.roll_hist.aging = Aging::SlidingK; c.roll_hist.sliding_k = 256;

    c.xyz_hist.bins = 512; c.xyz_hist.range_min = -1.5; c.xyz_hist.range_max = 1.5;
    c.xyz_hist.aging = Aging::SlidingK; c.xyz_hist.sliding_k = 256;
    return c;
}

// One conjugated hand-eye window: base motion A (rotation `theta` about unit `axis` +
// a small translation), sensor report B = X^-1 o A o X — the exact measurement model the
// sim/inject_calib pin. Returns the observe() status.
Status feed_window(Phase2Calibrator& cal, SourceId id, const SE3& X,
                   const Vec3& axis, Scalar theta, Scalar rot_noise = 0,
                   std::mt19937_64* rng = nullptr) {
    SE3 A;
    A.R = so3::exp(axis.normalized() * theta);
    A.t = Vec3(0.5, 0.05, -0.02);                       // any sizable base translation
    SE3 B = se3::compose(se3::compose(se3::inverse(X), A), X);
    if (rot_noise > Scalar(0) && rng != nullptr) {
        // Independent white rotation noise on BOTH reports (the per-window composition of
        // per-step noise — see the noise-sanity test for the scaling).
        std::normal_distribution<Scalar> nd(0.0, rot_noise);
        const Vec3 na(nd(*rng), nd(*rng), nd(*rng));
        const Vec3 nb(nd(*rng), nd(*rng), nd(*rng));
        A.R = A.R * so3::exp(na);
        B.R = B.R * so3::exp(nb);
    }
    const Scalar dt = Scalar(0.5);                      // 0.5 s windows
    const Vec3 omega = so3::log(A.R) / dt;              // > turn gate for theta >= 0.15
    SourceId ids[1] = { id };
    SE3      rep[1] = { B };
    return cal.observe(1, ids, rep, A, omega);
}

// Feed `count` windows. multiaxis cycles through 3 genuinely distinct rotation axes
// (yaw-, pitch- and a mixed-axis turn); yaw-only keeps every window about +z (the planar
// ground regime — BBw must stay rank 1). theta varies window-to-window (excitation
// variety, all above the turn gate).
void feed_stream(Phase2Calibrator& cal, SourceId id, const SE3& X, int count,
                 bool multiaxis, Scalar rot_noise = 0, unsigned seed = 42) {
    std::mt19937_64 rng(seed);
    const Vec3 axes[3] = { Vec3(0, 0, 1), Vec3(0, 1, 0.2), Vec3(0.3, -0.4, 0.85) };
    for (int k = 0; k < count; ++k) {
        const Vec3 axis = multiaxis ? axes[k % 3] : Vec3(0, 0, 1);
        const Scalar theta = Scalar(0.25) + Scalar(0.10) * std::sin(Scalar(k) * Scalar(0.7));
        feed_window(cal, id, X, axis, theta, rot_noise, rot_noise > 0 ? &rng : nullptr);
    }
}

// --- estimator-level helpers (mirror test_calib_feedback) ---------------------------
// TURN-ONLY multi-axis trajectory — the EuRoC/drone regime the slice exists for: NO
// straight segments, so Phase-1 (straight-gated yaw/pitch) STARVES and the existing
// ext/roll path can never fix a yaw/pitch prior error. Only the rot3d publish can.
Trajectory multiaxis_turnonly_traj() {
    Trajectory tr;
    Vec6 turnA; turnA << 2.0, 0, 0, 0, 0.35,  0.6;   // yaw+pitch rate (axis 1)
    Vec6 turnB; turnB << 2.0, 0, 0, 0, -0.35, 0.6;   // flipped pitch (axis 2)
    for (int rep = 0; rep < 3; ++rep) {
        tr.add_segment(turnA, 1.6);
        tr.add_segment(turnB, 1.6);
    }
    return tr;
}

void set_hists(Config& c) {
    c.vote_weight = VoteWeight::One;
    c.straight_omega_max = 0.05;
    c.straight_trans_min = 0.02;
    c.turn_omega_min     = 0.20;
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

const CalibSnapshot* snap(const Result& r, SourceId id) {
    for (int i = 0; i < r.source_count; ++i)
        if (r.calib[i].id == id) return &r.calib[i];
    return nullptr;
}

// MIXED straight + multi-axis-turn trajectory: BOTH calibration regimes fire, so
// Phase-1 ext/roll AND rot3d can commit simultaneously (the precedence-stability rig).
Trajectory mixed_straight_turn_traj(int reps) {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0,     0;
    Vec6 turnA;    turnA    << 2.0, 0, 0, 0, 0.35,  0.6;
    Vec6 turnB;    turnB    << 2.0, 0, 0, 0, -0.35, 0.6;
    for (int rep = 0; rep < reps; ++rep) {
        tr.add_segment(straight, 1.5);
        tr.add_segment(turnA, 1.6);
        tr.add_segment(turnB, 1.6);
    }
    return tr;
}

// An ISource that switches between two underlying sources at t_switch (keyed on the
// query window END) — used to plant a mid-run MOUNT DRIFT so the rot3d commit can be
// driven through a genuine confidence DROP (full histogram, lost concentration).
class SwitchSource : public ISource {
public:
    SwitchSource(SourceId id, const ISource* before, const ISource* after,
                 Timestamp t_switch)
        : id_(id), before_(before), after_(after), t_switch_(t_switch) {}
    SourceId id() const override { return id_; }
    Expected<Delta> query(Timestamp t0, Timestamp t1) const override {
        return (t1 < t_switch_) ? before_->query(t0, t1) : after_->query(t0, t1);
    }
private:
    SourceId  id_;
    const ISource* before_;
    const ISource* after_;
    Timestamp t_switch_;
};
} // namespace

// ===========================================================================
// 1. Clean multi-axis recovery + the composition/sign convention pins
// ===========================================================================
TEST_CASE("rot3d: clean multi-axis stream recovers a planted full rotation < 1e-3 rad") {
    // EuRoC-magnitude planted mount: yaw 8 / pitch 5 / roll 4 deg.
    const SE3 X = make_extrinsic(8 * kPi / 180, 5 * kPi / 180, 4 * kPi / 180,
                                 Vec3(0.10, -0.05, 0.20));
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);     // IDENTITY prior (wrong-ish)

    feed_stream(cal, 1, X, 200, /*multiaxis=*/true);

    CHECK(cal.rot3d_observable(1));
    CHECK(cal.rot3d_vote_count(1) > Scalar(50));
    const Mat3 R = cal.rot3d(1);
    const Scalar err = so3::log(R * X.R.transpose()).norm();
    INFO("rot3d recovery err = " << err << " rad");
    CHECK(err < 1e-3);

    // CONVENTION PIN (HANDOFF §5): on clean conjugated data B = X^-1 o A o X the recovered
    // rotation satisfies a = R·b EXACTLY (a = log R_A, b = log R_B) — the sensor->base map.
    {
        SE3 A; A.R = so3::exp(Vec3(0.2, 0.5, 0.8).normalized() * Scalar(0.3));
        A.t = Vec3(0.4, 0, 0);
        const SE3 B = se3::compose(se3::compose(se3::inverse(X), A), X);
        const Vec3 a = so3::log(A.R);
        const Vec3 b = so3::log(B.R);
        CHECK((a - R * b).norm() < 1e-3);               // a = R_X · b (within mode resolution)
        CHECK((a - X.R * b).norm() < 1e-12);            // ...and exactly for the planted X
    }
}

TEST_CASE("rot3d: a planted +yaw reads +yaw (sign convention vs parse_extrinsic)") {
    // PURE +0.3 rad yaw mount. The recovered euler-zyx yaw must be +0.3 (not -0.3) —
    // pinning the composition direction against config_loader's "yaw pitch roll" parse
    // and tools/inject_calib.py.
    const Scalar yaw_t = 0.3;
    const SE3 X = make_extrinsic(yaw_t, 0, 0, Vec3::Zero());
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);

    feed_stream(cal, 1, X, 120, /*multiaxis=*/true);

    const Mat3 R = cal.rot3d(1);
    const Scalar yaw_rec = std::atan2(R(1, 0), R(0, 0));
    INFO("recovered yaw = " << yaw_rec);
    CHECK(near_abs(yaw_rec, yaw_t, 2e-3));
    CHECK(yaw_rec > 0.25);                              // genuinely POSITIVE (the sign pin)
}

TEST_CASE("rot3d: votes ~ 0 at prior == truth (fixed point)") {
    const SE3 X = make_extrinsic(0.14, -0.09, 0.07, Vec3(0.1, 0.0, -0.1));
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, X) == Status::Ok);         // prior == truth -> basepoint == X.R

    feed_stream(cal, 1, X, 150, /*multiaxis=*/true);

    // The residual mode sits at ~0 (votes deposit delta-phi = log(R_hat·R_bp^T) ~ 0), so
    // the readout equals the basepoint == truth, and the channels are CONCENTRATED.
    const Mat3 R = cal.rot3d(1);
    CHECK(so3::log(R * X.R.transpose()).norm() < 1e-3);
    CHECK(cal.rot3d_confidence(1) > Scalar(0.5));
}

// ===========================================================================
// 2. Observability self-test (LOAD-BEARING — never weaken)
// ===========================================================================
TEST_CASE("rot3d observability: yaw-only (rank-1) NEVER votes; multi-axis converges") {
    const SE3 X = make_extrinsic(0.20, -0.12, 0.30, Vec3(0.2, -0.1, 0.1));

    // Yaw-only stream: every window's rotation axis is +z -> BBw rank 1 -> the lambda_mid
    // gate must NEVER open: no votes, conf 0, readout frozen at the basepoint (the planar
    // ground regime stays on Phase-1 + roll — rot3d is inert there by construction).
    {
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        feed_stream(cal, 1, X, 200, /*multiaxis=*/false);

        CHECK_FALSE(cal.rot3d_observable(1));
        CHECK(cal.rot3d_vote_count(1) == doctest::Approx(0.0));
        CHECK(cal.rot3d_confidence(1) == doctest::Approx(0.0));
        // Readout frozen at the (identity-prior) basepoint — NOT pulled toward X.
        CHECK(so3::log(cal.rot3d(1)).norm() < 1e-12);
    }

    // Same planted mount, multi-axis stream: the gate opens and the solve converges.
    {
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        feed_stream(cal, 1, X, 200, /*multiaxis=*/true);

        CHECK(cal.rot3d_observable(1));
        CHECK(cal.rot3d_vote_count(1) > Scalar(50));
        CHECK(so3::log(cal.rot3d(1) * X.R.transpose()).norm() < 1e-3);
        CHECK(cal.rot3d_confidence(1) > Scalar(0.5));
    }
}

// ===========================================================================
// 3. Contractive re-anchor (calibrator level)
// ===========================================================================
TEST_CASE("rot3d re-anchor: large wrong basepoint -> recover -> re-anchor -> residual ~ 0") {
    const SE3 X = make_extrinsic(0.18, -0.10, 0.25, Vec3(0.3, 0.1, -0.2));
    // LARGE wrong prior basepoint: ~30 deg off in yaw + a pitch/roll error.
    const SE3 Xp = make_extrinsic(0.18 + 0.50, -0.10 + 0.15, 0.25 - 0.20, Vec3::Zero());
    REQUIRE(so3::log(Xp.R * X.R.transpose()).norm() > Scalar(0.4));

    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, Xp) == Status::Ok);

    // The Wahba solve is data-driven (basepoint-independent), so even from the wrong
    // basepoint the composed readout lands on truth (votes are just large residuals).
    feed_stream(cal, 1, X, 150, /*multiaxis=*/true);
    const Mat3 R1 = cal.rot3d(1);
    CHECK(so3::log(R1 * X.R.transpose()).norm() < 5e-3);

    // RE-ANCHOR (the estimator's rising-edge sequence): basepoint <- recovered, channels
    // reset. Mw/BBw are kept (basepoint-independent), so the very next window is already
    // observable and the residual re-votes ~ 0 around the new basepoint — contractive.
    REQUIRE(cal.set_rot3d_basepoint(1, R1) == Status::Ok);
    cal.reset_rot3d(1);
    CHECK(cal.rot3d_vote_count(1) == doctest::Approx(0.0));
    CHECK(so3::log(cal.rot3d(1) * R1.transpose()).norm() < 1e-12);   // fallback = new bp
    CHECK(cal.rot3d_observable(1));                                  // BBw survived

    feed_stream(cal, 1, X, 50, /*multiaxis=*/true, 0, /*seed=*/7);
    CHECK(cal.rot3d_vote_count(1) > Scalar(20));
    const Mat3 R2 = cal.rot3d(1);
    CHECK(so3::log(R2 * X.R.transpose()).norm() < 1e-3);             // still truth
    CHECK(so3::log(R2 * R1.transpose()).norm() < 5e-3);              // residual ~ 0
    CHECK(cal.rot3d_confidence(1) > Scalar(0.5));                    // re-concentrated
}

// ===========================================================================
// Lever coupling (review fix, MAJOR-1): once the two-axis gate is open, the
// TRANSLATION rows use the running Kabsch R̂ as R_X — so on turn-only motion
// (Phase-1 starved, R_yp frozen at the prior) the lever still recovers. This is
// the EuRoC regime and the spec §2 "lever coupling win": with the stale
// R_yp ∘ Rx(roll) row the ~|t|·θ bias (~0.33 rad × ~0.5 m here) would dwarf the
// planted lever and this pin fails — the R̂-driven row is load-bearing.
// ===========================================================================
TEST_CASE("rot3d lever coupling: turn-only planted rotation+lever — lever recovers "
          "despite Phase-1 never firing") {
    // Planted full mount: a large rotation (so the stale-row bias would be gross)
    // AND a lever; IDENTITY prior, so both start wrong and only rot3d can fix R.
    const SE3 X = make_extrinsic(0.25, -0.18, 0.12, Vec3(0.20, -0.15, 0.10));
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);

    // Turn-only conjugated stream (no straight windows ever reach Phase 2 anyway, and
    // set_yaw_pitch is never called — exactly the starved-Phase-1 bootstrap path).
    feed_stream(cal, 1, X, 250, /*multiaxis=*/true);

    CHECK(cal.rot3d_observable(1));
    CHECK(so3::log(cal.rot3d(1) * X.R.transpose()).norm() < 1e-3);   // rotation recovered
    const Vec3 lever = cal.lever_arm(1);
    INFO("lever err = " << (lever - X.t).norm() << " m  (planted |t| = " << X.t.norm() << ")");
    CHECK((lever - X.t).norm() < 0.02);                              // lever recovered too
    CHECK(cal.translation_confidence(1) > Scalar(0.5));
}

// ===========================================================================
// Combo vote weighting (review NIT): the λ_mid/λ_max gate is weight-ratio-
// invariant and the Wahba solve is mass-weighted — pin that rot3d behaves under
// the shipped real-data weighting (votes carry ‖ω‖ × Σ-confidence MASS, so
// rot3d_vote_count is a mass, mirroring the documented commit_min_votes
// semantics of the other DOFs).
// ===========================================================================
TEST_CASE("rot3d under Combo vote weighting: solve converges; yaw-only gate stays shut") {
    Config c = make_cfg();
    c.vote_weight = VoteWeight::Combo;
    const SE3 X = make_extrinsic(0.18, -0.11, 0.21, Vec3(0.10, 0.20, -0.10));
    const Vec3 axes[3] = { Vec3(0, 0, 1), Vec3(0, 1, 0.2), Vec3(0.3, -0.4, 0.85) };

    auto feed_combo = [&](Phase2Calibrator& cal, bool multiaxis) {
        for (int k = 0; k < 220; ++k) {
            const Vec3 axis = multiaxis ? axes[k % 3] : Vec3(0, 0, 1);
            const Scalar theta = Scalar(0.25) + Scalar(0.10) * std::sin(Scalar(k) * Scalar(0.7));
            SE3 A; A.R = so3::exp(axis.normalized() * theta); A.t = Vec3(0.5, 0.05, -0.02);
            const SE3 B = se3::compose(se3::compose(se3::inverse(X), A), X);
            const Scalar dt = Scalar(0.5);
            const Vec3 omega = so3::log(A.R) / dt;
            // Varying per-window source Σ-confidence — genuinely NON-UNIFORM vote mass.
            const Scalar conf = Scalar(0.4) + Scalar(0.6) *
                                (Scalar(0.5) + Scalar(0.5) * std::sin(Scalar(k) * Scalar(1.3)));
            SourceId ids[1] = { 1 };
            SE3      rep[1] = { B };
            Scalar   cf[1]  = { conf };
            cal.observe(1, ids, rep, A, omega, cf);
        }
    };

    // Multi-axis: converges under mass-weighted votes (loose bound — clean data).
    {
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(c, 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        feed_combo(cal, /*multiaxis=*/true);
        CHECK(cal.rot3d_observable(1));
        CHECK(cal.rot3d_vote_count(1) > Scalar(0));   // a MASS under Combo, not a count
        const Scalar err = so3::log(cal.rot3d(1) * X.R.transpose()).norm();
        INFO("Combo-weighted recovery err = " << err << " rad");
        CHECK(err < 5e-3);
        CHECK(cal.rot3d_confidence(1) > Scalar(0.5));
    }
    // Yaw-only: the gate is a RATIO test on BBw, so mass-weighting cannot fake a
    // second axis — channels stay empty (the observability self-test holds under Combo).
    {
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(c, 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        feed_combo(cal, /*multiaxis=*/false);
        CHECK_FALSE(cal.rot3d_observable(1));
        CHECK(cal.rot3d_vote_count(1) == doctest::Approx(0.0));
        CHECK(cal.rot3d_confidence(1) == doctest::Approx(0.0));
    }
}

// ===========================================================================
// 6. Noise sanity: 2 mrad/step synthetic noise; the rank gate stays honest
// ===========================================================================
TEST_CASE("rot3d noise sanity: 2 mrad/step recovery bounded; yaw-only gate stays shut") {
    // 0.5 s windows at 50 Hz = 25 steps; 2 mrad/step white rotation noise composes to
    // ~2e-3*sqrt(25) = 10 mrad per window, applied to BOTH A and B. The prototype's raw
    // batch Wahba read ~1.7 deg (0.03 rad) at 2-5 mrad/step on real EuRoC windows; the
    // DOCUMENTED LOOSE BOUND here is 0.06 rad (~3.4 deg) — histogram mode + decay add
    // robustness, but this test only pins the order of magnitude, not precision.
    const Scalar win_noise = Scalar(0.002) * std::sqrt(Scalar(25));
    const SE3 X = make_extrinsic(0.14, -0.087, 0.07, Vec3(0.1, 0.0, 0.2));

    {
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        feed_stream(cal, 1, X, 300, /*multiaxis=*/true, win_noise, /*seed=*/123);

        CHECK(cal.rot3d_observable(1));
        CHECK(cal.rot3d_vote_count(1) > Scalar(100));
        const Scalar err = so3::log(cal.rot3d(1) * X.R.transpose()).norm();
        INFO("noisy recovery err = " << err << " rad");
        CHECK(err < 0.06);
    }

    // GATE HONESTY: yaw-only + the same noise must NOT fake a second axis. The off-axis
    // noise mass in BBw is ~(0.01/0.3)^2 ~ 1e-3 of the about-axis mass — well below the
    // 1e-2 lambda_mid floor, so the channels must stay EMPTY.
    {
        auto calp = make_calib(); Phase2Calibrator& cal = *calp;
        REQUIRE(cal.configure(make_cfg(), 0) == Status::Ok);
        REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);
        feed_stream(cal, 1, X, 300, /*multiaxis=*/false, win_noise, /*seed=*/321);

        CHECK_FALSE(cal.rot3d_observable(1));
        CHECK(cal.rot3d_vote_count(1) == doctest::Approx(0.0));
        CHECK(cal.rot3d_confidence(1) == doctest::Approx(0.0));
    }
}

// ===========================================================================
// Decay aging: Mw/BBw are decay-aged per accepted window (so3_hist.decay_gamma), so the
// solve TRACKS a mount drift instead of averaging the old and new mounts forever.
// ===========================================================================
TEST_CASE("rot3d decay: the Wahba accumulator forgets an old mount (tracks drift)") {
    const SE3 X_old = make_extrinsic(0.30, -0.10, 0.05, Vec3::Zero());
    const SE3 X_new = make_extrinsic(0.05,  0.12, -0.08, Vec3::Zero());
    REQUIRE(so3::log(X_old.R * X_new.R.transpose()).norm() > Scalar(0.3));

    Config c = make_cfg();
    c.so3_hist.aging       = Aging::Decay;     // rot3d channels age too
    c.so3_hist.decay_gamma = 0.9;              // Mw/BBw decay per accepted window
    auto calp = make_calib(); Phase2Calibrator& cal = *calp;
    REQUIRE(cal.configure(c, 0) == Status::Ok);
    REQUIRE(cal.set_prior(1, SE3{}) == Status::Ok);

    // 150 windows at the OLD mount, then the mount DRIFTS: 150 windows at the NEW one.
    feed_stream(cal, 1, X_old, 150, /*multiaxis=*/true);
    feed_stream(cal, 1, X_new, 150, /*multiaxis=*/true, 0, /*seed=*/9);

    // With gamma=0.9 the old mass is down-weighted to ~1e-7 after 150 windows, so the
    // solve sits on the NEW mount. WITHOUT the Mw/BBw decay the accumulator would average
    // the two mounts (~0.17 rad off each) and this pin fails — the decay is load-bearing.
    const Scalar err_new = so3::log(cal.rot3d(1) * X_new.R.transpose()).norm();
    INFO("err vs new mount = " << err_new);
    CHECK(err_new < 0.02);
}

// ===========================================================================
// Determinism (calibrator level)
// ===========================================================================
TEST_CASE("rot3d determinism: identical input -> identical estimate") {
    const SE3 X = make_extrinsic(0.2, -0.1, 0.3, Vec3(0.4, -0.3, 0.15));
    auto run_once = [&]() {
        auto calp = make_calib();
        calp->configure(make_cfg(), 0);
        calp->set_prior(1, SE3{});
        feed_stream(*calp, 1, X, 150, true);
        const Vec3 r = so3::log(calp->rot3d(1));
        return std::vector<Scalar>{ r.x(), r.y(), r.z(),
                                    calp->rot3d_confidence(1),
                                    calp->rot3d_vote_count(1) };
    };
    const std::vector<Scalar> a = run_once();
    const std::vector<Scalar> b = run_once();
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(a[i] == b[i]);
}

// ===========================================================================
// 2b/3b. Estimator level: commit + bootstrap convergence on multi-axis motion,
// and (with persistence, item 7) the committed flag/value round-trip.
// ===========================================================================
TEST_CASE("rot3d estimator: commits on turn-only multi-axis motion, converges from a wrong "
          "prior that the existing path CANNOT fix, and round-trips through persistence") {
    const Trajectory tr = multiaxis_turnonly_traj();
    // Planted TRUE mount of source 3; the PRIOR rotation is ~15 deg off per axis (well
    // outside the small-deviation regime — the Slice-8 bootstrap shape). With NO straight
    // regime Phase-1 never votes, so the published rotation reaching truth REQUIRES the
    // rot3d publish precedence (the existing R_yp ∘ roll path can only correct roll).
    const SE3 X_true  = make_extrinsic(0.25, -0.18, 0.12, Vec3(0.20, -0.15, 0.10));
    const SE3 X_prior = make_extrinsic(0.25 + 0.26, -0.18 + 0.26, 0.12 - 0.20,
                                       Vec3(0.20, -0.15, 0.10));
    const Scalar prior_rerr = so3::log(X_true.R.transpose() * X_prior.R).norm();
    REQUIRE(prior_rerr > 0.35);

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X = X_true;                  // sources 0,1,2 identity (consensus anchor)
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = X_prior;

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
    cfg.rot3d_enabled = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    // COMMITTED + CONVERGED: the full 3-DOF rotation (incl. roll) walked from the large
    // wrong prior onto the planted mount, and the snapshot's extrinsic carries it (the
    // publish-precedence path: prior_extrinsic.R <- rot3d).
    CHECK(rig.estimator().rot3d_committed(3));
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    // Snapshot commit flag (review fix, MINOR): the published rotation is rot3d-driven
    // here while Phase-1 (straight-gated) never fired — extrinsic_committed must still
    // read TRUE (the rot3d commit is OR-ed in), so a consumer keying on the flag does
    // not misread the fully-committed published R as uncommitted.
    CHECK(cs->extrinsic_committed);
    const Scalar rerr = so3::log(X_true.R.transpose() * cs->extrinsic.R).norm();
    INFO("prior rot err=" << prior_rerr << "  converged rot err=" << rerr);
    CHECK(rerr < 0.10);
    CHECK(rerr < prior_rerr * Scalar(0.4));    // genuinely moved toward truth
    // The fused trajectory was not destabilized by the feedback.
    Scalar te, re; rig.max_error(te, re);
    CHECK(te < 0.5);
    CHECK(re < 0.5);

    // --- Persistence round-trip (acceptance item 7) ---------------------------------
    // Serialize the committed state; restore into a FRESH estimator with the same rig;
    // the rot3d commit flag must survive (held through the empty refilling channels by
    // the re-fill hysteresis) and the blob must re-serialize BYTE-IDENTICAL (the flag +
    // the committed value — the restored prior_extrinsic — round-trip losslessly).
    unsigned char blob[8192];
    const Expected<int> wr = rig.estimator().serialize(blob, sizeof(blob));
    REQUIRE(wr.ok());

    // PUBLISH PIN: the blob persists the FUSION prior_extrinsic, so source 3's persisted
    // rotation must be the rot3d-committed one (≈ X_true.R, NOT the wrong config prior) —
    // this directly pins the prior_extrinsic[i].R <- rot3d(id) publish, independent of the
    // snapshot path. Parse the documented payload schema (estimator.cpp serialize()).
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
            for (int b = 0; b < 6; ++b) (void)pr.get_bool();   // 5 + rot3d commit flags
            (void)pr.get_f64(); (void)pr.get_f64(); (void)pr.get_f64();   // EMA triple
            (void)pr.get_i32();      // resid_n
            if (rid == 3) {
                found = true;
                CHECK(so3::log(X_true.R.transpose() * Xp.R).norm() < 0.10);
                CHECK(so3::log(X_prior.R.transpose() * Xp.R).norm() > 0.25);  // moved off prior
            }
        }
        REQUIRE(found);
        REQUIRE_FALSE(pr.underflow);
    }

    Estimator fresh;
    REQUIRE(fresh.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(fresh.add_source(sp.get()) == Status::Ok);
    CHECK_FALSE(fresh.rot3d_committed(3));                    // cold before restore
    REQUIRE(fresh.deserialize(blob, wr.value()) == Status::Ok);
    CHECK(fresh.rot3d_committed(3));                          // restored + held

    unsigned char blob2[8192];
    const Expected<int> wr2 = fresh.serialize(blob2, sizeof(blob2));
    REQUIRE(wr2.ok());
    REQUIRE(wr2.value() == wr.value());
    bool same = true;
    for (int i = 0; i < wr.value(); ++i) {
        if (blob[i] != blob2[i]) { same = false; break; }
    }
    CHECK(same);
}

// ===========================================================================
// Cold start under ReferenceOnly (review fix, MAJOR-2): on a never-straight rig
// (turn-only — the EuRoC/drone regime) Phase-1 starves, so ext_committed can never
// latch; a rot3d-committed source must STILL join the fusion median (participates()
// keys on EITHER rotation commit). Without the fix the non-reference source dead-
// reckons behind the reference forever and the slice's win never reaches fusion.
// ===========================================================================
TEST_CASE("rot3d cold start: ReferenceOnly + turn-only rig — source joins the median "
          "after the rot3d commit (and never joins with the knob off)") {
    const Trajectory tr = multiaxis_turnonly_traj();
    const SE3 X1 = make_extrinsic(0.20, -0.12, 0.15, Vec3(0.20, -0.10, 0.05));

    for (int pass = 0; pass < 2; ++pass) {
        const bool rot3d_on = (pass == 1);

        std::vector<SourceParams> planted(2);
        planted[0].id = 0;
        planted[1].id = 1;
        planted[1].X = X1;
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

        std::vector<SensorConfig> sensors(2);
        sensors[0].id = 0;
        sensors[1].id = 1;
        sensors[1].prior_extrinsic = X1;   // prior == truth: the GAP is participates(), not
                                           // convergence — rot3d votes ≈ 0 and commits fast.
        Config cfg;
        cfg.max_sources = 2; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
        cfg.timesync_enabled = false; cfg.cold_start = ColdStart::ReferenceOnly;
        set_hists(cfg);
        cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
        cfg.rot3d_enabled = rot3d_on;
        cfg.sensors = sensors.data(); cfg.sensor_count = 2;

        Rig rig;
        rig.set_trajectory(tr);
        REQUIRE(rig.init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
        const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 200);

        const Result& r = rig.estimator().latest();
        REQUIRE(r.source_count == 2);
        REQUIRE(r.health[1].id == 1);
        if (rot3d_on) {
            // Committed (turn-only, prior == truth) -> PARTICIPATES: a nonzero fusion
            // weight on the final fused step proves the source joined the median.
            CHECK(rig.estimator().rot3d_committed(1));
            CHECK(r.health[1].weight > Scalar(0));
            // ...and the join is genuinely gated on the commit: before any commit the
            // recorded weights are 0 (scan the early fused records).
            bool early_weight_zero = true;
            int  seen = 0;
            for (const auto& rec : rig.records()) {
                if (!rec.fused) continue;
                if (++seen > 30) break;                 // the first ~30 fused steps are
                if (rec.result.health[1].weight > Scalar(0)) early_weight_zero = false;
            }                                           // well before commit_min_votes=60
            CHECK(early_weight_zero);
        } else {
            // Knob off on the same never-straight rig: ext never commits -> the source
            // NEVER joins (the pre-fix behaviour, now pinned as the off-path).
            CHECK_FALSE(rig.estimator().rot3d_committed(1));
            CHECK(r.health[1].weight == Scalar(0));
        }
    }
}

// ===========================================================================
// Estimator rising-edge re-anchor mutation pin (review fix, MINOR): a prior rotation
// error BEYOND the so3_hist range (1.2 rad yaw vs the ±0.8 rad channels). The δφ votes
// CLAMP into the edge bin, so the first commit publishes a clamped (≈0.8 rad) partial
// correction — only the rising-edge re-anchor (basepoint <- committed value + channel
// reset) lets the remaining ≈0.4 rad residual re-vote IN range and converge. Deleting
// the set_rot3d_basepoint/reset_rot3d pair in apply_calib_feedback leaves the published
// rotation stuck ≈0.4 rad off truth and this pin fails (mirrors the Slice-8 large-error
// ext walk).
// ===========================================================================
TEST_CASE("rot3d estimator bootstrap: prior error beyond the histogram range converges "
          "only through the rising-edge re-anchor") {
    const Trajectory tr = multiaxis_turnonly_traj();
    const SE3 X_true  = make_extrinsic(0.15, -0.10, 0.12, Vec3(0.20, -0.15, 0.10));
    const SE3 X_prior = make_extrinsic(0.15 + 1.20, -0.10, 0.12, X_true.t);   // yaw +1.2 rad
    const Scalar prior_rerr = so3::log(X_true.R.transpose() * X_prior.R).norm();
    REQUIRE(prior_rerr > 1.1);                       // genuinely beyond the ±0.8 range

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X = X_true;
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = X_prior;

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
    cfg.rot3d_enabled = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    CHECK(rig.estimator().rot3d_committed(3));
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    const Scalar rerr = so3::log(X_true.R.transpose() * cs->extrinsic.R).norm();
    INFO("prior rot err=" << prior_rerr << "  converged rot err=" << rerr);
    // The clamped single-commit value sits ≈0.4 rad off truth; converging well below
    // that REQUIRES the contractive re-anchor walk.
    CHECK(rerr < 0.15);
    CHECK(rerr < prior_rerr * Scalar(0.15));
}

// ===========================================================================
// Pre-gate lever pollution (real-data follow-up, EuRoC euroc_calib_rot3d.ini): lever
// rows deposited BEFORE the BBw two-axis gate opens are built with the WRONG row
// rotation R_X (on turn-only motion the R_yp ∘ Rx(roll) composition is frozen at the
// prior — Phase-1 starves), and the lever normal equations ata_/atb_/rows_ have NO
// aging, so that early pollution sits in the cumulative LS forever; the xyz channels
// vote the RUNNING ridge solve, so the polluted cumulative solution drags the mode
// even after gate-open rows turn clean (on EuRoC: rotation recovered to 1.6e-7 rad and
// lever y to 1.4 mm, but lever x stuck 3.5 cm off — truth 0.30 m, read 0.2653 m).
// Rig: a LONG yaw-only (rank-1 — gate shut) turning phase deposits hundreds of
// polluted rows, then multi-axis turns open the gate and rot3d converges + commits.
// The fix under test: the rot3d COMMIT RISING EDGE in apply_calib_feedback() also
// reset_lever()s — drop the polluted accumulator + xyz histograms (rows/votes cast at
// the wrong rotation are stale, exactly the so3-reset-on-re-anchor rationale); the
// rows rebuild quickly from gate-open R̂-driven windows and the lever lands on truth.
// MUTATION GUARD: deleting the reset_lever call in apply_calib_feedback leaves the
// published lever dragged by the polluted cumulative solve and the bound below fails.
// ===========================================================================
TEST_CASE("rot3d lever pollution: pre-gate yaw-only rows are dropped on the commit "
          "rising edge — lever converges despite early wrong-R_X pollution") {
    // 8 s of PURE yaw turning (turn-gated -> lever rows deposit; rank-1 -> rot3d gate
    // shut -> every row uses the wrong prior rotation), then the multi-axis turns.
    Trajectory tr;
    Vec6 yawturn; yawturn << 2.0, 0, 0, 0, 0,     0.6;
    Vec6 turnA;   turnA   << 2.0, 0, 0, 0, 0.35,  0.6;
    Vec6 turnB;   turnB   << 2.0, 0, 0, 0, -0.35, 0.6;
    tr.add_segment(yawturn, 8.0);
    for (int rep = 0; rep < 3; ++rep) {
        tr.add_segment(turnA, 1.6);
        tr.add_segment(turnB, 1.6);
    }

    // Planted mount: rotation ≥ 0.3 rad + a lever; IDENTITY prior (both start wrong,
    // and with no straight regime only rot3d can ever fix the rotation).
    const SE3 X_true = make_extrinsic(0.25, -0.18, 0.12, Vec3(0.20, -0.15, 0.10));
    REQUIRE(so3::log(X_true.R).norm() > Scalar(0.30));

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X = X_true;                  // sources 0,1,2 identity (consensus anchor)
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = SE3{};     // identity prior — rotation AND lever wrong

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
    cfg.rot3d_enabled = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 500);

    // The rotation committed + converged (the rising edge genuinely fired).
    REQUIRE(rig.estimator().rot3d_committed(3));
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 3);
    REQUIRE(cs != nullptr);
    const Scalar rerr = so3::log(X_true.R.transpose() * cs->extrinsic.R).norm();
    INFO("converged rot err = " << rerr << " rad");
    CHECK(rerr < 0.05);

    // THE PIN: the published lever lands on truth INCLUDING the early-window pollution
    // (~400 polluted rows). Without the rising-edge reset_lever the cumulative solve
    // stays dragged several cm (the EuRoC lever-x failure) and this bound fails.
    const Scalar lerr = (cs->extrinsic.t - X_true.t).norm();
    INFO("lever err = " << lerr << " m  (planted |t| = " << X_true.t.norm() << ")");
    CHECK(lerr < 0.02);
}

// ===========================================================================
// Precedence stability (review fix, MINOR): MIXED straight + multi-axis-turn run where
// Phase-1 ext/roll AND rot3d BOTH commit. Pins: (a) ext commits first (straight regime)
// and rot3d on top — both latch; (b) while rot3d is committed the published rotation is
// STABLE (no ext/roll-vs-rot3d publish thrash) and converged; (c) a genuine rot3d DROP
// (a planted mid-run mount drift deconcentrates the full channels) falls back CLEANLY —
// the published rotation stays orthonormal, fusion keeps running, and the fused error
// stays bounded throughout.
// ===========================================================================
TEST_CASE("rot3d precedence: mixed straight+turn — both paths commit without thrash; "
          "clean fallback when rot3d drops") {
    const int  reps = 7;
    const Trajectory tr = mixed_straight_turn_traj(reps);
    const Scalar rep_s = 4.7;
    const Scalar t_switch_s = 3 * rep_s;             // mount drifts after 3 clean reps
    const SE3 X_true  = make_extrinsic(0.10, -0.07, 0.09, Vec3(0.15, -0.10, 0.05));
    const SE3 X_drift = make_extrinsic(0.10 + 0.35, -0.07, 0.09 - 0.25, X_true.t);
    const SE3 X_prior = make_extrinsic(0.10 + 0.08, -0.07 - 0.06, 0.09, X_true.t);

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[3].X = X_true;
    SourceParams drifted = planted[3];
    drifted.X = X_drift;

    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (int i = 0; i < 3; ++i) srcs.emplace_back(new SyntheticSource(tr, planted[i]));
    SyntheticSource src3_before(tr, planted[3]);
    SyntheticSource src3_after(tr, drifted);
    SwitchSource sw3(3, &src3_before, &src3_after,
                     static_cast<Timestamp>(t_switch_s * 1e9));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[3].prior_extrinsic = X_prior;

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    // Decay aging so BOTH the channels and Mw/BBw track the drift (a SlidingK histogram
    // only evicts on new votes; decay lets concentration genuinely fall during the walk).
    // gamma=0.98 saturates the total vote MASS at ~50, so commit_min_votes=30.
    cfg.so3_hist.aging = Aging::Decay;   cfg.so3_hist.decay_gamma = 0.98;
    cfg.roll_hist.aging = Aging::Decay;  cfg.roll_hist.decay_gamma = 0.98;
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 30;
    cfg.rot3d_enabled = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Estimator est;
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(est.add_source(sp.get()) == Status::Ok);
    REQUIRE(est.add_source(&sw3) == Status::Ok);

    const Scalar tick_hz = 50.0;
    const Timestamp tick = static_cast<Timestamp>(1e9 / tick_hz);
    const Timestamp t0   = static_cast<Timestamp>(0.2 * 1e9);
    const Timestamp tend = static_cast<Timestamp>((tr.duration_s() - 0.1) * 1e9);
    const Timestamp t_switch = static_cast<Timestamp>(t_switch_s * 1e9);
    const Timestamp window_ns = static_cast<Timestamp>(cfg.window_s * 1e9);

    bool saw_ext_only      = false;   // Phase-1 ext committed BEFORE any rot3d commit
    bool saw_rot3d         = false;
    bool saw_drop          = false;   // rot3d went committed -> uncommitted post-drift
    int  clean_rises       = 0;       // rot3d rising edges during the CLEAN phase
    int  clean_drops       = 0;       // rot3d drops during the CLEAN phase (must be 0)
    int  jumps_committed   = 0;       // published-R steps > 0.05 rad while committed (clean)
    bool prev_rc           = false;
    bool have_prev_R       = false;
    Mat3 prev_R = Mat3::Identity();
    Scalar clean_end_rerr  = Scalar(-1);
    Scalar max_terr = 0, max_rerr = 0;
    bool have_anchor = false;
    SE3  anchor_inv;
    int  fuses = 0;

    for (Timestamp now = t0; now <= tend; now += tick) {
        const Status st = est.step(now);
        const bool rc = est.rot3d_committed(3);
        const bool clean_phase = (now < t_switch);
        if (rc && !prev_rc && clean_phase) ++clean_rises;
        if (!rc && prev_rc) {
            if (clean_phase) ++clean_drops; else saw_drop = true;
        }
        prev_rc = rc;
        if (!ok(st)) { have_prev_R = false; continue; }
        ++fuses;

        const Result& r = est.latest();
        const CalibSnapshot* cs = snap(r, 3);
        REQUIRE(cs != nullptr);

        // (a) Phase-1's OWN commit fires first (straight regime) — observable as the
        // snapshot flag set while rot3d is still uncommitted.
        if (cs->extrinsic_committed && !rc && !saw_rot3d) saw_ext_only = true;
        if (rc) saw_rot3d = true;

        // (c) The published rotation is ALWAYS sane (orthonormal, finite) — committed,
        // fallen-back, or mid-walk.
        const Mat3& R = cs->extrinsic.R;
        REQUIRE(R.allFinite());
        REQUIRE((R.transpose() * R - Mat3::Identity()).norm() < 1e-9);

        // (b) No publish thrash while committed on CLEAN data: consecutive published
        // rotations move smoothly (the rising-edge step itself is excluded — committed
        // on both ends of the pair).
        if (clean_phase && rc && have_prev_R) {
            if (so3::log(R * prev_R.transpose()).norm() > Scalar(0.05)) ++jumps_committed;
        }
        prev_R = R; have_prev_R = rc;
        if (clean_phase && rc) {
            clean_end_rerr = so3::log(X_true.R.transpose() * R).norm();
        }

        // Fused-vs-GT error (the rig's anchoring, inlined — manual loop so the per-step
        // commit flags are pollable).
        const Timestamp frontier = r.frontier.stamp;
        if (!have_anchor) {
            anchor_inv = se3::inverse(tr.pose(frontier - window_ns));
            have_anchor = true;
        }
        const SE3 gt = se3::compose(anchor_inv, tr.pose(frontier));
        const SE3 e  = se3::compose(se3::inverse(gt), r.frontier.pose);
        max_terr = std::max(max_terr, e.t.norm());
        max_rerr = std::max(max_rerr, so3::log(e.R).norm());
    }

    REQUIRE(fuses > 500);
    // (a) Both paths committed — ext first (straight), rot3d on top (multi-axis turns).
    CHECK(saw_ext_only);
    CHECK(saw_rot3d);
    // (b) No thrash on clean data: EXACTLY ONE rising edge, ZERO drops, smooth publish.
    CHECK(clean_rises == 1);
    CHECK(clean_drops == 0);
    CHECK(jumps_committed == 0);
    // ...and the committed published rotation converged onto the planted mount.
    REQUIRE(clean_end_rerr >= Scalar(0));
    CHECK(clean_end_rerr < 0.05);
    // (c) The planted drift deconcentrated the FULL channels -> a genuine drop, and the
    // run survived it (fusion never destabilized; the median carries 3 clean sources).
    CHECK(saw_drop);
    INFO("max fused terr=" << max_terr << " rerr=" << max_rerr);
    CHECK(max_terr < 1.0);
    CHECK(max_rerr < 0.6);
}

// ===========================================================================
// Estimator commit under Combo weighting (review NIT): rot3d_vote_count is a vote MASS
// under Combo (‖ω‖ × Σ-confidence per vote), so commit_min_votes is a mass threshold —
// pin that the commit machinery still latches on a turn-only rig with the shipped
// real-data weighting (loose: just the latch + a sane published rotation).
// ===========================================================================
TEST_CASE("rot3d estimator commit under Combo vote weighting") {
    const Trajectory tr = multiaxis_turnonly_traj();
    const SE3 X1 = make_extrinsic(0.20, -0.12, 0.15, Vec3(0.20, -0.10, 0.05));

    std::vector<SourceParams> planted(4);
    for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
    planted[1].X = X1;
    std::vector<std::unique_ptr<SyntheticSource>> srcs;
    for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));

    std::vector<SensorConfig> sensors(4);
    for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
    sensors[1].prior_extrinsic = X1;                 // prior == truth (commit, not walk)

    Config cfg;
    cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
    cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
    set_hists(cfg);
    cfg.vote_weight = VoteWeight::Combo;             // mass-weighted votes
    cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 60;
    cfg.rot3d_enabled = true;
    cfg.sensors = sensors.data(); cfg.sensor_count = 4;

    Rig rig;
    rig.set_trajectory(tr);
    REQUIRE(rig.init(cfg) == Status::Ok);
    for (auto& sp : srcs) REQUIRE(rig.add_source(sp.get()) == Status::Ok);
    const int fuses = rig.run(0.2, tr.duration_s() - 0.1, 50.0);
    REQUIRE(fuses > 200);

    CHECK(rig.estimator().rot3d_committed(1));
    const CalibSnapshot* cs = snap(rig.estimator().latest(), 1);
    REQUIRE(cs != nullptr);
    CHECK(cs->extrinsic_committed);
    CHECK(so3::log(X1.R.transpose() * cs->extrinsic.R).norm() < 0.05);
}

// ===========================================================================
// 4/5. Publish precedence + default-off: planar (yaw-only) motion with the flag ON is
// BIT-IDENTICAL to the flag OFF (rot3d never votes/commits there -> the existing path
// publishes as today). This is also the BBw-gate mutation guard: removing the gate lets
// planar windows vote/commit and the exact-equality pins below fail.
// ===========================================================================
TEST_CASE("rot3d planar inertness: yaw-only run with rot3d_enabled is bit-identical to off") {
    Trajectory tr;
    Vec6 straight; straight << 2.0, 0, 0, 0, 0, 0;
    Vec6 yawturn;  yawturn  << 2.0, 0, 0, 0, 0, 0.6;       // PURE yaw rate (planar)
    for (int rep = 0; rep < 2; ++rep) {
        tr.add_segment(straight, 1.5);
        tr.add_segment(yawturn,  2.0);
    }

    const SE3 X1 = make_extrinsic(0.15, -0.10, 0.20, Vec3(0.25, -0.15, 0.10));

    auto run = [&](bool rot3d_on) {
        std::vector<SourceParams> planted(4);
        for (int i = 0; i < 4; ++i) planted[i].id = static_cast<SourceId>(i);
        planted[1].X = X1;
        std::vector<std::unique_ptr<SyntheticSource>> srcs;
        for (const auto& sp : planted) srcs.emplace_back(new SyntheticSource(tr, sp));
        std::vector<SensorConfig> sensors(4);
        for (int i = 0; i < 4; ++i) sensors[i].id = static_cast<SourceId>(i);
        sensors[1].prior_extrinsic = X1;
        Config cfg;
        cfg.max_sources = 4; cfg.fusion_delay_s = 0.05; cfg.window_s = 0.10;
        cfg.timesync_enabled = false; cfg.cold_start = ColdStart::MedianFromStart;
        set_hists(cfg);
        cfg.commit_concentration = 0.5; cfg.commit_drop = 0.3; cfg.commit_min_votes = 40;
        cfg.rot3d_enabled = rot3d_on;
        cfg.sensors = sensors.data(); cfg.sensor_count = 4;
        auto rig = std::unique_ptr<Rig>(new Rig());
        rig->set_trajectory(tr);
        REQUIRE(rig->init(cfg) == Status::Ok);
        for (auto& sp : srcs) REQUIRE(rig->add_source(sp.get()) == Status::Ok);
        const int fuses = rig->run(0.2, tr.duration_s() - 0.1, 50.0);
        REQUIRE(fuses > 100);
        // Yaw-only: rot3d must never commit, regardless of the flag.
        CHECK_FALSE(rig->estimator().rot3d_committed(1));
        return rig;
    };

    auto rig_off = run(false);
    auto rig_on  = run(true);

    // EXACT (bit-identical) equality of the published outputs: the frontier pose and the
    // per-source calibration snapshot. Any behavioral leak of the enabled-but-unobservable
    // rot3d path (a planar vote, a commit, a publish) breaks these exact pins.
    const Result& a = rig_off->estimator().latest();
    const Result& b = rig_on->estimator().latest();
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) CHECK(a.frontier.pose.R(r, c) == b.frontier.pose.R(r, c));
        CHECK(a.frontier.pose.t(r) == b.frontier.pose.t(r));
    }
    REQUIRE(a.source_count == b.source_count);
    for (int i = 0; i < a.source_count; ++i) {
        const CalibSnapshot& ca = a.calib[i];
        const CalibSnapshot& cb = b.calib[i];
        CHECK(ca.id == cb.id);
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) CHECK(ca.extrinsic.R(r, c) == cb.extrinsic.R(r, c));
        CHECK(ca.extrinsic.t.x() == cb.extrinsic.t.x());
        CHECK(ca.extrinsic.t.y() == cb.extrinsic.t.y());
        CHECK(ca.extrinsic.t.z() == cb.extrinsic.t.z());
        CHECK(ca.scale == cb.scale);
        CHECK(ca.extrinsic_committed == cb.extrinsic_committed);
        CHECK(ca.scale_committed == cb.scale_committed);
        CHECK(ca.translation_committed == cb.translation_committed);
    }
}

// ===========================================================================
// 7. Config-hash covers the flag; old format (v1) blobs reject
// ===========================================================================
TEST_CASE("rot3d persistence guards: config-hash flip rejects; format v1 rejects") {
    // A minimal valid rig config (no sources needed — the framing checks run first).
    std::vector<SensorConfig> sensors(2);
    sensors[0].id = 0; sensors[1].id = 1;
    Config cfg;
    cfg.max_sources = 2;
    set_hists(cfg);
    cfg.sensors = sensors.data(); cfg.sensor_count = 2;
    cfg.rot3d_enabled = false;

    auto est = std::unique_ptr<Estimator>(new Estimator());
    REQUIRE(est->init(cfg) == Status::Ok);
    unsigned char blob[4096];
    const Expected<int> wr = est->serialize(blob, sizeof(blob));
    REQUIRE(wr.ok());

    // (a) CONFIG-HASH FLIP: the same rig with rot3d_enabled=true must REJECT the blob
    // written with the flag off (InvalidConfig) — the flag is calibration-shaping.
    {
        Config cfg_on = cfg;
        cfg_on.rot3d_enabled = true;
        auto est_on = std::unique_ptr<Estimator>(new Estimator());
        REQUIRE(est_on->init(cfg_on) == Status::Ok);
        CHECK(est_on->deserialize(blob, wr.value()) == Status::InvalidConfig);
    }

    // (b) FORMAT v1 REJECT: a blob stamped with the OLD format version (1) must reject
    // with VersionMismatch (version is checked BEFORE the checksum — the documented
    // order), exactly the "old blobs reject -> cold start" contract of the v2 bump.
    {
        REQUIRE(persist::kFormatVersion == 2u);          // this slice's bump (pin it)
        unsigned char old_blob[4096];
        for (int i = 0; i < wr.value(); ++i) old_blob[i] = blob[i];
        old_blob[4] = 1u;                                 // version word (LE) -> 1
        old_blob[5] = old_blob[6] = old_blob[7] = 0u;
        CHECK(est->deserialize(old_blob, wr.value()) == Status::VersionMismatch);
    }
}
