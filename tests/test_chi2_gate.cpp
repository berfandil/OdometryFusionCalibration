// test_chi2_gate.cpp — Slice 11b: per-DOF Mahalanobis chi-square gate (Eskf::chi2_gate).
//
// The single-scalar Mahalanobis gate (one cfg.mahalanobis_chi2 applied regardless of the
// measurement DOF n) gated a 3-DOF position fix and a 6-DOF pose fix at DIFFERENT statistical
// confidence (the base is tuned at the n=3 ~95% quantile, so n=6 only saw ~80%). chi2_gate()
// scales that single base by the chi-square quantile ratio q[n]/q[3] so EVERY n in 1..6 gates
// at the SAME confidence level. This file proves:
//   (1) n==3 returns the base unchanged (back-compat: default 9.0 / loosened 100.0 stay identical
//       to the dim=3 fixes shipping now, since the ratio at n=3 is exactly 1.0).
//   (2) the threshold is strictly increasing in n (a higher-DOF fix is allowed a larger NIS at
//       the same confidence).
//   (3) n is clamped to [1,6] (0 -> n=1 value, >6 -> n=6 value).
//   (4) base <= 0 returns base (validate() forbids it; the helper must not produce NaN/garbage).
//   (5) a behavioral end-to-end check through Eskf::update: a dim=6 measurement whose NIS sits
//       BETWEEN the n=3 base threshold and the per-n(6) threshold is REJECTED under the flat base
//       but ACCEPTED under chi2_gate(base, 6) — the per-n gate is correctly less conservative for
//       higher DOF at fixed confidence.
//   (6) a numeric spot-check of one ratio value.
#include <doctest/doctest.h>

#include "ofc/core/absolute_ref.hpp"
#include "ofc/core/eskf.hpp"
#include "ofc/core/lie.hpp"

#include <cmath>

using namespace ofc;

namespace {
// The 0.95-confidence chi-square quantiles the gate is built on, DOF 1..6 (index 0 unused).
// Mirrors the const table in eskf.cpp so the tests pin the exact expected scaling.
constexpr Scalar kQ[7] = {Scalar(0),
                          Scalar(3.841459), Scalar(5.991465), Scalar(7.814728),
                          Scalar(9.487729), Scalar(11.070498), Scalar(12.591587)};
} // namespace

// (1) n==3 is the identity (back-compat anchor).
TEST_CASE("chi2_gate: n==3 returns the base unchanged (back-compat)") {
    CHECK(Eskf::chi2_gate(9.0, 3) == 9.0);       // the default gate
    CHECK(Eskf::chi2_gate(100.0, 3) == 100.0);   // the loosened drift-removal test value
    CHECK(Eskf::chi2_gate(1.0, 3) == 1.0);
}

// (2) strictly increasing in n for a fixed base > 0.
TEST_CASE("chi2_gate: threshold is strictly increasing in n (1..6) at fixed base") {
    const Scalar base = 9.0;
    Scalar prev = Eskf::chi2_gate(base, 1);
    for (int n = 2; n <= 6; ++n) {
        const Scalar cur = Eskf::chi2_gate(base, n);
        CHECK(cur > prev);
        prev = cur;
    }
    // Each value matches base * q[n]/q[3] (the documented scaling).
    for (int n = 1; n <= 6; ++n) {
        CHECK(Eskf::chi2_gate(base, n) == doctest::Approx(base * kQ[n] / kQ[3]).epsilon(1e-12));
    }
}

// (3) n is clamped to [1,6].
TEST_CASE("chi2_gate: n is clamped to [1,6] (0 -> n=1, >6 -> n=6)") {
    const Scalar base = 9.0;
    CHECK(Eskf::chi2_gate(base, 0)   == Eskf::chi2_gate(base, 1));   // below range -> n=1
    CHECK(Eskf::chi2_gate(base, -5)  == Eskf::chi2_gate(base, 1));
    CHECK(Eskf::chi2_gate(base, 7)   == Eskf::chi2_gate(base, 6));   // above range -> n=6
    CHECK(Eskf::chi2_gate(base, 100) == Eskf::chi2_gate(base, 6));
}

// (4) base <= 0 passthrough (the validate()-forbidden disabled path; never NaN/garbage).
TEST_CASE("chi2_gate: base <= 0 returns base unchanged") {
    CHECK(Eskf::chi2_gate(0.0, 6)  == 0.0);
    CHECK(Eskf::chi2_gate(-1.0, 6) == -1.0);
    CHECK(Eskf::chi2_gate(0.0, 3)  == 0.0);
    CHECK(Eskf::chi2_gate(-2.5, 1) == -2.5);
}

// (6) numeric spot-check: scaling the n=3 quantile to n=6 lands on the n=6 quantile.
TEST_CASE("chi2_gate: numeric ratio spot-check (q3 scaled to n=6 == q6)") {
    CHECK(Eskf::chi2_gate(7.814728, 6) == doctest::Approx(12.591587).epsilon(1e-6));
    CHECK(Eskf::chi2_gate(7.814728, 1) == doctest::Approx(3.841459).epsilon(1e-6));
}

// (5) behavioral: a 6-DOF fix whose NIS is between the flat base and the per-n(6) threshold is
// REJECTED at the flat base but ACCEPTED at chi2_gate(base, 6). We construct P, H, R and the
// residual directly so we control S and thus the NIS exactly.
TEST_CASE("chi2_gate: a dim=6 fix between base and per-n(6) is rejected flat, accepted per-n") {
    const Scalar base = 9.0;
    const Scalar gate6 = Eskf::chi2_gate(base, 6);   // ~ 9.0 * 12.591587/7.814728 ~= 14.50
    REQUIRE(gate6 > base);

    // Build a full-rank dim=6 pose measurement with S = I6 (H = I, P = 0, R = I6) so the NIS is
    // exactly the squared residual norm. Pick a residual whose NIS lands strictly BETWEEN base
    // and gate6 -> rejected under base, accepted under gate6.
    const Scalar target_nis = Scalar(0.5) * (base + gate6);   // midpoint, strictly in (base, gate6)
    REQUIRE(target_nis > base);
    REQUIRE(target_nis < gate6);

    // Helper that runs one update on a fresh filter with the constructed measurement.
    auto run = [&](Scalar threshold) -> bool {
        Eskf f;
        f.init(SE3{}, Mat12::Zero());   // P_pose = 0 -> S = H P H^T + R = R = I6
        Measurement m;
        m.dim = 6;
        m.H.setZero();
        m.H.block<6, 6>(0, 0) = Mat6::Identity();   // H = [I6 | 0], observes the full pose
        m.R = Mat6::Identity();                     // R = I6 -> S = I6
        // residual with squared norm == target_nis (all mass in the first component).
        m.residual.setZero();
        m.residual(0) = std::sqrt(target_nis);
        return f.update(m, threshold);
    };

    // Sanity: with P=0 and S=I6 the NIS is exactly target_nis.
    {
        Eskf f;
        f.init(SE3{}, Mat12::Zero());
        Measurement m;
        m.dim = 6;
        m.H.setZero();
        m.H.block<6, 6>(0, 0) = Mat6::Identity();
        m.R = Mat6::Identity();
        m.residual.setZero();
        m.residual(0) = std::sqrt(target_nis);
        f.update(m, 1e9);   // accept so last_nis is the computed NIS
        CHECK(f.last_nis() == doctest::Approx(target_nis).epsilon(1e-9));
    }

    // Flat base (the OLD single-scalar behavior): NIS > base -> REJECTED.
    CHECK_FALSE(run(base));
    // Per-n(6) gate (the NEW behavior): NIS <= gate6 -> ACCEPTED.
    CHECK(run(gate6));
}
