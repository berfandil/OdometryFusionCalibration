# Slice 6 review — Phase-1 calibration (yaw/pitch + scale) @ 7d10477

Scope: `calibration.{hpp,cpp}`, `estimator.cpp`/`result.hpp` diffs, `test_calib_phase1.cpp`, against DESIGN §3/§6, DECISIONS D5/D11/D20, CONFIG §6/§7/§8, `histogram.hpp`, oracle `synthetic_source.cpp`. Geometry verified self-consistent: `g_obs = R_Xprior·R_Xtrue⁻¹·(±e_x)`, fold→fwd hemisphere, `δR=rot(e_x→g)`, `f=δR·e_x=g_unit`, `extrinsic.R=δR·R_Xprior`; prior==true ⇒ δR≈I; scale residual ⇒ `cs.scale=prior_scale·residual`. Matches the test's `expected_forward_axis` and oracle `B=X⁻¹∘A∘X; B.t*=scale`.

src/core/calibration.cpp:286-294: MAJOR - antiparallel branch builds a 180° δR via so3::exp(perp*π); rotation_between is then logged by so3::log, which is singular at θ=π (s=θ/(2·sinθ) → ÷0 → huge/NaN, lie.cpp:49). Triggerable with fold OFF or a sideways mount whose g_obs≈-e_x. Guard: skip/clamp the vote when c<0 (cos(e_x,g)<0) instead of voting a π-rotation log.
src/core/calibration.cpp:444: MAJOR - reverse-fold tests g_obs.x<0 against FIXED +e_x, not the consensus sign as DESIGN §6/leaf-default says ("fold into consensus/prior hemisphere"). For a near-sideways mount (forward≈±e_y ⇒ g_obs.x≈0) the x-sign is pure noise, so fwd/reverse samples split arbitrarily across antipodal y peaks — the header's "robust even when fused sign is near zero" claim does not hold off the forward axis. Fold against the per-step base forward sign (sign of fused_trans·prior-fwd projection), or document Phase-1 as forward-mount-only.
src/core/calibration.cpp:444: MINOR - fold uses strict `< 0`, so g_obs.x==0 (exactly perpendicular) is kept un-negated; a degenerate but harmless tie. Use a deterministic tie-break tied to a stable component.
include/ofc/core/calibration.hpp:19,76 / src/core/calibration.cpp:443: MINOR - header/DESIGN say "180°-off second peak" with fold off, but so3_hist range (test ±0.8 rad) clamps the π-norm reverse log into the boundary bins (non-circular clamp, histogram.cpp:95-98) — the "second peak" is edge-clamped mass, not a real bin at π. Doc-only; the test's conf_on>conf_off still holds.
include/ofc/core/calibration.hpp:154: NIT - comment says "prior_R_[slot] is R_Xprior" but the member is `prior_[slot]` (full SE3); stale name from a rename.
include/ofc/core/calibration.hpp:83: NIT - observe() doc references `fused_v` (linear part) but the signature takes `fused_trans`; no `fused_v` parameter exists.
src/core/estimator.cpp:621: MINOR - fused_omega=so3::log(med.value.R)/dt is an angular SPEED gated against straight_omega_max, but fused_trans=med.t is a per-step DISPLACEMENT gated against straight_trans_min; the two gate operands have different time-normalization (rad/s vs m). Consistent with the header and the unit test, but tie straight_trans_min to a speed or document that the displacement threshold is cadence-dependent.

CLEAN categories (verified, no finding):
Observability self-test — LOAD-BEARING assertion is strong: voted==0, vote_count==0, extrinsic_confidence==0, forward_axis==e_x@1e-9, scale==1@1e-9. A weak "runs" pass is impossible here.
Convergence test — pins recovered axis to planted @1e-2 and scale @2e-2 (not just "confidence rose"); a transposed/inverted extrinsic or scale-vs-residual confusion would fail. Reference stays pinned (scale==1, axis==e_x).
Scale guards — ‖B_ref.t‖ guarded by have_ref=(ref_mag>kUnitEps); per-source ‖B.t‖<kUnitEps skipped; ref/scale_calib==false/unvoted pinned at 1.
Gate (D5) — votes only when ‖fused_ω‖<ε ∧ ‖fused_t‖>δ; ref_cross_check folds the ref's own delta to +e_x, requires >30°-tol forward, returns NotReady when ref absent; cross-check on/off both tested.
Estimator wiring — per-DOF confidences added to CalibSnapshot non-breaking; `confidence` stays the Slice-5 time-offset; no feedback into fusion (Slice 8); fusion path unchanged.
Strict core / determinism — all histograms configured at init; observe()/readouts heap-free, bounded by source_count_/3; double math; SlidingK fixed-order; determinism test asserts exact equality.
Scale de-scale wiring — estimator feeds B_corr (B.t/prior_scale) matching pre-median de-scale; calib1.scale recovers RESIDUAL; cs.scale=prior_scale·residual correct.

TOTAL: 7 findings (0 CRITICAL, 2 MAJOR, 3 MINOR, 2 NIT) + 7 clean categories.
