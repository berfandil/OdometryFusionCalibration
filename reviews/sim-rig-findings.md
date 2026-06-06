# Simulation rig (ground-truth oracle) — code review

Commit `4c186b7` on `main`. Files: `sim/include/ofc_sim/{trajectory,synthetic_source,rig}.hpp`, `sim/src/{trajectory,synthetic_source,rig}.cpp`, `tests/test_sim.cpp`, `sim/CMakeLists.txt`, root + `tests/CMakeLists.txt`. Contract: DESIGN §6/§10, DECISIONS D20/D21/D24, `include/ofc/core/{types,source,estimator,lie}.hpp`, `src/core/estimator.cpp`. No code changed. Built + ran: 17 cases / 1579 assertions pass.

One finding per line. severity ∈ {CRITICAL, MAJOR, MINOR, NIT}.

## Measurement model — frame-align inverse, scale, offset
sim/src/synthetic_source.cpp:45: CLEAN - B = X^-1 o A o X is the exact inverse of the estimator's A = X o B o X^-1 (src/core/estimator.cpp:196); identity X reproduces A and planted X round-trips (pinned, test_sim.cpp:143-167). X = sensor->base matches D21.
sim/src/synthetic_source.cpp:47: MAJOR - the planted `scale` multiplies B.t but the estimator (src/core/estimator.cpp) reads `prior_extrinsic` only and NEVER de-scales by `prior_scale`. So "estimator-with-prior == planted recovers true base motion" is FALSE whenever scale != 1: the rig end-to-end claim only holds for X, not scale. Pre-median scale-correction is unimplemented (DESIGN §2 defers it). Document this gap in synthetic_source.hpp and add a `scale != 1` rig case once the estimator de-scales, else the oracle silently over-claims.
sim/src/synthetic_source.cpp:34: MAJOR - `time_offset_s` sign convention is undocumented in the contract (CONFIG.md `prior_time_offset_s` gives no sign meaning; the param comment says "leads/lags" — both). The oracle DEFINES the offset that Slice-5 xcorr must invert: positive off shifts the SAMPLED window later (ts=t0+off, reads a later trajectory slice). Pin the sign in DESIGN §6 / CONFIG so Slice-5 recovers it with the right sign; an oracle with an unstated convention can validate a sign-flipped time-sync as "correct".
sim/src/synthetic_source.cpp:47: MINOR - scale is applied to B.t AFTER the X-conjugation (scale of the sensor-frame delta, post hand-eye), consistent with D20 "per-source translation scale". The composition order is fine but undocumented re: whether scale precedes or follows the mount; state it.

## Noise model + covariance
sim/src/synthetic_source.cpp:112: MAJOR - `eps << nt(gen), nt(gen), nt(gen), nr(gen), nr(gen), nr(gen)` — the relative evaluation order of the six generator calls in an Eigen comma-initializer is unspecified pre-C++17, and this target is C++14 (sim/CMakeLists.txt:14). Within one toolchain it is stable (replay tests pass) but the noise draw is NOT portable across compilers, undermining the §10 byte-stable golden-regression requirement. Draw into a local `double a=nt(gen); double b=nt(gen); ...` in explicit statement order, then assign.
sim/src/synthetic_source.cpp:121: MINOR - `modeled_cov_(B)` is computed from the NOISY B (line 114 already reassigned B), so the reported sigma uses the perturbed dist/angle, not the clean dist/angle actually used to draw the noise (lines 102-106). For small noise negligible, but the cov is then not exactly the injected sigma. Compute cov from the clean delta (move modeled_cov_ before the perturb, or pass clean dist/angle).
sim/src/synthetic_source.cpp:121: MINOR - on an OUTLIER window the cov is synthesized from the gross-wrong delta → huge dist/angle → large variance → sigma_confidence down-weights the outlier in the median. The "median rejects the outlier" rig test (test_sim.cpp:397) thus partially passes via Σ-downweighting, not pure geometric-median rejection. Report a fixed/normal cov on outlier windows so the test isolates median robustness.
sim/src/synthetic_source.cpp:114: CLEAN - noise is zero-mean (normal_distribution mean 0), right-perturbed in the body tangent B o exp(eps), in [trans;rot] order matching Delta.cov and types.hpp Vec6 = [v;omega].
sim/src/synthetic_source.cpp:108: MINOR - `std::normal_distribution`/`mt19937_64` output is not specified by the standard across stdlib implementations; golden output is byte-stable only on a fixed stdlib. Acceptable for the dev toolchain but note it as a golden-regression caveat (§10).

## Determinism
sim/src/synthetic_source.cpp:17: CLEAN - mix_seed is a pure function of (seed, t0, t1) via a splitmix-style mix; noise on a window is independent of call sequence → replay-invariant to tick alignment, as claimed. No map iteration, no unseeded RNG, no time/IO. (One residual cross-toolchain caveat: the comma-init order MAJOR above.)

## Trajectory — integration + presets
sim/src/trajectory.cpp:31: CLEAN - pose(t) = cum_pose_[k] o exp(xi*tau) with cumulative boundary poses precomputed by the same se3::exp the estimator integrates with → continuity across boundaries is exact by construction (pinned, test_sim.cpp:41). Clamp before t0 / after end is correct (no motion extrapolation).
sim/src/trajectory.cpp:88: CLEAN - twist(t) returns the active segment's body [v;omega]; matches a finite-difference log(p0^-1 p1)/2h of pose(t) away from boundaries (pinned 1e-3, test_sim.cpp:60).
sim/src/trajectory.cpp:130: CLEAN - straight(): omega=0, v>0 (pinned). turning(): wz>0 → ||omega||>0 (pinned). Both realize their regimes per DESIGN §6.
sim/src/trajectory.cpp:144: CLEAN - omega_varying() alternates straight / left / hard / right with a 0.5*seg straight gap → ||omega||(t) has a distinctive non-constant shape (min<0.05, max>0.5 pinned, test_sim.cpp:99) usable for the Slice-5 cross-correlation.
sim/src/trajectory.cpp:80: MINOR - twist_s uses `rel >= total_s_` → zero at exactly t=end, but pose_s clamps to the final pose; the twist is discontinuous (drops to 0) one tick before pose stops moving. Harmless for queries but means twist(end) != the integrand of pose near end. Note the half-open convention.
sim/src/trajectory.cpp:64: NIT - segment_at_ is a linear scan per query; fine for short sim trajectories (documented), but the rig run() calls pose() twice per tick over many ticks → O(ticks*segments). Acceptable; no action.

## Rig — gauge anchoring + error metric
sim/src/rig.cpp:65: CLEAN - GT anchored at pose(frontier - window_ns), which equals the estimator's bootstrap q0 = t1 - window_s (src/core/estimator.cpp:167) exactly (same stamp, same rounding). The estimator predicts from identity over [q0,t1], so its origin IS the GT pose at that point — principled, NOT a fudge; later fuses chain gap-free from the same origin. The "0.2 m error" symptom was the genuine window-length offset, correctly fixed here.
sim/src/rig.cpp:78: CLEAN - pose_error compares fused.frontier.pose and gt_frontier both in the anchored odom frame; trans = ||dt||, rot = ||log(gt.R^T fused.R)||. Consistent frame; zero iff fused == GT.
sim/src/rig.cpp:61: MINOR - the anchor is captured from the FIRST fused record's frontier; if the very first step's source coverage differs across sources/offsets, the anchor stamp shifts. In the current rig (offset=0 sources) this is stable, but with planted offsets the per-source windows no longer share the frontier and the single-anchor assumption is untested. Add a planted-offset rig case (and pin which frontier defines the anchor).

## Injection
sim/src/synthetic_source.cpp:85: CLEAN - dropout returns Status::NoData; outlier integrates a fixed gross body twist; early-return on dropout before noise; outlier vs clean windows do not leak (pinned, test_sim.cpp:234-265).
sim/src/synthetic_source.cpp:71: MINOR - fault membership is decided by the offset-shifted MIDPOINT only → a query straddling a window edge is all-or-nothing (no partial-window fault). Acceptable for an oracle but undocumented; a window boundary landing mid-query silently flips the whole delta.

## CMake
CMakeLists.txt:57: CLEAN - sim built when OFC_BUILD_SIM OR OFC_BUILD_TESTS; ofc_tests links ofc_sim (tests/CMakeLists.txt:18); sim is a relaxed edge with no strict warning flags forced (sim/CMakeLists.txt:13). Matches D24 / §9.

## Test quality
tests/test_sim.cpp:127: CLEAN - analytically PINS: identity reproduces GT delta (1e-9); planted-X conjugation + A-recovery (1e-9); scale on translation only (1e-9); offset-shifted window AND not-equal-to-unshifted (offset observable); noise zero-mean/bounded/bit-deterministic + different-seed; dropout NoData; outlier gross; buffer-backed == analytic; rig tracks GT; deterministic replay WITH noise; outlier rejected. These would catch a transposed extrinsic, a no-op offset, and non-determinism.
tests/test_sim.cpp:314: MAJOR - NO end-to-end rig case plants scale != 1 or time_offset != 0: every rig SourceParams uses the defaults. Combined with the estimator not de-scaling, a WRONG scale handling (or a sign-flipped offset) passes the whole rig suite. Add rig cases with planted scale and planted offset once the estimator consumes them.
tests/test_sim.cpp:201: MINOR - the noise test asserts zero-MEAN only implicitly (single draw bounded < 0.2); it never averages many windows to assert the empirical mean ≈ 0. A non-zero-mean noise (a planted bias) would pass. Average eps over many windows and assert ||mean|| small.
tests/test_sim.cpp:143: MINOR - planted-X test uses one window on turning(); does not assert the conjugation under translation-only motion where X.t matters most distinctly from X.R. Add a straight-motion planted-X case to disambiguate a transposed vs inverted X.t.

TOTAL: 27 findings (0 CRITICAL, 5 MAJOR, 11 MINOR, 2 NIT, 9 CLEAN)
