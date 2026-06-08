# Review findings — commit 1142e41 (true interior robust median + covariance recalibration)

Reviewer: independent deep review. Method: full diff read + clean build (`dev.ps1 -Task clean` then build) +
ctest 100% green (unit 211 cases / 8995 assertions, adapters green) + targeted throwaway experiments
(old-solver reconstruction, scale_calib probe) since reverted. Working tree clean; nothing committed.

Format: `path:line: <emoji> <severity>: problem. fix.`

## Re-baselines — was any test WEAKENED to mask a regression? (highest priority)

reviews/slice-median-fix-findings.md:0: ✅ NONE FOUND. No re-baseline masks a real regression. Each weakened bound was traced to a real, measured, understood mechanism (verified empirically below); the observability/safety INTENT is preserved in every case.

tests/test_sim.cpp:448: ✅ CLEAN (root-cause REAL, verified): flipping `scale_calib=false` is a legitimate isolation, not a mask. PROBE: I temporarily set `scale_calib=true` in `make_rig_config`, rebuilt, ran — `max_te` jumped to 0.167714 (exactly the cited ~0.17) on BOTH "tracks GT >=3 sources" AND "tracks GT scale!=1"; identical magnitude across both confirms it is the scale-calib feedback corrupting prior_scale, not the planted scales. calibration.cpp is untouched by this commit, so the spurious ~0.984 commit is a PRE-EXISTING latent bug newly EXPOSED by the correct interior median (the old pinning median masked it by returning the clean reference vertex). This is "latent bug newly exposed", not "the fix degraded scale calibration."

tests/test_sim.cpp:448: ✅ CLEAN (no shared-helper coverage lost): all 5 `make_rig_config` users (lines 483, 522, 573, 633, 673) set `prior_scale == planted`, so NONE of them ever exercised ONLINE scale CALIBRATION — they test de-scaling-via-prior_scale (estimator.cpp:976, independent of the `scale_calib` flag) + fusion tracking. The "scale != 1" test (501) still de-scales {1.2,0.8,1.5} correctly with the calibrator off (verified: it passes). Dedicated online-scale-calibration coverage is intact and GREEN in test_calib_feedback.cpp:186 (wrong prior 1.0 -> commits ~1.18), test_calib_phase1/phase2, test_smoother_calib, test_weights. No coverage silently lost.

tests/test_weights.cpp:351: ✅ CLEAN (mechanism REAL, verified): reliability threshold 0.7x->0.5x clean. MEASURED biased rel=0.652511 (clean=1.0, floor=0.2). The interior median now BLENDS the biased source into the consensus, adding zero-mean residual scatter -> 0.65. Intent (D17: constant bias is KEPT not collapsed to floor) holds decisively: 0.65 >> floor+0.3. Not a relaxation-to-pass.

tests/test_weights.cpp:369: ✅ CLEAN (mechanism REAL, verified): scale-recovery >1.1->1.05. MEASURED converged scale=1.07232 (planted 1.2). Genuine partial recovery from consensus feedback contamination (the median blend the frame-align loop feeds back partially masks the residual). New bounds keep a DIRECTIONAL lower bound (>1.05, moved off 1.0 toward truth) + a NEW upper bound (<1.25) — tighter in one axis, not just loosened. Observability intent (bias routed to calibrator, scale moves toward truth) preserved. NOTE: the test comment says "the consensus the calibrator votes AGAINST is pulled toward the biased source"; the scale vote is actually `bn/ref_mag` vs the REFERENCE (calibration.cpp:295), not vs the median — the contamination is via the rising-edge prior_scale feedback loop, not a direct median vote. Imprecise wording, correct conclusion.

tests/test_validation.cpp:445: ✅ CLEAN (re-measured): NEES band [4.0,5.6], chi2 99% interval on the mean = [5.94277, 6.05754], measured ensemble-mean = 4.8189 (N=24180). 4.82 < 5.94 so `CHECK_FALSE(truly_consistent)` is HONEST; band brackets 4.82 with ~0.8 margin each side; upper 5.6 < DOF=6 keeps never-overconfident.

tests/test_cov_calibration.cpp:184: ✅ CLEAN (re-measured, all 3 trajectories): nees_traj=4.8189, turning=2.00051, straight=4.84321 — every one < overconf_cap=5.5 and < DOF. q_scale=1 -> 3.38005 (the D4 sign-flip is real: smaller q_scale -> larger NEES). NIS re-measured = 2.72697 (DOF=3, band [1.5,3.5]).

tests/test_validation.cpp:781: ✅ CLEAN (golden genuinely unchanged, verified): noise-free golden passes byte-identical on the new code — all kGoldenPose samples within kAbsTol, c1->scale=1.08281 holds to 1e-4, fused-vs-GT max_te/max_re < 1e-9. Its `vote_weight=One` config (golden_config does NOT disable scale_calib) tracks to machine precision even with the calibrator on + the interior median, confirming the implementer's claim. No regen needed.

## Median fix correctness

src/core/median.cpp:128: ✅ CLEAN: off-vertex weighted-mean init (translation arithmetic mean + one Karcher rotation step about highest-weight R) is correct, deterministic, allocation-free (fixed-size Eigen locals; so3::log/exp are heap-free per lie.cpp). n==1 (87) / n==2 (96) closed forms untouched. VZ guard `if (d<=eps) continue` (156) correctly excludes a coincident self-term; the all-coincident edge -> usum==0 -> converged + returns the mean init (correct).

tests/test_median.cpp:165: ✅ CLEAN (convergence-flag edge is a legit fix, verified): PROBE on the `converges` input — tol=1e-4 reaches the step in 25 iters, tol=1e-6 in 77, tol=1e-8 in 132. The old test (tol=1e-8, max_iters=20) demanded a step unreachable within its cap once the solver genuinely iterates off-vertex (132 > 20). Moving to production tol=1e-6 / max_iters=100 (77<100) is a real fix, not papering over non-convergence. The VALUE converges to (1,0,0) to ~1e-7 at every tol.

## High-weight-outlier guard (new) — is it NON-VACUOUS?

tests/test_median.cpp:201: ✅ CLEAN (NON-VACUOUS, proven by reconstruction): I rebuilt the OLD pinning solver (vertex init, no VZ guard) — it returns t=(50.00000,-40.00000,7.00000), the outlier VERBATIM, so the test's `CHECK(t.x ~= 1.0)` would FAIL against it. The NEW solver returns t=(1.01,-0.0004,0.0001), iters=61 — rejects the outlier. The outlier carries weight 2.0 but the inliers are the weighted MAJORITY (4x1.0=4.0 > 2.0), so rejecting it is correct (a >50%-mass point would legitimately BE the median). WCET case at the production cap (max_iters=10): iters=10, conv=false, t.x=1.02879 — flag not latched (Weiszfeld linear tail) but VALUE correct, far from the outlier. Guard + interior-not-pinned assertions are meaningful.

## /n_eff removal completeness

include/ofc/core/eskf.hpp:205: ✅ CLEAN: adaptive_q back to 3-arg everywhere — decl (eskf.hpp:205), defn (eskf.cpp:368), all callers (estimator.cpp:1292/1294, test_eskf.cpp). No 4-arg caller anywhere (incl. adapters); clean build links.

include/ofc/core/config.hpp:107: ✅ CLEAN: `adaptive_q_source_reduction` removed from config.hpp; zero references in include/ src/ adapters/. `chi2_gate` (per-n gate, eskf.cpp/test_chi2_gate.cpp) untouched (not in diff).

src/core/estimator.cpp:156: ✅ CLEAN (hash-neutral CONFIRMED): config_hash never hashed the removed knob — the global section (156-194) hashes no q_* fields; the per-sensor section (201-217) hashes scale_calib/bias_states/etc. but never adaptive_q_source_reduction. validate() never referenced it. Removal does not change any persisted config hash.

## WCET / strict-core

src/core/median.cpp:129: ✅ CLEAN: off-vertex init is one bounded O(n) pre-loop with only fixed-size Eigen locals (no heap, no new allocations); main loop still hard-capped by max_iters (10 in production). Strict-core preserved.

## Scope / housekeeping

.deps_cache/: ✅ CLEAN (not committed): commit 1142e41 contains only the 11 expected files; `.deps_cache/` is untracked (`??`) in the working tree and was NOT in the commit. MINOR housekeeping: it is not in .gitignore (build/ and *.obj are), so it could be `git add .`-ed by accident later — consider adding `/.deps_cache/` to .gitignore.

reviews/slice-median-fix-findings.md:0: ✅ No scope creep / no stray TODO-FIXME-HACK introduced in the diff. The 11 changed files all serve D3/D4.

## Documentation drift (the commit changed code+tests but not the canonical specs)

CONFIG.md:40: ⚠️ MAJOR: CONFIG.md still documents the REMOVED `adaptive_q_source_reduction` knob as a live bool config field (full row, default true, formula `q_scale·spread²/max(1,n_eff)`). The field no longer exists in the code — a user setting it per this reference would not compile / would be silently ignored by the loader. Delete the row.

CONFIG.md:39: ⚠️ MAJOR: CONFIG.md q_scale row says default **0.5** with the old sweep rationale (worst-case NEES ~2.9, q_scale=0.2 overconfident). Code default is now 0.7 with the OPPOSITE sign (smaller q_scale -> larger NEES, 0.5 is overconfident). Also CONFIG.md:150 says "now 0.5". Update to 0.7 + the D4 sign/rationale.

DECISIONS.md:28: ⚠️ MAJOR: D4 impl-note describes the now-REVERTED approach A as current — "the adaptive Q now divides the spread term by n_eff", the `adaptive_q_source_reduction` knob, and `q_scale 1.0 -> 0.5`. This is the canonical decision record and the production comments cite "DECISIONS D4" as authoritative; it now contradicts the code. Update to record the D3 median fix + D4 /n_eff removal + q_scale 0.7.

DESIGN.md:143: ⚠️ MINOR: DESIGN §10 NEES/NIS status still describes approach A `/n_eff` as the live fix (q_scale 0.5, NEES band [1.50,2.70], the knob-OFF guard [0.45,1.00]) — all of which this commit deleted. Production comments reference "DESIGN §4". Refresh to the new regime (NEES ~4.8, band [4.0,5.6], no knob-OFF guard).

tests/test_cov_calibration.cpp:172: ⚠️ NIT: comment says "1.6 leaves headroom" but the actual constant `pessimism_lo = 1.3` (175). The per-trajectory lower-bound loop with 1.3 would NOT trip on a q_scale=1 regression (turning ~1.4 > 1.3); the real q_scale=1 regression guard is the ratio CHECK at line 201 (nees_new > nees_old*1.25). Comment overstates what the loop bound catches. Fix the comment (or raise pessimism_lo to 1.6 to match the prose — 1.6 < 2.0 measured turning, still passes).

reviews/slice-median-fix-findings.md:0: ⚠️ MINOR: two genuine pre-existing latent bugs newly EXPOSED by this fix are documented only in test comments, not surfaced to ISSUES.md: (1) the cold-start scale-calib spurious commit (~0.984 = 63/64 coarse-histogram/bootstrap artifact when priors already equal truth — test_sim.cpp:435); (2) the calibrator consensus-contamination under-recovery (scale 1.07 vs 1.2 because the median blends the source-under-calibration — test_weights.cpp:369). The implementer flagged both inline; per the project's surface-to-orchestrator workflow they should be logged as open issues (the second has a candidate fix: exclude the source-under-calibration from its consensus / use ReferenceOnly cold-start).
