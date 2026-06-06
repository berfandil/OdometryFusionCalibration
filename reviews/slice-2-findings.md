# Slice 2 (median fusion) review — commit 7828e4d. Strict code review; severity ∈ {CRITICAL, MAJOR, MINOR, NIT}.

## Median (median.cpp / median.hpp)
src/core/median.cpp:85: MINOR - n==2/n>=3 spread is computed from the RAW `weights` array, but the median itself uses `w_of()` (uniform fallback when all weights<=0). When all weights<=0 the median is the uniform geodesic/IRLS solution yet spread returns 0, understating Q. Pass the same effective weights (uniform-corrected) to weighted_rms_spread.
src/core/median.cpp:85: MINOR - weighted_rms_spread reads raw `weights[i]` with no non-negativity guard, so a negative weight subtracts from acc_d2/acc_w and can yield a NaN/garbage spread. Clamp w to >=0 inside the RMS like w_of does.
src/core/median.cpp:131: NIT - convergence step_norm mixes rotation (rad) and lambda*translation under one tol, so `tol` is not a pure-radian threshold; document that tol is in split-metric units (matches solver intent, just unstated).
src/core/median.cpp:107: NIT - n>=3 has no early-out when all inputs are identical; it runs to a tol-triggered break on iter 1 (fine), but spread==0 inputs still cost one full pass. Acceptable; noted only for WCET accounting.

## ESKF (eskf.cpp / eskf.hpp)
src/core/eskf.cpp:62: MINOR - F's twist block is zeroed and Qmap is block-diagonal, so the pose<->twist cross-covariance is forced to 0 every predict even though twist=log(delta)/dt is a deterministic function of the same delta that moves the pose (they ARE correlated). Tip extrapolation (pose+twist) is mildly over/under-confident as a result. Acceptable for Slice 2 but flag: a future correct model couples them via the readout Jacobian.
src/core/eskf.cpp:55: MINOR - twist = log(delta)/dt_eff but the published twist.cov = q_pose/dt^2 only (the fresh-readout noise); it ignores that delta's own posterior pose covariance also feeds the readout. Consistent with "re-read each step" intent; record as a modeling approximation, not a bug.
src/core/eskf.cpp:50: NIT - F = blkdiag(Ad(delta^-1),0) for T<-T∘delta with right error eta (T_true=T∘Exp(eta) => eta+=Ad(delta^-1)eta) is CORRECT; ordering [trans;rot] matches lie.cpp adjoint. Verified, no action.
src/core/eskf.cpp:74: NIT - internal state_.stamp accumulates secs_to_ns(dt_eff) from 0 and diverges from the real frontier t1 (estimator overwrites frontier.stamp=t1 and tip.stamp=now, so it's masked). Either drive the eskf stamp from the estimator or drop it to avoid a latent inconsistency.

## Estimator (estimator.cpp)
src/core/estimator.cpp:212: MAJOR - predict() is always called with dt=cfg.window_s regardless of the actual spacing between consecutive frontier t1 values. Nothing enforces that the caller steps exactly one window per tick; if tick cadence != window_s the integrator overlaps/gaps and the fused pose drifts. Either derive dt from (t1 - last_t1) or document+validate that step cadence must equal window_s.
src/core/estimator.cpp:200: MINOR - adaptive_q hardcodes q_scale=1.0 and q_floor=1e-6 (six axes); non-adaptive branch hardcodes 1e-4. CONFIG §3 defines q_scale/q_floor as Config fields. KNOWN: to be fixed in the follow-up slice — recorded per review scope.
src/core/estimator.cpp:48: MINOR - sigma_confidence uses 1/(mean_diag+eps) over all 6 diagonal entries, mixing translation (m^2) and rotation (rad^2) variances into one scalar; a source with tiny rad var but huge m var gets a misleading confidence. Acceptable placeholder (reliability EMA is Slice 9) but note the unit-mixing.
src/core/estimator.cpp:201: NIT - adaptive Q is isotropic (same spread^2 on all 6 axes) though spread folds rotation+lambda*translation together; trans and rot Q are not separable. Minor modeling choice, consistent with split-metric spread.
src/core/estimator.cpp:190: CLEAN - the residual loop re-walks with an independent k incremented only on in_window, correctly matching the packed aligned[] indexing vs health[i] slot indexing. No off-by-one.
src/core/estimator.cpp:36: CLEAN - kMaxSourcesCap=32 matches Result::health[32]/calib[32]; add_source caps on both cfg.max_sources and the cap; validate() rejects max_sources>32. Bounds consistent.

## Strict core (no-heap / bounded / status)
src/core/estimator.cpp:131: CLEAN - step() allocates nothing (all scratch in fixed Impl arrays; median uses caller arrays; eskf uses fixed Eigen); the sole heap is the allocate-once Impl at init(). Loops bounded by source_count and weiszfeld_max_iters. Status returns correct (NotInitialized/NotReady/Ok).

## Determinism
src/core/median.cpp:98: CLEAN - tie-break in the highest-weight start uses strict '>' (keeps first max); accumulation order is fixed by index; double, no fast-math. test_fusion replay asserts bit-identical cov. No ordering/uninitialized-read nondeterminism found.

## Config
include/ofc/core/config.hpp:74: MINOR - confidence_blend added correctly (default 0.5, CONFIG §4) but validate() does not range-check it to [0,1]; likewise weiszfeld_tol/eps, fusion_delay_s, weight_floor/cap, tip_cov_inflation are unvalidated (TODO at estimator.cpp:257 acknowledges). Add the missing bounds.

## Test quality
tests/test_eskf.cpp:90: MAJOR - covariance test only asserts symmetric + PSD + trace>tr0; it never pins P to a closed-form value, so a WRONG F (e.g. Ad(delta) instead of Ad(delta^-1), or an identity twist block) would still pass. Add a one-step analytic check of P = F P0 F^T + Q against a hand-computed matrix.
tests/test_fusion.cpp:178: MINOR - outlier-rejection test pins the fused pose to GT within 5e-2 but never asserts the outlier source's health.residual is large / its weight low; it confirms the OUTCOME but not that rejection (not luck/averaging) caused it. Add a residual/weight assertion on the biased source.
tests/test_median.cpp:50: MINOR - no test for n==0 (identity, spread 0) or all-zero/negative weights (uniform fallback); the degenerate-count and weight-guard paths in solve() are unexercised. Add cases.
tests/test_fusion.cpp:100: NIT - GT is integrated with the same const-twist composition the buffers use, so the test validates the pipeline against its own integrator, not an independent closed form; tolerance 1e-2 hides small model error. Acceptable for a tracer bullet.
tests/test_eskf.cpp:60: CLEAN - twist readout, pose composition, multi-step accumulation, adaptive_q floor/scale, and tip const-velocity are pinned analytically. Good coverage of the mean path.

TOTAL: 0 CRITICAL, 2 MAJOR, 9 MINOR, 6 NIT (+5 CLEAN categories noted).
