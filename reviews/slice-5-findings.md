# Slice 5 (time-sync) review — commit e094920

Scope: `timesync.hpp/.cpp`, `test_timesync.cpp`, `estimator.cpp` + `config.hpp` diffs, against DESIGN §2/§6, DECISIONS D16/D20/D21, CONFIG §5/§8/§9, `histogram.hpp`. Verified against the oracle `synthetic_source.cpp::base_delta`.

## Sign convention & xcorr (load-bearing) — CLEAN
- Oracle: query `[t0,t1]` reads trajectory `[t0+off,t1+off]` ⇒ `omega_src(t)=omega_ref(t+off)` (positive off = source clock ahead). Matches header lines 15-26 and D21 line 104.
- xcorr aligns `ref[i]` vs `src[i-lag]`; with `src[g]=ref[g+p]`, `p=off/dt`, peak sits at `lag=p` ⇒ recovered `=lag*dt=off`, positive→positive. CORRECT.
- Estimator removes it via `query(q0-off, t1-off)` (estimator.cpp:287) ⇒ source's internal +off lands back on `[q0,t1]`. Cancels, does not double. CORRECT.
- Four metrics all normalized to maximize (cost metrics negated; NCC zero-mean/variance-normalized, `denom<kRatioEps` guard at timesync.cpp:239; Ratio denom `+kRatioEps` at :253). CORRECT.

## Findings
src/core/estimator.cpp:153,363: MAJOR - Commit gate checks only `confidence >= commit_concentration`; DESIGN §2/§6 require commit gated ALSO on `votes >= commit_min_votes` (N_min) AND hysteresis (re-open below `commit_drop`). With SlidingK=64 a few sharp votes commit instantly; a transient spurious peak can drive fusion. Add an N_min check (`hists_[s].total() >= commit_min_votes`) and drop/re-open hysteresis state, or document the deferral in DESIGN.
src/core/estimator.cpp:424: MINOR - `validate()` never range-checks `excitation_min_var` (>=0), `max_lag_s` (0,2], or `match_metric`. Only enforced inside `timesync.configure()` and only when `timesync_enabled`; a caller validating config standalone (or with timesync off) misses them. Add the §5 checks to `validate()`.
tests/test_timesync.cpp:233-266: MINOR - The sim e2e PINS the sign (`est>0`) but NOT sub-sample accuracy: tol is `1.5*dt` (15 ms on a 30 ms plant) and `omega_varying()` is piecewise-constant (step edges → flat-topped xcorr, weak parabola). Sub-sample is only pinned on the smooth `shape()` unit test (:135, `0.4*dt`). Acceptable split, but add one sim e2e assertion at sub-sample tol on a smooth-‖ω‖ trajectory to guard the refine end-to-end.
tests/test_timesync.cpp:158-174,135-156: MINOR - Negative-lag recovery is exercised under L2 only; fractional under L2+NCC only. A sign/refine regression in L1 or Ratio for negative/fractional lags would pass. Loop all four metrics over a negative and a fractional plant (the integer test already covers all four).
src/core/timesync.cpp:296-302: NIT - Parabolic refine guards `|denom|>kParabEps` but not `denom<0` (concave). On a flat-topped plateau (`s0` ties a neighbour) `denom` can be >0 and δ points toward a minimum; bounded only by the ±0.5 clamp. Harmless given the clamp, but a `denom<0` precondition would make the concavity assumption explicit.
src/core/estimator.cpp:259-282: NIT - ‖ω‖ is sampled at the source's RAW timeline with the offset deliberately NOT applied (correct — that skew is the signal). The query interval is `[t1-h, t1]` (one sample_dt), so the per-tick ‖ω‖ is a backward-difference at the frontier; on a piecewise-constant trajectory the window straddling a step yields a ramp, not the instantaneous magnitude. Fine for xcorr (both channels see the same smoothing), worth a comment.
include/ofc/core/timesync.hpp:152-156: NIT - Default member initializers (sample_dt_=0.02, dt_ns_=2e7, max_lag_samp_=5) are dead — `configure()` always overwrites them before any use, and every accessor short-circuits on `!configured_`. Cosmetic only.

## Confirmed-correct categories (no action)
- Excitation gate: variance>min on BOTH ref+src overlaps (timesync.cpp:284-286); `straight()` (zero ω) and flat unit signal both correctly gated (tests :179, :268). CORRECT.
- ‖ω‖ extraction `‖log(ΔR)‖/h` extrinsic-invariant, sub-interval consistent, common-grid resample with latest-wins overwrite + carry-forward gap fill + bounded channel-restart on > kMaxSamples jump (timesync.cpp:123-174). Sound and bounded.
- Histogram reuse: offset_hist range ±max_lag, mode/confidence/empty→0 special-case (timesync.cpp:325-338). CORRECT.
- Estimator wiring: offset applied only when confident else prior (effective_offset); `timesync_enabled==false` ⇒ exact Slice-2 (test :382 passes prior through, conf 0, uncommitted); offset+confidence surfaced in CalibSnapshot (:355-364). CORRECT.
- Common-offset unobservability: reference pinned at 0, all others relative to it — handled gracefully (no absolute anchor assumed).
- Strict core: fixed-capacity rings, bounded lag scan `[-max_lag,+max_lag]`, bounded loops, Status returns, no heap after configure(); new `excitation_min_var` field added per CONFIG §5. CORRECT.
- Determinism: replay test (:293) asserts bit-identical; double math, deterministic loops. CORRECT.

TOTAL: 7 findings (0 CRITICAL, 1 MAJOR, 3 MINOR, 3 NIT)
