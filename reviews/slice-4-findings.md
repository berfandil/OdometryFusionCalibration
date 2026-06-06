# Slice 4 (histogram primitive) — code review

Commit `e6b5dc2` on `main`. Files: `include/ofc/core/histogram.hpp`, `src/core/histogram.cpp`, `tests/test_histogram.cpp`, `include/ofc/core/config.hpp` (HistogramConfig). Contract: DESIGN §6, DECISIONS D9, CONFIG §8. No code changed.

One finding per line. severity ∈ {CRITICAL, MAJOR, MINOR, NIT}.

## Aging — SlidingK (float drift)
src/core/histogram.cpp:148: MAJOR - SlidingK forgetting subtracts stored deposits via apply_vote(.., -1); float addition is non-associative so individual bins_[] accumulate drift and can settle slightly negative over many evictions. No per-bin clamp/rebuild exists; only total_ is patched (and only at ring_count_==0). Periodically rebuild total_=sum(bins_) and clamp each bin's tiny-negative residual to 0 (e.g. when |bins_[i]| < kEmptyEps), or rebuild from the ring.
src/core/histogram.cpp:194: MINOR - total_ desync guard (`if (ring_count_ == 0) total_ = 0`) is effectively dead in steady saturated operation: once ring_count_ reaches sliding_k_ it never returns to 0, so accumulated total_ float noise is never scrubbed. Clamp total_ to 0 whenever total_ <= kEmptyEps, or recompute total_ from bins_ each evict.

## Confidence (doc vs impl)
src/core/histogram.cpp:252: MINOR - confidence = (peak + two neighbours)/total. DESIGN §6 line 102 states "Confidence = peak concentration (peak mass ÷ total)" — i.e. peak bin only. Impl matches the review brief/header but contradicts DESIGN.md. Reconcile DESIGN.md §6 wording with the implemented 3-bin concentration (or vice versa).

## Binning / wrap
src/core/histogram.cpp:87: MINOR - non-circular fold() clamps inclusively to [range_min, range_max] while the binning math assumes a half-open [range_min, range_max); a value exactly == range_max yields f = nbins-0.5 (top edge of the last bin) rather than being rejected/treated symmetrically with the bottom edge. Behaviour is in-bounds (lround/floor both resolve into bin nbins-1) but the half-open/closed asymmetry is undocumented. Note it, or clamp the top to range_max - tiny.
src/core/histogram.cpp:93: CLEAN - wrap_bin / mass_at neighbour indexing: circular uses i%nbins_ with negative fix-up, non-circular clamps (wrap_bin) or returns 0 (mass_at). No OOB, no seam off-by-one; min bins=4 prevents peak±1 neighbour aliasing/double-count in confidence.

## Vote split
src/core/histogram.cpp:131: CLEAN - linear split fl=floor(f), frac=f-fl, w_lo=w(1-frac)+w_hi=w·frac sums to weight; wrap-aware at the seam (lo=-1→nbins-1, hi=nbins→0); single-bin path (vote_split==false) deposits full weight at lround(f) nearest bin.

## Sub-bin peak
src/core/histogram.cpp:236: CLEAN - parabolic offset 0.5(hm-hp)/(hm-2h0+hp) with kParabEps denom guard, clamp to [-0.5,0.5], fold() of the circular result; boundary fallback to bin center on a missing non-circular neighbour is the documented, acceptable behaviour.

## Aging — Decay
src/core/histogram.cpp:164: CLEAN - every add scales all bins_ and total_ by decay_gamma before deposit; mode tracks a shifting distribution; steady-state total bounded by weight/(1-gamma).

## Empty / degenerate
src/core/histogram.cpp:218: CLEAN - total_<=eps → mode = midpoint (range_min+range_max)/2, peak_bin -1, confidence 0; all-equal and single-vote paths well-defined.

## Strict core
src/core/histogram.cpp:45: CLEAN - no heap after (or in) configure(); fixed bins_[4096]/ring_[4096]; loops bounded by nbins_ or sliding window; Status returns; configure range-checks bins∈[4,kMaxBins], range_max>range_min, gamma∈(0,1), sliding_k∈[1,kMaxSlidingK].

## Determinism
src/core/histogram.cpp:30: CLEAN - all double math, no RNG, no time/IO; results are deterministic for a given add() sequence.

## Test quality
tests/test_histogram.cpp:273: MINOR - the circular wrap-seam test pins only dist_to_seam<0.1 and |m|>3.0; it does NOT pin an analytic sub-bin offset / split fraction AT the seam. Add a seam case with planted neighbour masses asserting the exact parabolic offset and the wrap-aware split fractions.
tests/test_histogram.cpp:230: MINOR - SlidingK tests run only ~45 adds (K=20) and ~4 adds (K=3); the negative-bin / total-desync drift over many evictions (the MAJOR above) is never exercised. Add a long SlidingK run (e.g. >10*K mixed-value adds) asserting all bins_ >= 0 (or > -eps) and total_ == sum(bins_).
tests/test_histogram.cpp:348: NIT - the [0,1] confidence stress loop uses the default Decay config, not SlidingK, so it does not cover SlidingK confidence/total over many evictions either.
tests/test_histogram.cpp:381: CLEAN - weight<=0 (0 and -1) ignored is tested; add-before-configure no-op tested; vote-split fractions, parabolic offset, decay A→B migration, SlidingK saturation, circular neighbour wrap, and boundary missing-neighbour are all pinned to analytic values.

TOTAL: 16 findings (0 CRITICAL, 1 MAJOR, 6 MINOR, 1 NIT, 8 CLEAN)
