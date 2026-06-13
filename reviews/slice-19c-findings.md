# Slice 19c Review: GPS-course heading-drift monitor (split policy layer c)

## Findings

include/ofc/core/heading_monitor.hpp:130: CRITICAL: block MEAN (not median) for staircase reduction breaks robustness. On multi-hour urban canyon (e.g. urban12 multipath), a single bad 60-s block with +30deg outlier pair will shift the mean slope estimate. The spec (1.3) specifies "block medians"; the prototype uses median. A short near-linear 60-s block is NOT guaranteedly free of outliers (multipath spikes, sensor glitches, wheel slip steps can all fit in 60s). The comment claims "mean and median coincide to estimator tolerance" — this is unvalidated and contradicts the event-channel design itself (which exists precisely because clamped deviations DO occur). Fix: use median (per-block SORT is one tiny alloc/free in finalize_block; HeadingMonitor is not a WCET-critical path like the fusion loop).

src/core/heading_monitor.cpp:88-98: MEDIUM: slope reservoir overflow behavior under-documented. When full, the OLDEST slope sample (slope_head position) is evicted. The spec comment (hpp lines 52-56) says "the most recent blocks dominate the baseline-weighted median" — this is only TRUE if baseline-weighted MEDIAN correctly ignores old samples with small baseline weights. On a sustained 10-minute GPS denial, the post-denial pool will have 256 short-baseline slopes (120-125s, tiny weights in the median). This is CORRECT (the median handles it), but add explicit comment: "old short-baseline slopes are naturally down-weighted by the baseline-weighted median, so FIFO eviction is safe."

src/core/heading_monitor.cpp:28-40 (weighted_median): MINOR: the algorithm terminates at the first index where `acc >= half`. For tied weights, this picks the FIRST (smallest-value) slope, not the upper-median. This is consistent with the prototype and correct for slope estimation, but document it: "first-half termination = lower median for ties."

src/core/heading_monitor.cpp:265: CRITICAL: zero-drift handling creates boosting discontinuity. Line 265 returns `cap` when `!(s_i > 0)`, meaning a source with EXACTLY zero drift gets FULL cap while a source with tiny drift gets proportionally less. The boost() is designed as a RANKING (min_score / score_i), so zero-drift should be handled the same way. This breaks the ratio principle: a source at 1e-9 rad/s gets boosted 1e9× less than a source at 1e-18 rad/s. On KAIST this will not surface (no source hits zero), but on noiseless simulators or future multi-FOG systems it will cause chatter. Fix: clamp denominator to max(best_score, 1e-4 rad/s) to keep the ratio well-behaved, or return 1.0 for zero-drift (not special-cased to cap).

src/core/heading_monitor.cpp:77-114: MINOR: block ring management is correct but the saturation transition (when block_count reaches kSlopeReservoir and wraps) lacks an assertion. The code is safe (block_count always increments until cap, then the ring wraps), but add a comment: "block_count is the water level; once it saturates at kSlopeReservoir, the ring wraps and oldest block is evicted FIFO."

src/core/heading_monitor.cpp:149-219 (submit_fix): MINOR: the comment on lines 174-175 is misleading. It says "Per-source forward + yaw increments", but these are differences (unsigned). The cross-val gate line 186 compares SPEEDS (magnitudes), and the forward-motion gate line 187 checks the sign separately. Clarify the comment: "cross-val compares GPS speed vs odometry speed; forward-only gating is separate below."

tests/test_heading_monitor.cpp:486: MEDIUM/RISK: byte-identical pin is WEAK. It runs the estimator with `heading_monitor=false`, so the OFF code path is exercised. However, the TRUE safety case should verify that boost=1.0 (the abstain state: monitor ON but < 120s data) reproduces pre-monitor fusion BYTE-IDENTICALLY. The current test only validates that the monitor doesn't CORRUPT data when OFF, not that abstain mode is byte-identical. Fix: add a second test case with 150+ fixes (forms >= 1 baseline), run with monitor ON (boosts will be non-neutral), then run the SAME config with monitor OFF, verify fusion results byte-match and boosts are exactly 1.0. This is ~50 lines and critical for the "byte-identical OFF" contract.

tests/test_heading_monitor.cpp:355-380 (item 3, telescoping): the test injects turn-correlated slip as extra yaw MUTATIONS (lines 366-369). This simulates the slip but does NOT validate one-sided slip (injection on turn-in, no injection on turn-out). The test is valid and catches the drifter, but a stronger test would inject slip ONLY on the inbound turn and verify the bridged pair still catches it. Current test is conservative (both-sided slip is easier to detect than one-sided); add a comment: "two-sided slip is conservative; one-sided would be a stronger check."

tests/test_heading_monitor.cpp:261-286 (item 1, headline): the assertion `CHECK(on.boost[0] == doctest::Approx(10.0))` pins the boost to the hardcoded cap=10.0. A more robust test would verify the RATIO of weights (clean vs. drifter in the rotation channel). If boost_max changes to 5.0, this PIN will fail even though the ranking (and the fusion effect) is correct. Recommendation: also CHECK that boost[0]/boost[2] > 5.0 (the ratio is what matters for the headline).

## Tests Assessment

- **Item 1 (headline)**: Validates that boost moves fusion via R(0,1). Valid, but PIN is hard-coded to cap=10.0.
- **Item 2 (anchor gates)**: Each gate + positive control. Valid and comprehensive.
- **Item 3 (telescoping)**: Turn-bridging catches drifter. Valid but conservative (two-sided slip only).
- **Item 4 (denial freeze)**: Gap > pairgap closes segment, scores hold exactly. Valid and direct.
- **Item 5 (abstain)**: < 3 sources / no blocks -> boost 1.0. Valid.
- **Item 6 (chatter guard)**: One-fix anomaly doesn't flip; persistent flip does. Valid and strong.
- **Item 7 (byte-identical)**: WEAK. Pins OFF path, not abstain (ON, no boosts yet) path. See critical finding above.
- **Item 8 (config-hash)**: Serialize with two configs, cross-reject. Valid.
- **Item 9 (cov-cal)**: Translation weights match (monitor=ON) vs (monitor=OFF). Valid.

## Verdict

**REQUEST CHANGES** (two critical, one high-risk)

The implementation is WELL-STRUCTURED and the fusion integration is CLEAN. However:

### Critical (must fix before merge):
1. **Block MEAN vs. MEDIAN**: Spec says median; code uses mean. On clean KAIST data this will not surface, but on a multi-hour urban canyon with multipath, a single bad 60-s block will skew the mean slope. The event-channel design itself contradicts the claim that mean and median coincide (event-channel clamping exists because deviations DO occur in short blocks). **FIX**: Use median; it's one per-block sort and HeadingMonitor is not WCET-critical.

2. **Zero-drift boost discontinuity**: Line 265 gives zero-drift sources the full cap (10.0) while tiny-drift sources get proportionally less. This breaks the ratio principle of ranking. The fix is simple: clamp the denominator to max(best_score, 1e-4 rad/s) so the ratio stays well-behaved. Affects noiseless simulators and future multi-FOG systems.

### High-risk (strongly recommended):
3. **Byte-identical pin incomplete**: The test only validates that monitor=OFF doesn't corrupt data. It should ALSO validate that monitor=ON with abstaining boosts (< 120s data, so boosts=1.0) reproduces pre-monitor fusion byte-identically. Add a 150+ fix case; this is ~50 lines and closes the safety loop on the byte-identical contract.

### Minor (nice to have):
- Clarify block-ring saturation logic with an assertion or comment.
- Strengthen the byte-identical test headline to check ratio, not hard-coded boost value.

After fixes: **READY FOR MERGE**.
