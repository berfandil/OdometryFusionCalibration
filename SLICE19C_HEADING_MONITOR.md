# Slice 19c — GPS-course heading-drift monitor (split policy layer c: auto-discovery of the heading-grade source)

**Goal**: replace the hand-tuned `rot_weight_prior = 10` in the urban recipe with online auto-discovery: rank the per-source heading quality against GPS course-over-ground and boost the best source's rotation-channel weight. Weighting needs only a RANKING (relative drift rates), not an absolute bias estimate — this sidesteps the Slice-18/D28 observability boundary (position-only GPS cannot observe per-source yaw bias absolutely; source-vs-course drift RATES are directly comparable).

**Prototype (validated 2026-06-12, `tools/proto_heading_monitor.py`, `28c7d07`)**: on KAIST urban07/12/17 the monitor auto-discovers the FOG source in **79/79 scored sweep configurations (zero wrong top-1; 2 abstain)**; full 3-way ordering matches the GT reference on 07/17 everywhere (junk-vs-junk order on urban12 flips in 8/25 settings — irrelevant to the weight rule, both clip to the floor). FOG drift estimates land within ~1–3 deg/h of GT. Discovery latency = first GPS-contiguous segment with two clean 60-s blocks ≥120 s apart (~3 min clean at-speed GPS; urban12's stop-and-go canyon start delays it to ~1300 s). GPS-denial probe: 522 s mask → 0 pairs inside, scores held bit-exact, ranking unchanged.

---

## 1. Design (prototype-proven mechanics → strict core)

### 1.1 Monitor input + anchor gating

The monitor consumes (a) applied/evaluated GPS position fixes in the odom frame (already flowing through the correction path) and (b) each source's aligned yaw increments (already computed for fusion). Course is compared via DELTAS (Δcourse vs Δyaw between anchor pairs) so any fixed frame/mounting yaw offset cancels — no `odom_from_enu_yaw` dependency.

**Anchor gates** (a fix pair becomes a course anchor only if ALL hold; constants fixed in namespace, prototype-swept):
- fix spacing `dt ∈ (0.05, 3] s`; horizontal GPS speed `3 ≤ v ≤ 30 m/s`;
- `|v_gps − v_odo| < max(3 m/s, 0.5·v)` — cross-validate course speed against the sources' own translation (kills multipath jumps);
- net forward motion (reverse flips course 180°);
- `|yaw rate| < 3°/s` over the fix interval (chord-vs-heading transient hits −30..−45° in turns).

**Load-bearing prototype lesson**: anchors are gated to straight-at-speed driving, but residual pairs must BRIDGE turns — wheel-yaw drift is turn-correlated; gating turns out of pair interiors hides exactly the junk source (falsified three prototype designs on urban17). Therefore: **never hard-drop a pair after anchoring** — every drop breaks the telescoping of the cumulative residual and injects a permanent offset.

### 1.2 Residual + consensus split

Per anchor pair, per source `i`: `r_i = wrap(Δcourse − Δyaw_i)`. Common-mode course error is removed by the cross-source median: `m = median_i(r_i)`; the per-source deviation `dev_i = r_i − m` splits into a **rate channel** `m + clamp(dev_i, ±5°)` accumulated cumulatively, and the clamp excess into a per-source **event channel** (catches slip steps, e.g. urban17's 89.5° wheel slip, which any magnitude gate would reject). Requires **n ≥ 3 sources** (median consensus degenerates below; monitor inert otherwise).

### 1.3 Drift estimator (slope, not endpoint sum)

The cumulative residual is a staircase with steps at canyon entries/denials → cross-segment slopes are poison. Estimator: within GPS-CONTIGUOUS segments only (a gap > `pairgap` 60 s closes the segment), form 60-s block medians; drift slope = baseline-weighted median of pairwise block-median slopes with baseline ≥ 120 s. Score_i = |slope_i| + event_rate_i (deg/h). All state online + fixed-capacity: running accumulators, one open block, a bounded slope pool per source (fixed-size reservoir; document the bound).

### 1.4 Weight rule + estimator wiring

`monitor_boost_i = clip(boost_max · min_j score_j / score_i, 1, boost_max)`, `boost_max = 10` (knob). On KAIST this reproduces the hand recipe with no config: fog→10.0, wheel→1.0, wheelimu→1.0–1.1.

- Wiring: `w_rot_i = clamp(w_base·rel_rot_i) · rot_weight_prior_i · monitor_boost_i` — the boost slots beside the config prior (same position outside the clamp, Slice-19 contract). Boost is 1.0 for every source until the monitor has ≥ 2 sources scored (abstain-don't-guess: insufficient data ⇒ all boosts 1.0).
- Boost updates are SLOW (once per completed block at most) + hysteresis on rank changes (a rank flip must persist ≥ 2 consecutive evaluations before boosts move) — fusion weights must not chatter.
- `Config::heading_monitor` (bool, default **false**) — OFF byte-identical (coupled AND split paths). Requires `split_median=true` to have any effect (boost feeds the rotation channel only); flag-on with coupled median = validate() error (no silent no-op).
- `Config::heading_monitor_boost_max` (double, default 10, ≥1). Both knobs in the config-hash. Monitor state NOT persisted (re-warms; same stance as 19b per-channel reliability).
- Diagnostics: per-source score + boost + scored flag in `SourceHealth` (or CalibSnapshot — implementer picks the natural home, reports back).

### 1.5 Honest limits (documented, not solved here)

- Absolute floor ~5–15 deg/h: multiple sources all better than the floor are mutually indistinguishable (KAIST has one heading-grade source; the multi-FOG case is unproven) — boosts then spread between 1 and boost_max harmlessly.
- Drift accrued during a GPS denial longer than `pairgap` is permanently invisible (deltas telescope across the gap by design).
- Course needs motion: a rig that never exceeds v_min never scores (boosts stay 1.0 — graceful).

## 2. Acceptance

Unit (TDD):
1. Sim observability-style headline: 3 sources, one with planted yaw-rate drift + GPS fixes → drifter's score ≫ clean sources', boost ranks clean source(s) up; drift removed from the planted source's boost only (ranking, not absolute recovery).
2. Anchor gates: stopped/slow/turn/reverse/multipath-jump fixes produce no anchors (each gate exercised).
3. Telescoping: a turn between anchors is bridged (pair not dropped), and the turn-correlated drifter is still caught.
4. Denial freeze: a fix gap > pairgap closes the segment, scores hold exactly, no cross-segment slope forms.
5. Abstain: <3 sources, or no scored blocks → all boosts exactly 1.0.
6. Chatter guard: a one-evaluation rank flip does not move boosts; a persistent flip does.
7. Default-off byte-identical (split-ON without monitor unchanged vs Slice-19b; coupled untouched); validate() rejects monitor-without-split.
8. Config-hash flips on both knobs; loader keys.
9. Split-ON cov-cal guard + observability self-tests stay green (boost is a weight, the spread→Q path is unchanged structurally).

Gate: `scripts/dev.ps1 -Task test` fully green.

Real-data (orchestrator, post-merge — KAIST urban07/12/17, urban recipe WITHOUT `rot_weight_prior`, monitor ON):
- urban12 full-drive heading rms ≈ the manual-prior result (≈5.4°, FOG floor 4.6°); local rot p50 in the same band (~0.004 rad).
- urban07/17 no regression vs the 19b baseline.
- Reported boosts ≈ {fog≈10, others≈1}; discovery latency consistent with the prototype (~3 min clean GPS; urban12 late start tolerated).

## 3. Status

- [ ] Implemented (TDD, gate green, committed)
- [ ] Reviewed + findings fixed
- [ ] Real-data validation
- [ ] Docs updated (DECISIONS D29 / CONFIG / ISSUES) — orchestrator
