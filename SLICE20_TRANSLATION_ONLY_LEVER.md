# Slice 20 — Translation-only source lever (velocity/Doppler sensors)

**Goal**: let a TRANSLATION-ONLY source (Doppler radar, optical-flow, any sensor that measures its own velocity but not its rotation) calibrate its LEVER ARM. Today such a source gets lever conf 0 (never commits) AND, worse, commits a CONFIDENTLY-WRONG rotation extrinsic. This slice adds a `translation_only` source flag that pins the rotation extrinsic to the prior and runs the existing (already-correct) hand-eye lever LS with that trusted `R_X` — recovering the observable lever axes.

Investigated + prototyped first (the 15-19 pattern). Evidence: `tools/proto_radar_lever.py` on real nuScenes radar.

---

## 1. Why it fails today (investigation, `src/core/calibration.cpp`)

The Phase-2 lever row is already the right equation and already gates on the CONSENSUS turn, not the source's rotation:
- Row: `(R_A - I) t_X = R_X t_B - t_A` (calibration.cpp:1025/1103) — `R_A,t_A` = fused/consensus motion, `t_B` = source translation. Never uses the source's own rotation `R_B`.
- Turn gate: `||fused_omega|| > turn_omega_min` (calibration.cpp:895) + per-window `||log R_A|| >= kRotRowMin` (calibration.cpp:928) — BOTH on the base, not the source.

The single break point is **`R_X`** (calibration.cpp:1030): `R_X = rot3d_row ? Rhat : R_yp*Rx(roll)`.
- A heading-blind source (`R_B = I`) -> rot3d's `BBw = sum w*log(R_B)log(R_B)^T = 0` -> the two-axis gate never opens -> `rot3d_row = false`.
- So `R_X = R_yp * Rx(roll)`, where `R_yp` is Phase-1's yaw/pitch. Phase-1 tries to fit yaw/pitch from the source's (absent / noise) rotation and **commits garbage** (prototype/probe: a front radar committed `rx = -0.984 rad` at conf 0.95). A polluted `R_X` -> polluted lever rows -> lever conf 0 (and a wrong rotation extrinsic fed back to fusion = a real footgun).

So the fix is NOT a new estimator: it is to STOP corrupting `R_X` for a declared translation-only source — keep `R_X = prior` and let the existing lever LS run.

## 2. Prototype evidence (`tools/proto_radar_lever.py`, real nuScenes radar)

The linear estimator `(R_A-I)t_X = R_X t_B - t_A` accumulated over turn windows (GT-base or wheel/imu-consensus-base, true `R_X`):
- Synthetic control (planted `t_B`): all 5 levers recovered to ~1e-13 cm -> math + code exact.
- Real radar: **lateral lever recovers to ~5-25 cm** (back_right 4 cm, front_right 11 cm, front 19 cm, back_left 25 cm) on the clean radars.
- Along-track lever: **~2 m error — a sensor SNR wall** (per-window `t_B` carries ~0.5 m along-track Doppler drift; the turn-induced lever signal is ~0.65 m; SNR ~1.26). Accumulation does NOT shrink it (systematic, not averageable). NOT a math/data-quantity gap.
- Falsification: wrong `R_X` (+20 deg) -> 113-179 cm (R_X-correctness is load-bearing); z eigenvalue ~120x smaller than x,y (planar null) -> pinned at prior; a zero-lever source recovers ~0.

## 3. Design

### 3.1 Config
- `SensorConfig::translation_only` (bool, default **false**). Semantics: this source's per-step rotation increment is uninformative (`R_B ~ I`); calibrate its lever + scale from translation, take its rotation extrinsic AS GIVEN (the prior). Loader key (per-sensor) `translation_only`. Joins the persistence config-hash. `validate()`: a `translation_only` source must still have a valid `prior_extrinsic` (its rotation is trusted, not recovered).

### 3.2 Calibrator behaviour when `translation_only` (per source)
- **Rotation extrinsic: PINNED to the prior.** Skip Phase-1 yaw/pitch voting, skip rot3d accumulation/vote, skip the per-window roll recovery FOR THIS SOURCE. `R_yp` stays = the prior yaw/pitch; `roll = 0`. So the lever row uses `R_X = prior_rotation` (correct). The rotation-extrinsic snapshot is the prior, `extrinsic_committed = false` (honestly: not recovered, trusted). This ALSO removes the garbage-rotation-commit footgun.
- **Lever LS: unchanged**, now fed the clean `R_X = prior`. Existing turn gate (consensus), existing obs gate (`AtA(c,c) >= kAxisInfoMin*max_diag`) for the planar-null (z) axis.
- **Scale: KEEP** (straight-regime magnitude ratio is valid for a velocity source). The joint lever+scale path (17b) may run; if it gates on rot3d (which is off here) confirm it degrades to the 3-unknown lever (R_X from prior) — do not require rot3d for a translation_only source's lever.
- **Never commit a wrong extrinsic — residual-gated lever confidence (the along-track guard).** The along-track axis is well-CONDITIONED but sensor-noise-BIASED, so the conditioning obs gate alone would commit a ~2 m-wrong value. Add a per-source normalized LS residual `rho = ||M t_X_hat - rhs|| / ||rhs||` (running, decay-aged) and scale the lever vote weight / confidence down as `rho` rises (a high-residual fit = the source's `t_B` is too noisy for that geometry). Pin the threshold from the prototype so the lateral lever (low residual) commits and the along-track (high residual) is WITHHELD (conf below the commit gate). Constant, namespaced; document. (If a single scalar residual cannot separate lateral-good from along-track-bad, fall back to per-axis residual attribution — decide during impl from the prototype numbers.)

### 3.3 What stays untouched
- Sources without the flag: byte-identical (the flag default false; all new code behind `if (translation_only)`).
- Fusion, the split median, the heading monitor, persistence format (hash gains the flag -> pre-20 blobs cold-start by design).
- The observability self-tests (each DOF converges in its regime) — add a translation_only case, do not weaken existing ones.

## 4. Acceptance

Unit (TDD):
1. **Headline (sim)**: a planted translation-only source (`R_B = I`, true lever, base turns) -> the existing-calibrator path gives lever conf 0 / garbage rotation; with `translation_only=true` the OBSERVABLE lever axes recover (to the sim tolerance) and the rotation extrinsic stays = prior (not garbage). Pin both the recovery AND the no-garbage-rotation.
2. **Footgun fix**: a heading-blind source with `translation_only=true` never commits a rotation extrinsic off its prior (the probe's `rx=-0.984` case must not occur).
3. **Residual gate**: a source whose `t_B` is clean on the lateral axis but noisy along-track commits ONLY the lateral lever (along-track withheld via the residual gate) — never commits the wrong along-track value.
4. **z / planar null**: yaw-only turns leave z unobservable -> withheld (existing obs gate), x,y observable.
5. **Default-off byte-identical**: `translation_only=false` reproduces current behaviour exactly (golden + existing tests untouched; exact-equality pin).
6. **Scale preserved**: a translation_only source still recovers its straight-regime scale.
7. **Loader key + config-hash flip**; `validate()` rejects nonsense.
8. Observability self-tests (incl. a new translation_only case) green; cov-cal band + golden untouched.

Gate: `scripts/dev.ps1 -Task test` fully green.

Real-data (orchestrator, post-merge):
- nuScenes mini concat (5 radars `translation_only=true`, true mounts as prior, lever perturbed +0.20 m): the LATERAL lever recovers toward truth (target ~the prototype's 5-25 cm), the along-track lever is WITHHELD (not committed wrong), no garbage radar rotation commits. Compare to the pre-flag run (all conf 0 / garbage rotation).
- Sources without the flag on the same run unchanged.

## 4b. Slice 20b — per-axis lever commit (the radar unlock)

**Problem**: the lever commits as a whole 3-vector (`lever_committed[i]` gated on `translation_confidence` = min over axes). For a translation_only radar the clean LATERAL axis can never commit because the observed-but-NOISY along-track axis drags the min below the gate (min-over-observed got conf to 0.19-0.55, still < 0.6). The axes are independently observable — commit them independently.

**Design (per-axis commit, translation_only sources ONLY — non-flagged stay whole-lever, byte-identical)**:
- **Calibrator accessors** (per-axis values already exist in `xyz_[3*s+c]`): add `translation_confidence_axis(id, c)` = `xyz_[3*s+c].confidence()` and `xyz_axis_mass(id, c)` = `xyz_[3*s+c].total()`. `lever_arm(id)` already returns the per-axis 3-vector.
- **Estimator**: `bool lever_committed_xyz[kMaxSourcesCap][3]`. For a translation_only source, replace the whole-lever gate with three per-axis `commit_gate_reanchor(prev_xyz[c], conf_axis(c), mass_axis(c), commit_concentration, commit_min_votes)` (same hysteresis, per axis). On a committed axis c, publish `prior_extrinsic[i].t[c] = lever_arm(id)[c]` (the other axes stay at prior). `lever_committed[i] = OR(lever_committed_xyz[i][*])` (drives `CalibSnapshot.translation_committed` + the existing median-join / serialize). Non-flagged sources: the existing whole-lever path UNCHANGED.
- **Diagnostics**: add `CalibSnapshot.translation_committed_xyz[3]` (per-axis flags; whole `translation_committed` = OR). The replay summary can show which axes committed.
- **Persistence**: serialize the 3 per-axis flags (format bump v4; pre-v4 blobs cold-start). For a non-flagged source the 3 mirror the whole flag.
- **Outcome (radar)**: lateral axis (clean, `c >= 0.6`) COMMITS + folds into the prior; along-track (noisy, `c < 0.6`) + null z (empty) held at prior — never a wrong commit.

**Acceptance (added to the 8)**:
9. Per-axis headline (sim): a translation_only source clean on the lateral axis, noisy on along-track, null on z -> ONLY the lateral axis commits (folded into prior); along-track + z held at prior; per-axis flags correct; the whole `translation_committed` = OR.
10. Non-flagged byte-identical (whole-lever path untouched; exact-equality pin).
11. Persistence round-trips the per-axis flags (format v4; pre-v4 cold-start).
12. Real-data: nuScenes radar LATERAL lever commits (if its concentration clears the gate) with along-track/z held — or an honest report if per-window radar noise keeps even the lateral below the gate.

## 5. Status
- [x] Investigated (calibration.cpp lever path) + prototyped (`tools/proto_radar_lever.py`, `8df5de8`) — lateral lever recovers, along-track is a sensor SNR wall, `R_X`-from-prior is load-bearing.
- [x] Implemented (TDD, gate green, committed `629d379`): `translation_only` flag pins rotation to prior (skip Phase-1 direction vote / rot3d / roll), lever LS runs with `R_X = prior`, scale kept. Along-track guard = the existing concentration gate withholds (no residual gate added). 8 unit cases.
- [x] Reviewed (`reviews/slice-20-review.md`: APPROVE, clean) — byte-identical default-off verified incl. the scale-vote-drop equivalence the implementer self-caught.
- [x] Real-data validation (nuScenes mini concat, 5 radars `translation_only=true`, lever perturbed +0.20 m) — see below.
- [x] Docs updated (DECISIONS D30 / CONFIG / ISSUES) — orchestrator

### Real-data result (nuScenes radar, 2026-06-15) — footgun FIXED, committed lever NOT delivered (honest)

- **Footgun fix CONFIRMED on real radar**: every radar's `extr_rot = [0, 0, prior_yaw]`, conf 0 — exactly the prior, no spurious roll/pitch, nothing committed off prior. The pre-flag garbage (`rx -0.984` at conf 0.95) is GONE. Rotation is trusted-at-prior as designed. This is the primary, verified win + a real safety fix (a wrong rotation extrinsic was being fed back to fusion).
- **The radar LEVER does NOT commit on real data** (conf 0 on all 5), and the readout OVERSHOOTS the lateral axis (src3 y: perturbed 1.0 -> readout 0.516 vs true 0.8; src6 y: -0.418 -> -0.710 vs -0.618) — noisy, NOT the ~cm clean recovery the prototype showed. Two causes: (1) `translation_confidence = min(cx,cy,cz)` and z is unobservable on planar -> the min is pinned at 0 -> the lever can never commit on planar driving (a structural min-over-unobservable-axis issue, not specific to this flag); (2) the PRODUCTION per-window / consensus-base estimator is materially noisier than the prototype's idealized GT-base BATCH LS -> the radar along-track Doppler noise dominates the per-window lever vote, so even the lateral axis does not concentrate to the commit gate.
- **Lesson (the project's recurring one)**: the `proto_radar_lever.py` clean lateral recovery used GT-base + a single batch LS over all windows; the real estimator uses the consensus base + per-window histogram votes + the concentration commit gate. Do NOT claim a real-data win from an idealized prototype — the full production path revealed the gap. The flag + footgun fix are solid; committed lever recovery on real radar is NOT achieved.

### Follow-ups (for committed radar lever recovery)
1. **`min`-over-OBSERVED-axes lever confidence** — DONE (`kAxisObservedFrac=0.1`; for `translation_only` sources `translation_confidence` is the min over axes whose histogram reached >=10% of the most-voted axis's mass, excluding the empty planar-null z). **Result (nuScenes radar re-validation)**: lever confidence rose from **0 -> 0.19-0.55** (the z-min pin is gone; sim clean planar now COMMITS — pinned in `test_calib_translation_only.cpp`). But **still no commit on real radar** (all < the 0.6 gate; src2 closest 0.553): the observed-but-NOISY along-track x drags the joint `min(cx,cy)` below the gate. Necessary fix, NOT sufficient — as predicted.
2. **Per-axis lever commit** (the real unlock): commit each lever axis when ITS concentration clears the gate, independent of the others, so the clean LATERAL axis commits while the noisy along-track + the null z are each held at the prior. A larger change (the commit + re-anchor feedback are whole-lever today). This is what would actually deliver a committed radar lateral lever.
3. **Lower per-window radar noise** — better Doppler ego-velocity (tighter RANSAC / scale handling), longer windows, or a batch lever-only solve at the deeper calibration frontier instead of per-window histogram votes. The along-track axis remains a sensor SNR wall regardless (prototype-confirmed).
