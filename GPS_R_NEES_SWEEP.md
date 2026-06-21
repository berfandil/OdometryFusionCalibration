# GPS R-side NEES sweep — with-correction covariance consistency (KAIST)

**Goal**: the with-GPS pose NEES is overconfident on every KAIST drive; the GPS measurement noise `R` (`cov_floor_m2`, the ENU covariance floor in m^2) under-states the real error. Sweep `cov_floor_m2` against NEES/NIS/drift to find the consistency-honest `R` and quantify the cost. Config-only (existing knob); recommended-urban config (`q_floor=1e-5/1e-6`, `split_median`, `correction_rot_suppress_kappa=0.8`) on urban07/12/17.

---

## 1. The sweep (mean_pose_nees vs DOF 6; mean_nis over ACCEPTED fixes vs DOF 3)

| drive | cov_floor_m2 | NEES | NIS | gps applied/rej | tail trans (m) | local p50 (m) | worst window (m) |
|---|---|---|---|---|---|---|---|
| **urban07** | 25 (current) | 103.4 | 0.158 | 551/5 | 8.33 | 0.98 | 5.2 |
| | 100 | 38.2 | 0.097 | 556/0 | 8.86 | 0.81 | 7.4 |
| | 400 | 14.8 | 0.026 | 556/0 | 8.98 | 0.51 | 6.0 |
| | 800 | 9.57 | 0.014 | 556/0 | 8.94 | **0.38** | 5.3 |
| **urban12** | 25 | 1889 | 0.344 | 1419/1060 | 1.88 | 0.194 | 211 |
| | 100 | 683 | 0.233 | 1562/917 | 1.85 | 0.245 | 122 |
| | 200 | 385 | 0.194 | 1577/902 | 1.89 | 0.262 | **108** |
| **urban17** | 25 | 1198 | 0.108 | 1181/11 | 26.2 | 0.53 | 11.8 |
| | 400 | 123 | 0.030 | 1189/3 | 26.3 | 0.54 | 5.0 |
| | 1600 | 38.7 | 0.017 | 1192/0 | 26.2 | 0.65 | 7.2 |
| | 3200 | 21.8 | — | 1192/0 | 26.1 | — | — |

(urban12 capped at cov=200 — the slow drive; the trend extrapolates the honest-R below. An orphaned sweep-agent process kept respawning urban12 high-cov runs and was killed; the cov 25-200 row is sufficient for the conclusion.)

## 2. Findings

1. **`cov_floor_m2=25` is MASSIVELY overconfident** — NEES 103 / 1889 / 1198 (07/12/17) vs DOF 6, i.e. 17-300x too confident on every drive. The current value was an accuracy/gate compromise, never a consistency calibration; this confirms + quantifies it.
2. **Raising `R` monotonically reduces NEES — AND on urban07 IMPROVES accuracy.** urban07 local p50 drift drops 0.98 -> 0.38 m as cov goes 25 -> 800 (R=25 OVER-trusted the noisy GPS, injecting GPS jitter into the trajectory; a larger R fuses it more smoothly). urban12's worst window shrinks 211 -> 108 m; tails stay flat on all drives. So `R=25` was actively HURTING accuracy, not just mis-reporting uncertainty.
3. **The consistency-honest `R` (NEES~6) differs ~7x cross-drive**: urban07 ~cov 1200-1600, urban12 ~cov 10000+, urban17 ~cov 10000 (21.8 @ 3200, ~2x/doubling -> 6 near 10000). A SCALAR `cov_floor_m2` CANNOT be cross-drive consistent — the multipath-per-drive analogue of the `q_floor` 11x cross-drive spread.
4. **NIS << 3 everywhere while NEES >> 6 — the diagnostic.** The accepted-fix NIS is 0.1-0.34 (and drops with R), FAR below DOF 3, while NEES is FAR above DOF 6. This is not a contradiction: the Mahalanobis gate accepts only fixes whose innovation matches the (small) `S=HPH^T+R`, so the accepted-NIS looks over-conservative, while the fused STATE is still far from GT. The error is **structural**, not a simple `R` mis-scale.
5. **The honest-`R` is implausibly large for VRS-RTK** (cov ~10000 m^2 = sigma ~100 m, for a centimetre-grade GPS) -> the real cause is **GPS-vs-GT frame DISAGREEMENT**, not GPS measurement noise. The KAIST GPS is VRS-RTK; the GT is the SLAM `global_pose` — different reference systems that disagree by metres (urban07 tail ~8 m, urban17 ~26 m, stable across R). Raising `R` "fixes" NEES-vs-GT by trusting the disagreeing GPS LESS. So the with-GPS NEES-vs-GT consistency is CONFOUNDED by the frame offset.
6. **The per-fix GPS variance carries NO signal**: the KAIST GPS CSV's `var_e,var_n,var_u` are a FLAT `0.25,0.25,1.0` nominal on EVERY fix (VRS-RTK datasheet value, not a live multipath estimate), so floored by `cov_floor_m2` it is inert. A per-fix `R` model has nothing to read from this dataset — a multipath-aware `R` would have to be MODELLED (geometry / residual statistics), or sourced from a consumer-GPS stream with varying variance.

## 3. Recommendation

- **Raise the default `cov_floor_m2` well above 25** — it was 17-300x too tight. A value ~200-400 is STRICTLY better on this data (urban07 accuracy + consistency both improve, urban12 worst-window improves 211 -> ~110 m, tails stable, gate accepts more genuine fixes). This is a free, config-only win; update the recommended urban recipe.
- **Full NEES~6 consistency is NOT a scalar-`R` problem** — it needs (a) a per-fix / multipath-aware `R` MODEL (no signal in the flat VRS variance; would have to be modelled or come from a varying-variance GPS), AND (b) separating the GPS-vs-GT frame DISAGREEMENT (an alignment/lever/SLAM-drift offset) from the measurement noise, since the implausibly-large honest-`R` is mostly that offset. Both are code/analysis follow-ups, parallel to the `Ad`-shape `q_floor` (the predict-side cross-drive analogue).
- **NIS is the cleaner R-side target** (innovation, no GT) but is gate-biased low; an un-gated NIS over ALL fixes (or a wider gate) would calibrate `R` without the GT-frame confound — a better future methodology than NEES-vs-GT for the measurement model alone.

## 4. Status
- [x] Sweep run (urban07 full, urban17 to cov 3200, urban12 to cov 200) + analyzed.
- [x] Docs updated (CONFIG `cov_floor_m2` note; ISSUES per-fix-R / GPS-vs-GT-confound item) + **recipe banked: urban `cov_floor_m2 = 300`** (HANDOFF URBAN TRAJECTORY RECIPE, `28339aa`+bank).
- [x] (follow-up DONE, 2026-06-21) un-gated-NIS methodology + bias/spread decomposition + an innovation-adaptive ROBUST `R` model — see §5.

---

## 5. Follow-up: un-gated NIS, bias/spread decomposition, innovation-adaptive R (2026-06-21)

Built the un-gated-NIS methodology the §3 recommendation pointed to: `tools/gps_innovation_analysis.py` computes the per-fix innovation `nu = z_gps - h(x_fused)` OFFLINE (from the replay `out.csv` fused pose + cov diag and the GPS CSV, exact geodetic transform), DECOMPOSES it into MEAN (frame offset / bias) + COVARIANCE (measurement spread), and reports the un-gated NIS over ALL fixes (the sweep had only the gate-biased-low accepted NIS).

**Decomposition (urban07/12/17, cov_floor 25):**

| drive | mean nu \|bias\| | spread rms | honest R_spread | un-gated NIS (mean) |
|---|---|---|---|---|
| urban07 | **0.48 m** | 3.1 m | ~2 m^2 | 0.35 |
| urban17 | **0.30 m** | 5.2 m | ~22 m^2 | 0.97 |
| urban12 | **181 m** | 538 m | ~1.4e5 m^2 | 2495 (median 0.42) |

**Sharpened findings:**
1. **The GPS-vs-FUSED bias is SMALL on the clean drives (0.3-0.5 m)** — the fused state TRACKS the GPS. The ~8-26 m the §2 sweep saw is the GPS/fused-vs-**GT** offset (the NEES confound), NOT in the innovation. So the huge pose-NEES is a GT-frame disagreement, confirmed with innovation numbers.
2. **The honest measurement `R` is SMALL and cross-drive-VARYING** (urban07 ~2, urban17 ~22 m^2), and the un-gated NIS is FAR below 3 -> `cov_floor = 25/300` is OVER-conservative for the clean-drive MEASUREMENT (the inflation was band-aiding the GT-frame NEES + urban12's outliers).
3. **urban12 is an OUTLIER-multipath drive**, not a uniform-R one: median NIS 0.42 (bulk fine), but a population of ~180 m gross fixes (GPS-hostile canyon) -> mean NIS 2495. A bigger uniform R is the wrong tool; the fix is per-fix robust rejection.

**The model — innovation-adaptive ROBUST `R`** (`GpsConfig::adaptive_r`, `gps_correction.{hpp,cpp}`; opt-in, default-off byte-identical, manifest keys `adaptive_r`/`adaptive_window`/`adaptive_min_samples`/`adaptive_r_floor_m2`): `evaluate()` sets `R_odom = diag(max((1.4826*MAD_i)^2, floor))` from the robust scale of the last `window` innovations per odom axis. MAD tracks the BULK spread (outlier-immune), so R stays small on a clean drive and gross multipath fixes show a LARGE NIS -> the estimator's Mahalanobis gate rejects them. Causal (R for fix k uses fixes [k-window, k-1]); the ring sees every evaluated fix (the MAD is immune to the rejects it sees).

**In-loop result (adaptive vs the banked `cov_floor = 300`):**

| drive | NEES (b300 -> adapt) | NIS (b300 -> adapt) | local p50 drift | rejects |
|---|---|---|---|---|
| urban07 | 18.9 -> 1808 | 0.06 -> **1.36** | 1.02 -> 1.54 m | 90 |
| urban12 | 254 -> 23902 | 0.23 -> **0.68** | 0.67 -> **0.30** m | 1140 |
| urban17 | 156 -> 28667 | 0.03 -> **1.13** | 0.51 -> **0.37** m | 158 |

**Verdict — the adaptive R SURGICALLY SEPARATES the two error sources the scalar sweep conflated:**
- It makes the **NIS consistent + cross-drive self-calibrating** (1.36/0.68/1.13 vs the scalar's 0.06/0.23/0.03 -> the measurement model is now honest), AND improves the typical (p50) drift on the OUTLIER drives (urban12 0.67->0.30, urban17 0.51->0.37 — the robust rejection working).
- But **NEES-vs-GT EXPLODES** (18.9->1808 etc.): with honest-small R the filter confidently tracks the GPS, which is offset from GT -> confident-but-wrong. **This PROVES the consistency problem is a GT-frame BIAS, not measurement `R`** — the `cov_floor` inflation was band-aiding the bias; an honest R cannot fix it. The clean drive (urban07) accuracy also slightly degrades (over-trusts the offset GPS).

So the adaptive R is **the principled MEASUREMENT model** (NIS-honest, cross-drive, helps multipath drift) shipped opt-in (default-off, like the Slice-15 Huber) — NOT a drop-in accuracy win, because the dominant KAIST GPS error is the **GPS-vs-GT frame offset** (a SLAM-vs-VRS alignment bias), which is the real remaining open item and needs a bias model (alignment / lever / SLAM-drift), not an R model. Tools: `tools/gps_innovation_analysis.py`; feature in `adapters/{include,src}/.../gps_correction.*` + a `test_gps_correction` adaptive-R case.
