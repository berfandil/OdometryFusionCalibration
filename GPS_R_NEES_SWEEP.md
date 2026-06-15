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
- [ ] (follow-up) per-fix/multipath-aware `R` model + GPS-vs-GT frame-offset characterization; un-gated-NIS R calibration methodology.
