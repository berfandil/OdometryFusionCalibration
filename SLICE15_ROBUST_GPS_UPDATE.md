# Slice 15 — Robust correction update (large-innovation safety)

Status: `[ ]` DESIGN (not implemented). Design-first; no code until accepted.
Deps: Slice 13 (replay harness + local metric), Slice 14 (covariance calibration).

## 1. Problem

On the KAIST `urban12` drive (41 min, dense urban, GPS multipath) the fused estimate
diverges to **4214 m** global tail / pose-NEES 60k. The new GT-anchored local metric
(`7c71b63`) deconstructs that divergence:

| run    | global tail | local p50 | local mean | local max |
|--------|-------------|-----------|------------|-----------|
| urban07| 7.9 m       | 4.68 m    | 5.51 m     | 31.9 m    |
| urban17| 24.8 m      | 0.23 m    | 1.35 m     | 18.2 m    |
| urban12| **4214 m**  | **0.17 m**| 9.54 m     | **856 m** |

`urban12`'s **median 10 s window is 0.17 m — the best of the three.** The km-scale
divergence is NOT pervasive: mean/p50 = 56× skew, and `max = 856 m` in a single 10 s
window. The odometry+fusion is locally excellent; a **handful of catastrophic GPS-correction
events** corrupt the global frame, and the continuous filter then carries that corruption
forever. **Fix target is narrow: bound the damage a single large-innovation correction can do.**

## 2. Mechanism (grounded in `src/core/eskf.cpp` `update()`)

A GPS fix is a dim-3 position measurement. `update()` (eskf.cpp:139):

```
S  = H P H^T + R                 // innovation cov                       (eskf.cpp:161)
d2 = r^T S^-1 r                  // NIS (Mahalanobis), gate on this      (eskf.cpp:166)
K  = P H^T S^-1                  // 12x6 gain                            (eskf.cpp:172-176)
dx = K r                         // 12-vector error correction           (eskf.cpp:179)
pose <- pose o exp(dx[0..5])     // rotation kick = exp(dx[3..5])         (eskf.cpp:187)
```

**Why heading is corrupted even with lever = 0:** with `lever = 0`, `H`'s rotation columns
are zero (the lever-arm coupling `−R[lever]×` vanishes — that was the `7e98dcc` fix). But
`K`'s **rotation rows are still nonzero**, fed by the **trans↔rot cross-covariance** block of
`P` (`P[3:6, 0:3]`, built by the SE(3)-Ad predict). So a large position residual `r` injects a
large `dx[3..5]` → `exp(dx[3..5])` rotates the heading. One bad fix (multipath jump, or a fix
applied while the global frame is already offset) → large `r` → large heading kick → the
continuous filter diverges. This is the residual we must bound.

The Mahalanobis gate (eskf.cpp:170) is a *binary* defense and the wrong tool here: with the
overconfident `P`, `S ≈ R` is small, so `urban12` already rejects 2383/2479 fixes — and the few
it accepts are the blowups. We need a **graded** defense, not a harder binary gate.

## 3. Design options

| # | lever | mechanism | principled | scope |
|---|-------|-----------|------------|-------|
| A | Huber / M-estimator | down-weight gain smoothly as NIS grows (≡ inflate R for large `r`) | ✅ textbook | core `update()` |
| B | cap injected `dx[3..5]` | hard-clamp rotation-correction magnitude per update | ❌ ad-hoc, biases | core `update()` |
| C | decouple pos-fix from rot | zero/shrink `K`'s rotation rows (or `P` trans-rot cross-cov) for a pure position fix | ⚠️ discards real info | core `update()` |
| D | realistic `R` | raise GPS `cov_floor_m2` so `S` reflects true urban GPS error | ✅ root cause | config only (exists) |

**A vs C:** C throws away that position fixes legitimately inform heading over time
(along/cross-track). A keeps that signal in the normal regime and only attenuates the *outlier*
fixes. A wins. **B** is a blunt safety clamp — keep it only as an optional belt-and-suspenders
behind its own flag. **D** is the documented root cause (Slice 14: `cov_floor_m2` 0→100 dropped
NEES 687→87) and is *complementary* to A, not redundant: D fixes the covariance *magnitude*
(NEES, and unblocks the gate so good fixes are accepted); A bounds the *residual-spike injection*
(the 856 m windows). Ship A+D together.

## 4. Chosen approach — Huber-robust gain (A) + realistic R (D)

**Robust weight (IGG/Huber down-weighting), in the generic `update()` so every `ICorrection`
benefits:**

```
d2  = r^T S^-1 r                          // existing NIS
dbar = sqrt(d2 / n)                        // RMS Mahalanobis per active DOF
w   = (dbar <= kappa) ? 1 : kappa / dbar   // Huber weight in (0,1]
R'  = R / w                                // inflate noise for outliers  (w<1 -> R'>R)
S'  = H P H^T + R' ;  K' = P H^T S'^-1 ;  dx = K' r
```

- `w = 1` in the normal regime → **bit-identical to today** (no regression on clean data).
- For an outlier (`dbar > kappa`) the gain shrinks ∝ `w`, smoothly bounding `dx` (and thus the
  heading kick) instead of a binary accept/reject. As `dbar → ∞`, influence → 0.
- Inflating `R` (not just scaling `dx`) keeps the **Joseph-form covariance consistent** with the
  gain actually applied (`state_.cov` update at eskf.cpp:194-195 uses `R'`, `K'`). One extra 6×6
  LDLT solve only when `w < 1`.
- Apply the SAME path in `update_aug()` (eskf.cpp:305, the bias-augmented update).

**Config:** add `Config::correction_robust_kappa` (Huber threshold in per-DOF σ units).
`0` (default) = **disabled** → exact current behavior (opt-in, like the gate / `cov_floor`).
Expose via `ConfigLoader` `[global]` + manifest. Optional separate `correction_max_rot_rad`
(lever B clamp; default 0 = off).

`D` needs no code — set `cov_floor_m2` to a realistic urban value in the KAIST manifests.

## 5. Acceptance test (the local metric makes this objective)

Pass criteria, all three KAIST drives, full length (no truncated slices — the metric would
spike `max` on a hidden blowup, so a slice can't lie this time):

1. **urban12** `local max` 856 m → **≤ ~50 m** (toward its p50 0.17 m); global tail bounded to
   tens of m; `gps_applied` recovers well above 96/2479 (robust + realistic R accepts good fixes
   AND attenuates bad ones).
2. **urban07 / urban17** `local p50` and `mean` **not regressed** (±10 %).
3. **Sim** `test_cov_calibration` NEES still ≈ 4.82 / never overconfident; `test_gps`/e2e GPS still
   removes drift and rejects the planted outlier. (`kappa = 0` path must be bit-identical → these
   are unchanged unless the sim configs opt in.)
4. New unit test: a planted huge GPS innovation → with `kappa` set, injected `‖dx[3..5]‖` is
   bounded and shrinks monotonically as the residual grows; with `kappa = 0` the old (large) kick
   is reproduced.

## 6. Risks / rollback

- **Strict-core edit** (the most sensitive path). Mitigate: `kappa = 0` is exact-current, so the
  change is provably inert until opted in; Joseph form keeps PSD with the scaled gain.
- **Masking vs fixing:** A could paper over a still-wrong `P`. Mitigate: ship D (realistic R) so
  NEES is addressed at the root, not just the spike.
- **Tuning `kappa`:** too low → rejects legitimate corrections (drift); too high → no protection.
  Start `kappa ≈ 3` (3σ per DOF), sweep against the acceptance metric.
- **Slow iteration:** full `urban12` ≈ minutes/run. Accept; the metric gives a clear signal.
- Rollback: set `kappa = 0` (one config line) — no revert needed.

## 7. Out of scope

Distance-windowed metric (discarded). Source-fusion robustness (the median already handles
gross-outlier *sources* — D3 fix). NHC / `Ad`-shape predict-Q (Slice 14 second-order). The
`urban12` post-event recovery dynamics beyond bounding the injection.

## 8. OUTCOME (urban12 sweep, full drives)

| variant | tail m | local max m | NEES | accepted-NIS | gps_applied |
|---------|--------|-------------|------|--------------|-------------|
| baseline (kappa 0, floor 0, gate 9) | 4214 | 856 | 59907 | 1.49 | 96 |
| **D: floor 25, gate 9** (realistic R) | **3075** | **476** | **225** | **0.88** | 188 |
| A: kappa 3, floor 0, gate 9 | 4214 | 856 | 59907 | — | 96 (inert) |
| open gate (kappa 0, floor 25, gate 100) | 4605 | 1951 | 66747 | 38.2 | 298 |
| open gate + Huber (kappa 2, floor 25, gate 100) | 3055 | 2347 | 7534 | 20.8 | 171 |

**D (realistic GPS `R`) works and is the keeper** — config-only, big consistency win (NEES 59907→225),
local max halved, accepted fixes consistent (NIS 0.88). **A (Huber) is the wrong lever for urban12**:
under the tight gate its trigger band (`d2 > kappa²·n`) lies *inside* the gate reject band (`d2 > gate`)
→ inert; opening the gate to let it fire admits genuine multipath (admitted-NIS 38 ≫ 3) → net WORSE.
**The Mahalanobis gate is correctly rejecting the bad fixes** — they are not good fixes wrongly dropped.
A is KEPT as an opt-in robustness primitive (loose-gate / non-GPS corrections, helped NEES 66k→7.5k when
the gate was open), NOT as the urban12 fix. Lever B (rot clamp) NOT built — A's outcome made it moot.

**urban12 remains unsolved** and is reframed as a *localize-the-bad-window* problem (which event makes the
one 476 m window), not a robust-update gain problem. Lesson: Huber/gain levers were the 4th thing to
help-but-not-fix urban12 — stop adding gain knobs; instrument the bad window instead.

## 9. Open questions

- `kappa` per-DOF RMS (`sqrt(d2/n)`) vs raw `sqrt(d2)` threshold — RMS makes `kappa` dimension-agnostic; confirm.
- Should B (rot clamp) ship at all, or is A alone sufficient on `urban12`? Decide empirically after A.
- Realistic `cov_floor_m2` value for KAIST VRS-RTK vs consumer (per-source, from the data).
