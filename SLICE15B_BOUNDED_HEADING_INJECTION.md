# Slice 15b — Bounded heading injection from position-only corrections (lever C)

Status: `[ ]` DESIGN (not implemented). Design-first; no code until accepted.
Deps: Slice 13 (local metric), Slice 15 (robust update + the urban12 diagnosis below).

## 1. Problem — pinned to one timestamp

`tools/diagnose_window.py` on the keeper config (`urban12_k0f25`: realistic R, tight gate) localized
the urban12 divergence to a single applied GPS fix at **t = 1929.2 s**:

```
t=1928.8s  trans_err=648.93m  rot_err=0.975rad
t=1929.2s  trans_err=597.53m  rot_err=2.162rad   <-- APPLIED, NIS=8.45
```

In ONE step an applied dim-3 *position* fix pulls translation 51 m closer but **rotates heading
+1.19 rad (68 deg)**. Heading is now 2.16 rad wrong -> dead-reckoning drives the wrong way ->
every later fix is a genuine ~800 m outlier (NIS 780+) -> correctly gated out -> 3 km runaway.

Upstream trigger: a **522 s GPS-denied coast** (last applied fix at t=1398.9 s) let DR drift
position to ~650 m and heading to ~1 rad; the returning fix then delivered the fatal kick. Earlier
windows (t=574, 584, 1396, each applied=1) show the SAME mechanism in miniature — position fixes
intermittently kick heading, usually recovered; the post-coast one was fatal.

## 2. Mechanism

A position fix has `H = [R | 0 | 0]` — its rotation columns are zero (lever=0). Yet in `update()`:

```
K  = P H^T S^-1                  // 12x6 gain
dx = K r                         // dx[3:6] = rotation correction
pose <- pose o exp(dx[0:6])      // exp(dx[3:6]) rotates heading
```

`K`'s rotation rows (3..5) are nonzero because of the **P trans-rot cross-covariance** `P[3:6,0:3]`
(built by the SE(3)-Ad predict). So a position residual injects heading via cross-cov. This is fine
under small-residual linearization, but a **650 m residual is far outside linear validity** — the
"position informs heading" coupling becomes garbage and over-rotates.

The damaging fix has **NIS 8.45 -> dbar = sqrt(8.45/3) = 1.68**, *just under* the gate (9, dbar
ceiling 1.73). So:
- The Mahalanobis gate cannot stop it (it is in-gate by design).
- Slice-15 Huber (A) cannot stop it at any practical kappa: kappa<1.68 is needed to fire, and even
  then A scales the WHOLE gain — it would shrink the *useful* 51 m position pull along with the
  heading kick. A is the wrong tool: the translation correction is good; only the rotation is fatal.

## 3. Design options (lever C = suppress the rotation a position fix injects)

| # | option | mechanism | pro | con |
|---|--------|-----------|-----|-----|
| C1 | always zero K rotation rows for a position-only fix | dx[3:6] := 0 always | kills over-rotation entirely | forbids ALL position->heading info (loses slow heading observability) |
| C2 | zero P trans-rot cross-cov in the gain | P'[3:6,0:3]=0 for this update | same effect, cov-flavored | gain/Joseph consistency care; still all-or-nothing |
| C3 | cap injected rotation magnitude | clamp ‖dx[3:6]‖ <= theta_max | dead-simple | ad-hoc, biases, Joseph inconsistent |
| C4 | **residual-gated rotation suppression** | suppress dx[3:6] ONLY when the position residual is large | keeps the legit info path in the normal regime; kills only the pathological kick | one threshold to tune |

C4 is the principled choice: the position->heading coupling is valid only under small-residual
linearization, so disable it exactly when that assumption breaks (large residual), keep it otherwise.
C1/C2 throw away legitimate slow heading observability; C3 is a blunt clamp.

## 4. Chosen approach — C4: residual-gated rotation-row suppression

In `update()` / `update_aug()`, for a correction that does NOT directly observe rotation
(`H` rotation columns ~ 0 — detect `H.middleCols(3,3).norm() < eps`):

```
// after computing the gain K and the per-DOF RMS Mahalanobis dbar = sqrt(d2/n):
if (rotation_unobserving && dbar > rot_suppress_kappa) {
    w_rot = rot_suppress_kappa / dbar;      // in (0,1): smooth, not a hard cut
    K.block(rotation_rows) *= w_rot;        // attenuate ONLY the rotation-error gain rows
}
dx = K r;                                   // translation rows unchanged -> position still pulled
// Joseph uses the SAME modified K  -> covariance stays consistent (rotation block barely reduced)
```

- Translation correction is untouched (the useful 51 m pull survives).
- Rotation injection is smoothly down-weighted as the residual grows; for the t=1929 fix
  (dbar 1.68, say kappa_rot=0.8) w_rot=0.48 -> 68 deg kick becomes ~33 deg; pushing kappa_rot lower
  (e.g. 0.3 -> w_rot 0.18 -> ~12 deg) bounds it further. Tune against the acceptance metric.
- `correction_rot_suppress_kappa` config knob; 0 = DISABLED -> bit-identical to today (opt-in).
- Joseph form with the modified K stays PSD (K scaling preserves the (I-KH)P(I-KH)^T + KRK^T form).

This is ORTHOGONAL to Slice-15 Huber (A): A scales the whole gain on the truly-huge outliers; C4
selectively protects heading on in-gate-but-large position residuals. They compose.

## 5. Acceptance test (objective via the local metric + the diagnosis tool)

1. **The t=1929 s kick is bounded**: with C4 on, `diagnose_window.py` shows the applied fix rotates
   heading by << 68 deg (target < ~10 deg) while still pulling translation.
2. **urban12 divergence bounded**: global tail far below 3 km; `local max` 476 -> materially lower;
   ideally the death spiral is averted (heading never runs to 2 rad).
3. **urban07 / urban17 not regressed**: `local p50` within +/-10% (the normal-regime position->heading
   path is untouched there — their fixes are small-residual).
4. **Sim unchanged at kappa=0**: `test_cov_calibration` / `test_gps` bit-identical (opt-in).
5. New unit test: a position fix with a planted trans-rot cross-cov and a LARGE residual injects a
   bounded yaw with C4 on (and the old large yaw with C4 off); a SMALL-residual fix is unaffected
   (rotation still corrected) — proving the residual gate, not a blanket cut.

## 6. Risks

- Strict-core edit (same sensitive path as Slice 15). Mitigate: kappa=0 inert default; Joseph PSD.
- Over-suppression could slow legitimate heading correction from position in GPS-sparse driving.
  Mitigate: the residual gate keeps small-residual fixes fully effective; only large residuals are
  attenuated.
- May not FULLY save urban12: heading was already ~1 rad wrong from earlier kicks + the 522 s coast
  before the fatal fix. C4 should convert a 3 km runaway into a bounded (large-but-recoverable)
  error; full recovery may also need the upstream issue (see below). Acceptance #2 measures this.

## 7. Out of scope (separate work)

- **The 522 s heading drift during the GPS coast** (DR/fusion holding heading ~1 rad over 8.7 min of
  dense urban) — a derived-source / FOG-IMU heading-hold quality issue, not a correction-update bug.
- Huber (A, Slice 15) — kept as-is; C4 composes with it.
- Lever B (hard rot clamp) — superseded by C4's smooth residual-gated form.

## 8. OUTCOME — ACCEPTANCE MET (`ba895b6`)

Full-drive sweep, keeper base (realistic R `cov_floor_m2=25`, default gate). C4 = `correction_rot_suppress_kappa`.

| urban12 | tail m | max m | rms m | local max m | NEES | gps_applied |
|---------|--------|-------|-------|-------------|------|-------------|
| baseline (no R, no C4) | 4214 | 4596 | 2648 | 856 | 59907 | 96 |
| D only (cov25) | 3075 | — | — | 476 | 225 | 188 |
| **D + C4 kappa=0.8** | **1.99** | 409 | 120 | 422 | 371 | **1375** |
| D + C4 kappa=0.5 | 1.99 | 304 | 106 | 313 | 367 | 1383 |

**C4 is the first lever that fixes urban12.** Tail 4214 -> **1.99 m**; the death spiral is averted
(gps_applied 96 -> 1375 — heading no longer runs away, so fixes stay in-gate and keep correcting).
`diagnose_window.py` confirms the t=1929 s 68 deg kick is GONE (rot_err there 2.0 -> 0.08) and the
whole arc recovers to ~2 m by the end. A mid-drive transient excursion remains (max ~300-400 m, one
window) but now RECOVERS instead of carrying forever — fully flattening it needs the upstream 522 s
GPS-coast heading drift (Sec. 7, out of scope).

**No regression on urban07/17** — C4 is *inert* there: their accepted fixes have dbar ~0.18-0.55,
below kappa=0.8, so C4 never fires. (urban07 cov25+C4: local p50 4.68->2.37, NEES 388->28; urban17:
NEES 4438->321, local p50 0.23->0.35 — the latter a `cov_floor` (D) tradeoff, not C4.) Sim
bit-identical at kappa=0 (full suite green). Validated on FULL drives — the local metric guards
against a hidden blowup, so no truncated-slice over-claim.

**Recommended config**: `cov_floor_m2=25` + `correction_rot_suppress_kappa=0.8` (safe default, clear
of urban07/17 dbar). kappa=0.5 lowers the urban12 transient max (304 vs 409) but risks touching
urban07 (dbar 0.55). `cov_floor_m2` is an independent per-deploy consistency/accuracy knob.

## 9. Open questions

- `rot_suppress_kappa` value: start ~0.8 (just below the in-gate dbar ceiling 1.73) and sweep down.
- Detect "rotation-unobserving" by `H` rotation-column norm vs a config flag on the correction type?
  Norm test is automatic and correct for any ICorrection; confirm no false positives for 6-DOF pose
  fixes (those legitimately observe rotation -> H rotation cols nonzero -> not suppressed).
- Should suppression scale to 0 (full cut) past a second, larger residual threshold, or asymptote via
  the 1/dbar law? The smooth 1/dbar form is simplest; revisit if urban12 needs a harder cut.
