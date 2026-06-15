# Slice 20 Review — Translation-Only Source Lever

**Commit**: 629d379

## Findings

### 1. Byte-Identical Default-Off Path

src/core/calibration.cpp:275–276: VERIFIED CORRECT
- When translation_only=false, `if (!trans_only_[slot])` guard ensures phase1_direction_vote() is called.
- On guard drop (return false), caller does `continue` past scale vote — byte-identical to pre-20 inline `continue` paths.
- When translation_only=true, direction block skipped, scale vote always reached.
- Implementer's self-review correctly identified and fixed initial bug (return void → return bool).

### 2. Rotation Extrinsic Pinned to Prior

src/core/calibration.cpp:1003–1007 (R_yp forcing): VERIFIED CORRECT
- R_yp = trans_only ? prior_[slot].R : (ryp_set_[slot] ? ryp_[slot] : yaw_pitch_of(prior_[slot].R))
- When translation_only=true, R_yp forced to FULL prior (including roll), so R_X = R_yp * Rx(0) = prior exactly.

src/core/calibration.cpp:1012–1016 (Roll voting skipped): VERIFIED CORRECT
- roll = trans_only ? Scalar(0) : best_roll(R_yp, R_A, R_B)
- Roll vote skipped; roll stays 0 in lever row.
- if (!trans_only) roll_[slot].add(roll, w) — roll histogram untouched.

src/core/calibration.cpp:1040–1044 (rot3d block gated): VERIFIED CORRECT
- if (rot3d_enabled_ && !trans_only) — rot3d Wahba accumulation/vote skipped.
- rot3d_row stays false → R_X pinned at prior, not corrupted by Wahba.

src/core/calibration.cpp:1294–1299 (extrinsic() readout): VERIFIED CORRECT
- R_yp = trans_only_[s] ? prior_[s].R : (ryp_set_[s] ? ryp_[s] : prior_[s].R)
- Readout also forces R to prior for translation_only sources.

### 3. Footgun Fix — Rotation Never Commits Off Prior

src/core/estimator.cpp:765–768 (ext_committed forcing): VERIFIED CORRECT
- ext_committed[i] = !tonly && commit_gate_reanchor(...)
- Yaw/pitch commit forced false for translation_only sources.

src/core/estimator.cpp:774–777 (roll_committed forcing): VERIFIED CORRECT
- roll_committed[i] = !tonly && commit_gate_reanchor(...)
- Roll commit forced false.

src/core/estimator.cpp:795–802 (rot3d_committed forcing): VERIFIED CORRECT
- if (cfg.rot3d_enabled && !tonly) { rot3d_committed[i] = ... } else { rot3d_committed[i] = false; }
- rot3d commit forced false. Belt-and-suspenders: calibrators skip votes AND estimator forces commits false.

### 4. Persistence / Restore Path (CONFIG-HASH GUARDS)

src/core/estimator.cpp:268–274 (config_hash includes translation_only): VERIFIED CORRECT
- w.put_bool(sc.translation_only) included in per-sensor hash block.
- Blob written with translation_only=true has different hash than identical rig with translation_only=false.
- Hash mismatch rejects restore at deserialization.
- Prevents stale blob from non-translation_only rig restoring rotation into translation_only slot.

Restore orthonormality guard (line 2319–2327): APPROPRIATE
- Existing guard validates all restored rotations are orthonormal (RᵀR ≈ I, det > 0).
- Applied uniformly; no special case needed (config-hash already blocks cross-flag restores).

### 5. Participates() and ReferenceOnly Cold-Start

src/core/estimator.cpp:701–706 (participates carve-out): VERIFIED CORRECT
- if (trans_only[i]) return lever_committed[i]
- Without this, ReferenceOnly would require committed rotation — but both ext_committed and rot3d_committed pinned false.
- With carve-out, translation_only sources join median once lever commits (the DOF they actually recover).
- Rotation IS trusted "by declaration (the prior)" — correct logic.

### 6. Validate() Orthonormality Check

src/core/estimator.cpp:2511–2524 (validate check): VERIFIED CORRECT
- Gated on if (cfg.sensors[i].translation_only)
- Checks RᵀR ≈ I (tolerance 1e-6) and det > 0.
- Rejects at init; matches restore path tolerance (line 2322).

### 7. Scale Calibration Preserved

src/core/calibration.cpp:280–282 (scale vote): VERIFIED CORRECT
- Scale vote independent of translation_only flag.
- Runs for all sources with scale_calib_[slot] = true.
- Phase-1 and Phase-2 scale both untouched.

### 8. Along-Track Withheld (Concentration Gate Suffices)

Test case 3: VERIFIED VIA TEST
- Under 0.5 m/window along-track noise, per-channel concentration collapses.
- translation_confidence falls below commit gate (0.6).
- No new residual gate added (spec decision: concentration already withholds).

### 9. Config-Hash and Loader Key

adapters/src/config_loader.cpp:408–414: VERIFIED CORRECT
- Key translation_only=<bool> parsed.
- Default false (byte-identical).
- Comment documents Slice 20 intent and prior-rotation requirement.

### 10. Strict-Core Compliance

Fixed-capacity arrays:
- trans_only_[kMaxSources] in Phase1Calibrator (calibration.hpp:258)
- trans_only_[kMaxSources] in Phase2Calibrator (calibration.hpp:680)
- trans_only[kMaxSourcesCap] in Estimator::Impl (estimator.cpp:456)
- All zero-initialized on ensure_slot/init.

No heap post-init: VERIFIED
Bounded WCET: VERIFIED (validate orthonormality = fixed 3×3 op)
Status returns: VERIFIED

### 11. Tests (test_calib_translation_only.cpp, 8 cases)

Test 1 (Headline): VERIFIED — flag OFF lever err > 0.02m; flag ON err < 2e-3m + rotation stays = prior (dev < 1e-12 rad)
Test 2 (Footgun rot3d enabled): VERIFIED — rot3d never votes (count=0), rotation stays = prior (dev < 1e-12 rad)
Test 3 (Along-track noisy): VERIFIED — translation_confidence collapses below gate under drift
Test 4 (Planar null): VERIFIED — yaw-only: x,y observable, z withheld by conditioning gate
Test 5 (Scale preserved): VERIFIED — Phase-1 scale votes despite direction skipped
Test 6 (Estimator integration): VERIFIED — radar lever recovers (err < 0.04m), rotation never leaves prior (dev < 1e-9 rad), extrinsic_committed=false
Test 7 (Validate + hash flip): VERIFIED — garbage prior rejected when flag=true; hash flip blocks cross-flag restore
Test 8 (Observability): VERIFIED — no turns: lever unvoted; multi-axis: all 3 axes observable

Byte-identical pin: Existing test suite untouched (CMakeLists.txt only).

## Verdict

APPROVE

Implementation correct, comprehensive, meets all 8 spec acceptance criteria. phase1_direction_vote extraction properly fixed to preserve byte-identity on default-off path. Rotation extrinsic provably pinned to prior by multiple layers (calibrator skips votes, estimator forces commits false, readout pins rotation, config-hash rejects cross-flag restores).
