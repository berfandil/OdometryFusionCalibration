# Slice 16 — Centroid sub-bin readout (calibration precision floor)

**Goal**: committed calibration precision toward the user target (rotation < 0.1° = 1.745e-3 rad), by removing the dominant precision floor in `Histogram1D::mode()`.

---

## 1. Problem + evidence (2026-06-10 investigation)

The Slice-15c committed-yaw error (~0.23° vs the 0.10° uncommitted estimate) was hypothesized to be commit-dynamics or vote-weighting. A discriminator experiment on the KAIST `calib_test` rig (urban07 wheel ×3, injected yaw on src2) ruled both out and found the real cause.

**Discriminator (truth yaw = 0.174533 rad):**

| Run | rz read | err |
|---|---|---|
| `vote_weight=one`, commit ON | 0.17849 | +0.229° |
| `vote_weight=one`, commit OFF | 0.172707 | −0.103° |
| `vote_weight=combo`, commit OFF | 0.172741 | −0.101° |

Weighting is irrelevant (one ≈ combo to 3e-5 rad). The committed/uncommitted gap is NOT a feedback defect either — see below.

**Root cause — parabolic sub-bin bias.** `so3_hist` = 64 bins over [−1,1] rad → bin width 0.03125 rad = 1.79°. All target residuals live deep inside one bin, so precision rides entirely on the sub-bin readout. `mode()` does parabolic interpolation over (peak−1, peak, peak+1). For a point mass vote-split between two bins at fraction f above the peak-bin center, the parabola reads `0.5f/(2−3f)` bins instead of `f` — at f=0.25 it reads 0.1: a systematic **pull toward the peak-bin center of up to ~70–80% of the offset**, worst around quarter-bin positions (~0.3–0.45° worst case at this bin width).

**Falsification probes (uncommitted, combo; manifest `kaist_run/calib_binc.ini` / `calib_q18.ini`):**

| Truth (rad) | sub-bin position | read | err |
|---|---|---|---|
| 0.171875 (= exact bin-37 center) | 0.00 | 0.171953 | **+0.0045°** |
| 0.174533 (original) | +0.084 bins | 0.172707 | −0.103° |
| 0.180 | +0.26 bins | 0.175399 | **−0.264°** |

Error is a deterministic function of the truth's position within its bin — quantization bias, not noise. The "0.10°" Slice-15c number was luck (10° ≈ a bin center).

**Same bias in every DOF of the same run** (~70–76% pull toward the peak-bin center):
- scale: truth 1.1, peak-bin center 1.1015625, read 1.10112;
- time-offset: truth −0.05 s, peak-bin center −0.046875, read −0.0478.

**Why commit made it look worse**: the rising-edge re-anchor moves the histogram basepoint, so the post-reset residual lands at a different (uncontrolled) sub-bin phase → a different deterministic bias (+0.23° instead of −0.10°). Commit dynamics are otherwise healthy (lever even improved under commit: [0.509, 0.203] vs [0.545, 0.071] uncommitted, truth [0.5, 0.2]).

## 2. Design

**Mechanism**: add an opt-in centroid sub-bin readout to `Histogram1D` — the mass-weighted mean over the peak bin and its two immediate neighbors:

```
mode_centroid = Σ_{i∈{p−1,p,p+1}} m_i · c_i / Σ m_i
```

- **Exact for a vote-split point mass** at any sub-bin position (the linear vote split is itself the sub-bin encoding: masses are linear in the vote's position, and their 2-bin centroid reconstructs it exactly — the parabola does not).
- For a concentrated noisy distribution (our regime: `confidence()` ≈ 0.99 means ≤~1% of mass outside peak±1) the 3-bin truncation bias is ≤ ~0.02 bins ≈ 0.04° at the so3 bin width — under target.
- Window = peak±1 deliberately matches `confidence()`'s 3-bin concentration notion.

**Placement**: new field `bool subbin_centroid = false` in `HistogramConfig` (honored only when `subbin == true`; `subbin == false` still returns the bin center). Default OFF → every existing config byte-identical; the sim and the tuned trust apparatus are untouched.

**Edge handling**:
- Circular: wrap-aware neighbor indices; compute neighbor centers as `center(p) ± width` (unwrapped), fold the centroid back into range. A seam-straddling peak must read continuously.
- Non-circular boundary peak (p=0 or p=nbins−1): missing neighbor contributes zero mass and is excluded from both sums — centroid degrades gracefully (no parabola-style fall-back-to-center cliff).
- Empty histogram: unchanged (midpoint).

**Config plumbing**:
- `ConfigLoader`: one `[global]` key `subbin_centroid = true|false` applied to all five calibration histogram configs (`so3_hist`, `roll_hist`, `xyz_hist`, `scale_hist`, `offset_hist`).
- **Persistence config-hash**: `subbin_centroid` is calibration-shaping → include it in the FNV config-hash stream (a flipped flag must reject a stale warm-restore). This changes the hash for ALL configs vs prior builds → old persisted files reject with `VersionMismatch`-style cold start. Intentional; note in DECISIONS.
- Doc-affecting: CONFIG.md (HistogramConfig table + loader key), DECISIONS.md (this design), DESIGN.md (histogram primitive description).

**Out of scope (fallback if centroid empirically insufficient)**: zoom-on-commit re-anchor (narrow the histogram range post-commit for a finer bin width). Not built now — centroid alone is expected to clear the target; zoom adds range-vs-noise-spread coupling and re-open complexity.

## 3. Acceptance

Unit (TDD; `tests/test_histogram.cpp` + loader + persistence tests):
1. Centroid: single split vote at sub-bin positions {0, 0.1, 0.25, 0.4, 0.5} of a bin → `mode()` returns the exact value (≤1e-12). Pin that Parabolic reads the f=0.25 case with ≥0.1-bin error (documents the rationale).
2. Symmetric multi-vote spread → centroid ≈ true mean (≪ parabola error).
3. Circular seam: peak straddling the wrap reads continuously.
4. Boundary peak (non-circular): no NaN, graceful one-sided centroid.
5. Default-off: `subbin_centroid=false` → bit-identical `mode()` vs current (existing tests stay green untouched).
6. ConfigLoader parses the key; round-trips into all five histogram configs.
7. Config-hash flips when the flag flips (persistence rejects cross-flag restores).

Gate: `scripts/dev.ps1 -Task test` green; cov-cal NEES unchanged (~4.82, band [4.0, 5.6]) — default is OFF so byte-identical.

Real-data validation (orchestrator, post-merge):
- Rerun the three probes (`calib_combo_nocommit.ini`, `calib_binc.ini`, `calib_q18.ini`) with `subbin_centroid=true` → ALL yaw reads within **<0.1°** (expect ≤0.02°); scale err ≤5e-4; time-offset err ≤1 ms.
- Rerun `calib_vw1_commit.ini` with the flag → **committed** yaw err <0.1° (the actual Slice-16 goal), scale/offset committed and tightened.
- No degradation of lever recovery.

## 4. Status

- [x] Implemented (TDD, gate green, committed) — `c237bfb`; unit 214→221 cases, adapters 39
- [x] Reviewed (`reviews/slice-16-findings.md`: APPROVE, 1 MINOR + 3 NIT) + all findings fixed — `00978c8`
- [x] Real-data validation — see table below
- [x] Docs updated (CONFIG/DECISIONS D25/DESIGN/ISSUES Slice 16) — orchestrator

### Real-data validation results (KAIST urban07 calib_test, `subbin_centroid=true`, 2026-06-10)

| Run | truth yaw (rad) | read | err BEFORE | err AFTER |
|---|---|---|---|---|
| combo, no commit | 0.174533 | 0.174791 | −0.103° | **+0.0148°** |
| bin-center probe | 0.171875 | 0.172133 | +0.0045° | +0.0148° |
| quarter-bin probe | 0.180 | 0.180258 | −0.264° | +0.0148° |
| **one, COMMITTED** | 0.174533 | 0.174711 | +0.229° | **+0.0102°** |

Error is now sub-bin-phase-independent (identical +2.6e-4 rad across probes — the quantization bias is gone; the residual is a small common systematic). Other DOFs, committed run: scale 1.10008 vs 1.1 (err ~8e-5, was 1.8e-3); time-offset −0.0500026 vs −0.05 (err **2.6 µs**, was 2.2 ms); lever [0.5007, 0.1975] vs [0.5, 0.2] (~2.5 mm, translation target met on the observable axes). **Committed yaw 0.0102° — 10× under the 0.1° rotation target** on this calibration DOF. Manifests: `kaist_run/calib_*_cen.ini`.
