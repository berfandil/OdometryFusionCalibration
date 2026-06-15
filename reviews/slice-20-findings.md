# Slice 20 review findings — translation-only source lever

Reviewer: caveman:cavecrew-reviewer (diff review) + orchestrator self-review.

## MAJOR (found in self-review, FIXED before commit)

- src/core/calibration.cpp (Phase1Calibrator::observe / phase1_direction_vote): 🔴 MAJOR:
  the `phase1_direction_vote` extraction initially used `return` (void) on the two direction
  guards (degenerate fold `gn < kUnitEps`; the π-singularity `g_unit.dot(kFwd) < 0`), and the
  caller then voted scale unconditionally afterward. In the PRE-20 inline code those guards were
  `continue`s that DROPPED THE WHOLE SAMPLE — including its scale vote. So a non-translation_only
  source's reverse sample (fold off) or near-π folded direction would have voted scale where the
  original voted nothing — a default-OFF behavioral regression (the existing tests passed only
  because their trajectories never hit those guards). FIX: `phase1_direction_vote` now returns
  `bool` (false = a guard dropped the whole sample), and the caller does
  `if (!phase1_direction_vote(...)) continue;` for a non-translation_only source — byte-identical
  to the pre-20 `continue`s. A translation_only source skips the direction block entirely and
  always reaches the scale vote (its direction is irrelevant). Re-verified green.

## Reviewed clean (no change)

- CORRECTNESS: R_X pinned to the FULL prior rotation (R_yp = prior_[slot].R, roll = 0, rot3d
  block gated `&& !trans_only`) -> the lever row uses R_X = prior exactly. Phase1 skips ONLY the
  direction vote, keeps scale. Estimator forces ext/roll/rot3d commits false for a
  translation_only source (`!tonly && commit_gate_reanchor(...)`); snapshot extrinsic_committed
  stays false. Phase2::extrinsic() pins R = prior for the flag.
- BYTE-IDENTICAL default-off: all new behavior behind `if (translation_only)` /
  `if (trans_only)` / `!trans_only` / `!tonly`; full existing suite (golden + cov-cal) untouched
  and green.
- STRICT CORE: fixed-capacity `trans_only_[kMaxSources]` (both calibrators) +
  `trans_only[kMaxSourcesCap]` (estimator), zero-initialized in ensure_slot/init; no heap; the
  validate() orthonormality check is a fixed 3x3 op; bounded loops.
- config_hash: `translation_only` hashed in the per-sensor block (a flip rejects a stale
  restore; pre-20 blobs cold-start). validate(): a translation_only source's prior rotation must
  be orthonormal (RᵀR ≈ I, det > 0), gated on the flag — a non-translation_only source tolerates
  a rough prior (the recovered-rotation path).
- participates(): a translation_only source never commits a rotation, so under ReferenceOnly
  cold-start it joins the median on `lever_committed` (its rotation is trusted at the prior). It
  is still fed calibration before joining (EVERY covered source feeds calibration, line ~1371),
  so no cold-start deadlock.
- persistence restore: a self-consistent blob from a translation_only config always has
  ext/roll/rot3d committed = false (the feedback forced them); restoring those false flags is
  correct, and the first post-restore step re-forces them false anyway. Config-hash guarantees
  flag agreement across restore.

## Along-track guard decision (measured)

The existing histogram-CONCENTRATION gate ALREADY withholds the noisy/along-track axis — NO new
residual gate was added (SLICE20 path (b), simpler). See the test
`trans-only along-track: concentration gate already withholds the noisy axis`: under 0.5 m/window
along-track sensor drift the per-channel concentration collapses and the MIN-over-channels
`translation_confidence` -> ~0, FAR below `commit_concentration` (0.6), so the wrong along-track
value is never committed. The clean control (no noise) concentrates and commits the lever.
