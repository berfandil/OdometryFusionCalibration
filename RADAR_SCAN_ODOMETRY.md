# Radar scan-matching odometry (descriptor-based) — a project side-part

A radar odometry front-end that recovers **full rigid motion (rotation + translation)** from sparse radar detections by matching rotation/translation-INVARIANT per-point descriptors across consecutive scans. User-contributed algorithm (2026-06-15).

**Why it matters here**: the Doppler ego-velocity front-end (`tools/nuscenes_to_csv.py`) gives a radar TRANSLATION only (`R_B = I`, heading-blind) — which is exactly why the radar could not self-calibrate its rotation extrinsic and the hand-eye lever stalled (Slice 20/20b: `translation_only`, lever blocked by per-window noise + missing rotation). This algorithm gives the radar a real frame-to-frame `R` AND `t`, so the radar becomes an ordinary full-6DOF (here full-3DOF planar) source: rotation extrinsic + lever both observable through the existing calibrators, no `translation_only` flag. It is the principled radar-front-end fix flagged at the end of Slice 20b ("committing a real radar lever is now a radar-front-end odometry-quality problem").

---

## 1. The 2D algorithm

Inputs per step: the CURRENT radar detection set and the PREVIOUS one (both already static-filtered, step A); the reference sensor's odometry; the prior extrinsics of the reference sensor (`X_ref` = base_from_ref) and the radar (`X_radar` = base_from_radar). Assumes close-enough calibration on both sensors and good-enough reference odometry.

### A. Static-object filtering (reject movers)
The descriptor geometry must come from STATIONARY world points, so movers are removed first:
1. From the reference odometry over `[t_prev, t_curr]` (its own frame delta `B_ref`), form the base motion `A = X_ref ∘ B_ref ∘ X_ref⁻¹`, then the radar's expected ego-motion `B_radar = X_radar⁻¹ ∘ A ∘ X_radar`. Its instantaneous twist `(v, ω) = log(B_radar)/dt`.
2. For each detection `i` at sensor-frame position `r_i` with measured (radial) velocity `ṽ_i`, the velocity a STATIC world point would show is `v_static(r_i) = −(v + ω × r_i)` (the ego-motion seen from a fixed point; `ω × r_i` is the rotational/lever term). Keep `i` iff the measured radial velocity matches the static prediction within a threshold: `| proj_radial(ṽ_i) − proj_radial(v_static(r_i)) | < τ_v`. Detections that disagree are moving objects → dropped.
   - This is the same static/dynamic split as the Doppler-RANSAC front-end, but it uses the REFERENCE ODOMETRY as the prior instead of self-fitting — cheaper + more robust when the reference odometry is trustworthy.
3. The PREVIOUS frame was filtered the same way in its own turn (so both sets are static-only). The very first radar frame has no usable predecessor → skipped (harmless).

### B. Per-detection descriptor
For each static detection `i` in a frame, the descriptor is the SORTED multiset of Euclidean distances to every other static detection in the SAME frame:
```
desc(i) = sort({ ‖r_i − r_j‖ : j ≠ i })
```
The multiset of inter-point distances of a point is invariant under any rigid motion (rotation + translation) — so the same world point keeps (approximately) the same descriptor across two scans, regardless of how the radar moved.

### C. Cross-frame matching (linear, sorted)
Match a current detection `i` to a previous detection `k` by comparing `desc(i)` with `desc(k)`. The two descriptors have DIFFERENT lengths (different detection counts) and shared distances are not at the same indices — but BOTH are sorted, so the comparison is a single LINEAR merge walk: advance two pointers, count distance pairs `(a ∈ desc(i), b ∈ desc(k))` with `|a − b| < τ_d` as "agreeing" (advance the smaller; a matched pair advances both). The similarity score = the agreeing count (optionally normalized by `min(len_i, len_k)`).
- **Partial-overlap robustness (refinement)**: between consecutive scans the static-point sets only PARTIALLY overlap (points enter/leave the field of view), so a common point's neighbor set — hence its descriptor — shifts. The tolerance-merge COUNT (not exact-equality) handles this: a true match still agrees on the MANY distances to the common neighbors even if some distances (to the non-shared points) do not align.
- Accept `i↔k` as a correspondence when it is a MUTUAL best match and the similarity clears a floor (`≥ ρ·min(len_i, len_k)` agreeing distances).

### D. Motion recovery
The accepted correspondences `{(r_i^curr, r_k^prev)}` are the same world points seen at two times, so the radar's frame-to-frame motion `B_radar = (R, t)` is the rigid transform with `r^prev ≈ R·r^curr + t`. Solve by Umeyama/Kabsch (rotation + translation, NO scale), wrapped in **RANSAC** (refinement: sample minimal pairs — 2 points in 2D — score inliers by residual, refit on inliers) to reject the descriptor mismatches that survive C. Output `B_radar` as the radar source's odometry increment over `[t_prev, t_curr]` — now with a real `R`, not identity.

**Complexity**: building all descriptors is `O(N² log N)` (N detections/frame, ~50–200 for ARS408); matching is `O(N²)` merge walks of `O(N)` each → `O(N³)` worst case, trivial at radar N. Bounded, no heap concerns for an offline front-end.

---

## 2. Extension to 3D (straightforward)

The algorithm is dimension-agnostic except the motion solve and the static-filter velocity:
- **Descriptor (B)** — unchanged: `‖r_i − r_j‖` is the 3D Euclidean distance; the sorted-distance multiset is rigid-invariant in 3D exactly as in 2D. Matching (C) is byte-identical (still a linear merge of sorted scalar lists).
- **Motion recovery (D)** — Umeyama/Kabsch in 3D → `SE(3)` (`R ∈ SO(3)`, `t ∈ ℝ³`); RANSAC minimal sample = 3 non-collinear pairs. The only code change is the dimension of the point matrices.
- **Static filter (A)** — needs the 3D ego-velocity `−(v + ω × r_i)` with `r_i ∈ ℝ³` and a 3D measured velocity. This is the ONE part that needs a **3D radar** (elevation + `vz`): the nuScenes Continental ARS408 is 2D (all `z = 0`, velocity `vx,vy` only — verified), so on it the z-lever and pitch/roll stay unobservable from the radar itself. A 4D/imaging radar (elevation + Doppler in 3D) drops straight in: same descriptor, same matching, 3D Kabsch, 3D static filter.

So 3D is a mechanical lift (point dimension + Kabsch dimension); the descriptor/matching core — the clever part — is unchanged.

---

## 3. Fit into the project

- **Side-part deliverable**: `tools/radar_scan_odometry.py` (relaxed-edge, the `tools/` prototype tier) — a 2D radar odometry front-end with a 3D-ready structure (point dimension parameterized). Reads nuScenes radar scans + the reference (CAN) odometry + extrinsics, emits a radar increment-form CSV WITH rotation (the project's standard schema), as a drop-in alternative to the Doppler front-end in `nuscenes_to_csv.py`.
- **Validation**: per-frame recovered `B_radar` vs the GT frame-delta (does the scan-match track the true radar motion, and does it now carry a real heading?); then optionally feed it through the calibrator WITHOUT `translation_only` and check whether the radar rotation extrinsic + lever become observable (the unblock).
- **Not core (yet)**: the radar front-end is a sensor adapter, not strict-core fusion math; it lives in `tools/`/`adapters/` like the other converters. If it proves out, a hardened `adapters/` version is the productionization step.

## 4. Status
- [x] Algorithm documented (this file) + 3D notes — user-contributed, orchestrator-captured.
- [x] Implemented (`tools/radar_scan_odometry.py`, `61b79f3`; 2D, 3D-ready; both refinements — tolerance-merge matching + RANSAC). Core math proven EXACT (synthetic Kabsch/RANSAC round-trips a planted 3 deg to 1e-9). Two real bugs found+fixed: the CAN base delta integrated ZERO samples per 75 ms radar window (CAN is ~500 ms) -> use the instantaneous nearest-in-time CAN twist held over dt; + a plausibility gate rejecting RANSAC-blunder yaw/translation beyond the CAN prior + margin.
- [x] Validated on nuScenes radar (scene-0061/0103, all 5 radars):
  - **Translation tracks** (median 0.5-1.0 m/frame ≈ the ~0.9 m GT step; step-length corr ~0.4).
  - **Rotation NOT usable** on sparse 2D ARS408: the true 75 ms inter-scan yaw (~0.2 deg) sits BELOW the ~0.8 deg Kabsch rotation-noise floor on ~9 noisy inlier points -> cumulative recovered heading scatters (-53 to +23 deg vs GT net +98 deg), per-frame yaw correlation ~0. Match stats: dense radars 44-84 static / 20-49 corr / 9-33 inliers / 8-45% identity-fallback; sparse corner radars worse (33-86% fallback).
- [x] (stretch) Fed through the calibrator as a non-`translation_only` radar: the rotation channel is STRUCTURALLY unblocked — radar `extr_rot` confidence becomes NONZERO (0.11-0.46) vs the Doppler `translation_only` baseline's structural conf 0. BUT the noisy per-step rotation drives the recovered extrinsic to ~1 rad nonsense, collapses the radar scale, and wrecks fused drift (max 801 m, NEES 78). The heading-blind wall is removed structurally; the limiter is now front-end rotation SNR.

**Verdict**: the algorithm is correct and gives usable radar TRANSLATION + a structural rotation unblock, but per-step ROTATION is below the noise floor on sparse 2D 13 Hz radar. **Next levers** (rotation, in order): (1) **longer matching baseline** — match scan k to k−N so the accumulated yaw (N x ~0.2 deg) rises above the ~0.8 deg Kabsch floor while the FOV still overlaps (the most promising, cheap fix; trades overlap for signal); (2) temporal smoothing of the recovered R; (3) a 4D/imaging radar (§2: more points + elevation -> 3D + a far better-conditioned Kabsch). Translation is already usable as a drop-in radar source today. Committing a real radar LEVER remains gated on this front-end rotation quality (Slice 20b).
