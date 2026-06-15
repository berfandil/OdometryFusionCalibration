# Slice 20b Review — Per-Axis Lever Commit

**Commit**: c85442e
**Spec**: SLICE20_TRANSLATION_ONLY_LEVER.md section "4b"
**Reviewer**: Claude Code

## Summary

Per-axis lever commit for `translation_only` sources. Each lever axis gates independently when its concentration clears the threshold, enabling clean axes to commit while noisy/unobservable axes stay at prior. Non-flagged sources remain byte-identical. Persistence v3→v4 with per-axis flags appended after scale2_committed.

## Findings

### 1. Persistence Version Guard ✓

**File**: include/ofc/core/persistence.hpp:73, src/core/estimator.cpp:2319
- Version bumped 3→4. Deserialize checks `version != kFormatVersion` before any state touched.
- **Test**: test_calib_translation_only.cpp L955-966 stamps v3 on v4 blob, confirms VersionMismatch rejection.
- **No OOB**: version check pre-commit, so truncated v3 blobs cannot cause read beyond bounds.

### 2. Per-Axis Flags Schema Position ✓

**Files**: serialize (estimator.cpp:2256-2258), deserialize (estimator.cpp:2375-2377)
- Schema: 5 base + rot3d + scale2 + 3 per-axis levers = 10 flags.
- **Test pins**: test_calib_translation_only.cpp L889, test_calib_rot3d.cpp L326 expect 10 flags.
- **Blob parse test**: L878-905 manually parses blob, confirms flags[7,8,9] at correct positions.

### 3. Whole Flag Re-Derived as OR on Restore ✓

**File**: estimator.cpp:2437
- Stored `rc.lever_c` discarded; whole = `rc.lever_c_xyz[0] || rc.lever_c_xyz[1] || rc.lever_c_xyz[2]`.
- **Invariant**: `lever_committed = OR(lever_committed_xyz)` provably consistent after restore.
- **Test**: L927-931 confirms restored flags round-trip and OR maintained through refill.

### 4. Byte-Identical for Non-Flagged ✓

**File**: estimator.cpp:797-814
- Whole-lever path unchanged (else branch L808-813).
- Per-axis array mirrors whole flag (L812-813, all 3 = whole).
- Per-axis publish: all 3 true ⟺ old `t = lever_arm(id)` (L887-891).
- **Test**: L785-797 pins per-axis == whole at every step for non-flagged rig.
- **Byte-identical proof**: L945-950 round-trip serialization is identical.

### 5. Per-Axis Gate Correctness ✓

**File**: estimator.cpp:799-804
- Each axis gates on ITS OWN confidence + mass via commit_gate_reanchor (same hysteresis, per-axis).
- **Accessors** (calibration.cpp:1347-1361): `translation_confidence_axis(id,c)` → `xyz_[3*s+c].confidence()`, `xyz_axis_mass(id,c)` → `.total()`.
- Guards: configured, axis ∈ [0,2], valid slot → return 0 for out-of-range (safe).
- Empty channel → 0 (histogram's own convention).
- **Test**: L653-675 measures per-axis confidence (cx≥gate, cy<gate, cz=0); L663-665 measures mass.

### 6. Per-Axis Publish ✓

**File**: estimator.cpp:887-891
- Only committed axes folded: `if (lever_committed_xyz[i][c]) prior_extrinsic[i].t[c] = calib2.lever_arm(id)[c]`.
- Uncommitted axes stay at prior.
- **Test**: L872-909 parses blob, confirms x moved off prior, y,z exactly at prior.

### 7. Reference Source ✓

**File**: estimator.cpp:749-756
- Reference: all commit flags false (including per-axis L752-753).
- Gauge stays fixed; per-axis array consistent.

### 8. Init() ✓

**File**: estimator.cpp:1040-1042
- Per-axis levers initialized false (consistent with whole flag L1039).
- No heap, fixed-capacity loop.

### 9. Snapshot Diagnostics ✓

**File**: include/ofc/core/result.hpp:93
- `translation_committed_xyz[3]` defaults false. Comment explains semantics.
- **Test**: L725-732 confirms per-axis split published; L793-795 confirms mirror for non-flagged.

### 10. rot3d + reset_lever ✓

**File**: estimator.cpp:907-923
- On rot3d commit RISING EDGE, lever rows/histograms reset.
- Per-axis flags held through refill via commit_gate_reanchor hysteresis.
- **Contract**: no double-fold, no drift. Lever reset post-rot3d, per-axis flags respected in gate.

### 11. Test Coverage ✓

**Test 9** (per-axis headline): cx=0.61 commits, cy=0.17 withheld, cz=empty; whole-lever min would NOT commit.

**Test 10** (non-flagged byte-identical): every snapshot pins xyz==whole; byte-identical round-trip blob1==blob2.

**Test 11** (persistence): per-axis split confirmed, v4 constant pinned, blob parse flags at [7,8,9], folded axis off prior, withheld axes exactly at prior, restored flags held through refill, re-serialize byte-identical, v3-stamp rejected VersionMismatch.

### 12. Strict-Core ✓

- No heap post-init: `lever_committed_xyz[kMaxSourcesCap][3]` stack-resident.
- Bounds-checked axis access: guards `axis<0||axis>2` → return 0 (safe).
- No FP undefined: lever_arm safe 3-vector access; publish loop [0,1,2] in bounds.

---

## Verdict

APPROVE — Implementation correct, spec-compliant, thoroughly tested. Persistence v3→v4 clean. Per-axis gates independent. Byte-identical guarantee pinned. No strict-core violations. All 3 new test cases (9-11) pass with high assertion count (70252 + adapters 46/772).

