# Slice 13 CSV-Replay Code Review (commit ddcf436)

## Summary
Real-data ingestion scaffolding: hand-rolled CSV parser, GT-track reader, offline replay/eval harness, and CLI. Dep-free. Extensive tests verify round-trips, equivalence, GT eval, GPS correction, and parse robustness.

RESULT: Clean. No blocking issues. Ship as-is.

---

## Findings by Area

### CSV Parser Robustness

adapters/src/csv_source.cpp (lines 22-114)

- Comments/blanks/headers: # and // stripped; blank lines skipped; leading non-numeric header row tolerated. Correct.
- Delimiters: Comma OR whitespace; empty fields dropped. Correct.
- Numeric parsing: Both parse_double and parse_i64 catch exceptions, check used > 0, verify no trailing garbage, return false on failure. stod/stoll throws caught and mapped to Status::InvalidConfig. Correct.
- Covariance optional: Checked explicitly; absent = zero diagonal. Correct.
- Finite check: std::isfinite(d) rejects NaN/inf. Correct.
- Empty file: rows.size() < 2 returns Status::NoData. Test verified (line 505-507). Correct.
- Non-increasing timestamps: Checked after each row; returns InvalidConfig with line number. Correct.

No issues found.

---

### Quaternion Handling

adapters/src/csv_source.cpp (lines 119-140); test_csv_source.cpp (line 169-184)

- quat_to_mat3: Norm-squared threshold 1e-12 catches degenerate quats; isfinite(n2) catches NaN. Degenerate -> identity + ok=false -> parser error. Normalizes before matrix conversion. Correct.
- mat3_to_quat: Constructs from matrix, normalizes, canonical hemisphere (w >= 0). Correct.
- Round-trip test: Non-trivial rotation written and read back; query matches expectation to 1e-9. Correct.
- Equivalence tolerance 1e-6: Justified by window interpolation + 200 Hz CSV sampling. Tight but honest.

No issues found.

---

### Form Semantics

adapters/src/csv_source.cpp (L145-277); test_csv_source.cpp (L71-166)

- Absolute: Cumulative pose -> push_absolute(). Test: query([t0,t1]) = pose(t0)^-1 o pose(t1). Correct.
- Increment: Per-step relative; first row identity seed -> push_increment() for k>=1. Test: integration verified. Correct.
- Twist: Body v,w -> push_twist(). Test: constant twist = exp(xi*(t1-t0)). Correct.
- Covariance diagonal [trans;rot]: Applied correctly via diag6(). Correct.
- Form tag + ctor override: Tag sets unless force_form=true (ctor wins). Tested. Correct.

No issues found.

---

### NEES/NIS/Drift Math

adapters/src/replay_harness.cpp (L21-29, L129-131); test_csv_source.cpp (L332-346)

- NEES formula: e = se3::log(T_est^-1 o T_gt), nees = e^T (P_pp)^-1 e. EXACT match to test_validation.cpp. Correct.
- GT anchor: Anchored at gt_pose(frontier - window_s). Matches sim Rig. Correct.
- Drift metrics: trans_err_m and rot_err_rad correctly computed. Correct.
- NIS: Read from CorrectionDiag.last_nis on steps with corr_applied > 0. Correct.
- Warmup skip: NEES accumulation after warmup_steps (default 20). Correct.
- Hand-check NEES: Computed and verified to epsilon 1e-9. Correct.

No issues found.

---

### Equivalence Test

test_csv_source.cpp (L192-268)

- Setup: Analytic vs CSV path with same constant-twist scenario. Both use same SourceBuffer/Estimator.
- Comparison: For each (stamp, pose) in analytic, find match in CSV, check se3_close(a,b,1e-6).
- Timeline: 50 Hz steps 0.2-2.0s; CSV spans 2.1s at 5ms grid. Bracketing guaranteed.
- Tolerance: Interpolation error justified by SourceBuffer window + 200 Hz sampling ~ 1e-6m. Honest.
- Genuine equivalence: Real CSV == in-memory (same paths), not a tautology. Correct.

No issues found.

---

### Timeline Construction

replay_harness.cpp (L53-88)

- Regular grid: Tick at 1/tick_rate_hz over union of source spans. Correct.
- Fusion delay: GPS-fix frontier ticks = fix.stamp + delay_ns added to grid. No aliasing. Correct.
- Deduping: Sort + unique merges duplicates. Correct.
- Source span union: Min/Max over all sources. Correct.
- Empty timeline guard: Returns Status::NoData if no ticks. Correct.
- GPS submission: Before step(), submit fixes whose stamp <= frontier. Cadence matches e2e test. Correct.

No issues found.

---

### GpsCorrection Reuse

tools/ofc_replay/main.cpp (L322-328, L342-343); replay_harness.hpp (L90-92)

- Config wiring: GpsCorrection(man.gps_cfg) from manifest [gps] block. Correct.
- Ownership: Caller-owned; ptrs/refs in ReplayInputs. Lifetimes safe. Correct.
- Datum lazy-latch: First fix latches if datum_lat/lon/alt omitted. Contract allows. Correct.
- No core changes: GpsCorrection is existing adapter (Slice 11b); only used here. Correct.

No issues found.

---

### Manifest Parser

tools/ofc_replay/main.cpp (L136-223)

- Section routing: [global], [sensor.N], [gps], [gt], [replay] recognized. Unknown sections error. Correct.
- Extension consumption: csv= and form= consumed; other keys passed to ConfigLoader. Correct.
- GPS block: Keys parsed into GpsConfig; unknown keys error. Correct.
- GT block: Only csv= allowed. Correct.
- Replay block: tail_window_s, warmup_steps, out parsed. Correct.
- ConfigLoader reuse: cfg_out fed to ConfigLoader; Config validated normally. Correct.
- Source matching: Loop over source_entries, load sources[si], match to cfg.sensors[slot] with bounds check. Safe. Correct.

No issues found.

---

### Dep-Free / Relaxed-Edge / CLI Guarded

adapters/CMakeLists.txt (L62-69); adapters/src/*.cpp; tools/ofc_replay/main.cpp

- No strict-core changes: Only adapters/ + tools/ modified. Correct.
- Quat helpers adapter-local: Not exported to core. Correct.
- CLI guarded: add_executable(ofc_replay) in adapters block, only when OFC_BUILD_ADAPTERS=ON. Correct.
- Relaxed-edge compliance: std::string, std::vector, std::fstream, exceptions documented. Correct.
- Exception safety: stod/stoll throws caught at source; parse funcs return false, no escape. Correct.

No issues found.

---

## Verdict

No bugs. No risks. No nits. Ship as-is.

Correct for ship:
- CSV parser: comments, blanks, headers, delimiters, exception handling all solid.
- Quaternion: normalization, degeneracy detection, canonical hemisphere, round-trip accurate.
- Form semantics: absolute/increment/twist fed correctly; covariance applied.
- NEES/NIS/drift: exact match to sim conventions; GT anchoring + 6-DOF math sound.
- Equivalence: genuine (same paths), tolerance honest.
- Timeline: regular grid + GPS-fix ticks, no aliasing, fusion-delay aware.
- GpsCorrection: wiring correct, lifetimes safe.
- Manifest: routing correct, ConfigLoader reused, bounds-safe.
- Dep-free/guarded: no core changes, adapter-local helpers, CLI guarded.
- Exception handling: errors caught early, Status returned.

Test coverage verified:
- Round-trip x3 forms.
- Quaternion round-trip.
- Equivalence (50+ steps).
- GT eval (drift/NEES hand-checked).
- GPS via CSV (drift reduction).
- Parse robustness (all edge cases).

Safe.
