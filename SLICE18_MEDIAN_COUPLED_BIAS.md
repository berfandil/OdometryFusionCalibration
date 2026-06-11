# Slice 18 — Median-coupled multi-source bias states (11b Option B)

**Goal**: the fused trajectory's heading on real data is structurally worse than its best input, and the urban12 coast residual is the visible symptom. Learn each source's slowly-varying body-twist bias online (observable through GPS), de-bias every source BEFORE the median — so the consensus stops tracking the middle drifter and coast drift collapses toward the best-source floor.

---

## 1. Evidence (2026-06-11, `tools/diag_heading_drift.py` + `urban12_c4k08_out.csv`)

Per-source heading drift vs GT, urban12 (truth from `global_pose`):

| source | full-drive heading err | rate | 522 s coast drift |
|---|---|---|---|
| wheel (diff yaw) | −51.8° | −75.2°/h | −5.2° |
| wheel+FOG | **+4.6°** | +6.7°/h | **+0.76°** |
| wheel+IMU | −16.3° (max 29.7°) | −23.6°/h | **+23.9°** |

The FUSED heading (current recommended config) wanders 5–30° all drive, spiking 61° — **worse than the best input everywhere, not just coasts**. Mechanism: per window, the geometric median tracks the MIDDLE yaw increment, and the middle drifter here is wheel+IMU. Slice-9 reliability cannot fix this: the middle source has the LOWEST residual-vs-consensus by construction, and `SourceHealth.bias` is an unsigned magnitude. The median structurally cannot express "trust FOG's yaw, wheel's forward speed".

De-biasing each source before the median CAN: wheel's −75°/h is a near-constant yaw-rate bias (one random-walk state nails it); IMU's is slowly varying (tracked by `bias_process_noise`); FOG learns ≈0. The median of de-biased sources then sits near the FOG-quality consensus, and a coast inherits learned biases instead of raw drift.

## 2. Design (the deferred 11b Option B, now unblocked)

The ISSUES Slice-11b DESIGN NOTE + UPDATE govern; key derivation pre-validated in spike 5a and re-validated by the median fix:

- **State**: `[pose(6); twist(6); bias_1(6) … bias_k(6)]` for the k sources with `bias_states=true`, k ≤ `kMaxBiasSources = 4` (compile-time; augmented dim ≤ 36). Option A is the k=1 sole-driver special case — its tests and behavior must stay green; the Mat18 machinery generalizes to the capped fixed size (implementer may unify A onto the generalized path ONLY if every existing Option-A test stays green unmodified).
- **De-bias point**: each biased source's SOURCE-FRAME delta, before frame-align: `B_i' = B_i ∘ exp(−b_i·dt)`; then the usual `A_i = X_i∘B_i'∘X_i⁻¹` → median. (Bias is a body-twist bias; the Option-A de-bias mechanic per source.)
- **Coupling Jacobian — CORRECTED 2026-06-11 (the spike-5a scalar form is FALSIFIED post-median-fix)**: spike 5a's `J = −dt·ω_i·Ad(X_i)` was FD-verified against the OLD pinning median (degenerate ω = (1,0,…) — an input passthrough, for which a scalar weight is exact). Against the TRUE interior Weiszfeld median it is 38–98% wrong (FD-falsified, `tests/test_multi_bias.cpp`): the median's input sensitivity carries a **radial projector** — perturbing an input ALONG its offset from the median doesn't move the median — which no scalar can express. The EXACT first-order influence (implicit function theorem at the converged fixed point, FD-verified worst 1.5%, typical <1%):
  - `Ω_i = M⁻¹·u_i·P_i`, `P_i = I₆ − ξ_i(Wξ_i)ᵀ/d_i²`, `M = Σ_j u_j P_j`, `u_i = w_i/d_i`, `ξ_i = [t_i−t_m; log(R_mᵀR_i)]`, `W = diag(λI₃, I₃)`; `Σ_i Ω_i = I` exactly (ω_i·I is its isotropic approximation).
  - `J_{pose,b_i} = blkdiag(R_mᵀ, I)·Ω_i·blkdiag(R_{A_i}, I)·(−dt·Ad(X_i))` — the sign/frame factor `−dt·Ad(X_i)` from spike 5a survives (n=1 pin exact).
  - Needs only the per-source `split_distance` d_i + fusion weights (estimator-side, no median API change) + one 6×6 solve per step. Edge cases subsumed: absent source → u_i=0 → Ω_i=0 (frozen); sole driver → Ω=I (Option A exact); weight-dominant vertex → Ω_i→I. **n=2**: the solver interpolates with FIXED weights → influence is `(w_i/Σw)·I` (the d-based form is wrong there).
  - Honest caveat: near a weight-dominant vertex the median is non-smooth — any first-order J has a shrinking linearization radius; Ω_i degrades gracefully (→ I) and the EKF tolerates it.
- **Predict**: median over de-biased aligned deltas → consensus drives the ESKF exactly as today; augmented `F` adds the per-source pose↔bias blocks above + per-bias random walk; `Q_bias,i = bias_process_noise_i·dt·I₆`.
- **Update**: `update_aug` generalized to the active augmented dim (the per-`n` χ² gate, Joseph, robust kappa, C4 rot-suppression all carry over unchanged — they act on the pose rows).
- **Observability**: biases observable ONLY with an absolute ref (self-test: no-ref ⇒ `bias_observable` ≈ 0 for all, planted bias NOT recovered); with GPS, a planted constant bias on one source is recovered while a clean co-source's bias stays ≈ 0 (separation test — the ω-weighted coupling is what distinguishes them).
- **Config**: `Config::multi_bias_enabled` (default **false** — the existing multi-`bias_states` `InvalidConfig` guard and every Option-A behavior stay byte-identical). Loader `[global] multi_bias_enabled`. Config-hash. `CalibSnapshot.bias`/`bias_observable` extend per-source (the fields exist; today only the single bias slot fills them).
- **Calibrator-consensus note (ISSUES Slice-14 contamination item)**: de-biased deltas make the consensus cleaner, which the calibrators consume — strictly an improvement, but the calib observe path must receive the SAME de-biased deltas fusion used (consistency; document the single de-bias site).
- **Strict core**: fixed `Mat<36>` worst-case allocated once; bounded; no heap.

**REJECTED**:
- Lightweight signed-EMA de-bias outside the filter — no covariance bookkeeping; the filter would double-count or fight it; violates the honest-Σ discipline that the whole Slice-14 saga enforces.
- Per-DOF median weighting — changes the fusion primitive itself (D3); far more invasive than removing the bias upstream.
- Marking FOG as a "heading-priority" source — config hack, doesn't generalize, leaves the bias in every other consumer.

## 3. Acceptance

Unit (TDD):
1. Jacobian FD-check test: the analytic `Ω_i`-based block matches finite differences of the de-biased-median pipeline at randomized states — DONE (`def0cee`, `tests/test_multi_bias.cpp`: exact block <5% worst, scalar-form falsification pinned, n=1/n=2 regime pins, FD two-epsilon self-consistency).
2. Multi-source separation: 3-source sim + GPS, planted constant twist-bias on source A only → A's bias recovered (<10% err), B/C's stay ≈ 0; fused drift with the fix ≪ without.
3. No-ref observability self-test: planted bias NOT recovered, `bias_observable` ≈ 0 (never weaken).
4. Coast scenario: GPS-rich learn phase, then a GPS-denied stretch → de-biased fused heading drift ≪ biased baseline (the urban12 shape, in sim).
5. Influence edge cases: source dropped from a window (Ω=0, no coupling), single participant (Ω=I ≡ Option A), n=2 fixed-weight interpolation ((w_i/Σw)·I), VZ-guard coincident vertex, weight-dominant vertex (Ω→I, no blow-up).
6. Option A regression: every existing 11b test green unmodified; `multi_bias_enabled=false` byte-identical (exact pin).
7. Loader key; config-hash flip; NEES guard untouched (cov-cal band [4.0, 5.6] — run it with the flag ON too: the augmented filter must not corrupt the 12-DOF marginals).

Gate: `scripts/dev.ps1 -Task test` green.

Real-data (orchestrator):
- urban12, recommended urban config + `multi_bias_enabled` + `bias_states` on wheel/wheel+IMU (FOG too — its bias should learn ≈0): mid-drive transient max (currently ~300–400 m) **< 100 m**, fused heading error trajectory materially below the current 5–30° band, tail stays ~2 m, NEES sane.
- urban07/17: no regression (drift + NEES).

## 4. Status

- [ ] Implemented (TDD, gate green, committed)
- [ ] Reviewed (`reviews/slice-18-findings.md`) + findings fixed
- [ ] Real-data validation table filled in
- [ ] Docs updated (CONFIG/DECISIONS/DESIGN/ISSUES) — orchestrator
