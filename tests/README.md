# tests — unit / observability / consistency (relaxed edge)

Layered validation (D24). Determinism (caller-pumped, lock-free core) makes golden regression cheap.

## Layers
1. **Unit** — per block: Lie ops vs analytic, Weiszfeld median, histogram mode + aging, hand-eye split solve, ‖ω‖ cross-correlation, ESKF update.
2. **Observability self-tests** — assert each DOF converges *only* in its regime:
   - yaw/pitch + scale converge under straight motion, frozen otherwise;
   - roll + xyz lever-arm converge under turning, frozen under pure-straight;
   - time-offset converges only when ‖ω‖ varies.
   These regression-guard the observability spine.
3. **NEES / NIS consistency** — Monte-Carlo over sim runs: the published covariance must be neither over- nor under-confident.
4. **Golden regression** — recorded/sim input replayed → byte-stable `Result` output.

Driven by the `sim/` rig; registered with CTest via the root `OFC_BUILD_TESTS` option.
