# sim — simulation rig (relaxed edge)

The validation backbone (D24): it is the only place **ground truth is known**, so the only place calibration *correctness* can be checked.

## Responsibilities (ISSUES.md Slice 14)
- Generate a parameterized GT trajectory (straight / turning / mixed regimes, configurable excitation).
- Spawn N synthetic `ISource`s with **planted** extrinsics, scale, time-offset, noise, and injectable outliers / dropout.
- Drive `ofc::Estimator::step()` deterministically and record `Result` streams.
- Provide scenario presets for the observability self-tests (pure-straight, pure-turning, ‖ω‖-varying).

Used by `tests/` for convergence, observability, and NEES/NIS checks.
