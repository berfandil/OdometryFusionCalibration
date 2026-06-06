# src/core — strict core implementation

**Profile (AUTOSAR C++14):** no heap after `init()`, no exceptions, no RTTI, bounded WCET, fixed-capacity from `Config`, `double`/no-fast-math. See [`../../DESIGN.md`](../../DESIGN.md) §9.

## Modules (by roadmap slice)
| File | Slice | Responsibility |
|---|---|---|
| `estimator.cpp` | facade | `init/step/latest`, validate, persistence entry points |
| `lie.cpp` | 0 | in-house SO(3)/SE(3) exp/log/adjoint/compose, split metric |
| `buffer.cpp` | 1 | per-source ring buffer + `delta(t0,t1)` query |
| `median.cpp` | 2 | Weiszfeld split-metric geometric median |
| `eskf.cpp` | 2 | predict on median, adaptive Q, absolute-ref update, tip |
| `histogram.cpp` | 4 | fixed-bin vote/age/mode primitive (shared) |
| `timesync.cpp` | 5 | ‖ω‖ cross-correlation → offset histogram |
| `calibration.cpp` | 6–8 | Phase 1 (so(3)+scale), Phase 2 (hand-eye), commit/feedback |
| `lifecycle.cpp` | 3 | INIT→WARMUP→DEGRADED→NOMINAL |

All allocation is owned by `Estimator::Impl`, sized from `Config` at `init()`.
