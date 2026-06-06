# OdometryFusionCalibration

A lightweight framework for **odometry fusion** and **online calibration**.

Fuses N odometry sources into one robust motion estimate (a geometric-median-driven error-state integrator) and self-calibrates their relative extrinsics, scale, and time-offsets online — from odometry disagreement alone, no external reference required.

## Documentation
- [`DESIGN.md`](./DESIGN.md) — architecture specification.
- [`DECISIONS.md`](./DECISIONS.md) — decision log (chosen / rejected / why).
- [`CONFIG.md`](./CONFIG.md) — full configuration reference.
- [`ISSUES.md`](./ISSUES.md) — implementation roadmap (tracer-bullet slices).

## Layout
```
include/ofc/core/   public core headers   — STRICT (no heap post-init, no exceptions, bounded WCET, C++14/AUTOSAR)
src/core/           core implementation    — STRICT
adapters/           middleware integrations — relaxed (ROS, config parsers, threading, file persistence)
sim/                simulation rig         — relaxed
tests/              unit / observability / consistency — relaxed
```
The **core** is platform-, middleware-, and allocation-discipline-strict; everything that never reaches the ECU (adapters, sim, tests) is relaxed. See [`DESIGN.md`](./DESIGN.md) §9.

## Status
Scaffold. Interfaces and structure are in place; implementation follows [`ISSUES.md`](./ISSUES.md), starting at Slice 0. **Not yet building** — requires Eigen and the core implementation.

## Build (planned)
```sh
cmake -B build -DOFC_BUILD_TESTS=ON
cmake --build build
```
