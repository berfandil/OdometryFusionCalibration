// ofc/core/smoother.hpp — per-sensor fixed-lag RTS twist smoother (Slice 10, D18).
//
// THE USER INSIGHT (D18, load-bearing). Calibration tolerates latency, so it can run at
// a DEEPER frontier (now − delay_fast − L) than the causal fusion frontier. There it is
// fed body twists that have been smoothed by a TWO-SIDED (forward + backward) pass over a
// fixed lag window L — a fixed-lag Rauch-Tung-Striebel (RTS) smoother. Two-sided ≈
// ZERO-PHASE: unlike a one-sided causal filter (EMA / forward-KF) it does NOT lag a
// time-varying signal, so the calibration histogram peaks SHARPEN (lower input variance)
// WITHOUT shifting (no calibration bias). Fusion stays causal + byte-identical; only the
// CALIBRATION feed sees the smoothed twist.
//
// THE MODEL (per source, D18 "CV twist ESKF, random-walk accel"). An augmented
// constant-velocity (CV) error-state KF on the body twist:
//   state  x = [ w ; a ]   ∈ ℝ¹²    (w = body twist ℝ⁶, a = twist-rate "accel" ℝ⁶)
//   transition (Δt)   F = [[ I, Δt·I ], [ 0, I ]]      (w advances by a·Δt; CV)
//   process noise     random-walk on a, strength kf_process_noise = q:
//       Q(Δt) = q · [[ Δt³/3·I,  Δt²/2·I ],
//                    [ Δt²/2·I,  Δt·I    ]]            (continuous white-accel discretized)
//   measurement       z = w (the measured twist log(B_corr)/dt), H = [ I, 0 ], R = r·I.
// The "twist" is the body twist log(delta)/dt consistent with the rest of the codebase
// (eskf.cpp twist readout, median fused_omega/fused_trans). w[0..2] = linear v, w[3..5] =
// angular ω, matching Vec6 [v; omega].
//
// FIXED-LAG RTS (the two-sided pass). We keep a ring of the last L+1 filtered states +
// their one-step prediction + the F used, where L = lag_steps. Each push():
//   (1) FORWARD: KF predict (F, Q) then measurement update with z — append to the ring.
//   (2) BACKWARD: once the ring holds L+1 entries, run the RTS recursion from the newest
//       back to the OLDEST entry:
//         C_k       = P_k Fᵀ_{k+1} (P^pred_{k+1})⁻¹           (smoother gain)
//         x^s_k     = x_k + C_k (x^s_{k+1} − x^pred_{k+1})
//         P^s_k     = P_k + C_k (P^s_{k+1} − P^pred_{k+1}) C_kᵀ
//       and EMIT the oldest entry's smoothed twist w (= the sample L steps in the PAST).
//   Before the ring fills (fewer than L+1 pushes for a slot) the slot is NOT ready: ready()
//   returns false and smoothed() falls back to the latest filtered twist (pass-through), so
//   a caller that ignores ready() degrades to the causal estimate rather than garbage.
//
// LAG FORMULA. lag_steps is supplied by the estimator as round(calib_lag_s · tick_rate_hz),
// CLAMPED to [1, kMaxLag] (a compile-time cap sizing the fixed rings — no post-configure
// heap). L = 0 would be a degenerate "no smoothing" pass-through; we clamp UP to 1.
//
// STRICT CORE: all storage is fixed-capacity, sized at configure() (max_sources × (L+1)
// ring of 12-vectors + 12×12 covariances); push()/smoothed() allocate nothing; bounded
// loops (the RTS recursion is ≤ L+1 iterations); Status returns; no exceptions; double.
#ifndef OFC_CORE_SMOOTHER_HPP
#define OFC_CORE_SMOOTHER_HPP

#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

namespace ofc {

class TwistSmoother {
public:
    // Compile-time caps sizing the fixed rings (strict core: no post-configure heap).
    // kMaxSources matches the Result / calibrator arrays; kMaxLag bounds L+1 ring depth.
    static constexpr int kMaxSources = 32;
    static constexpr int kMaxLag     = 64;   // max lag_steps; ring depth is lag_steps+1

    TwistSmoother() = default;

    // Configure the smoother once (preallocates the per-source rings). After this no heap
    // is allocated.
    //   max_sources  — number of source slots used (1..kMaxSources).
    //   lag_steps    — the fixed lag L in STEPS (clamped to [1, kMaxLag]); the smoothed
    //                  output trails the latest push by L steps (the deeper frontier).
    //   process_noise— random-walk accel strength q (kf_process_noise); must be > 0.
    //   meas_noise   — measurement noise r on the twist (default 1.0; > 0). Only the RATIO
    //                  q/r shapes the smoothing; r is a convenience knob.
    // Returns OutOfRange for an out-of-range max_sources / non-positive noise.
    Status configure(int max_sources, int lag_steps,
                     Scalar process_noise, Scalar meas_noise = Scalar(1));

    // Drop all per-source state (keeps the configuration + storage). Every slot returns to
    // not-ready; the next push() restarts the filter.
    void reset();

    // The clamped lag in steps (>= 1) the smoother actually runs at. 0 if unconfigured.
    int lag_steps() const { return configured_ ? lag_ : 0; }

    // Push one measured body twist for source `slot` over an interval of length `dt` (s).
    //   slot       — source index in [0, max_sources).
    //   twist_meas — the measured body twist z = log(B_corr)/dt ([v; omega]).
    //   dt         — the step interval in seconds (> 0; a non-positive dt is clamped to a
    //                tiny positive value, matching eskf.cpp's dt guard).
    // Runs the forward KF step + (once the ring is full) the backward RTS pass, advancing
    // the slot's emitted (lag-L) smoothed twist. Heap-free; bounded.
    // Returns OutOfRange for a bad slot, NotInitialized if unconfigured, else Ok.
    Status push(int slot, const Vec6& twist_meas, Scalar dt);

    // True once source `slot` has received at least lag_steps+1 pushes — i.e. its emitted
    // smoothed twist is the genuine TWO-SIDED estimate of the sample L steps in the past.
    // Before that the slot is warming up (smoothed() is a causal pass-through fallback).
    bool ready(int slot) const;

    // The smoothed body twist for source `slot` at the EMITTED (lag-L-old) timestamp. Once
    // ready(), this is the RTS two-sided estimate (lower variance, zero phase). Before
    // ready() it falls back to the latest filtered twist (causal pass-through). Zero for an
    // unconfigured / never-pushed / out-of-range slot.
    Vec6 smoothed(int slot) const;

    // The refined covariance (12×12, [w; a]) of the emitted smoothed state — the RTS
    // posterior P^s of the lag-L-old entry. Identity for a not-ready / bad slot. Cheap (it
    // is already maintained by the backward pass); exposed for the calibrators' vote weight.
    Mat12 smoothed_cov(int slot) const;

    bool configured() const { return configured_; }

private:
    // Build the CV transition F(dt) and process-noise Q(dt) for a step of length dt.
    static Mat12 transition(Scalar dt);
    static Mat12 process_noise(Scalar dt, Scalar q);

    bool   configured_  = false;
    int    max_sources_ = 0;
    int    lag_         = 1;        // clamped lag L (>= 1); ring depth = lag_ + 1
    Scalar q_           = Scalar(1);
    Scalar r_           = Scalar(1);

    // Per-source fixed-lag ring. depth = lag_ + 1 entries, indexed by a head/count. Each
    // entry stores the FORWARD filtered posterior (x, P), the one-step PREDICTION that
    // produced it (x_pred, P_pred), and the transition F used — all needed by the backward
    // RTS recursion. Sized at the compile-time cap so configure() is the sole allocator.
    static constexpr int kRing = kMaxLag + 1;
    struct Slot {
        Vec6   emitted     = Vec6::Zero();          // last emitted (lag-L) smoothed twist
        Mat12  emitted_cov = Mat12::Identity();     // its RTS covariance
        bool   has_filt    = false;                 // at least one forward step taken
        int    count       = 0;                     // pushes seen (caps the ring fill)
        // Ring storage (logical order: 0 = oldest, count-1 = newest).
        Eigen::Matrix<Scalar, 12, 1> x[kRing];      // filtered posterior mean
        Mat12 P[kRing];                             // filtered posterior cov
        Eigen::Matrix<Scalar, 12, 1> xp[kRing];     // one-step prediction mean
        Mat12 Pp[kRing];                            // one-step prediction cov
        Mat12 F[kRing];                             // transition that produced entry k
    };
    Slot slots_[kMaxSources];
};

} // namespace ofc
#endif // OFC_CORE_SMOOTHER_HPP
