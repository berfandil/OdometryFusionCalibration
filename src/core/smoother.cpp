// ofc/core/smoother.cpp — per-sensor fixed-lag RTS twist smoother (Slice 10, D18).
//
// See smoother.hpp for the model + math. In brief, per source we run an augmented
// constant-velocity (CV) KF on the body twist:
//   x = [ w ; a ]      (twist w ∈ ℝ⁶, twist-rate a ∈ ℝ⁶)
//   F(dt) = [[ I, dt·I ], [ 0, I ]]
//   Q(dt) = q · [[ dt³/3·I, dt²/2·I ], [ dt²/2·I, dt·I ]]   (white-accel discretization)
//   z = w,  H = [ I, 0 ],  R = r·I
// FORWARD: predict then measurement-update each push; store (x,P) + the prediction
// (x_pred,P_pred) + F in a fixed-lag ring of depth L+1. BACKWARD: once the ring is full,
// the RTS recursion
//   C_k   = P_k Fᵀ_{k+1} (P_pred_{k+1})⁻¹
//   x^s_k = x_k + C_k (x^s_{k+1} − x_pred_{k+1})
//   P^s_k = P_k + C_k (P^s_{k+1} − P_pred_{k+1}) C_kᵀ
// runs from the newest entry back to the OLDEST, whose smoothed twist (the lag-L sample)
// is EMITTED. The pass is two-sided -> near zero-phase (uses future measurements to
// de-lag a time-varying signal), unlike a forward-only causal filter.
//
// STRICT CORE: fixed-size Eigen only (no heap after configure()); bounded loops (the
// recursion is ≤ L+1 iterations); no exceptions; covariances symmetrized; double.
#include "ofc/core/smoother.hpp"

#include <Eigen/LU>          // inverse() for the small 12x12 RTS gain solve

#include <algorithm>

namespace ofc {

namespace {
constexpr Scalar kMinDt = Scalar(1e-9);   // dt guard (matches eskf.cpp)

using Vec12 = Eigen::Matrix<Scalar, 12, 1>;

Mat12 symmetrize(const Mat12& M) { return Scalar(0.5) * (M + M.transpose()); }
} // namespace

Mat12 TwistSmoother::transition(Scalar dt) {
    Mat12 F = Mat12::Identity();
    F.block<6, 6>(0, 6) = dt * Mat6::Identity();   // w_{k} = w_{k-1} + dt * a_{k-1}
    return F;
}

Mat12 TwistSmoother::process_noise(Scalar dt, Scalar q) {
    // Continuous white-acceleration model discretized over dt (the textbook CV/CA Q):
    //   Q = q * [[ dt^3/3 I, dt^2/2 I ], [ dt^2/2 I, dt I ]].
    const Scalar dt2 = dt * dt;
    const Scalar dt3 = dt2 * dt;
    const Scalar a = q * dt3 / Scalar(3);
    const Scalar b = q * dt2 / Scalar(2);
    const Scalar c = q * dt;
    Mat12 Q = Mat12::Zero();
    Q.block<6, 6>(0, 0) = a * Mat6::Identity();
    Q.block<6, 6>(0, 6) = b * Mat6::Identity();
    Q.block<6, 6>(6, 0) = b * Mat6::Identity();
    Q.block<6, 6>(6, 6) = c * Mat6::Identity();
    return Q;
}

Status TwistSmoother::configure(int max_sources, int lag_steps,
                                Scalar process_noise, Scalar meas_noise) {
    if (max_sources < 1 || max_sources > kMaxSources) return Status::OutOfRange;
    if (!(process_noise > Scalar(0)))                 return Status::OutOfRange;
    if (!(meas_noise > Scalar(0)))                    return Status::OutOfRange;

    max_sources_ = max_sources;
    // Clamp the lag to [1, kMaxLag]: L=0 is a degenerate no-smoothing pass (clamp UP);
    // beyond the compile-time ring cap clamp DOWN (the rings are sized at kMaxLag).
    lag_ = std::max(1, std::min(kMaxLag, lag_steps));
    q_   = process_noise;
    r_   = meas_noise;
    reset();
    configured_ = true;
    return Status::Ok;
}

void TwistSmoother::reset() {
    for (int i = 0; i < kMaxSources; ++i) {
        Slot& sl = slots_[i];
        sl.emitted     = Vec6::Zero();
        sl.emitted_cov = Mat12::Identity();
        sl.has_filt    = false;
        sl.count       = 0;
        // Ring buffers are overwritten on push before use; no need to clear them all.
    }
}

Status TwistSmoother::push(int slot, const Vec6& twist_meas, Scalar dt) {
    if (!configured_)                       return Status::NotInitialized;
    if (slot < 0 || slot >= max_sources_)   return Status::OutOfRange;
    const Scalar dt_eff = (dt > kMinDt) ? dt : kMinDt;

    Slot& sl = slots_[slot];
    const int depth = lag_ + 1;             // ring capacity for this configuration

    // --- FORWARD KF step -----------------------------------------------------
    const Mat12 F = transition(dt_eff);
    const Mat12 Q = process_noise(dt_eff, q_);

    Vec12 x_prev;
    Mat12 P_prev;
    if (!sl.has_filt) {
        // Initialize from the first measurement: twist = z, accel = 0; a generous prior on
        // the unobserved accel so the first few steps adapt quickly to a ramp.
        x_prev.setZero();
        x_prev.head<6>() = twist_meas;
        P_prev = Mat12::Identity();
        P_prev.block<6, 6>(0, 0) = r_ * Mat6::Identity();     // twist prior ~ meas noise
        P_prev.block<6, 6>(6, 6) = Scalar(1e3) * Mat6::Identity();  // accel: weak prior
    } else {
        // Newest filtered entry is at logical index count-1 (capped at depth-1).
        const int newest = std::min(sl.count, depth) - 1;
        x_prev = sl.x[newest];
        P_prev = sl.P[newest];
    }

    // Predict.
    const Vec12 x_pred = F * x_prev;
    const Mat12 P_pred = symmetrize(F * P_prev * F.transpose() + Q);

    // Measurement update: z = w = H x, H = [I6, 0]. Innovation on the top-6 block only.
    //   S = H P_pred H^T + R = P_pred[0:6,0:6] + r I   (6x6)
    //   K = P_pred H^T S^-1                            (12x6)
    //   x = x_pred + K (z - H x_pred)
    //   P = (I - K H) P_pred  (then symmetrize)
    const Mat6 S = P_pred.block<6, 6>(0, 0) + r_ * Mat6::Identity();
    // P_pred H^T = the first 6 COLUMNS of P_pred (12x6).
    const Eigen::Matrix<Scalar, 12, 6> PHt = P_pred.leftCols<6>();
    const Eigen::Matrix<Scalar, 12, 6> K   = PHt * S.inverse();   // 12x6, S small & SPD
    const Vec6  innov = twist_meas - x_pred.head<6>();
    const Vec12 x_upd = x_pred + K * innov;
    // Joseph-free but PSD-safe enough: (I - K H) P_pred. K H is 12x12 with only the first 6
    // columns populated (H selects the top block).
    Mat12 KH = Mat12::Zero();
    KH.leftCols<6>() = K;                      // (K H) places K into the first 6 columns
    const Mat12 P_upd = symmetrize((Mat12::Identity() - KH) * P_pred);

    // --- append to the fixed-lag ring (logical order: 0 oldest .. newest) -----
    if (sl.count < depth) {
        // Ring not yet full: append at the next slot.
        const int idx = sl.count;
        sl.x[idx]  = x_upd;
        sl.P[idx]  = P_upd;
        sl.xp[idx] = x_pred;
        sl.Pp[idx] = P_pred;
        sl.F[idx]  = F;
    } else {
        // Ring full: shift everything one step older (drop the oldest), append at the end.
        // depth <= kRing so this bounded shift is cheap + heap-free (strict core).
        for (int k = 0; k < depth - 1; ++k) {
            sl.x[k]  = sl.x[k + 1];
            sl.P[k]  = sl.P[k + 1];
            sl.xp[k] = sl.xp[k + 1];
            sl.Pp[k] = sl.Pp[k + 1];
            sl.F[k]  = sl.F[k + 1];
        }
        const int last = depth - 1;
        sl.x[last]  = x_upd;
        sl.P[last]  = P_upd;
        sl.xp[last] = x_pred;
        sl.Pp[last] = P_pred;
        sl.F[last]  = F;
    }
    sl.has_filt = true;
    if (sl.count < depth) ++sl.count;   // saturate count at the ring depth

    // --- BACKWARD RTS pass (only once the ring is full) -----------------------
    // The smoothed estimate of the OLDEST entry (logical 0 = the lag-L sample) is what we
    // emit. The recursion needs the entry AHEAD's smoothed (x^s, P^s) and the entry AHEAD's
    // PREDICTION (x_pred, P_pred) + the F that produced the entry ahead.
    if (sl.count >= depth) {
        const int last = depth - 1;
        Vec12 xs = sl.x[last];          // newest entry: smoothed == filtered
        Mat12 Ps = sl.P[last];
        for (int k = last - 1; k >= 0; --k) {
            // Smoother gain C_k = P_k F_{k+1}^T (P_pred_{k+1})^{-1}.
            const Mat12& Pk    = sl.P[k];
            const Mat12& Fk1   = sl.F[k + 1];        // transition that produced entry k+1
            const Mat12& Pp1   = sl.Pp[k + 1];       // prediction cov of entry k+1
            const Vec12& xp1   = sl.xp[k + 1];       // prediction mean of entry k+1
            const Mat12  C     = Pk * Fk1.transpose() * Pp1.inverse();
            const Vec12  xs_k  = sl.x[k] + C * (xs - xp1);
            const Mat12  Ps_k  = symmetrize(Pk + C * (Ps - Pp1) * C.transpose());
            xs = xs_k;
            Ps = Ps_k;
        }
        // xs / Ps now hold the smoothed state of the OLDEST entry (the lag-L sample).
        sl.emitted     = xs.head<6>();
        sl.emitted_cov = Ps;
    } else {
        // Not yet full: causal pass-through fallback (latest filtered twist).
        sl.emitted     = x_upd.head<6>();
        sl.emitted_cov = P_upd;
    }
    return Status::Ok;
}

bool TwistSmoother::ready(int slot) const {
    if (!configured_ || slot < 0 || slot >= max_sources_) return false;
    return slots_[slot].count >= lag_ + 1;
}

Vec6 TwistSmoother::smoothed(int slot) const {
    if (!configured_ || slot < 0 || slot >= max_sources_) return Vec6::Zero();
    return slots_[slot].emitted;
}

Mat12 TwistSmoother::smoothed_cov(int slot) const {
    if (!configured_ || slot < 0 || slot >= max_sources_) return Mat12::Identity();
    if (!slots_[slot].has_filt)                            return Mat12::Identity();
    return slots_[slot].emitted_cov;
}

} // namespace ofc
