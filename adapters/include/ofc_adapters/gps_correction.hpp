// ofc_adapters/gps_correction.hpp — the PRODUCTION GPS absolute-reference adapter (D22,
// Slice 11b/13). Turns real GPS fixes (WGS-84 geodetic lat/lon/alt + ENU position covariance)
// into the core's dim=3 POSITION correction (ofc::ICorrection) — the SAME ESKF update path the
// sim ref (ofc::sim::SyntheticAbsoluteRef) exercises — so a real GPS stream removes the
// fused-odometry drift. The estimator Mahalanobis-gates + applies the measurement at the frontier
// (see estimator.cpp's correction loop; the per-`n` chi2_gate scales the n=3 base for this dim).
//
// RELAXED EDGE. This is an adapter, never shipped to the ECU: std / heap / exceptions are fine.
// It depends only on the PUBLIC ofc/core headers (via ofc_core) + Eigen.
//
// MEASUREMENT MODEL (position fix, dim=3; mirrors the sim ref, + an antenna lever arm). The
// adapter converts geodetic -> ECEF -> ENU about a tangent-plane datum, then ENU -> the
// estimator's odom frame via a rigid alignment, yielding the measured ANTENNA position z. The
// predicted antenna position is h(x) = x.pose.t + x.pose.R * lever_arm_base (the base position
// plus the rotated lever arm). With the ESKF's right-error tangent T_true = T o exp([rho;theta]):
//       h(eta) = t + R*rho + R*(I + [theta]x)*l   (first order)
//   =>  H = [ R | -R*[l]x | 0_3x6 ]   (3x12),   residual = z (-) h(x).
// With lever = 0 this reduces to the sim ref's H = [R | 0 | 0]. R_odom rotates the ENU position
// covariance into odom: R_odom = R_align*(cov_enu + cov_floor_m2*I3)*R_align^T.
//
// SCOPE: position-only (dim=3). A full-pose (dim=6) or orientation fix is a noted extension —
// same machinery, a different h/H — NOT required here.
#ifndef OFC_ADAPTERS_GPS_CORRECTION_HPP
#define OFC_ADAPTERS_GPS_CORRECTION_HPP

#include "ofc/core/absolute_ref.hpp"
#include "ofc/core/types.hpp"

#include <vector>

namespace ofc {
namespace adapters {

// One GPS fix: geodetic position + ENU position covariance + stamp. (Relaxed edge.)
struct GpsFix {
    double    lat_deg = 0.0;     // WGS-84 latitude  (degrees)
    double    lon_deg = 0.0;     // WGS-84 longitude (degrees)
    double    alt_m   = 0.0;     // ellipsoidal height (m)
    Mat3      cov_enu = Mat3::Identity();  // ENU position covariance (m^2). diag(sx^2,sy^2,su^2) typical.
    Timestamp stamp   = 0;       // fix time (ns), same clock as Estimator::step()
    bool      valid   = true;    // a no-fix / bad-fix sentinel -> evaluate() yields nothing
};

struct GpsConfig {
    // Datum (ENU tangent-plane origin). If has_datum is false, the FIRST valid fix is latched
    // as the datum (origin = first fix; common "local origin" pattern, mirrors the sim ref's
    // lazy anchor). If true, datum_lat/lon/alt are used.
    bool   has_datum     = false;
    double datum_lat_deg = 0.0;
    double datum_lon_deg = 0.0;
    double datum_alt_m   = 0.0;

    // Rigid alignment ENU -> odom frame (the estimator's pose frame). z_odom = odom_from_enu o p_enu.
    // Default identity = "GPS ENU already aligned to odom" (the integrator supplies the real
    // alignment, e.g. from an initial heading). FOOTGUN (like the sim ref's window_s): a WRONG
    // rotation here biases EVERY residual by a constant rotation of the position error — the drift
    // test would still partially improve, masking the misconfig. The caller owns this alignment.
    SE3 odom_from_enu = SE3();   // identity

    // GPS antenna position in the BASE frame (lever arm), meters. Default zero (antenna at the
    // base origin). A nonzero lever makes the measured antenna position depend on the base
    // orientation; the Jacobian carries the -R*[l]x rotation-coupling block.
    Vec3 lever_arm_base = Vec3::Zero();

    // Optional R inflation / floor (m^2 added to the ENU cov diagonal before rotating into odom)
    // for robustness against an over-optimistic receiver covariance. Default 0.
    double cov_floor_m2 = 0.0;

    // INNOVATION-ADAPTIVE ROBUST R (off by default; GPS_R_NEES_SWEEP.md follow-up). A scalar
    // cov_floor cannot be cross-drive consistent: the honest GPS measurement R varies ~2 -> 22 m^2
    // across KAIST drives, and a multipath drive (urban12) carries a population of gross outliers.
    // When adaptive_r is true, evaluate() sets R from the ROBUST scale (MAD) of the last
    // `adaptive_window` innovations (z - h) per odom axis instead of the cov_enu+floor path:
    //   R_odom = diag( max( (1.4826*MAD_i)^2, adaptive_r_floor_m2 ) ).
    // MAD tracks the BULK spread (outlier-immune), so R stays small on a clean drive (consistent
    // NIS) and gross multipath fixes show a LARGE NIS -> the estimator's Mahalanobis gate rejects
    // them, instead of a uniform cov_floor inflation that over-distrusts the GOOD fixes too. The
    // ring fills causally (R for fix k uses fixes [k-window, k-1]); until `adaptive_min_samples`
    // innovations are seen it falls back to the cov_enu+floor path. Relaxed-edge: heap ring.
    bool   adaptive_r           = false;
    int    adaptive_window      = 60;     // innovations in the robust-scale window
    int    adaptive_min_samples = 10;     // fall back to cov_floor until this many seen
    double adaptive_r_floor_m2  = 0.25;   // R floor (m^2) — the best-case receiver noise
};

// A production ICorrection that emits a dim=3 GPS position fix at the current frontier state.
//
// USAGE: construct with a GpsConfig, register via Estimator::add_correction(&gps), and call
// submit_fix(...) from the GPS subscriber/integrator thread context (single-threaded into the
// core: submit between steps, like the sim ref returns a fix per evaluate()). The adapter holds
// the most recent UNCONSUMED valid fix; evaluate() (called inside step()) emits it once.
class GpsCorrection : public ICorrection {
public:
    explicit GpsCorrection(const GpsConfig& cfg);

    // Push the latest GPS fix (called by the integrator/subscriber). The adapter holds the most
    // recent UNCONSUMED valid fix; evaluate() emits it once. A newer fix REPLACES an unconsumed
    // one (latest-wins; GPS at a few Hz vs the step rate). An invalid (valid=false) fix is
    // ignored. The FIRST valid fix latches the lazy datum when !cfg.has_datum.
    void submit_fix(const GpsFix& fix);

    // ICorrection: if a fresh (unconsumed, valid) fix is pending, build the dim=3 position
    // Measurement at the current frontier state x and mark the fix consumed (return true); else
    // return false (no fix available this step). Stamp policy: out.stamp = the FIX's stamp (the
    // real sensor time), NOT x.stamp — a deviation from the sim ref (which stamps x.stamp because
    // its "fix" is sampled AT the frontier). A real GPS fix has its own acquisition time; the
    // estimator's per-step correction loop overwrites x.stamp = t1 for evaluate(), and uses
    // out.stamp only as a diagnostic / for the measurement record.
    bool evaluate(const State& x, Measurement& out) const override;

    // Diagnostics.
    bool has_datum() const;            // datum latched yet?
    Vec3 last_enu() const;             // last computed ENU position (m), for diagnostics/tests

private:
    // Convert a geodetic fix to the ENU position about the (latched) datum. Requires has_datum_.
    Vec3 fix_to_enu(const GpsFix& fix) const;
    // Latch the datum from a fix's geodetic position (and precompute the ECEF origin + rotation).
    void latch_datum(double lat_deg, double lon_deg, double alt_m) const;

    GpsConfig cfg_;

    // mutable lazy datum + pending-fix state: evaluate() is const (ICorrection contract) but the
    // datum / pending fix are established on the fly (relaxed-edge, mirrors the sim ref's mutable
    // lazy anchor).
    mutable bool      have_datum_   = false;
    mutable Vec3      datum_ecef_   = Vec3::Zero();   // ECEF of the datum origin (m)
    mutable Mat3      R_ecef_to_enu_ = Mat3::Identity(); // ECEF->ENU rotation at the datum

    mutable bool      have_pending_ = false;          // an unconsumed valid fix is waiting
    mutable GpsFix    pending_;                       // the most recent unconsumed valid fix
    mutable Vec3      last_enu_     = Vec3::Zero();    // last computed ENU (diagnostics)

    // Adaptive-R: a causal ring of the last innovations (odom-frame z - h) + a helper that returns
    // the robust per-axis variance (MAD^2) over the ring, or false if fewer than min_samples.
    mutable std::vector<Vec3> innov_ring_;            // most-recent at the back; trimmed to window
    bool adaptive_var(Vec3& var_out) const;
};

} // namespace adapters
} // namespace ofc
#endif // OFC_ADAPTERS_GPS_CORRECTION_HPP
