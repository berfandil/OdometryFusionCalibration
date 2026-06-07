// ofc_adapters/gps_correction.cpp — see gps_correction.hpp for the model.
//
// WGS-84 geodetic -> ECEF -> ENU -> odom, then a dim=3 position Measurement (residual + H + R).
// RELAXED EDGE (adapter): std/<cmath> fine.
#include "ofc_adapters/gps_correction.hpp"

#include "ofc/core/lie.hpp"   // so3::hat (skew of the lever arm)

#include <cmath>

namespace ofc {
namespace adapters {

namespace {

// WGS-84 ellipsoid constants.
constexpr double kWgs84A  = 6378137.0;                 // semi-major axis (m)
constexpr double kWgs84F  = 1.0 / 298.257223563;       // flattening
const     double kWgs84E2 = kWgs84F * (2.0 - kWgs84F); // first eccentricity squared
const     double kDeg2Rad = std::atan2(0.0, -1.0) / 180.0; // pi/180

// Geodetic (lat,lon in radians; h in m) -> ECEF (m). Standard WGS-84 forward transform.
Vec3 geodetic_to_ecef(double lat, double lon, double h) {
    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);
    const double sin_lon = std::sin(lon);
    const double cos_lon = std::cos(lon);
    const double N = kWgs84A / std::sqrt(1.0 - kWgs84E2 * sin_lat * sin_lat);
    Vec3 ecef;
    ecef.x() = (N + h) * cos_lat * cos_lon;
    ecef.y() = (N + h) * cos_lat * sin_lon;
    ecef.z() = (N * (1.0 - kWgs84E2) + h) * sin_lat;
    return ecef;
}

// The ECEF->ENU rotation at a datum (lat0=phi, lon0=lambda, both radians). Rows are the local
// East / North / Up axes expressed in ECEF:
//   [ -sinλ        cosλ        0    ]   (East)
//   [ -sinφcosλ  -sinφsinλ   cosφ  ]   (North)
//   [  cosφcosλ   cosφsinλ   sinφ  ]   (Up)
Mat3 ecef_to_enu_rotation(double lat0, double lon0) {
    const double sin_phi = std::sin(lat0);
    const double cos_phi = std::cos(lat0);
    const double sin_lam = std::sin(lon0);
    const double cos_lam = std::cos(lon0);
    Mat3 R;
    R << -sin_lam,            cos_lam,            0.0,
         -sin_phi * cos_lam, -sin_phi * sin_lam,  cos_phi,
          cos_phi * cos_lam,  cos_phi * sin_lam,  sin_phi;
    return R;
}

} // namespace

GpsCorrection::GpsCorrection(const GpsConfig& cfg) : cfg_(cfg) {
    // An explicit datum is latched up front so geodetic->ENU works before the first submit_fix
    // and last_enu() is meaningful immediately. A lazy datum (!has_datum) is latched on the first
    // valid submit_fix instead.
    if (cfg_.has_datum) {
        latch_datum(cfg_.datum_lat_deg, cfg_.datum_lon_deg, cfg_.datum_alt_m);
    }
}

void GpsCorrection::latch_datum(double lat_deg, double lon_deg, double alt_m) const {
    const double lat = lat_deg * kDeg2Rad;
    const double lon = lon_deg * kDeg2Rad;
    datum_ecef_     = geodetic_to_ecef(lat, lon, alt_m);
    R_ecef_to_enu_  = ecef_to_enu_rotation(lat, lon);
    have_datum_     = true;
}

Vec3 GpsCorrection::fix_to_enu(const GpsFix& fix) const {
    const double lat = fix.lat_deg * kDeg2Rad;
    const double lon = fix.lon_deg * kDeg2Rad;
    const Vec3 ecef  = geodetic_to_ecef(lat, lon, fix.alt_m);
    return R_ecef_to_enu_ * (ecef - datum_ecef_);
}

void GpsCorrection::submit_fix(const GpsFix& fix) {
    if (!fix.valid) return;   // a no-fix / bad-fix sentinel is ignored

    // Lazy datum latch: the FIRST valid fix sets the tangent-plane origin when no datum was
    // configured (mirrors the sim ref's lazy anchor).
    if (!have_datum_) {
        latch_datum(fix.lat_deg, fix.lon_deg, fix.alt_m);
    }

    // Latest-wins: a newer fix replaces an unconsumed one. Refresh last_enu_ for diagnostics.
    pending_      = fix;
    have_pending_ = true;
    last_enu_     = fix_to_enu(fix);
}

bool GpsCorrection::evaluate(const State& x, Measurement& out) const {
    if (!have_pending_ || !have_datum_) return false;   // no fresh fix this step

    // --- geodetic -> ENU -> odom : the measured ANTENNA position z in the odom frame ----------
    const Vec3 p_enu = fix_to_enu(pending_);
    last_enu_        = p_enu;
    const Vec3 z     = cfg_.odom_from_enu.R * p_enu + cfg_.odom_from_enu.t;

    // --- predicted antenna position h(x) = t + R*lever -----------------------------------------
    const Vec3& l = cfg_.lever_arm_base;
    const Vec3  h = x.pose.t + x.pose.R * l;

    // --- residual = z (-) h(x) -----------------------------------------------------------------
    out.dim = 3;
    out.residual.setZero();
    out.residual.head<3>() = z - h;

    // --- Jacobian H = [ R | -R*[l]x | 0_3x6 ]  (3x12) ------------------------------------------
    // Right-error tangent: h(eta) = t + R*rho + R*(I+[theta]x)*l => dh/drho = R, dh/dtheta = -R*[l]x
    // (since [theta]x * l = -[l]x * theta). With l = 0 this collapses to the sim ref's [R | 0 | 0].
    out.H.setZero();
    out.H.block<3, 3>(0, 0) = x.pose.R;
    out.H.block<3, 3>(0, 3) = -x.pose.R * so3::hat(l);

    // --- noise R_odom = R_align*(cov_enu + floor*I3)*R_align^T ----------------------------------
    const Mat3 cov_enu_floored = pending_.cov_enu + cfg_.cov_floor_m2 * Mat3::Identity();
    const Mat3 R_align         = cfg_.odom_from_enu.R;
    out.R.setZero();
    out.R.block<3, 3>(0, 0) = R_align * cov_enu_floored * R_align.transpose();

    // Stamp policy: the FIX's own acquisition time (a real sensor stamp), not x.stamp. See header.
    out.stamp = pending_.stamp;

    have_pending_ = false;   // consumed: one emit per submitted fix
    return true;
}

bool GpsCorrection::has_datum() const { return have_datum_; }

Vec3 GpsCorrection::last_enu() const { return last_enu_; }

} // namespace adapters
} // namespace ofc
