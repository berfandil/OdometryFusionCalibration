// ofc_adapters/config_loader.cpp — see config_loader.hpp for the format + ownership contract.
//
// RELAXED EDGE: std::string / std::vector / heap / exceptions fine. No exception escapes a
// public method (parse-time std::stod throws are caught -> InvalidConfig).
#include "ofc_adapters/config_loader.hpp"

#include "ofc/core/types.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ofc {
namespace adapters {

namespace {

// --- small string helpers ----------------------------------------------------------------
std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

bool parse_bool(const std::string& v, bool& out) {
    const std::string s = lower(trim(v));
    if (s == "true" || s == "1" || s == "on" || s == "yes")  { out = true;  return true; }
    if (s == "false" || s == "0" || s == "off" || s == "no") { out = false; return true; }
    return false;
}

bool parse_double(const std::string& v, Scalar& out) {
    try {
        std::size_t used = 0;
        const double d = std::stod(trim(v), &used);
        if (used == 0) return false;
        // reject trailing garbage after the number (e.g. "1.0abc")
        const std::string rest = trim(trim(v).substr(used));
        if (!rest.empty()) return false;
        out = static_cast<Scalar>(d);
        return true;
    } catch (...) {
        return false;
    }
}

// Parse a base-10 int. On failure returns false; if the value was a well-formed number that
// overflows `long`, *overflow (when non-null) is set true so the caller can report it distinctly.
bool parse_int(const std::string& v, long& out, bool* overflow = nullptr) {
    if (overflow != nullptr) *overflow = false;
    try {
        std::size_t used = 0;
        const long n = std::stol(trim(v), &used, 10);
        if (used == 0) return false;
        const std::string rest = trim(trim(v).substr(used));
        if (!rest.empty()) return false;
        out = n;
        return true;
    } catch (const std::out_of_range&) {
        if (overflow != nullptr) *overflow = true;
        return false;
    } catch (...) {
        return false;
    }
}

// 6 whitespace-separated numbers: yaw pitch roll x y z.
bool parse_extrinsic(const std::string& v, SE3& out) {
    std::istringstream is(v);
    Scalar a[6];
    for (int i = 0; i < 6; ++i) {
        double d;
        if (!(is >> d)) return false;
        a[i] = static_cast<Scalar>(d);
    }
    double extra;
    if (is >> extra) return false;   // too many numbers

    const Scalar yaw = a[0], pitch = a[1], roll = a[2];
    Mat3 Rz, Ry, Rx;
    Rz << std::cos(yaw), -std::sin(yaw), 0,
          std::sin(yaw),  std::cos(yaw), 0,
          0, 0, 1;
    Ry << std::cos(pitch), 0, std::sin(pitch),
          0, 1, 0,
          -std::sin(pitch), 0, std::cos(pitch);
    Rx << 1, 0, 0,
          0, std::cos(roll), -std::sin(roll),
          0, std::sin(roll),  std::cos(roll);
    out.R = Rz * Ry * Rx;
    out.t = Vec3(a[3], a[4], a[5]);
    return true;
}

// q_floor: 1 number (applied to ALL 6 axes) OR 6 numbers (per-axis, [trans; rot] order).
bool parse_q_floor(const std::string& v, Scalar out[6]) {
    std::istringstream is(v);
    double vals[6];
    int n = 0;
    double d;
    while (n < 6 && (is >> d)) vals[n++] = d;
    double extra;
    if (is >> extra) return false;                       // > 6 numbers
    if (n == 1) { for (int i = 0; i < 6; ++i) out[i] = static_cast<Scalar>(vals[0]); return true; }
    if (n == 6) { for (int i = 0; i < 6; ++i) out[i] = static_cast<Scalar>(vals[i]); return true; }
    return false;                                        // must be exactly 1 or 6 numbers
}

const char* kKnobDoc =
    "ofc config (key=value / INI). '#' or ';' comment. Sections: [global], [sensor.N].\n"
    "[global]\n"
    "  max_sources=<int>            reference_sensor_id=<int>\n"
    "  window_s=<f>   fusion_delay_s=<f>   tick_rate_hz=<f>   calib_lag_s=<f>\n"
    "  cold_start=reference_only|median_from_start\n"
    "  vote_weight=one|confidence|rotation|combo\n"
    "  commit_concentration=<f>     commit_min_votes=<int>   commit_drop=<f>\n"
    "  straight_omega_max=<f>  straight_trans_min=<f>  turn_omega_min=<f>\n"
    "  timesync_enabled=<bool>\n"
    "  rot3d_enabled=<bool>     (turn-regime FULL rotation extrinsic, Slice 17; default off)\n"
    "  joint_lever_scale=<bool> (turn-regime joint lever+scale 4-unknown LS, Slice 17b; default off)\n"
    "  multi_bias_enabled=<bool> (median-coupled multi-source bias states, Slice 18; default off)\n"
    "  multi_bias_cov0=<f>      (multi-bias per-DOF prior variance seed, > 0; default 0.04)\n"
    "  subbin_centroid=<bool>   (centroid sub-bin readout, all 5 calib histograms)\n"
    "  adaptive_q=<bool>   q_scale=<f>   q_floor=<f | 6 f's [trans;rot]>\n"
    "  mahalanobis_chi2=<f>   correction_robust_kappa=<f>   correction_rot_suppress_kappa=<f>  (0=off)\n"
    "[sensor.N]\n"
    "  id=<int>   prior_extrinsic=<yaw pitch roll x y z>   prior_scale=<f>\n"
    "  prior_time_offset_s=<f>   weight_prior=<f>\n"
    "  per_sensor_kf=<bool>   scale_calib=<bool>   is_reference=<bool>\n"
    "  bias_states=<bool>     bias_process_noise=<f | 6 f's [v;omega]>   (per-source body-twist\n"
    "      bias, Slice 11b/18; 6-number per-DOF form needs multi_bias_enabled — a DOF rate of\n"
    "      exactly 0 PINS that bias DOF at zero)\n";

} // namespace

const char* ConfigLoader::knob_doc() { return kKnobDoc; }

Status ConfigLoader::parse(const std::string& text) {
    cfg_ = Config{};
    sensors_.clear();
    error_.clear();

    auto fail = [&](const std::string& msg, int line_no, const std::string& line) -> Status {
        std::ostringstream os;
        os << "line " << line_no << ": " << msg << "  [" << line << "]";
        error_ = os.str();
        return Status::InvalidConfig;
    };

    // Parse an int field, distinguishing a well-formed-but-out-of-range value (NIT) from generic
    // junk in the message. Returns true on success (out set); on failure leaves out untouched and
    // sets *why to the appropriate diagnostic for the caller to pass to fail().
    auto parse_int_field = [](const std::string& v, long& out, std::string& why) -> bool {
        bool overflow = false;
        if (parse_int(v, out, &overflow)) return true;
        why = overflow ? "int value out of range" : "expected int";
        return false;
    };

    enum class Scope { Global, Sensor };
    Scope scope = Scope::Global;
    int   cur_sensor = -1;   // index into sensors_

    std::istringstream in(text);
    std::string raw;
    int line_no = 0;
    while (std::getline(in, raw)) {
        ++line_no;
        // strip comments
        std::string line = raw;
        const std::size_t hash = line.find_first_of("#;");
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = trim(line);
        if (line.empty()) continue;

        // section header
        if (line.front() == '[') {
            if (line.back() != ']') return fail("malformed section header", line_no, raw);
            const std::string sec = lower(trim(line.substr(1, line.size() - 2)));
            if (sec == "global" || sec.empty()) {
                scope = Scope::Global;
            } else if (sec.compare(0, 7, "sensor.") == 0) {
                long idx;
                if (!parse_int(sec.substr(7), idx) || idx < 0)
                    return fail("bad sensor section index", line_no, raw);
                // sensors are stored densely in first-seen order; the index must match the
                // next slot (so [sensor.0] then [sensor.1] ...). Re-entering an existing slot
                // is allowed (append more keys); skipping ahead is an error.
                if (idx == static_cast<long>(sensors_.size())) {
                    sensors_.push_back(SensorConfig{});
                    cur_sensor = static_cast<int>(idx);
                } else if (idx < static_cast<long>(sensors_.size())) {
                    cur_sensor = static_cast<int>(idx);   // re-open
                } else {
                    return fail("sensor index skips a slot (expected sequential)", line_no, raw);
                }
                scope = Scope::Sensor;
            } else {
                return fail("unknown section", line_no, raw);
            }
            continue;
        }

        // key = value
        const std::size_t eq = line.find('=');
        if (eq == std::string::npos) return fail("expected key = value", line_no, raw);
        const std::string key = lower(trim(line.substr(0, eq)));
        const std::string val = trim(line.substr(eq + 1));
        if (key.empty()) return fail("empty key", line_no, raw);

        if (scope == Scope::Global) {
            long   iv;
            Scalar dv;
            bool   bv;
            std::string why;
            if (key == "max_sources") {
                if (!parse_int_field(val, iv, why)) return fail(why, line_no, raw);
                cfg_.max_sources = static_cast<int>(iv);
            } else if (key == "reference_sensor_id") {
                if (!parse_int_field(val, iv, why)) return fail(why, line_no, raw);
                if (iv < 0) return fail("expected non-negative int", line_no, raw);
                cfg_.reference_sensor_id = static_cast<SourceId>(iv);
            } else if (key == "window_s") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.window_s = dv;
            } else if (key == "fusion_delay_s") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.fusion_delay_s = dv;
            } else if (key == "tick_rate_hz") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.tick_rate_hz = dv;
            } else if (key == "calib_lag_s") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.calib_lag_s = dv;
            } else if (key == "cold_start") {
                const std::string lv = lower(val);
                if (lv == "reference_only")       cfg_.cold_start = ColdStart::ReferenceOnly;
                else if (lv == "median_from_start") cfg_.cold_start = ColdStart::MedianFromStart;
                else return fail("cold_start: reference_only|median_from_start", line_no, raw);
            } else if (key == "vote_weight") {
                const std::string lv = lower(val);
                if (lv == "one")             cfg_.vote_weight = VoteWeight::One;
                else if (lv == "confidence") cfg_.vote_weight = VoteWeight::Confidence;
                else if (lv == "rotation")   cfg_.vote_weight = VoteWeight::Rotation;
                else if (lv == "combo")      cfg_.vote_weight = VoteWeight::Combo;
                else return fail("vote_weight: one|confidence|rotation|combo", line_no, raw);
            } else if (key == "commit_concentration") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.commit_concentration = dv;
            } else if (key == "commit_min_votes") {
                if (!parse_int_field(val, iv, why)) return fail(why, line_no, raw);
                cfg_.commit_min_votes = static_cast<int>(iv);
            } else if (key == "commit_drop") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.commit_drop = dv;
            } else if (key == "straight_omega_max") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.straight_omega_max = dv;
            } else if (key == "straight_trans_min") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.straight_trans_min = dv;
            } else if (key == "turn_omega_min") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.turn_omega_min = dv;
            } else if (key == "timesync_enabled") {
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                cfg_.timesync_enabled = bv;
            } else if (key == "rot3d_enabled") {
                // Slice 17: turn-regime FULL rotation extrinsic (axis-correspondence
                // hand-eye). Default off = byte-identical legacy behavior.
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                cfg_.rot3d_enabled = bv;
            } else if (key == "joint_lever_scale") {
                // Slice 17b: turn-regime JOINT lever+scale (4-unknown hand-eye LS).
                // Default off = byte-identical 3-unknown lever numerics.
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                cfg_.joint_lever_scale = bv;
            } else if (key == "multi_bias_enabled") {
                // Slice 18: median-coupled multi-source bias states (11b Option B).
                // Default off = byte-identical legacy behavior (Option A unchanged).
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                cfg_.multi_bias_enabled = bv;
            } else if (key == "multi_bias_cov0") {
                // Slice 18 review/B2 (MAJOR-3): the per-DOF multi-bias prior variance seed
                // (> 0; validate()). Rig-dependent stability knob; default 0.04.
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.multi_bias_cov0 = dv;
            } else if (key == "subbin_centroid") {
                // Slice 16: one switch, applied to ALL FIVE calibration histograms.
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                cfg_.so3_hist.subbin_centroid    = bv;
                cfg_.roll_hist.subbin_centroid   = bv;
                cfg_.xyz_hist.subbin_centroid    = bv;
                cfg_.scale_hist.subbin_centroid  = bv;
                cfg_.offset_hist.subbin_centroid = bv;
            } else if (key == "adaptive_q") {
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                cfg_.adaptive_q = bv;
            } else if (key == "q_scale") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.q_scale = dv;
            } else if (key == "q_floor") {
                if (!parse_q_floor(val, cfg_.q_floor))
                    return fail("q_floor: 1 number (all axes) or 6 numbers [trans;rot]",
                                line_no, raw);
            } else if (key == "mahalanobis_chi2") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.mahalanobis_chi2 = dv;
            } else if (key == "correction_robust_kappa") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.correction_robust_kappa = dv;
            } else if (key == "correction_rot_suppress_kappa") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                cfg_.correction_rot_suppress_kappa = dv;
            } else {
                return fail("unknown global key '" + key + "'", line_no, raw);
            }
        } else {  // Scope::Sensor
            SensorConfig& sc = sensors_[static_cast<std::size_t>(cur_sensor)];
            long   iv;
            Scalar dv;
            bool   bv;
            std::string why;
            if (key == "id") {
                if (!parse_int_field(val, iv, why)) return fail(why, line_no, raw);
                if (iv < 0) return fail("expected non-negative int", line_no, raw);
                sc.id = static_cast<SourceId>(iv);
            } else if (key == "prior_extrinsic") {
                if (!parse_extrinsic(val, sc.prior_extrinsic))
                    return fail("prior_extrinsic: yaw pitch roll x y z (6 numbers)", line_no, raw);
            } else if (key == "prior_scale") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                sc.prior_scale = dv;
            } else if (key == "prior_time_offset_s") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                sc.prior_time_offset_s = dv;
            } else if (key == "weight_prior") {
                if (!parse_double(val, dv)) return fail("expected number", line_no, raw);
                sc.weight_prior = dv;
            } else if (key == "per_sensor_kf") {
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                sc.per_sensor_kf = bv;
            } else if (key == "scale_calib") {
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                sc.scale_calib = bv;
            } else if (key == "is_reference") {
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                sc.is_reference = bv;
            } else if (key == "bias_states") {
                // Slice 11b/18: this source carries a body-twist bias state (Option A when
                // [global] multi_bias_enabled is off — single source only; Option B when on).
                if (!parse_bool(val, bv)) return fail("expected bool", line_no, raw);
                sc.bias_states = bv;
            } else if (key == "bias_process_noise") {
                // Slice 11b/18: the bias random-walk process-noise rate(s) (>= 0 per DOF;
                // validate()). ONE number = uniform across the 6 DOFs (the legacy form);
                // SIX numbers = per-DOF [v; omega] order (Slice 18 review/B2 — a DOF rate
                // of exactly 0 PINS that bias DOF under multi_bias_enabled). Reuses the
                // q_floor 1-or-6 parse.
                Scalar pn[6];
                if (!parse_q_floor(val, pn))
                    return fail("expected 1 or 6 numbers", line_no, raw);
                for (int d = 0; d < 6; ++d) sc.bias_process_noise[d] = pn[d];
            } else {
                return fail("unknown sensor key '" + key + "'", line_no, raw);
            }
        }
    }

    // Cross-section semantic checks the core validate() does not yet make (estimator.cpp's
    // validate() has a "TODO: per-sensor checks" and only range-checks reference_sensor_id):
    //
    // (a) Duplicate sensor id. The [sensor.N] section index is just a dense storage slot; the
    // bound id is `id=`. Two sections can set the same id, which parses + validates Ok yet
    // silently mis-binds at add_source/restore time. Reject it loudly here.
    for (std::size_t i = 0; i < sensors_.size(); ++i) {
        for (std::size_t j = i + 1; j < sensors_.size(); ++j) {
            if (sensors_[i].id == sensors_[j].id) {
                std::ostringstream os;
                os << "duplicate sensor id " << static_cast<long>(sensors_[i].id)
                   << " (sections [sensor." << i << "] and [sensor." << j << "])";
                error_ = os.str();
                return Status::InvalidConfig;
            }
        }
    }

    // (b) is_reference (per-sensor) vs reference_sensor_id (global) consistency. We treat
    // reference_sensor_id as the source of truth: AT MOST ONE sensor may carry is_reference=true,
    // and if one does its id must equal reference_sensor_id — otherwise the two disagree silently.
    {
        int  ref_flag_count = 0;
        long flagged_id = -1;
        for (const SensorConfig& sc : sensors_) {
            if (sc.is_reference) { ++ref_flag_count; flagged_id = static_cast<long>(sc.id); }
        }
        if (ref_flag_count > 1) {
            error_ = "more than one sensor flagged is_reference=true";
            return Status::InvalidConfig;
        }
        if (ref_flag_count == 1 &&
            static_cast<SourceId>(flagged_id) != cfg_.reference_sensor_id) {
            std::ostringstream os;
            os << "is_reference=true on sensor id " << flagged_id
               << " disagrees with reference_sensor_id=" << static_cast<long>(cfg_.reference_sensor_id);
            error_ = os.str();
            return Status::InvalidConfig;
        }
    }

    // Point Config at the owned sensor storage (the ownership contract).
    cfg_.sensors      = sensors_.empty() ? nullptr : sensors_.data();
    cfg_.sensor_count = static_cast<int>(sensors_.size());

    // Run the core validator — the adapter's job is to BUILD a struct the core accepts.
    const Status vs = validate(cfg_);
    if (!ok(vs)) {
        std::ostringstream os;
        os << "validate() rejected the parsed config (status " << static_cast<int>(vs) << ")";
        error_ = os.str();
        return vs;
    }
    return Status::Ok;
}

Status ConfigLoader::parse_file(const std::string& path) {
    std::string text;
    try {
        std::ifstream f(path, std::ios::binary);
        if (!f) { error_ = "could not open file: " + path; return Status::NoData; }
        std::ostringstream ss;
        ss << f.rdbuf();
        text = ss.str();
    } catch (...) {
        error_ = "IO error reading file: " + path;
        return Status::CorruptData;
    }
    return parse(text);
}

Status load_config(const std::string& text, Config& out_cfg,
                   std::vector<SensorConfig>& out_sensors, std::string* out_error) {
    ConfigLoader loader;
    const Status st = loader.parse(text);
    out_sensors = loader.sensors();                 // copy the owned storage to the caller's
    out_cfg = loader.config();
    out_cfg.sensors = out_sensors.empty() ? nullptr : out_sensors.data();  // re-point at caller's
    out_cfg.sensor_count = static_cast<int>(out_sensors.size());
    if (out_error != nullptr) *out_error = loader.error();
    return st;
}

} // namespace adapters
} // namespace ofc
