// tools/ofc_replay/main.cpp — the standalone real-data replay CLI (Slice 13, real-data
// scaffolding). The USER runs this on real CSVs: it reads a MANIFEST (the config_loader INI text
// format, EXTENDED with the CSV file paths + the GPS datum/lever/alignment + the GT file path),
// loads the CsvSources / GPS fixes / GT track, runs the ReplayHarness, and writes a results CSV +
// a human-readable summary.
//
// RELAXED EDGE (a tool): std iostream / fstream / heap / exceptions all fine. Built only when
// OFC_BUILD_ADAPTERS is ON (guarded in CMake), so the default build is unaffected.
//
// =====================================================================================
// MANIFEST FORMAT (a SUPERSET of the config_loader format):
//
//   [global]                     # any config_loader [global] knob (window_s, tick_rate_hz, ...)
//   tick_rate_hz = 50
//   fusion_delay_s = 0.05
//   window_s = 0.10
//   mahalanobis_chi2 = 9.0        # (any Config knob the loader supports)
//   max_sources = 4
//   reference_sensor_id = 0
//
//   [sensor.0]                    # any config_loader [sensor.N] knob (id, prior_extrinsic, ...)
//   id = 0
//   is_reference = true
//   csv = data/odom0.csv          # <-- EXTENSION: the CSV file for THIS source
//   form = increment              # <-- EXTENSION: absolute|increment|twist (else the file tag)
//
//   [sensor.1]
//   id = 1
//   csv = data/odom1.csv
//   form = twist
//
//   [gps]                         # <-- EXTENSION (optional): the GPS absolute ref
//   csv = data/gps.csv            # rows: t_ns, lat_deg, lon_deg, alt_m [, sx2, sy2, su2]
//   datum_lat_deg = 47.0          # omit the 3 datum_* to lazily latch the first fix
//   datum_lon_deg = 8.0
//   datum_alt_m = 400.0
//   odom_from_enu_yaw = 0.0       # ENU->odom yaw (rad); pitch/roll optional
//   lever_x = 0.0  lever_y = 0.0  lever_z = 0.0   # antenna lever in the base frame (m)
//   cov_floor_m2 = 0.0
//
//   [gt]                          # <-- EXTENSION (optional): the GT track (eval only)
//   csv = data/gt.csv             # rows: t_ns, x, y, z, qw, qx, qy, qz
//
//   [replay]                      # <-- EXTENSION (optional): harness knobs
//   tail_window_s = 1.0
//   warmup_steps = 20
//   local_batch_len = 0           # GT-anchored relative-pose-error window (num fused-gt records);
//                                 # 0 = off. Same value across recordings -> length-fair drift.
//   out = results.csv             # output results-CSV path (default: replay_results.csv)
//
// GPS CSV schema (one fix per row; `#`/`//`/blank skipped; whitespace OR comma delimited):
//   t_ns, lat_deg, lon_deg, alt_m [, var_e, var_n, var_u]   (ENU position variances m^2; absent
//   -> a 1 m^2 isotropic default).
// =====================================================================================
#include "ofc_adapters/config_loader.hpp"
#include "ofc_adapters/csv_source.hpp"
#include "ofc_adapters/gps_correction.hpp"
#include "ofc_adapters/replay_harness.hpp"

#include "ofc/core/config.hpp"
#include "ofc/core/lie.hpp"
#include "ofc/core/status.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace ofc;
using namespace ofc::adapters;

namespace {

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
std::string strip_comment(const std::string& line) {
    std::size_t cut = line.size();
    const std::size_t h = line.find('#');  if (h != std::string::npos && h < cut) cut = h;
    const std::size_t s = line.find("//"); if (s != std::string::npos && s < cut) cut = s;
    return line.substr(0, cut);
}

bool to_double(const std::string& v, double& out) {
    try { std::size_t u = 0; out = std::stod(v, &u); return u != 0 && trim(v.substr(u)).empty(); }
    catch (...) { return false; }
}
bool to_i64(const std::string& v, long long& out) {
    try { std::size_t u = 0; out = std::stoll(v, &u, 10); return u != 0 && trim(v.substr(u)).empty(); }
    catch (...) { return false; }
}

std::vector<std::string> split_fields(const std::string& line) {
    std::vector<std::string> out; std::string cur;
    auto flush = [&]() { const std::string t = trim(cur); if (!t.empty()) out.push_back(t); cur.clear(); };
    for (char c : line) {
        if (c == ',' || std::isspace(static_cast<unsigned char>(c))) flush();
        else cur.push_back(c);
    }
    flush();
    return out;
}

// One per-source manifest entry (the CSV path + an optional explicit form).
struct SourceEntry {
    SourceId    id        = 0;
    std::string csv;
    OdomForm    form      = OdomForm::Increment;
    bool        has_form  = false;
};

// The parsed manifest: a config text (only the config_loader-known lines), per-source CSV paths,
// the optional GPS block, the optional GT path, and the replay knobs.
struct Manifest {
    std::string                config_text;     // re-fed to ConfigLoader (sans our extensions)
    std::map<int, SourceEntry> source_entries;  // keyed by [sensor.N] slot index
    bool        has_gps = false;  GpsConfig gps_cfg;  std::string gps_csv;
    bool        has_gt  = false;  std::string gt_csv;
    Scalar      tail_window_s   = Scalar(1.0);
    int         warmup_steps    = 20;
    int         local_batch_len = 0;
    std::string out_path        = "replay_results.csv";
};

// Parse the manifest text. EXTENSION keys (csv/form inside [sensor.N]; the whole [gps]/[gt]/
// [replay] sections) are consumed HERE; everything else is passed through to `config_text` so the
// existing ConfigLoader validates the core Config. Returns "" on success, else an error string.
std::string parse_manifest(const std::string& text, Manifest& m) {
    std::ostringstream cfg_out;
    std::istringstream in(text);
    std::string raw;
    int line_no = 0;
    enum class Sec { Global, Sensor, Gps, Gt, Replay };
    Sec sec = Sec::Global;
    int cur_sensor = -1;

    auto err = [&](const std::string& msg) {
        std::ostringstream os; os << "manifest line " << line_no << ": " << msg; return os.str();
    };

    while (std::getline(in, raw)) {
        ++line_no;
        const std::string line = trim(strip_comment(raw));
        if (line.empty()) continue;

        if (line.front() == '[') {
            if (line.back() != ']') return err("malformed section header");
            const std::string s = lower(trim(line.substr(1, line.size() - 2)));
            if (s == "gps")        { sec = Sec::Gps;    m.has_gps = true; continue; }
            if (s == "gt")         { sec = Sec::Gt;     m.has_gt  = true; continue; }
            if (s == "replay")     { sec = Sec::Replay; continue; }
            if (s == "global" || s.empty()) { sec = Sec::Global; cfg_out << line << "\n"; continue; }
            if (s.compare(0, 7, "sensor.") == 0) {
                long long idx;
                if (!to_i64(s.substr(7), idx) || idx < 0) return err("bad sensor section index");
                sec = Sec::Sensor; cur_sensor = static_cast<int>(idx);
                m.source_entries[cur_sensor];   // ensure the entry exists
                cfg_out << line << "\n";        // the section header goes to the config too
                continue;
            }
            return err("unknown section '" + s + "'");
        }

        const std::size_t eq = line.find('=');
        if (eq == std::string::npos) return err("expected key = value");
        const std::string key = lower(trim(line.substr(0, eq)));
        const std::string val = trim(line.substr(eq + 1));

        if (sec == Sec::Sensor) {
            SourceEntry& e = m.source_entries[cur_sensor];
            if (key == "csv") { e.csv = val; continue; }
            if (key == "form") {
                const std::string lv = lower(val);
                if (lv == "absolute")       { e.form = OdomForm::AbsolutePose; e.has_form = true; }
                else if (lv == "increment") { e.form = OdomForm::Increment;    e.has_form = true; }
                else if (lv == "twist")     { e.form = OdomForm::Twist;        e.has_form = true; }
                else return err("form must be absolute|increment|twist");
                continue;
            }
            // capture the bound id (for matching the loaded SensorConfig later)
            if (key == "id") { long long iv; if (to_i64(val, iv) && iv >= 0) e.id = static_cast<SourceId>(iv); }
            cfg_out << line << "\n";   // any other sensor key is a real Config knob
            continue;
        }
        if (sec == Sec::Global) { cfg_out << line << "\n"; continue; }
        if (sec == Sec::Gps) {
            double d;
            if (key == "csv")            { m.gps_csv = val; }
            else if (key == "datum_lat_deg" && to_double(val, d)) { m.gps_cfg.has_datum = true; m.gps_cfg.datum_lat_deg = d; }
            else if (key == "datum_lon_deg" && to_double(val, d)) { m.gps_cfg.datum_lon_deg = d; }
            else if (key == "datum_alt_m"   && to_double(val, d)) { m.gps_cfg.datum_alt_m = d; }
            else if (key == "odom_from_enu_yaw"   && to_double(val, d)) { /* set below */ m.gps_cfg.odom_from_enu.R = so3::exp(Vec3(0,0,d)) * m.gps_cfg.odom_from_enu.R; }
            else if (key == "lever_x" && to_double(val, d)) { m.gps_cfg.lever_arm_base.x() = d; }
            else if (key == "lever_y" && to_double(val, d)) { m.gps_cfg.lever_arm_base.y() = d; }
            else if (key == "lever_z" && to_double(val, d)) { m.gps_cfg.lever_arm_base.z() = d; }
            else if (key == "cov_floor_m2" && to_double(val, d)) { m.gps_cfg.cov_floor_m2 = d; }
            else return err("unknown [gps] key '" + key + "'");
            continue;
        }
        if (sec == Sec::Gt) {
            if (key == "csv") { m.gt_csv = val; continue; }
            return err("unknown [gt] key '" + key + "'");
        }
        if (sec == Sec::Replay) {
            double d; long long iv;
            if (key == "tail_window_s" && to_double(val, d)) { m.tail_window_s = static_cast<Scalar>(d); }
            else if (key == "warmup_steps" && to_i64(val, iv)) { m.warmup_steps = static_cast<int>(iv); }
            else if (key == "local_batch_len" && to_i64(val, iv)) { m.local_batch_len = static_cast<int>(iv); }
            else if (key == "out") { m.out_path = val; }
            else return err("unknown [replay] key '" + key + "'");
            continue;
        }
    }
    m.config_text = cfg_out.str();
    return "";
}

// Load the GPS fix CSV (t_ns, lat, lon, alt [, var_e, var_n, var_u]).
std::string load_gps_fixes(const std::string& path, std::vector<GpsFix>& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return "could not open GPS csv: " + path;
    std::string raw;
    int line_no = 0;
    bool seen_data = false;
    while (std::getline(f, raw)) {
        ++line_no;
        const std::string body = trim(strip_comment(raw));
        if (body.empty()) continue;
        std::vector<std::string> fld = split_fields(body);
        if (fld.empty()) continue;
        long long t;
        if (!to_i64(fld[0], t)) { if (!seen_data) continue; /* header */ return "gps line " + std::to_string(line_no) + ": bad timestamp"; }
        if (fld.size() != 4 && fld.size() != 7)
            return "gps line " + std::to_string(line_no) + ": expected 4 or 7 columns";
        double lat, lon, alt;
        if (!to_double(fld[1], lat) || !to_double(fld[2], lon) || !to_double(fld[3], alt))
            return "gps line " + std::to_string(line_no) + ": bad lat/lon/alt";
        GpsFix fx;
        fx.lat_deg = lat; fx.lon_deg = lon; fx.alt_m = alt; fx.stamp = static_cast<Timestamp>(t);
        fx.valid = true;
        if (fld.size() == 7) {
            double ve, vn, vu;
            if (!to_double(fld[4], ve) || !to_double(fld[5], vn) || !to_double(fld[6], vu))
                return "gps line " + std::to_string(line_no) + ": bad covariance";
            fx.cov_enu = Vec3(ve, vn, vu).asDiagonal();
        } else {
            fx.cov_enu = Mat3::Identity();   // 1 m^2 isotropic default
        }
        out.push_back(fx);
        seen_data = true;
    }
    if (out.empty()) return "gps csv had no fixes: " + path;
    return "";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: ofc_replay <manifest.ini> [results_out.csv]\n"
                     "  reads the manifest, replays the CsvSources (+ optional GPS/GT) through the\n"
                     "  Estimator, writes the results CSV + a summary to stdout.\n";
        return 2;
    }
    const std::string manifest_path = argv[1];

    std::string manifest_text;
    {
        std::ifstream f(manifest_path, std::ios::binary);
        if (!f) { std::cerr << "error: cannot open manifest: " << manifest_path << "\n"; return 1; }
        std::ostringstream ss; ss << f.rdbuf(); manifest_text = ss.str();
    }

    Manifest man;
    const std::string mperr = parse_manifest(manifest_text, man);
    if (!mperr.empty()) { std::cerr << "error: " << mperr << "\n"; return 1; }
    if (argc >= 3) man.out_path = argv[2];

    // Build the core Config via the existing loader (it owns the sensor storage).
    ConfigLoader loader;
    const Status cs = loader.parse(man.config_text);
    if (!ok(cs)) { std::cerr << "error: config: " << loader.error() << "\n"; return 1; }
    const Config& cfg = loader.config();

    // Load each source CSV. Match the [sensor.N] slot to its SensorConfig (same dense order the
    // loader uses) so the modeled-cov knobs flow through.
    std::vector<CsvSource> sources(man.source_entries.size());
    std::vector<CsvSource*> source_ptrs;
    int si = 0;
    for (const auto& kv : man.source_entries) {
        const int slot = kv.first;
        const SourceEntry& e = kv.second;
        if (e.csv.empty()) { std::cerr << "error: [sensor." << slot << "] has no csv= path\n"; return 1; }
        CsvSourceConfig sc;
        sc.id = e.id;
        sc.combine = cfg.conf_combine;
        sc.confidence_blend = cfg.confidence_blend;
        if (slot < cfg.sensor_count && cfg.sensors != nullptr) {
            const SensorConfig& sen = cfg.sensors[slot];
            sc.id = sen.id;
            sc.native_confidence     = sen.native_confidence;
            sc.modeled_noise_per_m   = sen.modeled_noise_per_m;
            sc.modeled_noise_per_rad = sen.modeled_noise_per_rad;
        }
        if (e.has_form) { sc.form = e.form; sc.force_form = true; }
        const Status ls = sources[si].load_file(e.csv, sc);
        if (!ok(ls)) { std::cerr << "error: source [sensor." << slot << "] csv '" << e.csv
                                 << "': " << sources[si].error() << "\n"; return 1; }
        source_ptrs.push_back(&sources[si]);
        ++si;
    }
    if (source_ptrs.empty()) { std::cerr << "error: no sources in the manifest\n"; return 1; }

    // Optional GPS.
    GpsCorrection gps(man.gps_cfg);
    std::vector<GpsFix> gps_fixes;
    if (man.has_gps) {
        if (man.gps_csv.empty()) { std::cerr << "error: [gps] has no csv= path\n"; return 1; }
        const std::string ge = load_gps_fixes(man.gps_csv, gps_fixes);
        if (!ge.empty()) { std::cerr << "error: " << ge << "\n"; return 1; }
    }

    // Optional GT track.
    CsvGtTrack gt;
    if (man.has_gt) {
        if (man.gt_csv.empty()) { std::cerr << "error: [gt] has no csv= path\n"; return 1; }
        const Status gs = gt.load_file(man.gt_csv);
        if (!ok(gs)) { std::cerr << "error: gt csv '" << man.gt_csv << "': " << gt.error() << "\n"; return 1; }
    }

    // Run the harness.
    ReplayInputs in;
    in.cfg = &cfg;
    in.sources = source_ptrs;
    if (man.has_gps) { in.gps = &gps; in.gps_fixes = gps_fixes; }
    if (man.has_gt)  { in.gt = &gt; }
    in.tail_window_s   = man.tail_window_s;
    in.warmup_steps    = man.warmup_steps;
    in.local_batch_len = man.local_batch_len;

    ReplayHarness h;
    const Status rs = h.run(in);
    if (!ok(rs)) { std::cerr << "error: replay: " << h.error() << "\n"; return 1; }

    const Status ws = h.write_results_csv_file(man.out_path);
    if (!ok(ws)) { std::cerr << "error: could not write results to " << man.out_path << "\n"; return 1; }

    std::cout << "wrote " << man.out_path << "\n";
    h.write_summary(std::cout);
    return 0;
}
