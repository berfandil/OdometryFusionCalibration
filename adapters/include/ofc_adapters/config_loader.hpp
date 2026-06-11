// ofc_adapters/config_loader.hpp — dependency-free text -> Config (Slice 13, D19).
//
// RELAXED EDGE. DESIGN §9 / D19: "File parsing (YAML/JSON/ROS-param) is an adapter that
// BUILDS the struct." The strict core only ever receives a validated `Config`. This adapter
// parses a SIMPLE key=value / INI text format — NO external parser (the project keeps deps to
// Eigen + doctest), so this is a tiny hand-rolled lexer, not a YAML/JSON library.
//
// FORMAT (a documented SUBSET of CONFIG.md knobs):
//   * `#` or `;` begin a comment (to end of line). Blank lines ignored. Whitespace trimmed.
//   * `key = value` sets a knob. Booleans: true/false/1/0/on/off/yes/no (case-insensitive).
//   * `[sensor.<N>]` opens a per-sensor section; subsequent keys set that SensorConfig. <N> is
//     the 0-based slot index (sensors are stored densely in section order; the index is for
//     readability + duplicate detection, not a sparse array).
//   * `[global]` (or no section) sets top-level Config knobs.
//
// SUPPORTED TOP-LEVEL KEYS (see kKnobDoc in the .cpp for the authoritative list):
//   max_sources, reference_sensor_id, window_s, fusion_delay_s, tick_rate_hz, calib_lag_s,
//   cold_start (reference_only|median_from_start), commit_concentration, commit_min_votes,
//   commit_drop, straight_omega_max, straight_trans_min, turn_omega_min, timesync_enabled,
//   rot3d_enabled (bool — turn-regime FULL rotation extrinsic, Slice 17; default off),
//   joint_lever_scale (bool — turn-regime joint lever+scale 4-unknown LS, Slice 17b; default off),
//   subbin_centroid (bool -> all five calibration histograms' HistogramConfig::subbin_centroid),
//   adaptive_q (bool), q_scale (f), q_floor (1 number -> all 6 axes, or 6 numbers [trans;rot]).
// SUPPORTED PER-SENSOR KEYS ([sensor.N]):
//   id, prior_extrinsic (yaw pitch roll x y z — 6 numbers, ZYX euler + translation),
//   prior_scale, prior_time_offset_s, weight_prior, per_sensor_kf, scale_calib,
//   is_reference.
//
// UNKNOWN KEYS / SECTIONS -> a hard error (Status::InvalidConfig) with the offending line in
// `error()` — typo-loud rather than silently dropping a misspelled knob. (Documented choice;
// switch to skip-with-warning by relaxing parse() if a forward-compat use-case appears.)
//
// OWNERSHIP (mirrors Config::sensors being a CALLER-OWNED pointer): the loader OWNS the
// std::vector<SensorConfig> backing Config::sensors. The returned Config's `sensors` pointer
// is valid ONLY for the lifetime of the ConfigLoader that produced it (and is invalidated by a
// re-parse). The caller holds the ConfigLoader for as long as it uses the Config — exactly the
// "per-sensor SensorConfig storage must outlive the Config" contract.
#ifndef OFC_ADAPTERS_CONFIG_LOADER_HPP
#define OFC_ADAPTERS_CONFIG_LOADER_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/status.hpp"

#include <string>
#include <vector>

namespace ofc {
namespace adapters {

class ConfigLoader {
public:
    ConfigLoader() = default;

    // Parse `text` into the internal Config + the owned sensor vector, then run the core
    // validate(). On success the built Config (with `sensors` pointing at the owned vector) is
    // available via config(); on failure config() is the partially-built state and error()
    // describes the first problem.
    //   Ok            — parsed + validate() == Ok.
    //   InvalidConfig — a syntax error, an unknown key/section, a malformed value, or a
    //                   duplicate sensor slot. error() holds the message + offending line.
    //   (any validate() status) — the text parsed but the core rejected the config; error()
    //                   notes it failed validate().
    Status parse(const std::string& text);

    // Convenience: read a file then parse() it. NotInitialized-class IO failures map to
    // CorruptData (the file could not be read); a missing file -> NoData.
    Status parse_file(const std::string& path);

    // The built Config. Its `sensors` pointer is owned by THIS loader (see ownership note in
    // the header doc) — valid until the loader is destroyed or re-parsed.
    const Config& config() const { return cfg_; }

    // The owned per-sensor storage (also reachable via config().sensors). Exposed for callers
    // that want to tweak a sensor before re-pointing Config::sensors.
    const std::vector<SensorConfig>& sensors() const { return sensors_; }

    // Human-readable description of the first parse/validate failure ("" on success).
    const std::string& error() const { return error_; }

    // The supported-knob reference text (the authoritative key list), for `--help`-style
    // surfacing. Pure documentation; no parsing side effects.
    static const char* knob_doc();

private:
    Config                    cfg_;
    std::vector<SensorConfig> sensors_;
    std::string               error_;
};

// Free-function convenience wrapper. Parses `text` into `out_cfg`, taking ownership of the
// sensor storage in `out_sensors` (the caller MUST keep `out_sensors` alive as long as it uses
// `out_cfg`, and must repoint out_cfg.sensors at out_sensors.data() — done here for them).
Status load_config(const std::string& text, Config& out_cfg,
                   std::vector<SensorConfig>& out_sensors, std::string* out_error = nullptr);

} // namespace adapters
} // namespace ofc
#endif // OFC_ADAPTERS_CONFIG_LOADER_HPP
