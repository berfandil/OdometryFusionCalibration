// ofc_adapters/csv_source.hpp — REAL-DATA INGESTION (Slice 13, real-data scaffolding, D7/D19).
//
// RELAXED EDGE (adapters/). Never ships to the ECU as the strict core does: std containers /
// std::fstream / heap / exceptions are all fine here. Depends only on the PUBLIC ofc/core
// headers (via ofc_core) + Eigen. NO csv library — the CSV is parsed by hand (the project keeps
// deps to Eigen + doctest), matching the dependency-free hand-rolled style of config_loader.
//
// WHY: today the system is exercised only by the SIM rig (synthetic sources from ground truth).
// To run on a REAL dataset (KITTI/EuRoC/...) we need (1) an ISource that reads real odometry from
// a file and answers the core's delta(t0,t1) contract, (2) a GT-track reader, and (3) a
// replay+eval harness (see replay_harness.hpp) that pumps the Estimator over a timeline and
// computes drift/NEES/NIS vs GT. The user converts any public dataset to the CSV schema below
// with a small script.
//
// =====================================================================================
// CSV SCHEMA (all numbers are decimal; `#`, `//`, and blank lines are skipped; a leading header
// row of non-numeric column names is tolerated; whitespace is trimmed; the delimiter is a comma
// OR any run of whitespace). The native form is selected by a header TAG line `# form: <form>`
// (anywhere before data, case-insensitive: absolute|increment|twist) OR by the ctor `form` arg;
// the ctor arg wins if both are present and disagree.
//
//   absolute / increment:
//     t_ns, x, y, z, qw, qx, qy, qz [, var_tx, var_ty, var_tz, var_rx, var_ry, var_rz]
//       * t_ns   — sample time, nanoseconds (int64), STRICTLY increasing.
//       * x,y,z  — translation (m). For `absolute` this is the cumulative pose in the source's
//                  own integrated frame; for `increment` it is the per-step relative motion from
//                  the previous row (the FIRST increment row is ignored as the identity seed).
//       * qw,qx,qy,qz — unit quaternion (w-first), normalized on read.
//       * the optional 6 trailing numbers are the DIAGONAL of the per-step increment covariance,
//         in [trans; rot] order (matching Delta::cov). Absent -> the buffer's modeled cov path
//         (identity native, combined per Config::conf_combine) is used.
//
//   twist:
//     t_ns, vx, vy, vz, wx, wy, wz [, var_tx, var_ty, var_tz, var_rx, var_ry, var_rz]
//       * vx,vy,vz — body linear velocity (m/s); wx,wy,wz — body angular velocity (rad/s).
//       * the optional 6 vars are the per-SECOND twist covariance diagonal in [trans; rot] order.
//
//   GT track (CsvGtTrack — eval only, absolute pose):
//     t_ns, x, y, z, qw, qx, qy, qz       (no covariance columns; extra columns are ignored)
// =====================================================================================
//
// QUATERNION CONVENTION: (qw, qx, qy, qz), w-first, normalized on read; the rotation matrix is
// Eigen::Quaternion(w,x,y,z).toRotationMatrix(). A degenerate (near-zero-norm) quaternion is a
// parse error. The inverse (Mat3 -> quaternion) is provided as a small adapter helper for the
// harness's pose logging — quat<->Mat3 lives HERE in the adapter, NOT in the strict core.
#ifndef OFC_ADAPTERS_CSV_SOURCE_HPP
#define OFC_ADAPTERS_CSV_SOURCE_HPP

#include "ofc/core/buffer.hpp"
#include "ofc/core/config.hpp"
#include "ofc/core/source.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include <string>
#include <vector>

namespace ofc {
namespace adapters {

// --- quaternion <-> rotation helpers (adapter-local; NOT in the strict core) -------------------
// (qw, qx, qy, qz) w-first. quat_to_mat3 normalizes; a near-zero-norm input returns identity and
// sets ok=false (the parser treats that as a malformed row). mat3_to_quat returns w-first too.
Mat3 quat_to_mat3(Scalar qw, Scalar qx, Scalar qy, Scalar qz, bool* ok = nullptr);
void mat3_to_quat(const Mat3& R, Scalar& qw, Scalar& qx, Scalar& qy, Scalar& qz);

// Parameters for a CsvSource (everything the SourceBuffer needs that does not come from the file).
struct CsvSourceConfig {
    SourceId    id            = 0;          // the source id (the core keys calibration on this)
    OdomForm    form          = OdomForm::Increment;  // overridden by `# form:` unless force_form
    bool        force_form    = false;      // true => `form` wins over any `# form:` tag
    int         capacity      = 0;          // SourceBuffer capacity; 0 => sized to the row count + 2
    ConfCombine combine       = ConfCombine::Sum;  // native+modeled combine (CONFIG §4)
    Scalar      confidence_blend = Scalar(0.5);    // blend factor when combine == Weighted

    // The per-sensor knobs the buffer's modeled-cov path uses (CONFIG §9). Defaults match
    // SensorConfig; the manifest/CLI fills these from the parsed Config so a CsvSource's modeled
    // cov matches the estimator's view of the sensor.
    bool   native_confidence      = true;
    Scalar modeled_noise_per_m    = Scalar(0.01);
    Scalar modeled_noise_per_rad  = Scalar(0.01);
};

// An ISource backed by a CSV file. Loads the file into a SourceBuffer (REUSING the same buffer
// the sim's SyntheticSource::build_buffer uses) and answers query(t0,t1) straight from it — so
// the CSV ingestion path is EQUIVALENT to the in-memory path (proved by the equivalence test).
class CsvSource : public ISource {
public:
    CsvSource() = default;

    // Parse `text` (the whole CSV) into the internal buffer. On success the source answers
    // query(); on failure error() describes the first problem and the buffer is left empty.
    //   Ok            — parsed + at least 2 rows ingested.
    //   InvalidConfig — a malformed row / bad quaternion / wrong column count / a `# form:` tag
    //                   disagreeing with a forced form. error() holds the message + line number.
    //   NoData        — fewer than 2 usable data rows (the buffer cannot answer a delta).
    Status load(const std::string& text, const CsvSourceConfig& cfg);

    // Convenience: read a file then load() it. A missing/unreadable file -> NoData / CorruptData.
    Status load_file(const std::string& path, const CsvSourceConfig& cfg);

    // ISource.
    SourceId id() const override { return cfg_.id; }
    Expected<Delta> query(Timestamp t0, Timestamp t1) const override;

    // Stamp range of the loaded data (0/0 if empty) — the harness merges these into the timeline.
    Timestamp oldest() const { return buf_.oldest(); }
    Timestamp newest() const { return buf_.newest(); }
    int       row_count() const { return buf_.size(); }
    bool      loaded() const { return loaded_; }

    const std::string& error() const { return error_; }
    OdomForm           form()  const { return cfg_.form; }

private:
    SourceBuffer    buf_;
    CsvSourceConfig cfg_;
    bool            loaded_ = false;
    std::string     error_;
};

// A ground-truth absolute-pose track (EVAL ONLY — not an ISource). Loads `t_ns,x,y,z,qw..qz`
// rows and answers pose_at(t) by LINEAR SE(3) interpolation (se3::interpolate) between the two
// bracketing samples; clamps to the nearest end outside the data range. The harness compares the
// fused frontier against this.
class CsvGtTrack {
public:
    CsvGtTrack() = default;

    // Parse `text` into the internal sorted sample list.
    //   Ok            — >= 1 sample.
    //   InvalidConfig — a malformed row / bad quaternion / out-of-order stamp. error() set.
    //   NoData        — no samples.
    Status load(const std::string& text);
    Status load_file(const std::string& path);

    // Interpolated GT pose at `t` (clamped at the ends). Precondition: loaded() (else identity).
    SE3 pose_at(Timestamp t) const;

    Timestamp oldest() const { return samples_.empty() ? Timestamp(0) : samples_.front().stamp; }
    Timestamp newest() const { return samples_.empty() ? Timestamp(0) : samples_.back().stamp; }
    int       size()   const { return static_cast<int>(samples_.size()); }
    bool      loaded() const { return !samples_.empty(); }
    const std::string& error() const { return error_; }

private:
    struct Sample { Timestamp stamp = 0; SE3 pose; };
    std::vector<Sample> samples_;
    std::string         error_;
};

} // namespace adapters
} // namespace ofc
#endif // OFC_ADAPTERS_CSV_SOURCE_HPP
