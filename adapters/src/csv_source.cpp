// ofc_adapters/csv_source.cpp — see csv_source.hpp for the schema + conventions.
//
// RELAXED EDGE: std::string / std::vector / std::fstream / heap / exceptions fine (no exception
// escapes a public method — std::stod/stoll throws are caught and mapped to InvalidConfig). The
// CSV is parsed BY HAND (no csv library), mirroring config_loader's hand-rolled lexer style.
#include "ofc_adapters/csv_source.hpp"

#include "ofc/core/lie.hpp"

#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ofc {
namespace adapters {

namespace {

// --- small string helpers (same conventions as config_loader) ----------------------------------
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

// Strip a `#`/`//` comment to end of line. (A bare `#` at column 0 is a full-line comment; a
// trailing comment after data is also honored. `;` is NOT a comment here — CSVs may legitimately
// carry it, and the format tag uses `#`.)
std::string strip_comment(const std::string& line) {
    std::size_t cut = line.size();
    const std::size_t h = line.find('#');
    if (h != std::string::npos && h < cut) cut = h;
    const std::size_t s = line.find("//");
    if (s != std::string::npos && s < cut) cut = s;
    return line.substr(0, cut);
}

// Split on commas OR runs of whitespace (so both `1,2,3` and `1 2 3` and `1, 2 ,3` work). Empty
// fields (e.g. a trailing comma) are dropped.
std::vector<std::string> split_fields(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    auto flush = [&]() {
        const std::string t = trim(cur);
        if (!t.empty()) out.push_back(t);
        cur.clear();
    };
    for (char c : line) {
        if (c == ',' || std::isspace(static_cast<unsigned char>(c))) flush();
        else cur.push_back(c);
    }
    flush();
    return out;
}

bool parse_double(const std::string& v, Scalar& out) {
    try {
        std::size_t used = 0;
        const double d = std::stod(v, &used);
        if (used == 0) return false;
        if (!trim(v.substr(used)).empty()) return false;   // trailing garbage
        out = static_cast<Scalar>(d);
        return std::isfinite(d);
    } catch (...) {
        return false;
    }
}

bool parse_i64(const std::string& v, std::int64_t& out) {
    try {
        std::size_t used = 0;
        const long long n = std::stoll(v, &used, 10);
        if (used == 0) return false;
        if (!trim(v.substr(used)).empty()) return false;
        out = static_cast<std::int64_t>(n);
        return true;
    } catch (...) {
        return false;
    }
}

// A line with a leading `# form: <tag>` sets the native form. Returns true + sets `form` if the
// line is a form tag (anywhere in the comment-stripped raw line's `#` body); false otherwise.
bool parse_form_tag(const std::string& raw, OdomForm& form) {
    const std::size_t h = raw.find('#');
    if (h == std::string::npos) return false;
    std::string body = lower(trim(raw.substr(h + 1)));
    const std::string key = "form:";
    const std::size_t k = body.find(key);
    if (k == std::string::npos) return false;
    const std::string tag = trim(body.substr(k + key.size()));
    if (tag.compare(0, 8, "absolute") == 0)      { form = OdomForm::AbsolutePose; return true; }
    if (tag.compare(0, 9, "increment") == 0)     { form = OdomForm::Increment;    return true; }
    if (tag.compare(0, 5, "twist") == 0)         { form = OdomForm::Twist;        return true; }
    return false;
}

// Diagonal Mat6 from 6 [trans;rot] variances.
Mat6 diag6(const Scalar v[6]) {
    Mat6 m = Mat6::Zero();
    for (int i = 0; i < 6; ++i) m(i, i) = v[i];
    return m;
}

} // namespace

// --- quaternion helpers (adapter-local; NOT in the strict core) --------------------------------
Mat3 quat_to_mat3(Scalar qw, Scalar qx, Scalar qy, Scalar qz, bool* ok) {
    const Scalar n2 = qw * qw + qx * qx + qy * qy + qz * qz;
    if (!(n2 > Scalar(1e-12)) || !std::isfinite(n2)) {
        if (ok) *ok = false;
        return Mat3::Identity();
    }
    if (ok) *ok = true;
    // Eigen::Quaternion ctor is (w, x, y, z); normalized() handles the unit constraint.
    Eigen::Quaternion<Scalar> q(qw, qx, qy, qz);
    q.normalize();
    return q.toRotationMatrix();
}

void mat3_to_quat(const Mat3& R, Scalar& qw, Scalar& qx, Scalar& qy, Scalar& qz) {
    Eigen::Quaternion<Scalar> q(R);
    q.normalize();
    if (q.w() < Scalar(0)) {           // canonical hemisphere (w >= 0) so logs are stable
        qw = -q.w(); qx = -q.x(); qy = -q.y(); qz = -q.z();
    } else {
        qw = q.w(); qx = q.x(); qy = q.y(); qz = q.z();
    }
}

// ================================================================================================
// CsvSource
// ================================================================================================
Status CsvSource::load(const std::string& text, const CsvSourceConfig& cfg) {
    cfg_    = cfg;
    loaded_ = false;
    error_.clear();
    buf_ = SourceBuffer{};   // reset any previous load

    auto fail = [&](const std::string& msg, int line_no, const std::string& line) -> Status {
        std::ostringstream os;
        os << "line " << line_no << ": " << msg << "  [" << trim(line) << "]";
        error_ = os.str();
        return Status::InvalidConfig;
    };

    // --- pass 1: scan for the form tag + collect the data rows --------------------------------
    OdomForm form = cfg.form;
    bool form_from_tag = false;
    std::vector<std::vector<std::string>> rows;
    std::vector<int> row_lines;

    std::istringstream in(text);
    std::string raw;
    int line_no = 0;
    while (std::getline(in, raw)) {
        ++line_no;
        // A `# form:` tag is honored before stripping.
        OdomForm tag_form;
        if (parse_form_tag(raw, tag_form)) {
            if (cfg.force_form) {
                if (tag_form != cfg.form)
                    return fail("`# form:` tag disagrees with the forced ctor form", line_no, raw);
            } else {
                form = tag_form;
                form_from_tag = true;
            }
            continue;
        }
        const std::string body = trim(strip_comment(raw));
        if (body.empty()) continue;
        std::vector<std::string> f = split_fields(body);
        if (f.empty()) continue;
        // Tolerate a header row: if the FIRST field is non-numeric, skip the row.
        std::int64_t probe;
        if (!parse_i64(f[0], probe)) {
            // Only skip a *leading* header (no rows yet); a non-numeric stamp later is an error.
            if (rows.empty()) continue;
            return fail("non-numeric timestamp", line_no, raw);
        }
        rows.push_back(std::move(f));
        row_lines.push_back(line_no);
    }
    (void)form_from_tag;

    if (rows.size() < 2) {
        error_ = "CSV has fewer than 2 usable data rows";
        return Status::NoData;
    }

    cfg_.form = form;   // record the resolved form

    // --- configure the buffer ------------------------------------------------------------------
    SensorConfig sc;
    sc.id                   = cfg.id;
    sc.native_confidence    = cfg.native_confidence;
    sc.modeled_noise_per_m  = cfg.modeled_noise_per_m;
    sc.modeled_noise_per_rad= cfg.modeled_noise_per_rad;
    const int cap = (cfg.capacity > 0) ? cfg.capacity
                                       : (static_cast<int>(rows.size()) + 2);
    const Status cs = buf_.configure(sc, form, cap, cfg.combine, cfg.confidence_blend);
    if (!ok(cs)) {
        error_ = "SourceBuffer.configure rejected the parameters";
        return cs;
    }

    // Column expectations: absolute/increment need 8 (t + xyz + quat) + optional 6; twist needs
    // 7 (t + 6) + optional 6.
    const bool is_twist = (form == OdomForm::Twist);
    const int  base_cols = is_twist ? 7 : 8;

    // --- pass 2: ingest -----------------------------------------------------------------------
    for (std::size_t r = 0; r < rows.size(); ++r) {
        const std::vector<std::string>& f = rows[r];
        const int ln = row_lines[r];
        const std::string joined = [&]() { std::string s; for (auto& x : f) { s += x; s += ' '; } return s; }();
        if (static_cast<int>(f.size()) != base_cols && static_cast<int>(f.size()) != base_cols + 6)
            return fail("wrong column count (expected " + std::to_string(base_cols) + " or " +
                            std::to_string(base_cols + 6) + ", got " + std::to_string(f.size()) + ")",
                        ln, joined);

        std::int64_t t_ns;
        if (!parse_i64(f[0], t_ns)) return fail("bad timestamp", ln, joined);

        // Optional covariance diagonal.
        const bool have_cov = (static_cast<int>(f.size()) == base_cols + 6);
        Scalar cv[6] = {0, 0, 0, 0, 0, 0};
        if (have_cov) {
            for (int i = 0; i < 6; ++i)
                if (!parse_double(f[base_cols + i], cv[i]))
                    return fail("bad covariance value", ln, joined);
        }
        const Mat6 cov = have_cov ? diag6(cv) : Mat6::Zero();

        if (is_twist) {
            Scalar v[6];
            for (int i = 0; i < 6; ++i)
                if (!parse_double(f[1 + i], v[i])) return fail("bad twist value", ln, joined);
            Vec6 xi; xi << v[0], v[1], v[2], v[3], v[4], v[5];
            const Status ps = buf_.push_twist(static_cast<Timestamp>(t_ns), xi, cov);
            if (!ok(ps)) return fail(ps == Status::OutOfRange
                                         ? "non-increasing timestamp" : "buffer push failed",
                                     ln, joined);
        } else {
            Scalar p[3], q[4];
            for (int i = 0; i < 3; ++i)
                if (!parse_double(f[1 + i], p[i])) return fail("bad position value", ln, joined);
            for (int i = 0; i < 4; ++i)
                if (!parse_double(f[4 + i], q[i])) return fail("bad quaternion value", ln, joined);
            bool qok = false;
            const Mat3 R = quat_to_mat3(q[0], q[1], q[2], q[3], &qok);
            if (!qok) return fail("degenerate (near-zero) quaternion", ln, joined);
            SE3 pose; pose.R = R; pose.t = Vec3(p[0], p[1], p[2]);
            Status ps;
            if (form == OdomForm::AbsolutePose)
                ps = buf_.push_absolute(static_cast<Timestamp>(t_ns), pose, cov);
            else  // Increment
                ps = buf_.push_increment(static_cast<Timestamp>(t_ns), pose, cov);
            if (!ok(ps)) return fail(ps == Status::OutOfRange
                                         ? "non-increasing timestamp" : "buffer push failed",
                                     ln, joined);
        }
    }

    loaded_ = true;
    return Status::Ok;
}

Status CsvSource::load_file(const std::string& path, const CsvSourceConfig& cfg) {
    std::string text;
    try {
        std::ifstream fcheck(path, std::ios::binary);
        if (!fcheck) { error_ = "could not open file: " + path; return Status::NoData; }
        std::ostringstream ss;
        ss << fcheck.rdbuf();
        text = ss.str();
    } catch (...) {
        error_ = "IO error reading file: " + path;
        return Status::CorruptData;
    }
    return load(text, cfg);
}

Expected<Delta> CsvSource::query(Timestamp t0, Timestamp t1) const {
    if (!loaded_) return Expected<Delta>(Status::NoData);
    return buf_.query(t0, t1);
}

// ================================================================================================
// CsvGtTrack
// ================================================================================================
Status CsvGtTrack::load(const std::string& text) {
    samples_.clear();
    error_.clear();

    auto fail = [&](const std::string& msg, int line_no, const std::string& line) -> Status {
        std::ostringstream os;
        os << "line " << line_no << ": " << msg << "  [" << trim(line) << "]";
        error_ = os.str();
        return Status::InvalidConfig;
    };

    std::istringstream in(text);
    std::string raw;
    int line_no = 0;
    while (std::getline(in, raw)) {
        ++line_no;
        const std::string body = trim(strip_comment(raw));
        if (body.empty()) continue;
        std::vector<std::string> f = split_fields(body);
        if (f.empty()) continue;
        std::int64_t t_ns;
        if (!parse_i64(f[0], t_ns)) {
            if (samples_.empty()) continue;   // tolerate a leading header row
            return fail("non-numeric timestamp", line_no, raw);
        }
        if (f.size() < 8) return fail("GT row needs at least 8 columns (t,x,y,z,qw,qx,qy,qz)",
                                      line_no, raw);
        Scalar p[3], q[4];
        for (int i = 0; i < 3; ++i)
            if (!parse_double(f[1 + i], p[i])) return fail("bad position value", line_no, raw);
        for (int i = 0; i < 4; ++i)
            if (!parse_double(f[4 + i], q[i])) return fail("bad quaternion value", line_no, raw);
        bool qok = false;
        const Mat3 R = quat_to_mat3(q[0], q[1], q[2], q[3], &qok);
        if (!qok) return fail("degenerate (near-zero) quaternion", line_no, raw);
        if (!samples_.empty() && static_cast<Timestamp>(t_ns) <= samples_.back().stamp)
            return fail("non-increasing GT timestamp", line_no, raw);
        Sample s;
        s.stamp  = static_cast<Timestamp>(t_ns);
        s.pose.R = R;
        s.pose.t = Vec3(p[0], p[1], p[2]);
        samples_.push_back(s);
    }

    if (samples_.empty()) {
        error_ = "GT track has no usable samples";
        return Status::NoData;
    }
    return Status::Ok;
}

Status CsvGtTrack::load_file(const std::string& path) {
    std::string text;
    try {
        std::ifstream fcheck(path, std::ios::binary);
        if (!fcheck) { error_ = "could not open file: " + path; return Status::NoData; }
        std::ostringstream ss;
        ss << fcheck.rdbuf();
        text = ss.str();
    } catch (...) {
        error_ = "IO error reading file: " + path;
        return Status::CorruptData;
    }
    return load(text);
}

SE3 CsvGtTrack::pose_at(Timestamp t) const {
    if (samples_.empty()) return SE3{};
    if (t <= samples_.front().stamp) return samples_.front().pose;
    if (t >= samples_.back().stamp)  return samples_.back().pose;
    // Binary search for the first sample with stamp >= t (bracketing pair is [hi-1, hi]).
    std::size_t lo = 0, hi = samples_.size();
    while (lo < hi) {
        const std::size_t mid = lo + (hi - lo) / 2;
        if (samples_[mid].stamp < t) lo = mid + 1;
        else                          hi = mid;
    }
    // lo is the first index with stamp >= t; lo >= 1 here (t > front, t < back).
    const Sample& a = samples_[lo - 1];
    const Sample& b = samples_[lo];
    if (b.stamp == t) return b.pose;
    const Scalar span = static_cast<Scalar>(b.stamp - a.stamp);
    const Scalar u = (span > Scalar(0)) ? static_cast<Scalar>(t - a.stamp) / span : Scalar(0);
    return se3::interpolate(a.pose, b.pose, u);
}

} // namespace adapters
} // namespace ofc
