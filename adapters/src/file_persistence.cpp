// ofc_adapters/file_persistence.cpp — see file_persistence.hpp for the contract.
//
// RELAXED EDGE: std::fstream / std::vector / heap / exceptions are fine. The one hard rule
// is that NO exception escapes a public method — IO failures are caught and mapped to Status.
#include "ofc_adapters/file_persistence.hpp"

#include "ofc/core/persistence.hpp"   // public core: blob layout + checksum, for the validity check

#include <cstring>
#include <fstream>
#include <ios>
#include <vector>

namespace ofc {
namespace adapters {

namespace {

constexpr int kSeqBytes = 8;   // uint64 seq
constexpr int kLenBytes = 4;   // uint32 blob_len
constexpr int kHdrBytes = kSeqBytes + kLenBytes;

void put_u64_le(unsigned char* b, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) b[i] = static_cast<unsigned char>((v >> (8 * i)) & 0xFFu);
}
void put_u32_le(unsigned char* b, std::uint32_t v) {
    for (int i = 0; i < 4; ++i) b[i] = static_cast<unsigned char>((v >> (8 * i)) & 0xFFu);
}
std::uint64_t get_u64_le(const unsigned char* b) {
    std::uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(b[i]) << (8 * i);
    return v;
}
std::uint32_t get_u32_le(const unsigned char* b) {
    std::uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v |= static_cast<std::uint32_t>(b[i]) << (8 * i);
    return v;
}

// One parsed on-disk record.
struct Record {
    bool                       present = false;   // file existed + held a full header
    std::uint64_t              seq     = 0;
    std::vector<unsigned char> blob;              // exactly blob_len bytes (may be torn-short)
};

// Read one record file. A missing file, a file too short to hold the header, or any IO
// exception -> present=false (treated as "no record"). A header that claims more blob bytes
// than are actually on disk (a torn write) yields present=true with a SHORT blob; the core
// deserialize then rejects it on length/checksum, which is exactly the fallback we want.
Record read_record(const std::string& path) {
    Record rec;
    try {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) return rec;
        const std::streamoff sz = f.tellg();
        if (sz < static_cast<std::streamoff>(kHdrBytes)) return rec;  // no full header

        f.seekg(0);
        unsigned char hdr[kHdrBytes];
        f.read(reinterpret_cast<char*>(hdr), kHdrBytes);
        if (!f) return rec;

        rec.seq = get_u64_le(hdr);
        const std::uint32_t claimed = get_u32_le(hdr + kSeqBytes);
        const std::streamoff on_disk = sz - kHdrBytes;
        // Read whatever is actually present, up to the claimed length. If the file is short
        // (torn write) we hand the truncated blob to the core, which rejects it.
        std::size_t to_read = static_cast<std::size_t>(claimed);
        if (static_cast<std::streamoff>(to_read) > on_disk)
            to_read = static_cast<std::size_t>(on_disk);
        rec.blob.resize(to_read);
        if (to_read > 0) {
            f.read(reinterpret_cast<char*>(rec.blob.data()),
                   static_cast<std::streamsize>(to_read));
            if (!f) { rec.blob.resize(static_cast<std::size_t>(f.gcount())); }
        }
        rec.present = true;
    } catch (...) {
        return Record{};   // any IO exception -> no record
    }
    return rec;
}

// Lightweight, NON-DESTRUCTIVE validity check on a candidate OFCP blob — the exact subset of the
// core deserialize() framing guards that decide whether load() would ACCEPT this record:
// magic -> version -> exact length -> checksum -> config-hash. It deliberately does NOT touch an
// Estimator (so save() can test a candidate without a scratch estimator and without mutating the
// live one): instead it compares the candidate's config-hash to `want_cfg_hash`, the config-hash
// of the estimator we are about to save (extracted from the blob we just serialized). A blob that
// passes here is one the core deserialize() would accept for the SAME rig — i.e. the file load()
// would return — so save() must PRESERVE it and overwrite the other. A torn/truncated higher-seq
// file (intact 12-byte adapter header, short body) fails the length+checksum here, so it is NOT
// treated as the last-good file. Mirrors ofc/core/persistence.hpp + estimator.cpp deserialize().
bool blob_is_valid(const std::vector<unsigned char>& blob, std::uint64_t want_cfg_hash) {
    const int len = static_cast<int>(blob.size());
    const int header_len = 4 + 4 + 8 + 4;   // magic + version + config_hash + payload_len
    if (len < header_len + 4) return false;  // too short for header + trailing checksum
    const unsigned char* b = blob.data();
    // magic
    for (int i = 0; i < 4; ++i) {
        if (b[i] != persist::kMagic[i]) return false;
    }
    // version (LE u32 at offset 4)
    std::uint32_t version = 0;
    for (int i = 0; i < 4; ++i) version |= static_cast<std::uint32_t>(b[4 + i]) << (8 * i);
    if (version != persist::kFormatVersion) return false;
    // config_hash (LE u64 at offset 8) must match the rig we are saving
    std::uint64_t stored_hash = 0;
    for (int i = 0; i < 8; ++i) stored_hash |= static_cast<std::uint64_t>(b[8 + i]) << (8 * i);
    if (stored_hash != want_cfg_hash) return false;
    // payload_len (LE u32 at offset 16) -> exact total length
    std::uint32_t payload_len = 0;
    for (int i = 0; i < 4; ++i) payload_len |= static_cast<std::uint32_t>(b[16 + i]) << (8 * i);
    constexpr std::uint32_t kMaxPayloadLen = 0x40000000u;   // mirrors deserialize()'s cap
    if (payload_len > kMaxPayloadLen) return false;
    const long long need = static_cast<long long>(header_len) +
                           static_cast<long long>(payload_len) + 4;
    if (static_cast<long long>(len) != need) return false;   // truncated / torn / extra
    // checksum over header+payload (everything before the trailing 4-byte checksum word)
    const int covered = header_len + static_cast<int>(payload_len);
    const std::uint32_t want = persist::fnv1a32(b, covered);
    std::uint32_t got = 0;
    for (int i = 0; i < 4; ++i) got |= static_cast<std::uint32_t>(b[covered + i]) << (8 * i);
    return got == want;
}

// Write { seq, blob_len, blob } to `path`, truncating, then flush + close. Returns false on
// any IO failure/exception. (A real ECU adapter would fsync the fd here; std::ofstream::flush
// + dtor close is the std-portable analogue — it pushes the bytes to the OS before we flip.)
bool write_record(const std::string& path, std::uint64_t seq,
                  const unsigned char* blob, int n) {
    try {
        std::ofstream f(path, std::ios::binary | std::ios::trunc);
        if (!f) return false;
        unsigned char hdr[kHdrBytes];
        put_u64_le(hdr, seq);
        put_u32_le(hdr + kSeqBytes, static_cast<std::uint32_t>(n));
        f.write(reinterpret_cast<const char*>(hdr), kHdrBytes);
        if (n > 0) f.write(reinterpret_cast<const char*>(blob), n);
        f.flush();
        if (!f) return false;
        f.close();
        return static_cast<bool>(f);
    } catch (...) {
        return false;
    }
}

} // namespace

FilePersistence::FilePersistence(std::string path_a, std::string path_b, int max_blob_bytes)
    : path_a_(std::move(path_a)),
      path_b_(std::move(path_b)),
      max_blob_bytes_(max_blob_bytes > 0 ? max_blob_bytes : kDefaultMaxBlobBytes),
      last_seq_(0) {}

Status FilePersistence::save(const Estimator& est) {
    // Serialize into the staging buffer first — a serialize failure must NOT touch a file.
    std::vector<unsigned char> buf(static_cast<std::size_t>(max_blob_bytes_), 0);
    const Expected<int> ser = est.serialize(buf.data(), max_blob_bytes_);
    if (!ser.ok()) return ser.status();   // NotInitialized / CapacityExceeded straight through
    const int n = ser.value();

    // The config-hash of the estimator we are saving lives at offset 8 (LE u64) of the blob we
    // just serialized — it identifies the rig a candidate must match to be load()-acceptable.
    std::uint64_t want_cfg_hash = 0;
    for (int i = 0; i < 8; ++i)
        want_cfg_hash |= static_cast<std::uint64_t>(buf[8 + i]) << (8 * i);

    // Pick the INACTIVE target by VALIDITY, not by the raw on-disk seq word. The file to
    // PRESERVE is the one load() would return = the highest-seq candidate whose blob actually
    // passes the core framing guards (a complete, checksum + config-hash-valid OFCP blob). A
    // file with an intact 12-byte adapter header but a TORN body reads as present=true with a
    // (possibly higher) seq, yet is NOT valid — overwriting it is correct; overwriting the
    // genuinely-good lower-seq file (the old seq-word-only logic) could lose BOTH copies on a
    // crash during this very write. So: overwrite the INVALID/missing file; if both are valid,
    // overwrite the LOWER-seq one; if neither is valid, overwrite either (no good state to lose).
    const Record ra = read_record(path_a_);
    const Record rb = read_record(path_b_);
    const std::uint64_t sa = ra.present ? ra.seq : 0;
    const std::uint64_t sb = rb.present ? rb.seq : 0;
    const bool va = ra.present && blob_is_valid(ra.blob, want_cfg_hash);
    const bool vb = rb.present && blob_is_valid(rb.blob, want_cfg_hash);

    bool write_to_a;
    if (va && vb)       write_to_a = (sa <= sb);   // both good -> overwrite the older (preserve newest)
    else if (va)        write_to_a = false;        // only A good -> overwrite B (preserve A)
    else if (vb)        write_to_a = true;         // only B good -> overwrite A (preserve B)
    else                write_to_a = true;         // neither good -> overwrite either (-> A)

    // next_seq is one past the highest seq PRESENT on disk (valid or not) so load() — which
    // picks the highest-seq VALID record — chooses our fresh write next, even past a torn
    // higher-seq residue that we are about to overwrite or that sits on the preserved file.
    const std::uint64_t next_seq = (sa > sb ? sa : sb) + 1;

    const std::string& target = write_to_a ? path_a_ : path_b_;
    if (!write_record(target, next_seq, buf.data(), n)) return Status::CorruptData;

    last_seq_ = next_seq;
    return Status::Ok;
}

Status FilePersistence::load(Estimator& est) {
    const Record ra = read_record(path_a_);
    const Record rb = read_record(path_b_);

    // Order the two candidates by DESCENDING seq, then try the core deserialize on each;
    // the first the core ACCEPTS wins (so a corrupt higher-seq file falls back to the lower).
    const Record* hi = &ra;
    const Record* lo = &rb;
    if (rb.present && (!ra.present || rb.seq > ra.seq)) { hi = &rb; lo = &ra; }

    const Record* order[2] = { hi, lo };
    Status last_err = Status::NoData;
    for (const Record* c : order) {
        if (!c->present) continue;
        const Status st = est.deserialize(c->blob.data(), static_cast<int>(c->blob.size()));
        if (ok(st)) { last_seq_ = c->seq; return Status::Ok; }
        // NotInitialized means the caller's estimator is unusable — surface it directly rather
        // than masking it as NoData (deserialize would fail identically on the other file).
        if (st == Status::NotInitialized) return st;
        last_err = st;
    }
    (void)last_err;
    return Status::NoData;
}

} // namespace adapters
} // namespace ofc
