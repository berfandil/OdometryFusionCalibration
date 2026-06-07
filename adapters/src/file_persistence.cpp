// ofc_adapters/file_persistence.cpp — see file_persistence.hpp for the contract.
//
// RELAXED EDGE: std::fstream / std::vector / heap / exceptions are fine. The one hard rule
// is that NO exception escapes a public method — IO failures are caught and mapped to Status.
#include "ofc_adapters/file_persistence.hpp"

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

    // Pick the INACTIVE target stateless-ly: read both seqs, write to the LOWER-seq (or
    // missing) file with seq = max+1. A missing file counts as seq 0 (so it is chosen first).
    const Record ra = read_record(path_a_);
    const Record rb = read_record(path_b_);
    const std::uint64_t sa = ra.present ? ra.seq : 0;
    const std::uint64_t sb = rb.present ? rb.seq : 0;

    // Inactive = the file we will OVERWRITE. Prefer a missing file; else the lower seq; ties
    // -> A (deterministic). next_seq is one past the highest seq present.
    bool write_to_a;
    if (!ra.present)      write_to_a = true;
    else if (!rb.present) write_to_a = false;
    else                  write_to_a = (sa <= sb);   // overwrite the older (lower-seq) file
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
