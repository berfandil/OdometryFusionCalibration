// ofc/core/persistence.hpp — warm-restart serialization-format primitives (Slice 12, D23).
//
// STRICT CORE. These are the explicit-endian, padding-free byte read/write helpers and the
// FNV-1a hash the Estimator's serialize()/deserialize() use to lay down a STABLE, portable
// on-the-wire format into a CALLER-OWNED fixed buffer. No heap, no exceptions, no file I/O,
// bounded — everything operates on (buf, cap, cursor) with an explicit overflow guard. File
// I/O + the double-buffer ping-pong live RELAXED-EDGE (an adapter / the tests), on top of
// these core primitives.
//
// WHY EXPLICIT BYTES (not memcpy of a struct). A struct memcpy bakes in the compiler's
// padding + the host's endianness, so a blob written by one build would not deserialize on
// another. We instead serialize every field little-endian, byte-by-byte for the integers and
// via the IEEE-754 bit pattern for doubles (the codebase is `double` everywhere, no
// fast-math, so the bit pattern is the value). The result is a fixed, documented layout that
// round-trips bit-exact and is stable across builds of the same ISA family.
//
// THE BLOB LAYOUT (see estimator.cpp serialize() for the payload schema):
//   [ "OFCP" magic : 4 bytes                                            ]
//   [ uint32 format_version (little-endian)                             ]
//   [ uint64 config_hash    (little-endian, FNV-1a over the rig config) ]
//   [ uint32 payload_len    (little-endian, bytes of payload below)     ]
//   [ payload : payload_len bytes (the committed calibration state)     ]
//   [ uint32 checksum       (little-endian, FNV-1a over header+payload) ]
// The checksum covers the magic+version+config_hash+payload_len+payload (everything BEFORE
// the trailing checksum word) so a single flipped byte anywhere is detected on read.
#ifndef OFC_CORE_PERSISTENCE_HPP
#define OFC_CORE_PERSISTENCE_HPP

#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include <cstdint>
#include <cstring>

namespace ofc {
namespace persist {

// On-the-wire format identity. Bump kFormatVersion on any incompatible payload-schema change
// (deserialize rejects a mismatch with Status::VersionMismatch).
constexpr unsigned char kMagic[4]      = { 'O', 'F', 'C', 'P' };
constexpr std::uint32_t kFormatVersion = 1u;

// FNV-1a (64-bit) over a byte span — the config_hash hash family. Deterministic, no deps,
// strict-core (a plain loop, no heap). Stream more bytes by threading the returned state back
// in as `h` (fnv1a64(more, n, fnv1a64(first, m))).
constexpr std::uint64_t kFnvOffset64 = 1469598103934665603ull;
constexpr std::uint64_t kFnvPrime64  = 1099511628211ull;

inline std::uint64_t fnv1a64(const void* data, int n, std::uint64_t h = kFnvOffset64) {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (int i = 0; i < n; ++i) {
        h ^= static_cast<std::uint64_t>(p[i]);
        h *= kFnvPrime64;
    }
    return h;
}

// FNV-1a (32-bit) — the trailing blob checksum. Same construction, 32-bit width.
constexpr std::uint32_t kFnvOffset32 = 2166136261u;
constexpr std::uint32_t kFnvPrime32  = 16777619u;

inline std::uint32_t fnv1a32(const void* data, int n) {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    std::uint32_t h = kFnvOffset32;
    for (int i = 0; i < n; ++i) {
        h ^= static_cast<std::uint32_t>(p[i]);
        h *= kFnvPrime32;
    }
    return h;
}

// ---- Bounded little-endian writer over a caller-owned fixed buffer --------------------
// Every put_* checks capacity and trips `overflow` (then becomes a no-op) instead of writing
// out of bounds — the strict-core "no exceptions, status-by-flag" idiom. The caller checks
// .overflow once at the end and returns CapacityExceeded.
struct Writer {
    unsigned char* buf = nullptr;
    int            cap = 0;
    int            pos = 0;
    bool           overflow = false;

    Writer(unsigned char* b, int c) : buf(b), cap(c) {}

    bool fits(int n) const { return pos + n <= cap; }

    void put_bytes(const void* data, int n) {
        if (!fits(n)) { overflow = true; return; }
        std::memcpy(buf + pos, data, static_cast<std::size_t>(n));
        pos += n;
    }
    void put_u32(std::uint32_t v) {
        unsigned char b[4];
        for (int i = 0; i < 4; ++i) b[i] = static_cast<unsigned char>((v >> (8 * i)) & 0xFFu);
        put_bytes(b, 4);
    }
    void put_u64(std::uint64_t v) {
        unsigned char b[8];
        for (int i = 0; i < 8; ++i) b[i] = static_cast<unsigned char>((v >> (8 * i)) & 0xFFu);
        put_bytes(b, 8);
    }
    void put_i32(std::int32_t v) { put_u32(static_cast<std::uint32_t>(v)); }
    void put_u8(unsigned char v) { put_bytes(&v, 1); }
    void put_bool(bool v)        { put_u8(v ? 1u : 0u); }
    // Doubles by their IEEE-754 bit pattern (the codebase is `double`, no fast-math, so the
    // pattern IS the value and round-trips bit-exact).
    void put_f64(Scalar v) {
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        put_u64(bits);
    }
    // An SE3 as 9 rotation doubles (row-major) + 3 translation doubles (explicit, no padding).
    void put_se3(const SE3& X) {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) put_f64(X.R(r, c));
        for (int i = 0; i < 3; ++i) put_f64(X.t(i));
    }
};

// ---- Bounded little-endian reader over a caller-owned fixed buffer --------------------
// Mirror of Writer. Each get_* checks remaining length and trips `underflow` (returning a
// zero value) rather than reading past `len`; the caller checks .underflow once and rejects
// the blob as CorruptData (a truncated / too-short blob).
struct Reader {
    const unsigned char* buf = nullptr;
    int                  len = 0;
    int                  pos = 0;
    bool                 underflow = false;

    Reader(const unsigned char* b, int l) : buf(b), len(l) {}

    bool has(int n) const { return pos + n <= len; }

    void get_bytes(void* out, int n) {
        if (!has(n)) { underflow = true; std::memset(out, 0, static_cast<std::size_t>(n)); return; }
        std::memcpy(out, buf + pos, static_cast<std::size_t>(n));
        pos += n;
    }
    std::uint32_t get_u32() {
        unsigned char b[4] = {0,0,0,0};
        get_bytes(b, 4);
        std::uint32_t v = 0;
        for (int i = 0; i < 4; ++i) v |= static_cast<std::uint32_t>(b[i]) << (8 * i);
        return v;
    }
    std::uint64_t get_u64() {
        unsigned char b[8] = {0,0,0,0,0,0,0,0};
        get_bytes(b, 8);
        std::uint64_t v = 0;
        for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(b[i]) << (8 * i);
        return v;
    }
    std::int32_t  get_i32()  { return static_cast<std::int32_t>(get_u32()); }
    unsigned char get_u8()   { unsigned char v = 0; get_bytes(&v, 1); return v; }
    bool          get_bool() { return get_u8() != 0u; }
    Scalar        get_f64() {
        std::uint64_t bits = get_u64();
        Scalar v;
        std::memcpy(&v, &bits, sizeof(v));
        return v;
    }
    SE3 get_se3() {
        SE3 X;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) X.R(r, c) = get_f64();
        for (int i = 0; i < 3; ++i) X.t(i) = get_f64();
        return X;
    }
};

} // namespace persist
} // namespace ofc
#endif // OFC_CORE_PERSISTENCE_HPP
