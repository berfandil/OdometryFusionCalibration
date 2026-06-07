// ofc_adapters/file_persistence.hpp — production file-persistence adapter (Slice 13, D23).
//
// RELAXED EDGE. The strict core's Estimator::serialize()/deserialize() touch ONLY a
// caller-owned fixed buffer (no heap / file IO / exceptions — see ofc/core/persistence.hpp).
// This adapter wraps those primitives in the production DOUBLE-BUFFER PING-PONG over two
// files A/B. save() always writes to the file load() would REJECT (the invalid / older / missing
// one) and never the highest-seq VALID file, so a crash mid-write leaves at most that
// rejected-anyway file torn and load() always recovers the last good state.
//
// THE ON-DISK RECORD (one per file): a small header then the core blob.
//   [ uint64 seq        : 8 bytes, little-endian — monotonic write sequence    ]
//   [ uint32 blob_len   : 4 bytes, little-endian — bytes of core blob below     ]
//   [ blob : blob_len bytes — the Estimator::serialize() "OFCP" blob            ]
// The core blob is itself versioned + checksummed + config-hash-guarded, so a torn /
// truncated / corrupt inactive file is rejected by the core deserialize() guards and load()
// falls back to the other file. blob_len lets load() detect a short write before handing the
// (possibly garbage) tail to the core.
//
// STATELESS ACROSS RESTARTS. save() reads BOTH files and picks the INACTIVE target by VALIDITY,
// not by the raw on-disk seq word: it checks whether each file holds a complete, checksum +
// config-hash-valid blob for the rig being saved (the same framing guards load()/the core apply),
// and OVERWRITES the file load() would reject — the invalid / torn / missing one; if both are
// valid it overwrites the LOWER-seq one. It NEVER overwrites the highest-seq VALID file (the one
// load() would return), then writes seq = max(seq)+1. This handles the crash residue the whole
// design exists to survive: a higher-seq file with an intact 12-byte adapter header but a TORN
// body reads as present with a higher seq YET is invalid, so it is the overwrite target and the
// genuinely-good lower-seq file is preserved (the seq-word-only target choice could clobber the
// good file and lose BOTH copies on a crash during that write). A fresh process resumes the
// ping-pong correctly without remembering which file it wrote last. load() reads both, tries the
// HIGHER-seq blob first, and falls back to the other if the higher one fails the core guards.
//
// EXCEPTIONS DO NOT ESCAPE. std::fstream IO can throw / fail; every public method catches and
// maps failures to a Status. The adapter never lets an exception cross its own API boundary.
#ifndef OFC_ADAPTERS_FILE_PERSISTENCE_HPP
#define OFC_ADAPTERS_FILE_PERSISTENCE_HPP

#include "ofc/core/estimator.hpp"
#include "ofc/core/status.hpp"

#include <cstdint>
#include <string>

namespace ofc {
namespace adapters {

// Default core-blob staging buffer size (bytes). The "OFCP" blob is small (header + a few
// per-source SE3s + checksum); 4 KiB is generous. save() returns CapacityExceeded if the
// blob does not fit (bump max_blob_bytes via the ctor for a very large rig).
constexpr int kDefaultMaxBlobBytes = 4096;

class FilePersistence {
public:
    // path_a / path_b are the two ping-pong files. They need not exist yet; save() creates
    // them. max_blob_bytes sizes the internal serialize staging buffer (must hold the core
    // blob); the default fits any realistic rig.
    FilePersistence(std::string path_a, std::string path_b,
                    int max_blob_bytes = kDefaultMaxBlobBytes);

    // Serialize `est` into the INACTIVE file (the one load() would reject — the invalid / torn /
    // older / missing one; never the highest-seq VALID file) with the next seq, flushing +
    // closing before returning so a crash leaves at most that rejected-anyway file torn.
    //   Ok               — written + flushed.
    //   NotInitialized   — est.serialize() reports the estimator was never init()'d.
    //   CapacityExceeded — the blob did not fit the staging buffer.
    //   CorruptData      — the file could not be opened / written (IO failure; mapped here
    //                      because the on-disk record is then unusable).
    Status save(const Estimator& est);

    // Read BOTH files, pick the HIGHEST seq whose blob the core ACCEPTS (deserialize == Ok),
    // falling back to the other.
    //   Ok          — a valid record was restored into `est`.
    //   NoData      — neither file holds a record the core accepts (missing / torn / corrupt
    //                 / config-mismatch on both).
    //   NotInitialized — `est` was not init()'d (the core deserialize rejects it).
    Status load(Estimator& est);

    // The seq of the record load() chose, or save() last wrote (0 before any IO). Lets a
    // caller observe the ping-pong (used by the crash-safety test to assert which file won).
    std::uint64_t last_seq() const { return last_seq_; }

private:
    std::string   path_a_;
    std::string   path_b_;
    int           max_blob_bytes_;
    std::uint64_t last_seq_;
};

} // namespace adapters
} // namespace ofc
#endif // OFC_ADAPTERS_FILE_PERSISTENCE_HPP
