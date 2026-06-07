// Adapters Slice 13: FilePersistence double-buffer ping-pong (productionizes the Slice-12
// test-side logic into ofc_adapters). Covers the done-criteria: round-trip restore, crash-mid-
// write fallback to the last good file, and config-change rejection.
#include <doctest/doctest.h>

#include "ofc_adapters/file_persistence.hpp"

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/status.hpp"

#include "adapter_test_util.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace ofc;
using namespace adptest;

namespace {

std::string temp_path(const char* tag) {
    const char* base = std::getenv("TEMP");
    if (base == nullptr) base = std::getenv("TMP");
    std::string dir = (base != nullptr) ? base : ".";
    return dir + "/ofc_adapt_persist_" + tag;
}

// Init an estimator over the rig + step it a few times so it has real state to serialize.
void make_and_run(Estimator& est, const Config& cfg, std::vector<TwistSource>& srcs) {
    REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);
    const Timestamp tick = secs_to_ns(0.02);
    for (Timestamp now = secs_to_ns(0.2); now <= secs_to_ns(1.5); now += tick) (void)est.step(now);
}

} // namespace

TEST_CASE("FilePersistence round-trip: save then load restores into a fresh estimator") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);

    const std::string pa = temp_path("rt_A.bin");
    const std::string pb = temp_path("rt_B.bin");
    std::remove(pa.c_str());
    std::remove(pb.c_str());

    // (1) Run + save.
    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; make_and_run(est, cfg, srcs);

    adapters::FilePersistence db(pa, pb);
    CHECK(db.save(est) == Status::Ok);
    CHECK(db.last_seq() == 1);                       // first write -> seq 1

    // The serialized blob of the SAME estimator (to compare a re-serialize after load).
    std::vector<unsigned char> ref(1024, 0);
    const Expected<int> r = est.serialize(ref.data(), static_cast<int>(ref.size()));
    REQUIRE(r.ok());
    ref.resize(r.value());

    // (2) Load into a FRESH estimator with the same config.
    std::vector<TwistSource> srcs2; make_rig_sources(srcs2);
    Estimator est2;
    REQUIRE(est2.init(cfg) == Status::Ok);
    for (auto& s : srcs2) REQUIRE(est2.add_source(&s) == Status::Ok);

    adapters::FilePersistence db2(pa, pb);
    CHECK(db2.load(est2) == Status::Ok);
    CHECK(db2.last_seq() == 1);

    // The restored estimator re-serializes to the SAME blob (idempotent restore).
    std::vector<unsigned char> got(1024, 0);
    const Expected<int> r2 = est2.serialize(got.data(), static_cast<int>(got.size()));
    REQUIRE(r2.ok());
    got.resize(r2.value());
    REQUIRE(got.size() == ref.size());
    CHECK(std::memcmp(got.data(), ref.data(), ref.size()) == 0);

    std::remove(pa.c_str());
    std::remove(pb.c_str());
}

TEST_CASE("FilePersistence ping-pong: alternates files + picks the inactive (lower-seq) target") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);

    const std::string pa = temp_path("pp_A.bin");
    const std::string pb = temp_path("pp_B.bin");
    std::remove(pa.c_str());
    std::remove(pb.c_str());

    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; make_and_run(est, cfg, srcs);

    adapters::FilePersistence db(pa, pb);
    CHECK(db.save(est) == Status::Ok); CHECK(db.last_seq() == 1);   // A (both missing -> A)
    CHECK(db.save(est) == Status::Ok); CHECK(db.last_seq() == 2);   // B (A=1,B=missing -> B)
    CHECK(db.save(est) == Status::Ok); CHECK(db.last_seq() == 3);   // A (overwrite lower seq)
    CHECK(db.save(est) == Status::Ok); CHECK(db.last_seq() == 4);   // B

    // A fresh adapter (stateless across "restart") still loads the highest seq.
    adapters::FilePersistence db_restart(pa, pb);
    std::vector<TwistSource> s2; make_rig_sources(s2);
    Estimator e2; REQUIRE(e2.init(cfg) == Status::Ok);
    for (auto& s : s2) REQUIRE(e2.add_source(&s) == Status::Ok);
    CHECK(db_restart.load(e2) == Status::Ok);
    CHECK(db_restart.last_seq() == 4);            // chose the highest seq, no in-memory state

    // And a save AFTER restart continues the ping-pong from the on-disk seqs (writes seq 5).
    CHECK(db_restart.save(e2) == Status::Ok);
    CHECK(db_restart.last_seq() == 5);

    std::remove(pa.c_str());
    std::remove(pb.c_str());
}

TEST_CASE("FilePersistence crash-mid-write: a torn higher-seq file falls back to the last good") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);

    const std::string pa = temp_path("crash_A.bin");
    const std::string pb = temp_path("crash_B.bin");
    std::remove(pa.c_str());
    std::remove(pb.c_str());

    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; make_and_run(est, cfg, srcs);

    // (1) A clean save -> file A, seq 1.
    adapters::FilePersistence db(pa, pb);
    REQUIRE(db.save(est) == Status::Ok);
    REQUIRE(db.last_seq() == 1);

    // (2) Simulate a crash mid-write to B: write a record with a HIGHER seq (2) but a TORN
    // (truncated) payload by hand. The header's blob_len claims the full length; the file is
    // short, so the core deserialize rejects it and load() must fall back to A.
    {
        // Read A's good blob back out (skip the 12-byte adapter header) to make a half-blob.
        std::ifstream fa(pa, std::ios::binary | std::ios::ate);
        REQUIRE(static_cast<bool>(fa));
        const std::streamoff sz = fa.tellg();
        fa.seekg(12);   // skip seq(8)+len(4)
        std::vector<unsigned char> good(static_cast<std::size_t>(sz) - 12);
        fa.read(reinterpret_cast<char*>(good.data()), static_cast<std::streamsize>(good.size()));

        std::ofstream fb(pb, std::ios::binary | std::ios::trunc);
        std::uint64_t seq = 2;                       // HIGHER than A's 1
        std::uint32_t len = static_cast<std::uint32_t>(good.size());  // claim full length...
        fb.write(reinterpret_cast<const char*>(&seq), 8);
        fb.write(reinterpret_cast<const char*>(&len), 4);
        fb.write(reinterpret_cast<const char*>(good.data()),
                 static_cast<std::streamsize>(good.size() / 2));       // ...but write only half
        fb.flush();
    }

    // load() sees B's higher seq first, the core rejects the torn blob, it falls back to A.
    std::vector<TwistSource> s2; make_rig_sources(s2);
    Estimator e2; REQUIRE(e2.init(cfg) == Status::Ok);
    for (auto& s : s2) REQUIRE(e2.add_source(&s) == Status::Ok);
    adapters::FilePersistence db2(pa, pb);
    CHECK(db2.load(e2) == Status::Ok);
    CHECK(db2.last_seq() == 1);                       // chose A (the last good), not torn B

    // (3) A clean save now targets the INACTIVE file by VALIDITY: A=1 is the only VALID record
    // (B=2 has an intact header but a torn body), so save() PRESERVES A and overwrites B with
    // next_seq = max(1,2)+1 = 3. (See the dedicated "overwrites the TORN higher-seq file"
    // test below for the crash-safety assertion that A is left untouched.)
    REQUIRE(db2.save(e2) == Status::Ok);
    CHECK(db2.last_seq() == 3);
    // A subsequent load now returns the clean seq-3 record.
    std::vector<TwistSource> s3; make_rig_sources(s3);
    Estimator e3; REQUIRE(e3.init(cfg) == Status::Ok);
    for (auto& s : s3) REQUIRE(e3.add_source(&s) == Status::Ok);
    adapters::FilePersistence db3(pa, pb);
    CHECK(db3.load(e3) == Status::Ok);
    CHECK(db3.last_seq() == 3);

    std::remove(pa.c_str());
    std::remove(pb.c_str());
}

namespace {

// Read the full on-disk record bytes (adapter header + body) of `path`.
std::vector<unsigned char> read_all(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) return {};
    const std::streamoff sz = f.tellg();
    f.seekg(0);
    std::vector<unsigned char> v(static_cast<std::size_t>(sz));
    if (!v.empty()) f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size()));
    return v;
}

// Hand-write a TORN record into `path`: seq + a blob_len that CLAIMS `claim_len` bytes, but only
// `body_bytes` of the (valid) body `good` actually written -> an intact 12-byte header over a
// short body (the realistic crash-mid-write residue).
void write_torn(const std::string& path, std::uint64_t seq,
                const std::vector<unsigned char>& good, std::size_t body_bytes) {
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    std::uint32_t len = static_cast<std::uint32_t>(good.size());   // claim the full length
    f.write(reinterpret_cast<const char*>(&seq), 8);
    f.write(reinterpret_cast<const char*>(&len), 4);
    if (body_bytes > 0)
        f.write(reinterpret_cast<const char*>(good.data()),
                static_cast<std::streamsize>(body_bytes));         // ... but write only a prefix
    f.flush();
}

} // namespace

// CRITICAL (Slice-13 review): save() must NOT overwrite the genuinely last-good file when the
// OTHER file has a HIGHER seq but a TORN body. Pre-fix, save() picked the overwrite target by the
// raw on-disk seq word alone (overwrite the lower seq), so it clobbered the VALID lower-seq file
// A and a crash during THAT write would lose BOTH copies. Post-fix, save() picks by VALIDITY:
// A (valid, seq 1) is preserved and the torn B (seq 2) is overwritten with seq 3. This test FAILS
// on the pre-fix code (which would leave A holding seq 3 and lose the original good A blob).
TEST_CASE("FilePersistence crash-safety: save() overwrites the TORN higher-seq file, not the good one") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);

    const std::string pa = temp_path("torn_A.bin");
    const std::string pb = temp_path("torn_B.bin");
    std::remove(pa.c_str());
    std::remove(pb.c_str());

    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; make_and_run(est, cfg, srcs);

    // (1) A clean save -> A holds a VALID seq-1 record. Snapshot A's exact bytes to prove later
    // that save() left A byte-for-byte untouched.
    adapters::FilePersistence db(pa, pb);
    REQUIRE(db.save(est) == Status::Ok);
    REQUIRE(db.last_seq() == 1);
    const std::vector<unsigned char> a_before = read_all(pa);
    REQUIRE_FALSE(a_before.empty());

    // (2) Crash-mid-write residue in B: a HIGHER seq (2) with an intact 12-byte header but a TORN
    // (half-written) body. read_record() reports it present with seq 2; the core would reject it.
    {
        std::vector<unsigned char> good(a_before.begin() + 12, a_before.end());  // A's valid blob
        REQUIRE(good.size() >= 2);
        write_torn(pb, /*seq=*/2, good, /*body_bytes=*/good.size() / 2);          // only half
    }

    // (3) save() must PRESERVE the valid A and overwrite the torn B with the next seq (3).
    adapters::FilePersistence db2(pa, pb);
    REQUIRE(db2.save(est) == Status::Ok);
    CHECK(db2.last_seq() == 3);

    // A is byte-for-byte the ORIGINAL valid seq-1 record (NOT overwritten).
    const std::vector<unsigned char> a_after = read_all(pa);
    REQUIRE(a_after.size() == a_before.size());
    CHECK(std::memcmp(a_after.data(), a_before.data(), a_before.size()) == 0);

    // B now holds the fresh seq-3 record, and a load picks it (highest-seq valid).
    {
        std::vector<TwistSource> sl; make_rig_sources(sl);
        Estimator el; REQUIRE(el.init(cfg) == Status::Ok);
        for (auto& s : sl) REQUIRE(el.add_source(&s) == Status::Ok);
        adapters::FilePersistence dl(pa, pb);
        CHECK(dl.load(el) == Status::Ok);
        CHECK(dl.last_seq() == 3);                    // the new B record
    }

    // (4) Simulate a crash DURING that seq-3 write to B: tear B's body again. The invariant holds
    // only because save() preserved A — load() falls back to A's good seq-1 state.
    {
        const std::vector<unsigned char> b_now = read_all(pb);
        REQUIRE(b_now.size() > 12);
        std::vector<unsigned char> b_body(b_now.begin() + 12, b_now.end());
        write_torn(pb, /*seq=*/3, b_body, /*body_bytes=*/b_body.size() / 2);   // tear seq-3 B
    }
    {
        std::vector<TwistSource> sr; make_rig_sources(sr);
        Estimator er; REQUIRE(er.init(cfg) == Status::Ok);
        for (auto& s : sr) REQUIRE(er.add_source(&s) == Status::Ok);
        adapters::FilePersistence dr(pa, pb);
        CHECK(dr.load(er) == Status::Ok);            // recovered, NOT NoData
        CHECK(dr.last_seq() == 1);                    // fell back to the preserved good A
    }

    std::remove(pa.c_str());
    std::remove(pb.c_str());
}

TEST_CASE("FilePersistence config-change: a blob from a different rig is rejected") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);

    const std::string pa = temp_path("cfg_A.bin");
    const std::string pb = temp_path("cfg_B.bin");
    std::remove(pa.c_str());
    std::remove(pb.c_str());

    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; make_and_run(est, cfg, srcs);
    adapters::FilePersistence db(pa, pb);
    REQUIRE(db.save(est) == Status::Ok);

    // A DIFFERENT rig (moved per-sensor prior_scale) -> different config_hash -> the core
    // deserialize rejects BOTH candidates -> load() returns NoData.
    std::vector<SensorConfig> sensors2;
    Config cfg2 = make_rig_cfg(sensors2);
    sensors2[1].prior_scale = 1.5;                   // was 1.0 -> config_hash differs
    cfg2.sensors = sensors2.data();

    std::vector<TwistSource> s2; make_rig_sources(s2);
    Estimator e2; REQUIRE(e2.init(cfg2) == Status::Ok);
    for (auto& s : s2) REQUIRE(e2.add_source(&s) == Status::Ok);
    adapters::FilePersistence db2(pa, pb);
    CHECK(db2.load(e2) == Status::NoData);           // both blobs rejected by the config guard

    std::remove(pa.c_str());
    std::remove(pb.c_str());
}

TEST_CASE("FilePersistence error paths: missing files load NoData; tiny buffer CapacityExceeded") {
    std::vector<SensorConfig> sensors;
    const Config cfg = make_rig_cfg(sensors);

    const std::string pa = temp_path("err_A.bin");
    const std::string pb = temp_path("err_B.bin");
    std::remove(pa.c_str());
    std::remove(pb.c_str());

    // Loading when NEITHER file exists -> NoData (not a crash).
    std::vector<TwistSource> srcs; make_rig_sources(srcs);
    Estimator est; REQUIRE(est.init(cfg) == Status::Ok);
    for (auto& s : srcs) REQUIRE(est.add_source(&s) == Status::Ok);
    adapters::FilePersistence db(pa, pb);
    CHECK(db.load(est) == Status::NoData);

    // A too-small staging buffer -> save() surfaces the core CapacityExceeded (no file written).
    adapters::FilePersistence tiny(pa, pb, /*max_blob_bytes=*/8);
    CHECK(tiny.save(est) == Status::CapacityExceeded);
    {
        std::ifstream fa(pa, std::ios::binary);
        CHECK_FALSE(static_cast<bool>(fa));          // nothing written on the serialize failure
    }

    // serialize/load before init -> NotInitialized surfaces through the adapter.
    Estimator un;
    adapters::FilePersistence db_un(pa, pb);
    CHECK(db_un.save(un) == Status::NotInitialized);

    std::remove(pa.c_str());
    std::remove(pb.c_str());
}
