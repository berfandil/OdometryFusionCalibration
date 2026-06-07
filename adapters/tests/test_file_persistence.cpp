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

    // (3) A clean save now targets the INACTIVE file. A=1 (valid), B=2 (torn but seq 2 reads as
    // the higher), so next_seq = 3 written to A (the lower valid... actually A=1 < B=2 -> A).
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
