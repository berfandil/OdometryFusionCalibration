// ofc_adapters/threaded_estimator.hpp — threading wrapper for the caller-pumped core (Slice 13, D14).
//
// RELAXED EDGE. The strict core Estimator is SINGLE-THREADED + caller-pumped (D14): the
// consumer owns the thread and pumps step(). This adapter provides the "library owns the
// thread" ergonomic that middleware (ROS spinners, DDS) expect, WITHOUT touching the core's
// threading model: exactly ONE worker thread ever calls into the Estimator, so the core is
// still only ever pumped single-threaded. A mutex guards only the PUBLISHED Result copy that
// the worker hands to consumer threads — never the Estimator itself.
//
// FEED MODEL: timestamps are fed via a queue (the caller injects the clock), NOT a wall-clock
// timer — so the wrapper is deterministic + testable (the same stamps produce the same final
// snapshot as a single-threaded reference run). The worker pops a stamp, calls step(stamp),
// and on a successful step publishes a COPY of latest() under the snapshot mutex.
//
// LIFECYCLE: construct (not running) -> start() spawns the worker -> submit()/submit_batch()
// enqueue stamps -> snapshot() reads a tear-free Result copy from any thread -> stop() (or the
// dtor) signals drain+exit and JOINS the worker. drain_and_stop() blocks until the queue is
// empty and the last step published, then stops (the deterministic "process these N stamps
// then stop" mode the determinism test uses).
//
// EXCEPTIONS DO NOT ESCAPE: the worker swallows any exception from step()/latest() (the
// strict core does not throw, but Eigen under MSVC can in pathological cases) so a stray throw
// cannot terminate the process from the worker thread; the public API returns Status.
#ifndef OFC_ADAPTERS_THREADED_ESTIMATOR_HPP
#define OFC_ADAPTERS_THREADED_ESTIMATOR_HPP

#include "ofc/core/estimator.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <mutex>
#include <thread>

namespace ofc {
namespace adapters {

class ThreadedEstimator {
public:
    // Wraps an ALREADY-init()'d estimator (caller owns it; it must outlive this wrapper and
    // must NOT be touched by any other thread once start() is called). add_source/init are the
    // caller's job BEFORE construction — those are setup-time, single-threaded operations.
    explicit ThreadedEstimator(Estimator& est);
    ~ThreadedEstimator();   // stops + joins the worker

    ThreadedEstimator(const ThreadedEstimator&) = delete;
    ThreadedEstimator& operator=(const ThreadedEstimator&) = delete;

    // Spawn the worker thread. Idempotent: a second start() while running is a no-op (Ok).
    Status start();

    // Enqueue one / many timestamps for the worker to step() through (in order). Stamps fed
    // after stop()/drain_and_stop() are ignored (Rejected). Submitting before start() just
    // buffers them; start() then processes them.
    Status submit(Timestamp now);
    Status submit_batch(const Timestamp* stamps, int count);

    // A tear-free COPY of the latest published Result (the worker copies latest() under the
    // snapshot mutex after each successful step). Before the first successful step this is the
    // default-constructed Result (Phase::Init). Safe to call from any thread, any time.
    Result snapshot() const;

    // How many stamps the worker has successfully step()'d (Ok return) since start. Lets a
    // consumer poll progress without racing on the Result. Monotonic.
    std::uint64_t steps_done() const;

    // Block until the queue is EMPTY (every submitted stamp has been step()'d), then stop +
    // join. The deterministic batch mode: submit_batch(...) then drain_and_stop() leaves the
    // final snapshot equal to a single-threaded run over the same stamps. Returns Ok.
    Status drain_and_stop();

    // Signal the worker to finish the CURRENT step then exit, and JOIN it (does NOT wait for
    // the queue to drain — pending stamps are dropped). Idempotent. Called by the dtor.
    Status stop();

    bool running() const;

private:
    void worker_loop_();

    Estimator& est_;   // not owned

    // --- Work queue + lifecycle (queue mutex) ---
    mutable std::mutex      q_mtx_;
    std::condition_variable q_cv_;
    std::deque<Timestamp>   queue_;
    bool                    started_ = false;
    bool                    stop_now_ = false;     // stop(): exit without draining
    bool                    drain_then_stop_ = false;  // drain_and_stop(): exit when empty
    std::uint64_t           steps_done_ = 0;

    // --- Published snapshot (separate mutex so a consumer's snapshot() read never blocks the
    // worker on the queue mutex, and vice-versa) ---
    mutable std::mutex      snap_mtx_;
    Result                  snapshot_;

    std::thread             worker_;
};

} // namespace adapters
} // namespace ofc
#endif // OFC_ADAPTERS_THREADED_ESTIMATOR_HPP
