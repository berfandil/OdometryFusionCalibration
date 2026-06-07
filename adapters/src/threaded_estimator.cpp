// ofc_adapters/threaded_estimator.cpp — see threaded_estimator.hpp for the contract.
//
// RELAXED EDGE: std::thread / std::mutex / std::condition_variable / heap / exceptions fine.
// INVARIANT: exactly ONE thread (the worker) ever calls into est_; the snapshot mutex guards
// only the published Result copy. No exception escapes a public method or terminates the
// worker thread.
#include "ofc_adapters/threaded_estimator.hpp"

namespace ofc {
namespace adapters {

ThreadedEstimator::ThreadedEstimator(Estimator& est) : est_(est) {}

ThreadedEstimator::~ThreadedEstimator() {
    stop();   // signal + join (no-op if never started / already stopped)
}

Status ThreadedEstimator::start() {
    std::lock_guard<std::mutex> lk(q_mtx_);
    if (started_) return Status::Ok;            // idempotent
    started_   = true;
    stop_now_  = false;
    drain_then_stop_ = false;
    worker_ = std::thread(&ThreadedEstimator::worker_loop_, this);
    return Status::Ok;
}

Status ThreadedEstimator::submit(Timestamp now) {
    {
        std::lock_guard<std::mutex> lk(q_mtx_);
        if (stop_now_ || drain_then_stop_) return Status::Rejected;   // shutting down
        queue_.push_back(now);
    }
    q_cv_.notify_one();
    return Status::Ok;
}

Status ThreadedEstimator::submit_batch(const Timestamp* stamps, int count) {
    if (stamps == nullptr || count < 0) return Status::OutOfRange;
    {
        std::lock_guard<std::mutex> lk(q_mtx_);
        if (stop_now_ || drain_then_stop_) return Status::Rejected;
        for (int i = 0; i < count; ++i) queue_.push_back(stamps[i]);
    }
    q_cv_.notify_one();
    return Status::Ok;
}

Result ThreadedEstimator::snapshot() const {
    std::lock_guard<std::mutex> lk(snap_mtx_);
    return snapshot_;   // a COPY — tear-free
}

std::uint64_t ThreadedEstimator::steps_done() const {
    std::lock_guard<std::mutex> lk(q_mtx_);
    return steps_done_;
}

bool ThreadedEstimator::running() const {
    std::lock_guard<std::mutex> lk(q_mtx_);
    return started_ && !stop_now_;
}

Status ThreadedEstimator::drain_and_stop() {
    bool already_stopped;
    {
        std::lock_guard<std::mutex> lk(q_mtx_);
        // A prior stop() (which DROPS pending stamps) must stay honored: a later
        // drain_and_stop() must not resurrect those dropped stamps via the inline-drain branch
        // below. Capture it now so we can skip the inline drain when stop already happened.
        already_stopped = stop_now_;
        if (!started_) {
            // Never started: drain synchronously here (the caller may have submitted stamps
            // before start()). Process them inline, then mark stopped.
            // Fall through to inline drain below.
        }
        drain_then_stop_ = true;
    }
    q_cv_.notify_all();

    if (worker_.joinable()) {
        worker_.join();
    } else if (already_stopped) {
        // A prior stop() already dropped the pending queue's claim to be processed; honor the
        // drop. Clear the leftover stamps so a later snapshot/observation does not imply they
        // ran, then fall through to mark stopped. (No worker exists to drain them in any case.)
        std::lock_guard<std::mutex> lk(q_mtx_);
        queue_.clear();
    } else {
        // No worker (start() was never called): drain the queue inline on this thread so the
        // contract ("the final snapshot equals a single-threaded run over the stamps") still
        // holds. Still single-threaded into est_ (this is the only caller).
        for (;;) {
            Timestamp now;
            {
                std::lock_guard<std::mutex> lk(q_mtx_);
                if (queue_.empty()) break;
                now = queue_.front();
                queue_.pop_front();
            }
            try {
                if (ok(est_.step(now))) {
                    Result snap = est_.latest();
                    {
                        std::lock_guard<std::mutex> sl(snap_mtx_);
                        snapshot_ = snap;
                    }
                    std::lock_guard<std::mutex> lk(q_mtx_);
                    ++steps_done_;
                }
            } catch (...) {
                // swallow — a bad step must not abort the drain
            }
        }
    }

    std::lock_guard<std::mutex> lk(q_mtx_);
    stop_now_ = true;
    started_  = false;
    return Status::Ok;
}

Status ThreadedEstimator::stop() {
    {
        std::lock_guard<std::mutex> lk(q_mtx_);
        if (!started_ && !worker_.joinable()) { stop_now_ = true; return Status::Ok; }
        stop_now_ = true;
    }
    q_cv_.notify_all();
    if (worker_.joinable()) worker_.join();
    std::lock_guard<std::mutex> lk(q_mtx_);
    started_ = false;
    return Status::Ok;
}

void ThreadedEstimator::worker_loop_() {
    for (;;) {
        Timestamp now;
        {
            std::unique_lock<std::mutex> lk(q_mtx_);
            q_cv_.wait(lk, [&] {
                return stop_now_ || !queue_.empty() ||
                       (drain_then_stop_ && queue_.empty());
            });
            if (stop_now_) return;                                // immediate stop
            if (queue_.empty()) {
                if (drain_then_stop_) return;                     // drained -> exit
                continue;                                         // spurious wake
            }
            now = queue_.front();
            queue_.pop_front();
        }

        // Single-threaded into est_ (only this thread ever calls it). Publish a snapshot of
        // latest() under the snapshot mutex after a successful step.
        try {
            const Status st = est_.step(now);
            if (ok(st)) {
                Result snap = est_.latest();        // copy under no lock (worker-exclusive est_)
                {
                    std::lock_guard<std::mutex> sl(snap_mtx_);
                    snapshot_ = snap;               // tear-free publish
                }
                std::lock_guard<std::mutex> lk(q_mtx_);
                ++steps_done_;
            }
        } catch (...) {
            // A throw must not terminate the worker thread / process. Swallow + continue.
        }
    }
}

} // namespace adapters
} // namespace ofc
