// ofc/core/status.hpp — error handling for the strict core (no exceptions).
// AUTOSAR C++14: no std::optional / std::variant (C++17) → hand-rolled here.
#ifndef OFC_CORE_STATUS_HPP
#define OFC_CORE_STATUS_HPP

#include <utility>

namespace ofc {

enum class Status {
    Ok = 0,
    InvalidConfig,
    OutOfRange,
    CapacityExceeded,
    NotInitialized,
    NoData,
    NotReady,
    Unobservable,
    Rejected,        // e.g. Mahalanobis gate
};

constexpr bool ok(Status s) { return s == Status::Ok; }

// Minimal Optional (no dynamic alloc, trivially destructible payloads).
template <typename T>
class Optional {
public:
    Optional() : has_(false) {}
    explicit Optional(const T& v) : has_(true), val_(v) {}
    bool has_value() const { return has_; }
    const T& value() const { return val_; }      // precondition: has_value()
    T value_or(const T& d) const { return has_ ? val_ : d; }
private:
    bool has_;
    T    val_{};
};

// Minimal Expected: a value or a Status. No exceptions, no heap.
template <typename T>
class Expected {
public:
    Expected(const T& v) : status_(Status::Ok), val_(v) {}
    Expected(Status s)   : status_(s), val_() {}
    bool ok() const { return status_ == Status::Ok; }
    Status status() const { return status_; }
    const T& value() const { return val_; }       // precondition: ok()
private:
    Status status_;
    T      val_{};
};

} // namespace ofc
#endif // OFC_CORE_STATUS_HPP
