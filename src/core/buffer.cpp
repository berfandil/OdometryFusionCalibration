// ofc/core/buffer.cpp — per-source ring buffer + uniform delta query (Slice 1).
//
// STRICT CORE: storage is sized once in configure(); the push helpers and query()
// allocate nothing. No exceptions; status-code returns. Every scan is bounded by
// the fixed capacity. double math, no fast-math.
//
// Normalization (D7, DESIGN §8): each native form is converted at ingest to a
// {stamp, cum_pose, incr_cov} entry. cum_pose is the cumulative pose in the
// source's own integrated frame; incr_cov is the covariance of the increment from
// the previous entry, in [trans; rot] order (matching Delta::cov).
//
// Covariance frame convention: incr_cov lives in the source's body-increment
// tangent frame and is accumulated additively across the query window (with the two
// end increments scaled by their temporal overlap fraction). This is a simple,
// consistent first-order model adequate for the median weighting this feeds; a
// fuller adjoint-transported accumulation is deferred to a later slice.
#include "ofc/core/buffer.hpp"

#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>

namespace ofc {

namespace {
// Seconds per Timestamp tick (Timestamp is nanoseconds — see types.hpp).
constexpr Scalar kNanosToSec = Scalar(1e-9);

Scalar dt_seconds(Timestamp a, Timestamp b) {
    return static_cast<Scalar>(b - a) * kNanosToSec;
}
} // namespace

Status SourceBuffer::configure(const SensorConfig& sensor, OdomForm form,
                               int capacity, ConfCombine combine,
                               Scalar confidence_blend) {
    if (capacity < 2) return Status::OutOfRange;
    switch (combine) {
        case ConfCombine::NativeOnly:
        case ConfCombine::ModeledOnly:
        case ConfCombine::Sum:
        case ConfCombine::Max:
        case ConfCombine::Weighted:
            break;
        default:
            return Status::InvalidConfig;
    }

    // Allocate once; never grow after (re-configure with a larger capacity is the
    // only path that reallocates, and it happens before any post-init steady state).
    if (capacity > capacity_) {
        buf_.assign(static_cast<std::size_t>(capacity), BufferEntry{});
    }
    capacity_   = capacity;
    head_       = 0;
    count_      = 0;
    configured_ = true;
    form_       = form;
    sensor_     = sensor;
    combine_    = combine;
    blend_      = std::max(Scalar(0), std::min(Scalar(1), confidence_blend));
    return Status::Ok;
}

void SourceBuffer::reset() {
    head_  = 0;
    count_ = 0;
}

const BufferEntry& SourceBuffer::at(int logical) const {
    const int phys = (head_ + logical) % capacity_;
    return buf_[static_cast<std::size_t>(phys)];
}

Status SourceBuffer::push_twist(Timestamp t, const Vec6& xi, const Mat6& cov) {
    if (!configured_)            return Status::NotInitialized;
    if (form_ != OdomForm::Twist) return Status::InvalidConfig;
    if (count_ > 0 && t <= newest()) return Status::OutOfRange;

    BufferEntry e;
    e.stamp = t;
    if (count_ == 0) {
        // First sample seeds the integrated frame at identity; no preceding
        // increment, so its covariance is zero.
        e.cum_pose = SE3{};
        e.incr_cov = Mat6::Zero();
    } else {
        const BufferEntry& prev = at(count_ - 1);
        const Scalar dt = dt_seconds(prev.stamp, t);
        const SE3 increment = se3::exp(xi * dt);
        e.cum_pose = se3::compose(prev.cum_pose, increment);
        e.incr_cov = cov * dt;   // integrate the per-second twist covariance
    }

    const int phys = (head_ + count_) % capacity_;
    buf_[static_cast<std::size_t>(phys)] = e;
    if (count_ < capacity_) {
        ++count_;
    } else {
        head_ = (head_ + 1) % capacity_;   // overwrote the oldest; advance head
    }
    return Status::Ok;
}

Status SourceBuffer::push_increment(Timestamp t, const SE3& incr, const Mat6& cov) {
    if (!configured_)                 return Status::NotInitialized;
    if (form_ != OdomForm::Increment) return Status::InvalidConfig;
    if (count_ > 0 && t <= newest())  return Status::OutOfRange;

    BufferEntry e;
    e.stamp = t;
    if (count_ == 0) {
        e.cum_pose = SE3{};
        e.incr_cov = Mat6::Zero();
    } else {
        const BufferEntry& prev = at(count_ - 1);
        e.cum_pose = se3::compose(prev.cum_pose, incr);
        e.incr_cov = cov;
    }

    const int phys = (head_ + count_) % capacity_;
    buf_[static_cast<std::size_t>(phys)] = e;
    if (count_ < capacity_) ++count_;
    else                    head_ = (head_ + 1) % capacity_;
    return Status::Ok;
}

Status SourceBuffer::push_absolute(Timestamp t, const SE3& pose, const Mat6& cov) {
    if (!configured_)                    return Status::NotInitialized;
    if (form_ != OdomForm::AbsolutePose) return Status::InvalidConfig;
    if (count_ > 0 && t <= newest())     return Status::OutOfRange;

    BufferEntry e;
    e.stamp    = t;
    e.cum_pose = pose;   // absolute pose is already the cumulative pose
    // Increment covariance choice (documented): absolute-pose providers report an
    // absolute pose covariance, not an increment covariance, and the difference of
    // two absolute covariances is not guaranteed PSD. We therefore take the supplied
    // covariance as the per-step increment covariance directly. The first entry has
    // no preceding increment, so its incr_cov is zero.
    e.incr_cov = (count_ == 0) ? Mat6::Zero() : cov;

    const int phys = (head_ + count_) % capacity_;
    buf_[static_cast<std::size_t>(phys)] = e;
    if (count_ < capacity_) ++count_;
    else                    head_ = (head_ + 1) % capacity_;
    return Status::Ok;
}

Timestamp SourceBuffer::oldest() const {
    return (count_ == 0) ? Timestamp(0) : at(0).stamp;
}

Timestamp SourceBuffer::newest() const {
    return (count_ == 0) ? Timestamp(0) : at(count_ - 1).stamp;
}

bool SourceBuffer::covers(Timestamp t0, Timestamp t1) const {
    if (t1 <= t0)     return false;
    if (count_ < 2)   return false;
    return t0 >= oldest() && t1 <= newest();
}

// Smallest logical index whose stamp >= t (bounded linear scan; capacity-bounded).
int SourceBuffer::bracket_upper(Timestamp t) const {
    for (int i = 0; i < count_; ++i) {
        if (at(i).stamp >= t) return i;
    }
    return count_;   // t is past the newest entry
}

SE3 SourceBuffer::pose_at(Timestamp t) const {
    const int hi = bracket_upper(t);   // first entry with stamp >= t
    if (hi == 0)        return at(0).cum_pose;             // at/just before oldest
    if (hi >= count_)   return at(count_ - 1).cum_pose;    // at/just after newest
    const BufferEntry& lo = at(hi - 1);
    const BufferEntry& up = at(hi);
    if (up.stamp == t)  return up.cum_pose;
    const Scalar span = dt_seconds(lo.stamp, up.stamp);
    const Scalar u = (span > Scalar(0)) ? dt_seconds(lo.stamp, t) / span : Scalar(0);
    return se3::interpolate(lo.cum_pose, up.cum_pose, u);
}

Mat6 SourceBuffer::accumulate_native_cov(Timestamp t0, Timestamp t1) const {
    Mat6 sum = Mat6::Zero();
    // Each entry i (i >= 1) carries the increment covariance over (stamp[i-1], stamp[i]].
    // Add the portion of each increment overlapping [t0, t1], scaled by the temporal
    // overlap fraction so the two end increments are partially counted.
    for (int i = 1; i < count_; ++i) {
        const Timestamp a = at(i - 1).stamp;
        const Timestamp b = at(i).stamp;
        const Timestamp lo = std::max(a, t0);
        const Timestamp hi = std::min(b, t1);
        if (hi <= lo) continue;                 // no overlap
        const Scalar full = dt_seconds(a, b);
        if (full <= Scalar(0)) continue;
        const Scalar frac = dt_seconds(lo, hi) / full;
        sum += at(i).incr_cov * frac;
    }
    return sum;
}

Mat6 SourceBuffer::modeled_cov(const SE3& motion) const {
    // Synthesize a diagonal cov from the motion magnitude and the per-sensor noise
    // coefficients (CONFIG §9). Translation axes scale with ||t||, rotation axes with
    // the rotation angle. A tiny floor keeps the matrix strictly positive-definite so
    // it can never make a combined covariance singular.
    const Scalar trans_mag = motion.t.norm();
    const Scalar rot_mag   = so3::log(motion.R).norm();
    const Scalar st = sensor_.modeled_noise_per_m   * trans_mag;
    const Scalar sr = sensor_.modeled_noise_per_rad * rot_mag;
    constexpr Scalar kFloor = Scalar(1e-12);
    Mat6 m = Mat6::Zero();
    const Scalar tvar = st * st + kFloor;
    const Scalar rvar = sr * sr + kFloor;
    m(0, 0) = tvar; m(1, 1) = tvar; m(2, 2) = tvar;
    m(3, 3) = rvar; m(4, 4) = rvar; m(5, 5) = rvar;
    return m;
}

Mat6 SourceBuffer::combine_cov(const Mat6& native, const Mat6& modeled,
                               bool have_native) const {
    // Per D7: no usable native covariance (provider disabled it, or none supplied)
    // -> fall back to the modeled covariance.
    if (!have_native || !sensor_.native_confidence) return modeled;

    switch (combine_) {
        case ConfCombine::NativeOnly:  return native;
        case ConfCombine::ModeledOnly: return modeled;
        case ConfCombine::Sum:         return native + modeled;
        case ConfCombine::Max:         return native.cwiseMax(modeled);
        case ConfCombine::Weighted:
            return blend_ * native + (Scalar(1) - blend_) * modeled;
    }
    return modeled;   // unreachable; defensive default
}

Expected<Delta> SourceBuffer::query(Timestamp t0, Timestamp t1) const {
    if (t1 <= t0)            return Expected<Delta>(Status::OutOfRange);
    if (!covers(t0, t1))     return Expected<Delta>(Status::NoData);

    const SE3 p0 = pose_at(t0);
    const SE3 p1 = pose_at(t1);

    Delta d;
    d.t0     = t0;
    d.t1     = t1;
    d.motion = se3::compose(se3::inverse(p0), p1);

    const Mat6 native  = accumulate_native_cov(t0, t1);
    const bool have_native = !native.isZero(Scalar(0));
    const Mat6 modeled = modeled_cov(d.motion);
    d.cov = combine_cov(native, modeled, have_native);

    return Expected<Delta>(d);
}

} // namespace ofc
