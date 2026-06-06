// ofc/core/buffer.hpp — per-source fixed-capacity ring buffer + uniform delta query.
// Slice 1 (D7, DESIGN §8): ingest one of three native odometry forms, normalize at
// ingest, and answer delta(t0,t1) -> (SE3 motion, Sigma 6x6) with in-window
// interpolation and the native-(+)-modeled covariance combine.
//
// STRICT CORE: storage is sized once in configure() and never reallocated; query()
// and the push helpers allocate nothing. No exceptions; status-code returns; bounded
// loops (every scan is bounded by the fixed capacity).
#ifndef OFC_CORE_BUFFER_HPP
#define OFC_CORE_BUFFER_HPP

#include "ofc/core/config.hpp"
#include "ofc/core/status.hpp"
#include "ofc/core/types.hpp"

#include <vector>

namespace ofc {

// The native form a source reports. Each is converted to a normalized entry at
// ingest (see SourceBuffer push helpers).
enum class OdomForm { Twist, Increment, AbsolutePose };

// A normalized stream entry. cum_pose is the cumulative pose in the source's own
// integrated frame (first entry = identity for relative forms, or the reported
// pose for AbsolutePose); incr_cov is the covariance of the increment from the
// previous entry, in [trans; rot] order (matching Delta::cov). First entry's
// incr_cov is zero (no preceding increment).
struct BufferEntry {
    Timestamp stamp    = Timestamp(0);
    SE3       cum_pose = SE3{};          // identity by default
    Mat6      incr_cov = Mat6::Zero();   // cov of increment (prev -> this)
};

// Fixed-capacity ring buffer for one odometry source.
//
// Covariance convention (this slice, all three native forms): the covariance
// supplied to a push helper is treated as the *per-step increment* covariance
// (the covariance of the motion from the previous entry to this one), expressed in
// the source's body-increment tangent frame, in [trans; rot] order (matching
// Delta::cov). In particular, push_absolute does NOT difference two absolute-pose
// covariances (the difference of two absolute Sigmas is not guaranteed PSD); it
// takes the supplied Sigma directly as the increment covariance. query() then
// accumulates these increments additively across the window. This is a simple,
// consistent first-order model; a fuller adjoint-transported model is deferred (it
// is not needed for the median weighting this feeds).
class SourceBuffer {
public:
    SourceBuffer() = default;

    // Preallocate to 'capacity' entries and bind the per-sensor knobs used by the
    // covariance combine. No heap allocation occurs after this call. Re-calling
    // configure() re-binds and clears the buffer (reallocates only if capacity grew).
    // Returns OutOfRange for capacity < 2, InvalidConfig for an unknown combine.
    Status configure(const SensorConfig& sensor, OdomForm form, int capacity,
                     ConfCombine combine, Scalar confidence_blend = Scalar(0.5));

    // Drop all entries (keeps configuration + storage).
    void reset();

    // Push one sample in the configured native form. Pushing a different form than
    // configured returns InvalidConfig. Stamps must be strictly increasing; a stamp
    // <= the newest returns OutOfRange. Pushing beyond capacity drops the oldest.
    Status push_twist(Timestamp t, const Vec6& xi, const Mat6& cov);
    Status push_increment(Timestamp t, const SE3& incr, const Mat6& cov);
    Status push_absolute(Timestamp t, const SE3& pose, const Mat6& cov);

    // True iff [t0, t1] is fully bracketed by stored entries (and t1 > t0).
    bool covers(Timestamp t0, Timestamp t1) const;

    // delta(t0,t1): motion = inverse(T(t0)) o T(t1) with cum_pose interpolated at
    // the endpoints; cov = window-accumulated native incr_cov (fractional end
    // scaling) combined with a modeled cov synthesized from the motion magnitude.
    //   OutOfRange if t1 <= t0; NoData if the window is not covered.
    Expected<Delta> query(Timestamp t0, Timestamp t1) const;

    Timestamp oldest() const;   // stamp of the oldest entry (0 if empty)
    Timestamp newest() const;   // stamp of the newest entry (0 if empty)
    int       size() const { return count_; }
    int       capacity() const { return capacity_; }
    bool      configured() const { return configured_; }

private:
    const BufferEntry& at(int logical) const;  // logical 0 = oldest
    int  bracket_upper(Timestamp t) const;      // smallest logical idx with stamp >= t
    SE3  pose_at(Timestamp t) const;            // interpolated cum_pose at t (covered)
    Mat6 accumulate_native_cov(Timestamp t0, Timestamp t1) const;
    Mat6 modeled_cov(const SE3& motion) const;
    // Per D7: combine the native and modeled covariances per the configured rule.
    // When no usable native covariance is available (have_native == false) the
    // native operand is the *identity* (D7: "missing native -> assume identity"),
    // NOT modeled-only — so e.g. Sum still adds the modeled term to Identity.
    Mat6 combine_cov(const Mat6& native, const Mat6& modeled, bool have_native) const;

    std::vector<BufferEntry> buf_;          // sized in configure(); never grows after
    int          capacity_   = 0;
    int          head_       = 0;           // physical index of the oldest entry
    int          count_      = 0;
    bool         configured_ = false;
    bool         any_native_ = false;       // a non-zero native cov was pushed at least once
    OdomForm     form_       = OdomForm::Increment;
    SensorConfig sensor_{};
    ConfCombine  combine_    = ConfCombine::Sum;
    Scalar       blend_      = Scalar(0.5);
};

} // namespace ofc
#endif // OFC_CORE_BUFFER_HPP
