// ofc_adapters/ros_node.hpp — ROS node SKETCH (Slice 13 deliverable; the REAL node is 13b).
//
// ====================================================================================
//  DEFERRED TO SLICE 13b. This header is a DOCUMENTED, COMPILE-GUARDED SKETCH only.
//  It is NOT compiled into ofc_adapters (no ROS toolchain on the build box; the project
//  keeps deps to Eigen + doctest), it links nothing, and it registers no test. The whole
//  body is behind `#if defined(OFC_HAS_ROS)`, which is NEVER defined by this build — so
//  including the header is a no-op everywhere.
//
//  The ISSUES.md Slice-13 "Done when: ROS node round-trips on a recorded bag" criterion is
//  the 13b deliverable: a real <ros/ros.h> (or rclcpp) build, a node wiring the mappings
//  sketched below, and a recorded-bag round-trip test on a machine with ROS installed.
// ====================================================================================
//
// INTENDED MAPPING (the design this sketch records, for 13b to implement):
//
//   INPUTS  (ROS -> core)
//   ---------------------------------------------------------------------------
//   * Each odometry topic (nav_msgs/Odometry) or TF chain becomes an ISource adapter. The
//     subscriber callback pushes the incoming pose/twist increments into a per-source
//     SourceBuffer (the same ring buffer sim/ uses); ISource::query(t0,t1) answers the core's
//     uniform "relative SE(3) motion over [t0,t1] + covariance" contract from that buffer
//     (converting the native form — absolute pose differenced, or twist integrated — into a
//     body delta, attaching the message's covariance as Delta.cov).
//   * An absolute-fix topic (GPS NavSatFix / a map-match PoseWithCovariance) becomes an
//     ICorrection adapter feeding add_correction() (Slice 11 path).
//   * The clock comes from the message stamps (or /clock under sim time), fed as the
//     Timestamp into step() — exactly the timestamp queue the ThreadedEstimator wrapper takes,
//     so the ROS node = a ThreadedEstimator driven by the subscriber-stamp stream.
//
//   OUTPUTS (core -> ROS), once per published step from Result:
//   ---------------------------------------------------------------------------
//   * Result.frontier (the trustworthy causal state) -> nav_msgs/Odometry on ~odometry, with
//     Result.frontier.cov mapped into the 6x6 pose+twist covariance; the base->odom TF.
//   * Result.tip (predict-only extrapolation to now) -> an optional low-latency Odometry on
//     ~odometry_tip when Result.tip_valid (tip_cov_inflation already applied by the core).
//   * Result.calib[i] (per-source CalibSnapshot: extrinsic, scale, time_offset + the per-DOF
//     confidences/commit flags) -> a custom calibration msg (or static TF per committed
//     extrinsic) on ~calibration, so downstream nodes see the online-calibrated mounts.
//   * Result.phase / readiness / Result.correction (CorrectionDiag) -> a diagnostics msg.
//
//   THREADING: the core stays caller-pumped (D14). The node owns a ThreadedEstimator
//   (ofc_adapters/threaded_estimator.hpp): ROS subscriber callbacks submit() timestamps; the
//   wrapper's single worker thread pumps step(); a ROS timer publishes snapshot() at the
//   output rate. The ROS spinner threads never touch the core directly.
//
//   CONFIG: ros::NodeHandle params -> the same key=value text the ConfigLoader
//   (ofc_adapters/config_loader.hpp) parses (or build the Config + the owned
//   std::vector<SensorConfig> directly from the param tree), then Estimator::init(cfg).
//
//   PERSISTENCE: a ~save_path param wires FilePersistence (ofc_adapters/file_persistence.hpp)
//   — load() on startup for a warm restart, save() on a timer / on shutdown.
#ifndef OFC_ADAPTERS_ROS_NODE_HPP
#define OFC_ADAPTERS_ROS_NODE_HPP

// The entire node is compiled ONLY when a ROS toolchain defines OFC_HAS_ROS. This build never
// defines it, so nothing below is compiled here — the real node + its recorded-bag round-trip
// is Slice 13b.
#if defined(OFC_HAS_ROS)

#include "ofc/core/config.hpp"
#include "ofc/core/estimator.hpp"
#include "ofc/core/result.hpp"
#include "ofc/core/source.hpp"

#include "ofc_adapters/config_loader.hpp"
#include "ofc_adapters/file_persistence.hpp"
#include "ofc_adapters/threaded_estimator.hpp"

// ---- 13b will include the real ROS headers here, e.g.: ------------------------------------
//   #include <ros/ros.h>
//   #include <nav_msgs/Odometry.h>
//   #include <tf2_ros/transform_broadcaster.h>
//   #include <sensor_msgs/NavSatFix.h>
// (ROS 2: rclcpp / rclcpp::Node, nav_msgs/msg/odometry.hpp, tf2_ros, etc.)

#include <memory>
#include <vector>

namespace ofc {
namespace adapters {

// SKETCH ONLY — signatures showing the intended shape; bodies are 13b. None of this is
// instantiated by the current build (the whole TU is behind OFC_HAS_ROS).
class RosOdometrySource;  // : public ISource — wraps a nav_msgs/Odometry subscriber + buffer
class RosAbsoluteFix;     // : public ICorrection — wraps a GPS / map-match subscriber

class OfcRosNode {
public:
    // 13b: read params -> Config (via ConfigLoader / the param tree), init the Estimator,
    // create one RosOdometrySource per configured odom topic + any RosAbsoluteFix, optionally
    // FilePersistence::load() for a warm restart, then wrap the Estimator in a
    // ThreadedEstimator.
    explicit OfcRosNode(/* ros::NodeHandle& nh, ros::NodeHandle& pnh */);

    // 13b: subscriber callbacks push into the per-source buffers + submit() the stamp to the
    // ThreadedEstimator; an output timer publishes snapshot() as Odometry/TF/calib/diag msgs;
    // a save timer calls FilePersistence::save().
    void spin();

private:
    // Estimator             est_;
    // ConfigLoader          loader_;
    // std::vector<std::unique_ptr<RosOdometrySource>> sources_;
    // std::unique_ptr<FilePersistence>                persist_;
    // std::unique_ptr<ThreadedEstimator>              threaded_;
};

} // namespace adapters
} // namespace ofc

#endif // OFC_HAS_ROS
#endif // OFC_ADAPTERS_ROS_NODE_HPP
