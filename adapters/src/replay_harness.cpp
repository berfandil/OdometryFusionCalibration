// ofc_adapters/replay_harness.cpp — see replay_harness.hpp for the design.
//
// RELAXED EDGE. The real-data analogue of sim/src/rig.cpp: pump Estimator::step() over a merged
// timeline, record each published step, compute drift/NEES/NIS vs an optional GT track.
#include "ofc_adapters/replay_harness.hpp"

#include "ofc/core/lie.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace ofc {
namespace adapters {

namespace {
constexpr Scalar kNanosPerSec = Scalar(1e9);
Timestamp secs_to_ns(Scalar s) { return static_cast<Timestamp>(std::llround(s * kNanosPerSec)); }

// 6-DOF pose NEES (the EXACT convention of tests/test_validation.cpp pose_nees):
//   e = se3::log(T_est^-1 o T_gt)  ([trans;rot]),  nees = e^T (P_pp)^-1 e, P_pp = cov(0..5).
Scalar pose_nees(const SE3& est, const SE3& gt, const Mat12& cov12) {
    const SE3 err_T   = se3::compose(se3::inverse(est), gt);
    const Vec6 e      = se3::log(err_T);
    const Mat6 Ppp    = cov12.block<6, 6>(0, 0);
    const Vec6 Pinv_e = Ppp.ldlt().solve(e);
    return e.dot(Pinv_e);
}
} // namespace

Status ReplayHarness::run(const ReplayInputs& in) {
    records_.clear();
    summary_ = RunSummary{};
    error_.clear();

    if (in.cfg == nullptr)       { error_ = "null config"; return Status::InvalidConfig; }
    if (in.sources.empty())      { error_ = "no sources";  return Status::InvalidConfig; }
    for (CsvSource* s : in.sources)
        if (s == nullptr || !s->loaded()) { error_ = "a source is null/unloaded"; return Status::InvalidConfig; }

    const Status is = est_.init(*in.cfg);
    if (!ok(is)) { error_ = "estimator init failed"; return is; }
    for (CsvSource* s : in.sources) {
        const Status as = est_.add_source(s);
        if (!ok(as)) { error_ = "add_source failed"; return as; }
    }
    if (in.gps != nullptr) {
        const Status cs = est_.add_correction(in.gps);
        if (!ok(cs)) { error_ = "add_correction failed"; return cs; }
    }

    // --- build the merged tick timeline ---------------------------------------------------------
    // We step on a REGULAR grid at the configured tick_rate over the union of the source data
    // spans, then fold in the GPS-fix frontier ticks. A regular grid keeps the estimator's
    // fixed-cadence assumptions intact while the per-source SourceBuffers interpolate the (possibly
    // irregular) real samples to the queried windows — so an irregular real dataset still replays.
    std::vector<Timestamp> ticks;
    Timestamp span_lo = 0, span_hi = 0;
    bool have_span = false;
    for (CsvSource* s : in.sources) {
        if (s->row_count() < 2) continue;
        const Timestamp o = s->oldest(), n = s->newest();
        if (!have_span) { span_lo = o; span_hi = n; have_span = true; }
        else { span_lo = std::min(span_lo, o); span_hi = std::max(span_hi, n); }
    }
    if (!have_span) { error_ = "no source data span"; return Status::NoData; }

    // Tick grid at the configured rate over the data span. The estimator's frontier is
    // now - fusion_delay, so ticks before span_lo + fusion_delay simply do not fuse; that is fine.
    const Scalar rate = (in.cfg->tick_rate_hz > Scalar(0)) ? in.cfg->tick_rate_hz : Scalar(50.0);
    const Timestamp tick = secs_to_ns(Scalar(1) / rate);
    if (tick <= 0) { error_ = "bad tick rate"; return Status::InvalidConfig; }
    for (Timestamp now = span_lo; now <= span_hi; now += tick) ticks.push_back(now);
    // Fold in the GPS fix stamps (so a fix is never missed by grid aliasing) + their frontier
    // ticks (a fix is consumed at the step whose frontier reaches it -> tick = fix_stamp + delay).
    const Timestamp delay_ns = secs_to_ns(in.cfg->fusion_delay_s);
    for (const GpsFix& fx : in.gps_fixes) {
        ticks.push_back(fx.stamp + delay_ns);
    }
    std::sort(ticks.begin(), ticks.end());
    ticks.erase(std::unique(ticks.begin(), ticks.end()), ticks.end());
    if (ticks.empty()) { error_ = "empty timeline"; return Status::NoData; }

    // GPS fixes sorted by stamp (caller is asked to pre-sort; sort defensively).
    std::vector<GpsFix> fixes = in.gps_fixes;
    std::sort(fixes.begin(), fixes.end(),
              [](const GpsFix& a, const GpsFix& b) { return a.stamp < b.stamp; });
    std::size_t next_fix = 0;

    // --- GT anchoring (DESIGN §7; same as the sim Rig) -----------------------------------------
    const Timestamp window_ns = secs_to_ns(in.cfg->window_s);
    bool have_anchor = false;
    SE3  gt_anchor_inv;
    const bool has_gt = (in.gt != nullptr && in.gt->loaded());

    // --- pump ----------------------------------------------------------------------------------
    int fuses = 0;
    for (Timestamp now : ticks) {
        // Submit any GPS fix whose stamp the upcoming frontier (now - delay) has reached. The fix
        // is consumed INSIDE this step() at the frontier, matching the e2e GPS test's cadence.
        if (in.gps != nullptr) {
            const Timestamp frontier = now - delay_ns;
            while (next_fix < fixes.size() && fixes[next_fix].stamp <= frontier) {
                in.gps->submit_fix(fixes[next_fix]);
                ++next_fix;
            }
        }

        const Status st = est_.step(now);
        ReplayRecord rec;
        rec.now   = now;
        rec.fused = ok(st);
        if (!ok(st)) { records_.push_back(rec); continue; }

        rec.result = est_.latest();
        ++fuses;

        if (has_gt) {
            const Timestamp frontier = rec.result.frontier.stamp;
            const SE3 gt_abs = in.gt->pose_at(frontier);
            if (!have_anchor) {
                gt_anchor_inv = se3::inverse(in.gt->pose_at(frontier - window_ns));
                have_anchor   = true;
            }
            rec.has_gt  = true;
            rec.gt_pose = se3::compose(gt_anchor_inv, gt_abs);
            const SE3& fused = rec.result.frontier.pose;
            rec.trans_err_m = (fused.t - rec.gt_pose.t).norm();
            rec.rot_err_rad = so3::log(rec.gt_pose.R.transpose() * fused.R).norm();
            rec.pose_nees   = pose_nees(fused, rec.gt_pose, rec.result.frontier.cov);
        }

        records_.push_back(rec);
    }

    // --- aggregate -----------------------------------------------------------------------------
    summary_.steps  = static_cast<int>(records_.size());
    summary_.has_gt = has_gt;
    Timestamp last_fused_now = 0;
    for (const ReplayRecord& r : records_) if (r.fused) last_fused_now = r.now;
    const Timestamp tail_start = last_fused_now - secs_to_ns(in.tail_window_s);

    Scalar sum_te2 = 0, sum_re2 = 0;
    long   drift_n = 0;
    Scalar tail_te = 0, tail_re = 0; long tail_n = 0;
    Scalar sum_nees = 0; long nees_n = 0;
    Scalar sum_nis  = 0; long nis_n  = 0;
    int    fused_seen = 0;
    // Post-warmup fused-with-GT poses, in order, for the local relative-pose-error windows.
    std::vector<SE3> seg_est, seg_gt;

    for (const ReplayRecord& r : records_) {
        if (!r.fused) continue;
        ++summary_.fused_steps;
        ++fused_seen;
        const CorrectionDiag& cd = r.result.correction;
        summary_.gps_evaluated += cd.corr_evaluated;
        summary_.gps_applied   += cd.corr_applied;
        summary_.gps_rejected  += cd.corr_rejected;
        if (cd.corr_applied > 0) { sum_nis += cd.last_nis; ++nis_n; }

        if (r.has_gt) {
            sum_te2 += r.trans_err_m * r.trans_err_m;
            sum_re2 += r.rot_err_rad * r.rot_err_rad;
            ++drift_n;
            if (r.trans_err_m > summary_.max_trans_m) summary_.max_trans_m = r.trans_err_m;
            if (r.rot_err_rad > summary_.max_rot_rad)  summary_.max_rot_rad = r.rot_err_rad;
            if (r.now >= tail_start) { tail_te += r.trans_err_m; tail_re += r.rot_err_rad; ++tail_n; }
            // NEES + local windows after the warmup transient (mirrors the sim validation harness).
            if (fused_seen > in.warmup_steps) {
                sum_nees += r.pose_nees; ++nees_n;
                seg_est.push_back(r.result.frontier.pose);
                seg_gt.push_back(r.gt_pose);
            }
        }
    }

    if (drift_n > 0) {
        summary_.rms_trans_m = std::sqrt(sum_te2 / static_cast<Scalar>(drift_n));
        summary_.rms_rot_rad = std::sqrt(sum_re2 / static_cast<Scalar>(drift_n));
    }
    if (tail_n > 0) {
        summary_.tail_trans_m = tail_te / static_cast<Scalar>(tail_n);
        summary_.tail_rot_rad = tail_re / static_cast<Scalar>(tail_n);
    }
    summary_.nees_count     = nees_n;
    summary_.mean_pose_nees = (nees_n > 0) ? sum_nees / static_cast<Scalar>(nees_n) : Scalar(0);
    summary_.nis_count      = nis_n;
    summary_.mean_nis       = (nis_n > 0) ? sum_nis / static_cast<Scalar>(nis_n) : Scalar(0);

    // --- local (GT-anchored fixed-window) relative-pose error ----------------------------------
    // Tile the post-warmup fused-with-GT poses into non-overlapping windows of local_batch_len.
    // Per window: rel_est = est0^-1 o estL, rel_gt = gt0^-1 o gtL, error E = rel_gt^-1 o rel_est
    // (est's windowed motion expressed in the GT-start frame -> the window "starts from GT").
    // Length-independent: each window measures only intra-window drift, so the aggregate is
    // comparable across recordings of different total length. (KITTI-style segment metric.)
    const int L = in.local_batch_len;
    const int seg_n = static_cast<int>(seg_est.size());
    if (L > 0 && seg_n >= L) {
        const int nb = seg_n / L;
        std::vector<Scalar> wt; wt.reserve(nb);
        std::vector<Scalar> wr; wr.reserve(nb);
        Scalar sum_t = 0, sum_r = 0, max_t = 0, max_r = 0;
        for (int b = 0; b < nb; ++b) {
            const int i0 = b * L, i1 = b * L + L - 1;
            const SE3 rel_est = se3::compose(se3::inverse(seg_est[i0]), seg_est[i1]);
            const SE3 rel_gt  = se3::compose(se3::inverse(seg_gt[i0]),  seg_gt[i1]);
            const SE3 E       = se3::compose(se3::inverse(rel_gt), rel_est);
            const Scalar te = E.t.norm();
            const Scalar re = so3::log(E.R).norm();
            wt.push_back(te); wr.push_back(re);
            sum_t += te; sum_r += re;
            if (te > max_t) max_t = te;
            if (re > max_r) max_r = re;
        }
        std::sort(wt.begin(), wt.end());
        std::sort(wr.begin(), wr.end());
        summary_.local_batch_len    = L;
        summary_.local_batch_count  = nb;
        summary_.local_dropped      = seg_n - nb * L;
        summary_.local_mean_trans_m = sum_t / static_cast<Scalar>(nb);
        summary_.local_mean_rot_rad = sum_r / static_cast<Scalar>(nb);
        summary_.local_med_trans_m  = wt[nb / 2];
        summary_.local_med_rot_rad  = wr[nb / 2];
        summary_.local_max_trans_m  = max_t;
        summary_.local_max_rot_rad  = max_r;
    }

    if (fuses == 0) { error_ = "no steps fused"; return Status::NoData; }
    return Status::Ok;
}

Status ReplayHarness::write_results_csv(std::ostream& os) const {
    if (records_.empty()) { return Status::InvalidConfig; }

    os << std::setprecision(12);
    os << "# ofc replay results (one row per fused step)\n";
    os << "now_ns,frontier_ns,phase,readiness,"
          "x,y,z,qw,qx,qy,qz,"
          "p00,p11,p22,p33,p44,p55,"
          "tip_valid,tip_x,tip_y,tip_z,"
          "corr_evaluated,corr_applied,corr_rejected,last_nis,"
          "has_gt,gt_x,gt_y,gt_z,trans_err_m,rot_err_rad,pose_nees,"
          "source_count";
    // A fixed-width per-source block header (up to 8 sources covers the default max_sources).
    os << ",[per-source: id,scale,off_s,conf,extr_conf...]\n";

    for (const ReplayRecord& r : records_) {
        if (!r.fused) continue;
        const State& fr = r.result.frontier;
        Scalar qw, qx, qy, qz;
        mat3_to_quat(fr.pose.R, qw, qx, qy, qz);
        os << r.now << ',' << fr.stamp << ',' << static_cast<int>(r.result.phase) << ','
           << r.result.readiness << ','
           << fr.pose.t.x() << ',' << fr.pose.t.y() << ',' << fr.pose.t.z() << ','
           << qw << ',' << qx << ',' << qy << ',' << qz;
        for (int i = 0; i < 6; ++i) os << ',' << fr.cov(i, i);
        os << ',' << (r.result.tip_valid ? 1 : 0) << ','
           << r.result.tip.pose.t.x() << ',' << r.result.tip.pose.t.y() << ','
           << r.result.tip.pose.t.z();
        const CorrectionDiag& cd = r.result.correction;
        os << ',' << cd.corr_evaluated << ',' << cd.corr_applied << ',' << cd.corr_rejected
           << ',' << cd.last_nis;
        os << ',' << (r.has_gt ? 1 : 0) << ','
           << r.gt_pose.t.x() << ',' << r.gt_pose.t.y() << ',' << r.gt_pose.t.z() << ','
           << r.trans_err_m << ',' << r.rot_err_rad << ',' << r.pose_nees;
        os << ',' << r.result.source_count;
        for (int s = 0; s < r.result.source_count; ++s) {
            const CalibSnapshot& c = r.result.calib[s];
            os << ',' << static_cast<int>(c.id) << ',' << c.scale << ',' << c.time_offset_s
               << ',' << c.confidence << ',' << c.extrinsic_confidence;
        }
        os << '\n';
    }

    // Trailing summary block.
    os << "# summary:\n";
    write_summary(os);
    return Status::Ok;
}

Status ReplayHarness::write_results_csv_file(const std::string& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) return Status::NoData;
    return write_results_csv(f);
}

void ReplayHarness::write_summary(std::ostream& os) const {
    const RunSummary& s = summary_;
    os << "# steps=" << s.steps << " fused=" << s.fused_steps
       << " gps_eval=" << s.gps_evaluated << " gps_applied=" << s.gps_applied
       << " gps_rejected=" << s.gps_rejected << '\n';
    if (s.has_gt) {
        os << "# drift: max_trans_m=" << s.max_trans_m << " max_rot_rad=" << s.max_rot_rad
           << " rms_trans_m=" << s.rms_trans_m << " rms_rot_rad=" << s.rms_rot_rad
           << " tail_trans_m=" << s.tail_trans_m << " tail_rot_rad=" << s.tail_rot_rad << '\n';
        os << "# consistency: mean_pose_nees=" << s.mean_pose_nees << " (N=" << s.nees_count
           << ", target~DOF=6)";
        if (s.nis_count > 0)
            os << " mean_nis=" << s.mean_nis << " (N=" << s.nis_count << ", target~DOF=3)";
        os << '\n';
        if (s.local_batch_count > 0) {
            os << "# local (GT-anchored " << s.local_batch_len << "-step windows, N="
               << s.local_batch_count;
            if (s.local_dropped > 0) os << " dropped_tail=" << s.local_dropped;
            os << "): trans_m mean=" << s.local_mean_trans_m << " p50=" << s.local_med_trans_m
               << " max=" << s.local_max_trans_m
               << " ; rot_rad mean=" << s.local_mean_rot_rad << " p50=" << s.local_med_rot_rad
               << " max=" << s.local_max_rot_rad << '\n';
        }
    } else {
        os << "# (no GT track: drift/NEES not computed)\n";
        if (s.nis_count > 0)
            os << "# mean_nis=" << s.mean_nis << " (N=" << s.nis_count << ", target~DOF=3)\n";
    }

    // Final per-source calibration snapshot (the last fused step) — what online calibration
    // converged to. extrinsic rotation is reported as so3::log(R) = [rx,ry,rz] (rad; rz is yaw),
    // translation as the lever [x,y,z] (m). Confidences/commit flags expose the observability gate.
    const ReplayRecord* last = nullptr;
    for (auto it = records_.rbegin(); it != records_.rend(); ++it)
        if (it->fused) { last = &(*it); break; }
    if (last != nullptr && last->result.source_count > 0) {
        os << "# final calib (per source, at last fused step):\n";
        for (int i = 0; i < last->result.source_count; ++i) {
            const CalibSnapshot& c = last->result.calib[i];
            const Vec3 rlog = so3::log(c.extrinsic.R);
            os << "#   src" << static_cast<int>(c.id)
               << " scale=" << c.scale << " (conf " << c.scale_confidence
               << (c.scale_committed ? ",committed" : "") << ")"
               << " time_offset_s=" << c.time_offset_s << " (conf " << c.confidence
               << (c.committed ? ",committed" : "") << ")"
               << " extr_rot[rx,ry,rz]=[" << rlog.x() << "," << rlog.y() << "," << rlog.z()
               << "] (conf " << c.extrinsic_confidence
               << (c.extrinsic_committed ? ",committed" : "") << ")"
               << " lever[x,y,z]=[" << c.extrinsic.t.x() << "," << c.extrinsic.t.y() << ","
               << c.extrinsic.t.z() << "] (conf " << c.translation_confidence
               << (c.translation_committed ? ",committed" : "") << ")\n";
        }
    }
}

} // namespace adapters
} // namespace ofc
