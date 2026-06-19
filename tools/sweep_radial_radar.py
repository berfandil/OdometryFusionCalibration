"""Sweep radar scan-match params (tau_d, num_points, tau_inlier) on one RADIal run; report the
per-step yaw correlation vs lidar GT for each. Gate off (margins huge) so we measure raw geometry.
"""
import os, sys, io, contextlib
import radial_to_csv as r2c

run = sys.argv[1]

def trial(**kw):
    a = r2c.parse_args([run, '--hz', '4.1', '--tau-v', '1000', '--yaw-margin-deg', '180',
                        '--trans-margin-m', '100'])
    for k, v in kw.items():
        setattr(a, k, v)
    with contextlib.redirect_stdout(io.StringIO()):
        rep, stats, perframe = r2c.run(run, run, a)
    return rep

print('%-40s | %5s %6s %6s | %8s %8s %8s' %
      ('config', 'corr', 'inlier', 'rotmed', 'cumEst', 'cumGT', 'rotp90'))
configs = [
    dict(num_points=300, tau_d=1.0, tau_inlier=1.0, baseline=1),
    dict(num_points=300, tau_d=1.5, tau_inlier=1.5, baseline=1),
    dict(num_points=300, tau_d=2.0, tau_inlier=2.0, baseline=1),
    dict(num_points=300, tau_d=1.0, tau_inlier=1.0, baseline=2),
    dict(num_points=300, tau_d=1.0, tau_inlier=1.0, baseline=4),
    dict(num_points=300, tau_d=1.5, tau_inlier=1.5, baseline=8),
]
for c in configs:
    r = trial(**c)
    label = 'np=%d td=%.1f ti=%.1f N=%d' % (c['num_points'], c['tau_d'], c['tau_inlier'], c['baseline'])
    print('%-40s | %5.2f %6.1f %6.2f | %8.1f %8.1f %8.2f' %
          (label, r['yaw_corr'], r['avg_inl'], r['rot_med'],
           r['cum_est_yaw_deg'], r['cum_gt_yaw_deg'], r['rot_p90']))
