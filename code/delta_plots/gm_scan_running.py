"""gm_scan_running.py -- per-set chi^2(g_M) minima of the parameter-free prediction with the running Delta parameters (irun):
local parabola through g_M = 2.87, 2.97, 3.07 (the cached runs of theory_dp.Theory, each the complete prediction at its own
Pref, as in the figures' g_M band), n = 1, 150 < E <= 350 MeV; curvature error only.  Comparison values for the default
(notes/delta_plots.md, chi^2(g_M) scan of s05_delta.py): MAMI 75 2.91(2), MAMI 90 2.98(2), MAMI 131 3.00(5), MAMI 136 2.98(3),
A2 S2z 2.89(13), LEGS dsigma 3.02, LEGS S3 beyond 3.5.   usage: DP_IRUN=2 python gm_scan_running.py"""
import os, sys
IRUN = int(os.environ.get('DP_IRUN', '2'))
os.environ['DP_V2OPTS'] = 'pimassphys=t, nmassphys=t, irun=%d' % IRUN
os.environ['DP_MAIN'] = 'pred'
os.environ['DP_EXTRA'] = 'novc'       # no extra running run (2026-10-17 pipeline switch): the whole Theory is the running amplitude here
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import config as C, data_dp as D, theory_dp as TD, make_figs as MF

df = D.build(); pts = MF.build_points(df)
seen = set(); P = []
for p in pts:
    if p[0] not in seen:
        seen.add(p[0]); P.append(p)
ix = {p[0]: i for i, p in enumerate(P)}
T = TD.Theory(P, tag='dp_irun%d' % IRUN, verbose=False)
gm = T.gM_runs()
curves = {C.GM_LO: gm[C.GM_LO], C.GM0: T.main(), C.GM_HI: gm[C.GM_HI]}
s = df[(df.E > C.E_DELTA_MIN) & (df.E <= C.E_DELTA_MAX) & df.obs.isin(['CS', 'CSCM', 'S3', 'S2X', 'S2Z'])].copy()
s['set'] = s.dataset.replace({'Hallin1993_TabII': 'Hallin1993_TabII(D)', 'Hallin1993_TabII_D': 'Hallin1993_TabII(D)'})
rows = []
for ds in [k for k, _ in MF.PAPERROWS]:
    q = s[s.set == ds]
    ii = np.array([ix[p] for p in q.pid]); y = q.value.values; e = q.err.values
    g = np.array(sorted(curves)); c2 = np.array([np.sum(((y - curves[gv][ii]) / e) ** 2) for gv in g])
    a, b, c = np.polyfit(g, c2, 2)
    gmin = -b / (2 * a) if a > 0 else np.nan
    err = 1 / np.sqrt(a) if a > 0 else np.nan
    rows.append(dict(set=ds, N=len(q), chi2_287=c2[0], chi2_297=c2[1], chi2_307=c2[2], gM_min=gmin, err=err, chi2_min=np.polyval([a, b, c], gmin) if a > 0 else np.nan))
t = pd.DataFrame(rows)
out = os.path.join(C.RESULTS, 'running', 'gM_scan_irun%d.csv' % IRUN); os.makedirs(os.path.dirname(out), exist_ok=True)
t.to_csv(out, index=False)
print('irun=%d: per-set chi^2(gM), n = 1, parabola through 2.87/2.97/3.07 (minimum outside [2.6, 3.4] = not constrained)' % IRUN)
for r in t.itertuples():
    print('  %-22s N=%2d chi2(2.87/2.97/3.07) = %7.1f %7.1f %7.1f  gM_min = %s' % (
        r.set, r.N, r.chi2_287, r.chi2_297, r.chi2_307, ('%.2f(%.2f) [chi2 %.1f]' % (r.gM_min, r.err, r.chi2_min)) if 2.6 < r.gM_min < 3.4 else ('%.2f (outside)' % r.gM_min if np.isfinite(r.gM_min) else 'no minimum')))
