"""chi2_running.py -- chi^2 per point of the Delta-region data (150 < E <= 350 MeV) for the parameter-free NNLO prediction with
the running Delta parameters of code_v2 (&v2 irun = 1, 2, 3; work/notes/delta_running2.md), side by side with the paper numbers
(results/chi2_table_E350.csv, main curve = prediction, irun = 0), plus the data/curve excess on the rising flank at backward
angles (180-320 MeV; the statistics of notes/delta_plots.md).  Nothing of the paper pipeline (figures, results/*.csv) is
overwritten: outputs go to results/running/.  Two evaluations per scheme: at the own static polarisabilities of the running
amplitude (driver Pref: irun=2 shifts gamma_E1E1 by -0.05 and gamma_M1M1 by -0.37; irun=3 none) and at the static
polarisabilities of the default amplitude (contact terms restore the NNLO values).  The band ('incl. th.') is the LMP15-type
envelope computed with the running amplitude in the same way as for the paper.
usage: [DP_IRUN=2] python chi2_running.py   (~5 min per scheme from an empty cache)"""
import os, sys, json
IRUN = int(os.environ.get('DP_IRUN', '2'))
os.environ['DP_V2OPTS'] = 'pimassphys=t, nmassphys=t, irun=%d' % IRUN
os.environ['DP_MAIN'] = 'pred'
os.environ['DP_EXTRA'] = 'novc'       # no extra running run (2026-10-17 pipeline switch): the whole Theory is the running amplitude here
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import pandas as pd
import config as C
import data_dp as D
import theory_dp as TD
import make_figs as MF

OUTDIR = os.path.join(C.RESULTS, 'running')
os.makedirs(OUTDIR, exist_ok=True)
df = D.build()
pts = MF.build_points(df)
seen = set(); P = []
for p in pts:
    if p[0] not in seen:
        seen.add(p[0]); P.append(p)
ix = {p[0]: i for i, p in enumerate(P)}
print('irun=%d, &v2 options: %s, %d points' % (IRUN, C.v2opts(True), len(P)))
T = TD.Theory(P, tag='dp_irun%d' % IRUN)
val = T.validate()
for k, v in val.items():
    print('  validation %-28s %s' % (k, v))
info = json.load(open(os.path.join(C.RESULTS, 'run_info.json')))
P0 = np.array(info['Pref'], float)                 # static polarisabilities of the default prediction (paper)
print('  Pref(irun=%d) = %s' % (IRUN, np.round(T.Pref, 4)))
print('  Pref(default) = %s' % np.round(P0, 4))
cv_own = T.curves()
T.P_main = lambda: P0.copy()
cv_P0 = T.curves()

paper = pd.read_csv(os.path.join(C.RESULTS, 'chi2_table_E350.csv')).set_index('set')
rows = {}
for lab, cv in (('own', cv_own), ('P0', cv_P0)):
    t = MF.chi2_rows(df, cv, ix, C.E_DELTA_MIN, C.E_DELTA_MAX).set_index('set')
    rows[lab] = t
sets = [s for s in MF.SETORDER if s in paper.index]
tab = []
for s in sets:
    r = dict(set=s, label=MF.SETNAME.get(s, s), N=int(paper.loc[s, 'N']), floated=bool(paper.loc[s, 'norm'] > 0),
             paper_n1=paper.loc[s, 'chi2N_n1'], paper_nfit=paper.loc[s, 'chi2N_nfit'], paper_n=paper.loc[s, 'n'], paper_th=paper.loc[s, 'chi2N_th'])
    for lab in ('own', 'P0'):
        q = rows[lab].loc[s]
        r.update({lab + '_n1': q.chi2N_n1, lab + '_nfit': q.chi2N_nfit, lab + '_n': q.n, lab + '_th': q.chi2N_th})
    tab.append(r)
tab = pd.DataFrame(tab)
tab.to_csv(os.path.join(OUTDIR, 'chi2_table_irun%d.csv' % IRUN), index=False)

# ---- excess on the rising flank, 180-320 MeV, backward angles (as notes/delta_plots.md)
dwt = pd.read_csv(os.path.join(C.RESULTS, 'data_with_theory.csv'))            # paper: th_main = prediction
dd = df[df.obs.isin(['CS', 'CSCM'])].copy()
ii = np.array([ix[p] for p in dd.pid])
dd['th_own'] = cv_own['main'][ii]; dd['th_P0'] = cv_P0['main'][ii]
dd = dd.merge(dwt[['pid', 'th_main']], on='pid', how='left')
FLANK = [('Hunger1997_75', 'MAMI 97 75 cm', None), ('Hunger1997_90', 'MAMI 97 90 cm', None), ('Wissmann1999', 'MAMI 99 130.7 lab', None),
         ('Camen2002', 'MAMI 02 136.2 lab', None), ('LEGS2001_CS', 'LEGS dsigma, theta_cm >= 100', 100.0),
         ('Hallin1993_TabII', 'SAL 93, theta_cm >= 100', 100.0), ('Hallin1993_TabII_D', 'SAL 93 (D), theta_cm >= 100', 100.0)]
exc = []
for ds, lab, thmin in FLANK:
    m = (dd.dataset == ds) & (dd.E >= 180) & (dd.E <= 320)
    if thmin is not None:
        m &= dd.theta_cm >= thmin
    q = dd[m]
    if len(q) == 0:
        continue
    r = dict(set=lab, N=len(q))
    for k in ('th_main', 'th_own', 'th_P0'):
        x = q.value / q[k] - 1
        r[k.replace('th_', '') + '_mean'] = 100 * x.mean(); r[k.replace('th_', '') + '_median'] = 100 * x.median()
    exc.append(r)
exc = pd.DataFrame(exc)
exc.to_csv(os.path.join(OUTDIR, 'excess_irun%d.csv' % IRUN), index=False)

L = ['# Delta-region chi^2 per point: paper prediction (irun=0) vs running irun=%d' % IRUN, '',
     'Pref(irun=%d) = %s; Pref(default) = %s' % (IRUN, np.round(T.Pref, 4).tolist(), np.round(P0, 4).tolist()),
     "'own' = running amplitude at its own static polarisabilities; 'P0' = at the default static polarisabilities (contact terms).",
     'Columns: chi^2/N with n = 1 / n fitted / n / incl. the LMP15-type band.', '',
     '| data set | N | paper n=1 | own n=1 | P0 n=1 | paper nfit | own nfit | P0 nfit | paper n | own n | P0 n | paper th | own th | P0 th |',
     '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for r in tab.itertuples():
    f = lambda v: ('%.2f' % v) if r.floated else '--'
    g = lambda v: ('%.3f' % v) if r.floated else '1'
    L.append('| %s | %d | %.2f | %.2f | %.2f | %s | %s | %s | %s | %s | %s | %.2f | %.2f | %.2f |' % (
        r.label, r.N, r.paper_n1, r.own_n1, r.P0_n1, f(r.paper_nfit), f(r.own_nfit), f(r.P0_nfit), g(r.paper_n), g(r.own_n), g(r.P0_n),
        r.paper_th, r.own_th, r.P0_th))
L += ['', 'Totals (sum of chi^2 over the ten paper sets, n = 1 / n fitted / incl. band):']
pr = [s for s, _ in MF.PAPERROWS]
tt = tab.set_index('set').loc[pr]
for lab in ('paper', 'own', 'P0'):
    L.append('* %-5s %.1f / %.1f / %.1f' % (lab, (tt[lab + '_n1'] * tt.N).sum(), (tt[lab + '_nfit'] * tt.N).sum(), (tt[lab + '_th'] * tt.N).sum()))
L += ['', '## Excess of the data over the curve on the rising flank, 180-320 MeV, backward angles (data/curve - 1, %)', '',
      '| set | N | paper mean | paper median | own mean | own median | P0 mean | P0 median |', '|---|---|---|---|---|---|---|---|']
for r in exc.itertuples():
    L.append('| %s | %d | %+.1f | %+.1f | %+.1f | %+.1f | %+.1f | %+.1f |' % (r.set, r.N, r.main_mean, r.main_median, r.own_mean, r.own_median, r.P0_mean, r.P0_median))
open(os.path.join(OUTDIR, 'chi2_table_irun%d.md' % IRUN), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
