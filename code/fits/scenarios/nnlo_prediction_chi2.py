"""chi^2 of the parameter-free NNLO prediction (v2: physical masses, P = Pref) for the data sets below threshold,
normalisations fixed (n=1), point errors stat (+) p2p only, on the MAIN database (pre-2015 + A2 2022 minus MAIN_DROP,
the 3-sigma outlier; the HIgS 2022 sets are NOT part of the main database since 2026-10-17 and are listed separately
as a test); the sets containing a dropped point are listed again with all points for comparison.
Writes scenarios/nnlo_prediction_chi2.txt"""
import os, sys, numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(here, '..', 'lib')); os.chdir(os.path.join(here, '..'))
from common import *
df = load_db(); T = outputs(df, 'NNLO')
s_all = select(df, list(df.dataset.unique()), Emax=150.0)
s_all = s_all[s_all.theta_lab > 0.5]
s_all = s_all[~s_all.dataset.isin(HIGS)]                   # main-database sets (HIgS listed separately below)
s = s_all[~s_all.pid.isin(MAIN_DROP)]
s_higs = select(df, HIGS, Emax=150.0)
pred_all = T.predict_array(PREF_NNLO)
lines = ['Pref = ' + ' '.join('%.3f' % x for x in PREF_NNLO),
         'MAIN database: pre-2015 + A2 2022 minus MAIN_DROP = %s (HIgS 2022 excluded since 2026-10-17)' % (list(MAIN_DROP),)]
tot = 0; ntot = 0


def chi2_of(g):
    idx = [T.index[p] for p in g.pid]
    r = (g.value.values - pred_all[idx]) / g.err_pt.values
    return float(np.sum(r ** 2))


for ds, g in s.groupby('dataset'):
    c2 = chi2_of(g); tot += c2; ntot += len(g)
    lines.append(f'{ds:26s} N={len(g):3d} chi2={c2:7.1f}')
lines.append(f'TOTAL (main database) N={ntot} chi2={tot:.1f}')
dropped_sets = sorted(set(s_all[s_all.pid.isin(MAIN_DROP)].dataset))
for ds in dropped_sets:
    g = s_all[s_all.dataset == ds]
    lines.append(f'{ds:26s} N={len(g):3d} chi2={chi2_of(g):7.1f}   [all points incl. the 3-sigma outlier]')
th = 0
for ds, g in s_higs.groupby('dataset'):
    c2 = chi2_of(g); th += c2
    lines.append(f'{ds:26s} N={len(g):3d} chi2={c2:7.1f}   [test: HIgS 2022, not in the main database]')
lines.append(f'TOTAL (main database + HIgS) N={ntot + len(s_higs)} chi2={tot + th:.1f}')
open('scenarios/nnlo_prediction_chi2.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
