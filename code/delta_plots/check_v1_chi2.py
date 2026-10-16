"""
check_v1_chi2.py -- regression check of the chi^2 machinery: with the v1 executable (work/code) and the
main-fit values of fits_L.json (alpha = 11.6937, beta = 2.8793) the per-set chi^2/N (normalisation fitted)
must reproduce fit_report Sect. 5.1, column 'main-fit alpha, beta (spin NNLO)'.
Writes results/check_v1_chi2.md.   usage: python check_v1_chi2.py
"""
import os
os.environ['AMPL_EXE'] = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'code', 'execs', 'ampl.exe')
os.environ['DP_ALPHA'] = '11.6937'
os.environ['DP_BETA'] = '2.8793'
import numpy as np
import config as C
import data_dp as D
import theory_dp as TD
import make_figs as MF

df = D.build()
dd = df[df.obs.isin(['CS', 'CSCM', 'S3', 'S2X', 'S2Z'])]
P = D.driver_points(dd)
T = TD.Theory(P, tag='v1check', verbose=False)
assert not T.hasVC
main = T.main()
cv = dict(main=main, lo=main, hi=main)
ix = {p[0]: i for i, p in enumerate(P)}
t = MF.chi2_rows(df, cv, ix, C.E_DELTA_MIN, C.E_DELTA_MAX)
t['fit_report'] = [MF.FITREPORT.get(x, np.nan) for x in t.set]
L = ['| set | N | chi2/N (n fitted), this pipeline with v1 | fit_report Sect. 5.1 |', '|---|---|---|---|']
for r in t.itertuples():
    L.append('| %s | %d | %.2f | %s |' % (r.set, r.N, r.chi2N_nfit, '' if np.isnan(r.fit_report) else '%.2f' % r.fit_report))
open(os.path.join(C.RESULTS, 'check_v1_chi2.md'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
