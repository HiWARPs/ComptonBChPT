"""Decomposition of the change of the NNLO fit (F2, main database, truncation covariance) caused by the switch of the
pre-2001 data to Baranov et al. 2001, Table 2 (2026-10-19).  Starting from the new database, the changes are undone
one at a time (cumulatively): (1) MacGibbon untagged full error matrix -> diagonal; (2) Baranov 74, 150 deg,
86.3 MeV error 1.0 -> 2.0 nb/sr; (3) the normalisation uncertainties and the Zieger error model -> the previous ones.
The remaining difference to the previous fit (results backed up in work/_to_delete/pre_baranov2001_2026-10-19/) is due
to the corrected energies (Hyman, Pugh at 50 deg), the Bernardini values and the rounding of the Table-2 values.
Run from work/fits (after `make`).  Writes results/baranov2001_decomposition.json."""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import common
import s03_fits_L as S

df = S.df
out = {}


def fit(tag):
    r = S.do_fit(tag, S.ALL, 'F2', thcov=True, store=False, decompose=False)
    out[tag] = dict(N=r['N'], chi2=r['chi2'], dof=r['dof'], alpha=r['P'][0], beta=r['P'][1], dalpha=r['dP'][0], dbeta=r['dP'][1])


fit('new (Baranov 2001, MacGibbon full matrix)')
common.MG_COV = {k: (v if k[1] == k[2] else 0.0) for k, v in common.MG_COV.items()}
fit('(1) MacGibbon untagged diagonal')
m = (df.dataset == 'Baranov1974_75_th150') & (abs(df.E - 86.3) < 0.01)
df.loc[m, ['err_stat', 'err_pt']] = 2.0
fit('(2) + Baranov 74 86.3 MeV error 2.0')
OLD = {'Oxley1958': 0.10, 'Hyman1959': 0.10, 'Goldansky1960': 0.08, 'Bernardini1960': 0.10, 'Pugh1957': 0.10,
       'Baranov1974_75_th90': 0.08, 'Baranov1974_75_th150': 0.08, 'Federspiel1991': 0.02, 'Hallin1993_TabII': 0.036,
       'MacGibbon1995_tagged': 0.04, 'Zieger1992': 0.0}
for k, v in OLD.items():
    df.loc[df.dataset == k, 'norm'] = v
z = df.dataset == 'Zieger1992'                    # previous Zieger error model: stat (+) 6.4% / 4.3% per point, not floated
df.loc[z, 'err_pt'] = np.sqrt(df.loc[z, 'err_stat'] ** 2 + (np.where(df.loc[z, 'E'] < 100, 0.064, 0.043) * df.loc[z, 'value']) ** 2)
fit('(3) + previous normalisations, Zieger')
PREV = os.path.join(common.FITS, '..', '_to_delete', 'pre_baranov2001_2026-10-19', 'fits_results', 'fits_L.json')
if os.path.exists(PREV):                          # backup of the previous results (project only, not in the package)
    prev = json.load(open(PREV))['F2_d_all_thB']
    out['previous database (KLP18 energies/values)'] = dict(N=prev['N'], chi2=prev['chi2'], dof=prev['dof'], alpha=prev['P'][0],
                                                             beta=prev['P'][1], dalpha=prev['dP'][0], dbeta=prev['dP'][1])
for k, v in out.items():
    print('%-46s N=%d chi2=%6.1f/%d  alpha=%.2f(%.2f) beta=%.2f(%.2f)' % (k, v['N'], v['chi2'], v['dof'], v['alpha'], v['dalpha'], v['beta'], v['dbeta']))
common.save_json(out, 'baranov2001_decomposition.json')
