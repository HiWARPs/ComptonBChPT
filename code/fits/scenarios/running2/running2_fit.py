"""running2_fit.py -- effect of the running Delta parameters (code_v2 &v2 irun = 2 doubly subtracted, 3 = + slopes at s0;
work/notes/delta_running2.md) on the NNLO fit of the paper: F2_d_all_thB (main database, GHP forward points to 145 MeV as data
with correlated nuisances, truncation covariance with the default settings, Lambda = 750 MeV), repeated with the base amplitude
'running2' / 'running3' (fits library configs; physical masses and vertex corrections as in the pipeline).
Truncation covariance (scheme B, LEC-matched): the orders Born+pi0 and NLO contain no Delta pole and are unchanged by the
running; y3 is the running amplitude, all three orders evaluated at the same six static polarisabilities P0 of the fit.
Variants per scheme: (i) P0 = own Pref of the running amplitude (irun=2: gamma_M1M1 = 2.56; irun=3: = default), covariance
recomputed; (ii) P0 = default Pref (contact terms restore the NNLO gammas), covariance recomputed; (iii) as (i) with the
default-configuration covariance; (iv) without theory covariance.  Reproduction of results/fits_L.json for the default first.
Run: .venv/bin/python running2_fit.py  (~2 min)  -> running2_fit.{json,md}"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
FITS = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(FITS, 'scripts')); sys.path.insert(0, os.path.join(FITS, 'lib'))
os.chdir(FITS)
import numpy as np
import pandas as pd
import bchpt
import theory as TH
from common import load_db, outputs, PREF_NNLO, RES
import s03_fits_L as S

df = S.df
T0 = S.T
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
thcov_std = S.thcov_for


def conv(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))


def make_thcov(cfg, Pmatch):
    """truncation covariance with y3 = running amplitude cfg, all orders at the static polarisabilities Pmatch"""
    def order_values(dfx, s):
        T3 = outputs(dfx, cfg)
        idx = T3.subset(s['pid'].values)
        y3 = T3.predict_array(np.asarray(Pmatch, float), idx)
        y0 = outputs(dfx, 'BornPi0').predict_array(np.asarray(Pmatch, float), idx)
        y2 = outputs(dfx, 'NLO').predict_array(np.asarray(Pmatch, float), idx)
        return y0, y2, y3

    def thcov(s, dfx, tag, mode='sqrtQ', Lam=750.0, cbar='rms', corr=True, lE=25.0, lth=40.0):
        y0, y2, y3 = order_values(dfx, s)
        tr = TH.truncation(y0, y2, y3, s['dobs'].values, s['E'].values, mode=mode, Lam=Lam, cbar=cbar)
        if corr:
            Cm = TH.covariance(tr['sigma'], s['dobs'].values, s['E'].values, s['theta_cm'].values, lE=lE, lth=lth)
        else:
            Cm = np.diag(tr['sigma'] ** 2)
        tr['y0'], tr['y2'], tr['y3'] = y0, y2, y3
        return Cm, tr
    return thcov


def strip(r):
    keep = ['name', 'kind', 'N', 'free', 'P', 'dP', 'covP', 'chi2', 'dof', 'p', 'chi2_data', 'chi2_norm', 'chi2_con', 'chi2_shift',
            'perset', 'g0', 'gpi', 'apb', 'dapb', 'amb', 'damb', 'rho_ab', 'dP_stat', 'dP_statnorm', 'dP_fwd', 'dP_norm', 'lam_shift']
    return {k: r[k] for k in keep if k in r}


fits = {}
print('=== default: reproduction of fits_L.json F2_d_all_thB ===')
r0 = S.do_fit('F2_d_all_thB', S.ALL, 'F2', thcov=True, store=False)
ref = FL['F2_d_all_thB']
print('    max |dP| = %.1e, dchi2 = %.1e' % (np.max(np.abs(np.array(r0['P']) - np.array(ref['P']))), r0['chi2'] - ref['chi2']))
fits['F2_d_all_thB'] = strip(r0)
fits['F2_d_all'] = strip(S.do_fit('F2_d_all', S.ALL, 'F2', store=False))
P0d = PREF_NNLO.copy()
for cfg in ('running2', 'running3'):
    Tx = outputs(df, cfg)
    assert list(Tx.ids) == list(T0.ids)
    print('=== %s: Pref = %s (default %s) ===' % (cfg, np.round(Tx.Pref, 4), np.round(P0d, 4)))
    fits[cfg + '_Pref'] = Tx.Pref.tolist()
    # (i) own Pref, covariance recomputed
    S.thcov_for = make_thcov(cfg, Tx.Pref)
    fits['F2_d_all_thB_%s_own' % cfg] = strip(S.do_fit('F2_d_all_thB_%s_own' % cfg, S.ALL, 'F2', thcov=True, TT=Tx, P0=Tx.Pref.copy(), store=False))
    # (ii) default P0 (contact terms restore the NNLO gammas), covariance recomputed at P0
    S.thcov_for = make_thcov(cfg, P0d)
    fits['F2_d_all_thB_%s_P0' % cfg] = strip(S.do_fit('F2_d_all_thB_%s_P0' % cfg, S.ALL, 'F2', thcov=True, TT=Tx, P0=P0d.copy(), store=False))
    # (iii) own Pref with the default-configuration covariance
    S.thcov_for = thcov_std
    fits['F2_d_all_thB_%s_own_covdef' % cfg] = strip(S.do_fit('F2_d_all_thB_%s_own_covdef' % cfg, S.ALL, 'F2', thcov=True, TT=Tx, P0=Tx.Pref.copy(), store=False))
    # (iv) no theory covariance
    fits['F2_d_all_%s_own' % cfg] = strip(S.do_fit('F2_d_all_%s_own' % cfg, S.ALL, 'F2', TT=Tx, P0=Tx.Pref.copy(), store=False))
    fits['F2_d_all_%s_P0' % cfg] = strip(S.do_fit('F2_d_all_%s_P0' % cfg, S.ALL, 'F2', TT=Tx, P0=P0d.copy(), store=False))
    # size of the truncation-covariance change
    s_ = S.rows_for(df, S.ALL, 'F2')
    C0, tr0 = thcov_std(s_, df, 'db'); C1, tr1 = make_thcov(cfg, Tx.Pref)(s_, df, 'db')
    sig0, sig1 = np.sqrt(np.diag(C0)), np.sqrt(np.diag(C1))
    fits[cfg + '_trunc_sigma_change'] = dict(max_rel=float(np.max(np.abs(sig1 / sig0 - 1))), mean_rel=float(np.mean(sig1 / sig0 - 1)))
    print('    truncation sigma_th with %s: max rel. change %.2e, mean %.2e' % (cfg, fits[cfg + '_trunc_sigma_change']['max_rel'],
                                                                                 fits[cfg + '_trunc_sigma_change']['mean_rel']))
json.dump(fits, open(os.path.join(HERE, 'running2_fit.json'), 'w'), indent=1, default=conv)

L = ['# Running Delta parameters: NNLO fit (F2_d_all_thB type)', '',
     'Produced by `running2_fit.py`.  alpha, beta [1e-4 fm^3]; gamma_M1M1 held at P0[3]; chi^2/dof.',
     "Suffixes: `_own` = P0 = own Pref of the running amplitude, covariance recomputed; `_P0` = P0 = default Pref (contact terms",
     "restore the NNLO gammas), covariance recomputed; `_own_covdef` = default-configuration covariance; `F2_d_all_*` = no theory covariance.", '',
     '| fit | N | alpha_E1 | beta_M1 | alpha+beta | alpha-beta | gamma_M1M1 (fixed) | chi^2/dof | d alpha | d beta | d chi^2 |',
     '|---|---|---|---|---|---|---|---|---|---|---|']
order = ['F2_d_all_thB'] + ['F2_d_all_thB_%s_%s' % (c, v) for c in ('running2', 'running3') for v in ('own', 'P0', 'own_covdef')] + \
        ['F2_d_all'] + ['F2_d_all_%s_%s' % (c, v) for c in ('running2', 'running3') for v in ('own', 'P0')]
for k in order:
    r = fits[k]; b = fits['F2_d_all_thB'] if 'thB' in k else fits['F2_d_all']
    L.append('| %s | %d | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f | %.2f/%d | %+.3f | %+.3f | %+.2f |' % (
        k, r['N'], r['P'][0], r['dP'][0], r['P'][1], r['dP'][1], r['apb'], r['dapb'], r['amb'], r['damb'], r['P'][3],
        r['chi2'], r['dof'], r['P'][0] - b['P'][0], r['P'][1] - b['P'][1], r['chi2'] - b['chi2']))
for cfg in ('running2', 'running3'):
    t = fits[cfg + '_trunc_sigma_change']
    L.append('')
    L.append('%s: truncation sigma_th recomputed: max relative change %.1e (mean %.1e).' % (cfg, t['max_rel'], t['mean_rel']))
open(os.path.join(HERE, 'running2_fit.md'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
