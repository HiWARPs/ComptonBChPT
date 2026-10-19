"""The main BChPT fit (scenarios/combined_fit.py: alpha_E1, beta_M1 and R_EM up to 360 MeV, Delta region weighted by
0.01) with the four spin polarisabilities varied as well, within the uncertainties of the NNLO prediction (author
request, 2026-10-19): Gaussian priors centred at the NNLO values (physical masses) with the widths of
scenarios/spin_pol_uncertainties.py ('reexamined'; variant: the published LMP15 widths 0.8, 1.5, 0.2, 0.3), the four
priors independent.  Variants: R_EM fixed; below threshold only (Delta weight 1e-8, R_EM fixed); spin pols at NNLO
(= the current main fit).  Writes results/combined_fit_spin.json."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import combined_fit as CF
from common import build_fit, save_json, gamma0, gammapi, PREF_NNLO, RES

U = json.load(open(os.path.join(RES, 'spin_pol_uncertainties.json')))
SIG = {'reexamined': np.array(U['reexamined'][2:6]), 'LMP15': np.array([0.8, 1.5, 0.2, 0.3])}
NAMES = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
OUT = {}


def fit(name, prior='reexamined', ge_free=True, w=CF.W_DELTA, Emax=360.0):
    sD, sth = CF.rows_delta(Emax)
    s = pd.concat([CF.sL, sD])
    C = np.zeros((len(s), len(s))); C[:len(CF.sL), :len(CF.sL)] = CF.C_L
    e2 = sD.err_pt.values ** 2
    C[len(CF.sL):, len(CF.sL):] = np.diag(((e2 + sth ** 2) / w) - e2)
    idx = CF.T0.subset(s.pid.values)

    def tf(P, ex):
        O = np.array([CF.G[g].predict_array(P, idx) for g in CF.GE_GRID])
        return CubicSpline(CF.grid, O, axis=0)(ex.get('gE', CF.GE0))
    if prior is None:
        free, cons = NAMES[:2], ()
    else:
        free = NAMES
        cons = [(np.eye(6)[k], PREF_NNLO[k], SIG[prior][k - 2], 'prior_' + NAMES[k]) for k in range(2, 6)]
    R = build_fit(s, CF.T0, PREF_NNLO, free, constraints=cons, thcov=C, theory_fn=tf,
                  extra=dict(gE=CF.GE0) if ge_free else {}).run(decompose=False)
    g = R.extra.get('gE', (CF.GE0, 0.0))
    P, cov = np.array(R.P), np.array(R.covP)
    lin = lambda c: (float(c @ P), float(np.sqrt(c @ cov @ c)))
    chiL = sum(v['chi2'] for k, v in R.perset.items() if k in set(CF.sL.dataset))
    out = dict(name=name, prior=prior, w=w, P=P, dP=R.dP, covP=cov, chi2=R.chi2, dof=R.dof, gE=g,
               REM=(CF.REM_PER_GE * g[0], abs(CF.REM_PER_GE) * g[1]), g0=lin(np.r_[0, 0, -1, -1, -1, -1.]),
               gpi=lin(np.r_[0, 0, -1, 1, -1, 1.]), apb=lin(np.r_[1, 1, 0, 0, 0, 0.]), amb=lin(np.r_[1, -1, 0, 0, 0, 0.]),
               chi2_below150_sets=chiL, rho_ab=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])),
               chi2_priors=float(sum(((P[k] - PREF_NNLO[k]) / SIG[prior][k - 2]) ** 2 for k in range(2, 6))) if prior else 0.0)
    OUT[name] = out
    print('%-22s a=%.2f(%.2f) b=%.2f(%.2f) a+b=%.2f(%.2f) a-b=%.2f(%.2f) | gEE=%.2f(%.2f) gMM=%.2f(%.2f) gEM=%.2f(%.2f) '
          'gME=%.2f(%.2f) g0=%.2f(%.2f) gpi=%.2f(%.2f) | REM=%.2f(%.2f)%% chi2=%.1f/%d (sets<150 %.1f, priors %.1f)'
          % (name, P[0], R.dP[0], P[1], R.dP[1], *out['apb'], *out['amb'], P[2], R.dP[2], P[3], R.dP[3], P[4], R.dP[4],
             P[5], R.dP[5], *out['g0'], *out['gpi'], *out['REM'], R.chi2, R.dof, chiL, out['chi2_priors']))
    return out


if __name__ == '__main__':
    print('NNLO spin pols', np.round(PREF_NNLO[2:], 3), 'gamma0 %.2f gammapi %.2f' % (gamma0(PREF_NNLO), gammapi(PREF_NNLO)))
    print('prior widths', {k: np.round(v, 3).tolist() for k, v in SIG.items()})
    fit('spinNNLO', prior=None)
    fit('spin_reex')
    fit('spin_LMP15', prior='LMP15')
    fit('spin_reex_REMfixed', ge_free=False)
    fit('spin_reex_w0.03', w=0.03)
    fit('spin_reex_below150', ge_free=False, w=1e-8)
    save_json(OUT, 'combined_fit_spin.json')
