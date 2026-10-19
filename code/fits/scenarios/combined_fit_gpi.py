"""The main BChPT fit (scenarios/combined_fit.py) with the backward spin polarisability gamma_pi relaxed as well
(author question, 2026-10-19): gamma_pi free without a prior, shifted along U_GPI = (gE1E1 - 1/2, gM1M1 + 1/2) per unit,
so that gamma_0 and gE1M2, gM1E2 stay at their NNLO values (the direction of s10_fig_pred_S2.D_GPI).
Variants: as the main fit (w = 0.01, R_EM free); R_EM fixed; Delta weights 0.03 and 1; below threshold only (Delta weight 1e-8, R_EM fixed);
and, for comparison, the main fit with gamma_pi at NNLO.  Writes results/combined_fit_gpi.json."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import combined_fit as CF
from common import build_fit, save_json, gamma0, gammapi, PREF_NNLO

U_GPI = np.array([0, 0, -0.5, 0.5, 0, 0])
assert abs(gammapi(U_GPI) - 1) < 1e-12 and abs(gamma0(U_GPI)) < 1e-12
OUT = {}


def fit(name, ge_free=True, gpi_free=True, w=CF.W_DELTA, Emax=360.0):
    sD, sth = CF.rows_delta(Emax)
    s = pd.concat([CF.sL, sD])
    C = np.zeros((len(s), len(s))); C[:len(CF.sL), :len(CF.sL)] = CF.C_L
    e2 = sD.err_pt.values ** 2
    C[len(CF.sL):, len(CF.sL):] = np.diag(((e2 + sth ** 2) / w) - e2)
    idx = CF.T0.subset(s.pid.values)
    ex0 = {}
    if ge_free:
        ex0['gE'] = CF.GE0
    if gpi_free:
        ex0['dgpi'] = 0.0

    def tf(P, ex):
        Pe = P + ex.get('dgpi', 0.0) * U_GPI
        O = np.array([CF.G[g].predict_array(Pe, idx) for g in CF.GE_GRID])
        return CubicSpline(CF.grid, O, axis=0)(ex.get('gE', CF.GE0))
    R = build_fit(s, CF.T0, PREF_NNLO, ['alphaE1', 'betaM1'], thcov=C, theory_fn=tf, extra=ex0).run(decompose=False)
    g = R.extra.get('gE', (CF.GE0, 0.0)); dg = R.extra.get('dgpi', (0.0, 0.0))
    P = np.array(R.P) + dg[0] * U_GPI
    chiL = sum(v['chi2'] for k, v in R.perset.items() if k in set(CF.sL.dataset))
    out = dict(name=name, w=w, P=P, dP=R.dP, covP=R.covP, chi2=R.chi2, dof=R.dof, gE=g, REM=(CF.REM_PER_GE * g[0],
               abs(CF.REM_PER_GE) * g[1]), gpi=(gammapi(P), dg[1]), g0=gamma0(P), chi2_below150_sets=chiL,
               rho_ab=float(R.covP[0][1] / np.sqrt(R.covP[0][0] * R.covP[1][1])))
    OUT[name] = out
    print('%-22s a=%.2f(%.2f) b=%.2f(%.2f) a+b=%.2f a-b=%.2f gpi=%.2f(%.2f) REM=%.2f(%.2f)%% chi2=%.1f/%d '
          '(sets below 150: %.1f) gE=%.4f' % (name, P[0], R.dP[0], P[1], R.dP[1], P[0] + P[1], P[0] - P[1], *out['gpi'],
                                       *out['REM'], R.chi2, R.dof, chiL, g[0]))
    return out


if __name__ == '__main__':
    print('NNLO: gamma_pi = %.2f, gamma_0 = %.2f' % (gammapi(PREF_NNLO), gamma0(PREF_NNLO)))
    fit('main_gpiNNLO', gpi_free=False)
    fit('main_gpifree')
    fit('REMfixed_gpifree', ge_free=False)
    fit('w0.03_gpiNNLO', gpi_free=False, w=0.03)
    fit('w0.03_gpifree', w=0.03)
    fit('w1_gpiNNLO', gpi_free=False, w=1.0)
    fit('w1_gpifree', w=1.0)
    fit('below150_gpifree', ge_free=False, w=1e-8)
    fit('below150_gpiNNLO', ge_free=False, gpi_free=False, w=1e-8)
    save_json(OUT, 'combined_fit_gpi.json')
