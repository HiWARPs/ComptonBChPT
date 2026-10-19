"""Combined NNLO fit of alpha_E1, beta_M1 and R_EM to the data up to 360 MeV (the main BChPT fit of the paper since
2026-10-19, author decision): the main database below the pion-production threshold (angular data E <= 150 MeV minus
MAIN_DROP and the forward points, with the truncation covariance, exactly as F2_d_all_thB) and the Delta-region data
150 < E <= 360 MeV of the ten sets of Table "delta" (A2 Sigma_2x, Sigma_2z; LEGS Sigma_3 and dsigma; MAMI 75, 90, 131,
136 deg; SAL; TAPS 157-164 MeV).  The Delta-region points enter with the uncertainty band of the NNLO amplitude
there (the LMP15-type envelope with delta = omega/Lambda_chi of work/delta_plots, results/data_with_theory.csv, half
width (hi - lo)/2) added in quadrature, as in the "+th." column of Table "delta", and the region below 150 MeV carries
the weight: the Delta-region chi^2 is multiplied by W_DELTA = 0.01 (author decision; with W = 1 the Delta region
dominates and drives alpha, beta to 7.9, 6.5 through the omega^2 growth of the contact terms).
Free: alpha_E1, beta_M1 and the E2/M1 coupling g_E, R_EM = (M_D - M_N)/(2(M_D + M_N)) g_E/g_M (default g_E = -1.0,
R_EM = -2.28%; cubic-spline interpolation over a grid of complete NNLO runs, the static polarisabilities held fixed),
unconstrained in the main fit (author decision; variant: within +-15% of the default, applied as a bound);
g_M = 2.97.  The four spin polarisabilities are varied within the uncertainties of the NNLO prediction (author request,
2026-10-19): independent Gaussian priors at the NNLO values with the re-examined widths of
scenarios/spin_pol_uncertainties.py (results/spin_pol_uncertainties.json, 'reexamined'); variant 'combined_spinNNLO':
spin polarisabilities fixed at NNLO (the main fit until then).
Variants: R_EM fixed at the default; R_EM within +-15%; Delta weight 0.03 and 1; E <= 350 MeV.
The fits of Table "fits" (below 150 MeV, g_E at its default) are unchanged.
Run from work/fits after `make` and work/delta_plots/make_all.sh.  Writes results/combined_fit.json."""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from common import outputs, build_fit, save_json, gamma0, gammapi, PREF_NNLO, FITS, RES
import s03_fits_L as S

GE0, GE_WIN = -1.0, 0.15
W_DELTA = 0.01
GE_GRID = [-3.5, -3.0, -2.6, -2.3, -2.0, -1.8, -1.6, -1.4, -1.3, -1.2, -1.15, -1.1, -1.05, -1.0, -0.95, -0.9, -0.85, -0.8, -0.7, -0.6]
REM_PER_GE = 100 * (1232.0 - 938.272) / (2 * (1232.0 + 938.272)) / 2.97      # R_EM [%] per unit g_E at g_M = 2.97
DSETS = ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_S3', 'LEGS2001_CS', 'Hunger1997_75', 'Hunger1997_90',
         'Wissmann1999', 'Camen2002', 'Hallin1993_TabII', 'Hallin1993_TabII_D', 'OlmosdeLeon2001']
NAMES = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
SIG_SPIN = np.array(json.load(open(os.path.join(RES, 'spin_pol_uncertainties.json')))['reexamined'][2:6])
BAND = pd.read_csv(os.path.join(FITS, '..', 'delta_plots', 'results', 'data_with_theory.csv')).set_index('pid')

df = S.df
G = {g: outputs(df, 'gE%.2f' % g) for g in GE_GRID}
T0 = G[GE0]
grid = np.array(GE_GRID)
sL = S.rows_for(df, S.ALL, 'F2')                               # main database + forward points (as F2_d_all_thB)
C_L, _ = S.thcov_for(sL, df, 'db')
OUT = {}


def rows_delta(Emax):
    m = df.dataset.isin(DSETS) & (df.E > 150) & (df.E <= Emax)
    s = df[m].copy()
    sth = 0.5 * (BAND.loc[s.pid, 'th_hi'].values - BAND.loc[s.pid, 'th_lo'].values)
    assert np.all(np.isfinite(sth)), 'theory band missing for some Delta-region points'
    return s, sth


def fit(name, ge='free', band=True, Emax=360.0, w=W_DELTA, spin=True):
    sD, sth = rows_delta(Emax)
    s = pd.concat([sL, sD])
    C = np.zeros((len(s), len(s))); C[:len(sL), :len(sL)] = C_L
    e2 = sD.err_pt.values ** 2                                 # build_fit puts err_pt^2 on the diagonal
    C[len(sL):, len(sL):] = np.diag(((e2 + (sth ** 2 if band else 0.0)) / w) - e2)   # Delta chi^2 weighted by w
    idx = T0.subset(s.pid.values)
    free = NAMES if spin else NAMES[:2]
    cons = [(np.eye(6)[k], PREF_NNLO[k], SIG_SPIN[k - 2], 'prior_' + NAMES[k]) for k in range(2, 6)] if spin else ()

    def tf(P, ex):
        O = np.array([G[g].predict_array(P, idx) for g in GE_GRID])
        return CubicSpline(grid, O, axis=0)(ex['gE'])

    def run(gfix=None):
        if gfix is None:
            F = build_fit(s, T0, PREF_NNLO, free, constraints=cons, thcov=C, theory_fn=tf, extra=dict(gE=GE0))
        else:
            F = build_fit(s, T0, PREF_NNLO, free, constraints=cons, thcov=C, theory_fn=lambda P, ex: tf(P, dict(gE=gfix)))
        return F.run(decompose=False)
    edge = None
    if ge == 'fixed':
        R = run(GE0); g = (GE0, 0.0)
    else:
        R = run(); g = R.extra['gE']
        lo, hi = GE0 * (1 + GE_WIN), GE0 * (1 - GE_WIN)                 # -1.15 ... -0.85
        if ge == 'window' and not (lo <= g[0] <= hi):
            edge = lo if g[0] < lo else hi
            g_free = g
            R = run(edge); g = (edge, 0.0)
    P, cov = np.array(R.P), np.array(R.covP)
    der = {k: (float(np.dot(v, P)), float(np.sqrt(np.dot(v, cov @ np.array(v, float)))))
           for k, v in {'apb': [1, 1, 0, 0, 0, 0], 'amb': [1, -1, 0, 0, 0, 0], 'g0': [0, 0, -1, -1, -1, -1],
                        'gpi': [0, 0, -1, 1, -1, 1]}.items()}
    chiL = sum(v['chi2'] for k, v in R.perset.items() if k in set(sL.dataset))
    out = dict(name=name, spin=spin, N=len(s), N_L=len(sL), N_D=len(sD), Emax=Emax, band=band, w=w, P=P, dP=R.dP, covP=cov, chi2=R.chi2,
               dof=R.dof, p=R.pvalue, chi2_data=R.chi2_data, perset=R.perset, der=der, gE=g,
               REM=(REM_PER_GE * g[0], abs(REM_PER_GE) * g[1]), edge=edge, chi2_below150_sets=chiL,
               rho_ab=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])))
    if edge is not None:
        out['gE_free'] = g_free; out['REM_free'] = (REM_PER_GE * g_free[0], abs(REM_PER_GE) * g_free[1])
    OUT[name] = out
    print('%-24s N=%d (%d+%d) a=%.2f(%.2f) b=%.2f(%.2f) a+b=%.2f(%.2f) a-b=%.2f(%.2f) spin=%s g0=%.2f(%.2f) gpi=%.2f(%.2f) '
          'gE=%.3f(%.3f) REM=%.2f(%.2f)%%%s chi2=%.1f/%d' % (name, len(s), len(sL), len(sD), P[0], R.dP[0], P[1], R.dP[1],
                            *der['apb'], *der['amb'], np.round(P[2:], 2).tolist(), *der['g0'], *der['gpi'],
                            g[0], g[1], *out['REM'], (' [edge; free: REM=%.2f(%.2f)%%]' % out['REM_free']) if edge is not None else '',
                            R.chi2, R.dof))
    return out


if __name__ == '__main__':
    fit('combined')                                  # MAIN: alpha, beta, R_EM free, spin pols within NNLO priors,
    fit('combined_spinNNLO', spin=False)             #   E <= 360 MeV, Delta band, weight 0.01
    fit('combined_REMfixed', ge='fixed')
    fit('combined_REMwindow', ge='window')
    fit('combined_w0.03', w=0.03)
    fit('combined_w1', w=1.0)
    fit('combined_E350', Emax=350.0)
    for k in ('combined', 'combined_spinNNLO', 'combined_REMfixed'):
        r = OUT[k]
        print('per set (%s):' % k, {ds: round(v['chi2'], 1) for ds, v in r['perset'].items()})
    save_json(OUT, 'combined_fit.json')
