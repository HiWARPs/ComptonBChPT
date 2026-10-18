"""
forward.py -- error model for the forward 'data points' dsigma/dOmega_lab(0) and Sigma_2z(0).

The points are derived from the empirical forward amplitudes Re f, Re g (db_forward_amp.csv):
  FWD_SOURCE = '2026' (default since 2026-10-18): the 2026 update of the ForwardRCS repository (github.com/HiWARPs/
     ForwardRCS; work/data/raw/ForwardRCS2026, built by work/data/scripts/build_db_forward.py);
  FWD_SOURCE = 'ghp2016' (environment variable): Gryniuk, Hagelstein, Pascalutsa, PRD 92 (2015) 074031;
     PRD 94 (2016) 034043 (db_forward_{amp,obs}_ghp2016.csv = the files used until 2026-10-18).
    CS(0)  = K (f^2 + g^2),           K = 389.379/389.379^2*1000 = 2.5682 nb per (mub GeV)^2
    S2z(0) = -2 f g / (f^2 + g^2)     (sigma_3/2 - sigma_1/2 convention = driver S2Z convention, checked)
Their uncertainties stem from the same fits of sigma_tot and Delta sigma at all energies and are therefore
strongly correlated.  Model ('nuisance'):
    delta f_i = sqrt(df_i^2 - r_f,i^2) * lam_f + (uncorrelated part r_f,i),   same for g with lam_g,
    r_i = err_unc column of db_forward_amp.csv (2026: remainder after the leading exp+fit eigenmode (+) systematic
          mode), or, if the column is absent (GHP 2016), r = 0.005/sqrt(3) mub GeV (uniform rounding of the 2-decimal
          table values),
lam_f, lam_g ~ N(0,1) common to all energies (two nuisance parameters).  The uncorrelated part gives a 2x2
covariance block (CS_i, S2Z_i) per energy.  Variant 'uncorr': diagonal errors as listed in db_forward_obs.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'data')
K = 389379.0 / 389.379 ** 2
R_ROUND = 0.005 / np.sqrt(3)
SOURCE = os.environ.get('FWD_SOURCE', '2026')
assert SOURCE in ('2026', 'ghp2016'), SOURCE
SUFFIX = '' if SOURCE == '2026' else '_ghp2016'          # db_forward_{amp,obs}<SUFFIX>.csv


def amplitudes(unc=False):
    """E, f, df, g, dg (df, dg = total errors); unc=True also returns the uncorrelated parts rf, rg"""
    d = pd.read_csv(os.path.join(DATA, 'db_forward_amp%s.csv' % SUFFIX))
    if 'err_unc' not in d:
        d['err_unc'] = R_ROUND
    f = d[d.obs == 'FWD_REF'].set_index('E_lab_MeV')
    g = d[d.obs == 'FWD_REG'].set_index('E_lab_MeV')
    E = f.index.values.astype(float)
    out = (E, f['value'].values, f['stat_err'].values, g.loc[E, 'value'].values, g.loc[E, 'stat_err'].values)
    return out + (f['err_unc'].values, g.loc[E, 'err_unc'].values) if unc else out


def obs(fv, gv):
    return K * (fv ** 2 + gv ** 2), -2 * fv * gv / (fv ** 2 + gv ** 2)


def derivs(fv, gv):
    s = fv ** 2 + gv ** 2
    dCS_df, dCS_dg = 2 * K * fv, 2 * K * gv
    dS_df = -2 * gv / s + 4 * fv ** 2 * gv / s ** 2
    dS_dg = -2 * fv / s + 4 * fv * gv ** 2 / s ** 2
    return dCS_df, dCS_dg, dS_df, dS_dg


def model(df_rows, mode='nuisance'):
    """df_rows: fit-db rows (in fit order) -- returns (cov_block (n,n), shifts (2,n) or None, idx_fwd)
    for the forward rows among df_rows (positions in df_rows). Non-forward rows untouched."""
    E, f, df, g, dg, rf, rg = amplitudes(unc=True)
    pos = {}
    for j, r in enumerate(df_rows.itertuples()):
        if r.dataset in ('GHP_fwd_CS', 'GHP_fwd_S2Z'):
            pos[(r.dataset, float(r.E))] = j
    n = len(df_rows)
    cov = np.zeros((n, n))
    shifts = np.zeros((2, n))
    a, b, c, d = derivs(f, g)
    for i, e in enumerate(E):
        jc = pos.get(('GHP_fwd_CS', e)); js = pos.get(('GHP_fwd_S2Z', e))
        if mode == 'uncorr':
            for j in (jc, js):
                if j is not None:
                    cov[j, j] = df_rows.iloc[j]['err_pt'] ** 2
            continue
        # uncorrelated part (GHP: rounding): J diag(rf^2, rg^2) J^T
        Ju = {}
        if jc is not None:
            Ju[jc] = np.array([a[i] * rf[i], b[i] * rg[i]])
        if js is not None:
            Ju[js] = np.array([c[i] * rf[i], d[i] * rg[i]])
        for j1, v1 in Ju.items():
            for j2, v2 in Ju.items():
                cov[j1, j2] += v1 @ v2
        cf = np.sqrt(max(df[i] ** 2 - rf[i] ** 2, 0)); cg = np.sqrt(max(dg[i] ** 2 - rg[i] ** 2, 0))
        if jc is not None:
            shifts[0, jc] = a[i] * cf; shifts[1, jc] = b[i] * cg
        if js is not None:
            shifts[0, js] = c[i] * cf; shifts[1, js] = d[i] * cg
    idx = sorted(pos.values())
    if mode == 'uncorr':
        return cov, None, idx
    return cov, shifts, idx


def check_table():
    """compare the error propagation with the errors listed in db_forward_obs.csv"""
    E, f, df, g, dg = amplitudes()
    cs, s2 = obs(f, g)
    a, b, c, d = derivs(f, g)
    ecs = np.sqrt((a * df) ** 2 + (b * dg) ** 2)
    es = np.sqrt((c * df) ** 2 + (d * dg) ** 2)
    return pd.DataFrame(dict(E=E, CS=cs, dCS=ecs, S2Z=s2, dS2Z=es))
