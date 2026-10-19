"""forward_amp.py -- the forward sum rules (forward SRs) above the pion-production threshold at the level of the
amplitudes, for the main BChPT fit (author request 2026-10-19: forward SRs up to 360 MeV, fitted as amplitudes).

Data: Re f, Im f, Re g, Im g of the ForwardRCS 2026 update at E = 160, 165, ..., 360 MeV (5-MeV grid; load() selects the fit grid) with the complete exp+fit
covariance and the signed systematic modes (data selection of f, PWA model of g) -- work/data/scripts/
export_forward_amp_grid.py -> work/data/raw/ForwardRCS2026/forward_amp_grid.json (mub GeV).  The one-sided elastic-channel
shift of f is not included (as below threshold).
BChPT: f, g extracted EXACTLY from the fit driver's quadratic forms of dsigma/dOmega_lab(0) and Sigma_2z(0) (the 'exact'
route of work/pwa/scenarios/nnlo_forward.py, same formulas); in the forward direction the contact terms give exactly
delta f = delta(alpha+beta) nu^2 and delta g = delta(gamma_0) nu^3 (lab nu), so f, g depend on the static
polarisabilities only through alpha+beta and gamma_0.
Theory band: the LMP15-type envelope of work/delta_plots (shift sets piN, piDelta, vertex corrections, their sum, each by
+-delta = omega_cm/Lambda_chi, Lambda_chi = 750 MeV) applied to the amplitudes, which are linear in the mechanisms:
sigma_X = delta max_e |X_e|; the mechanism amplitudes from differences of complete driver runs at their own static values
(NLO - BornPi0 = piN loops; NNLO_noDpole - NLO = piDelta loops; NNLO - NNLO without vertex corrections = VC)."""
import os, json
import numpy as np
import bchpt
from fitter import C_BALDIN, C_GAMMA0

HERE = os.path.dirname(os.path.abspath(__file__))
GRID_JSON = os.path.join(HERE, '..', '..', 'data', 'raw', 'ForwardRCS2026', 'forward_amp_grid.json')
UBGEV_PER_FM = 1973.27                    # 1 fm = 1973.27 mub GeV
K_NB = 1.0e7                              # nb/sr per fm^2
HBARC, MP = 197.33, 938.272                 # MeV fm (the Fortran value), MeV
E_PI0 = 144.7                             # pi0 p threshold (code masses)
LAMBDA_CHI = 750.0
PIDX = {p: i for i, p in enumerate([(a, b) for a in range(7) for b in range(a, 7)])}
COMP = ['Re f', 'Im f', 'Re g', 'Im g']


def load(emax=360.0, step=10.0):
    """empirical amplitudes on the grid (E <= emax, every `step` MeV): E, y (Re f, Im f, Re g, Im g blocks), covariance"""
    J = json.load(open(GRID_JSON))
    E = np.array(J['E'])                      # 5-MeV grid 160-360 MeV
    keep = (E <= emax + 1e-9) & (np.abs(((E - E[0]) / step) - np.round((E - E[0]) / step)) < 1e-9)
    n = len(E)
    sel2 = np.r_[keep, keep]
    Cf = np.array(J['cov_f']) + np.outer(J['sel_f'], J['sel_f'])
    Cg = np.array(J['cov_g']) + np.outer(J['pwa_g'], J['pwa_g'])
    Cf, Cg = Cf[np.ix_(sel2, sel2)], Cg[np.ix_(sel2, sel2)]
    m = keep.sum()
    C = np.zeros((4 * m, 4 * m)); C[:2 * m, :2 * m] = Cf; C[2 * m:, 2 * m:] = Cg
    y = np.r_[np.array(J['re_f'])[keep], np.array(J['im_f'])[keep], np.array(J['re_g'])[keep], np.array(J['im_g'])[keep]]
    return E[keep], y, C


def points(E, tag='FSR'):
    return [('%s_CS_%03d' % (tag, i), float(e), 0.0, 'CS', 'LAB') for i, e in enumerate(E)] + \
           [('%s_S2Z_%03d' % (tag, i), float(e), 0.0, 'S2Z', 'LAB') for i, e in enumerate(E)]


def nu_fm(E):
    return np.asarray(E, float) / HBARC


def fg(T, E):
    """f, g [mub GeV, complex] at T.Pref from a driver output with the rows of points(E) (CS rows, then S2Z rows)"""
    E = np.asarray(E, float); n = len(E)
    assert np.allclose(T.E[:n], E) and np.allclose(T.E[n:], E) and all(T.obs[:n] == 'CS') and all(T.obs[n:] == 'S2Z')
    CN, SN, SD = T.CN[:n], T.CN[n:], T.CD[n:]
    c = {k: CN[:, PIDX[k]] for k in ((0, 0), (0, 1), (1, 1), (0, 3), (3, 3))}
    u, w = np.sqrt(c[(1, 1)] / K_NB), np.sqrt(c[(3, 3)] / K_NB)          # |df/dalpha|, |dg/dgamma| (fm per unit)
    nu = nu_fm(E)
    assert np.allclose(u, nu ** 2 * 1e-4, rtol=1e-6) and np.allclose(w, nu ** 3 * 1e-4, rtol=1e-6), 'contact-term check'
    Ref = c[(0, 1)] / (2 * K_NB * u)
    Reg = -c[(0, 3)] / (2 * K_NB * w)
    Refg = -0.5 * (SN[:, 0] / SD[:, 0]) * c[(0, 0)] / K_NB
    S = c[(0, 0)] / K_NB - Ref ** 2 - Reg ** 2                         # Im f^2 + Im g^2
    Pp = Refg - Ref * Reg                                               # Im f Im g
    disc = S ** 2 - 4 * Pp ** 2
    above = E >= E_PI0
    tot = c[(0, 0)] / K_NB                                              # |f|^2 + |g|^2 (scale of the rounding)
    assert np.all(disc[above] > -1e-3 * S[above] ** 2 - 1e-9 * tot[above] ** 2), 'inconsistent imaginary parts'
    Imf = np.sqrt(np.maximum((S + np.sqrt(np.maximum(disc, 0))) / 2, 0))   # |Im f| > |Im g|, Im f = nu sigma_T/4pi > 0
    Img = np.where(Imf > 0, Pp / np.where(Imf > 0, Imf, 1.0), 0.0)
    f = Ref + 1j * np.where(above, Imf, 0.0); g = Reg + 1j * np.where(above, Img, 0.0)
    rec = K_NB * (abs(f) ** 2 + abs(g) ** 2)
    assert np.allclose(rec, T.predict_array(T.Pref)[:n], rtol=1e-5), 'dsigma(0) reconstruction'
    return f * UBGEV_PER_FM, g * UBGEV_PER_FM


def vec(f, g):
    return np.r_[f.real, f.imag, g.real, g.imag]


def shift(E, P, Pref):
    """contact-term shift of (Re f, Im f, Re g, Im g) for static polarisabilities P relative to Pref (mub GeV)"""
    nu = nu_fm(E)
    df = (C_BALDIN @ (np.asarray(P) - Pref)) * nu ** 2 * 1e-4 * UBGEV_PER_FM
    dg = (C_GAMMA0 @ (np.asarray(P) - Pref)) * nu ** 3 * 1e-4 * UBGEV_PER_FM
    z = np.zeros_like(nu)
    return np.r_[df, z, dg, z]


def band(E, tag='fsr'):
    """half width of the truncation envelope of (Re f, Im f, Re g, Im g) [mub GeV] and the mechanism amplitudes"""
    pts = points(E)
    R = {c: vec(*fg(bchpt.run(pts, c, name='%s_%s' % (tag, c if isinstance(c, str) else 'noVC'), verbose=False), E))
         for c in ('NNLO', 'NLO', 'BornPi0', 'NNLO_noDpole')}
    R['noVC'] = vec(*fg(bchpt.run(pts, dict(dvertex=False), name='%s_noVC' % tag, verbose=False), E))
    mech = dict(piN=R['NLO'] - R['BornPi0'], piD=R['NNLO_noDpole'] - R['NLO'], VC=R['NNLO'] - R['noVC'])
    mech['sum'] = mech['piN'] + mech['piD'] + mech['VC']
    Ecm = np.asarray(E) * MP / np.sqrt(MP ** 2 + 2 * MP * np.asarray(E))
    d = np.tile(Ecm / LAMBDA_CHI, 4)
    return d * np.max(np.abs(np.array(list(mech.values()))), axis=0), mech, R['NNLO']
