"""
theory.py -- Bayesian/EFT truncation-error estimate from the order-by-order results of the code.

Orders (delta counting, LMP2015 Table 1; delta^2 ~ m_pi/Lambda, delta ~ m_pi/Delta):
    LO   = O(p^2)       : Born graphs without anomalous magnetic moment (Thomson/Dirac)          y0
    NLO  = O(p^3)       : + Pauli (kappa) Born terms + pi0 pole + piN loops                        y2  (relative delta^2)
    NNLO = O(p^{7/2})   : + Delta pole + piDelta loops                                              y3  (relative delta^3)
Expansion: y_NNLO = y_ref (c0 + c2 delta^2 + c3 delta^3),  c2 = (y_NLO - y_LO)/(y_ref delta^2),
c3 = (y_NNLO - y_NLO)/(y_ref delta^3).  Truncation error (all omitted orders delta^4, delta^5, ..., summed in
quadrature, Furnstahl et al. PRC 92 (2015) 024005; Melendez et al. PRC 100 (2019) 044001):
    sigma_th = y_ref * cbar * delta^4 / sqrt(1 - delta^2),   cbar^2 = (c2^2 + c3^2)/2  ('rms', default)
                                                               cbar = max(|c2|,|c3|)    ('max')
Expansion parameter (default 'sqrtQ'): delta(E) = sqrt(max(m_pi, E_lab)/Lambda_b), Lambda_b = 650 MeV
(delta = 0.463 at and below E = m_pi); alternative 'fixed': delta = m_pi/Delta = 0.48.
Reference scale: y_ref = |y_LO| for cross sections (LO Born), y_ref = 1 for asymmetries (natural size; the LO
Sigma_2x,z are accidentally small).
Two schemes for the order-by-order input (y0, y2, y3):
  'A' (naive):   y0 = LO (Born, kappa=0), y2 = NLO, y3 = NNLO, each with its own static polarisabilities.
                 The LO->NLO->NNLO steps are then dominated by the static polarisabilities, i.e. by exactly the
                 O(p^4) LEC directions that are FITTED -> double counting; overestimates the error of a fit.
  'B' (matched, default for fits): all orders evaluated with the SAME six static polarisabilities (those of
                 NNLO) via the O(p^4) contact terms; Born (with kappa) + pi0 pole treated as exact:
                 y0 = Born+pi0+contact, y2 = NLO amplitude + contact, y3 = NNLO.  The differences then measure
                 only the energy dependence beyond the static polarisabilities (dispersive/quadrupole effects,
                 cusp, Delta tail), which the fit parameters cannot absorb.
Correlations (for the theory covariance in chi^2): Gaussian kernel in (E_lab, theta_cm) within the same
observable family, rho = exp(-dE^2/2 lE^2 - dth^2/2 lth^2), lE = 25 MeV, lth = 40 deg (default, from a GP
marginal-likelihood fit of the scheme-B coefficients, s04_theory_lengths.py); 0 between families.
"""
import numpy as np

MPI = 139.57
ASYM = ('S3', 'S2X', 'S2Z')
FAMILY = {'CS': 'CS', 'CSCM': 'CS', 'SPAR': 'SPAR', 'SPERP': 'SPERP', 'S3': 'S3', 'S2X': 'S2X', 'S2Z': 'S2Z'}


def delta_of(E, mode='sqrtQ', Lam=650.0, fixed=0.48):
    E = np.asarray(E, float)
    if mode == 'fixed':
        return np.full_like(E, fixed)
    return np.sqrt(np.maximum(MPI, E) / Lam)


def truncation(y0, y2, y3, obs, E, mode='sqrtQ', Lam=650.0, cbar='rms', fixed=0.48):
    obs = np.asarray(obs)
    d = delta_of(E, mode, Lam, fixed)
    yref = np.where(np.isin(obs, ASYM), 1.0, np.abs(y0))
    c2 = (y2 - y0) / (yref * d ** 2)
    c3 = (y3 - y2) / (yref * d ** 3)
    if cbar == 'max':
        cb = np.maximum(np.abs(c2), np.abs(c3))
    else:
        cb = np.sqrt(0.5 * (c2 ** 2 + c3 ** 2))
    sig = yref * cb * d ** 4 / np.sqrt(1 - d ** 2)
    return dict(sigma=sig, c2=c2, c3=c3, cbar=cb, delta=d, yref=yref)


def covariance(sig, obs, E, theta_cm, lE=25.0, lth=40.0):
    fam = np.array([FAMILY[o] for o in obs])
    E = np.asarray(E, float); th = np.asarray(theta_cm, float)
    same = fam[:, None] == fam[None, :]
    rho = np.exp(-0.5 * ((E[:, None] - E[None, :]) / lE) ** 2 - 0.5 * ((th[:, None] - th[None, :]) / lth) ** 2) * same
    return rho * np.outer(sig, sig)
