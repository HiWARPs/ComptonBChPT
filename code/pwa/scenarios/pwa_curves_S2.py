"""Curves of the PWA fit of the paper (scenarios/main_fits_legacynorm.py, fit 'MAIN' -- the same fit and the same
evaluation as scenarios/pwa_curves.py for the paper's Figs. cs, sigma3, higs) for the double-polarisation asymmetries
    S2X  Sigma_2x   (circularly polarised beam, target spin in the scattering plane, perpendicular to the beam)
    S2Z  Sigma_2z   (target spin along the beam)
in the conventions of the A2 data and of the Fortran driver (compton.S2_from_phis, validated in validation/validate_S2.py),
at the panel energies of paper/figs/fig_pred_S2_abs.pdf (work/fits/scripts/s11_fig_pred_S2_abs.py): E = 100, 120, 140 MeV,
lab angles 0-180 deg in 2.5-deg steps.  Writes results/pwa_curves_S2.csv with the columns
    obs, E, theta_lab, value, err
value = Model.curve(E, theta, ib, p, obs) = exact Born + pi0 amplitude plus the non-Born multipoles (l = 1 fitted, l = 2 from
NNLO BChPT) of the bin nearest to E (bins 99, 117, 143 MeV; energy_mode='actual', as pwa_curves.py); err = linear
propagation of the fit covariance of the multipole parameters (cov_L_MAIN.npy, central differences in the parameters;
the observable is a ratio of quadratic forms, so this is essentially exact), statistical only (no 'method' variants).
Also printed (not written): the difference with the 'main' fit of make_figures.fits() (same data, the pre-1961 sets with
fixed normalisation), and the asymmetries at the bin centres.
Run from anywhere (~3 s):  .venv/bin/python work/pwa/scenarios/pwa_curves_S2.py"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
os.chdir(PWA)
sys.path.insert(0, PWA)
sys.path.insert(0, HERE)
import update_fits as UF                      # noqa: E402
import main_fits_legacynorm as MFL            # noqa: E402
import compton as C                           # noqa: E402

OUT = os.path.join(PWA, 'results', 'pwa_curves_S2.csv')
ENERGIES = [100.0, 120.0, 140.0]
TH = np.arange(0.0, 180.01, 2.5)
OBS = ['S2X', 'S2Z']


def nearest_bin(m, E):
    return int(np.argmin(abs(m.bins - E)))


def curve_err_fd(m, E, th, ib, p, cov, obs, h=1e-3):
    """value and 1-sigma error of Model.curve(E, th, ib, p, obs): numerical Jacobian w.r.t. the multipole parameters
    p[:npm] (central differences, step h), propagated with cov[:npm, :npm]."""
    npm = m.ng + m.nr * m.nb
    v0 = m.curve(E, th, ib, p, obs=obs)
    J = np.zeros((len(th), npm))
    for j in range(npm):
        pp, pm_ = p.copy(), p.copy(); pp[j] += h; pm_[j] -= h
        J[:, j] = (m.curve(E, th, ib, pp, obs=obs) - m.curve(E, th, ib, pm_, obs=obs)) / (2 * h)
    err = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', J, cov[:npm, :npm], J), 0))
    return v0, err


def main():
    d = MFL.load_leg()
    m = UF.build_model(UF.select(d), UF.BINS_A, **MFL.MAIN)
    p = np.load('results/update/p_L_MAIN.npy')
    cov = np.load('results/update/cov_L_MAIN.npy')
    J = json.load(open('results/main_fits_legacynorm.json'))['MAIN']
    assert p.shape == (m.npar,) and cov.shape == (m.npar, m.npar), 'parameter layout differs from the saved fit'
    chi2 = m.chi2_data(p)
    print('PWA main fit (paper): N = %d, chi2 = %.4f (json %.4f), npar = %d' % (int(m.mask.sum()), chi2, J['chi2'], m.npar))
    assert abs(chi2 - J['chi2']) < 1e-6, 'saved parameters do not reproduce the chi2 of the main fit'
    # consistency of the helper with the model: S3 and cs_lab of Model.curve are unchanged by the new keys (regression)
    k0 = C.kin(100.0, TH, 'lab')
    o = C.obs_from_phis(C.phis_from_M(C.born_amplitudes(100.0, np.degrees(k0['thcm']), 'cm')), k0)
    assert set(['cs_lab', 'S3', 'cs_par_lab', 'cs_perp_lab', 'S2X', 'S2Z']) <= set(o)
    # the same fit without the legacy normalisations floated (make_figures.fits()['main']) for comparison
    from main_fits import MAIN as MAIN0
    m0 = UF.build_model(UF.select(UF.load()), UF.BINS_A, **MAIN0)
    r0 = m0.fit()
    print("make_figures 'main' fit: N = %d, chi2 = %.4f" % (r0.npts, r0.chi2))
    rows = []
    for obs in OBS:
        for E in ENERGIES:
            ib = nearest_bin(m, E)
            v, err = curve_err_fd(m, E, TH, ib, p, cov, obs)
            v0 = m0.curve(E, TH, nearest_bin(m0, E), r0.p, obs=obs)
            print('  %-4s E = %5.1f MeV -> bin %3g MeV: value %.3f .. %.3f, fit error %.4f .. %.4f, max |main(make_figures) - paper| = %.4f'
                  % (obs, E, m.bins[ib], v.min(), v.max(), err.min(), err.max(), np.abs(v0 - v).max()))
            if obs == 'S2X':
                assert abs(v[0]) < 1e-12 and abs(v[-1]) < 1e-12, 'Sigma_2x must vanish at 0 and 180 deg'
            rows += [dict(obs=obs, E=E, theta_lab=float(t), value=float(x), err=float(e)) for t, x, e in zip(TH, v, err)]
    out = pd.DataFrame(rows, columns=['obs', 'E', 'theta_lab', 'value', 'err'])
    out.to_csv(OUT, index=False, float_format='%.6g')
    print('written', OUT, '(%d rows)' % len(out))
    # bin-centre values for reference
    for obs in OBS:
        for E in ENERGIES:
            ib = nearest_bin(m, E)
            vb = m.curve(float(m.bins[ib]), np.array([0.0, 60.0, 90.0, 120.0, 150.0, 180.0]), ib, p, obs=obs)
            print('  %-4s at the bin centre %3g MeV, theta = 0, 60, 90, 120, 150, 180: %s' % (obs, m.bins[ib], np.round(vb, 4)))


if __name__ == '__main__':
    main()
