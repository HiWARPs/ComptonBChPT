"""Bare curves of the PWA fit of the paper (scenarios/main_fits_legacynorm.py, fit 'MAIN') on the angular grids
of the fit figures of work/fits/scripts/s06_figures.py (paper Figs. cs, sigma3, higs).  Writes
results/pwa_curves.csv with the columns obs, E, theta_lab, value:
    obs = CS     unpolarised lab cross section dsigma/dOmega_lab [nb/sr]
          S3     beam asymmetry Sigma_3 (lab)
          SPAR   dsigma_parallel/dOmega_lab [nb/sr]   (linearly polarised beam, photon polarisation in the scattering plane)
          SPERP  dsigma_perp/dOmega_lab    [nb/sr]
i.e. the names of the BChPT fit database (work/fits/db/fitdb.csv: CS, S3, SPAR, SPERP = PWA database CS, S3, CS_PAR, CS_PERP).

The values are the THEORY curve of the PWA fit WITHOUT any normalisation factor (like the NNLO curves of the figures,
which are compared with the data without the fitted normalisations): the fitted normalisation N_k of a data set
multiplies the theory only inside the chi^2, model = N_k * theory (pwa_model.Model.residuals), so at the data points
    value(E, theta) * N_k  ==  'pred' of results/update/pulls_L_MAIN.csv
which is verified below for every figure energy (A2 cross sections and Sigma_3, HIGS circular and linear).

Energy within a bin.  The PWA parametrisation (KLP18, pwa_model.Model.multipoles_at) is the Born + pi0-pole amplitude
plus the l = 1 non-Born multipoles
    f_i(E) = pre_i(E) * [ static polarisabilities (alpha, beta, gamma's; E-dependent LO/NLO terms) + (E/M)^2 f^R_i ]
with the l = 2 multipoles fixed to NNLO BChPT.  Everything except the residual values f^R_i depends continuously on the
actual photon energy E (the kinematic prefactors x^n/sqrt(s), x = E/hbar c, the static-polarisability terms, the Born
amplitudes, the BChPT l = 2 multipoles, and the (E/M)^2 factor multiplying f^R_i).  The six f^R_i (f^R_2, f^R_3 are
free per bin in the chi^2 sum-rule mode of the main fit) are constants per energy bin; the main fit uses
energy_mode='actual' (update_fits.build_model), i.e. every data point is evaluated at its own energy with the f^R_i of
the bin whose centre is nearest to it (bins BINS_A = 59, 69, 79, 89, 99, 109, 117, 127, 135, 143, 150 MeV).  The curves
here are evaluated in exactly the same way: at the panel energy E, with the residuals of the bin nearest to E
(Model.curve(E, theta, ib, p); the same function as in make_figures.py).  Hence at a data energy the curve is the
PWA description of these data up to the normalisation; between the data energies the curve is a smooth interpolation
within a bin and may jump by the change of f^R_i at a bin boundary.  The panel-energy -> bin assignment is printed.
149.87 MeV lies in the highest bin (150 MeV, which contains the TAPS and SAL 93 data up to 150 MeV), so it is covered.
The fig1b energies of s06 (58.94, 79.19, 108.42 MeV, older data) are added for the CS curves (bins 59, 79, 109).

Run from work/pwa (or anywhere: the script changes to work/pwa):  .venv/bin/python scenarios/pwa_curves.py   (~2 s).
The fit parameters are read from results/update/p_L_MAIN.npy (written by scenarios/main_fits_legacynorm.py); the chi^2
of the data at these parameters must reproduce results/main_fits_legacynorm.json['MAIN']['chi2'] (checked)."""
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
import main_fits_legacynorm as MFL            # noqa: E402  (the scenario: load_leg(), MAIN; main() not run on import)

OUT = os.path.join(PWA, 'results', 'pwa_curves.csv')
TH_FULL = np.arange(0.0, 175.1, 2.5)
TH_S3 = np.arange(20.0, 175.1, 2.5)
# (obs of the fit database, PWA observable key of Model.curve, energies, angular grid)
GRIDS = [('CS', 'cs_lab', [81.3, 92.2, 103.2, 113.6, 124.1, 135.3, 149.87, 58.94, 79.19, 108.42], TH_FULL),
         ('S3', 'S3', [83.4, 92.2, 108.5, 129.5], TH_S3),
         ('SPAR', 'cs_par_lab', [83.4], TH_FULL),
         ('SPERP', 'cs_perp_lab', [83.4], TH_FULL)]
# PWA database (dataset, obs) -> figure obs, for the check against the fitted values
CHECK = [('A2_2022_CS', 'CS', 'cs_lab'), ('A2_2022_S3', 'S3', 'S3'), ('HIGS_2022_circ', 'CS', 'cs_lab'),
         ('HIGS_2022_lin', 'CS_PAR', 'cs_par_lab'), ('HIGS_2022_lin', 'CS_PERP', 'cs_perp_lab')]


def nearest_bin(m, E):
    return int(np.argmin(abs(m.bins - E)))


def main():
    d = MFL.load_leg()
    m = UF.build_model(UF.select(d), UF.BINS_A, **MFL.MAIN)
    p = np.load('results/update/p_L_MAIN.npy')
    cov = np.load('results/update/cov_L_MAIN.npy')
    J = json.load(open('results/main_fits_legacynorm.json'))['MAIN']
    assert p.shape == (m.npar,) and cov.shape == (m.npar, m.npar), 'parameter layout differs from the saved fit'
    chi2 = m.chi2_data(p)
    print('PWA main fit: N = %d, chi2 = %.4f (json %.4f), npar = %d' % (int(m.mask.sum()), chi2, J['chi2'], m.npar))
    assert abs(chi2 - J['chi2']) < 1e-6, 'saved parameters do not reproduce the chi2 of the main fit'
    npm = m.ng + m.nr * m.nb
    norms = {g: 1.0 + m.gsig[j] * p[npm + j] for j, g in enumerate(m.ngroups)}
    for g, N in norms.items():
        assert abs(N - J['norms'][g]['N']) < 1e-9
    # ---- check: bare curve x N_k == fitted value 'pred' at every data point of the figure energies ----
    pulls = pd.read_csv('results/update/pulls_L_MAIN.csv')
    worst = 0.0
    for ds, obs, key in CHECK:
        s = pulls[(pulls.dataset == ds) & (pulls.obs == obs)]
        if len(s) == 0:        # HIGS 2022 is not part of the main database (2026-10-17): the curve is a PWA prediction (N = 1)
            print('  %-15s %-7s not in the main fit (HIgS excluded since 2026-10-17): curve = prediction of the PWA fit, no check' % (ds, obs))
            continue
        for E, sE in s.groupby('E'):
            ib = nearest_bin(m, E)
            assert np.all(sE.ibin.values == ib), 'bin assignment differs from the fit'
            v = m.curve(float(E), sE.theta_lab.values, ib, p, obs=key)
            N = norms.get(ds, 1.0)
            dev = np.abs(v * N - sE.pred.values) / np.abs(sE.pred.values)
            worst = max(worst, dev.max())
            print('  %-15s %-7s E = %6.2f MeV (bin %3g, N = %.4f): %2d pts, max |curve*N/pred - 1| = %.1e' % (
                ds, obs, E, m.bins[ib], N, len(sE), dev.max()))
    assert worst < 1e-8, 'bare curve times normalisation does not reproduce the fitted values'
    assert not pulls.dataset.isin(UF.HIGS).any(), 'HIGS data found in the main fit (excluded since 2026-10-17)'
    print('check passed: curve * N_k reproduces the fitted values at all data points of the main fit (max rel. dev. %.1e)' % worst)
    # ---- curves ----
    rows = []
    for obs, key, Es, th in GRIDS:
        for E in Es:
            ib = nearest_bin(m, E)
            v = m.curve(float(E), th, ib, p, obs=key)
            _, err = m.curve_err(float(E), th, ib, p, cov, obs=key)
            print('  %-5s E = %6.2f MeV -> bin %3g MeV (residuals f^R_i of this bin); value range %.3f .. %.3f, '
                  'fit error (stat.) %.3f .. %.3f' % (obs, E, m.bins[ib], v.min(), v.max(), err.min(), err.max()))
            rows += [dict(obs=obs, E=E, theta_lab=float(t), value=float(x)) for t, x in zip(th, v)]
    out = pd.DataFrame(rows, columns=['obs', 'E', 'theta_lab', 'value'])
    out.to_csv(OUT, index=False, float_format='%.6g')
    print('written', OUT, '(%d rows)' % len(out))


if __name__ == '__main__':
    main()
