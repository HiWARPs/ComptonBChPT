"""Total cross section of proton Compton scattering, sigma(gamma p -> gamma p), from the v2 BChPT amplitude,
integrating dsigma/dOmega_cm over the c.m. angle (32-point Gauss-Legendre), and its contribution to the Baldin
integral (1/2pi^2) int dnu sigma/nu^2.

Main result ('NNLO'): the parameter-free NNLO prediction, P = Pref of the driver run (= common.PREF_NNLO of
work/fits/lib: the code's physical-mass reference values, alpha = 11.415, beta = 3.729, spin pols -3.246, 2.939,
0.170, 1.133), with the band 'NNLO_lec' of the O(p^4) LECs of natural size, delta alpha_E1 = delta beta_M1 = +-DLEC
= +-0.7 (Lensky-Pascalutsa 2010, Eq. 25): central differences of the integrated cross section at Pref +- 0.7 e_k,
added in quadrature (exact: predict_array is the exact quadratic form of the amplitude in the LECs).
Also: the NNLO fit ('NNLO_fit': alpha, beta from work/fits/results/fits_L.json['F2_d_all_thB'], spin pols at NNLO)
with the band 'NNLO_fit_err' from the (alpha, beta) fit covariance; Born + pi0 pole ('BornPi0'); the prediction
without the vertex corrections ('NNLO_noVC'); Thomson.  Until 2026-10-16 the fit was the main result
(work/backups/sigma_tot_before_predmain.*).
Run:  python3 sigma_tot.py   (any cwd; writes sigma_tot.json here)
then  plot_sigma_tot.py  ->  paper/figs/fig_sigma_tot.pdf"""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import os, sys, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, (_WORK + '/fits/lib'))
os.environ.setdefault('FIT_PROFILE', 'v2')
import bchpt
DLEC = 0.7
# NNLO fit (F2 + truncation covariance, 3-sigma outlier removed): alpha, beta and their covariance are read from
# work/fits/results/fits_L.json ('F2_d_all_thB'); the values below are only the fallback if the file is missing
_FL = (_WORK + '/fits/results/fits_L.json')
try:
    _m = json.load(open(_FL))['F2_d_all_thB']
    A_FIT, B_FIT = float(_m['P'][0]), float(_m['P'][1])
    COV_AB = np.array(_m['covP'])[:2, :2]
    print('NNLO fit from fits_L.json: alpha = %.4f, beta = %.4f, rho = %.3f' % (A_FIT, B_FIT, COV_AB[0, 1] / np.sqrt(COV_AB[0, 0] * COV_AB[1, 1])))
except (OSError, KeyError):
    A_FIT, B_FIT = 11.672, 2.876
    COV_AB = np.array([[0.4668**2, -0.84*0.4668*0.4319], [-0.84*0.4668*0.4319, 0.4319**2]])
    print('WARNING: fits_L.json not found, using fallback alpha, beta')
E = np.concatenate([np.arange(5, 140, 5.0), np.arange(140, 160, 1.0), np.arange(160, 451, 5.0)])
x, w = np.polynomial.legendre.leggauss(32)            # cos(theta_cm)
th = np.degrees(np.arccos(x))
pts = [(f'p{i}_{j}', float(e), float(t), 'CSCM', 'CM') for i, e in enumerate(E) for j, t in enumerate(th)]


def integrate(out, P):
    """sigma_tot [mu b] on the grid E from the driver output at polarisabilities P"""
    vals = out.predict_array(P)                 # nb/sr
    sig = []
    for i, e in enumerate(E):
        v = np.array([vals[out.index[f'p{i}_{j}']] for j in range(len(th))])
        sig.append(2 * np.pi * np.sum(w * v) / 1000.0)   # mu b
    return np.array(sig)


res, info = {}, {}
for cfg in ('NNLO', 'BornPi0', 'NNLO_noVC'):
    c = cfg if cfg != 'NNLO_noVC' else dict(dvertex=False)
    out = bchpt.run(pts, c, name='sigtot_' + (cfg if isinstance(cfg, str) else 'x'), verbose=False)
    Pref = np.array(list(out.Pref))
    res[cfg] = integrate(out, Pref)              # prediction (BornPi0: no polarisabilities)
    if cfg == 'NNLO':
        info['Pref'] = Pref.tolist()
        # LEC band of the prediction: central differences at Pref +- DLEC e_k, in quadrature
        grads_lec = []
        for k in range(2):
            Pp, Pm = Pref.copy(), Pref.copy()
            Pp[k] += DLEC; Pm[k] -= DLEC
            grads_lec.append((integrate(out, Pp) - integrate(out, Pm)) / 2)
        res['NNLO_lec_da'], res['NNLO_lec_db'] = grads_lec
        res['NNLO_lec'] = np.sqrt(grads_lec[0]**2 + grads_lec[1]**2)
        # NNLO fit: alpha, beta of the fit, spin pols at NNLO; band from the (alpha, beta) fit covariance
        Pfit = Pref.copy(); Pfit[0], Pfit[1] = A_FIT, B_FIT
        res['NNLO_fit'] = integrate(out, Pfit)
        grads = []
        for k in range(2):
            Pp, Pm = Pfit.copy(), Pfit.copy()
            Pp[k] += 0.05; Pm[k] -= 0.05
            grads.append((integrate(out, Pp) - integrate(out, Pm)) / 0.1)
        G = np.vstack(grads).T
        res['NNLO_fit_err'] = np.sqrt(np.einsum('ij,jk,ik->i', G, COV_AB, G))
r2 = (1/137.036 * 197.327 / 938.272)**2 * 1e4   # Thomson r^2 in mu b (fm^2 -> mu b: 1 fm^2 = 1e4 mu b)
res['Thomson'] = np.full_like(E, 8*np.pi/3*r2)


# Baldin-type integral: (1/2pi^2) int dnu sigma/nu^2 ; sigma [mu b] -> fm^2 *1e-4 ; nu [MeV] -> /197.327 fm^-1
def baldin(sig, Emin, Emax):
    m = (E >= Emin) & (E <= Emax)
    nu = E[m] / 197.327                          # fm^-1
    s = sig[m] * 1e-4                            # fm^2
    return np.trapezoid(s / nu**2, nu) / (2 * np.pi**2) * 1e4   # in 1e-4 fm^3


RANGES = [(5, 150), (50, 150), (150, 450), (150, 350), (200, 450), (50, 450)]
summary = {}
for cfg in ('NNLO', 'NNLO_fit', 'BornPi0', 'NNLO_noVC'):
    summary[cfg] = {f'{a}-{b}': round(baldin(res[cfg], a, b), 4) for a, b in RANGES}
# LEC uncertainty of the integrals of the prediction (the integral is linear in sigma: same central differences)
summary['NNLO_lec'] = {f'{a}-{b}': round(float(np.hypot(baldin(res['NNLO_lec_da'], a, b), baldin(res['NNLO_lec_db'], a, b))), 4)
                       for a, b in RANGES}
summary['NNLO_minus_BornPi0'] = {f'{a}-{b}': round(baldin(res['NNLO'] - res['BornPi0'], a, b), 4) for a, b in [(5, 150), (150, 450), (5, 450)]}
summary['NNLO_fit_minus_BornPi0'] = {f'{a}-{b}': round(baldin(res['NNLO_fit'] - res['BornPi0'], a, b), 4) for a, b in [(5, 150), (150, 450), (5, 450)]}
for key in ('NNLO', 'NNLO_fit'):
    i = int(np.argmax(res[key]))
    summary['peak_' + key] = dict(E=float(E[i]), sigma_mub=float(res[key][i]),
                                  **({'lec': float(res['NNLO_lec'][i])} if key == 'NNLO' else {'fit_err': float(res['NNLO_fit_err'][i])}))
summary['prediction'] = dict(Pref=info['Pref'], dlec=DLEC, band='central differences at Pref +- dlec e_k (alpha_E1, beta_M1), in quadrature')
summary['nnlo_fit'] = dict(alphaE1=A_FIT, betaM1=B_FIT, cov=COV_AB.tolist(), source=_FL + " ['F2_d_all_thB']")
json.dump(dict(E=E.tolist(), **{k: v.tolist() for k, v in res.items()}, baldin=summary),
          open(os.path.join(HERE, 'sigma_tot.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
for e in (50, 100, 140, 150, 200, 250, 300, 320, 350, 400, 450):
    i = int(np.argmin(abs(E - e)))
    print(f"E={E[i]:5.0f}  NNLO pred={res['NNLO'][i]:.3f}+-{res['NNLO_lec'][i]:.3f}  fit={res['NNLO_fit'][i]:.3f}+-{res['NNLO_fit_err'][i]:.3f}"
          f"  BornPi0={res['BornPi0'][i]:.3f}  noVC={res['NNLO_noVC'][i]:.3f}  Thomson={res['Thomson'][i]:.3f} mub")
