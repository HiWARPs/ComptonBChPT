"""Total cross section of proton Compton scattering, sigma(gamma p -> gamma p), from the v2 BChPT amplitude
(main fit alpha, beta; spin pols NNLO), integrating dsigma/dOmega_cm over the c.m. angle (Gauss-Legendre),
and its contribution to the Baldin integral (1/2pi^2) int dnu sigma/nu^2.  Writes sigma_tot.json/.csv."""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import os, sys, json, numpy as np
sys.path.insert(0, (_WORK + '/fits/lib'))
os.environ.setdefault('FIT_PROFILE', 'v2')
import bchpt
# main BChPT fit (F2 + truncation covariance, 3-sigma outlier removed): alpha, beta and their covariance are read
# from work/fits/results/fits_L.json ('F2_d_all_thB'); the values below are only the fallback if the file is missing
_FL = (_WORK + '/fits/results/fits_L.json')
try:
    _m = json.load(open(_FL))['F2_d_all_thB']
    A_MAIN, B_MAIN = float(_m['P'][0]), float(_m['P'][1])
    COV_AB = np.array(_m['covP'])[:2, :2]
    print('main fit from fits_L.json: alpha = %.4f, beta = %.4f, rho = %.3f' % (A_MAIN, B_MAIN, COV_AB[0, 1] / np.sqrt(COV_AB[0, 0] * COV_AB[1, 1])))
except (OSError, KeyError):
    A_MAIN, B_MAIN = 11.672, 2.876
    COV_AB = np.array([[0.4668**2, -0.84*0.4668*0.4319], [-0.84*0.4668*0.4319, 0.4319**2]])
    print('WARNING: fits_L.json not found, using fallback alpha, beta')
E = np.concatenate([np.arange(5, 140, 5.0), np.arange(140, 160, 1.0), np.arange(160, 451, 5.0)])
x, w = np.polynomial.legendre.leggauss(32)            # cos(theta_cm)
th = np.degrees(np.arccos(x))
pts = [(f'p{i}_{j}', float(e), float(t), 'CSCM', 'CM') for i, e in enumerate(E) for j, t in enumerate(th)]
res = {}
P_main = None
for cfg in ('NNLO', 'BornPi0', 'NNLO_noVC'):
    c = cfg if cfg != 'NNLO_noVC' else dict(dvertex=False)
    out = bchpt.run(pts, c, name='sigtot_' + (cfg if isinstance(cfg, str) else 'x'), verbose=False)
    P = np.array(list(out.Pref))
    if cfg != 'BornPi0':
        P[0], P[1] = A_MAIN, B_MAIN             # main fit (v2, from fits_L.json)
    vals = out.predict_array(P)                 # nb/sr
    sig = []
    for i, e in enumerate(E):
        v = np.array([vals[out.index[f'p{i}_{j}']] for j in range(len(th))])
        sig.append(2 * np.pi * np.sum(w * v) / 1000.0)   # mu b
    res[cfg] = np.array(sig)
    if cfg == 'NNLO':
        # alpha, beta fit band: derivative wrt alpha, beta via exact quadratic forms
        cov = COV_AB
        grads = []
        for k in range(2):
            Pp = P.copy(); Pp[k] += 0.05; Pm = P.copy(); Pm[k] -= 0.05
            vp, vm = out.predict_array(Pp), out.predict_array(Pm)
            g = []
            for i, e in enumerate(E):
                a = np.array([vp[out.index[f'p{i}_{j}']] for j in range(len(th))]); b = np.array([vm[out.index[f'p{i}_{j}']] for j in range(len(th))])
                g.append(2 * np.pi * np.sum(w * (a - b) / 0.1) / 1000.0)
            grads.append(np.array(g))
        G = np.vstack(grads).T
        res['NNLO_err'] = np.sqrt(np.einsum('ij,jk,ik->i', G, cov, G))
r2 = (1/137.036 * 197.327 / 938.272)**2 * 1e4   # Thomson r^2 in mu b (fm^2 -> mu b: 1 fm^2 = 1e4 mu b)
res['Thomson'] = np.full_like(E, 8*np.pi/3*r2)
# Baldin-type integral: (1/2pi^2) int dnu sigma/nu^2 ; sigma [mu b] -> fm^2 *1e-4 ; nu [MeV] -> /197.327 fm^-1
def baldin(sig, Emin, Emax):
    m = (E >= Emin) & (E <= Emax)
    nu = E[m] / 197.327                          # fm^-1
    s = sig[m] * 1e-4                            # fm^2
    return np.trapezoid(s / nu**2, nu) / (2 * np.pi**2) * 1e4   # in 1e-4 fm^3
summary = {}
for cfg in ('NNLO', 'BornPi0', 'NNLO_noVC'):
    summary[cfg] = {f'{a}-{b}': round(baldin(res[cfg], a, b), 4) for a, b in [(5, 150), (50, 150), (150, 450), (150, 350), (200, 450), (50, 450)]}
    sig_nb = res[cfg] - res['BornPi0'] if cfg != 'BornPi0' else None
summary['NNLO_minus_BornPi0'] = {f'{a}-{b}': round(baldin(res['NNLO'] - res['BornPi0'], a, b), 4) for a, b in [(5, 150), (150, 450), (5, 450)]}
summary['peak'] = dict(E=float(E[np.argmax(res['NNLO'])]), sigma_mub=float(res['NNLO'].max()))
summary['main_fit'] = dict(alphaE1=A_MAIN, betaM1=B_MAIN, cov=COV_AB.tolist(), source=_FL + " ['F2_d_all_thB']")
json.dump(dict(E=E.tolist(), **{k: v.tolist() for k, v in res.items()}, baldin=summary), open('sigma_tot.json', 'w'), indent=1)
print(json.dumps(summary, indent=1))
for e in (50, 100, 140, 150, 200, 250, 300, 320, 350, 400, 450):
    i = int(np.argmin(abs(E - e))); print(f"E={E[i]:5.0f}  NNLO={res['NNLO'][i]:.3f}+-{res['NNLO_err'][i]:.3f}  BornPi0={res['BornPi0'][i]:.3f}  noVC={res['NNLO_noVC'][i]:.3f}  Thomson={res['Thomson'][i]:.3f} mub")
