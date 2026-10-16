"""Step 4: estimate the correlation lengths of the truncation-error GP from the order-by-order coefficient
functions c2(E,theta), c3(E,theta) of the LEC-matched scheme (B) on a kinematic grid (Melendez et al. 2019):
maximise the GP marginal likelihood of the observed coefficients (zero mean, variance cbar^2, squared-exponential
kernel in E_lab and theta_cm).  Output: results/theory_lengths.json"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from common import *  # noqa
import theory as TH

rows = []
for E in np.arange(50, 151, 7.5):
    for th in np.arange(30, 171, 10):
        for o in ['CS', 'S3']:
            thc = float(database.th_lab2cm(E, th))
            rows.append(dict(pid='grid_%s_%g_%g' % (o, E, th), dataset='grid', E=E, theta=th, theta_lab=th, theta_cm=thc,
                             dobs=o, frame='LAB', value=0, err_pt=1, norm=0, flags=''))
g = pd.DataFrame(rows)
out = {}
for o in ['CS', 'S3']:
    s = g[g.dobs == o].reset_index(drop=True)
    y0, y2, y3 = order_values(g, s, tag='grid')
    tr = TH.truncation(y0, y2, y3, s['dobs'].values, s['E'].values)
    X = np.c_[s.E.values, s.theta_cm.values]
    res = {}
    for cn, k in [('c2', 0), ('c3', 0), ('c2', 3), ('c3', 3), ('c2', 4), ('c3', 4)]:
        # optional stationarisation: divide out the power-law growth (E/m_pi)^k of the matched coefficients
        c = tr[cn] / (np.minimum(s.E.values, 139.57) / 139.57) ** k

        def nll(p):
            lE, lt = np.exp(p)
            d2 = ((X[:, None, 0] - X[None, :, 0]) / lE) ** 2 + ((X[:, None, 1] - X[None, :, 1]) / lt) ** 2
            K = np.exp(-0.5 * d2) + 1e-6 * np.eye(len(c))
            L = np.linalg.cholesky(K)
            a = np.linalg.solve(L, c)
            n = len(c)
            v2 = a @ a / n            # profiled variance cbar^2
            return 0.5 * n * np.log(v2) + np.sum(np.log(np.diag(L)))
        best = None
        for lE0 in [20, 50, 100]:
            for lt0 in [20, 50, 100]:
                r = minimize(nll, np.log([lE0, lt0]), method='Nelder-Mead')
                if best is None or r.fun < best.fun:
                    best = r
        lE, lt = np.exp(best.x)
        res['%s_k%d' % (cn, k)] = dict(lE=float(lE), lth=float(lt), rms=float(np.sqrt(np.mean(c ** 2))))
        print('%s %s k=%d: l_E = %.0f MeV, l_theta = %.0f deg, rms = %.3f' % (o, cn, k, lE, lt, np.sqrt(np.mean(c ** 2))))
    out[o] = res
save_json(out, 'theory_lengths.json')
