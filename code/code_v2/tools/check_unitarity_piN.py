#!/usr/bin/env python3
"""Unitarity check of the code's piN-loop Compton multipoles: Im f = k sum_c |M_c|^2 with the
tree-level (LO, kappa=0) photoproduction multipoles of born_multipoles.py (channels pi0 p, pi+ n,
isospin-symmetric pion mass), LMP2015 Eqs. (20)-(21)."""
import sys, numpy as np
from pwa_run import run_pwa
from born_multipoles import calibrated

def main(Es=(160, 200, 250, 300, 350), verbose=True):
    B = calibrated(kp=0.0, kn=0.0)
    r = run_pwa(list(Es), oneD=False, anomaly=False, delta=(False, False))
    M, m = 0.93827, 0.13957
    worst = 0.0
    for E in Es:
        s = M * M + 2 * M * E / 1000; W = np.sqrt(s)
        Eq = (s + m * m - M * M) / (2 * W); k = np.sqrt(Eq * Eq - m * m)
        a = B.multipoles(W, 'pi0p', 2); b = B.multipoles(W, 'pi+n', 2)
        S = lambda f: sum(f(x) for x in (a, b))
        pred = {'MM1+': k * S(lambda x: abs(x['M1+'])**2), 'EE2-': k * S(lambda x: abs(x['E1+'])**2),
                'ME1+': -k * S(lambda x: (x['E1+'] * np.conj(x['M1+'])).real),
                'EE1-': k * S(lambda x: abs(x['E0+'])**2), 'MM1-': k * S(lambda x: abs(x['M1-'])**2),
                'EE1+': k * S(lambda x: abs(x['E2-'])**2),
                'EM1+': k * S(lambda x: (x['E2-'] * np.conj(x['M2-'])).real),
                'MM2-': k * S(lambda x: abs(x['M2-'])**2)}
        d = r[float(E)]
        line = 'E=%3d' % E
        for key in ['MM1+', 'EE2-', 'ME1+', 'EE1-', 'MM1-', 'EE1+', 'EM1+', 'MM2-']:
            rel = d[key].imag / pred[key] - 1
            worst = max(worst, abs(rel))
            line += ' %s %+.1e' % (key, rel)
        if verbose:
            print(line)
    print('piN-loop unitarity: max |Im f_code/(k sum|M_Born|^2) - 1| = %.2e' % worst)
    return worst

if __name__ == '__main__':
    main()
