#!/usr/bin/env python3
"""Unitarity check of the piN loops with physical pion masses (&v2 pimassphys=t), separately for the
pi0 p and pi+ n channels (&v2 pichan=2/3) and for their sum:
   Im f_code = k_c sum |M_c|^2 ,  M_c = LO Born multipoles of gamma p -> pi0 p (m_pi0) / pi+ n (m_pi+)."""
import numpy as np
from pwa_run import run_pwa
from born_multipoles import calibrated
M = 0.93827; MC = 0.13957; MN = 0.13498
KEYS = ['MM1+', 'EE2-', 'ME1+', 'EE1-', 'MM1-', 'EE1+', 'EM1+']

def pred(B, W, m, ch):
    s = W * W
    Eq = (s + m * m - M * M) / (2 * W)
    if Eq <= m:
        return {k: 0.0 for k in KEYS}
    k = np.sqrt(Eq * Eq - m * m)
    a = B.multipoles(W, ch, 2)
    return {'MM1+': k * abs(a['M1+'])**2, 'EE2-': k * abs(a['E1+'])**2, 'ME1+': -k * (a['E1+'] * np.conj(a['M1+'])).real,
            'EE1-': k * abs(a['E0+'])**2, 'MM1-': k * abs(a['M1-'])**2, 'EE1+': k * abs(a['E2-'])**2,
            'EM1+': k * (a['E2-'] * np.conj(a['M2-'])).real}

def main(Es=(147.0, 150.0, 155.0, 200.0, 300.0)):
    B0 = calibrated(kp=0.0, kn=0.0, mpi=MN)
    Bc = calibrated(kp=0.0, kn=0.0, mpi=MC)
    worst = 0.0
    for chan, lab in [(2, 'pi0 p'), (3, 'pi+ n'), (0, 'sum')]:
        r = run_pwa(list(Es), oneD=False, anomaly=False, delta=(False, False),
                    v2={'pimassphys': '.true.', 'pichan': chan})
        for E in Es:
            W = np.sqrt(M * M + 2 * M * E / 1000)
            p0 = pred(B0, W, MN, 'pi0p'); pc = pred(Bc, W, MC, 'pi+n')
            P = {k: (p0[k] if chan in (0, 2) else 0) + (pc[k] if chan in (0, 3) else 0) for k in KEYS}
            d = r[float(E)]
            line = '%-6s E=%5.1f' % (lab, E)
            for k in KEYS:
                ref = P[k]; got = d[k].imag
                if abs(ref) < 1e-14:
                    rel = abs(got) / 1e-6     # should vanish below threshold (report in units of 1e-6 GeV^-1)
                    line += ' %s %8.1e(abs)' % (k, abs(got))
                else:
                    rel = got / ref - 1
                    worst = max(worst, abs(rel))
                    line += ' %s %+.1e' % (k, rel)
            print(line)
    print('channel unitarity: max |Im f/(k|M_Born|^2) - 1| = %.2e' % worst)
    return worst

if __name__ == '__main__':
    main()
