#!/usr/bin/env python3
"""Compare Im of the l=1 Compton multipoles of ampl.exe (v2, with and without the gamma N Delta
vertex corrections) with the curves of the LMP2015 draft figure (red solid / red dashed) and MAID."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figdigit'))
import multfig
from pwa_run import run_pwa
M = 0.93827; MPI = 0.13957

def code_curves(oms, v2):
    Es = []
    for om in oms:
        W = om / 1000 + np.sqrt((om / 1000) ** 2 + M * M)
        Es.append(round((W * W - M * M) / (2 * M) * 1000, 6))
    r = run_pwa(Es, v2=v2)
    out = {}
    for om, E in zip(oms, Es):
        for k, v in r[float(E)].items():
            out.setdefault(k, []).append(v.imag * MPI * 1e3)
    return {k: np.array(v) for k, v in out.items()}

def main():
    oms = np.arange(140.0, 291.0, 10.0)
    D = multfig.digitise()
    C0 = code_curves(oms, {'dvertex': '.false.'})
    C1 = code_curves(oms, {'dvertex': '.true.'})
    print('Im f (1e-3/m_pi) at omega_cm = 200, 230, 260, 290 MeV;  [code v2 | LMP2015 draft curve | MAID]')
    for name in ['MM1+', 'ME1+', 'EE1+', 'EE1-', 'MM1-', 'EM1+']:
        for tag, C, kind in [('noVC', C0, 'novc'), ('VC  ', C1, 'vc')]:
            om_l, v_l = D[name][kind]
            om_m, v_m = D[name]['maid']
            s = '%-5s %s' % (name, tag)
            dev = []
            for om in [200, 230, 260, 290]:
                c = np.interp(om, oms, C[name]); l = np.interp(om, om_l, v_l); m = np.interp(om, om_m, v_m)
                s += ' | %7.4f %7.4f %7.4f' % (c, l, m)
            for om in oms[(oms >= 150) & (oms <= 290)]:
                dev.append(np.interp(om, oms, C[name]) - np.interp(om, om_l, v_l))
            scale = np.max(np.abs(v_l))
            s += '   max|code-LMP|/max|LMP| = %.3f' % (np.max(np.abs(dev)) / scale)
            print(s)
    return C0, C1, D

if __name__ == '__main__':
    main()
