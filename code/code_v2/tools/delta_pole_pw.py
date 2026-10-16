#!/usr/bin/env python3
"""J=3/2+ projection of the code's s-channel Delta-pole Compton amplitude:
f_XY = N_XY(gM,gE)/D(W), D = s - MD^2 + i MD Gamma(W) (code's EDEM); returns F = -N (real)."""
import numpy as np
from pwa_run import run_pwa
SMM = 0.93827
HA = float(np.float32(2.85))
GA, FPI, MPI, MD = 1.27, 0.09221 / SMM, 0.13957 / SMM, 1.232 / SMM

def edem(W, Md=MD, hA=HA, fpi=FPI, mpi=MPI):
    if W < 1 + mpi:
        return 0.0
    Csq = (hA * Md / (8 * np.pi * fpi))**2
    what2 = ((W * W - 1 + mpi**2) / (2 * W))**2
    Eh = (W * W + 1 - mpi**2) / (2 * W)
    ph = np.sqrt(what2 - mpi**2)
    K0 = (Eh + 1) * ph**3 / W
    return 4 * np.pi / 3 * Csq * K0 / Md**2 * (W / Md)**2

def Dfun(W):
    s = W * W
    return s - MD**2 + 1j * MD * edem(W)        # M units

def F_quadratic(Es):
    """returns dict E -> {'MM': (a,b,c), 'EE':..., 'ME':...} with F = a gM^2 + b gM gE + c gE^2 (GeV^-1 * M^2-units D)"""
    res = {}
    runs = {}
    for (gM, gE) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]:
        runs[(gM, gE)] = run_pwa(list(Es), gM=gM, gE=gE, oneN=False, oneD=False, anomaly=False,
                                 nucl=(False, False), delta=(True, False), v2={'dvertex': '.false.'})
    for E in Es:
        W = np.sqrt(1 + 2 * E / 1000 / SMM)
        D = Dfun(W)
        out = {}
        for key, lab in [('MM', 'MM1+'), ('EE', 'EE2-'), ('ME', 'ME1+')]:
            v = {gg: -(runs[gg][float(E)][lab] * D) for gg in runs}
            a = v[(1.0, 0.0)]; c = v[(0.0, 1.0)]; b = v[(1.0, 1.0)] - a - c
            out[key] = (a, b, c)
        res[E] = out
    return res

if __name__ == '__main__':
    Es = [150, 200, 250, 300, 340]
    R = F_quadratic(Es)
    for E in Es:
        W = np.sqrt(1 + 2 * E / 1000 / SMM)
        MM, EE, ME = R[E]['MM'], R[E]['EE'], R[E]['ME']
        print('E=%d W=%.4f GeV' % (E, W * SMM))
        for k, v in R[E].items():
            print('   F_%s = %+.5e gM^2 %+.5e gMgE %+.5e gE^2   (max Im %.1e)' % (k, v[0].real, v[1].real, v[2].real, max(abs(x.imag) for x in v)))
        a = MM[0].real
        kap = -MM[1].real / (2 * a)
        print('   rank-1 checks: F_MM c/a vs kappa^2 : %.6f %.6f ; F_EE a,b = %.1e %.1e ; F_ME a = %.1e' %
              (MM[2].real / a, kap**2, EE[0].real, EE[1].real, ME[0].real))
        print('   kappa(W) = %.6f ;  (W-M)/(2(W+M)) = %.6f ;  (MD-M)/(2(MD+M)) = %.6f' %
              (kap, (W - 1) / (2 * (W + 1)), (MD - 1) / (2 * (MD + 1))))
        print('   F_ME^2 / (F_MM F_EE) at gM=2.97,gE=-1: %.6f' % (
            (ME[1].real * 2.97 * -1 + ME[2].real) ** 2 / ((a * 2.97**2 + MM[1].real * 2.97 * -1 + MM[2].real) * EE[2].real)))
