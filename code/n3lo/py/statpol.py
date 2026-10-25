"""statpol.py -- static dipole polarisabilities from A_1..A_6 (A_7 = A_8 = 0) via the L'vov amplitudes
(same conversion as work/code_v2/driver/static_pols.py: LMP2015 Eqs. A2, B3)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code_v2', 'driver'))
from static_pols import CAL, pols_from_a

def lvov(A6, nu, t):
    A = np.zeros(8, complex); A[:6] = A6
    return A @ CAL(nu, t)

def extract(Afun, fs=(0.5, 1.0, 2.0), nus=np.linspace(0.02, 0.12, 11), degs=(4, 5, 6)):
    """Afun(x, t) -> A_1..A_6; returns (alphaE1, betaM1, gE1E1, gM1M1, gE1M2, gM1E2) [1e-4 fm^3, fm^4] and spread"""
    rows = {}
    for f in fs:
        for nu in nus:
            t = -f * nu ** 2
            x = nu - t / 4            # nu = (s-u)/(4M) = x + t/4
            rows[(f, nu)] = lvov(Afun(x, t), nu, t).real
    a, sp = [], []
    for i in range(6):
        ests = []
        for f in fs:
            x2 = nus ** 2; y = np.array([rows[(f, nu)][i] for nu in nus])
            for d in degs: ests.append(np.polyfit(x2, y, d)[-1])
        a.append(np.median(ests)); sp.append(np.ptp(ests))
    p = pols_from_a(np.array(a)); dp = np.abs(pols_from_a(np.array(sp)))
    return p, dp

def g0pi(p):
    return -(p[2] + p[3] + p[4] + p[5]), p[3] + p[5] - p[2] - p[4]
