"""Independent check of the kappa_v-linear Baldin integral with work/code_v2/tools/born_multipoles.py"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code_v2', 'tools'))
import born_multipoles as bm
M = 0.93827
def sigma(B, W, n=24):
    xs, ws = np.polynomial.legendre.leggauss(n); tot = 0.
    for ch in ('pi0p', 'pi+n'):
        for x, w in zip(xs, ws):
            K = B.kin(W, x)
            s2 = 0.
            for e4 in (np.array([0, 1.0, 0, 0]), np.array([0, 0, 1.0, 0])):
                T = B.matrix(K, ch, e4); s2 += np.sum(abs(T)**2)
            tot += w * 2*np.pi * s2 / 4 * K['q'] / K['k'] / (8*np.pi*W)**2
    return tot   # GeV^-2
def baldin(kp, kn, mpi):
    B = bm.calibrated(kp=kp, kn=kn, mpi=mpi)
    zs, wz = np.polynomial.legendre.leggauss(48); zs = 0.5*(zs+1); wz = 0.5*wz
    nuth = mpi + mpi**2/(2*M); I = 0.
    for z, w in zip(zs, wz):
        sc = 2*mpi; nu = nuth + sc*z/(1-z)**2; jac = sc*(1+z)/(1-z)**3
        W = np.sqrt(M*M + 2*M*nu)
        I += w*jac*sigma(B, W)/nu**2/(2*np.pi**2)
    return I * 0.19733**3 * 1e4      # GeV^-3 -> 1e-4 fm^3
for mpi in (0.13957, 0.03*M, 0.003*M):
    d = 0.5
    kv = (baldin(d/2, -d/2, mpi) - baldin(-d/2, d/2, mpi)) / (2*d)
    ks = (baldin(d/2, d/2, mpi) - baldin(-d/2, -d/2, mpi)) / (2*d)
    print(f'mpi/M={mpi/M:.4f}: Baldin, linear in kappa_v: {kv:.5f}  in kappa_s: {ks:.5f}  (1e-4 fm^3 per unit)')
