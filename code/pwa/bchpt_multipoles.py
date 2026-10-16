"""Non-Born BChPT multipoles (NNLO, Lensky-McGovern-Pascalutsa 2015) from the Fortran code:
project the non-Born helicity amplitudes (all mechanisms except Born and pi0 pole) on partial waves
with an endpoint-avoiding 20-point Gauss rule (see validation), on a grid of lab energies.
Output: results/bchpt_multipoles.csv  (E_lab [MeV], multipoles in fm, real and imaginary parts)"""
import numpy as np
import pandas as pd
import compton as C
import fortran_runner as FR

KEYS = [('EE', 1, '-'), ('EE', 1, '+'), ('MM', 1, '-'), ('MM', 1, '+'), ('EM', 1, '+'), ('ME', 1, '+'),
        ('EE', 2, '-'), ('EE', 2, '+'), ('MM', 2, '-'), ('MM', 2, '+'), ('EM', 2, '+'), ('ME', 2, '+'),
        ('EE', 3, '-'), ('MM', 3, '-')]


def keyname(k):
    return f'f{k[0]}{k[1]}{"p" if k[2] == "+" else "m"}'


def compute(Egrid=None, pols=None):
    if Egrid is None:
        Egrid = np.arange(20.0, 170.01, 1.0)
    xg, wg = np.polynomial.legendre.leggauss(20)
    thg = np.degrees(np.arccos(xg))
    EE, TT = np.meshgrid(Egrid, thg, indexing='ij')
    M, _ = FR.run(EE.ravel(), TT.ravel(), 'CM', 'bchpt', pols)
    M = M.reshape(len(Egrid), len(thg), *M.shape[1:])
    rows = []
    for i, E in enumerate(Egrid):
        P = C.phis_from_M(M[i, :, 2:].sum(axis=1))
        f = C.multipoles_from_pw(C.pw_project(P, xg, wg, 9))
        r = {'E_lab': E}
        for k in KEYS:
            r[keyname(k) + '_re'] = np.real(f[k]); r[keyname(k) + '_im'] = np.imag(f[k])
        rows.append(r)
    return pd.DataFrame(rows)


class BChPT:
    """interpolator of BChPT non-Born multipoles (real parts) in fm."""
    def __init__(self, path='results/bchpt_multipoles.csv'):
        self.df = pd.read_csv(path)

    def __call__(self, key, E):
        return np.interp(E, self.df.E_lab.values, self.df[keyname(key) + '_re'].values)


if __name__ == '__main__':
    D = compute()
    D.to_csv('results/bchpt_multipoles.csv', index=False)
    print(D.iloc[::20, :8].to_string())
