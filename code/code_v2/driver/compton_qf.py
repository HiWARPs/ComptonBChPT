#!/usr/bin/env python3
"""
compton_qf.py -- read the fit-driver output of ampl.exe (run_points.sh) and evaluate / fit observables.

Every observable O at a kinematic point is an exact ratio of two quadratic forms in the
shifts d_k = P_k - Pref_k of the six static dipole polarisabilities

    P = (alphaE1, betaM1 [1e-4 fm^3];  gE1E1, gM1M1, gE1M2, gM1E2 [1e-4 fm^4])

    O(P) = N(d) / D(d),   N(d) = sum_{a<=b} cN_ab d_a d_b ,  d_0 = 1   (28 coefficients, same for D)

For cross sections D == 1.  The coefficients are computed once per point by ampl.exe; afterwards
any parameter set is evaluated with two (Npoints x 28) matrix-vector products, i.e. microseconds.

Usage as a module:
    from compton_qf import DriverOutput
    T = DriverOutput('test_points.out')
    T.predict(alphaE1=11.0, betaM1=3.0)          # unspecified parameters stay at Pref (NNLO)
    T.predict_array(P)                            # P: (6,) or (n,6) absolute polarisabilities

Command line:
    python3 compton_qf.py test_points.out                      # print values at NNLO and a shifted set
    python3 compton_qf.py test_points.out --demo-fit           # closure test of a chi^2 fit on pseudo-data
"""
import sys
import numpy as np

PNAMES = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
PAIRS = [(a, b) for a in range(7) for b in range(a, 7)]


class DriverOutput:
    def __init__(self, fn):
        self.fn = fn
        self.Pref = None
        self.Puser = None
        ids, obs, E, thl, thc, val, flag, cn, cd = [], [], [], [], [], [], [], [], []
        for line in open(fn):
            if line.startswith('#'):
                if line.startswith('# Pref'):
                    self.Pref = np.array([float(v) for v in line.split('=')[1].split()])
                elif line.startswith('# Puser'):
                    self.Puser = np.array([float(v) for v in line.split('=')[1].split()])
                continue
            p = line.split()
            if not p:
                continue
            ids.append(p[0]); obs.append(p[1])
            E.append(float(p[2])); thl.append(float(p[3])); thc.append(float(p[4]))
            val.append(float(p[5])); flag.append(int(p[6]))
            c = np.array([float(v) for v in p[7:]])
            if c.size != 56:
                raise ValueError('bad coefficient count in line: ' + p[0])
            cn.append(c[:28]); cd.append(c[28:])
        self.ids = np.array(ids); self.obs = np.array(obs)
        self.E = np.array(E); self.theta_lab = np.array(thl); self.theta_cm = np.array(thc)
        self.value_file = np.array(val); self.flag = np.array(flag)
        self.CN = np.array(cn); self.CD = np.array(cd)
        if self.Pref is None:
            raise ValueError('no "# Pref" header line found in ' + fn)

    @staticmethod
    def monomials(d):
        """d: (...,6) shifts -> (...,28) monomials d_a d_b (a<=b) with d_0 = 1."""
        d = np.atleast_2d(d)
        dd = np.concatenate([np.ones((d.shape[0], 1)), d], axis=1)
        return np.stack([dd[:, a] * dd[:, b] for a, b in PAIRS], axis=1)

    def predict_array(self, P):
        """P: absolute polarisabilities, shape (6,) or (n,6). Returns (npoints,) or (n, npoints)."""
        P = np.asarray(P, dtype=float)
        single = (P.ndim == 1)
        m = self.monomials(np.atleast_2d(P) - self.Pref)          # (n,28)
        out = (m @ self.CN.T) / (m @ self.CD.T)                   # (n,npoints)
        return out[0] if single else out

    def predict(self, **kw):
        P = self.Pref.copy()
        for k, v in kw.items():
            P[PNAMES.index(k)] = v
        return self.predict_array(P)

    def gradient(self, P, eps=1e-4):
        """numerical Jacobian dO/dP_k (exactly linear+quadratic, so central differences are exact up to rounding)."""
        P = np.asarray(P, float)
        J = np.empty((len(self.ids), 6))
        for k in range(6):
            dp = np.zeros(6); dp[k] = eps
            J[:, k] = (self.predict_array(P + dp) - self.predict_array(P - dp)) / (2 * eps)
        return J


def chi2_fit(T, data, sigma, free=('alphaE1', 'betaM1'), P0=None, norm=None, baldin=None):
    """Least-squares fit of selected polarisabilities.
    data, sigma : arrays aligned with T.ids
    free        : names of the free parameters (others fixed at P0, default Pref)
    norm        : optional dict {group_label: (mask, sigma_norm)} for floating normalisations of CS subsets
    baldin      : optional (value, error) Gaussian constraint on alphaE1+betaM1
    """
    from scipy.optimize import least_squares
    P0 = T.Pref.copy() if P0 is None else np.asarray(P0, float)
    idx = [PNAMES.index(f) for f in free]
    groups = list(norm.items()) if norm else []

    def resid(x):
        P = P0.copy(); P[idx] = x[:len(idx)]
        th = T.predict_array(P)
        r = []
        th = th.copy()
        for j, (lab, (mask, sn)) in enumerate(groups):
            th[mask] *= x[len(idx) + j]
        r.append((th - data) / sigma)
        for j, (lab, (mask, sn)) in enumerate(groups):
            r.append(np.atleast_1d((x[len(idx) + j] - 1.0) / sn))
        if baldin is not None:
            r.append(np.atleast_1d((P[0] + P[1] - baldin[0]) / baldin[1]))
        return np.concatenate(r)

    x0 = np.concatenate([P0[idx], np.ones(len(groups))])
    sol = least_squares(resid, x0, method='lm')
    J = sol.jac
    cov = np.linalg.inv(J.T @ J)
    chi2 = float(np.sum(sol.fun ** 2))
    return sol.x, np.sqrt(np.diag(cov)), cov, chi2


def _main():
    fn = sys.argv[1] if len(sys.argv) > 1 else 'test_points.out'
    T = DriverOutput(fn)
    print('Pref (reference = NNLO) :', np.round(T.Pref, 5))
    v0 = T.predict_array(T.Pref)
    v1 = T.predict(alphaE1=10.65, betaM1=3.15, gM1M1=2.2)
    chk = np.max(np.abs(T.predict_array(T.Puser) - T.value_file) / np.maximum(np.abs(T.value_file), 1e-12))
    print('max rel. diff between python evaluation at Puser and file value column: %.1e' % chk)
    print('%-10s %-5s %8s %8s %16s %16s' % ('id', 'obs', 'E_lab', 'th_lab', 'NNLO', 'a=10.65,b=3.15,gMM=2.2'))
    for i in range(len(T.ids)):
        print('%-10s %-5s %8.2f %8.2f %16.8g %16.8g' % (T.ids[i], T.obs[i], T.E[i], T.theta_lab[i], v0[i], v1[i]))
    if '--demo-fit' in sys.argv:
        rng = np.random.default_rng(1)
        Ptrue = T.Pref.copy(); Ptrue[0] = 10.6; Ptrue[1] = 3.2
        truth = T.predict_array(Ptrue)
        sig = np.where(np.char.startswith(T.obs.astype(str), 'S') & ~np.isin(T.obs, ['SPAR', 'SPERP']),
                       0.01, 0.03 * np.abs(truth))
        data = truth + sig * rng.standard_normal(truth.size)
        x, dx, cov, chi2 = chi2_fit(T, data, sig, free=('alphaE1', 'betaM1'))
        print('\ndemo fit (pseudo-data generated at alphaE1=10.6, betaM1=3.2):')
        print('  alphaE1 = %.3f +- %.3f   betaM1 = %.3f +- %.3f   chi2 = %.2f / %d dof'
              % (x[0], dx[0], x[1], dx[1], chi2, len(data) - 2))


if __name__ == '__main__':
    _main()
