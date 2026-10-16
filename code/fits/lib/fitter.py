"""
fitter.py -- chi^2 machinery for the polarisability fits.

chi^2 = (y - m)^T C^{-1} (y - m) + sum_k lam_k^2 + sum_c ((c.P - v_c)/s_c)^2

  m_i   = n_{k(i)} * t_i(P, extra) + sum_m lam_m * s_{m,i}         (model)
  n_k   = 1 + delta_k * lam_k                                       (floating normalisation of set k;
                                                                     penalty lam_k^2 = ((n_k-1)/delta_k)^2)
  s_m   = additive correlated shift vectors (e.g. GHP forward-amplitude nuisances), penalty lam_m^2
  C     = point-to-point covariance (diagonal stat (+) p2p syst, + uncorrelated forward blocks,
          + optionally the EFT truncation covariance)
  P     = P0 + A x  (free physics parameters x; P = 6 static polarisabilities)

Degrees of freedom (GMPF/MPG convention): dof = N_points + N_constraints - N_free - N_norm - N_shift.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2dist
from scipy.linalg import cho_factor, cho_solve, solve_triangular

PNAMES = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
# useful parameter directions (columns of A)
DIR = {
    'alphaE1': np.array([1, 0, 0, 0, 0, 0.]), 'betaM1': np.array([0, 1, 0, 0, 0, 0.]),
    'gE1E1': np.array([0, 0, 1, 0, 0, 0.]), 'gM1M1': np.array([0, 0, 0, 1, 0, 0.]),
    'gE1M2': np.array([0, 0, 0, 0, 1, 0.]), 'gM1E2': np.array([0, 0, 0, 0, 0, 1.]),
    # gamma_pi direction at fixed gamma_0: dgM1M1=+1/2, dgE1E1=-1/2  -> dgamma_pi=1, dgamma_0=0
    'gpi@g0': np.array([0, 0, -0.5, 0.5, 0, 0.]),
}
C_BALDIN = np.array([1, 1, 0, 0, 0, 0.])
C_GAMMA0 = np.array([0, 0, -1, -1, -1, -1.])
C_GAMMAPI = np.array([0, 0, -1, 1, -1, 1.])


def gamma0(P):
    return -(P[2] + P[3] + P[4] + P[5])


def gammapi(P):
    return P[3] + P[5] - P[2] - P[4]


class FitResult(dict):
    __getattr__ = dict.__getitem__


class Chi2Fit:
    def __init__(self, y, theory, C, datasets, norms, P0, free, constraints=(), shifts=None,
                 extra=None, extra_prior=None, fix_norms=None, labels=None):
        """
        y        : (N,) data
        theory   : callable theory(P, extra) -> (N,) predictions (no normalisation)
        C        : (N,N) covariance or (N,) variances
        datasets : (N,) dataset label per point
        norms    : dict dataset -> delta (fractional 1-sigma); 0 or missing => not floated
        P0       : (6,) starting / fixed polarisabilities
        free     : list of direction names (keys of DIR) or 6-vectors
        constraints : list of (cvec(6), value, sigma, name)
        shifts   : (Ns, N) additive nuisance shift vectors (1 sigma), or None
        extra    : dict name -> start value for non-linear extra parameters (e.g. gM); passed to theory
        extra_prior : dict name -> (mean, sigma) Gaussian prior (optional)
        fix_norms: dict dataset -> fixed lam (for error decomposition), or 'all' to fix all at given values
        """
        self.y = np.asarray(y, float)
        self.N = len(self.y)
        self.theory = theory
        C = np.asarray(C, float)
        if C.ndim == 1:
            self.diag = True
            self.sig = np.sqrt(C)
        else:
            self.diag = False
            self.L = np.linalg.cholesky(C)
            self.C = C
        self.datasets = np.asarray(datasets)
        self.norm_sets = [k for k in dict.fromkeys(self.datasets) if norms.get(k, 0) > 0]
        self.norm_delta = np.array([norms[k] for k in self.norm_sets])
        self.norm_idx = [np.where(self.datasets == k)[0] for k in self.norm_sets]
        self.P0 = np.asarray(P0, float)
        self.free = list(free)
        self.A = np.array([DIR[f] if isinstance(f, str) else np.asarray(f, float) for f in self.free]).T.reshape(6, -1)
        self.constraints = list(constraints)
        self.shifts = np.zeros((0, self.N)) if shifts is None else np.atleast_2d(shifts)
        self.extra_names = list(extra.keys()) if extra else []
        self.extra0 = np.array([extra[k] for k in self.extra_names]) if extra else np.zeros(0)
        self.extra_prior = extra_prior or {}
        self.fix_norms = fix_norms
        self.labels = labels

    # ---- parameter vector layout: [x (free), extra, lam_norm, lam_shift]
    def _split(self, q):
        nf, ne, nn, ns = len(self.free), len(self.extra_names), len(self.norm_sets), len(self.shifts)
        x = q[:nf]; e = q[nf:nf + ne]; ln = q[nf + ne:nf + ne + nn]; ls = q[nf + ne + nn:]
        return x, e, ln, ls

    def P_of(self, x):
        return self.P0 + self.A @ x

    def model(self, q, with_norm=True):
        x, e, ln, ls = self._split(q)
        P = self.P_of(x)
        ex = dict(zip(self.extra_names, e))
        t = np.array(self.theory(P, ex), float)
        nvec = np.ones(self.N)
        if with_norm:
            for j, idx in enumerate(self.norm_idx):
                nvec[idx] = 1.0 + self.norm_delta[j] * ln[j]
        m = nvec * t
        if len(ls):
            m = m + ls @ self.shifts
        return m, t, nvec, P

    def whiten(self, r):
        if self.diag:
            return r / self.sig
        return solve_triangular(self.L, r, lower=True)

    def residuals(self, q):
        x, e, ln, ls = self._split(q)
        m, t, nvec, P = self.model(q)
        r = [self.whiten(self.y - m)]
        r.append(ln)
        r.append(ls)
        for (c, v, s, nm) in self.constraints:
            r.append(np.atleast_1d((c @ P - v) / s))
        for k, (mu, s) in self.extra_prior.items():
            r.append(np.atleast_1d((e[self.extra_names.index(k)] - mu) / s))
        return np.concatenate(r)

    def fit(self, x0=None, fixed_mask=None, q0=None):
        nf, ne, nn, ns = len(self.free), len(self.extra_names), len(self.norm_sets), len(self.shifts)
        if q0 is None:
            q0 = np.concatenate([np.zeros(nf) if x0 is None else x0, self.extra0, np.zeros(nn), np.zeros(ns)])
        if fixed_mask is None:
            fixed_mask = np.zeros(len(q0), bool)
        fm = np.asarray(fixed_mask, bool)
        qfix = q0.copy()

        def res(qf):
            q = qfix.copy(); q[~fm] = qf
            return self.residuals(q)
        sol = least_squares(res, q0[~fm], method='lm', xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=20000)
        q = qfix.copy(); q[~fm] = sol.x
        # accurate Jacobian by central differences
        J = self._jac(res, sol.x)
        H = J.T @ J
        try:
            covf = np.linalg.inv(H)
        except np.linalg.LinAlgError:
            covf = np.linalg.pinv(H)
        cov = np.zeros((len(q), len(q)))
        ii = np.where(~fm)[0]
        cov[np.ix_(ii, ii)] = covf
        return q, cov, sol

    @staticmethod
    def _jac(res, x, h=None):
        x = np.asarray(x, float)
        f0 = res(x)
        J = np.empty((len(f0), len(x)))
        for j in range(len(x)):
            hj = 1e-5 * max(1.0, abs(x[j])) if h is None else h
            xp = x.copy(); xp[j] += hj
            xm = x.copy(); xm[j] -= hj
            J[:, j] = (res(xp) - res(xm)) / (2 * hj)
        return J

    def summarize(self, q, cov):
        x, e, ln, ls = self._split(q)
        m, t, nvec, P = self.model(q)
        rr = self.y - m
        rw = self.whiten(rr)
        chi2_data = float(rw @ rw)
        # per-point chi^2 contribution (sums to chi2_data): r_i (C^-1 r)_i
        if self.diag:
            contrib = (rr / self.sig) ** 2
        else:
            Cinv_r = cho_solve((self.L, True), rr)
            contrib = rr * Cinv_r
        chi2_norm = float(ln @ ln)
        chi2_shift = float(ls @ ls)
        chi2_con = 0.0
        for (c, v, s, nm) in self.constraints:
            chi2_con += float(((c @ P - v) / s) ** 2)
        for k, (mu, s) in self.extra_prior.items():
            chi2_con += float(((e[self.extra_names.index(k)] - mu) / s) ** 2)
        chi2 = chi2_data + chi2_norm + chi2_shift + chi2_con
        nf = len(self.free) + len(self.extra_names)
        dof = self.N + len(self.constraints) + len(self.extra_prior) - nf - len(self.norm_sets) - len(self.shifts)
        # parameter covariance of P
        covx = cov[:len(self.free), :len(self.free)]
        covP = self.A @ covx @ self.A.T
        perset = {}
        for k in dict.fromkeys(self.datasets):
            idx = np.where(self.datasets == k)[0]
            d = dict(N=len(idx), chi2=float(contrib[idx].sum()))
            if k in self.norm_sets:
                j = self.norm_sets.index(k)
                d['n'] = 1 + self.norm_delta[j] * ln[j]
                d['dn'] = self.norm_delta[j] * np.sqrt(cov[len(self.free) + len(self.extra_names) + j,
                                                          len(self.free) + len(self.extra_names) + j])
                d['chi2_norm'] = float(ln[j] ** 2)
            else:
                d['n'] = 1.0; d['dn'] = 0.0; d['chi2_norm'] = 0.0
            perset[k] = d
        ie = len(self.free)
        extra = {k: (e[i], np.sqrt(cov[ie + i, ie + i])) for i, k in enumerate(self.extra_names)}
        return FitResult(P=P, covP=covP, dP=np.sqrt(np.clip(np.diag(covP), 0, None)), x=x, covx=covx, extra=extra,
                         chi2=chi2, chi2_data=chi2_data, chi2_norm=chi2_norm, chi2_shift=chi2_shift, chi2_con=chi2_con,
                         dof=dof, pvalue=float(chi2dist.sf(chi2, dof)) if dof > 0 else np.nan,
                         perset=perset, model=m, theory=t, nvec=nvec, resid=rr, contrib=contrib, lam_shift=ls,
                         q=q, cov=cov)

    def run(self, decompose=True):
        """full fit; optionally error decomposition by fixing nuisance groups at their best-fit values."""
        q, cov, sol = self.fit()
        R = self.summarize(q, cov)
        if decompose:
            nf, ne, nn, ns = len(self.free), len(self.extra_names), len(self.norm_sets), len(self.shifts)
            n = len(q)
            # (1) all nuisances fixed -> "stat"
            fm = np.zeros(n, bool); fm[nf + ne:] = True
            q1, c1, _ = self.fit(fixed_mask=fm, q0=q.copy())
            R['covP_stat'] = self.A @ c1[:nf, :nf] @ self.A.T
            # (2) norms floated, shifts fixed -> stat+norm
            fm = np.zeros(n, bool); fm[nf + ne + nn:] = True
            q2, c2, _ = self.fit(fixed_mask=fm, q0=q.copy())
            R['covP_statnorm'] = self.A @ c2[:nf, :nf] @ self.A.T
        return R


def corr_from_cov(c):
    d = np.sqrt(np.diag(c))
    return c / np.outer(d, d)
