"""
Partial-wave model and fitter for proton Compton scattering below the pion threshold,
following Krupina, Lensky, Pascalutsa, PLB 782 (2018) 34 (KLP18):
  * Born + pi0-pole amplitudes exact (compton.born_amplitudes)
  * non-Born multipoles l=1 parametrised by static polarisabilities + residual functions f^R_i per energy bin
    (KLP18 Eq. nonBorn_multipoles_expansion), l=2 (EE,MM) fixed to BChPT or zero, higher neglected
  * forward sum rules (source sumrules.SOURCE: ForwardRCS 2026 update by default, GHP 2015/2016 with FWD_SOURCE=ghp2016):
    either exact linear constraints eliminating f^R_2, f^R_3 and alpha+beta, gamma_M1M1 (KLP18), or forward amplitudes
    as chi^2 points with correlated nuisances.  Defaults ab, g0, sig_ab, sig_g0 = sumrules.AB_SR, G0_SR, SIG_AB, SIG_G0
    (2026: 14.009(208), -1.025(102); ghp2016: 14.0(2), -0.929(105)); fwd_unc=None takes the uncorrelated part of the
    forward amplitudes per bin from fwd.split() (ghp2016: 0.0029 mub GeV, the value used until 2026-10-18).
Parameter vector p (all polarisabilities in 1e-4 fm^3 / fm^4, residuals in 1e-4 fm^3 (i=1..4), 1e-4 fm^4 (i=5,6)):
  global: D = alpha-beta, gEE, gEM, gME   [+ nf, ng in 'chi2' sum-rule mode]
  per bin: fR1, fR4, fR5, fR6            [+ fR2, fR3 in 'chi2' mode]
  per normalisation group: n_k = (N_k - 1)/sigma_k  (if floating normalisations are requested)
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import compton as C
from sumrules import ForwardAmps, HBARC_MEV, AB_SR, SIG_AB, G0_SR, SIG_G0

HBARC = HBARC_MEV
MN = C.MN_GEV * 1e3                   # MeV (same as Fortran / Born)
MFM = MN / HBARC                      # fm^-1
U = 1e-4                              # 1e-4 fm^n
L1KEYS = [('EE', 1, '+'), ('EE', 1, '-'), ('MM', 1, '+'), ('MM', 1, '-'), ('EM', 1, '+'), ('ME', 1, '+')]
L2KEYS = [('EE', 2, '-'), ('EE', 2, '+'), ('MM', 2, '-'), ('MM', 2, '+')]
ALLKEYS = L1KEYS + L2KEYS
WCS = np.array([0.5, 0.5, 1.0, 1.0, 0.5, 0.5])


def unit_phis(theta_cm_rad):
    """D[i, m, k]: contribution of a unit multipole m (in ALLKEYS order) to Phi_k at angle i."""
    th = np.atleast_1d(theta_cm_rad)
    D = np.zeros((len(th), len(ALLKEYS), 6))
    for m, key in enumerate(ALLKEYS):
        pw = C.pw_from_multipoles({key: 1.0}, 5)
        D[:, m, :] = np.real(C.phis_from_pw(pw, th))
    return D


class Model:
    def __init__(self, data, bins, l2='bchpt', energy_mode='actual', fix_gEM=0.0, sr_mode='constraint',
                 ab=AB_SR, g0=G0_SR, sig_ab=SIG_AB, sig_g0=SIG_G0, fwd=None, float_norm=False,
                 prior_w=(10., 10., 10., 10., 1., 1.), fwd_err='nuisance', fwd_unc=None,
                 bchpt=None, l2_scale=1.0, prior_center='zero', smooth_w=None, fwd_scale=1.0):
        self.d = data.reset_index(drop=True).copy()
        self.bins = np.array(sorted(bins), float)
        self.l2 = l2
        self.l2_scale = l2_scale
        self.energy_mode = energy_mode
        self.fix_gEM = fix_gEM
        self.sr_mode = sr_mode
        self.ab, self.g0, self.sig_ab, self.sig_g0 = ab, g0, sig_ab, sig_g0
        self.fwd = fwd if fwd is not None else ForwardAmps()
        self.float_norm = float_norm
        self.prior_w = np.array(prior_w, float)
        self.prior_center = prior_center
        self.smooth_w = None if smooth_w is None else np.array(smooth_w, float)
        self.fwd_err = fwd_err
        self.fwd_unc = fwd_unc
        self.fwd_scale = fwd_scale
        if bchpt is None and (l2 == 'bchpt' or prior_center == 'bchpt'):
            from bchpt_multipoles import BChPT
            bchpt = BChPT()
        self.bchpt = bchpt
        self._setup()

    # ---------------- parameter layout ----------------
    def _setup(self):
        d = self.d
        if 'bin' not in d:
            d['bin'] = self.bins[np.argmin(abs(d.E.values[:, None] - self.bins[None, :]), axis=1)]
        d['ibin'] = [int(np.argmin(abs(self.bins - b))) for b in d['bin']]
        self.gnames = ['D', 'gEE', 'gEM', 'gME']
        if self.fix_gEM is not None:
            self.gnames.remove('gEM')
        if self.sr_mode == 'chi2':
            self.gnames += ['nf', 'ng']
        self.rnames = ['fR1', 'fR4', 'fR5', 'fR6'] if self.sr_mode == 'constraint' else \
            ['fR1', 'fR2', 'fR3', 'fR4', 'fR5', 'fR6']
        self.ng = len(self.gnames)
        self.nr = len(self.rnames)
        self.nb = len(self.bins)
        if self.float_norm:
            grp = d['norm_group'].fillna('').values
            ok = (grp != '') & (d['norm_err'].fillna(0).values > 0)
            self.ngroups = sorted(set(grp[ok]))
        else:
            self.ngroups = []
        self.npar = self.ng + self.nr * self.nb + len(self.ngroups)
        self.names = self.gnames + [f'{r}@{b:g}' for b in self.bins for r in self.rnames] + \
            [f'norm:{g}' for g in self.ngroups]
        # kinematics
        E = d['bin'].values.astype(float) if self.energy_mode == 'bin' else d['E'].values.astype(float)
        self.Eeval = E
        k = C.kin(E, d['theta_lab'].values, 'lab')
        self.k = k
        self.thcm = k['thcm']
        PB = np.real(C.phis_from_M(C.born_amplitudes(E, np.degrees(self.thcm), 'cm')))
        self.PB = PB                                     # (N,6) fm
        self.Dm = unit_phis(self.thcm)                   # (N,10,6)
        self.lf = C.lab_factor(k)
        # affine map p -> multipoles at each point: f = A p + a  (A: N x 10 x npar)
        self._build_affine()
        # observables
        self.obs = d['obs'].values
        self.y = d['value'].values.astype(float)
        self.e = d['err'].values.astype(float)
        # normalisation groups
        self.gidx = np.full(len(d), -1)
        if self.ngroups:
            for j, g in enumerate(self.ngroups):
                self.gidx[d['norm_group'].values == g] = j
            self.gsig = np.array([d.loc[d['norm_group'] == g, 'norm_err'].iloc[0] for g in self.ngroups])
        # forward chi2 points (at bin centres)
        if self.sr_mode == 'chi2':
            self._build_forward()
        self.mask = np.ones(len(d), bool)
        # prior centres for the residual functions
        self.prior_mu = np.zeros((self.nb, self.nr))
        if self.prior_center == 'bchpt':
            fr = bchpt_residuals(self.bchpt, self.bins)
            for j, nm in enumerate(self.rnames):
                self.prior_mu[:, j] = fr[nm]

    def pidx(self, name, b=None):
        if b is None:
            return self.names.index(name)
        return self.ng + b * self.nr + self.rnames.index(name)

    def multipoles_at(self, E, ib, p):
        """non-Born multipoles (fm) at lab energy E (MeV) using residuals of bin ib, parameter vector p."""
        g = dict(zip(self.gnames, p[:self.ng]))
        gEM = self.fix_gEM if self.fix_gEM is not None else g['gEM']
        if self.sr_mode == 'chi2':
            S = self.ab + self.sig_ab * g['nf']
            G0 = self.g0 + self.sig_g0 * g['ng']
        else:
            S, G0 = self.ab, self.g0
        Dd = g['D']
        a = 0.5 * (S + Dd) * U; b = 0.5 * (S - Dd) * U
        gEE = g['gEE'] * U; gME = g['gME'] * U; gEM = gEM * U
        gMM = -G0 * U - gEE - gEM - gME
        r = {nm: p[self.ng + ib * self.nr + j] * U for j, nm in enumerate(self.rnames)}
        x = E / HBARC
        s = 1.0 + 2.0 * E / MN
        pre2 = x**2 / np.sqrt(s)
        pre3 = x**3 / np.sqrt(s)
        R = (x / MFM)**2
        f = {}
        f[('EE', 1, '+')] = pre2 * (a / 3 + x / 3 * ((-a + b) / MFM + gEE) + R * r['fR1'])
        f[('MM', 1, '-')] = pre2 * (b / 3 + x / 3 * ((-b + a) / MFM - 2 * gMM) + R * r['fR4'])
        f[('EM', 1, '+')] = pre3 * (gEM / 6 + x / 6 * ((-6 * gEM + 3 * gME + 3 * gMM) / (4 * MFM) - b / (8 * MFM**2)) + R * r['fR5'])
        f[('ME', 1, '+')] = pre3 * (gME / 6 + x / 6 * ((-6 * gME + 3 * gEM + 3 * gEE) / (4 * MFM) - a / (8 * MFM**2)) + R * r['fR6'])
        EE1m0 = pre2 * (a / 3 + x / 3 * ((-a + b) / MFM - 2 * gEE))
        MM1p0 = pre2 * (b / 3 + x / 3 * ((-b + a) / MFM + gMM))
        # l = 2
        for key in L2KEYS:
            f[key] = self.l2_scale * self.bchpt(key, E) if self.l2 == 'bchpt' else 0.0 * E
        if self.sr_mode == 'constraint':
            sq = np.sqrt(s)
            fsr = self.fwd.fbar(E)
            gsr = self.fwd.gbar(E)
            # f = sq*(EE1- + 2EE1+ + MM1- + 2MM1+ + 6EE2- + 9EE2+ + 6MM2- + 9MM2+)
            restf = (EE1m0 + 2 * f[('EE', 1, '+')] + f[('MM', 1, '-')] + 2 * MM1p0 + 6 * f[('EE', 2, '-')]
                     + 9 * f[('EE', 2, '+')] + 6 * f[('MM', 2, '-')] + 9 * f[('MM', 2, '+')])
            restg = (EE1m0 - f[('EE', 1, '+')] - 6 * f[('EM', 1, '+')] - 6 * f[('ME', 1, '+')] + f[('MM', 1, '-')] - MM1p0
                     + 3 * f[('EE', 2, '-')] - 3 * f[('EE', 2, '+')] + 3 * f[('MM', 2, '-')] - 3 * f[('MM', 2, '+')])
            Rf = fsr / sq - restf          # = c*(fR2 + 2 fR3)
            Rg = gsr / sq - restg          # = c*(fR2 - fR3)
            c = pre2 * R
            fR3 = (Rf - Rg) / (3 * c)
            fR2 = Rg / c + fR3
            f[('EE', 1, '-')] = EE1m0 + c * fR2
            f[('MM', 1, '+')] = MM1p0 + c * fR3
        else:
            f[('EE', 1, '-')] = EE1m0 + pre2 * R * r['fR2']
            f[('MM', 1, '+')] = MM1p0 + pre2 * R * r['fR3']
        return f

    def forward_amps(self, E, ib, p):
        f = self.multipoles_at(E, ib, p)
        sq = np.sqrt(1.0 + 2.0 * E / MN)
        fb = sq * (f[('EE', 1, '-')] + 2 * f[('EE', 1, '+')] + f[('MM', 1, '-')] + 2 * f[('MM', 1, '+')] + 6 * f[('EE', 2, '-')]
                   + 9 * f[('EE', 2, '+')] + 6 * f[('MM', 2, '-')] + 9 * f[('MM', 2, '+')])
        gb = sq * (f[('EE', 1, '-')] - f[('EE', 1, '+')] - 6 * f[('EM', 1, '+')] - 6 * f[('ME', 1, '+')] + f[('MM', 1, '-')]
                   - f[('MM', 1, '+')] + 3 * f[('EE', 2, '-')] - 3 * f[('EE', 2, '+')] + 3 * f[('MM', 2, '-')] - 3 * f[('MM', 2, '+')])
        return fb, gb

    def _affine_for(self, E, ib):
        """returns A (len(E),10,npar_model) and a (len(E),10) such that f = A p + a (exact: linear model)."""
        npm = self.ng + self.nr * self.nb
        p0 = np.zeros(npm)
        f0 = self.multipoles_at(E, ib, p0)
        a = np.stack([np.broadcast_to(f0[kk], np.shape(E)) for kk in ALLKEYS], -1)
        A = np.zeros(np.shape(E) + (len(ALLKEYS), npm))
        idx = list(range(self.ng))
        for j in idx:
            pj = p0.copy(); pj[j] = 1.0
            fj = self.multipoles_at(E, ib, pj)
            A[..., :, j] = np.stack([np.broadcast_to(fj[kk], np.shape(E)) for kk in ALLKEYS], -1) - a
        ibs = np.atleast_1d(ib)
        for b in np.unique(ibs):
            for r in range(self.nr):
                j = self.ng + b * self.nr + r
                pj = p0.copy(); pj[j] = 1.0
                fj = self.multipoles_at(E, ib, pj)
                A[..., :, j] = np.stack([np.broadcast_to(fj[kk], np.shape(E)) for kk in ALLKEYS], -1) - a
        return A, a

    def _build_affine(self):
        N = len(self.d)
        npm = self.ng + self.nr * self.nb
        self.A = np.zeros((N, len(ALLKEYS), npm))
        self.a = np.zeros((N, len(ALLKEYS)))
        for b in range(self.nb):
            sel = np.where(self.d['ibin'].values == b)[0]
            if len(sel) == 0:
                continue
            A, a = self._affine_for(self.Eeval[sel], b)
            self.A[sel] = A
            self.a[sel] = a
        # Phi = PB + Dm . f  =>  Phi = u + V p
        self.u = self.PB + np.einsum('imk,im->ik', self.Dm, self.a)
        self.V = np.einsum('imk,imj->ikj', self.Dm, self.A)      # N x 6 x npm

    def _build_forward(self):
        Eb = self.bins
        fb_t = self.fwd.fbar(Eb); gb_t = self.fwd.gbar(Eb)
        ef, eg = self.fwd.err(Eb)
        if self.fwd_err == 'nuisance':
            cf, cg, uf, ug = self.fwd.split(Eb)
            if self.fwd_unc is not None:
                uf = ug = np.full(len(Eb), self.fwd_unc)
        rows_A, rows_a, y, sig, nus = [], [], [], [], []
        npm = self.ng + self.nr * self.nb
        for b, E in enumerate(Eb):
            p0 = np.zeros(npm)
            f0, g0 = self.forward_amps(np.array([E]), b, p0)
            Af = np.zeros(npm); Ag = np.zeros(npm)
            for j in list(range(self.ng)) + [self.ng + b * self.nr + r for r in range(self.nr)]:
                pj = p0.copy(); pj[j] = 1
                fj, gj = self.forward_amps(np.array([E]), b, pj)
                Af[j] = fj[0] - f0[0]; Ag[j] = gj[0] - g0[0]
            # nuisance shift of the data: y + n * sigma  <=> model - n*sigma
            if self.fwd_err == 'nuisance':
                Af[self.names.index('nf')] -= cf[b]
                Ag[self.names.index('ng')] -= cg[b]
                sf, sg = uf[b], ug[b]
            else:
                sf, sg = ef[b] * self.fwd_scale, eg[b] * self.fwd_scale
            rows_A += [Af, Ag]; rows_a += [f0[0], g0[0]]; y += [fb_t[b], gb_t[b]]; sig += [sf, sg]
        self.FA = np.array(rows_A); self.Fa = np.array(rows_a); self.Fy = np.array(y); self.Fs = np.array(sig)

    # ---------------- observables ----------------
    def predict(self, p, deriv=False, idx=None):
        npm = self.ng + self.nr * self.nb
        pm = p[:npm]
        sl = slice(None) if idx is None else idx
        Phi = self.u[sl] + np.einsum('ikj,j->ik', self.V[sl], pm)
        cs = np.sum(WCS * Phi**2, -1)
        num = -((Phi[:, 0] + Phi[:, 4]) * Phi[:, 2] + (Phi[:, 1] - Phi[:, 5]) * Phi[:, 3])
        lf = self.lf[sl] * C.NB_PER_FM2
        obs = self.obs[sl]
        out = np.where(obs == 'S3', num / cs, np.where(obs == 'CS_PAR', (cs + num) * lf,
                       np.where(obs == 'CS_PERP', (cs - num) * lf, cs * lf)))
        if not deriv:
            return out, None
        V = self.V[sl]
        dcs = np.einsum('ik,ikj->ij', 2 * WCS * Phi, V)
        dnum = -(np.einsum('i,ij->ij', Phi[:, 2], V[:, 0] + V[:, 4]) + np.einsum('i,ij->ij', Phi[:, 0] + Phi[:, 4], V[:, 2])
                 + np.einsum('i,ij->ij', Phi[:, 3], V[:, 1] - V[:, 5]) + np.einsum('i,ij->ij', Phi[:, 1] - Phi[:, 5], V[:, 3]))
        J = np.where((obs == 'S3')[:, None], dnum / cs[:, None] - (num / cs**2)[:, None] * dcs,
                     np.where((obs == 'CS_PAR')[:, None], (dcs + dnum) * lf[:, None],
                              np.where((obs == 'CS_PERP')[:, None], (dcs - dnum) * lf[:, None], dcs * lf[:, None])))
        return out, J

    def residuals(self, p, with_prior=True, jac=False):
        npm = self.ng + self.nr * self.nb
        idx = np.where(self.mask)[0]
        pred, J = self.predict(p, deriv=jac, idx=idx)
        y, e = self.y[idx], self.e[idx]
        Jfull = None
        if self.ngroups:
            gi = self.gidx[idx]
            scale = np.ones(len(idx))
            has = gi >= 0
            scale[has] = 1.0 + self.gsig[gi[has]] * p[npm + gi[has]]
            r = (pred * scale - y) / e
            if jac:
                Jfull = np.zeros((len(idx), self.npar))
                Jfull[:, :npm] = J * (scale / e)[:, None]
                for j in range(len(self.ngroups)):
                    s = gi == j
                    Jfull[s, npm + j] = pred[s] * self.gsig[j] / e[s]
        else:
            r = (pred - y) / e
            if jac:
                Jfull = np.zeros((len(idx), self.npar))
                Jfull[:, :npm] = J / e[:, None]
        res = [r]
        jacs = [Jfull] if jac else []
        if self.sr_mode == 'chi2':
            fr = (self.FA @ p[:npm] + self.Fa - self.Fy) / self.Fs
            res.append(fr)
            if jac:
                Jf = np.zeros((len(fr), self.npar)); Jf[:, :npm] = self.FA / self.Fs[:, None]
                jacs.append(Jf)
            # nuisance penalties
            for nm in ('nf', 'ng'):
                j = self.names.index(nm)
                res.append(np.array([p[j]]))
                if jac:
                    Jn = np.zeros((1, self.npar)); Jn[0, j] = 1.0; jacs.append(Jn)
        if self.ngroups:
            res.append(p[npm:])
            if jac:
                Jn = np.zeros((len(self.ngroups), self.npar)); Jn[:, npm:] = np.eye(len(self.ngroups)); jacs.append(Jn)
        if with_prior:
            w = self.prior_w[[['fR1', 'fR2', 'fR3', 'fR4', 'fR5', 'fR6'].index(n) for n in self.rnames]]
            pr = (p[self.ng:npm].reshape(self.nb, self.nr) - self.prior_mu) / w
            res.append(pr.ravel())
            if jac:
                Jp = np.zeros((self.nb * self.nr, self.npar))
                Jp[:, self.ng:npm] = np.diag(np.tile(1 / w, self.nb))
                jacs.append(Jp)
            if self.smooth_w is not None and self.nb >= 3:
                # second divided differences of each residual function, scaled to a 10 MeV step
                sw = self.smooth_w[[['fR1', 'fR2', 'fR3', 'fR4', 'fR5', 'fR6'].index(n) for n in self.rnames]]
                E = self.bins
                rows = []
                for j in range(1, self.nb - 1):
                    h1, h2 = E[j] - E[j - 1], E[j + 1] - E[j]
                    c = 2.0 / (h1 * h2 * (h1 + h2)) * 100.0        # (10 MeV)^2 * f''
                    coef = {j - 1: c * h2, j: -c * (h1 + h2), j + 1: c * h1}
                    for r in range(self.nr):
                        row = np.zeros(self.npar)
                        for jj, cc in coef.items():
                            row[self.ng + jj * self.nr + r] = cc / sw[r]
                        rows.append(row)
                S = np.array(rows)
                res.append(S[:, :npm] @ p[:npm])
                if jac:
                    jacs.append(S)
        rr = np.concatenate(res)
        return (rr, np.vstack(jacs)) if jac else rr

    def chi2_data(self, p):
        idx = np.where(self.mask)[0]
        pred, _ = self.predict(p, idx=idx)
        npm = self.ng + self.nr * self.nb
        if self.ngroups:
            gi = self.gidx[idx]; scale = np.ones(len(idx)); has = gi >= 0
            scale[has] = 1.0 + self.gsig[gi[has]] * p[npm + gi[has]]
            pred = pred * scale
        return np.sum(((pred - self.y[idx]) / self.e[idx])**2)

    def pulls(self, p):
        pred, _ = self.predict(p)
        npm = self.ng + self.nr * self.nb
        if self.ngroups:
            gi = self.gidx; scale = np.ones(len(pred)); has = gi >= 0
            scale[has] = 1.0 + self.gsig[gi[has]] * p[npm + gi[has]]
            pred = pred * scale
        return (self.y - pred) / self.e, pred

    # ---------------- fitting ----------------
    def p_start(self):
        p = np.zeros(self.npar)
        st = {'D': 9.5, 'gEE': -3.3, 'gEM': 0.2, 'gME': 1.1, 'nf': 0.0, 'ng': 0.0}
        for j, n in enumerate(self.gnames):
            p[j] = st[n]
        return p

    def fit(self, p0=None, with_prior=True, fixed=None):
        """joint minimisation. fixed: dict name->value of parameters to hold fixed (by index)."""
        if p0 is None:
            p0 = self.p_start()
        free = np.ones(self.npar, bool)
        pf = p0.copy()
        if fixed:
            for nm, v in fixed.items():
                j = self.names.index(nm); free[j] = False; pf[j] = v

        def fun(q):
            pp = pf.copy(); pp[free] = q
            r, J = self.residuals(pp, with_prior, jac=True)
            return r, J[:, free]
        sol = least_squares(lambda q: fun(q)[0], pf[free], jac=lambda q: fun(q)[1], method='lm', xtol=1e-12, ftol=1e-12,
                            gtol=1e-12, max_nfev=2000)
        p = pf.copy(); p[free] = sol.x
        r, J = self.residuals(p, with_prior, jac=True)
        Jf = J[:, free]
        H = Jf.T @ Jf
        cov = np.zeros((self.npar, self.npar))
        cov[np.ix_(free, free)] = np.linalg.pinv(H)
        res = FitResult(self, p, cov, free, with_prior)
        return res

    def fit_globals_only(self, p, with_prior=False):
        """KLP18 stage 2: residual functions held fixed at the values in p; minimise over global parameters."""
        fixed = {n: p[self.names.index(n)] for n in self.names if n not in self.gnames}
        return self.fit(p0=p, with_prior=with_prior, fixed=fixed)

    def curve(self, E, theta_lab, ib, p, obs='cs_lab'):
        """model observable at lab energy E (MeV) and lab angles, multipoles from bin ib residuals."""
        th = np.atleast_1d(np.asarray(theta_lab, float))
        k = C.kin(E, th, 'lab')
        f = self.multipoles_at(np.array([float(E)]), ib, p)
        fd = {kk: float(np.asarray(v).ravel()[0]) for kk, v in f.items()}
        Pnb = np.real(C.phis_from_pw(C.pw_from_multipoles(fd, 5), k['thcm']))
        PB = np.real(C.phis_from_M(C.born_amplitudes(E, np.degrees(k['thcm']), 'cm')))
        o = C.obs_from_phis(PB + Pnb, k)
        return o[obs] if obs else o

    def curve_err(self, E, theta_lab, ib, p, cov, obs='cs_lab'):
        """model observable and its 1-sigma error (linear propagation of the full parameter covariance)."""
        th = np.atleast_1d(np.asarray(theta_lab, float))
        k = C.kin(E, th, 'lab')
        A, a = self._affine_for(np.array([float(E)]), ib)
        A, a = A[0], a[0]                                       # 10 x npm, 10
        Dm = unit_phis(k['thcm'])                                # N x 10 x 6
        PB = np.real(C.phis_from_M(C.born_amplitudes(E, np.degrees(k['thcm']), 'cm')))
        npm = A.shape[1]
        u = PB + np.einsum('imk,m->ik', Dm, a)
        V = np.einsum('imk,mj->ikj', Dm, A)
        Phi = u + np.einsum('ikj,j->ik', V, p[:npm])
        cs = np.sum(WCS * Phi**2, -1)
        num = -((Phi[:, 0] + Phi[:, 4]) * Phi[:, 2] + (Phi[:, 1] - Phi[:, 5]) * Phi[:, 3])
        dcs = np.einsum('ik,ikj->ij', 2 * WCS * Phi, V)
        dnum = -(Phi[:, 2, None] * (V[:, 0] + V[:, 4]) + (Phi[:, 0] + Phi[:, 4])[:, None] * V[:, 2]
                 + Phi[:, 3, None] * (V[:, 1] - V[:, 5]) + (Phi[:, 1] - Phi[:, 5])[:, None] * V[:, 3])
        lf = C.lab_factor(k) * C.NB_PER_FM2
        if obs == 'cs_lab':
            val, J = cs * lf, dcs * lf[:, None]
        elif obs == 'S3':
            val, J = num / cs, dnum / cs[:, None] - (num / cs**2)[:, None] * dcs
        elif obs == 'cs_par_lab':
            val, J = (cs + num) * lf, (dcs + dnum) * lf[:, None]
        elif obs == 'cs_perp_lab':
            val, J = (cs - num) * lf, (dcs - dnum) * lf[:, None]
        err = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', J, cov[:npm, :npm], J), 0))
        return val, err

    def multipoles_bin(self, ib, p, E=None):
        E = self.bins[ib] if E is None else E
        f = self.multipoles_at(np.array([float(E)]), ib, p)
        return {kk: float(np.asarray(v).ravel()[0]) for kk, v in f.items()}

    def multipole_cov(self, ib, cov, p, E=None):
        """covariance-propagated errors of the 10 non-Born multipoles at bin ib (linear map)."""
        E = self.bins[ib] if E is None else E
        A, a = self._affine_for(np.array([float(E)]), ib)
        A = A[0]                                   # 10 x npm
        npm = A.shape[1]
        C_ = cov[:npm, :npm]
        return np.sqrt(np.maximum(np.einsum('mj,jk,mk->m', A, C_, A), 0))

    # ---------------- derived quantities ----------------
    def derived(self, p):
        g = dict(zip(self.gnames, p[:self.ng]))
        gEM = self.fix_gEM if self.fix_gEM is not None else g['gEM']
        S = self.ab + (self.sig_ab * g['nf'] if 'nf' in g else 0.0)
        G0 = self.g0 + (self.sig_g0 * g['ng'] if 'ng' in g else 0.0)
        a = 0.5 * (S + g['D']); b = 0.5 * (S - g['D'])
        gMM = -G0 - g['gEE'] - gEM - g['gME']
        gpi = -g['gEE'] + gMM - gEM + g['gME']
        return dict(alphaE1=a, betaM1=b, alpha_m_beta=g['D'], alpha_p_beta=S, gE1E1=g['gEE'], gM1M1=gMM,
                    gE1M2=gEM, gM1E2=g['gME'], gamma0=G0, gamma_pi=gpi)


def bchpt_residuals(bchpt, Ebins, pols=None):
    """residual functions f^R_i(E) (1e-4 fm^3 / fm^4) implied by the NNLO BChPT non-Born l=1 multipoles
    and the BChPT static polarisabilities in the KLP18 parametrisation."""
    p = pols or dict(a=11.2461, b=3.87054, gEE=-3.34375, gMM=2.95934, gEM=0.185504, gME=1.14328)
    a, b, gEE, gMM, gEM, gME = [p[k] * U for k in ('a', 'b', 'gEE', 'gMM', 'gEM', 'gME')]
    out = {k: [] for k in ('fR1', 'fR2', 'fR3', 'fR4', 'fR5', 'fR6')}
    for E in np.atleast_1d(Ebins):
        x = E / HBARC; s = 1 + 2 * E / MN; p2 = x**2 / np.sqrt(s); p3 = x**3 / np.sqrt(s); R = (x / MFM)**2
        t = {'fR1': (('EE', 1, '+'), p2, a / 3 + x / 3 * ((-a + b) / MFM + gEE)),
             'fR2': (('EE', 1, '-'), p2, a / 3 + x / 3 * ((-a + b) / MFM - 2 * gEE)),
             'fR3': (('MM', 1, '+'), p2, b / 3 + x / 3 * ((-b + a) / MFM + gMM)),
             'fR4': (('MM', 1, '-'), p2, b / 3 + x / 3 * ((-b + a) / MFM - 2 * gMM)),
             'fR5': (('EM', 1, '+'), p3, gEM / 6 + x / 6 * ((-6 * gEM + 3 * gME + 3 * gMM) / (4 * MFM) - b / (8 * MFM**2))),
             'fR6': (('ME', 1, '+'), p3, gME / 6 + x / 6 * ((-6 * gME + 3 * gEM + 3 * gEE) / (4 * MFM) - a / (8 * MFM**2)))}
        for k, (key, pre, st) in t.items():
            out[k].append((bchpt(key, E) / pre - st) / R / U)
    return {k: np.array(v) for k, v in out.items()}


class FitResult:
    def __init__(self, model, p, cov, free, with_prior):
        self.m = model
        self.p = p
        self.cov = cov
        self.free = free
        self.chi2 = model.chi2_data(p)
        self.npts = int(model.mask.sum())
        self.chi2_total = float(np.sum(model.residuals(p, with_prior)**2))
        self.der = model.derived(p)
        # linear error propagation for derived quantities (numerical gradient, exact for linear maps)
        self.der_err = {}
        eps = 1e-6
        base = model.derived(p)
        G = {k: np.zeros(model.npar) for k in base}
        for j in range(model.ng):
            q = p.copy(); q[j] += eps
            dd = model.derived(q)
            for k in base:
                G[k][j] = (dd[k] - base[k]) / eps
        for k in base:
            self.der_err[k] = float(np.sqrt(max(G[k] @ cov @ G[k], 0)))
        self.G = G

    def corr(self, a, b):
        ga, gb = self.G[a], self.G[b]
        va, vb = ga @ self.cov @ ga, gb @ self.cov @ gb
        return float(ga @ self.cov @ gb / np.sqrt(va * vb)) if va > 0 and vb > 0 else np.nan

    def summary(self, label=''):
        d, e = self.der, self.der_err
        s = (f'{label:14s} chi2/pt={self.chi2/self.npts:5.3f} ({self.chi2:6.1f}/{self.npts})  '
             f'aE1={d["alphaE1"]:6.2f}({e["alphaE1"]:.2f}) bM1={d["betaM1"]:5.2f}({e["betaM1"]:.2f}) '
             f'gEE={d["gE1E1"]:5.2f}({e["gE1E1"]:.2f}) gMM={d["gM1M1"]:5.2f}({e["gM1M1"]:.2f}) '
             f'gEM={d["gE1M2"]:5.2f}({e["gE1M2"]:.2f}) gME={d["gM1E2"]:5.2f}({e["gM1E2"]:.2f}) '
             f'a-b={d["alpha_m_beta"]:5.2f}({e["alpha_m_beta"]:.2f}) gpi={d["gamma_pi"]:5.2f}({e["gamma_pi"]:.2f})')
        return s
