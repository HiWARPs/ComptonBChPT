#!/usr/bin/env python3
"""
Tree-level (Born) pion photoproduction gamma N -> pi N: CGLN amplitudes and multipoles,
computed numerically from explicit Dirac algebra (independent of the Fortran code).
Used for the gamma-N-Delta vertex corrections of ampl.exe v2 (Watson/unitarity construction)
and for the unitarity checks of the regression tests.

Graphs: nucleon s- and u-channel poles (charge + anomalous magnetic moment), t-channel pion
pole, Kroll-Ruderman contact (pseudovector coupling only). Relative signs of the graphs are
fixed by electromagnetic gauge invariance (checked numerically, k_mu J^mu = 0), the overall
sign by the low-energy theorems in the standard CGLN/MAID convention:
   E0+(gamma p -> pi+ n) ~ +27,  E0+(gamma p -> pi0 p) ~ -2.4  [1e-3/m_pi+]  (Born, with kappa).
Isospin: T^a = chi_f^+[ delta_a3 A(+) + 1/2[tau_a,tau_3] A(-) + tau_a A(0) ]chi_i,
   T(pi0 p) = T^3 = A(+)+A(0),  T(pi+ n) = (T^1 - i T^2)/sqrt2 = sqrt2 (A(-)+A(0)),
   A^{3/2} = A(+) - A(-) = T(pi0 p) - T(pi+ n)/sqrt2       (MAID normalisation)
   CG-normalised isospin-3/2 amplitude (sum_c |T_c|^2 = |T_3/2,CG|^2 + |T_1/2,CG|^2): sqrt(2/3) A^{3/2}.
CGLN: chi_f^+ F chi_i = T/(8 pi W) with ubar u = 2M spinors; dsigma/dOmega = (q/k)|chi^+ F chi|^2;
   F = i s.e F1 + s.qh s.(kh x e) F2 + i s.kh qh.e F3 + i s.qh qh.e F4 (c.m. frame, k along z).
Units: GeV, multipoles in GeV^-1 (mpi_unit() converts to 1e-3/m_pi+).
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
g0 = np.block([[I2, Z2], [Z2, -I2]])
GAM = [g0] + [np.block([[Z2, s], [-s, Z2]]) for s in (sx, sy, sz)]
G5 = np.block([[Z2, I2], [I2, Z2]])
METRIC = np.diag([1.0, -1.0, -1.0, -1.0])
I4 = np.eye(4, dtype=complex)
TAU = [sx, sy, sz]
QN = 0.5 * (I2 + sz)


def slash(p):
    return p[0] * GAM[0] - p[1] * GAM[1] - p[2] * GAM[2] - p[3] * GAM[3]


def mdot(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


SIGMN = [[0.5j * (GAM[m] @ GAM[n] - GAM[n] @ GAM[m]) for n in range(4)] for m in range(4)]


def spinor(p, M, chi):
    sp = p[1] * sx + p[2] * sy + p[3] * sz
    return np.concatenate([np.sqrt(p[0] + M) * chi, (sp @ chi) / np.sqrt(p[0] + M)])


def levi3(i, j, k):
    if len({i, j, k}) < 3:
        return 0
    return {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1}.get((i, j, k), -1)


class Born:
    """Born gamma N -> pi N amplitudes. Sign parameters fixed by gauge invariance (see check())."""

    def __init__(self, M=0.93827, mpi=0.13957, g=None, gA=1.27, fpi=0.09221, e=None,
                 kp=0.0, kn=0.0, coupling='PV', st=1.0, skr=1.0, overall=1.0):
        self.M, self.mpi = M, mpi
        self.g = g if g is not None else gA * M / fpi
        self.e = e if e is not None else np.sqrt(4 * np.pi / 137.036)
        self.kp, self.kn = kp, kn
        self.coupling = coupling
        self.st, self.skr, self.overall = st, skr, overall   # overall: global phase (+-1, +-i)

    def kin(self, W, cth):
        M, m = self.M, self.mpi
        s = W * W
        k = (s - M * M) / (2 * W)
        Eq = (s + m * m - M * M) / (2 * W)
        q = np.sqrt(max(Eq * Eq - m * m, 0.0))
        sth = np.sqrt(max(0.0, 1 - cth * cth))
        kv = np.array([k, 0, 0, k])
        qv = np.array([Eq, q * sth, 0, q * cth])
        pv = np.array([np.sqrt(k * k + M * M), 0, 0, -k])
        ppv = np.array([np.sqrt(q * q + M * M), -q * sth, 0, -q * cth])
        return dict(W=W, k=k, q=q, Eq=Eq, kv=kv, qv=qv, pv=pv, ppv=ppv, sth=sth, cth=cth)

    def vpi(self, pq):
        if self.coupling == 'PV':
            return (self.g / (2 * self.M)) * slash(pq) @ G5
        return self.g * G5

    def currents(self, K):
        """J[a][mu] : (2,2,4,4) isospin x Dirac operator, a = 0,1,2 (Cartesian)."""
        M = self.M
        kv, qv, pv = K['kv'], K['qv'], K['pv']
        kl = METRIC @ kv
        kS, kV = self.kp + self.kn, self.kp - self.kn
        Ss = np.linalg.inv(slash(pv + kv) - M * I4)
        Su = np.linalg.inv(slash(pv - qv) - M * I4)
        Vq = self.vpi(qv)
        tq = qv - kv
        Vt = self.vpi(tq)
        pole = 1.0 / (mdot(tq, tq) - self.mpi ** 2)
        J = []
        for a in range(3):
            Ja = np.zeros((4, 2, 2, 4, 4), dtype=complex)
            for mu in range(4):
                mag = sum(SIGMN[mu][nu] * kl[nu] for nu in range(4)) * (1j / (2 * M))
                # gamma NN vertex as isospin(2x2) x Dirac(4x4): Q gamma^mu + (kS+kV tau3)/2 i sigma k/2M
                Vg = np.zeros((2, 2, 4, 4), dtype=complex)
                for i in range(2):
                    for j in range(2):
                        Vg[i, j] = QN[i, j] * GAM[mu] + 0.5 * (kS * I2[i, j] + kV * sz[i, j]) * mag
                for i in range(2):
                    for j in range(2):
                        acc = np.zeros((4, 4), dtype=complex)
                        for m in range(2):
                            acc += TAU[a][i, m] * (Vq @ Ss @ Vg[m, j])      # s-channel
                            acc += (Vg[i, m] @ Su @ Vq) * TAU[a][m, j]      # u-channel
                        Ja[mu, i, j] += acc
                jt = (2 * qv - kv)[mu]
                for b in range(3):
                    ev = levi3(2, a, b)
                    if ev == 0:
                        continue
                    iso = (-1j * ev) * TAU[b]
                    Ja[mu] += self.st * np.einsum('ij,kl->ijkl', iso, Vt * (jt * pole))   # pion pole
                    if self.coupling == 'PV':
                        Ja[mu] += self.skr * np.einsum('ij,kl->ijkl', iso,
                                                       (self.g / (2 * M)) * GAM[mu] @ G5)  # Kroll-Ruderman
            J.append(Ja * (self.e * self.overall))
        return J

    @staticmethod
    def channel_op(J, mu, channel):
        if channel == 'pi0p':
            return J[2][mu][0, 0]
        if channel == 'pi+n':
            return (J[0][mu][1, 0] - 1j * J[1][mu][1, 0]) / np.sqrt(2)
        if channel == 'pi0n':
            return J[2][mu][1, 1]
        if channel == 'pi-p':
            return (J[0][mu][0, 1] + 1j * J[1][mu][0, 1]) / np.sqrt(2)
        raise ValueError(channel)

    def matrix(self, K, channel, epsvec):
        """2x2 matrix T_{fi} (Pauli spin basis) for polarisation 4-vector eps (upper index)."""
        J = self.currents(K)
        chis = [np.array([1, 0], dtype=complex), np.array([0, 1], dtype=complex)]
        out = np.zeros((2, 2), dtype=complex)
        epsl = METRIC @ epsvec
        Ops = [self.channel_op(J, mu, channel) for mu in range(4)]
        Op = sum(epsl[mu] * Ops[mu] for mu in range(4))
        for f, chf in enumerate(chis):
            ubf = spinor(K['ppv'], self.M, chf).conj() @ g0
            for i, chi in enumerate(chis):
                out[f, i] = ubf @ Op @ spinor(K['pv'], self.M, chi)
        return out

    def gauge_check(self, W=1.2, cth=0.3, channel='pi+n'):
        K = self.kin(W, cth)
        Tk = self.matrix(K, channel, K['kv'])
        Te = self.matrix(K, channel, np.array([0, 1.0, 0, 0]))
        return np.max(np.abs(Tk)) / np.max(np.abs(Te))

    def cgln(self, W, cth, channel):
        K = self.kin(W, cth)
        sth = K['sth']
        qh = np.array([sth, 0, cth])
        kh = np.array([0, 0, 1.0])

        def sd(v):
            return v[0] * sx + v[1] * sy + v[2] * sz
        A, b = [], []
        for e4, e3 in [(np.array([0, 1.0, 0, 0]), np.array([1.0, 0, 0])),
                       (np.array([0, 0, 1.0, 0]), np.array([0, 1.0, 0]))]:
            Fm = self.matrix(K, channel, e4) / (8 * np.pi * W)
            Bs = [1j * sd(e3), sd(qh) @ sd(np.cross(kh, e3)), 1j * sd(kh) * np.dot(qh, e3),
                  1j * sd(qh) * np.dot(qh, e3)]
            for i in range(2):
                for j in range(2):
                    A.append([B[i, j] for B in Bs])
                    b.append(Fm[i, j])
        A, b = np.array(A), np.array(b)
        F = np.linalg.lstsq(A, b, rcond=None)[0]
        resid = np.linalg.norm(A @ F - b) / max(np.linalg.norm(b), 1e-300)
        return F, K, resid

    def multipoles(self, W, channel, lmax=2, ngauss=40):
        xs, ws = np.polynomial.legendre.leggauss(ngauss)
        Fs = np.array([self.cgln(W, x, channel)[0] for x in xs])
        P = lambda L: (np.polynomial.legendre.legval(xs, [0] * L + [1]) if L >= 0 else 0 * xs)
        out = {}
        for l in range(0, lmax + 1):
            out['E%d+' % l] = np.sum(ws * (P(l) * Fs[:, 0] - P(l + 1) * Fs[:, 1]
                                           + l / (2 * l + 1.0) * (P(l - 1) - P(l + 1)) * Fs[:, 2]
                                           + (l + 1) / (2 * l + 3.0) * (P(l) - P(l + 2)) * Fs[:, 3])) / (2 * (l + 1))
            out['M%d+' % l] = 0.0 if l == 0 else np.sum(ws * (P(l) * Fs[:, 0] - P(l + 1) * Fs[:, 1]
                                           - 1.0 / (2 * l + 1) * (P(l - 1) - P(l + 1)) * Fs[:, 2])) / (2 * (l + 1))
            if l >= 1:
                out['E%d-' % l] = 0.0 if l == 1 else np.sum(ws * (P(l) * Fs[:, 0] - P(l - 1) * Fs[:, 1]
                                               - (l + 1) / (2 * l + 1.0) * (P(l - 1) - P(l + 1)) * Fs[:, 2]
                                               + l / (2 * l - 1.0) * (P(l) - P(l - 2)) * Fs[:, 3])) / (2 * l)
                out['M%d-' % l] = np.sum(ws * (-P(l) * Fs[:, 0] + P(l - 1) * Fs[:, 1]
                                               + 1.0 / (2 * l + 1.0) * (P(l - 1) - P(l + 1)) * Fs[:, 2])) / (2 * l)
        return out

    def mult32(self, W, lmax=1, ngauss=40):
        a = self.multipoles(W, 'pi0p', lmax, ngauss)
        b = self.multipoles(W, 'pi+n', lmax, ngauss)
        return {k: a[k] - b[k] / np.sqrt(2) for k in a}, a, b


def calibrated(**kw):
    """Born object with the graph signs fixed by gauge invariance and the overall sign by the LET
    (E0+(pi+ n) > 0)."""
    best = None
    for st in (1.0, -1.0):
        for skr in (1.0, -1.0):
            B = Born(st=st, skr=skr, **kw)
            gi = max(B.gauge_check(1.2, 0.3, 'pi+n'), B.gauge_check(1.25, -0.6, 'pi+n'))
            if best is None or gi < best[0]:
                best = (gi, st, skr)
    gi, st, skr = best
    B = Born(st=st, skr=skr, **kw)
    W = B.M + B.mpi + 1e-4
    e0 = B.multipoles(W, 'pi+n', lmax=0)['E0+']
    ph = np.conj(e0) / abs(e0)             # global phase convention (Feynman-rule factors of i)
    B.overall = complex(np.round(ph.real), np.round(ph.imag))
    assert abs(abs(B.overall) - 1) < 1e-12
    B.gauge_violation = gi
    return B


def mpi_unit(mpi=0.13957):
    return 1e-3 / mpi


def main():
    u = mpi_unit()
    full = calibrated(kp=1.793, kn=-1.913)
    print('gauge-invariance residual |k.J|/|e.J| = %.1e  (t-sign %+d, KR-sign %+d, overall phase %s)' %
          (full.gauge_violation, full.st, full.skr, full.overall))
    for ch, mp in [('pi+n', 0.13957), ('pi0p', 0.13498)]:
        B = calibrated(kp=1.793, kn=-1.913, mpi=mp)
        W = B.M + mp + 1e-5
        print('%-5s E0+(threshold, Born with kappa) = %8.3f [1e-3/mpi+]' % (ch, B.multipoles(W, ch, 0)['E0+'].real / u))
    B = calibrated(kp=1.793, kn=-1.913, mpi=0.13498)
    W = B.M + B.mpi + 2e-3
    mm = B.multipoles(W, 'pi0p', 1)
    K = B.kin(W, 0)
    q = K['q']
    P1 = 3 * mm['E1+'] + mm['M1+'] - mm['M1-']
    P2 = 3 * mm['E1+'] - mm['M1+'] + mm['M1-']
    P3 = 2 * mm['M1+'] + mm['M1-']
    print('pi0 p Born: P1/q=%7.3f P2/q=%7.3f P3/q=%7.3f [1e-3/mpi+^2]' %
          tuple(x.real / q / u * 0.13957 for x in (P1, P2, P3)))
    lo = calibrated(kp=0.0, kn=0.0)
    ps = calibrated(kp=0.0, kn=0.0, coupling='PS')
    m1 = lo.multipoles(1.2, 'pi+n', 2)
    m2 = ps.multipoles(1.2, 'pi+n', 2)
    print('PS vs PV at kappa=0 (pi+ n, W=1.2): max |diff| = %.1e of %.1e' %
          (max(abs(m1[k] - m2[k]) for k in m1), max(abs(m1[k]) for k in m1)))
    print('\nLO Born (kappa=0), isospin 3/2 (MAID normalisation) [1e-3/mpi+]:')
    for W in [1.10, 1.15, 1.20, 1.232, 1.26, 1.30]:
        m32, a, b = lo.mult32(W)
        print('  W=%5.3f  M1+(3/2)=%8.3f  E1+(3/2)=%8.3f   | pi0p M1+ %8.3f E1+ %7.3f  pi+n M1+ %8.3f E1+ %7.3f' %
              (W, m32['M1+'].real / u, m32['E1+'].real / u, a['M1+'].real / u, a['E1+'].real / u,
               b['M1+'].real / u, b['E1+'].real / u))


if __name__ == '__main__':
    main()
