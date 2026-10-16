#!/usr/bin/env python3
"""
Python prototype / cross-check of the 'drunning' option (LMP2015 Sect. 2.3, App. D):
  * running Re gM(s), Re gE(s) from the once-subtracted dispersion relation (D5) of the
    imaginary parts of the gamma N Delta vertex loops (pv2006_vertex.py);
  * running Delta mass M_Delta(s) = M_Delta + Re sigma_R(s) and field renormalisation
    1 + Re tau_R(s) from the piN and piDelta self-energy loops (D1-D4), Xi = 0, scale M_N.
Compared with LMP2015 Fig. 5 at read-off points (main()).  Units: GeV.
"""
import numpy as np
import mpmath as mp
import pv2006_vertex as PV

MN, MPI, MD, GA, HA, FPI = 0.93827, 0.13957, 1.232, 1.27, 2.85, 0.09221


def im_g(W, offshell=True):
    if W <= MN + MPI + 1e-9:
        return 0.0, 0.0
    a = PV.im_gME(W, MN=MN, MD=MD, offshell=offshell, Q2t=mp.mpf('1e-20'), mpi=MPI, gA=GA, hA=HA, fpi=FPI)
    return a[0], a[1]


def re_g_running(svals, offshell=True, c=1.0, n=400):
    """Re g(s) - g(MD^2) = (s-MD^2)/pi P int ds' Im g(s')/((s'-MD^2)(s'-s)); returns arrays dgM, dgE"""
    sth = (MN + MPI) ** 2
    y, wy = np.polynomial.legendre.leggauss(n)
    y = 0.5 * (y + 1); wy = 0.5 * wy
    # map y in (0,1) -> s' = sth + c*(y/(1-y))^2  (threshold behaviour ~ (s'-sth)^{3/2})
    u = y / (1 - y)
    sp = sth + c * u ** 2
    jac = c * 2 * u / (1 - y) ** 2
    img = np.array([im_g(np.sqrt(x), offshell) for x in sp])          # (n,2)
    MD2 = MD ** 2
    imgD = np.array(im_g(MD, offshell))

    def wfun(x, s_):            # w_x(s') = 1 at s'=x, ~1/s' at infinity
        return (x - sth + c) / (s_ - sth + c)
    out = []
    for s in svals:
        ims = np.array(im_g(np.sqrt(s), offshell)) if s > sth else np.zeros(2)
        G = np.zeros((n, 2))
        for k in range(2):
            if s > sth:
                G[:, k] += (img[:, k] - ims[k] * wfun(s, sp)) / (sp - s)
            else:
                G[:, k] += img[:, k] / (sp - s)
            G[:, k] -= (img[:, k] - imgD[k] * wfun(MD2, sp)) / (sp - MD2)
        F = (G * (wy * jac)[:, None]).sum(0)
        if s > sth:
            F -= ims * np.log((s - sth) / c)
        F += imgD * np.log((MD2 - sth) / c)
        out.append(F / np.pi)
    return np.array(out)


def sigma_tau(s, which='N', deriv=False, nx=200):
    """sigma(s), tau(s) (complex) of LMP2015 Eq. (D2), Xi=0, log scale mu = M_N."""
    if which == 'N':
        C, Mx = HA ** 2, MN
    else:
        C, Mx = (5.0 / 3 * 9.0 / 5 * GA) ** 2, MD
    pref = -C / (2 * (8 * np.pi * FPI) ** 2)
    # integrand in x, split at roots of M^2(x) = s x^2 - (s - mpi^2 + Mx^2) x + Mx^2 ... (M^2 = mpi^2 x + Mx^2 (1-x) - s x (1-x))
    a, b, cc = s, -(s - MPI ** 2 + Mx ** 2), Mx ** 2
    disc = b * b - 4 * a * cc
    pts = [0.0, 1.0]
    if disc > 0:
        for r in ((-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a)):
            if 0 < r < 1:
                pts.append(r)
    pts = sorted(pts)
    xs, ws = np.polynomial.legendre.leggauss(nx)
    sig = 0j; tau = 0j; dsig = 0j
    for x0, x1 in zip(pts[:-1], pts[1:]):
        x = 0.5 * (x1 - x0) * xs + 0.5 * (x1 + x0); w = 0.5 * (x1 - x0) * ws
        M2 = MPI ** 2 * x + Mx ** 2 * (1 - x) - s * x * (1 - x)
        L = np.log(np.abs(M2) / MN ** 2) - 1j * np.pi * (M2 < 0)
        sig += np.sum(w * (x * MD + Mx) * M2 * (-1 + L))
        tau += np.sum(w * x * M2 * (-1 + L))
        dsig += np.sum(w * (x * MD + Mx) * (-x * (1 - x)) * L)
    return pref * sig, pref * tau, pref * dsig


def running_mass_tau(svals):
    out = []
    MD2 = MD ** 2
    sD = {k: sigma_tau(MD2, k) for k in ('N', 'D')}
    for s in svals:
        sR = 0j; tR = 0j
        for k in ('N', 'D'):
            sg, tg, _ = sigma_tau(s, k)
            s0, t0, d0 = sD[k]
            sR += sg - s0.real - (s - MD2) * d0.real
            tR += tg - t0.real
        out.append((MD + sR.real, 1 + tR.real, -2 * sR.imag))
    return np.array(out)


def main():
    Ws = np.array([0.94, 1.0, 1.05, 1.1, 1.15, 1.2, 1.232, 1.26, 1.29])
    s = Ws ** 2
    for offs in (True, False):
        dg = re_g_running(s, offshell=offs)
        print('mapping M_D->W' if offs else 'on-shell mapping')
        for W, d in zip(Ws, dg):
            print('  W=%5.3f  gM=%7.4f  gE=%7.4f' % (W, 2.97 + d[0], -1.0 + d[1]))
    R = running_mass_tau(s)
    print('running mass, 1+Re tau_R, Gamma(App. D):')
    for W, r in zip(Ws, R):
        print('  W=%5.3f  M_D(s)/M_N=%.4f  1+Re tau_R=%.4f  Gamma=%.4f GeV' % (W, r[0] / MN, r[1], r[2]))


if __name__ == '__main__':
    main()
