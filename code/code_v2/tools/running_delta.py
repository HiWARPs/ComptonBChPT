#!/usr/bin/env python3
"""
Python prototype / cross-check of the 'drunning' option (LMP2015 Sect. 2.3, App. D) and of its
doubly-subtracted variant (&v2 irun=2, work/notes/delta_running2.md):
  * running Re gM(s), Re gE(s) from the dispersion relation of the imaginary parts of the
    gamma N Delta vertex loops (pv2006_vertex.py);
  * running Delta mass M_Delta(s) = M_Delta + Re sigma_R(s) and field renormalisation
    1 + Re tau_R(s) from the piN and piDelta self-energy loops (D1-D4), Xi = 0, scale M_N.

Subtraction schemes (nsub):
  1  once subtracted at s = M_D^2 (LMP2015 Eq. D5; sigma_R, tau_R of Eq. D4):
       Re g(s) = g + F1(s),  F1(s) = (s-M_D^2)/pi P int ds' Im g(s') / ((s'-M_D^2)(s'-s))
       sigma_R(s) = sigma(s) - Re sigma(M_D^2) - (s-M_D^2) Re sigma'(M_D^2),  tau_R(s) = tau(s) - Re tau(M_D^2)
  2  additionally subtracted at the static point s0 = M_N^2 (Re g(s0) = g, sigma_R(s0) = 0, tau_R(s0) = 0):
       Re g(s) = g + (s-M_D^2)(s-s0)/pi P int ds' Im g(s') / ((s'-s)(s'-M_D^2)(s'-s0))
               = g + F1(s) - (s-M_D^2)/(s0-M_D^2) F1(s0)                       [partial fractions]
       sigma_R(s) -> sigma_R(s) - c (s-M_D^2)^2,  c = sigma_R(s0)/(s0-M_D^2)^2
       tau_R(s)   -> tau_R(s) - c' (s-M_D^2),     c' = tau_R(s0)/(s0-M_D^2)
     (the twice-subtracted integral vanishes at both subtraction points, so the degree-1 polynomial of the
      general twice-subtracted representation is the constant g: Re g(M_D^2) = Re g(s0) = g).
  3  as 2, plus vanishing first derivative at s0 (Re g'(s0) = 0, sigma_R'(s0) = 0, tau_R'(s0) = 0):
       Re g(s) = g + F1(s) - (s-M_D^2) [a + b (s-s0)],  a = F1(s0)/(s0-M_D^2),  b = (F1'(s0) - a)/(s0-M_D^2),
       F1'(s0) = 1/pi int ds' Im g(s')/(s'-s0)^2 ;  sigma_R -> sigma_R - c (s-M_D^2)^2 - d (s-M_D^2)^2 (s-s0),
       d = sigma_R^(2)'(s0)/(s0-M_D^2)^2 ;  tau_R -> tau_R - c' (s-M_D^2) - d' (s-M_D^2)(s-s0), d' = tau_R^(2)'(s0)/(s0-M_D^2).
     The subtracted polynomials are analytic in s (counterterms of higher order).  Scheme 2 leaves alpha_E1,
     beta_M1 unchanged (they depend on the running functions only through their values at s0); the spin
     polarisabilities also depend on the slopes at s0, which scheme 3 removes.
Compared with LMP2015 Fig. 5 at read-off points (main()).  Units: GeV.
usage: running_delta.py [--table] [--check] [--fortran-compare FILE]
"""
import sys
import numpy as np
import mpmath as mp
import pv2006_vertex as PV

MN, MPI, MD, GA, HA, FPI = 0.93827, 0.13957, 1.232, 1.27, 2.85, 0.09221
S0 = MN ** 2
STH = (MN + MPI) ** 2


def im_g(W, offshell=True):
    if W <= MN + MPI + 1e-9:
        return 0.0, 0.0
    a = PV.im_gME(W, MN=MN, MD=MD, offshell=offshell, Q2t=mp.mpf('1e-20'), mpi=MPI, gA=GA, hA=HA, fpi=FPI)
    return a[0], a[1]


_NODES = {}


def nodes(offshell=True, c=1.0, n=400):
    """quadrature nodes s' = sth + c (y/(1-y))^2, y in (0,1) Gauss-Legendre (threshold behaviour ~ (s'-sth)^{3/2}),
    with the vertex Im parts tabulated on them (cached)."""
    key = (offshell, c, n)
    if key not in _NODES:
        y, wy = np.polynomial.legendre.leggauss(n)
        y = 0.5 * (y + 1); wy = 0.5 * wy
        u = y / (1 - y)
        sp = STH + c * u ** 2
        jac = c * 2 * u / (1 - y) ** 2
        img = np.array([im_g(np.sqrt(x), offshell) for x in sp])          # (n,2)
        imgD = np.array(im_g(MD, offshell))
        _NODES[key] = (sp, wy * jac, img, imgD)
    return _NODES[key]


def wfun(x, s_, c=1.0):
    """w_x(s') = 1 at s' = x, ~ 1/s' at infinity: P int_sth^inf ds' w_x(s')/(s'-x) = -log((x-sth)/c)"""
    return (x - STH + c) / (s_ - STH + c)


def F1(svals, offshell=True, c=1.0, n=400):
    """once-subtracted integral F1(s) = Re g(s) - g(M_D^2) = (s-M_D^2)/pi P int ds' Im g(s')/((s'-M_D^2)(s'-s));
    returns array (len(svals), 2) for (gM, gE)."""
    sp, w, img, imgD = nodes(offshell, c, n)
    MD2 = MD ** 2
    out = []
    for s in svals:
        ims = np.array(im_g(np.sqrt(s), offshell)) if s > STH else np.zeros(2)
        G = np.zeros((n, 2))
        for k in range(2):
            if s > STH:
                G[:, k] += (img[:, k] - ims[k] * wfun(s, sp, c)) / (sp - s)
            else:
                G[:, k] += img[:, k] / (sp - s)
            G[:, k] -= (img[:, k] - imgD[k] * wfun(MD2, sp, c)) / (sp - MD2)
        F = (G * w[:, None]).sum(0)
        if s > STH:
            F -= ims * np.log((s - STH) / c)
        F += imgD * np.log((MD2 - STH) / c)
        out.append(F / np.pi)
    return np.array(out)


def F1prime_s0(offshell=True, c=1.0, n=400):
    """dF1/ds at s0 (below threshold, no principal value): 1/pi int ds' Im g(s')/(s'-s0)^2"""
    sp, w, img, imgD = nodes(offshell, c, n)
    return (img * (w / (sp - S0) ** 2)[:, None]).sum(0) / np.pi


def re_g_running(svals, offshell=True, c=1.0, n=400, nsub=1):
    """Re g(s) - g for the subtraction scheme nsub (module doc); returns arrays dgM, dgE (shape (len(svals), 2))"""
    svals = np.asarray(svals, float)
    MD2 = MD ** 2
    out = F1(svals, offshell, c, n)
    if nsub == 1:
        return out
    f0 = F1([S0], offshell, c, n)[0]
    a = f0 / (S0 - MD2)
    if nsub == 2:
        return out - (svals - MD2)[:, None] * a[None, :]
    if nsub == 3:
        b = (F1prime_s0(offshell, c, n) - a) / (S0 - MD2)
        return out - (svals - MD2)[:, None] * (a[None, :] + (svals - S0)[:, None] * b[None, :])
    raise ValueError(nsub)


def re_g_running_direct2(svals, offshell=True, c=1.0, n=400):
    """twice-subtracted integral evaluated DIRECTLY (identity check of the partial-fraction form of nsub=2):
    F2(s) = (s-M_D^2)(s-s0)/pi P int ds' Im g(s') / ((s'-s)(s'-M_D^2)(s'-s0)),
    both principal values (s'=s and s'=M_D^2 lie on the cut for s > sth) by singularity subtraction with the
    weight w_x and the analytic integrals P int w_x/(s'-x) = -log((x-sth)/c)."""
    sp, w, img, imgD = nodes(offshell, c, n)
    MD2 = MD ** 2
    out = []
    for s in svals:
        ims = np.array(im_g(np.sqrt(s), offshell)) if s > STH else np.zeros(2)
        pref = (s - MD2) * (s - S0)
        rM = (s - S0) / (MD2 - S0)                   # residue factor at s' = M_D^2:  -rM Im(M_D^2)/(s'-M_D^2)
        G = pref * img / ((sp - s) * (sp - MD2) * (sp - S0))[:, None]
        if s > STH:
            G -= ims[None, :] * (wfun(s, sp, c) / (sp - s))[:, None]
        G += rM * imgD[None, :] * (wfun(MD2, sp, c) / (sp - MD2))[:, None]
        F = (G * w[:, None]).sum(0)
        if s > STH:
            F -= ims * np.log((s - STH) / c)
        F += rM * imgD * np.log((MD2 - STH) / c)
        out.append(F / np.pi)
    return np.array(out)


def sigma_tau(s, which='N', deriv=False, nx=96):
    """sigma(s), tau(s), dsigma/ds, dtau/ds (complex) of LMP2015 Eq. (D2), Xi=0, log scale mu = M_N."""
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
    sig = 0j; tau = 0j; dsig = 0j; dtau = 0j
    for x0, x1 in zip(pts[:-1], pts[1:]):
        if x1 - x0 < 1e-13:                       # coincident roots (exactly at threshold): 0 * log(0)
            continue
        # cubic endpoint map x = mid + half (3u - u^3)/2: the log singularities at the zeros of M^2(x) (interval
        # ends) converge like a smooth integrand (sigma'(M_D^2) to 1e-8 with 96 nodes; plain Gauss: 2e-5)
        half = 0.5 * (x1 - x0); mid = 0.5 * (x1 + x0)
        x = mid + half * (3 * xs - xs ** 3) / 2; w = half * 1.5 * (1 - xs ** 2) * ws
        M2 = MPI ** 2 * x + Mx ** 2 * (1 - x) - s * x * (1 - x)
        L = np.log(np.abs(M2) / MN ** 2) - 1j * np.pi * (M2 < 0)
        sig += np.sum(w * (x * MD + Mx) * M2 * (-1 + L))
        tau += np.sum(w * x * M2 * (-1 + L))
        dsig += np.sum(w * (x * MD + Mx) * (-x * (1 - x)) * L)      # d/ds [M2 (-1+L)] = dM2/ds L, dM2/ds = -x(1-x)
        dtau += np.sum(w * x * (-x * (1 - x)) * L)
    return pref * sig, pref * tau, pref * dsig, pref * dtau


def sigma_tau_R1(s):
    """once-subtracted (LMP2015 D4) Re sigma_R(s), Re tau_R(s) and their s-derivatives (sum of piN and piDelta)"""
    MD2 = MD ** 2
    sR = tR = dsR = dtR = 0.0
    for k in ('N', 'D'):
        sg, tg, dsg, dtg = sigma_tau(s, k)
        s0, t0, d0, _ = sigma_tau(MD2, k)
        sR += sg.real - s0.real - (s - MD2) * d0.real
        tR += tg.real - t0.real
        dsR += dsg.real - d0.real
        dtR += dtg.real
    return sR, tR, dsR, dtR


def running_mass_tau(svals, nsub=1):
    """[M_D(s), 1 + Re tau_R(s), Gamma(s)] for the subtraction scheme nsub"""
    MD2 = MD ** 2
    sR0, tR0, dsR0, dtR0 = sigma_tau_R1(S0)
    D = S0 - MD2
    c = sR0 / D ** 2; cp = tR0 / D                        # scheme 2: sigma_R(s0) = 0, tau_R(s0) = 0
    d = (dsR0 - 2 * c * D) / D ** 2; dp = (dtR0 - cp) / D  # scheme 3: + sigma_R'(s0) = 0, tau_R'(s0) = 0
    out = []
    for s in svals:
        sR, tR, _, _ = sigma_tau_R1(s)
        gam = -2 * sum(sigma_tau(s, k)[0].imag for k in ('N', 'D'))
        if nsub >= 2:
            sR -= c * (s - MD2) ** 2; tR -= cp * (s - MD2)
        if nsub == 3:
            sR -= d * (s - MD2) ** 2 * (s - S0); tR -= dp * (s - MD2) * (s - S0)
        out.append((MD + sR, 1 + tR, gam))
    return np.array(out)


# ----------------------------------------------------------------------------------------------------------
WGRID = [MN, MN + MPI, 1.1, 1.15, 1.2, 1.232, 1.3, 1.4]


def table(nsub, offshell=False, Ws=WGRID):
    """Re gM(s), Re gE(s), M_D(s)/M_N, 1 + Re tau_R(s) on the grid; nsub = 0 means no running"""
    s = np.array(Ws) ** 2
    if nsub == 0:
        return np.array([[W, 2.97, -1.0, MD / MN, 1.0] for W in Ws])
    dg = re_g_running(s, offshell=offshell, nsub=nsub)
    R = running_mass_tau(s, nsub=nsub)
    return np.array([[W, 2.97 + d[0], -1.0 + d[1], r[0] / MN, r[1]] for W, d, r in zip(Ws, dg, R)])


def checks(offshell=False):
    """consistency checks of scheme 2: partial-fraction form vs direct twice-subtracted integral; values at the
    subtraction points; once-subtracted limit; mass slope at the pole."""
    MD2 = MD ** 2
    s = np.array([S0, 1.05 ** 2, 1.15 ** 2, MD2, 1.3 ** 2, 1.6 ** 2])
    A = re_g_running(s, offshell=offshell, nsub=2)
    B = re_g_running_direct2(s, offshell=offshell)
    print('scheme 2: |partial-fraction form - direct twice-subtracted integral| max = %.2e' % np.max(np.abs(A - B)))
    print('scheme 2: Re g - g at s0 = %s, at M_D^2 = %s' % (A[0], A[3]))
    F = re_g_running(s, offshell=offshell, nsub=1)
    print('scheme 1: Re g - g at s0 = %s, at M_D^2 = %s' % (F[0], F[3]))
    print('scheme 2 - scheme 1 is linear in s: %s' % ((A - F) / (s - MD2)[:, None]))
    for ns in (1, 2, 3):
        R = running_mass_tau([S0, MD2, MD2 + 1e-4, MD2 - 1e-4], nsub=ns)
        slope = (R[2, 0] - R[3, 0]) / 2e-4
        print('scheme %d: M_D(s0)-M_D = %+.6f, M_D(M_D^2)-M_D = %+.2e, dM_D/ds at the pole = %+.2e, '
              '1+tau_R(s0) = %.6f, 1+tau_R(M_D^2) = %.6f' % (ns, R[0, 0] - MD, R[1, 0] - MD, slope, R[0, 1], R[1, 1]))
    h = 1e-4
    for ns in (2, 3):
        g = re_g_running([S0 - h, S0 + h], offshell=offshell, nsub=ns)
        R = running_mass_tau([S0 - h, S0 + h], nsub=ns)
        print('scheme %d: slopes at s0: dgM/ds = %+.5f, dgE/ds = %+.5f, dM_D/ds = %+.5f, dtau_R/ds = %+.5f [GeV^-2]'
              % (ns, (g[1, 0] - g[0, 0]) / (2 * h), (g[1, 1] - g[0, 1]) / (2 * h), (R[1, 0] - R[0, 0]) / (2 * h),
                 (R[1, 1] - R[0, 1]) / (2 * h)))


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
    args = sys.argv[1:]
    if not args:
        main()
    if '--check' in args:
        checks(offshell=False)
    if '--table' in args:
        for ns in (0, 1, 2, 3):
            print('scheme nsub=%d (on-shell mapping of the PV2006 Im parts = Fortran vxmode=0):' % ns)
            print('   W[GeV]    Re gM     Re gE    M_D(s)/M_N  1+Re tau_R')
            for row in table(ns):
                print('  %7.5f %9.5f %9.5f %11.6f %11.6f' % tuple(row))
    if '--fortran-compare' in args:
        # file written by tests (program tr.f90): columns W[M units], nsub, ReG1, ReG2, Mrun, ztau  (vxmode=0)
        fn = args[args.index('--fortran-compare') + 1]
        D = np.loadtxt(fn)
        worst = 0.0
        for row in D:
            W, ns = row[0] * MN, int(round(row[1]))
            t = table(ns, Ws=[W])[0]
            dev = np.abs(np.array([t[1], t[2], t[3], t[4]]) - row[2:6])
            worst = max(worst, dev.max())
            print('  W=%.5f nsub=%d  |Fortran - Python| gM %.1e gE %.1e M_D %.1e ztau %.1e' % (W, ns, *dev))
        print('FORTRAN-PYTHON max |diff| = %.2e' % worst)
