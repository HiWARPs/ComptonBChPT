"""check_bf30.py -- compare the virtual correction of radcorr.py (Denner-Dittmaier helicity amplitudes) with the
closed form of Brown & Feynman, Phys. Rev. 85 (1952) 231, Eq. (30) [+ Eq. (29): U^(1) = P(kappa,tau) + P(tau,kappa)],
dsigma = dsigma_KN [1 - (alpha/pi) U^(1)/U],  U = 4(1/kappa+1/tau)^2 - 4(1/kappa+1/tau) - (kappa/tau + tau/kappa) (Eq. 10),
kappa = -2 omega_1/m, tau = 2 omega_2/m (lab energies of the incoming/outgoing photon, electron at rest),
4 sinh^2 y = -(kappa+tau), h(y) = y^-1 int_0^y u coth u du, G_0(kappa) = -2/kappa int_{1-kappa}^1 ln(1-u) du/u;
real part: ln kappa -> ln|kappa|, G_0 -> its real part (BF App. W).  The photon mass lambda is in units of m.
Terms antisymmetric in (kappa, tau) cancel in the sum and are not needed.
"""
import os, sys
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import radcorr as rc

mp.mp.dps = 30


def h(y):
    y = mp.mpf(y)
    return mp.quad(lambda u: u * mp.coth(u) if u != 0 else mp.mpf(1), [0, y]) / y


def G0(k):
    """Re G_0(kappa) = (2/kappa) [Re L(1-kappa) - L(1)],  L(x) = -Li2(x)."""
    k = mp.mpf(k)
    L1k = -rc._re_li2(1 - k)
    L1 = -mp.pi ** 2 / 6
    return 2 / k * (L1k - L1)


def U_of(k, t):
    return 4 * (1 / k + 1 / t) ** 2 - 4 * (1 / k + 1 / t) - (k / t + t / k)


def P(k, t, lnlam):
    """symmetric part of BF Eq. (30) (real part)."""
    U = U_of(k, t)
    y = mp.asinh(mp.sqrt(-(k + t) / 4))
    c2 = mp.coth(2 * y); s2 = mp.sinh(2 * y); ch2 = mp.cosh(2 * y)
    lnk = mp.log(abs(k))
    val = (1 - 2 * y * c2) * lnlam * U
    val += -2 * y * c2 * (2 * h(y) - h(2 * y)) * U
    val += (-4 * y * s2 / (k * t) * (2 - ch2) + 2 * y * mp.coth(y)) * h(y)
    val += lnk * (4 * y * c2 * (4 / (k * t) * mp.cosh(y) ** 2 + (k - 6) / (2 * t) / ch2 + 4 / k ** 2 - 1 / k - t / (2 * k) - k / t - 1)
                  + 3 * t / (2 * k ** 2) + 3 * t / (2 * k) + 3 / t + 1 - 7 / (k * t) + 8 / k - 8 / k ** 2
                  + (2 * k - t ** 2 - k ** 2 * t) / (2 * k ** 2 * t * (k - 1))
                  - (2 * k ** 2 + t) / (2 * t * (k - 1) ** 2))
    val += y ** 2 / mp.sinh(y) ** 2 * (2 / k - 7 * k / 4 - 3 * t ** 2 / (4 * k))
    val += -4 * y * mp.tanh(y) * (mp.mpf(1) / 2 - 1 / k) + 4 * (1 / k + 1 / t) ** 2
    val += -12 / k - 3 * k / (2 * t) - 2 * k / t ** 2 + (k / t + mp.mpf(1) / 2) / (k - 1)
    val += G0(k) * (k ** 2 / t + t / k ** 2 + k / t + k + t / 2 + 2 / k - 3 / t - 1)
    return val


def delta_bf30(omega, theta_deg, M, lam):
    kin = rc.kin(omega, theta_deg, M)
    k = -2 * kin['omega'] / kin['M']; t = 2 * kin['omp'] / kin['M']
    lnlam = mp.log(mp.mpf(lam) / kin['M'])
    U = U_of(k, t)
    U1 = P(k, t, lnlam) + P(t, k, lnlam)
    return float(-rc.ALPHA / mp.pi * U1 / U), float(U)


if __name__ == '__main__':
    print('Normalisation check: BF Eq. (10) U(kappa,tau) should equal the Klein-Nishina bracket (omega_2/omega_1 + omega_1/omega_2 - sin^2 theta):')
    for (w, th, M) in [(2.62, 45., rc.M_E), (100., 90., rc.M_P)]:
        kin = rc.kin(w, th, M); k = -2 * kin['omega'] / kin['M']; t = 2 * kin['omp'] / kin['M']
        print('   U = %.6f   KN bracket = %.6f' % (float(U_of(k, t)), float(kin['omp'] / kin['omega'] + kin['omega'] / kin['omp'] - mp.sin(kin['theta']) ** 2)))
    print()
    print('| M | omega | theta_lab | lambda | delta_V (DD, radcorr) | delta_V (BF Eq. 30) | difference |')
    print('|---|---|---|---|---|---|---|')
    for (M, w, lam) in [(rc.M_E, 2.62, rc.M_E), (rc.M_E, 17.6, rc.M_E), (rc.M_E, 17.6, 0.01), (rc.M_E, 150., rc.M_E),
                        (rc.M_E, 1000., rc.M_E), (rc.M_P, 100., rc.M_P), (rc.M_P, 150., 1.0), (rc.M_P, 150., rc.M_P)]:
        for th in [0.01, 20., 45., 70., 90., 120., 150., 179.99]:
            d_dd = rc.delta_virtual(w, th, M, lam=lam)
            d_bf, U = delta_bf30(w, th, M, lam)
            print('| %.3f | %g | %g | %g | %+.6e | %+.6e | %+.1e |' % (M, w, th, lam, d_dd, d_bf, d_dd - d_bf))
