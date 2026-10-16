"""
radcorr.py -- O(alpha) QED (radiative) corrections to real Compton scattering on a point-like
spin-1/2 particle of mass M and charge e (on-shell renormalisation: physical mass and charge).

    dsigma = dsigma_Born * [1 + delta(omega, theta; Delta E)],
    delta  = delta_V(lambda) + delta_soft(lambda, Delta E),

* delta_V   : complete one-loop virtual correction, from the six independent one-loop helicity amplitudes
              of Denner & Dittmaier, Nucl. Phys. B 540 (1999) 58 [hep-ph/9805443], Eqs. (3.3)-(3.5)
              (exact in M, s, t; photon mass lambda as IR regulator).  Unpolarised:
              delta_V = sum_hel 2 Re(M0 M1*) / sum_hel |M0|^2.
* delta_soft: real emission of one photon with energy below Delta E in a chosen frame (default: LAB,
              target at rest), eikonal (soft-photon) approximation, integrated exactly with the photon-mass
              regulator ('t Hooft-Veltman integrals as given by Denner, Fortsch. Phys. 41 (1993) 307,
              Sect. 7.2 / arXiv:0709.1075).  Also available: brute-force 2d numerical integral (check) and
              the CM-frame closed form of Denner-Dittmaier Eq. (4.1) (check).
* The sum is independent of lambda (checked numerically).
* Linear photon polarisation: dsigma_par / dsigma_perp (E-vector in / perpendicular to the scattering
  plane) and the beam asymmetry Sigma_3 are built from the helicity amplitudes; the soft factor cancels
  in Sigma_3, so its correction is purely virtual.

Everything is evaluated with mpmath at high precision (the DD amplitudes contain 1/r^4, r = s_m cos(theta/2),
with cancellations near theta = 180 deg).

Units: MeV.  Angles: degrees.  Lab frame: target at rest, omega = incident photon energy.
Author: generated for the Compton-BChPT project, 2026-10.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 40

ALPHA = 1.0 / 137.035999084
M_P = 938.272088
M_E = 0.51099895
PI = mp.pi


# ----------------------------------------------------------------------------------------------------------
# kinematics
# ----------------------------------------------------------------------------------------------------------
def kin(omega, theta_deg, M=M_P):
    """lab kinematics -> Mandelstam variables and CM quantities (mpmath numbers)."""
    omega = mp.mpf(omega); M = mp.mpf(M)
    th = mp.mpf(theta_deg) * PI / 180
    c = mp.cos(th)
    omp = omega / (1 + omega / M * (1 - c))            # scattered photon energy (lab)
    s = M ** 2 + 2 * M * omega
    t = -2 * omega * omp * (1 - c)
    u = 2 * M ** 2 - s - t
    rs = mp.sqrt(s)
    Eg = (s - M ** 2) / (2 * rs)                        # CM photon energy
    Ep = (s + M ** 2) / (2 * rs)                        # CM target energy
    beta = Eg / Ep                                      # DD's beta
    cth = 1 + t / (2 * Eg ** 2)                         # CM cos(theta)
    # recoil in the lab
    Ep_lab = M + omega - omp
    p_lab = mp.sqrt(Ep_lab ** 2 - M ** 2)
    return dict(omega=omega, omp=omp, theta=th, s=s, t=t, u=u, M=M, Eg=Eg, Ep=Ep, beta=beta, cos_cm=cth,
                Ep_lab=Ep_lab, p_lab=p_lab, cos_lab=c)


def theta_cm_of_lab(omega, theta_lab_deg, M=M_P):
    k = kin(omega, theta_lab_deg, M)
    return float(mp.acos(k['cos_cm']) * 180 / PI)


def theta_lab_of_cm(omega, theta_cm_deg, M=M_P):
    """exact 2-body: cos_cm = (cos_L - b)/(1 - b cos_L), b = omega/(omega+M)  ->  invert."""
    b = omega / (omega + M)
    cc = np.cos(np.radians(theta_cm_deg))
    cl = (cc + b) / (1 + b * cc)
    return float(np.degrees(np.arccos(np.clip(cl, -1, 1))))


# ----------------------------------------------------------------------------------------------------------
# scalar one-loop integrals (DD Eq. (3.4)); only real parts are needed (Born amplitudes are real)
# ----------------------------------------------------------------------------------------------------------
def _li2(z):
    return mp.polylog(2, z)


def _re_li2(x):
    """Re Li2(x) for real x (continuous across the cut x > 1)."""
    x = mp.mpf(x)
    if x <= 1:
        return _li2(x)
    return PI ** 2 / 3 - mp.log(x) ** 2 / 2 - _li2(1 / x)


def scalar_integrals(s, t, u, M, lam):
    """returns dict of REAL parts of the DD scalar integrals, physical region s > M^2, t < 0."""
    m2 = M ** 2
    sm, um = s - m2, u - m2
    bt = mp.sqrt(1 - 4 * m2 / t)                       # t < 0 -> real, > 1
    xt = (bt - 1) / (bt + 1)                           # 0 < xt < 1
    lxt = mp.log(xt)
    L2 = mp.log(lam ** 2 / m2)
    out = {}
    out['Bp'] = -(1 + mp.log(lam / M)) / m2
    out['Bmm_t'] = bt * lxt
    # B_0m(v) = -(v_m/v) ln(-(v_m + i eps)/M^2): real part uses ln|v_m|
    out['B0m_s'] = -(sm / s) * mp.log(sm / m2)
    out['B0m_u'] = -(um / u) * mp.log(-um / m2)
    out['Cmmm_t'] = lxt ** 2 / (2 * t)
    # C_m0m(t): Re[2 Li2(1+xt) - 2 pi i ln(1+xt)] = 2 Re Li2(1+xt)
    out['Cm0m_t'] = (L2 * lxt + lxt ** 2 / 2 + 2 * _re_li2(1 + xt) - PI ** 2 / 2) / (bt * t)
    # C_0mm(v) = (1/v_m)[Li2(-(v_m+i eps)/M^2) + ln((v+i eps)/M^2) ln(-(v_m+i eps)/M^2)]
    #   v = s: Li2(-sm/M^2) real; ln(s/M^2) real; ln(-sm/M^2 - i eps) = ln(sm/M^2) - i pi -> Re = product of logs
    out['C0mm_s'] = (_re_li2(-sm / m2) + mp.log(s / m2) * mp.log(sm / m2)) / sm
    #   v = u: -um/M^2 > 0 (Li2 may need x > 1 branch -> real part); ln(u/M^2): u may be < 0 (electron, high E)
    lu = mp.log(abs(u) / m2)
    out['C0mm_u'] = (_re_li2(-um / m2) + lu * mp.log(-um / m2)) / um
    # D_m0mm(t, v) = 1/(bt t v_m) [2 ln xt ln(lam M/(-v_m - i eps)) - 2 Li2(1-xt) + 2 Li2(1+xt) - pi^2/2 - 2 pi i ln(1+xt)]
    common = -2 * _li2(1 - xt) + 2 * _re_li2(1 + xt) - PI ** 2 / 2
    out['D_ts'] = (2 * lxt * mp.log(lam * M / sm) + common) / (bt * t * sm)
    out['D_tu'] = (2 * lxt * mp.log(lam * M / (-um)) + common) / (bt * t * um)
    return out


# ----------------------------------------------------------------------------------------------------------
# helicity amplitudes (DD Eqs. (2.9) and (3.3)), upper-sign members:  A1=M(++++), A2=M(+-+-), A3=M(+--+),
# A4=M(+---)=M(++-+), A5=M(+++-)=M(+-++), A6=M(++--).  Lower-sign members: same for A1,A2,A5; opposite sign
# for A3,A4,A6 (both Born and one-loop; irrelevant for the unpolarised and linear-polarisation sums used here).
# ----------------------------------------------------------------------------------------------------------
def born_amps(s, t, u, M):
    m2 = M ** 2
    sm, um = s - m2, u - m2
    r = mp.sqrt(m2 ** 2 - s * u)
    mt = mp.sqrt(-t)
    pref = 8 * PI * ALPHA / (sm ** 2 * um)
    A1 = pref * r * (sm ** 2 + m2 * t)
    A2 = pref * r ** 3
    A3 = pref * M * s * t * mt
    A4 = pref * M * r ** 2 * mt
    A5 = pref * m2 * r * t
    A6 = pref * M ** 3 * t * mt
    return [A1, A2, A3, A4, A5, A6]


def loop_amps(s, t, u, M, lam):
    """real parts of the one-loop amplitudes M1 (DD Eq. (3.3)); photon mass lam (IR regulator)."""
    m2 = M ** 2
    m4 = m2 ** 2
    m6 = m2 ** 3
    sm, um, tm = s - m2, u - m2, t - 4 * m2
    r2 = m4 - s * u
    r = mp.sqrt(r2)
    r4 = r2 ** 2
    mt = mp.sqrt(-t)
    I = scalar_integrals(s, t, u, M, lam)
    Bp, Bmm, B0s, B0u = I['Bp'], I['Bmm_t'], I['B0m_s'], I['B0m_u']
    Cmmm, Cm0m, C0s, C0u = I['Cmmm_t'], I['Cm0m_t'], I['C0mm_s'], I['C0mm_u']
    Dts, Dtu = I['D_ts'], I['D_tu']
    a2 = ALPHA ** 2

    # M1(++++)
    B1 = 2 * a2 * r / sm * (
        (sm ** 2 * u + m2 * t ** 2 - 2 * m4 * t) / (s * um * u)
        + 4 * m2 * (sm ** 2 + m2 * t) / (sm * um) * Bp
        + 2 * t / tm * Bmm
        + m2 * (5 * s ** 2 - 2 * m2 * s + m4) / (s * sm ** 2) * B0s
        + (r2 * u + 2 * t * u ** 2 + 4 * m2 * u ** 2 + m4 * t) / (u * um ** 2) * B0u
        + (1 / r2) * (
            - t * (s ** 2 + t ** 2 + m2 * t - 9 * m4) * Cmmm
            + tm * (sm ** 2 - t ** 2 - m2 * t) * Cm0m
            - 2 * sm * (s ** 2 - m2 * (4 * s + u)) * C0s
            - 2 * um * (t ** 2 + 3 * m2 * s - 7 * m4) * C0u
            + (sm ** 3 * t - m2 * sm * (2 * sm ** 2 + sm * tm - t * tm) + 2 * m4 * t ** 2 - 4 * m6 * t) * Dts
            + (t ** 3 * um + m2 * (2 * t ** 3 + sm ** 2 * um - sm * um ** 2) - 2 * m4 * (sm ** 2 - 3 * um ** 2)
               - 4 * m6 * t) * Dtu))

    # M1(+-+-)
    B2 = 2 * a2 * r / sm * (
        - um / u
        + 4 * m2 * r2 / (sm * um) * Bp
        + 2 * t * (sm ** 2 + m2 * (3 * sm - um)) / (r2 * tm) * Bmm
        - (s + m2) * (sm ** 2 * (2 * sm + 3 * u) + 5 * m2 * sm * um - 2 * m4 * t) / (r2 * sm ** 2) * B0s
        - m2 * (3 * sm * u * um ** 2 - 2 * um ** 4 - m2 * (4 * r2 * u + um ** 3)) / (r2 * u * um ** 2) * B0u
        + (1 / r4) * (
            - t * (2 * r2 * s * tm + s ** 2 * sm ** 2 + m2 * (sm ** 3 - s * sm * um + sm * u * um - u * um ** 2)) * Cmmm
            + tm * (sm ** 3 * (sm + 2 * t) - m2 * t ** 2 * (4 * sm + t) - m4 * t ** 2) * Cm0m
            - 2 * sm * (s ** 2 * t * tm - m2 * s * sm * (3 * sm + 2 * t) + m4 * sm * um) * C0s
            + 2 * um * (r2 * s * u + m2 * u * (sm * t + um ** 2) + m4 * (2 * sm ** 2 + 3 * s * t)) * C0u
            + (s + m2) * (s ** 2 * t ** 3 + m2 * sm ** 2 * (2 * sm ** 2 + 3 * sm * t - 2 * t ** 2) - 5 * m4 * sm * t ** 2
                          - 2 * m6 * t ** 2) * Dts
            + (sm ** 2 * t * um ** 3 + m2 * um * (3 * t ** 3 * u + t ** 2 * um ** 2 - 5 * t * um ** 3 - 2 * um ** 4)
               - m4 * um ** 2 * (5 * t ** 2 + 6 * t * um - 4 * um ** 2) + 2 * m6 * t ** 2 * (t - 2 * u)) * Dtu))

    # M1(+--+)
    B3 = 4 * a2 * M * mt / sm * (
        sm / um
        + 2 * m2 * s * t / (sm * um) * Bp
        + s * (3 * s - m2) / sm ** 2 * B0s
        + (3 * s * u + m2 * s - 2 * m4) / um ** 2 * B0u
        - 2 * m2 * Cmmm
        - 2 * m2 * t / um * C0u
        - m2 * t * Dtu
        + (s * t - 3 * m2 * s + m4) * (Dts + Dtu))

    # M1(+---)
    B4 = 2 * a2 * M * mt / sm * (
        - r2 / (u * um)
        + 4 * m2 * r2 / (sm * um) * Bp
        + 2 * (s - u) / tm * Bmm
        + 4 * s ** 2 / sm ** 2 * B0s
        - (r2 * um + 2 * u * um ** 2 - 4 * m4 * u) / (u * um ** 2) * B0u
        + (1 / r2) * (
            - t * (2 * s * t + u * um - 7 * m2 * s - m4) * Cmmm
            + tm * (2 * sm ** 2 - um ** 2 - m2 * t) * Cm0m
            - 2 * sm * (r2 - 2 * s ** 2 - 2 * m2 * s) * C0s
            + 2 * um * (2 * s + t) * (u + m2) * C0u
            - (sm ** 3 * t - s ** 2 * t * tm - 4 * m2 * s * sm ** 2 - m2 * s * t ** 2) * Dts
            + (t ** 2 * u ** 2 + 2 * tm * um ** 3 + 4 * m2 * u * um ** 2 + m4 * t * tm) * Dtu))

    # M1(+++-)
    B5 = 2 * a2 * m2 * r / (sm * um) * (
        (s - u) / m2
        + 4 * m2 * t / sm * Bp
        + 4 * s * um / sm ** 2 * B0s
        + 2 * (2 * u - t) / um * B0u
        + (1 / r2) * (
            - 2 * um * (sm * (s - u) + 2 * m2 * t) * Cmmm
            - 2 * s * t * um * C0s
            - 2 * t * (r2 - m2 * um) * C0u
            - um * (s * t * (s + 2 * u) + 4 * m2 * r2 - 3 * m4 * t) * (Dts + Dtu)
            - t ** 2 * um * (s + m2) * Dtu))

    # M1(++--)
    B6 = 4 * a2 * M ** 3 * mt / sm ** 2 * (
        - sm / m2
        + 2 * m2 * t / um * Bp
        + (s + m2) / sm * B0s
        + sm * (3 * u - m2) / um ** 2 * B0u
        - 2 * sm * (Cmmm + C0s - C0u)
        - sm * (2 * s + um) * Dts
        - sm * (s + m2) * Dtu)
    return [B1, B2, B3, B4, B5, B6]


def hel_sums(A, B=None):
    """unpolarised sum D = sum_hel |M0|^2, linear-polarisation interference N = sum 2Re(M_+ M_-*)
    (initial photon helicity +/-), and their O(alpha) counterparts D1, N1 (B = Re M1)."""
    A1, A2, A3, A4, A5, A6 = A
    D0 = 2 * (A1 ** 2 + A2 ** 2 + A3 ** 2 + 2 * A4 ** 2 + 2 * A5 ** 2 + A6 ** 2)
    N0 = 4 * (A5 * (A1 + A2) + A4 * (A3 + A6))
    if B is None:
        return D0, N0
    B1, B2, B3, B4, B5, B6 = B
    D1 = 4 * (A1 * B1 + A2 * B2 + A3 * B3 + 2 * A4 * B4 + 2 * A5 * B5 + A6 * B6)
    N1 = 4 * (B5 * (A1 + A2) + A5 * (B1 + B2) + B4 * (A3 + A6) + A4 * (B3 + B6))
    return D0, N0, D1, N1


# ----------------------------------------------------------------------------------------------------------
# Born cross sections (checks)
# ----------------------------------------------------------------------------------------------------------
def klein_nishina_lab(omega, theta_deg, M=M_P, phi_deg=None):
    """dsigma/dOmega_lab [nb/sr] for a point Dirac particle. phi_deg: angle between the linear polarisation
    vector and the scattering plane (None = unpolarised)."""
    k = kin(omega, theta_deg, M)
    om, omp, c = k['omega'], k['omp'], k['cos_lab']
    hbarc2 = mp.mpf('0.389379372') * 1e6 * 1e3       # (hbar c)^2 in MeV^2 nb (0.3894 mb GeV^2)
    s2 = 1 - c ** 2
    if phi_deg is None:
        x = (omp / om + om / omp - s2) / 2
    else:
        x = (omp / om + om / omp - 2 * s2 * mp.cos(mp.radians(phi_deg)) ** 2) / 2
    return float(ALPHA ** 2 / M ** 2 * (omp / om) ** 2 * x * hbarc2)


def born_cs_lab(omega, theta_deg, M=M_P):
    """dsigma/dOmega_lab [nb/sr] from the helicity amplitudes (check of the transcription)."""
    k = kin(omega, theta_deg, M)
    A = born_amps(k['s'], k['t'], k['u'], k['M'])
    D0, N0 = hel_sums(A)
    dcm = D0 / 4 / (64 * PI ** 2 * k['s'])          # dsigma/dOmega_cm in MeV^-2
    # Jacobian dOmega_cm/dOmega_lab = (omp/omega)^2 * ... for a massless particle: dOmega_cm/dOmega_lab = (omega'_cm/omega'_lab)^2 ... use
    # dsigma/dOmega_lab = dsigma/dOmega_cm * (omp_lab/ Eg_cm)^2  (standard for photon scattering off a target at rest)
    dlab = dcm * (k['omp'] / k['Eg']) ** 2
    hbarc2 = mp.mpf('0.389379372') * 1e6 * 1e3
    return float(dlab * hbarc2)


# ----------------------------------------------------------------------------------------------------------
# virtual correction
# ----------------------------------------------------------------------------------------------------------
def virtual_parts(omega, theta_deg, M=M_P, lam=None):
    """returns (D0, N0, D1, N1) at lab kinematics; lam = photon mass (default M, i.e. ln(lam/M)=0)."""
    lam = mp.mpf(M) if lam is None else mp.mpf(lam)
    k = kin(omega, theta_deg, M)
    A = born_amps(k['s'], k['t'], k['u'], k['M'])
    B = loop_amps(k['s'], k['t'], k['u'], k['M'], lam)
    return hel_sums(A, B)


def delta_virtual(omega, theta_deg, M=M_P, lam=None):
    D0, N0, D1, N1 = virtual_parts(omega, theta_deg, M, lam)
    return float(D1 / D0)


# ----------------------------------------------------------------------------------------------------------
# soft-photon factor
# ----------------------------------------------------------------------------------------------------------
def _I_ii(p0, pabs, dE, lam):
    """Denner (7.13): I_ii = 2 pi { ln(4 dE^2/lam^2) + (p0/|p|) ln((p0-|p|)/(p0+|p|)) }"""
    if pabs == 0:
        return 2 * PI * (mp.log(4 * dE ** 2 / lam ** 2) - 2)
    return 2 * PI * (mp.log(4 * dE ** 2 / lam ** 2) + p0 / pabs * mp.log((p0 - pabs) / (p0 + pabs)))


def _I_ij(pi, pj, dE, lam):
    """'t Hooft-Veltman integral I_ij = int_{|k|<dE} d^3k/(2 w_k) 2 pi.pj/((pi.k)(pj.k)), Denner Eq. (7.12).
    pi, pj: 4-vectors (mp.matrix or lists) [E, px, py, pz]; both massive."""
    def dot(a, b):
        return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
    pi2, pj2, pij = dot(pi, pi), dot(pj, pj), dot(pi, pj)
    # alpha from alpha^2 pi^2 - 2 alpha pi.pj + pj^2 = 0 with (alpha pi0 - pj0)/pj0 > 0
    disc = mp.sqrt(pij ** 2 - pi2 * pj2)
    roots = [(pij + disc) / pi2, (pij - disc) / pi2]
    al = None
    for a in roots:
        if (a * pi[0] - pj[0]) / pj[0] > 0:
            al = a
    if al is None:
        raise ValueError('no valid alpha root')
    api = [al * x for x in pi]
    api2 = al ** 2 * pi2
    v = (api2 - pj2) / (2 * (api[0] - pj[0]))

    def F(u):
        u0 = u[0]
        uabs = mp.sqrt(u[1] ** 2 + u[2] ** 2 + u[3] ** 2)
        lg = mp.log((u0 - uabs) / (u0 + uabs))
        return (lg ** 2 / 4 + _re_li2(1 - (u0 + uabs) / v) + _re_li2(1 - (u0 - uabs) / v))
    pref = 4 * PI * al * pij / (api2 - pj2)
    return pref * (mp.log(api2 / pj2) / 2 * mp.log(4 * dE ** 2 / lam ** 2) + F(api) - F(pj))


def soft_lab(omega, theta_deg, dE, M=M_P, lam=None, method='analytic'):
    """delta_soft for a soft-photon cut |k| < dE in the LAB frame (target at rest).
    delta_soft = -(alpha/(2 pi^2)) [ I_pp/2 + I_p'p'/2 - I_pp' ]."""
    lam = mp.mpf(M) if lam is None else mp.mpf(lam)
    dE = mp.mpf(dE)
    k = kin(omega, theta_deg, M)
    M = k['M']
    if abs(mp.sin(k['theta'])) < mp.mpf('1e-12') and k['cos_lab'] > 0:
        return 0.0                                      # no recoil, no radiation
    Ep, pp = k['Ep_lab'], k['p_lab']
    if method == 'analytic':
        p = [M, 0, 0, 0]
        # recoil direction: choose z along p'
        q = [Ep, 0, 0, pp]
        Ipp = _I_ii(M, mp.mpf(0), dE, lam)
        Iqq = _I_ii(Ep, pp, dE, lam)
        Ipq = _I_ij(p, q, dE, lam)
        return float(-ALPHA / (2 * PI ** 2) * (Ipp / 2 + Iqq / 2 - Ipq))
    elif method == 'numeric':
        # brute force: p at rest -> integrand depends on |k| and cos(angle(k, p')) only
        def integrand(kk, c):
            wk = mp.sqrt(kk ** 2 + lam ** 2)
            pk = M * wk
            qk = Ep * wk - pp * kk * c
            F = 2 * (M * Ep) / (pk * qk) - M ** 2 / pk ** 2 - M ** 2 / qk ** 2    # -J.J (p.p' = M Ep)
            return kk ** 2 / (2 * wk) * F
        # radial integral split: [0, lam*50] and [lam*50, dE] (log scale) ; angular Gauss-Legendre
        def ang(kk):
            return mp.quad(lambda c: integrand(kk, c), [-1, 0, 1])
        val = mp.quad(ang, [0, lam, 10 * lam, 100 * lam, dE])
        val *= 2 * PI                                   # azimuth
        return float(ALPHA / (2 * PI ** 2) * val)       # (dsigma)_s = dsigma_0 * e^2/(2pi)^3 * int(-J.J)
    else:
        raise ValueError(method)


def soft_cm_DD(omega, theta_deg, dE_cm, M=M_P, lam=None):
    """Denner-Dittmaier Eq. (4.1): soft factor for a cut dE in the CM frame (check)."""
    lam = mp.mpf(M) if lam is None else mp.mpf(lam)
    k = kin(omega, theta_deg, M)
    s, t, M = k['s'], k['t'], k['M']
    m2 = M ** 2
    bt = mp.sqrt(1 - 4 * m2 / t)
    xt = (bt - 1) / (bt + 1)
    Ll = mp.log(lam ** 2 / (4 * mp.mpf(dE_cm) ** 2))
    val = (Ll + (s + m2) / (s - m2) * mp.log(s / m2)
           + (t - 2 * m2) / (bt * t) * (mp.log(xt) * Ll - mp.log(s / m2) ** 2 / 2
                                        - _re_li2(1 - s / m2 * xt) - _re_li2(1 - m2 / s * xt)
                                        + 2 * _re_li2(1 + xt) - PI ** 2 / 2))
    return float(ALPHA / PI * val)


def soft_frame(omega, theta_deg, dE, M=M_P, lam=None, frame='lab'):
    """general-frame soft factor via the 't Hooft-Veltman formula, frame = 'lab' or 'cm' (check vs DD)."""
    lam = mp.mpf(M) if lam is None else mp.mpf(lam)
    dE = mp.mpf(dE)
    k = kin(omega, theta_deg, M)
    M = k['M']
    if frame == 'lab':
        return soft_lab(omega, theta_deg, dE, M, lam)
    Eg, Ep, cth = k['Eg'], k['Ep'], k['cos_cm']
    sth = mp.sqrt(1 - cth ** 2)
    p = [Ep, 0, 0, Eg]                                  # DD: p = E_e(1,0,0,beta), beta = Eg/Ep -> |p| = Eg
    q = [Ep, Eg * sth, 0, Eg * cth]
    Ipp = _I_ii(Ep, Eg, dE, lam)
    Iqq = Ipp
    Ipq = _I_ij(p, q, dE, lam)
    return float(-ALPHA / (2 * PI ** 2) * (Ipp / 2 + Iqq / 2 - Ipq))


# ----------------------------------------------------------------------------------------------------------
# full correction factors
# ----------------------------------------------------------------------------------------------------------
def delta(omega, theta_deg, dE, M=M_P, lam=None, pol=None):
    """O(alpha) correction  dsigma/dsigma_Born - 1  for a soft-photon cut dE (lab).
    pol = None (unpolarised), 'par' (E in the scattering plane), 'perp'."""
    th = float(theta_deg)
    # exactly forward/backward: evaluate the smooth virtual part at a tiny offset (0/0 in the DD formulae)
    th_v = min(max(th, 1e-4), 180 - 1e-4)
    D0, N0, D1, N1 = virtual_parts(omega, th_v, M, lam)
    if pol is None:
        dv = D1 / D0
    elif pol == 'par':
        dv = (D1 + N1) / (D0 + N0)
    elif pol == 'perp':
        dv = (D1 - N1) / (D0 - N0)
    else:
        raise ValueError(pol)
    ds = 0.0 if th < 1e-9 else soft_lab(omega, th_v, dE, M, lam)
    return float(dv) + ds


def delta_virtual_only(omega, theta_deg, M=M_P, lam=None):
    th_v = min(max(float(theta_deg), 1e-4), 180 - 1e-4)
    return delta_virtual(omega, th_v, M, lam)


def sigma3(omega, theta_deg, M=M_P, lam=None, order=0):
    """beam asymmetry Sigma_3 = (dsigma_par - dsigma_perp)/(dsigma_par + dsigma_perp), point Dirac particle:
    order 0 (Born) or 1 (Born + virtual O(alpha); the soft factor cancels)."""
    th_v = min(max(float(theta_deg), 1e-4), 180 - 1e-4)
    D0, N0, D1, N1 = virtual_parts(omega, th_v, M, lam)
    if order == 0:
        return float(N0 / D0)
    return float((N0 + N1) / (D0 + D1))


def delta_sigma3(omega, theta_deg, M=M_P):
    """absolute O(alpha) shift of Sigma_3 (Born + one loop minus Born), lambda-independent."""
    return sigma3(omega, theta_deg, M, order=1) - sigma3(omega, theta_deg, M, order=0)


def delta_low_DD(theta_cm_deg, beta):
    """DD Eq. (3.6): low-energy (beta -> 0) virtual corrections M1/M0 per helicity amplitude (A1..A6 order)."""
    th = np.radians(theta_cm_deg)
    s2, c2 = np.sin(th / 2) ** 2, np.cos(th / 2) ** 2
    ab = ALPHA * beta / np.pi
    return np.array([ab * s2 / (c2 + beta * (2 - c2)), 0.0, ab, -ab / 2, ab / 2, -ab / s2])


# ----------------------------------------------------------------------------------------------------------
# helpers for tables
# ----------------------------------------------------------------------------------------------------------
def table(omegas=(60, 80, 100, 120, 140, 150), thetas=(0, 45, 90, 135, 180), dEs=(1.0, 3.0, 10.0, '3%'), M=M_P):
    """returns dict dE -> array[len(omegas), len(thetas)] of delta in percent."""
    out = {}
    for dE in dEs:
        arr = np.zeros((len(omegas), len(thetas)))
        for i, w in enumerate(omegas):
            de = 0.03 * w if dE == '3%' else dE
            for j, th in enumerate(thetas):
                arr[i, j] = 100 * delta(w, th, de, M)
        out[dE] = arr
    return out


if __name__ == '__main__':
    import sys
    w, th, dE = 100.0, 90.0, 3.0
    print('delta(100 MeV, 90 deg, dE=3 MeV) = %.5f %%' % (100 * delta(w, th, dE)))
