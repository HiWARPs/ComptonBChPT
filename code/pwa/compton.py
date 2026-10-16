"""
compton.py -- real Compton scattering off the proton: kinematics, exact Born(+pi0-pole)
helicity amplitudes, partial-wave / multipole machinery and observables.

Conventions (Hagelstein-Miskimen-Pascalutsa, PPNP 88 (2016) 29, Sec. 3 = "HMP16"):
  * c.m. helicity amplitudes T_{H'H}(omega,theta), H = lambda_gamma - lambda_N, normalised as
    Phi = T/(8 pi sqrt(s)); dsigma/dOmega_cm = 1/2|Phi1|^2+1/2|Phi2|^2+|Phi3|^2+|Phi4|^2+1/2|Phi5|^2+1/2|Phi6|^2
  * partial waves  T_{H'H} = sum_J (2J+1) T^J_{H'H} d^J_{H H'}(theta)    (Edmonds d-functions)
  * multipoles f^{l+-}_{rho' rho} from A^J_p, B^J_p, C^J_p   (HMP16 Eqs. 3.19, 3.21, 3.22)
The helicity-amplitude array used internally is  M[i,j,hp,h]  with exactly the index order and
phase conventions of the Fortran code ampl.f90 (i/j = final/initial photon helicity index,
0 -> +1, 1 -> -1; hp/h = final/initial nucleon helicity index, 0 -> +1/2, 1 -> -1/2), but
normalised to Phi (i.e. M = T/(8 pi sqrt s)).  All amplitudes are returned in fm.

Units: energies in MeV on input; internally everything in units of the nucleon mass M (as in the
Fortran) and converted to fm with hbarc/M.
"""
import numpy as np
from math import lgamma, exp

def factorial(n):
    return exp(lgamma(n + 1.0))

# ---------------- constants (identical to the Fortran code ampl.f90) ----------------
MN_GEV = 0.93827          # proton mass used by the Fortran (smm)
HBARC = 0.19733           # GeV fm (Fortran value)
ALPHA = 1.0 / 137.036
KAPPA_P = 1.793
GA = 1.27
FPI = 0.09221             # GeV
MPI0 = 0.13498            # GeV
FM_PER_M = HBARC / MN_GEV  # 1/M in fm
NB_PER_FM2 = 1.0e7

# polarisability conversions: 1e-4 fm^3 and 1e-4 fm^4 into M-units (alpha_em NOT included)
P3 = 1.0e-4 / FM_PER_M**3
P4 = 1.0e-4 / FM_PER_M**4


# ---------------- kinematics ----------------
def kin(Elab_MeV, theta_deg, frame='lab'):
    """Return dict of kinematic quantities (M units) for lab photon energy and angle."""
    E = np.asarray(Elab_MeV, float) * 1e-3 / MN_GEV
    th = np.radians(np.asarray(theta_deg, float))
    s = 1.0 + 2.0 * E
    W = np.sqrt(s)
    w = (s - 1.0) / (2 * W)          # c.m. photon energy
    EN = (s + 1.0) / (2 * W)         # c.m. nucleon energy
    if frame == 'lab':
        ct_lab = np.cos(th)
        ct = 1.0 - (1.0 + 2.0 * E) * (1.0 - ct_lab) / (1.0 + E * (1.0 - ct_lab))
        ct = np.clip(ct, -1, 1)
        thcm = np.arccos(ct)
        thlab = th
    else:
        thcm = th
        ct = np.cos(th)
        thlab = np.arccos(np.clip(1.0 - (1.0 - ct) / (1.0 + E * (1.0 + ct)), -1, 1))
    t = -2.0 * w**2 * (1.0 - np.cos(thcm))
    u = 2.0 - s - t
    nup = (1.0 - u) / 2.0            # lab energy of outgoing photon (M units)
    return dict(E=E, s=s, W=W, w=w, EN=EN, thcm=thcm, thlab=thlab, t=t, u=u, nu=E, nup=nup)


def theta_cm_from_lab(Elab_MeV, thlab_deg):
    return np.degrees(kin(Elab_MeV, thlab_deg, 'lab')['thcm'])


def lab_factor(k):
    """dsigma/dOmega_lab = lab_factor * dsigma/dOmega_cm  (HMP16 Eq. 3.11)."""
    return (k['nup'] / k['nu'])**2 * k['s']


# ---------------- Dirac algebra (Dirac representation, as in dirac.f90) ----------------
I2 = np.eye(2)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
Z2 = np.zeros((2, 2))
G0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
G1 = np.block([[Z2, SX], [-SX, Z2]])
G2 = np.block([[Z2, SY], [-SY, Z2]])
G3 = np.block([[Z2, SZ], [-SZ, Z2]])
G5 = np.block([[Z2, I2], [I2, Z2]]).astype(complex)
GAM = np.array([G0, G1, G2, G3])
METRIC = np.array([1.0, -1.0, -1.0, -1.0])
ONE4 = np.eye(4, dtype=complex)


def slash(a):
    """a: (...,4) contravariant components -> (...,4,4)."""
    a = np.asarray(a)
    return np.einsum('...m,mij->...ij', a * METRIC, GAM)


def mdot(a, b):
    return np.sum(np.asarray(a) * np.asarray(b) * METRIC, axis=-1)


def levi4():
    eps = np.zeros((4, 4, 4, 4))
    import itertools
    for p in itertools.permutations(range(4)):
        sgn = 1
        pl = list(p)
        for i in range(4):
            for j in range(i + 1, 4):
                if pl[i] > pl[j]:
                    sgn = -sgn
        eps[p] = sgn
    return eps          # eps^{0123}=+1 (upper indices), as lct4 in dirac.f90


EPS4 = levi4()


def eps_contract(a, b, c, d):
    """eps_{mu nu al be} a^mu b^nu c^al d^be with eps^{0123}=+1 (=> eps_{0123}=-1)."""
    am, bm, cm, dm = (np.asarray(x) * METRIC for x in (a, b, c, d))   # lower indices
    return np.einsum('...m,...n,...a,...b,mnab->...', am, bm, cm, dm, EPS4)


def spinor_u(p, theta, hel):
    """Helicity spinor of dirac.f90 u_(p,theta,0,hel,1) (M units), normalised to ubar u = 2M (=2)."""
    p = np.asarray(p, float)
    theta = np.asarray(theta, float)
    ep = np.sqrt(p * p + 1.0)
    comp = np.sqrt((ep + 1.0) / (2.0 * ep))
    if hel > 0:
        ih = 1
        c1, c2 = np.cos(theta / 2), np.sin(theta / 2)
    else:
        ih = -1
        c1, c2 = -np.sin(theta / 2), np.cos(theta / 2)
    r1 = comp
    r2 = ih * p * comp / (ep + 1.0)
    u = np.stack([r1 * c1, r1 * c2, r2 * c1, r2 * c2], axis=-1).astype(complex)
    return u * np.sqrt(2 * ep)[..., None]


def ubar(u):
    return np.conj(u) @ G0


# photon polarisation vectors (dirac.f90, ini_polar_real): eps_{+1}, eps_{-1} of the incoming photon
def eps_in():
    sq = 1 / np.sqrt(2)
    return np.array([[0, -sq, -1j * sq, 0], [0, sq, -1j * sq, 0]], complex)


def eps_out_conj(theta):
    """complex-conjugated polarisation vectors of the outgoing photon (lambda'=+1,-1)."""
    sq = 1 / np.sqrt(2)
    X, Z = np.sin(theta), np.cos(theta)
    z0 = np.zeros_like(X)
    e1 = np.stack([z0, -sq * Z, 1j * sq + z0, sq * X], axis=-1)
    e2 = np.stack([z0, sq * Z, 1j * sq + z0, -sq * X], axis=-1)
    return np.stack([e1, e2], axis=-2)        # (...,2,4)


# ---------------- exact Born + pi0-pole helicity amplitudes ----------------
def born_amplitudes(Elab_MeV, theta_deg, frame='lab', kappa=KAPPA_P, pi0=True, born=True,
                    gA=GA, fpi=FPI, mpi0=MPI0):
    """Return M[...,i,j,hp,h] = T/(8 pi sqrt s) in fm for the Born (Dirac+Pauli, nucleon s+u poles)
    and/or the pi0 t-channel pole.  Computed directly from Feynman rules with numerical Dirac algebra."""
    k = kin(Elab_MeV, theta_deg, frame)
    s, W, w, EN, th, t, u_ = k['s'], k['W'], k['w'], k['EN'], k['thcm'], k['t'], k['u']
    shape = np.broadcast(s, th).shape
    s = np.broadcast_to(s, shape); W = np.broadcast_to(W, shape); w = np.broadcast_to(w, shape)
    EN = np.broadcast_to(EN, shape); th = np.broadcast_to(th, shape)
    t = np.broadcast_to(t, shape); u_ = np.broadcast_to(u_, shape)
    z = np.zeros(shape)
    q = np.stack([w, z, z, w], -1)
    p = np.stack([EN, z, z, -w], -1)
    qp = np.stack([w, w * np.sin(th), z, w * np.cos(th)], -1)
    pp = np.stack([EN, -w * np.sin(th), z, -w * np.cos(th)], -1)
    e2 = 4 * np.pi * ALPHA
    epsin = eps_in()                      # (2,4)
    epsout = eps_out_conj(th)             # (...,2,4)
    Mout = np.zeros(shape + (2, 2, 2, 2), complex)
    # spinors
    uin = {hi: spinor_u(w, np.pi + 0 * th, hel) for hi, hel in ((0, 0.5), (1, -0.5))}
    uout = {hi: ubar(spinor_u(w, th + np.pi, hel)) for hi, hel in ((0, 0.5), (1, -0.5))}
    sl_q, sl_qp = slash(q), slash(qp)
    Sprop = (slash(p + q) + ONE4) / (s - 1.0)[..., None, None]
    Uprop = (slash(p - qp) + ONE4) / (u_ - 1.0)[..., None, None]

    def vertex(epsv, sl_k):
        se = slash(epsv)
        return se - (kappa / 4.0) * (se @ sl_k - sl_k @ se)

    for i in range(2):
        epo = epsout[..., i, :]
        for j in range(2):
            epi = np.broadcast_to(epsin[j], shape + (4,))
            amp = np.zeros(shape + (4, 4), complex)
            if born:
                Vout = vertex(epo, -sl_qp)
                Vin = vertex(epi, sl_q)
                amp = amp - e2 * (Vout @ Sprop @ Vin + Vin @ Uprop @ Vout)
            if pi0:
                gg = gA / (2 * np.pi * fpi / MN_GEV)**2          # M units
                coef = 4 * np.pi * ALPHA * gg / ((mpi0 / MN_GEV)**2 - t)
                epsc = eps_contract(epo, epi, q, qp)
                amp = amp + (coef * 1j * epsc)[..., None, None] * G5
            for hp in range(2):
                for h in range(2):
                    Mout[..., i, j, hp, h] = np.einsum('...a,...ab,...b->...', uout[hp], amp, uin[h])
    Mout = Mout / (8 * np.pi * W)[..., None, None, None, None]
    return Mout * FM_PER_M


# ---------------- helicity bookkeeping ----------------
# H -> (photon index, nucleon index)
HIDX = {0.5: (0, 0), -0.5: (1, 1), 1.5: (0, 1), -1.5: (1, 0)}


def T_HH(M, Hp, H):
    ig_p, in_p = HIDX[Hp]
    ig, iN = HIDX[H]
    return M[..., ig_p, ig, in_p, iN]


def phis_from_M(M):
    """Phi_1..6 (HMP16 Eq. 3.8) from the helicity array M (already normalised to T/(8 pi sqrt s))."""
    return np.stack([T_HH(M, -0.5, -0.5), T_HH(M, -0.5, 0.5), T_HH(M, -0.5, 1.5),
                     T_HH(M, -0.5, -1.5), T_HH(M, 1.5, 1.5), T_HH(M, -1.5, 1.5)], axis=-1)


def M_from_phis(P):
    """inverse of phis_from_M using the P,T relations of HMP16 Eq. (3.8)."""
    P = np.asarray(P)
    M = np.zeros(P.shape[:-1] + (2, 2, 2, 2), complex)
    rel = {(-0.5, -0.5): (0, 1), (0.5, 0.5): (0, 1),
           (-0.5, 0.5): (1, 1), (0.5, -0.5): (1, -1),
           (-0.5, 1.5): (2, 1), (0.5, -1.5): (2, 1), (1.5, -0.5): (2, 1), (-1.5, 0.5): (2, 1),
           (-0.5, -1.5): (3, 1), (0.5, 1.5): (3, -1), (1.5, 0.5): (3, 1), (-1.5, -0.5): (3, -1),
           (1.5, 1.5): (4, 1), (-1.5, -1.5): (4, 1),
           (-1.5, 1.5): (5, 1), (1.5, -1.5): (5, -1)}
    for (Hp, H), (k, sg) in rel.items():
        ig_p, in_p = HIDX[Hp]
        ig, iN = HIDX[H]
        M[..., ig_p, ig, in_p, iN] = sg * P[..., k]
    return M


# ---------------- Wigner d-functions (dfun of ampl.f90, Edmonds convention) ----------------
def dfun(j, mp, m, beta):
    beta = np.asarray(beta, float)
    c2, s2 = np.cos(beta / 2), np.sin(beta / 2)
    jmp, jm, jmp1, jm1, mpm = int(round(j + mp)), int(round(j - m)), int(round(j - mp)), int(round(j + m)), int(round(mp - m))
    lmin = max(0, mpm)
    lmax = min(jmp, jm)
    res = np.zeros_like(beta)
    for l in range(lmin, lmax + 1):
        res = res + (-1)**l * c2**(jm + jmp - 2 * l) * s2**(2 * l - mpm) / (
            factorial(jmp - l) * factorial(jm - l) * factorial(l) * factorial(l - mpm))
    return res * np.sqrt(factorial(jmp) * factorial(jmp1) * factorial(jm1) * factorial(jm))


# (H, H') of the d-function d^J_{H H'} for each Phi_k = T_{H'H}
PHI_HH = [(-0.5, -0.5), (0.5, -0.5), (1.5, -0.5), (-1.5, -0.5), (1.5, 1.5), (1.5, -1.5)]


def dmatrix(J, theta):
    """returns (...,6) d^J_{H H'}(theta) for the six Phi's (zero where |H|>J or |H'|>J)."""
    out = []
    for (H, Hp) in PHI_HH:
        if abs(H) > J or abs(Hp) > J:
            out.append(np.zeros_like(np.asarray(theta, float)))
        else:
            out.append(dfun(J, H, Hp, theta))
    return np.stack(out, -1)


def pw_project(phis_on_grid, xg, wg, Jmax2=21):
    """Phi^J_k = 1/2 int dcos Phi_k(theta) d^J_{H H'}(theta); grid = Gauss-Legendre in cos theta.
    phis_on_grid: (Nx,6) [or (...,Nx,6)].  Returns dict J -> (...,6)."""
    th = np.arccos(xg)
    res = {}
    for j2 in range(1, Jmax2 + 1, 2):
        J = j2 / 2
        d = dmatrix(J, th)                          # (Nx,6)
        res[J] = 0.5 * np.einsum('...xk,xk,x->...k', phis_on_grid, d, wg)
    return res


def phis_from_pw(pw, theta):
    out = 0
    for J, PJ in pw.items():
        out = out + (2 * J + 1) * PJ[..., None, :] * dmatrix(J, theta)
    return out


# ---------------- partial waves <-> multipoles (HMP16 Eqs. 3.19, 3.21, 3.22) ----------------
def ABC_from_pw(PJ):
    """PJ (...,6) = Phi^J_1..6 -> dict p -> (A,B,C) with p=+1,-1."""
    P1, P2, P3_, P4_, P5, P6 = [PJ[..., k] for k in range(6)]
    out = {}
    for p in (+1, -1):
        out[p] = (P1 + p * P2, -p * P3_ - P4_, P5 + p * P6)
    return out


def multipoles_from_pw(pw, lmax=None):
    """Return dict with keys ('EE',l,'-'),('MM',l,'-'),('EE',l,'+'),('MM',l,'+'),('EM',l,'+'),('ME',l,'+')."""
    f = {}
    for J, PJ in pw.items():
        abc = ABC_from_pw(PJ)
        if J == 0.5:
            f[('EE', 1, '-')] = 0.5 * abc[+1][0]
            f[('MM', 1, '-')] = 0.5 * abc[-1][0]
            continue
        c = 2.0 / (2 * J + 1)**2
        lm = int(round(J + 0.5))      # (J+1/2)-
        lp = int(round(J - 0.5))      # (J-1/2)+
        for rho, p in (('EE', +1), ('MM', -1)):
            A, B, C = abc[p]
            f[(rho, lm, '-')] = c * (A - 2 * np.sqrt((J - 0.5) / (J + 1.5)) * B + (J - 0.5) / (J + 1.5) * C)
            A, B, C = abc[-p]
            f[(rho, lp, '+')] = c * (A + 2 * np.sqrt((J + 1.5) / (J - 0.5)) * B + (J + 1.5) / (J - 0.5) * C)
        for rho, p in (('EM', +1), ('ME', -1)):
            A, B, C = abc[-p]
            f[(rho, lp, '+')] = c * (-A - 2.0 / np.sqrt((J + 1.5) * (J - 0.5)) * B + C)
    return f


def _mult_keys(J):
    lm = int(round(J + 0.5)); lp = int(round(J - 0.5))
    if J == 0.5:
        return [('EE', 1, '-'), ('MM', 1, '-')]
    return [('EE', lm, '-'), ('MM', lm, '-'), ('EE', lp, '+'), ('MM', lp, '+'), ('EM', lp, '+'), ('ME', lp, '+')]


_RINV = {}


def _pw_to_mult_matrix(J):
    """6x6 (or 2x2 for J=1/2) real matrix R with f_keys = R @ Phi^J (from multipoles_from_pw)."""
    keys = _mult_keys(J)
    n = 2 if J == 0.5 else 6
    R = np.zeros((len(keys), n))
    for k in range(n):
        e = np.zeros(6); e[k] = 1.0
        f = multipoles_from_pw({J: e})
        for a, key in enumerate(keys):
            R[a, k] = f[key]
    return keys, R


def pw_from_multipoles(f, Jmax2):
    """inverse of multipoles_from_pw: multipoles -> Phi^J_k (k=1..6), J=1/2..Jmax2/2 (missing = 0)."""
    pw = {}
    for j2 in range(1, Jmax2 + 1, 2):
        J = j2 / 2
        if J not in _RINV:
            keys, R = _pw_to_mult_matrix(J)
            _RINV[J] = (keys, np.linalg.inv(R))
        keys, Rinv = _RINV[J]
        fv = np.stack(np.broadcast_arrays(*[np.asarray(f.get(key, 0.0)) for key in keys]), -1)
        ph = np.einsum('kf,...f->...k', Rinv, fv)
        if J == 0.5:
            ph = np.concatenate([ph, np.zeros(ph.shape[:-1] + (4,), dtype=ph.dtype)], -1)
        pw[J] = ph
    return pw


# ---------------- observables ----------------
def obs_from_M(M, k):
    """M[...,i,j,hp,h] in fm; returns dict with cs_cm, cs_lab, cs_par_lab, cs_perp_lab (nb/sr), S3."""
    cs = 0.25 * np.sum(np.abs(M)**2, axis=(-4, -3, -2, -1))
    Mx = -(M[..., :, 0, :, :] - M[..., :, 1, :, :]) / np.sqrt(2)
    My = -1j * (M[..., :, 0, :, :] + M[..., :, 1, :, :]) / np.sqrt(2)
    spar = 0.5 * np.sum(np.abs(Mx)**2, axis=(-3, -2, -1))
    sperp = 0.5 * np.sum(np.abs(My)**2, axis=(-3, -2, -1))
    lf = lab_factor(k)
    out = dict(cs_cm=cs * NB_PER_FM2, cs_lab=cs * lf * NB_PER_FM2,
               cs_par_lab=spar * lf * NB_PER_FM2, cs_perp_lab=sperp * lf * NB_PER_FM2,
               S3=(spar - sperp) / (spar + sperp))
    return out


def obs_from_phis(P, k):
    """same, from Phi_1..6 (HMP16 Eqs. 3.10, 3.16)."""
    a = np.abs(P)**2
    cs = 0.5 * a[..., 0] + 0.5 * a[..., 1] + a[..., 2] + a[..., 3] + 0.5 * a[..., 4] + 0.5 * a[..., 5]
    S3num = -np.real((P[..., 0] + P[..., 4]) * np.conj(P[..., 2]) + (P[..., 1] - P[..., 5]) * np.conj(P[..., 3]))
    lf = lab_factor(k)
    return dict(cs_cm=cs * NB_PER_FM2, cs_lab=cs * lf * NB_PER_FM2, S3=S3num / cs,
                cs_par_lab=(cs + S3num) * lf * NB_PER_FM2, cs_perp_lab=(cs - S3num) * lf * NB_PER_FM2)


# ---------------- Fortran interface ----------------
def read_fortran_hel(path):
    """parse output of run_helpoints: returns list of ids and dict arrays.
    comps[c] -> M[N,i,j,hp,h] in Fortran normalisation (u^+u=1, includes alpha), M-units."""
    ids, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            ids.append(parts[0])
            rows.append([float(x) for x in parts[1:]])
    a = np.array(rows)
    head = a[:, :6]
    vals = a[:, 6:].reshape(len(a), 6, 2, 2, 2, 2, 2)   # c, h, hp, j, i, reim
    cplx = vals[..., 0] + 1j * vals[..., 1]
    M = np.transpose(cplx, (0, 1, 5, 4, 3, 2))          # N,c,i,j,hp,h
    return ids, dict(Elab=head[:, 0], thlab=head[:, 1], thcm=head[:, 2], s=head[:, 3], EN=head[:, 4],
                     q=head[:, 5], M=M)


def fortran_to_phi_norm(M, s, EN):
    """Fortran mgg -> T/(8 pi sqrt s) in fm:  Phi = (E_N/sqrt s) * mgg  (see notes)."""
    fac = (EN / np.sqrt(s))
    return M * fac.reshape(fac.shape + (1,) * (M.ndim - fac.ndim)) * FM_PER_M
