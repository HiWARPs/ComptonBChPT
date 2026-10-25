"""photoprod.py -- tree-level gamma p -> pi N (Born graphs of L^(1) + magnetic couplings), explicit Dirac matrices.
Independent of FORM and of the loop code; used for the forward sum-rule checks of the loop amplitudes:
  alpha+beta = 1/(2pi^2) int sigma_tot / nu^2,   gamma_0 = -1/(4 pi^2) int (sigma_3/2 - sigma_1/2) / nu^3.
Same Feynman rules as graphs.py (e included here):
  gamma N N: i e (Q eps/ - kappa/(2M) ag(eps,q));  pi N N: -(gA/2F) k/ g5 tau^a (k incoming);
  KR: -i e gA/(2F) eps^{3ab} eps/ g5 tau^b;  gamma pi pi: -e eps^{3cd}(k_in+k_out).eps.
Units M = 1.
"""
import numpy as np
I2 = np.eye(2); Z2 = np.zeros((2, 2))
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
g0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
gs = [np.block([[Z2, s], [-s, Z2]]) for s in (sx, sy, sz)]
G = [g0] + gs
g5 = np.block([[Z2, I2], [I2, Z2]]).astype(complex)
metric = np.diag([1., -1, -1, -1])
def sl(v): return sum(metric[m, m] * v[m] * G[m] for m in range(4))
def dot(a, b): return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
def ag(a, b): return sl(a) @ sl(b) - dot(a, b) * np.eye(4)
tau = [sx, sy, sz]
def eps3(a, b): return {(0, 1): 1, (1, 0): -1}.get((a, b), 0)

def spinor(p3, M, hel_axis_spin):
    """Dirac spinor, normalisation ubar u = 2M; chi = 2-spinor"""
    E = np.sqrt(M ** 2 + p3 @ p3)
    chi = hel_axis_spin
    return np.concatenate([np.sqrt(E + M) * chi, (sx * p3[0] + sy * p3[1] + sz * p3[2]) @ chi / np.sqrt(E + M)])
def bar(u): return u.conj() @ g0

def amplitudes(W, th, ph, c, kap_on=True):
    """M[a][(lam, hin, hout)] for pion isospin a (Cartesian), photon helicity lam=+-1 along z,
    nucleon spin projections hin (z), hout (z) of final nucleon; returns dict for 'lo' and 'kap' parts"""
    s = W ** 2; M = 1.0; m = np.sqrt(c['mpi2'])
    qv = (s - M ** 2) / (2 * W)
    kv = np.sqrt((s - (M + m) ** 2) * (s - (M - m) ** 2)) / (2 * W)
    q = np.array([qv, 0, 0, qv]); p = np.array([np.sqrt(M ** 2 + qv ** 2), 0, 0, -qv])
    n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    k = np.concatenate([[np.sqrt(m ** 2 + kv ** 2)], kv * n]); pp = p + q - k
    gA, F, e = c['gA'], c['F'], np.sqrt(4 * np.pi / 137.036)
    Q = np.diag([1., 0.]); K = np.diag([c['kp'], c['kn']])
    out = {}
    for lam in (1, -1):
        eps = np.array([0, -lam / np.sqrt(2), -1j / np.sqrt(2), 0])
        for part in ('lo', 'kap'):
            res = []
            for a in range(3):
                # isospin-Dirac amplitude as 2x2 (isospin) list of 4x4 (Dirac)
                def piNN_out(kk):      # outgoing pion momentum kk, isospin a
                    return [(gA / (2 * F)) * sl(kk) @ g5, tau[a]]
                terms = []
                Sps = sl(p + q) + M * np.eye(4); Spu = sl(p - k) + M * np.eye(4)
                if part == 'lo':
                    Vg = 1j * e * sl(eps); Viso = Q
                else:
                    Vg = 1j * e * (-1 / (2 * M)) * ag(eps, q); Viso = K
                v, ti = piNN_out(k)
                # s channel: piNN(out) S(p+q) gamma
                terms.append((v @ (1j * Sps / (dot(p + q, p + q) - M ** 2)) @ Vg, ti @ Viso))
                # u channel: gamma S(p-k) piNN(out)
                terms.append((Vg @ (1j * Spu / (dot(p - k, p - k) - M ** 2)) @ v, Viso @ ti))
                if part == 'lo':
                    # Kroll-Ruderman: -i e gA/(2F) eps^{3 a b} eps/ g5 tau^b, outgoing pion a
                    terms.append((-1j * e * gA / (2 * F) * sl(eps) @ g5, sum(eps3(a, b) * tau[b] for b in range(3))))
                    # pion pole: photon on pion: incoming pion c (momentum k-q) -> outgoing a (momentum k)
                    kt = k - q
                    for cc in range(3):
                        if eps3(cc, a) == 0: continue
                        vpp = -e * eps3(cc, a) * dot(kt + k, eps)
                        prop = 1j / (dot(kt, kt) - m ** 2)
                        # pion c absorbed... created at the nucleon vertex: outgoing momentum kt from the nucleon line
                        vN = (gA / (2 * F)) * sl(kt) @ g5
                        terms.append((vpp * prop * vN, tau[cc]))
                res.append(terms)
            out[(lam, part)] = (res, p, pp, qv, kv)
    return out

def helicity_xs(W, c, nth=48, nph=16):
    """sigma_1/2, sigma_3/2 (LO^2 and 2 Re LO* kap) summed over pi0 p and pi+ n, in M^-2 units"""
    s = W ** 2; M = 1.0
    xg, wg = np.polynomial.legendre.leggauss(nth)
    up = np.array([1, 0], complex); dn = np.array([0, 1], complex)
    res = {}
    for lam in (1,):
        sig = {('lo', 'h12'): 0., ('lo', 'h32'): 0., ('int', 'h12'): 0., ('int', 'h32'): 0.}
        for xc, wc in zip(xg, wg):
            th = np.arccos(xc)
            for ph in np.linspace(0, 2 * np.pi, nph, endpoint=False):
                amp = amplitudes(W, th, ph, c)
                (rl, p, pp, qv, kv) = amp[(lam, 'lo')]; (rk, *_ ) = amp[(lam, 'kap')]
                for hin, tag in ((dn, 'h12'), (up, 'h32')):
                    # target moving along -z: spin projection on z; photon helicity +1 along +z.
                    # total J_z = +1 + s_z: s_z=+1/2 -> 3/2, s_z=-1/2 -> 1/2
                    u = spinor(p[1:], M, hin)
                    for hout in (up, dn):
                        ub = bar(spinor(pp[1:], M, hout))
                        for a in range(3):
                            for Nf in range(2):          # final nucleon p (0) or n (1)
                                Ml = sum((ub @ D @ u) * T[Nf, 0] for D, T in rl[a])
                                Mk = sum((ub @ D @ u) * T[Nf, 0] for D, T in rk[a])
                                w = wc * (2 * np.pi / nph)
                                sig[('lo', tag)] += w * abs(Ml) ** 2
                                sig[('int', tag)] += w * 2 * (np.conj(Ml) * Mk).real
        flux = kv / (64 * np.pi ** 2 * s * qv)
        for k_ in sig: sig[k_] *= flux
        res = sig
    return res
