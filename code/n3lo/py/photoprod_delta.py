"""photoprod_delta.py -- tree-level gamma p -> pi N with Delta exchange (s and u channel), explicit 4x4 Dirac matrices
and an explicit Levi-Civita tensor (eps^{0123} = -1, i.e. eps_{0123} = +1).  Independent of FORM and of delta.py's
gamma5 representation of eps; same Lagrangian (see delta.py), e included:
  pi N Delta (Delta in p, pion in k):   c gamma^{mu nu lam} p_mu k_lam T^a,      c = hA/(2 F MD)
  pi Delta N (Delta out p', pion in k): -c gamma^{mu nu lam} p'_mu k_lam T^a+
  gamma N Delta (Delta in p, photon in q, pol e):  e [cM eps^{mu nu rho lam} p_mu q_rho e_lam + i cE g5 ((p.q) e^nu - (p.e) q^nu)] T^3
  gamma Delta N (Delta out p'):                    same with p -> p', T^3 -> T^3+
  Delta propagator  i N^{ab}(p)/(p^2 - MD^2),  N = -(p/+MD)[g - g g/3 - 2pp/(3MD^2) + (p g - g p)/(3MD)]   (4 dim)
Used with photoprod.amplitudes (Born, 'lo') for the sum rules of the gamma N Delta loop graphs (interference terms).
"""
import itertools, numpy as np
import photoprod as pp
from photoprod import sl, dot, g5, G, metric
import delta

def gam(mu):                     # gamma^mu (upper index)
    return G[mu]
def low(v): return metric @ v    # lower the index of a 4-vector
_LC = np.zeros((4, 4, 4, 4))
for perm in itertools.permutations(range(4)):
    _LC[perm] = -np.linalg.det(np.eye(4)[list(perm)])          # eps^{0123} = -1 (upper indices)

def g3_upper(a, b):
    """gamma^{mu nu lam} a_mu b_lam for each upper nu: list of 4 matrices (a, b: upper-index 4-vectors)"""
    A, B = sl(a), sl(b)
    return [(A @ gam(nu) @ B - B @ gam(nu) @ A) / 2 for nu in range(4)]

def rs_num(p, MD):
    """N^{al be}(p) as 4x4 array of 4x4 matrices (upper indices)"""
    g = metric
    S = sl(p) + MD * np.eye(4)
    out = np.empty((4, 4), object)
    for a in range(4):
        for b in range(4):
            br = g[a, b] * np.eye(4) - gam(a) @ gam(b) / 3 - 2 * p[a] * p[b] / (3 * MD ** 2) * np.eye(4) \
                 + (p[a] * gam(b) - gam(a) * p[b]) / (3 * MD)
            out[a, b] = -S @ br
    return out

def gnd_vertex(pD, q, e, cM, cE):
    """V^nu (upper nu) of the gamma N Delta vertex for Delta momentum pD, incoming photon q, polarisation e"""
    pl, ql, el = low(pD), low(q), low(e)
    V = []
    for nu in range(4):
        epsv = np.einsum('mrl,m,r,l->', _LC[:, nu, :, :], pl, ql, el)
        V.append(cM * epsv * np.eye(4) + 1j * cE * g5 * (dot(pD, q) * e[nu] - dot(pD, e) * q[nu]))
    return V

def contract(Vl, N, Vr):
    """sum_{a,b} Vl^a g_aa N^{ab} g_bb Vr^b  (Vl, Vr: lists of upper-index matrices)"""
    g = metric
    return sum(g[a, a] * g[b, b] * Vl[a] @ N[a, b] @ Vr[b] for a in range(4) for b in range(4))

def delta_amps(W, th, ph, c, which=('s', 'u'), lam=1):
    """Delta-exchange amplitudes [pion a][(Dirac, isospin 2x2)] for photon helicity lam along +z"""
    s = W ** 2; M = 1.0; m = np.sqrt(c['mpi2']); MD = c['MD']
    qv = (s - M ** 2) / (2 * W)
    kv = np.sqrt((s - (M + m) ** 2) * (s - (M - m) ** 2)) / (2 * W)
    q = np.array([qv, 0, 0, qv]); p = np.array([np.sqrt(M ** 2 + qv ** 2), 0, 0, -qv])
    n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    k = np.concatenate([[np.sqrt(m ** 2 + kv ** 2)], kv * n]); pp_ = p + q - k
    e = np.array([0, -lam / np.sqrt(2), -1j / np.sqrt(2), 0])
    ech = np.sqrt(4 * np.pi / 137.036)
    cpi = c['hA'] / (2 * c['F'] * MD)
    cM = 1.5 * c['gM'] / (1 + MD) * delta.SGND; cE = 1.5 * c['gE'] / (1 + MD) * delta.SGND   # see delta.SGND
    res = []
    for a in range(3):
        terms = []
        if 's' in which:      # gamma N -> Delta(p+q) -> pi(k) N
            P = p + q
            Vg = [ech * v for v in gnd_vertex(P, q, e, cM, cE)]                 # Deltabar <- N  (Delta out P): T3+
            Vp = [cpi * x for x in g3_upper(P, -k)]                            # Nbar <- Delta (Delta in P), pion in -k: T^a
            prop = 1j / (dot(P, P) - MD ** 2)
            terms.append((prop * contract(Vp, rs_num(P, MD), Vg), delta.T[a] @ delta.Td[2]))
        if 'u' in which:      # N -> pi(k) Delta(p-k), Delta + gamma -> N
            P = p - k
            Vp = [-cpi * x for x in g3_upper(P, -k)]                           # Deltabar <- N (Delta out P), pion in -k: T^a+
            Vg = [ech * v for v in gnd_vertex(P, q, e, cM, cE)]                 # Nbar <- Delta (Delta in P): T3
            prop = 1j / (dot(P, P) - MD ** 2)
            terms.append((prop * contract(Vg, rs_num(P, MD), Vp), delta.T[2] @ delta.Td[a]))
        res.append(terms)
    return res, p, pp_, qv, kv

def helicity_xs_int(W, c, which=('u',), nth=24, nph=8):
    """2 Re(Born_LO^* x Delta) in sigma_1/2, sigma_3/2 (pi0 p + pi+ n), M^-2 units"""
    s = W ** 2; M = 1.0
    xg, wg = np.polynomial.legendre.leggauss(nth)
    up = np.array([1, 0], complex); dn = np.array([0, 1], complex)
    sig = {'h12': 0., 'h32': 0.}
    for xc, wc in zip(xg, wg):
        th = np.arccos(xc)
        for ph in np.linspace(0, 2 * np.pi, nph, endpoint=False):
            amp = pp.amplitudes(W, th, ph, c)
            rl, p, pp_, qv, kv = amp[(1, 'lo')]
            rd, *_ = delta_amps(W, th, ph, c, which)
            for hin, tag in ((dn, 'h12'), (up, 'h32')):
                u = pp.spinor(p[1:], M, hin)
                for hout in (up, dn):
                    ub = pp.bar(pp.spinor(pp_[1:], M, hout))
                    for a in range(3):
                        for Nf in range(2):
                            Ml = sum((ub @ D @ u) * T[Nf, 0] for D, T in rl[a])
                            Md = sum((ub @ D @ u) * T[Nf, 0] for D, T in rd[a])
                            sig[tag] += wc * (2 * np.pi / nph) * 2 * (np.conj(Ml) * Md).real
    flux = kv / (64 * np.pi ** 2 * s * qv)
    return {k_: v * flux for k_, v in sig.items()}
