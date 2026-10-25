"""tchan.py -- O(p^4) graphs with a c_i pi pi N N vertex X and a closed pion loop (t-channel graphs).

L^(2) (Fettes et al. 2000):  c1 <chi_+> - c2/(8M^2) (<u.u> D D + h.c.) + c3/2 <u.u> + i c4/4 sigma [u,u] ...
Vertex rules (pions a: k1, b: k2 incoming; e removed; photon couplings from D_mu pi only):
  c1  pi pi N N        -4i c1 mpi^2/F^2 delta^{ab}
  c2  pi pi N N        -i c2/(M^2 F^2) delta^{ab} [(k1.p)(k2.p) + (k1.pp)(k2.pp)]
  c3  pi pi N N        -2i c3/F^2 delta^{ab} k1.k2
  c2  gamma pi pi N N  c2/(M^2 F^2) eps^{3ab} [((k2-k1).p)(eps.p) + ((k2-k1).pp)(eps.pp)]
  c3  gamma pi pi N N  2 c3/F^2 eps^{3ab} (k2-k1).eps
  c2  gg pi pi N N     2i c2/(M^2 F^2) P^{ab} [(ee.p)(ep.p) + (ee.pp)(ep.pp)],   P^{ab} = delta^{ab} - delta^{a3}delta^{b3}
  c3  gg pi pi N N     4i c3/F^2 P^{ab} ee.ep
c4 and the Weinberg-Tomozawa vertex give zero for the proton by isospin (eps^{abc} against a symmetric loop tensor).
Symmetry factors: triangle 1 (one orientation), two-vertex bubbles 1/2, tadpole 1/2.
"""
import numpy as np
from fractions import Fraction as Fr
from graphs import Mom, Tok, Graph, P_, Q_, QP_, L_, PP_, eps3

PHT = {'A': (Q_, 'ee'), 'B': (-QP_, 'ep')}
E3 = np.array([[eps3(i, j) for j in range(3)] for i in range(3)], float)
Pm = np.eye(3) - np.diag([0, 0, 1.])

def dotform(m, v):
    return m.form(lambda w: f'{w}.{v}')

def vertexX(ci, nph, k1, k2, phs):
    """(pref, FORM scalar text, isospin matrix V^{ab}) of the c_i vertex"""
    if nph == 0:
        if ci == 'c1': return -4j, 'mpi2', np.eye(3)
        if ci == 'c2':
            return -1j, f'({dotform(k1,"p")}*{dotform(k2,"p")}+{dotform(k1,"pp")}*{dotform(k2,"pp")})', np.eye(3)
        if ci == 'c3': return -2j, k1k2(k1, k2), np.eye(3)
    if nph == 1:
        (mom, pol), = [PHT[x] for x in phs]
        d = k2 - k1
        if ci == 'c2':
            return 1.0, f'({dotform(d,"p")}*p.{pol}+{dotform(d,"pp")}*pp.{pol})', E3
        if ci == 'c3':
            return 2.0, dotform(d, pol), E3
    if nph == 2:
        if ci == 'c2': return 2j, '(ee.p*ep.p+ee.pp*ep.pp)', Pm
        if ci == 'c3': return 4j, 'ee.ep', Pm
    return None

def k1k2(k1, k2):
    """FORM text of k1.k2 for two Mom's"""
    terms = []
    for a, ca in k1.items():
        for b, cb in k2.items():
            terms.append(f'{ca*cb}*{a}.{b}')
    return '(' + '+'.join(terms).replace('+-', '-') + ')'

def graphs():
    gs = []
    for ci in ('c1', 'c2', 'c3'):
        # (photons at X, photons on the loop in order)
        configs = [((), ('A', 'B')), ((), ('sg',))]
        if ci != 'c1':
            configs += [(('A',), ('B',)), (('B',), ('A',)), (('A', 'B'), ())]
        for atX, onloop in configs:
            g = Graph(f't_{ci}_X{"".join(atX)}_L{"".join(onloop)}')
            pref = 1 + 0j; scal = []; T = np.eye(3)
            k = L_; props = [k]; pref *= 1j
            for ph in onloop:
                if ph == 'sg':
                    scal.append('ee.ep'); pref *= 2j; T = T @ Pm; k = k + Q_ - QP_
                else:
                    mom, pol = PHT[ph]; kout = k + mom
                    scal.append((k + kout, pol)); pref *= -1; T = T @ E3; k = kout
                props.append(k); pref *= 1j
            k1, k2 = -L_, k                 # incoming momenta at X
            pv, sx, V = vertexX(ci, len(atX), k1, k2, atX)
            pref *= pv
            scal.append(sx)
            nlines = len(props)
            sym = 0.5 if nlines == 2 or nlines == 1 else 1.0
            g.toks = [Tok('I')]
            g.scal = scal
            g.loop = [(m, 'mpi2') for m in props]
            g.pref = pref * sym * (-1j) * 1j / (16 * np.pi ** 2)
            g.coup = {ci: 1, 'F': -2}
            iso = float(np.sum(V * T))
            g.iso = {'1': iso}
            g.kind = 'tch'
            if abs(iso) > 1e-14: gs.append(g)
    return gs

if __name__ == '__main__':
    for g in graphs():
        print(g.name, 'props', len(g.loop), 'iso', g.iso, 'pref*16pi^2', np.round(g.pref * 16 * np.pi ** 2, 6), g.scal)
