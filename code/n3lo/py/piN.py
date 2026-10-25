"""piN.py -- all one-loop graphs built on the pi N self-energy skeleton (two pi N N vertices V1, V2)
with the two photons attached in every possible way, plus their renormalisation subtraction terms.

Photon slots: 'in' (incoming nucleon leg, outside the loop), 'int' (internal nucleon), 'out' (outgoing leg),
'pi' (pion line), 'V1'/'V2' (Kroll-Ruderman at the emission/absorption vertex), 'sg' (gamma gamma pi pi seagull).
Nucleon-photon couplings: 'ch' (charge) or 'mag' (anomalous magnetic moment).
order 3: all nucleon-photon couplings 'ch' (LP2010);  order 4: exactly one 'mag' (the O(p^4) kappa insertions).
"""
import itertools, numpy as np
from fractions import Fraction as Fr
from graphs import Mom, Tok, Graph, P_, Q_, QP_, L_, tau, eps3

PH = {'A': dict(mom=Q_, pol='ee'), 'B': dict(mom=-QP_, pol='ep')}   # incoming momenta
Qm = np.diag([1., 0.]).astype(complex)
def Km(which): return np.diag([1., 0.] if which == 'p' else [0., 1.]).astype(complex)

def placements():
    """yield list of (photon, slot, position) with positions ordering photons inside the same slot"""
    slots = ['in', 'int', 'out', 'pi', 'V1', 'V2']
    for sA in slots:
        for sB in slots:
            if sA == sB:
                if sA in ('V1', 'V2'): continue
                yield [('A', sA, 0), ('B', sA, 1)]
                yield [('B', sA, 0), ('A', sA, 1)]
            else:
                yield [('A', sA, 0), ('B', sB, 0)]
    yield [('A', 'sg', 0), ('B', 'sg', 0)]

def classify(pl):
    s = sorted(x[1] for x in pl)
    ext = [x for x in s if x in ('in', 'out')]
    if len(ext) == 2:
        return 'R1' if ext == ['in', 'out'] else None     # both on one leg: external-leg correction
    if len(ext) == 1: return 'R2'
    return '1PI'

def build(pl, cpl, name, final=None, mode='compton'):
    """pl: placement; cpl: dict photon -> 'ch'/'mag' for nucleon-attached photons.
    Returns Graph with .iso_fn(kmat-of-mag) and bookkeeping for renormalisation."""
    g = Graph(name)
    seq = {s: sorted([x for x in pl if x[1] == s], key=lambda z: z[2]) for s in ('in', 'int', 'out', 'pi')}
    at = {s: [x for x in pl if x[1] == s] for s in ('V1', 'V2', 'sg')}
    toks = []                  # built right -> left, reversed at the end
    iso_ops = []               # right -> left: ('N', matrixfn) or ('V1'/'V2', kind)
    n = P_
    pref = 1 + 0j
    def nuc_photon(ph):
        nonlocal pref
        d = PH[ph]
        if cpl[ph] == 'ch':
            toks.append(Tok('V', d['pol'])); pref *= 1j
            iso_ops.append(('N', 'Q'))
        else:
            toks.append(Tok('AG', d['pol'], d['mom'])); pref *= -0.5j
            iso_ops.append(('N', 'K'))
    tree = []
    def treeprop(mom):
        nonlocal pref
        toks.append(Tok('S', mom)); pref *= 1j
        e = mom.ext()
        if e == P_ + Q_: tree.append('ixs')
        elif e == P_ - QP_: tree.append('iys')
        else: raise ValueError(mom)
    # incoming leg (at most one photon there)
    for ph, _, _ in seq['in']:
        nuc_photon(ph); n = n + PH[ph]['mom']; treeprop(n)
    # V1
    if at['V1']:
        ph = at['V1'][0][0]
        toks.append(Tok('V', PH[ph]['pol'], g5=True)); pref *= -0.5j      # -i gA/2F
        n = n + PH[ph]['mom']; iso_ops.append(('V1', 'KR'))
    else:
        toks.append(Tok('K', L_, g5=True)); pref *= 0.5                 # +gA/2F (outgoing pion l)
        iso_ops.append(('V1', 'PV'))
    n = n - L_
    loopN = []
    toks.append(Tok('S', n)); pref *= 1j; loopN.append(n)
    for ph, _, _ in seq['int']:
        nuc_photon(ph); n = n + PH[ph]['mom']
        toks.append(Tok('S', n)); pref *= 1j; loopN.append(n)
    # pion line
    k = L_; loopP = [k]; pscal = []; pion_iso = []
    pref *= 1j
    for ph, _, _ in seq['pi']:
        kout = k + PH[ph]['mom']
        pscal.append((k + kout, PH[ph]['pol'])); pref *= -1
        pion_iso.append('g'); k = kout; loopP.append(k); pref *= 1j
    if at['sg']:
        pscal.append('ee.ep'); pref *= 2j; pion_iso.append('sg')
        k = k + Q_ - QP_; loopP.append(k); pref *= 1j
    # V2
    if at['V2']:
        ph = at['V2'][0][0]
        toks.append(Tok('V', PH[ph]['pol'], g5=True)); pref *= -0.5j
        n = n + PH[ph]['mom']; iso_ops.append(('V2', 'KR'))
    else:
        toks.append(Tok('K', k, g5=True)); pref *= -0.5                 # -gA/2F (incoming pion k)
        iso_ops.append(('V2', 'PV'))
    n = n + k
    # outgoing leg (at most one photon there)
    for ph, _, _ in seq['out']:
        treeprop(n); nuc_photon(ph); n = n + PH[ph]['mom']
    assert n == (P_ + Q_ - QP_ if final is None else final), n
    g.toks = toks[::-1]
    g.scal = pscal
    g.loop = [(m, 'mpi2') for m in loopP] + [(m, '1') for m in loopN[::-1]]
    g.tree = tree
    g.pref = pref * (-1j) * 1j / (16 * np.pi ** 2)   # T = -i x graph;  loop measure i/(16 pi^2) x LoopTools
    g.coup = {'gA': 2, 'F': -2}
    g.iso_ops, g.pion_iso, g.cpl, g.pl = iso_ops, pion_iso, cpl, pl
    g.kind = classify(pl) if mode == 'compton' else mode
    g.mode = mode
    return g

def isospin(g, kmag):
    """<p| product of nucleon isospin operators |p> summed over pion indices; kmag = matrix of kappa"""
    tot = 0j
    for a in range(3):
        for b in range(3):
            # pion line tensor from index a (V1) to b (V2)
            P = np.eye(3)
            for x in g.pion_iso:
                if x == 'g': P = P @ np.array([[eps3(i, j) for j in range(3)] for i in range(3)])
                else: P = P @ (np.eye(3) - np.diag([0, 0, 1.]))
            pt = P[a, b]
            if pt == 0: continue
            Mt = np.eye(2, dtype=complex)        # accumulate right -> left: Mt = op @ Mt
            for kind, what in g.iso_ops:
                if kind == 'N':
                    op = Qm if what == 'Q' else kmag
                elif kind == 'V1':
                    op = tau[a] if what == 'PV' else sum(eps3(a, c) * tau[c] for c in range(3))
                else:
                    op = tau[b] if what == 'PV' else sum(eps3(b, c) * tau[c] for c in range(3))
                Mt = op @ Mt
            tot += pt * Mt[0, 0]
    return tot

def all_graphs(order):
    gs = []
    for pl in placements():
        kind = classify(pl)
        if kind is None: continue
        nucph = [ph for ph, s, _ in pl if s in ('in', 'int', 'out')]
        if order == 3:
            cpls = [{ph: 'ch' for ph in nucph}]
        else:
            cpls = []
            for m in nucph:
                c = {ph: 'ch' for ph in nucph}; c[m] = 'mag'; cpls.append(c)
        for c in cpls:
            tag = '_'.join(f'{ph}{s}{i}' + (c.get(ph, '')[:1] if s in ('in', 'int', 'out') else '') for ph, s, i in pl)
            g = build(pl, c, f'o{order}_{tag}')
            if order == 3:
                g.iso = {'1': isospin(g, None)}
            else:
                g.iso = {'kp': isospin(g, Km('p')), 'kn': isospin(g, Km('n'))}
            if all(abs(v) < 1e-14 for v in g.iso.values()): continue
            gs.append(g)
    return gs

if __name__ == '__main__':
    for o in (3, 4):
        gs = all_graphs(o)
        print(f'order {o}: {len(gs)} graphs')
        for g in gs:
            print(f'  {g.name:40s} {g.kind:4s} loop={len(g.loop)} tree={g.tree} iso={ {k: round(v.real, 6) for k, v in g.iso.items()} }')


# ------------------------------------------------------------ renormalisation pieces
def vertex_graph(slot, cpl, name):
    """on-shell gamma N N vertex loop, photon A (momentum q, q^2=-dl) at slot in ('int','pi','V1','V2');
    value Lambda = product of vertex factors, propagator i's and loop measure (no -i of T)"""
    pl = [('A', slot, 0)]
    g = build(pl, {'A': cpl}, name, final=P_ + Q_, mode='vertex')
    g.pref *= 1j
    return g

def self_graph(name):
    g = build([], {}, name, final=P_, mode='self')
    g.pref *= 1j
    return g

def vertex_iso(g, kmag):
    return isospin(g, kmag)

def sub_graph_R2(g):
    """tree graph replacing the vertex loop of an R2 graph by F1s ee/ + F2s ag(ee,k) (photon I inside)"""
    (O, sO, _), = [x for x in g.pl if x[1] in ('in', 'out')]
    (I, sI, _), = [x for x in g.pl if x[1] not in ('in', 'out')]
    h = Graph(g.name + '_sub')
    dO, dI = PH[O], PH[I]
    veff = Tok('VEFF', dI['pol'], dI['mom'])
    if g.cpl[O] == 'ch': vO = Tok('V', dO['pol']); fO = 1j
    else: vO = Tok('AG', dO['pol'], dO['mom']); fO = -0.5j
    if sO == 'in':
        h.toks = [veff, Tok('S', P_ + dO['mom']), vO]; h.tree = ['ixs' if O == 'A' else 'iys']
    else:
        h.toks = [vO, Tok('S', P_ + dI['mom']), veff]; h.tree = ['ixs' if I == 'A' else 'iys']
    h.pref = fO * 1j * (-1j)
    h.vkey = (sI, g.cpl.get(I, None) if sI == 'int' else None)
    h.Ocpl = g.cpl[O]
    return h

def sub_graph_R1(g):
    (O1, s1, _), = [x for x in g.pl if x[1] == 'in']
    (O2, s2, _), = [x for x in g.pl if x[1] == 'out']
    h = Graph(g.name + '_sub')
    ps = P_ + PH[O1]['mom']
    def vtx(O):
        d = PH[O]
        if g.cpl[O] == 'ch': return Tok('V', d['pol']), 1j
        return Tok('AG', d['pol'], d['mom']), -0.5j
    v1, f1 = vtx(O1); v2, f2 = vtx(O2)
    h.toks = [v2, Tok('S', ps), Tok('LINS', ps), Tok('S', ps), v1]
    h.tree = ['ixs', 'ixs'] if O1 == 'A' else ['iys', 'iys']
    h.pref = f1 * f2 * 1j * 1j * (-1j)
    h.cpls = (g.cpl[O1], g.cpl[O2])
    return h
