"""delta.py -- one-loop (and tree) Compton graphs with Delta(1232) lines, covariant BChPT (LP2010, PVY2007 couplings).

Lagrangian (e removed; code conventions of graphs.py: D_mu N = (d - i Q A) N, D_mu pi^a = d pi^a - eps^{3ab} A pi^b):
  L_Delta  = Dbar_mu (i gamma^{mu nu al} D_al - MD gamma^{mu nu}) D_nu,  D_al Delta = (d - i Q_D A) Delta
  L_piND   = (i hA/(2 F MD)) Nbar T^a gamma^{mu nu lam} (D_mu Delta_nu) (D_lam pi^a) + h.c.   (consistent coupling)
  L_gND    = (3 i gM/(2M(M+MD))) Nbar T^3 (D_mu Delta_nu) Ftilde^{mu nu}
           - (3 gE/(2M(M+MD))) Nbar T^3 gamma5 (D_mu Delta_nu) F^{mu nu} + h.c.,   Ftilde^{mu nu} = eps^{mu nu rho lam} d_rho A_lam
Feynman rules (k incoming photon momentum, polarisation e; pion momenta incoming at the vertex):
  Delta propagator   i N^{ab}(p)/(p^2-MD^2), N = D-dim RS numerator (dgraphs.DTok 'RS'; 'RSR' between two transverse vertices)
  gamma Delta Delta  i Q_D gamma^{mu nu al} e_al                         (mu: Deltabar index, nu: Delta index)
  pi N Delta         c gamma^{mu nu lam} p_mu k_lam T^a          (Delta in, p = its momentum; c = hA/(2F MD))
  pi Delta N         -c gamma^{mu nu lam} p'_mu k_lam T^a+       (Delta out, p' = its momentum)
    contacts: k_lam -> -i eps^{3 a b} e_lam (pion charge, b = pion line index);  p -> Q_D e (Delta charge)
  gamma N Delta      cM eps^{mu nu rho lam} p_mu k_rho e_lam T^3 + i cE gamma5 [(p.k) e^nu - (p.e) k^nu] T^3  (and T^3+ for
                     Delta out with p -> p'); cM,E = 3 gM,E/(2M(M+MD)); contacts p -> Q_D e (gamma gamma N Delta)
  eps^{mu nu rho lam} = SEPS * i gamma^{[mu nu rho lam]} gamma5   (4-dim identity, defines eps in D dims; SEPS fixed by the
                     gM gE interference of the Delta-pole graph, tests/t_dpole.py)
Isospin: T^a (2x4), T^a T^b+ = 2/3 delta^{ab} - i/3 eps^{abc} tau^c;  T^a Q_D - Q_N T^a = -i eps^{3ab} T^b.
"""
import itertools, numpy as np
from fractions import Fraction as Fr
from graphs import Mom, P_, Q_, QP_, L_, PP_, tau, eps3
from dgraphs import DTok, DGraph

SEPS = -1
SGND = -1        # overall sign of the gamma N Delta couplings relative to the Lagrangian as transcribed: fixed by
                 # the physical (constructive, Watson) interference of the Born and Delta M1+(3/2) amplitudes below
                 # the resonance and the one-loop Im gM, Im gE of the code (PV2006); only gA hA gM is physical
TREE4 = True     # Delta propagators outside the loop (Delta-pole graphs) in 4 dimensions (1/(D-1) -> 1/3)
PH = {'A': (Q_, 'ee'), 'B': (-QP_, 'ep')}
_s6, _s3 = np.sqrt(6), np.sqrt(3)
T = [np.array([[-_s3, 0, 1, 0], [0, -1, 0, _s3]], complex) / _s6,
     -1j * np.array([[_s3, 0, 1, 0], [0, 1, 0, _s3]], complex) / _s6,
     np.sqrt(2 / 3) * np.array([[0, 1, 0, 0], [0, 0, 1, 0]], complex)]
Td = [t.conj().T for t in T]
QN = np.diag([1., 0.]).astype(complex)
QD = np.diag([2., 1., 0., -1.]).astype(complex)
E3 = np.array([[eps3(i, j) for j in range(3)] for i in range(3)], complex)
Pm = np.eye(3) - np.diag([0, 0, 1.])
KVAL = dict(ch=1, chD=1, mag=2, gM=2, gE=2, gg=3, PV=1, KR=1, pND=1, cpi=1, cD=2, pipi=2, sg=2)

# ------------------------------------------------------------------ configurations
def configs():
    """photon placements: dict(in, int, out, pi: lists; V1, V2: photon or None; sg: bool)"""
    slots = ['in', 'int', 'out', 'pi', 'V1', 'V2']
    def mk(**kw):
        c = dict(**{'in': [], 'int': [], 'out': [], 'pi': []}, V1=None, V2=None, sg=False)
        for k, v in kw.items(): c[k] = v
        return c
    for sA in slots:
        for sB in slots:
            if sA == sB:
                if sA in ('V1', 'V2'): continue
                yield mk(**{sA: ['A', 'B']}); yield mk(**{sA: ['B', 'A']})
            else:
                yield mk(**{sA: 'A' if sA in ('V1', 'V2') else ['A'], sB: 'B' if sB in ('V1', 'V2') else ['B']})
    yield mk(sg=True)
    for w in ('in', 'int', 'out'):
        yield mk(**{w: ['AB']})

class Bad(Exception): pass

def build(cfg, types, cpl, mode='compton', final=None, name=None):
    """cfg: placement; types: tuple of segment types ('N'/'D') in the order in-tree, loop, out-tree;
    cpl: dict for baryon photon vertices / contacts, keys (slot, i) -> 'ch'|'gM'|'gE'|'cpi'|'cD'|'gg1'|'gg2'.
    Returns DGraph or raises Bad."""
    nin, nint, nout = len(cfg['in']), len(cfg['int']), len(cfg['out'])
    tin, tloop, tout = types[:nin], types[nin:nin + nint + 1], types[nin + nint + 1:]
    assert len(tout) == nout
    g = DGraph(name or 'dg')
    toks, pref, iso_ops = [], 1 + 0j, []
    cnt = dict(N=0, D=0, pi=0, k=0, gA=0, hA=0, gM=0, gE=0)
    tree, loopB = [], []
    state = dict(n=P_, open=None, idx=0, last_tr=None)
    dprops = []                          # [token-position, left-vertex-transverse, right-vertex-transverse]

    def newidx():
        state['idx'] += 1
        return f'a{state["idx"]}'
    def prop(typ, inloop):
        nonlocal pref
        n = state['n']
        if typ == 'N':
            toks.append(DTok('S', n)); pref *= 1j; cnt['N'] += 1
            if inloop: loopB.append((n, '1'))
            else: tree.append(_treeden(n, 'N'))
        else:
            iR = state['open_new']; iL = newidx()
            toks.append(DTok('RS', n, iL, iR)); pref *= 1j; cnt['D'] += 1
            dprops.append([len(toks) - 1, None, state['last_tr'], inloop])
            state['open'] = iL
            if inloop: loopB.append((n, 'MD2'))
            else: tree.append(_treeden(n, 'D'))
    def vtx_done(transverse_left):
        # mark the Delta propagator on the right of this vertex (if any) with this vertex's transversality
        if dprops and dprops[-1][1] is None and dprops[-1][0] == len(toks) - 2:
            dprops[-1][1] = transverse_left
    def photon_vertex(ph, X, Y, c):
        """baryon photon vertex, X: type on the right, Y: type on the left"""
        nonlocal pref
        k, e = PH[ph]
        n0 = state['n']; n1 = n0 + k
        if X == 'N' and Y == 'N':
            if c != 'ch': raise Bad
            toks.append(DTok('V', e)); pref *= 1j; iso_ops.append(('M', QN)); cnt['k'] += 1; tr = None
        elif X == 'D' and Y == 'D':
            if c != 'ch': raise Bad
            iL = newidx()
            toks.append(DTok('G3E', iL, state['open'], e)); pref *= 1j; iso_ops.append(('M', QD)); cnt['k'] += 1
            tr = False; state['open_new'] = iL
        elif X == 'N' and Y == 'D':
            if c not in ('gM', 'gE'): raise Bad
            iL = newidx()
            if c == 'gM': toks.append(DTok('EPS4', n1, iL, k, e, g5=True)); pref *= 1.5 * SEPS * 1j * SGND
            else: toks.append(DTok('SCI', n1, k, e, iL, g5=True)); pref *= 1.5j * SGND
            iso_ops.append(('M', Td[2])); cnt['k'] += 2; cnt[c] += 1; tr = True; state['open_new'] = iL
        else:
            if c not in ('gM', 'gE'): raise Bad
            if c == 'gM': toks.append(DTok('EPS4', n0, state['open'], k, e, g5=True)); pref *= 1.5 * SEPS * 1j * SGND
            else: toks.append(DTok('SCI', n0, k, e, state['open'], g5=True)); pref *= 1.5j * SGND
            iso_ops.append(('M', T[2])); cnt['k'] += 2; cnt[c] += 1; tr = True
        vtx_done(tr)
        state['n'] = n1; state['last_tr'] = tr
    def gg_vertex(X, Y, c):
        """gamma gamma N Delta contact (from D_mu Delta_nu in L_gND); c = 'gM1','gM2','gE1','gE2'
        (1: photon A at the Delta derivative, B in F; 2: the reverse)"""
        nonlocal pref
        (kA, eA), (kB, eB) = PH['A'], PH['B']
        coup, which = c[:2], c[2]
        ed, kf, ef = (eA, kB, eB) if which == '1' else (eB, kA, eA)
        if X == 'N' and Y == 'D':
            iL = newidx(); ix = iL; state['open_new'] = iL; op = QD @ Td[2]
        elif X == 'D' and Y == 'N':
            ix = state['open']; op = T[2] @ QD
        else: raise Bad
        if coup == 'gM': toks.append(DTok('EPS4', ed, ix, kf, ef, g5=True)); pref *= 1.5 * SEPS * 1j * SGND
        else: toks.append(DTok('SCI', ed, kf, ef, ix, g5=True)); pref *= 1.5j * SGND
        iso_ops.append(('M', op)); cnt['k'] += 3; cnt[coup] += 1
        vtx_done(False)
        state['n'] = state['n'] + kA + kB; state['last_tr'] = False
    def pion_vertex(which, X, Y, ph, c, kpi):
        """which: 'V1' (emission, outgoing pion l) or 'V2' (absorption, incoming pion kpi)"""
        nonlocal pref
        n0 = state['n']
        kin = -L_ if which == 'V1' else kpi
        dn = -L_ if which == 'V1' else kpi
        if ph is not None: dn = dn + PH[ph][0]
        n1 = n0 + dn
        e = PH[ph][1] if ph else None
        if X == 'N' and Y == 'N':
            if ph is None:
                toks.append(DTok('K', kin, g5=True)); pref *= -0.5          # -(gA/2F) k_in/ g5
                iso_ops.append((which, 'tau'))
            else:
                if c is not None: raise Bad
                toks.append(DTok('V', e, g5=True)); pref *= -0.5j
                iso_ops.append((which, 'tauKR'))
            cnt['k'] += 1; cnt['gA'] += 1; tr = None
        elif X == 'N' and Y == 'D':                                       # Deltabar <- N,  -c G3(p', nu, k_in) T+
            iL = newidx(); state['open_new'] = iL
            if ph is None:
                toks.append(DTok('G3', n1, iL, kin)); pref *= -0.5; iso_ops.append((which, 'Td')); tr = True; kk = 1
            elif c == 'cpi':
                toks.append(DTok('G3', n1, iL, e)); pref *= 0.5j; iso_ops.append((which, 'TdKR')); tr = True; kk = 1
            elif c == 'cD':
                toks.append(DTok('G3', e, iL, kin)); pref *= -0.5; iso_ops.append((which, 'QTd')); tr = False; kk = 2
            else: raise Bad
            cnt['k'] += kk; cnt['hA'] += 1
        elif X == 'D' and Y == 'N':                                       # Nbar <- Delta,  c G3(p, nu, k_in) T
            ix = state['open']
            if ph is None:
                toks.append(DTok('G3', n0, ix, kin)); pref *= 0.5; iso_ops.append((which, 'T')); tr = True; kk = 1
            elif c == 'cpi':
                toks.append(DTok('G3', n0, ix, e)); pref *= -0.5j; iso_ops.append((which, 'TKR')); tr = True; kk = 1
            elif c == 'cD':
                toks.append(DTok('G3', e, ix, kin)); pref *= 0.5; iso_ops.append((which, 'TQ')); tr = False; kk = 2
            else: raise Bad
            cnt['k'] += kk; cnt['hA'] += 1
        else: raise Bad
        vtx_done(tr)
        state['n'] = n1; state['last_tr'] = tr

    # ---------------- incoming leg
    X = 'N'
    for i, el in enumerate(cfg['in']):
        Y = tin[i]
        if el == 'AB': gg_vertex(X, Y, cpl[('in', i)])
        else: photon_vertex(el, X, Y, cpl.get(('in', i), 'ch'))
        prop(Y, False); X = Y
    # ---------------- V1
    pion_vertex('V1', X, tloop[0], cfg['V1'], cpl.get(('V1',)), None)
    prop(tloop[0], True); X = tloop[0]
    for i, el in enumerate(cfg['int']):
        Y = tloop[i + 1]
        if el == 'AB': gg_vertex(X, Y, cpl[('int', i)])
        else: photon_vertex(el, X, Y, cpl.get(('int', i), 'ch'))
        prop(Y, True); X = Y
    # ---------------- pion line
    kp = L_; loopP = [kp]; pscal = []; piso = []
    pref *= 1j; cnt['pi'] += 1
    for ph in cfg['pi']:
        kout = kp + PH[ph][0]
        pscal.append((kp + kout, PH[ph][1])); pref *= -1; piso.append('g'); cnt['k'] += 2
        kp = kout; loopP.append(kp); pref *= 1j; cnt['pi'] += 1
    if cfg['sg']:
        pscal.append('ee.ep'); pref *= 2j; piso.append('sg'); cnt['k'] += 2
        kp = kp + Q_ - QP_; loopP.append(kp); pref *= 1j; cnt['pi'] += 1
    # ---------------- V2
    Y = tout[0] if nout else 'N'
    pion_vertex('V2', X, Y, cfg['V2'], cpl.get(('V2',)), kp)
    X = Y
    for i, el in enumerate(cfg['out']):
        prop(X, False)
        Y = tout[i + 1] if i + 1 < nout else 'N'
        if el == 'AB': gg_vertex(X, Y, cpl[('out', i)])
        else: photon_vertex(el, X, Y, cpl.get(('out', i), 'ch'))
        X = Y
    if X != 'N': raise Bad
    if state['n'] != (PP_ if final is None else final): raise AssertionError(state['n'])
    # Delta propagators between two transverse vertices: reduced numerator
    for pos, trL, trR, inloop in dprops:
        tk = toks[pos]
        sfx = '' if inloop or not TREE4 else '4'
        if trL and trR: toks[pos] = DTok('RSR' + sfx, *tk.args)
        elif trL: toks[pos] = DTok('RSH' + sfx, *tk.args, 'L')
        elif trR: toks[pos] = DTok('RSH' + sfx, *tk.args, 'R')
        elif sfx: toks[pos] = DTok('RS' + sfx, *tk.args)
    g.toks = toks[::-1]
    g.scal = pscal
    g.loop = [(m, 'mpi2') for m in loopP] + [(m, s) for m, s in loopB[::-1]]
    g.tree = tree
    g.pref = pref * (-1j) * 1j / (16 * np.pi ** 2)
    g.coup = {'gA': cnt['gA'], 'hA': cnt['hA'], 'F': -(cnt['gA'] + cnt['hA']), 'MD': -cnt['hA'],
              'gM': cnt['gM'], 'gE': cnt['gE'], 'iMp': cnt['gM'] + cnt['gE']}
    g.coup = {k: v for k, v in g.coup.items() if v}
    g.order = 4 + cnt['k'] - 2 * cnt['pi'] - cnt['N'] - Fr(cnt['D'], 2)
    g.iso_ops, g.piso, g.cfg, g.types, g.cpl = iso_ops, piso, cfg, types, cpl
    g.nD = cnt['D']; g.ngND = cnt['gM'] + cnt['gE']
    g.mode = mode
    g.kind = classify(cfg) if mode == 'compton' else mode
    return g

def _treeden(n, typ):
    e = n.ext()
    if e == P_ + Q_: return 'ixs' if typ == 'N' else 'ixD'
    if e == P_ - QP_: return 'iys' if typ == 'N' else 'iyD'
    if e in (P_, PP_, P_ + Q_ - QP_):
        if typ == 'N': raise Bad          # external-leg self-energy
        return 'iMD'
    raise ValueError(n)

def classify(cfg):
    nin, nout = len(cfg['in']), len(cfg['out'])
    inside = cfg['int'] or cfg['pi'] or cfg['V1'] or cfg['V2'] or cfg['sg']
    if not nin and not nout: return '1PI'
    if nin and nout: return 'R1'
    if inside: return 'R2'
    return 'RLEG'

def isospin(g):
    """<p| ... |p> summed over the pion indices"""
    tot = 0j
    for a in range(3):
        for b in range(3):
            P = np.eye(3, dtype=complex)
            for x in g.piso:
                P = P @ (E3 if x == 'g' else Pm)
            pt = P[a, b]
            if pt == 0: continue
            Mt = np.eye(2, dtype=complex)
            for kind, what in g.iso_ops:
                if kind == 'M': op = what
                else:
                    j = a if kind == 'V1' else b
                    op = {'tau': lambda: tau[j], 'tauKR': lambda: sum(eps3(j, c) * tau[c] for c in range(3)),
                          'T': lambda: T[j], 'Td': lambda: Td[j],
                          'TKR': lambda: sum(eps3(c, j) * T[c] for c in range(3)),
                          'TdKR': lambda: sum(eps3(c, j) * Td[c] for c in range(3)),
                          'TQ': lambda: T[j] @ QD, 'QTd': lambda: QD @ Td[j]}[what]()
                Mt = op @ Mt
            tot += pt * Mt[0, 0]
    return tot

# ------------------------------------------------------------------ enumeration
def all_graphs(maxgND=1, minD=1, maxorder=None, exclude_rleg_N=True):
    """all one-loop graphs with >= minD Delta propagators and <= maxgND gamma N Delta couplings; no kappa"""
    out, seen = [], set()
    for cfg in configs():
        nin, nint, nout = len(cfg['in']), len(cfg['int']), len(cfg['out'])
        nseg = nin + nint + 1 + nout
        # coupling options per baryon vertex / contact
        keys = [('in', i) for i in range(nin)] + [('int', i) for i in range(nint)] + [('out', i) for i in range(nout)]
        vkeys = [k for k in ('V1', 'V2') if cfg[k]]
        for types in itertools.product('ND', repeat=nseg):
            if types.count('D') < minD: continue
            opts = []
            for key in keys:
                el = cfg[key[0]][key[1]]
                opts.append(['gM1', 'gM2', 'gE1', 'gE2'] if el == 'AB' else ['ch', 'gM', 'gE'])
            for key in vkeys:
                opts.append(['cpi', 'cD', None])
            for combo in itertools.product(*opts):
                cpl = {}
                for key, c in zip(keys, combo[:len(keys)]): cpl[key] = c
                for key, c in zip(vkeys, combo[len(keys):]): cpl[(key,)] = c
                try:
                    g = build(cfg, types, cpl)
                except Bad:
                    continue
                if g.ngND > maxgND: continue
                if maxorder is not None and g.order > maxorder: continue
                tag = _tag(cfg, types, cpl)
                if tag in seen: continue
                seen.add(tag)
                g.name = 'd_' + tag
                g.iso = {'1': isospin(g)}
                if abs(g.iso['1']) < 1e-14: continue
                out.append(g)
    return out

def _tag(cfg, types, cpl):
    s = []
    for w in ('in', 'int', 'out', 'pi'):
        if cfg[w]: s.append(w + ''.join(cfg[w]))
    for w in ('V1', 'V2'):
        if cfg[w]: s.append(w + cfg[w] + str(cpl.get((w,)))[:3])
    if cfg['sg']: s.append('sg')
    s.append(''.join(types))
    s.append(''.join(str(v)[:3] for k, v in sorted(cpl.items(), key=str) if v and k[0] in ('in', 'int', 'out')))
    return '_'.join(x for x in s if x)

if __name__ == '__main__':
    import collections
    gs = all_graphs()
    print(len(gs), 'graphs')
    by = collections.Counter((g.kind, str(g.order), g.ngND) for g in gs)
    for k, v in sorted(by.items()): print(k, v)

# ------------------------------------------------------------------ Delta-pole tree graphs
def dpole_graphs(c1, c2, tag=''):
    """tree Delta-pole graphs (s and u channel), coupling c1 ('gM'/'gE') at the incoming-photon vertex, c2 at the
    outgoing one; coup gM, gE powers, isospin T3 T3+ = 2/3 (proton).  Reproduces s_Delta2 of the code (t_dpole.py)."""
    from dgraphs import DTok, DGraph
    def vtx(c, pd, idx, k, e):
        if c == 'gM': return DTok('EPS4', pd, idx, k, e, g5=True), 1.5 * SEPS * 1j * SGND
        return DTok('SCI', pd, k, e, idx, g5=True), 1.5j * SGND
    gs = []
    pd = P_ + Q_
    vB, fB = vtx(c2, pd, 'a1', -QP_, 'ep'); vA, fA = vtx(c1, pd, 'a2', Q_, 'ee')
    g = DGraph(f'dpole{tag}_s_{c1}{c2}'); g.toks = [vB, DTok('RSR', pd, 'a1', 'a2'), vA]; g.tree = ['ixD']
    g.pref = fA * fB * 1j * (-1j); gs.append(g)
    pd = P_ - QP_
    vA, fA = vtx(c1, pd, 'a1', Q_, 'ee'); vB, fB = vtx(c2, pd, 'a2', -QP_, 'ep')
    g = DGraph(f'dpole{tag}_u_{c1}{c2}'); g.toks = [vA, DTok('RSR', pd, 'a1', 'a2'), vB]; g.tree = ['iyD']
    g.pref = fA * fB * 1j * (-1j); gs.append(g)
    for g in gs:
        g.coup = {c1: 1, c2: 1, 'iMp': 2} if c1 != c2 else {c1: 2, 'iMp': 2}
        g.iso = {'1': (T[2] @ Td[2])[0, 0]}
        g.kind = '1PI'; g.order = Fr(7, 2)
    return gs
