"""dassemble.py -- renormalised Compton amplitude of a set of Delta-line graphs (delta.py), Python check path.

On-shell renormalisation as for the piN loops (assemble.py): reducible graphs whose loop is a gamma N N vertex (R2) get
the subtraction F1(0) gamma + F2(0) magnetic, those with a nucleon self-energy (R1) Sigma(M) + (p/-M) Sigma'(M).
R2 graphs whose loop is a gamma N -> Delta vertex (Delta tree propagator next to the loop) are left in MS-bar here.
"""
import numpy as np, kin, delta, devalamp as dv, evalamp as ev, runform, formtools as ft
from graphs import Graph, Tok, P_, Q_, QP_
from fractions import Fraction as Fr

def _strip_legs(g):
    """configuration of the loop of an R2/R1 graph: (cfg with the inner photon called 'A' or no photon, loop types, cpl)"""
    cfg, types, cpl = g.cfg, g.types, g.cpl
    nin, nint = len(cfg['in']), len(cfg['int'])
    loop_types = types[nin:nin + nint + 1]
    legs = [('in', cfg['in']), ('out', cfg['out'])]
    inner = None
    for w in ('int', 'pi'):
        if cfg[w]: inner = cfg[w][0]
    for w in ('V1', 'V2'):
        if cfg[w]: inner = cfg[w]
    if cfg['sg']: inner = 'AB'
    ren = (lambda z: z) if inner in (None, 'A') else (lambda z: {'B': 'A', 'A': 'B'}.get(z, z))
    c = dict(**{'in': [], 'out': []}, int=[ren(z) for z in cfg['int']], pi=[ren(z) for z in cfg['pi']],
             V1=ren(cfg['V1']) if cfg['V1'] else None, V2=ren(cfg['V2']) if cfg['V2'] else None, sg=False)
    cp = {}
    for k, v in cpl.items():
        if k[0] == 'int': cp[k] = v
        if k[0] in ('V1', 'V2'): cp[k] = v
    tin = types[nin - 1] if nin else None
    tout = types[nin + nint + 1] if cfg['out'] else None
    return c, loop_types, cp, tin, tout

class DModel:
    def __init__(self, graphs, cs, nproc=8, iso_override=None, tag=''):
        """graphs: delta.py Compton graphs (1PI, R2, R1); iso_override: {kind: factor} multiplying the isospin factors
        (kinds 1PI, R2, R1, vertex, self; used to reproduce the 'corrected isospin' of the NNLO code)"""
        self.cs, self.iso_override = cs, iso_override
        if iso_override is not None:              # the override rescales g.iso: work on copies (graphs are shared)
            import copy
            graphs = [copy.copy(g) for g in graphs]
        self.graphs = graphs
        vert, selfe, subs = {}, {}, []
        for g in graphs:
            if g.kind == 'R2':
                c, lt, cp, tin, tout = _strip_legs(g)
                legt = tin if g.cfg['in'] else tout
                if legt != 'N': continue                                   # gamma N -> Delta vertex: MS-bar
                key = repr((c, lt, sorted(cp.items(), key=str)))
                if key not in vert:
                    v = delta.build(c, lt, cp, mode='vertex', final=P_ + Q_, name=f'dv{tag}_{len(vert)}')
                    v.iso = {'1': delta.isospin(v)}; v.pref *= 1j
                    vert[key] = v
                h = self._sub_R2(g); h.vkey = key; subs.append(h)
            elif g.kind == 'R1':
                c, lt, cp, tin, tout = _strip_legs(g)
                if tin != 'N' or tout != 'N': continue
                key = repr(lt)
                if key not in selfe:
                    c0 = dict(**{'in': [], 'int': [], 'out': [], 'pi': []}, V1=None, V2=None, sg=False)
                    s = delta.build(c0, lt, {}, mode='self', final=P_, name=f'ds{tag}_{len(selfe)}')
                    s.iso = {'1': delta.isospin(s)}; s.pref *= 1j
                    selfe[key] = s
                h = self._sub_R1(g); h.skey = key; subs.append(h)
        self.vert, self.selfe, self.subs = vert, selfe, subs
        loops = graphs + list(vert.values()) + list(selfe.values())
        self.dc = dv.compile_graphs(loops, nproc)
        res = runform.project_all(subs, nproc)
        self.sc = {h.name: ev.Compiled(h, res[h.name]) for h in subs}
        if iso_override is not None:          # dict kind -> factor multiplying the true isospin factor
            for c in list(self.dc.values()):
                c.g.iso = {'1': c.g.iso['1'] * iso_override[c.g.kind]}
        self.vc = {k: self.vconsts(self.dc[v.name]) for k, v in vert.items()}
        self.vcP = {k: self.vconsts(self.dc[v.name], pole=True) for k, v in vert.items()}
        self.sfc = {k: self.sconsts(self.dc[s.name]) for k, s in selfe.items()}
        self.sfcP = {k: self.sconsts(self.dc[s.name], pole=True) for k, s in selfe.items()}

    def _sub_R2(self, g):
        cfg = g.cfg
        sO = 'in' if cfg['in'] else 'out'
        O = (cfg['in'] or cfg['out'])[0]; I = 'B' if O == 'A' else 'A'
        h = Graph(g.name + '_sub')
        dO, dI = delta.PH[O], delta.PH[I]
        veff = Tok('VEFF', dI[1], dI[0]); vO = Tok('V', dO[1])
        if sO == 'in':
            h.toks = [veff, Tok('S', P_ + dO[0]), vO]; h.tree = ['ixs' if O == 'A' else 'iys']
        else:
            h.toks = [vO, Tok('S', P_ + dI[0]), veff]; h.tree = ['ixs' if I == 'A' else 'iys']
        h.pref = 1j * 1j * (-1j); h.coup = {}; h.iso = {'1': 1.0}
        return h

    def _sub_R1(self, g):
        O1 = g.cfg['in'][0]; O2 = g.cfg['out'][0]
        h = Graph(g.name + '_sub')
        ps = P_ + delta.PH[O1][0]
        h.toks = [Tok('V', delta.PH[O2][1]), Tok('S', ps), Tok('LINS', ps), Tok('S', ps), Tok('V', delta.PH[O1][1])]
        h.tree = ['ixs', 'ixs'] if O1 == 'A' else ['iys', 'iys']
        h.pref = 1j * 1j * 1j * 1j * (-1j); h.coup = {}; h.iso = {'1': 1.0}
        return h

    def vconsts(self, cv, dls=(1e-4, 2e-4, 3e-4), pole=False):      # quadratic extrapolation (evalamp.vertex_consts)
        out = []
        for dl in dls:
            en = dict(mpi2=self.cs['mpi2'], dl=dl, M=1.0, MD=self.cs['MD'])
            W = cv.raw(en, pole) * cv.coup(self.cs) * cv.isof(self.cs)
            A, B, C = np.linalg.solve(ev.vgram(en), W)
            out.append(np.array([A + 2 * B, B]))
        return 3 * out[0] - 3 * out[1] + out[2]

    def sconsts(self, cv, h=1e-4, pole=False):
        def ab(ss):
            en = dict(mpi2=self.cs['mpi2'], ss=ss, MD=self.cs['MD'])
            Y = cv.raw(en, pole) * cv.coup(self.cs) * cv.isof(self.cs)
            return Y[1] / ss, Y[0]
        a0, b0 = ab(1.0)
        (a1, b1), (a2, b2), (am1, bm1), (am2, bm2) = ab(1 + h), ab(1 + 2 * h), ab(1 - h), ab(1 - 2 * h)
        return a0 + b0, a0 + 2 * ((8 * (a1 - am1) - (a2 - am2) + 8 * (b1 - bm1) - (b2 - bm2)) / (12 * h))

    def V(self, x, t, pole=False, names=None):
        cs = self.cs
        en = dv.env(x, t, cs)
        V = np.zeros(6, complex); by = {}
        for g in self.graphs:
            if names is not None and g.name not in names: continue
            c = self.dc[g.name]
            v = c.raw(en, pole) * c.coup(cs) * c.isof(cs)
            V += v; by[g.name] = v
        for h in self.subs:
            if names is not None and h.name[:-4] not in names: continue
            c = self.sc[h.name]
            e = dict(en)
            if hasattr(h, 'vkey'):
                c1, c2 = (self.vcP if pole else self.vc)[h.vkey]; e.update(F1s=c1, F2s=c2)
            else:
                L0, L1 = (self.sfcP if pole else self.sfc)[h.skey]; e.update(L0s=L0, L1s=L1)
            f = -1.0
            v = c.raw(e, pole) * c.coup(cs) * f
            V += v; by[h.name] = v
        return V, by

    def A(self, x, t, pole=False, names=None):
        en = dv.env(x, t, self.cs)
        V, _ = self.V(x, t, pole, names)
        return np.linalg.solve(ev.gram(en), V)

    # ------------------------------------------------------------------ gauge invariance (ee -> q)
    def gi_compile(self, nproc=8):
        self.gdc = dv.compile_graphs(self.graphs, nproc, gi=True)
        res = runform.project_all(self.subs, nproc, gi=True)
        self.gsc = {h.name: [compile(ft.to_python(d.get(f'V{j}', '0')), f'W{j}', 'eval') for j in range(1, 9)]
                    for h in self.subs for d in [res[h.name]]}

    def GI(self, x, t):
        cs = self.cs
        en = dv.env(x, t, cs)
        W = np.zeros(8, complex); big = 0.
        for g in self.graphs:
            c = self.gdc[g.name]
            w = c.raw(en) * c.coup(cs) * self.dc[g.name].isof(cs)
            W += w; big = max(big, np.max(abs(w)))
        for h in self.subs:
            c = self.sc[h.name]
            T, U = c.tensors(en)
            e = dict(en); e['T'] = T; e['U'] = U
            if hasattr(h, 'vkey'):
                c1, c2 = self.vc[h.vkey]; e.update(F1s=c1, F2s=c2)
            else:
                L0, L1 = self.sfc[h.skey]; e.update(L0s=L0, L1s=L1)
            w = -np.array([eval(cc, e) for cc in self.gsc[h.name]]) * c.coup(cs)
            W += w; big = max(big, np.max(abs(w)))
        return W, big
