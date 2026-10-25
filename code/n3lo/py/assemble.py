"""assemble.py -- renormalised one-loop Compton amplitude from the generated graphs (Python check path).

On-mass-shell scheme of LP2010 Sect. 3.2: reducible graphs with Sigma_R = Sigma - Sigma(M) - (p/-M) Sigma'(M) and
Gamma_R = Gamma - [F1(0) gamma + F2(0) magnetic], i.e. physical mass, charge and magnetic moment in the Born graphs.
"""
import numpy as np
import piN, evalamp as ev
from graphs import Graph

class Model:
    def __init__(self, order, cs, nproc=8, extra_graphs=()):
        self.order, self.cs = order, cs
        comp = piN.all_graphs(order)
        subs = []
        for g in comp:
            if g.kind == 'R2':
                h = piN.sub_graph_R2(g); h.parent = g; subs.append(h)
            elif g.kind == 'R1':
                h = piN.sub_graph_R1(g); h.parent = g; subs.append(h)
        for h in subs:
            h.iso = {'1': 1.0}; h.coup = {}
        # vertex loops needed
        vkeys = sorted({h.vkey for h in subs if hasattr(h, 'vkey')}, key=str)
        verts = {}
        for (slot, c) in vkeys:
            v = piN.vertex_graph(slot, c, f'vtx_{slot}_{c}')
            if c == 'mag':
                v.iso = {'kp': piN.isospin(v, piN.Km('p')), 'kn': piN.isospin(v, piN.Km('n'))}
            else:
                v.iso = {'1': piN.isospin(v, None)}
            verts[(slot, c)] = v
        sg = piN.self_graph('self_piN'); sg.iso = {'1': piN.isospin(sg, None)}
        allg = comp + subs + list(verts.values()) + [sg] + list(extra_graphs)
        self.cmp = ev.compile_graphs(allg, nproc)
        self.vc = {k: ev.vertex_consts(self.cmp[v.name], cs) for k, v in verts.items()}
        self.L0, self.L1 = ev.self_consts(self.cmp[sg.name], cs)
        self.vcP = {k: ev.vertex_consts(self.cmp[v.name], cs, pole=True) for k, v in verts.items()}
        self.L0P, self.L1P = ev.self_consts(self.cmp[sg.name], cs, pole=True)
        self.subs = {h.name: h for h in subs}

    def sub_values(self, h, pole=False):
        """symbol values and coupling factor for a subtraction graph (overall minus sign included)"""
        cs = self.cs
        kp = cs['kp']
        if hasattr(h, 'vkey'):
            c1, c2 = (self.vcP if pole else self.vc)[h.vkey]
            fO = 1.0 if h.Ocpl == 'ch' else kp
            return dict(F1s=c1, F2s=c2), -fO
        c1, c2 = h.cpls
        f = (1.0 if c1 == 'ch' else kp) * (1.0 if c2 == 'ch' else kp)
        return (dict(L0s=self.L0P, L1s=self.L1P) if pole else dict(L0s=self.L0, L1s=self.L1)), -f

    def A(self, x, t, parts=False, pole=False):
        import kin
        cs = self.cs
        en = kin.env(x, t, cs['mpi2']); en['M'] = 1.0
        V = np.zeros(6, complex); byname = {}
        for name, c in self.cmp.items():
            if c.g.mode != 'compton': continue
            e = en
            f = c.coup(cs) * c.isof(cs)
            if name in self.subs:
                vals, s = self.sub_values(self.subs[name], pole); e = dict(en); e.update(vals); f *= s
            v = c.raw(e, pole) * f
            V += v
            if parts: byname[name] = v
        G = ev.gram(en)
        A = np.linalg.solve(G, V)
        return (A, V, G, byname) if parts else A

    def gi_compile(self, nproc=8):
        import runform, formtools as ft
        gs = [c.g for c in self.cmp.values() if c.g.mode == 'compton']
        res = runform.project_all(gs, nproc, gi=True)
        self.gic = {n: [compile(ft.to_python(d.get(f'V{j}', '0')), f'W{j}', 'eval') for j in range(1, 9)]
                    for n, d in res.items()}

    def GI(self, x, t):
        """projections of the amplitude with ee -> q (must vanish); returns total and the largest single graph"""
        import kin
        cs = self.cs
        en = kin.env(x, t, cs['mpi2']); en['M'] = 1.0
        W = np.zeros(8, complex); big = 0.
        for name, code in self.gic.items():
            c = self.cmp[name]
            T, U = c.tensors(en)
            e = dict(en); e['T'] = T; e['U'] = U
            f = c.coup(cs) * c.isof(cs)
            if name in self.subs:
                vals, sgn = self.sub_values(self.subs[name]); e.update(vals); f *= sgn
            w = np.array([eval(cc, e) for cc in code]) * f
            W += w; big = max(big, np.max(np.abs(w)))
        return W, big
