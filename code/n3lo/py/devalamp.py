"""devalamp.py -- FORM projections and numerical evaluation of graphs with Delta lines (topology-tagged tensors)."""
import os, re, sys, concurrent.futures as cf, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'lt'))
import ltpy, formtools as ft, kin, evalamp as ev, bhigh

MD_MEV = 1232.0
def dconsts(c=ev.LMP15, hA=2.85, gM=2.97, gE=-1.0, MD=MD_MEV):
    cs = ev.consts(c)
    md = MD / c['M']
    cs.update(hA=hA, gM=gM, gE=gE, MD=md, iMp=1 / (1 + md))
    return cs

def project(g, gi=False):
    body, _ = g.form_body(gi=gi)
    d = ft.run_form(g.name + ('_GI' if gi else ''), body, header='dcommon.h', bracket='TCFT,UCFT')
    stuck = [k for k, v in d.items() if 'STUCK' in v]
    if stuck:
        STUCK.append(g.name)
        if not ALLOW_STUCK: raise RuntimeError(f'{g.name}: rank reduction stuck in {stuck}')
    return g.name, d

STUCK, ALLOW_STUCK = [], True     # stuck terms: rank-4 bubbles, evaluated with bhigh.py

def project_all(graphs, nproc=8, gi=False):
    with cf.ThreadPoolExecutor(nproc) as ex:
        return dict(ex.map(lambda g: project(g, gi), graphs))

def _A(m2):
    fb, ub = ltpy.B(0., 0., m2)
    a0, u0 = m2 * fb['bb0'], m2 * ub['bb0']
    fin = {'aa0': a0, 'aa00': m2 * a0 / 4 + m2 ** 2 / 8, 'aa0000': m2 ** 2 * a0 / 24 + 5 * m2 ** 3 / 144}
    uv = {'aa0': u0, 'aa00': m2 * u0 / 4, 'aa0000': m2 ** 2 * u0 / 24}
    return fin, uv

class DCompiled:
    def __init__(self, g, d):
        self.g = g
        nproj = {'compton': 6, 'vertex': 3, 'self': 2}[g.mode] if not getattr(g, 'gi', False) else 8
        self.src = [ft.to_python(d.get(f'V{j}', '0')) for j in range(1, nproj + 1)]
        self.proj = [compile(s, f'V{j}', 'eval') for j, s in enumerate(self.src)]
        masks = set()
        for s in self.src:
            masks |= {int(m) for m in re.findall(r'[TU]\[\((\d+),\d+,\d+,\d+,\d+\)\]', s)}
        self.masks = sorted(masks)
        syms = sp.symbols('x t dl ss')
        self.args = {}
        for m in self.masks:
            kind, inv, ms = g.topo_args(m)
            self.args[m] = (kind, [sp.lambdify(syms, e, 'math') for e in inv], ms)
    def tensors(self, en):
        T, U = {}, {}
        for m in self.masks:
            kind, inv, ms = self.args[m]
            a = [f(en.get('x', 0), en.get('t', 0), en.get('dl', 0), en.get('ss', 0)) for f in inv]
            mm = [en['mpi2'] if s == 'mpi2' else (en['MD'] ** 2 if s == 'MD2' else 1.0) for s in ms]
            if kind == 'A': fin, uv = _A(mm[0])
            elif kind == 'B':
                fin, uv = ltpy.B(*a, *mm)
                for key in _keys(self.src, m):                 # ranks > 3: Feynman-parameter form (bhigh.py)
                    nm = ev.coef_name(kind, key)
                    if nm not in fin:
                        fin[nm], uv[nm] = bhigh.bcoef(key[0] // 2, key[1], a[0], mm[0], mm[1])
            elif kind == 'C': fin, uv = ltpy.C(*a, *mm)
            else: fin, uv = ltpy.D(*a, *mm)
            for key in _keys(self.src, m):
                nm = ev.coef_name(kind, key)
                T[(m,) + key] = fin[nm]; U[(m,) + key] = uv[nm]
        return T, U
    def raw(self, en, pole=False):
        T, U = self.tensors(en)
        if pole:
            T, U = U, {k: 0.0 for k in U}
        T.setdefault((0, 0, 0, 0), 1.0); U.setdefault((0, 0, 0, 0), 0.0)     # tree graphs: TCF(0,0,0,0) = 1
        e = dict(en); e['T'] = T; e['U'] = U; e['STUCK'] = 1.0   # STUCK only flags terms left at rank 4 (bubbles)
        return np.array([eval(c, e) for c in self.proj])
    def coup(self, cs):
        f = self.g.pref
        for k, v in self.g.coup.items(): f *= cs[k] ** v
        return f
    def isof(self, cs):
        return sum(v * (1.0 if k == '1' else cs[k]) for k, v in self.g.iso.items())

def _keys(srcs, m):
    ks = set()
    for s in srcs:
        for a in re.findall(r'[TU]\[\(' + str(m) + r',(\d+),(\d+),(\d+),(\d+)\)\]', s):
            ks.add(tuple(int(x) for x in a))
    return ks

def compile_graphs(graphs, nproc=8, gi=False):
    res = project_all(graphs, nproc, gi)
    out = {}
    for g in graphs:
        if gi: g.gi = True
        out[g.name] = DCompiled(g, res[g.name])
        if gi: g.gi = False
    return out

def env(x, t, cs, **kw):
    en = kin.env(x, t, cs['mpi2']); en['M'] = 1.0; en['MD'] = cs['MD']
    s, u = 1 + 2 * x, 1 + 2 * en['y']
    en.update(ixD=1 / (s - cs['MD'] ** 2), iyD=1 / (u - cs['MD'] ** 2), iMD=1 / (1 - cs['MD'] ** 2))
    en.update(kw)
    return en

def amplitude(comp, x, t, cs, extra=None, pole=False, names=None):
    en = env(x, t, cs)
    if extra: en.update(extra)
    V = np.zeros(6, complex)
    for nm, c in comp.items():
        if names is not None and nm not in names: continue
        V += c.raw(en, pole) * c.coup(cs) * c.isof(cs)
    return np.linalg.solve(ev.gram(en), V)
