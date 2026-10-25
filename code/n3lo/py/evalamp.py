"""evalamp.py -- numerical evaluation of the FORM-projected graphs (Python + LoopTools-dq), for checks.

Amplitudes A_1..A_6 in the basis of hel_from_A (ampl.f90), nucleon-mass units, Gam = T/e^2.
"""
import os, sys, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'lt'))
import ltpy, formtools as ft, kin
from graphs import Tok, Graph, Mom, P_, Q_, QP_

LMP15 = dict(M=938.27, mpi=139.57, F=92.21, gA=1.27, kp=1.793, kn=-1.913043)

def consts(c=LMP15):
    M = c['M']
    return dict(gA=c['gA'], F=c['F'] / M, mpi2=(c['mpi'] / M) ** 2, kp=c['kp'], kn=c['kn'],
                c1=c.get('c1', 0.) * M * 1e-3, c2=c.get('c2', 0.) * M * 1e-3,
                c3=c.get('c3', 0.) * M * 1e-3, c4=c.get('c4', 0.) * M * 1e-3)   # c_i in GeV^-1 -> M units

def coef_name(kind, key):
    n2, a, b, c = key
    pre = {'A': 'aa', 'B': 'bb', 'C': 'cc', 'D': 'dd'}[kind]
    idx = '00' * (n2 // 2) + '1' * a + '2' * b + '3' * c
    return pre + (idx if idx else '0')

class Compiled:
    def __init__(self, g, d):
        self.g = g
        self.proj = ft.compile_proj(d) if g.mode == 'compton' else \
            [compile(ft.to_python(d.get(f'V{j}', '0')), f'V{j}', 'eval') for j in range(1, 4 if g.mode == 'vertex' else 3)]
        self.kind, inv, self.ms = (g.lt_args() if g.loop else ('T', [], []))
        syms = sp.symbols('x t dl ss')
        self.inv = [sp.lambdify(syms, e, 'math') for e in inv]
    def tensors(self, en):
        if self.kind == 'T': return {(0, 0, 0, 0): 1.0}, {}
        args = [f(en.get('x', 0), en.get('t', 0), en.get('dl', 0), en.get('ss', 0)) for f in self.inv]
        ms = [en['mpi2'] if m == 'mpi2' else 1.0 for m in self.ms]
        if self.kind == 'A':
            fb, ub = ltpy.B(0., 0., ms[0])
            fin = {'aa0': ms[0] * fb['bb0']}; uv = {'aa0': ms[0] * ub['bb0']}
        elif self.kind == 'B': fin, uv = ltpy.B(*args, *ms)
        elif self.kind == 'C': fin, uv = ltpy.C(*args, *ms)
        elif self.kind == 'D': fin, uv = ltpy.D(*args, *ms)
        else: raise ValueError(self.kind)
        class LTdict(dict):
            def __init__(s, src, kind): s.src, s.kind = src, kind
            def __missing__(s, key): return s.src[coef_name(s.kind, key)]
        return LTdict(fin, self.kind), LTdict(uv, self.kind)
    def raw(self, en, pole=False):
        T, U = self.tensors(en)
        if pole and self.kind != 'T':
            class Z(dict):
                def __missing__(s, k): return 0.0
            T, U = U, Z()
        e = dict(en); e['T'] = T; e['U'] = U
        return np.array([eval(c, e) for c in self.proj])
    def coup(self, cs):
        f = self.g.pref
        for k, v in self.g.coup.items(): f *= cs[k] ** v
        return f
    def isof(self, cs):
        return sum(v * (1.0 if k == '1' else cs[k]) for k, v in self.g.iso.items())

def compile_graphs(graphs, nproc=8):
    import runform
    res = runform.project_all(graphs, nproc)
    return {g.name: Compiled(g, res[g.name]) for g in graphs}

_gram = None
def gram(en):
    global _gram
    if _gram is None:
        S = ['(-ep.ee)*gi_(1)', '(q.ep*qp.ee)*gi_(1)', '(-(g_(1,ep,ee)-ep.ee))', '(ep.ee*(g_(1,qp,q)-qp.q))',
             '(q.ep*(g_(1,qp,ee)-qp.ee)-qp.ee*(g_(1,q,ep)-q.ep))', '(q.ep*(g_(1,q,ee)-q.ee)-qp.ee*(g_(1,qp,ep)-qp.ep))']
        KS = '#procedure KSUB\n#endprocedure\n'
        _gram = [ft.compile_proj(ft.run_form(f'gram{k+1}', KS + f'L GAM = {S[k]};\n.sort\n#call Proj\n#call Finish\n'))
                 for k in range(6)]
    e = dict(en); e['T'] = {(0, 0, 0, 0): 1.0}; e['U'] = {}
    return np.array([[eval(_gram[k][j], e) for k in range(6)] for j in range(6)])

_vgram = None
def vgram(en):
    global _vgram
    if _vgram is None:
        R = ['g_(1,ee)', 'P.ee*gi_(1)', 'q.ee*gi_(1)']
        KS = '#procedure KSUB\n#endprocedure\n'
        _vgram = []
        for k in range(3):
            d = ft.run_form(f'vgram{k+1}', KS + f'L GAM = {R[k]};\n.sort\n#call VProj\n#call VFinish\n')
            _vgram.append([compile(ft.to_python(d[f'V{j}']), 'v', 'eval') for j in range(1, 4)])
    e = dict(en); e['T'] = {(0, 0, 0, 0): 1.0}; e['U'] = {}
    return np.array([[eval(_vgram[k][j], e) for k in range(3)] for j in range(3)])

def vertex_consts(cv, cs, dls=(1e-4, 2e-4, 3e-4), pole=False):
    """on-shell (c1, c2) of a compiled vertex-loop graph: Lambda ~ c1 ee/ + c2 ag(ee,q), q^2 -> 0 (quadratic
    extrapolation from dl = h, 2h, 3h, error O(h^3); the linear one from 2h, 4h left a residual nucleon pole ~1e-7)"""
    out = []
    for dl in dls:
        en = dict(mpi2=cs['mpi2'], dl=dl, M=1.0)
        W = cv.raw(en, pole) * cv.coup(cs) * cv.isof(cs)
        A, B, C = np.linalg.solve(vgram(en), W)
        out.append(np.array([A + 2 * B, B]))
    return 3 * out[0] - 3 * out[1] + out[2]

def self_consts(cv, cs, h=1e-4, pole=False):
    """L(M), L'(M) of a compiled self-energy graph, L(p) = a(p^2) p/ + b(p^2)"""
    def ab(ss):
        en = dict(mpi2=cs['mpi2'], ss=ss)
        Y = cv.raw(en, pole) * cv.coup(cs) * cv.isof(cs)
        return Y[1] / ss, Y[0]
    a0, b0 = ab(1.0)
    (a1, b1), (a2, b2), (am1, bm1), (am2, bm2) = ab(1 + h), ab(1 + 2 * h), ab(1 - h), ab(1 - 2 * h)
    da = (8 * (a1 - am1) - (a2 - am2)) / (12 * h); db = (8 * (b1 - bm1) - (b2 - bm2)) / (12 * h)   # O(h^4)
    return a0 + b0, a0 + 2 * (da + db)      # M = 1

def amplitude(comp, x, t, cs, extra=None):
    """sum of the compiled Compton graphs -> A_1..A_6"""
    en = kin.env(x, t, cs['mpi2']); en['M'] = 1.0
    if extra: en.update(extra)
    V = np.zeros(6, complex)
    for c in comp.values():
        if c.g.mode != 'compton': continue
        V += c.raw(en) * c.coup(cs) * c.isof(cs)
    return np.linalg.solve(gram(en), V), V
