"""graphs.py -- generator of the one-loop (and tree) Compton graphs of covariant BChPT for FORM.

Photon A: incoming, momentum q, polarisation ee.  Photon B: outgoing, momentum qp, polarisation ep.
Nucleon line read from the incoming nucleon p (right) to the outgoing nucleon pp (left).
Feynman rules (pseudovector form of L^(1), photon couplings via the covariant derivatives; e removed,
Gam = T/e^2 with iT = sum of graphs):
  N propagator  i(p/+M)/(p^2-M^2),  pi propagator  i delta^{ab}/(k^2-m^2)
  gamma N N     charge: i Q eps/ ;  magnetic (kappa): -i/(2M) kappa ag(eps, k_in)     [ag(a,b) = a/b/ - a.b]
  pi N N        -(gA/2F) k/ g5 tau^a  (k incoming pion momentum)
  KR gamma pi N N   -i (gA/2F) eps^{3ab} eps/ g5 tau^b   (pion a)
  gamma pi pi   -eps^{3cd} (k_in + k_out).eps     (pion c in, d out along the line)
  gamma gamma pi pi   2i (delta^{cd} - delta^{c3} delta^{d3}) ee.ep
Loop measure: int d^4l/(2pi)^4 = i/(16 pi^2) x LoopTools.
"""
from fractions import Fraction as Fr
import itertools, numpy as np

BASIS = ('l', 'p', 'q', 'qp')

class Mom(dict):
    """linear combination of l, p, q, qp (pp = p + q - qp)"""
    def __init__(self, **kw):
        super().__init__({k: Fr(v) for k, v in kw.items() if v != 0})
    def __add__(self, o):
        r = Mom(**self)
        for k, v in o.items():
            r[k] = r.get(k, Fr(0)) + v
            if r[k] == 0: del r[k]
        return r
    def __neg__(self): return Mom(**{k: -v for k, v in self.items()})
    def __sub__(self, o): return self + (-o)
    def __rmul__(self, c): return Mom(**{k: Fr(c) * v for k, v in self.items()})
    def subs_l(self, a, R):
        """l -> a*l + R"""
        c = self.get('l', Fr(0))
        r = Mom(**{k: v for k, v in self.items() if k != 'l'})
        if c != 0:
            r = r + Mom(l=c * a) + c * R
        return r
    def lcoef(self): return self.get('l', Fr(0))
    def ext(self): return Mom(**{k: v for k, v in self.items() if k != 'l'})
    def form(self, wrap):
        """FORM text of the linear combination; wrap(vec) gives the text of one vector"""
        if not self: return '0'
        out = ''
        for k in BASIS:
            if k in self:
                c = self[k]
                sgn = '+' if c > 0 else '-'
                a = abs(c)
                out += sgn + ('' if a == 1 else f'{a}*') + wrap(k)
        return '(' + out + ')'

P_ = Mom(p=1); Q_ = Mom(q=1); QP_ = Mom(qp=1); L_ = Mom(l=1); PP_ = Mom(p=1, q=1, qp=-1)

import sympy as _sp
_x, _t, _dl, _ss = _sp.symbols('x t dl ss')
KTAB = {'compton': {('p', 'p'): 1, ('q', 'q'): 0, ('qp', 'qp'): 0, ('p', 'q'): _x, ('p', 'qp'): _x + _t / 2,
                    ('q', 'qp'): -_t / 2},
        'vertex': {('p', 'p'): 1, ('q', 'q'): -_dl, ('p', 'q'): _dl / 2},
        'self': {('p', 'p'): _ss}}
def dot_sym(a, b, mode='compton'):
    """scalar product of two external momenta (no l), M=1 -> sympy expression"""
    import sympy as sp
    tab = KTAB[mode]
    r = 0
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            key = (k1, k2) if (k1, k2) in tab else (k2, k1)
            r += sp.Rational(v1.numerator, v1.denominator) * sp.Rational(v2.numerator, v2.denominator) * tab[key]
    return sp.expand(r)

# ---------------------------------------------------------------- Dirac-string tokens
class Tok:
    """one factor of the Dirac string.  kind: 'V' slashed vector (polarisation or external),
    'K' slashed momentum, 'S' nucleon numerator (mom/ + M), 'AG' ag(vec, mom), 'I' identity;
    g5: the factor is followed by gamma5"""
    def __init__(self, kind, a=None, b=None, g5=False, sign=1):
        self.kind, self.a, self.b, self.g5, self.sign = kind, a, b, g5, sign
    def subs_l(self, a, R):
        f = lambda m: m.subs_l(a, R) if isinstance(m, Mom) else m
        return Tok(self.kind, f(self.a), f(self.b), self.g5, self.sign)
    def flipped(self):
        """gamma5 X gamma5"""
        if self.kind in ('V', 'K'): return Tok(self.kind, self.a, self.b, self.g5, -self.sign)
        if self.kind == 'S': return Tok('Sf', self.a, self.b, self.g5, self.sign)
        if self.kind == 'Sf': return Tok('S', self.a, self.b, self.g5, self.sign)
        return Tok(self.kind, self.a, self.b, self.g5, self.sign)
    def form(self):
        s = '' if self.sign == 1 else '-'
        if self.kind == 'V': return f'{s}g_(1,{self.a})'
        if self.kind == 'K': return s + self.a.form(lambda v: f'g_(1,{v})')
        if self.kind == 'S': return s + '(' + self.a.form(lambda v: f'g_(1,{v})') + '+M*gi_(1))'
        if self.kind == 'Sf': return s + '(-' + self.a.form(lambda v: f'g_(1,{v})') + '+M*gi_(1))'
        if self.kind == 'AG':
            return s + '(' + self.b.form(lambda v: f'g_(1,{self.a},{v})') + '-' + self.b.form(lambda v: f'{self.a}.{v}') + ')'
        if self.kind == 'I': return s + 'gi_(1)'
        if self.kind == 'VEFF':
            return s + '(F1s*g_(1,' + self.a + ')+F2s*(' + self.b.form(lambda v: f'g_(1,{self.a},{v})') + '-' + \
                self.b.form(lambda v: f'{self.a}.{v}') + '))'
        if self.kind == 'LINS':
            return s + '(L0s*gi_(1)+L1s*(' + self.a.form(lambda v: f'g_(1,{v})') + '-M*gi_(1)))'
        raise ValueError(self.kind)

def remove_g5(toks):
    """(A g5) Y (B g5) = - A Ytilde B"""
    pos = [i for i, tk in enumerate(toks) if tk.g5]
    if not pos: return 1, toks
    assert len(pos) == 2, 'odd number of gamma5'
    i, j = pos
    new = []
    for k, tk in enumerate(toks):
        tk2 = Tok(tk.kind, tk.a, tk.b, False, tk.sign)
        if i < k < j: tk2 = tk2.flipped()
        new.append(tk2)
    return -1, new

# ---------------------------------------------------------------- graph container
class Graph:
    def __init__(self, name):
        self.name = name
        self.toks = []          # left -> right
        self.scal = []          # scalar factors: list of (Mom, vector-name) dot products or FORM text
        self.loop = []          # [(Mom, mass2 symbol)]
        self.tree = []          # tree denominators: 'ixs' / 'iys'
        self.pref = 1 + 0j      # numeric prefactor (i's, gA/2F etc. with gA=F=1; see self.coup)
        self.coup = {}          # symbolic coupling powers {'gA':2,'F':-2,...}
        self.iso = None         # isospin factor: dict {coupling-monomial: number}
        self.tags = {}
        self.mode = 'compton'    # 'compton' | 'vertex' (on-shell gamma N N) | 'self' (self-energy)
    def mul_coup(self, **kw):
        for k, v in kw.items(): self.coup[k] = self.coup.get(k, 0) + v

    def shifted(self):
        """move the loop momentum so that the first loop propagator is l; return (toks, scal, ks)"""
        if not self.loop:
            return self.toks, self.scal, []
        m0 = self.loop[0][0]
        s0 = m0.lcoef(); assert abs(s0) == 1
        a, R = Fr(1) / s0, -(Fr(1) / s0) * m0.ext()      # l = s0*(L - K0)
        toks = [tk.subs_l(a, R) for tk in self.toks]
        scal = [(it[0].subs_l(a, R), it[1]) if isinstance(it, tuple) else it for it in self.scal]
        ks = []
        for mom, m2 in self.loop:
            mm = mom.subs_l(a, R)
            c = mm.lcoef(); assert abs(c) == 1, mom
            ks.append(c * mm.ext())                       # (c L + R)^2 = (L + c R)^2
        assert not ks[0], ks[0]
        return toks, scal, ks[1:]

    def form_body(self, gi=False):
        toks, scal, ks = self.shifted()
        sgn, toks = remove_g5(toks)
        string = '*'.join(tk.form() for tk in toks) if toks else 'gi_(1)'
        sc = []
        for item in scal:
            if isinstance(item, tuple):
                m, v = item
                sc.append(m.form(lambda w: f'{w}.{v}'))
            else:
                sc.append(item)
        tree = '*'.join(self.tree)
        pieces = [x for x in [tree, '*'.join(sc), string] if x]
        expr = ('-' if sgn < 0 else '') + '*'.join(pieces)
        ksub = '\n'.join(f'  id k{i+1} = {k.form(lambda v: v)};' for i, k in enumerate(ks))
        ksub += '\n' + '\n'.join(f'  id k{i+1} = 0;' for i in range(len(ks), 3))
        pf = {'compton': '#call Proj\n#call Finish\n', 'vertex': '#call VProj\n#call VFinish\n',
              'self': '#call SFinish\n'}[self.mode]
        if gi:
            import re
            expr = re.sub(r'\bee\b', 'q', expr)
            pf = '#call GProj\n#call GFinish\n'
        return f'#procedure KSUB\n{ksub}\n#endprocedure\nL GAM = {expr};\n.sort\n{pf}', ks

    def lt_args(self):
        """LoopTools call (kind, [invariants as sympy], [masses]) for the loop"""
        import sympy as sp
        _, _, ks = self.shifted()
        kk = [Mom()] + ks
        ms = [m2 for _, m2 in self.loop]
        n = len(self.loop)
        sq = lambda m: dot_sym(m, m, self.mode)
        if n == 1: return ('A', [], ms)
        if n == 2: return ('B', [sq(kk[1])], ms)
        if n == 3: return ('C', [sq(kk[1]), sq(kk[2] - kk[1]), sq(kk[2])], ms)
        if n == 4: return ('D', [sq(kk[1]), sq(kk[2] - kk[1]), sq(kk[3] - kk[2]), sq(kk[3]), sq(kk[2]),
                                 sq(kk[3] - kk[1])], ms)
        raise ValueError(n)

# ---------------------------------------------------------------- isospin
tau = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
def eps3(a, b):
    """epsilon^{3ab}, a,b in 0..2"""
    if (a, b) == (0, 1): return 1
    if (a, b) == (1, 0): return -1
    return 0
