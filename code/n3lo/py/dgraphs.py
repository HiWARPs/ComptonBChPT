"""dgraphs.py -- Dirac-string tokens with Delta(1232) lines and the FORM program with loop-momentum reduction.

The RS propagator and the consistent couplings raise the tensor rank of the loop integrals beyond what LoopTools
provides (D: 5, C: 4, B: 3).  After the traces and polarisation sums the numerator depends on l only through
l.l, l.p, l.q, l.qp; terms whose rank exceeds the limit of their topology are reduced one scalar product at a time,
  l.l = DD(b) + m_b^2,   l.R_i = (DD(i) - DD(b) - R_i^2 + m_i^2 - m_b^2)/2,
(DD = inverse propagator, b = first propagator of the topology, R_i = K_i - K_b), and DD(i) cancels propagator i
(TOPO bitmask; the loop momentum is re-routed when the first propagator is cancelled).  Exact algebra, no Gram
determinants.  The remaining tensor integrals of each sub-topology are reduced as before (TCFT(mask, ...)).
"""
import itertools
from fractions import Fraction as Fr
from graphs import Mom, P_, Q_, QP_, L_, dot_sym

AUX = None          # dict Mom-key -> auxiliary vector name while a FORM body is generated (see DGraph.form_body)
AUX_ON = False      # auxiliary vectors for multi-term momenta: identical results, ~1.5x faster traces; needs the
                    # Vectors kv1..kv24 declaration and '#call AUXSUB' after tracen in dcommon.h (not active)

def _aux(m):
    """auxiliary vector for a multi-term momentum (fewer terms before the trace); None for +-single vectors"""
    if AUX is None or not isinstance(m, Mom) or not m: return None
    if len(m) == 1 and abs(list(m.values())[0]) == 1: return None
    key = tuple(sorted(m.items()))
    if key not in AUX: AUX[key] = (f'kv{len(AUX) + 1}', m)
    return AUX[key][0]

def sl(x):
    a = _aux(x)
    if a: return f'g_(1,{a})'
    if isinstance(x, Mom): return x.form(lambda v: f'g_(1,{v})')
    return f'g_(1,{x})'
def vidx(x, i):
    a = _aux(x)
    if a: return f'{a}({i})'
    if isinstance(x, Mom): return x.form(lambda v: f'{v}({i})')
    return f'{x}({i})'
def dot(x, y):
    ax, ay = _aux(x), _aux(y)
    if ax or ay:
        tx = ax if ax else None; ty = ay if ay else None
        if tx and ty: return f'{tx}.{ty}'
        if tx: return dot(Mom(**{tx: 1}), y) if isinstance(y, Mom) else f'{tx}.{y}'
        return dot(x, Mom(**{ty: 1})) if isinstance(x, Mom) else f'{x}.{ty}'
    if isinstance(x, Mom) and isinstance(y, Mom):
        t = [f'{ca*cb}*{a}.{b}' for a, ca in x.items() for b, cb in y.items()]
        return '(' + ('+'.join(t) if t else '0').replace('+-', '-') + ')'
    if isinstance(x, Mom): return x.form(lambda v: f'{v}.{y}')
    if isinstance(y, Mom): return y.form(lambda v: f'{x}.{v}')
    return f'{x}.{y}'

def _perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s

class DTok:
    """kind, args (Mom / vector name / index name), g5 (gamma5 at the right end of the factor)"""
    def __init__(self, kind, *args, g5=False):
        self.kind, self.args, self.g5 = kind, args, g5
    def subs_l(self, a, R):
        return DTok(self.kind, *[x.subs_l(a, R) if isinstance(x, Mom) else x for x in self.args], g5=self.g5)
    def nog5(self): return DTok(self.kind, *self.args, g5=False)
    def text(self, flip=False):
        k, a = self.kind, self.args
        sg = -1 if flip else 1
        if k == 'S':                       # p/ + M
            return f'({"-" if flip else ""}{sl(a[0])}+M*gi_(1))'
        if k in ('V', 'K'):
            return ('-' if flip else '') + sl(a[0])
        if k == 'AG':                      # a/ b/ - a.b (even)
            return f'({sl(a[0])}*{sl(a[1])}-{dot(a[0], a[1])})'
        if k in ('RS', 'RSR', 'RSH', 'RS4', 'RSR4', 'RSH4'):
            # D-dim RS numerator, -(p/+MD)[g - gg/(D-1) - (D-2)pp/((D-1)MD^2) + (p g - g p)/((D-1)MD)];
            # RSR: both neighbours transverse (p^i1, p^i2 terms drop); RSH: one side ('L'/'R');
            # suffix 4: tree propagator outside the loop, 4-dim (1/(D-1) -> 1/3; external line of the loop)
            if k.endswith('4'):
                return self.__class__(k[:-1], *a, g5=self.g5).text(flip).replace('iD1', '(1/3)')
            P, i1, i2 = a[:3]
            s1 = '-' if flip else ''
            s2 = '-' if flip else '+'
            br = f'd_({i1},{i2})-g_(1,{i1})*g_(1,{i2})*iD1'
            if k == 'RS':
                br += f'-(1-iD1)*{vidx(P, i1)}*{vidx(P, i2)}*MD^-2{s2}({vidx(P, i1)}*g_(1,{i2})-g_(1,{i1})*{vidx(P, i2)})*iD1*MD^-1'
            elif k == 'RSH':
                br += f'{s2}(-g_(1,{i1})*{vidx(P, i2)})*iD1*MD^-1' if a[3] == 'L' else f'{s2}({vidx(P, i1)}*g_(1,{i2}))*iD1*MD^-1'
            return f'(-({s1}{sl(P)}+MD*gi_(1))*({br}))'
        if k == 'G3':                      # gamma^{mu i lam} a_mu b_lam = (a/ g^i b/ - b/ g^i a/)/2  (odd)
            x, i, y = a
            return f'({"-" if flip else ""}({sl(x)}*g_(1,{i})*{sl(y)}-{sl(y)}*g_(1,{i})*{sl(x)})/2)'
        if k == 'G3E':                     # gamma^{i1 i2 al} e_al  (odd)
            i1, i2, e = a
            return f'({"-" if flip else ""}(g_(1,{i1})*g_(1,{i2})*g_(1,{e})-g_(1,{e})*g_(1,{i2})*g_(1,{i1}))/2)'
        if k == 'EPS4':                    # gamma^{[abcd]} (antisymmetrised, even), D-dim identity:
            A, B, C, Dd = a                # = abcd - (ab)cd + (ac)bd - (ad)bc - (bc)ad + (bd)ac - (cd)ab + (ab)(cd) - (ac)(bd) + (ad)(bc)
            if not EPS4_SHORT:
                terms = []
                for p in itertools.permutations(range(4)):
                    terms.append(('+' if _perm_sign(p) > 0 else '-') + '*'.join(sl(a[j]) for j in p))
                return '((' + ''.join(terms) + ')/24)'
            G = lambda *x: '*'.join(sl(y) for y in x)
            m = _metric
            return (f'({G(A, B, C, Dd)}-{m(A, B)}*{G(C, Dd)}+{m(A, C)}*{G(B, Dd)}-{m(A, Dd)}*{G(B, C)}'
                    f'-{m(B, C)}*{G(A, Dd)}+{m(B, Dd)}*{G(A, C)}-{m(C, Dd)}*{G(A, B)}'
                    f'+{m(A, B)}*{m(C, Dd)}-{m(A, C)}*{m(B, Dd)}+{m(A, Dd)}*{m(B, C)})')
        if k == 'SCI':                     # (P.K) e^i - (P.e) K^i   (c-number)
            P, K, e, i = a
            return f'({dot(P, K)}*{vidx(e, i)}-{dot(P, e)}*{vidx(K, i)})'
        if k == 'I': return 'gi_(1)'
        raise ValueError(k)

EPS4_SHORT = True

def _is_index(x): return isinstance(x, str) and x.startswith('a') and x[1:].isdigit()

def _metric(x, y):
    """g(x, y) for Mom / polarisation / Delta-index entries"""
    if _is_index(x) and _is_index(y): return f'd_({x},{y})'
    if _is_index(x): return vidx(y, x)
    if _is_index(y): return vidx(x, y)
    return dot(x, y)

def remove_g5(toks):
    """(A g5) Y (B g5) = A Ytilde Btilde, Xtilde = g5 X g5 (anticommuting gamma5)"""
    pos = [i for i, tk in enumerate(toks) if tk.g5]
    if not pos: return [(tk, False) for tk in toks]
    assert len(pos) == 2, 'odd number of gamma5'
    i, j = pos
    return [(tk.nog5(), i < k <= j) for k, tk in enumerate(toks)]

LIM = {4: 5, 3: 4, 2: 3, 1: 4}
RS4D = False          # True: 4-dim RS propagator (1/(D-1) -> 1/3); tested: not the scheme of the NNLO code
MODEVEC = {'compton': ['p', 'q', 'qp'], 'vertex': ['p', 'q'], 'self': ['p']}

def _vec(m, basis):
    return [Fr(m.get(b, 0)) for b in basis]

def _solve(cols, rhs):
    """least-squares-free exact solve: cols (list of vectors) spanning rhs; returns coefficients or None"""
    import sympy as sp
    A = sp.Matrix([[sp.Rational(c.numerator, c.denominator) for c in col] for col in cols]).T
    b = sp.Matrix([sp.Rational(r.numerator, r.denominator) for r in rhs])
    sol = A.gauss_jordan_solve(b)[0]
    return [Fr(int(sp.fraction(s)[0]), int(sp.fraction(s)[1])) for s in sol]

def masstext(m):
    return {'mpi2': 'mpi2', '1': 'M^2', 'MD2': 'MD^2'}[m]

class DGraph:
    def __init__(self, name):
        self.name = name
        self.toks = []; self.scal = []; self.loop = []; self.tree = []
        self.pref = 1 + 0j; self.coup = {}; self.iso = None; self.mode = 'compton'; self.tags = {}
    def mul_coup(self, **kw):
        for k, v in kw.items(): self.coup[k] = self.coup.get(k, 0) + v
    def shifted(self):
        if not self.loop: return self.toks, self.scal, []
        m0 = self.loop[0][0]
        s0 = m0.lcoef(); assert abs(s0) == 1
        a, R = Fr(1) / s0, -(Fr(1) / s0) * m0.ext()
        toks = [tk.subs_l(a, R) for tk in self.toks]
        scal = [(it[0].subs_l(a, R), it[1]) if isinstance(it, tuple) else it for it in self.scal]
        ks = []
        for mom, m2 in self.loop:
            mm = mom.subs_l(a, R); c = mm.lcoef(); assert abs(c) == 1, mom
            ks.append(c * mm.ext())
        assert not ks[0]
        return toks, scal, ks[1:]
    def Ks(self):
        _, _, ks = self.shifted()
        return [Mom()] + ks
    # ------------------------------------------------------------------ FORM
    def reduce_proc(self):
        """FORM procedure REDUCE (rank reduction) and KSUBT (LoopTools momenta per topology)"""
        n = len(self.loop)
        if n == 0:
            return '#procedure REDUCE\n#endprocedure\n#procedure KSUBT\n#endprocedure\n'
        K = self.Ks(); ms = [masstext(m) for _, m in self.loop]
        basis = MODEVEC[self.mode]
        red, canc, shifts, ksub = [], [], {}, []
        for mask in range(1, 2 ** n):
            mem = [i for i in range(n) if mask >> i & 1]
            b = mem[0]
            R = {i: K[i] - K[b] for i in mem}
            # independent R's and complement directions
            ind = []
            for i in mem[1:]:
                vecs = [_vec(R[j], basis) for j in ind] + [_vec(R[i], basis)]
                if _rank(vecs) > len(ind): ind.append(i)
            comp = []
            for v in basis:
                vecs = [_vec(R[j], basis) for j in ind] + [_vec(Mom(**{w: 1}), basis) for w in comp] + \
                       [_vec(Mom(**{v: 1}), basis)]
                if _rank(vecs) > len(ind) + len(comp): comp.append(v)
            rules = [f'      id,once l.l = DD({b})+{ms[b]};']
            for v in basis:
                if v in comp: continue
                cols = [_vec(R[j], basis) for j in ind] + [_vec(Mom(**{w: 1}), basis) for w in comp]
                c = _solve(cols, _vec(Mom(**{v: 1}), basis))
                terms = []
                for cj, j in zip(c[:len(ind)], ind):
                    if cj: terms.append(f'({cj})*(DD({j})-DD({b})-{dot(R[j], R[j])}+{ms[j]}-{ms[b]})/2')
                for cw, w in zip(c[len(ind):], comp):
                    if cw: terms.append(f'({cw})*l.{w}')
                rules.append(f'      id,once l.{v} = {"+".join(terms) if terms else "0"};')
            body = f'  if ( match(TOPO({mask})) );\n    if ( count(l,1) > {LIM[len(mem)]} );\n'
            first = True
            for r, v in zip(rules, ['l'] + [w for w in basis if w not in comp]):
                pat = 'l.l' if v == 'l' else f'l.{v}'
                body += f'      {"if" if first else "elseif"} ( match({pat}) );\n  {r}\n        $red = 1;\n'
                first = False
            body += '      else;\n        multiply STUCK;\n      endif;\n    endif;\n  endif;\n  id STUCK^n1?{>1} = STUCK;\n'
            red.append(body)
            for i in mem:
                rest = mask & ~(1 << i)
                if rest == 0:
                    canc.append(f'  id DD({i})*TOPO({mask}) = 0;')
                elif i == b:
                    b2 = min(j for j in range(n) if rest >> j & 1)
                    canc.append(f'  id DD({i})*TOPO({mask}) = TOPO({rest})*SHIFT({b},{b2});')
                    shifts[(b, b2)] = K[b2] - K[b]
                else:
                    canc.append(f'  id DD({i})*TOPO({mask}) = TOPO({rest});')
            kk = [R[i] for i in mem[1:]]
            ks = '\n'.join(f'    id k{j+1} = {kk[j].form(lambda v: v) if kk[j] else "0"};' if j < len(kk)
                           else f'    id k{j+1} = 0;' for j in range(3))
            ksub.append(f'  if ( match(TOPO({mask})) );\n{ks}\n  endif;')
        sh = '\n'.join(f'  if ( match(SHIFT({b},{b2})) );\n    Multiply replace_(l,l-{v.form(lambda w: w)});\n'
                       f'    id SHIFT({b},{b2}) = 1;\n  endif;' for (b, b2), v in shifts.items())
        kin = {'compton': '#call Kin', 'vertex': '#call KinV', 'self': '#call KinS'}[self.mode]
        nit = 20
        proc = '#procedure REDUCE\n' + f'  Multiply TOPO({2 ** n - 1});\n  .sort\n' + \
            f'#do it=1,{nit}\n  #$red = 0;\n' + ''.join(red) + '  ModuleOption maximum,$red;\n  .sort\n' + \
            '\n'.join(canc) + '\n  .sort\n' + sh + f'\n  {kin}\n  .sort\n  #if `$red\' == 0\n    #breakdo\n  #endif\n' + \
            '#enddo\n#endprocedure\n'
        proc += '#procedure KSUBT\n' + '\n'.join(ksub) + '\n#endprocedure\n'
        return proc

    def form_body(self, gi=False):
        global AUX
        toks, scal, ks = self.shifted()
        tl = remove_g5(toks)
        AUX = {} if AUX_ON else None
        try:
            string = '*'.join(tk.text(f) for tk, f in tl) if tl else 'gi_(1)'
            sc = [dot(it[0], it[1]) if isinstance(it, tuple) else it for it in scal]
            aux = dict(AUX.values()) if AUX_ON else {}
        finally:
            AUX = None
        rep = ','.join(f'{nm},{m.form(lambda v: v)[1:-1]}' for nm, m in aux.items())
        auxsub = '#procedure AUXSUB\n' + (f'  Multiply replace_({rep});\n' if aux else '') + '#endprocedure\n'
        pieces = [x for x in ['*'.join(self.tree), '*'.join(sc), string] if x]
        expr = '*'.join(pieces)
        pf = {'compton': '#call Proj\n#call FinishR\n', 'vertex': '#call VProj\n#call VFinishR\n',
              'self': '#call SFinishR\n'}[self.mode]
        if gi:
            import re
            expr = re.sub(r'\bee\b', 'q', expr)
            pf = '#call GProj\n#call GFinishR\n'
        pre = '#define RS4D\n' if RS4D else ''
        return pre + self.reduce_proc() + (auxsub if AUX_ON else '') + f'L GAM = {expr};\n.sort\n{pf}', ks

    def topo_args(self, mask):
        """LoopTools kind, invariants (sympy), masses for the sub-topology mask"""
        n = len(self.loop)
        K = self.Ks()
        mem = [i for i in range(n) if mask >> i & 1]
        kk = [K[i] - K[mem[0]] for i in mem]
        ms = [self.loop[i][1] for i in mem]
        sq = lambda m: dot_sym(m, m, self.mode)
        if len(mem) == 1: return ('A', [], ms)
        if len(mem) == 2: return ('B', [sq(kk[1])], ms)
        if len(mem) == 3: return ('C', [sq(kk[1]), sq(kk[2] - kk[1]), sq(kk[2])], ms)
        return ('D', [sq(kk[1]), sq(kk[2] - kk[1]), sq(kk[3] - kk[2]), sq(kk[3]), sq(kk[2]), sq(kk[3] - kk[1])], ms)

def _rank(vecs):
    import sympy as sp
    if not vecs: return 0
    return sp.Matrix([[sp.Rational(c.numerator, c.denominator) for c in v] for v in vecs]).rank()
