"""build_fortran.py -- generate work/n3lo/fortran/p4/*.F90: the O(p^4) piN-loop Compton amplitude for code_v2.

Content (all validated in the Python check path, tests/t_*.py):
  kappa insertions: piN.all_graphs(4) + their renormalisation subtraction graphs (on-shell scheme of LP2010)
  c1, c2, c3 graphs: tchan.graphs()
Output: SUBROUTINE CovOneLoopP4(x,t,A) -> A(1..6) in the basis of hel_from_A (A7=A8=0), nucleon-mass units,
alpha_em stripped (same normalisation as CovOneLoop*gafac).  SUBROUTINE p4_setup(gA,F,mpi2,kp,kn,c1,c2,c3)
must be called once (computes the per-graph factors and the on-shell renormalisation constants).
"""
import os, re, sys, concurrent.futures as cf, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import piN, tchan, fortgen, formtools as ft
from graphs import Graph, Tok
from evalamp import coef_name


def write_if_changed(path, text):
    """keep the timestamp of unchanged sources, so that make recompiles only what changed"""
    if not (os.path.exists(path) and open(path).read() == text):
        open(path, 'w').write(text)
OUT = os.path.join(HERE, '..', 'fortran', 'p4'); os.makedirs(OUT, exist_ok=True)

LTIDX = {}
for line in open(os.path.join(os.environ.get('LT', os.path.expanduser('~/.local/looptools-dq')), 'include', 'looptools.h')):
    m = re.match(r'#define ((?:aa|bb|cc|dd)[0-9]+) (\d+)', line)
    if m: LTIDX[m.group(1)] = int(m.group(2))

def lt_kind(g):
    return {0: 'T', 1: 'A', 2: 'B', 3: 'C', 4: 'D'}[len(g.loop)]

def fexpr(e):
    s = sp.fcode(sp.nsimplify(e), source_format='free', standard=95)
    return s.replace('d0', 'd0')

def graph_sub(g, nproj, mode):
    """Fortran subroutine for graph g: computes W(1..nproj) (complex), arguments x,t,dl,ss,F1s,F2s,L0s,L1s"""
    kind = lt_kind(g)
    code = fortgen.optimise(g, kind, nproj)
    used = sorted(set(re.findall(r'\b([TU](?:aa|bb|cc|dd)[0-9]+)\b', code)))
    nz = max([int(z) for z in re.findall(r'\bZ(\d+)_', code)] + [0])
    lines = [f'SUBROUTINE p4g_{g.name}(x,t,dl,ss,mpi2,F1s,F2s,L0s,L1s,W)',
             '  IMPLICIT NONE', '  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2',
             '  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s', f'  COMPLEX*16, INTENT(out) :: W({nproj})',
             '  REAL*8 ixs,iys,iPq,PP,p1,p2,p3,p4,p12,p23,ma(4)']
    if nz: lines.append(f'  COMPLEX*16 Z({nz})')
    if used: lines.append('  COMPLEX*16 ' + ','.join(used))
    lines += ['  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i']
    lines += ['  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t']
    if kind != 'T':
        _, inv, ms = g.lt_args()
        names = {'A': [], 'B': ['p1'], 'C': ['p1', 'p2', 'p12'], 'D': ['p1', 'p2', 'p3', 'p4', 'p12', 'p23']}[kind]
        for nm, e in zip(names, inv):
            lines.append(f'  {nm} = {fexpr(e)}')
        for i, m2 in enumerate(ms):
            lines.append(f'  ma({i+1}) = ' + ('mpi2' if m2 == 'mpi2' else '1.d0'))
        args = {'A': None, 'B': 'p1,ma(1),ma(2)', 'C': 'p1,p2,p12,ma(1),ma(2),ma(3)',
                'D': 'p1,p2,p3,p4,p12,p23,ma(1),ma(2),ma(3),ma(4)'}[kind]
        fn = {'B': 'B0i', 'C': 'C0i', 'D': 'D0i'}.get(kind)
        for u in used:
            nm = u[1:]
            off = '+1' if u[0] == 'U' else ''
            if kind == 'A':
                lines.append(f'  {u} = ma(1)*B0i({LTIDX["bb0"]}{off},0.d0,0.d0,ma(1))')
            else:
                lines.append(f'  {u} = {fn}({LTIDX[nm]}{off},{args})')
    body = code.replace('Z', 'Z')
    body = re.sub(r'\bZ(\d+)_', r'Z(\1)', body)
    lines.append(body)
    lines.append(f'END SUBROUTINE p4g_{g.name}\n')
    src = '\n'.join(lines)
    src = src.replace('\n     &', '&\n     &')         # fixed-form continuation -> free form (split tokens allowed)
    return src, used

def tree_graph(name, gam, mode='compton'):
    g = Graph(name); g.mode = mode
    g.toks = []; g.loop = []
    g._gam = gam
    return g

def all_graphs():
    comp = piN.all_graphs(4)
    subs, verts = [], {}
    for g in comp:
        if g.kind == 'R2':
            h = piN.sub_graph_R2(g); h.parent = g; subs.append(h)
            if h.vkey not in verts:
                slot, c = h.vkey
                v = piN.vertex_graph(slot, c, f'vtx_{slot}_{c}')
                v.iso = ({'kp': piN.isospin(v, piN.Km('p')), 'kn': piN.isospin(v, piN.Km('n'))} if c == 'mag'
                         else {'1': piN.isospin(v, None)})
                verts[h.vkey] = v
        elif g.kind == 'R1':
            h = piN.sub_graph_R1(g); h.parent = g; subs.append(h)
    sg = piN.self_graph('self_piN'); sg.iso = {'1': piN.isospin(sg, None)}
    return comp, subs, verts, sg, tchan.graphs()

def main(nproc=8):
    comp, subs, verts, sg, tg = all_graphs()
    jobs = [(g, 6, 'compton') for g in comp + subs + tg] + [(v, 3, 'vertex') for v in verts.values()] + [(sg, 2, 'self')]
    with cf.ThreadPoolExecutor(nproc) as ex:
        res = list(ex.map(lambda j: graph_sub(*j), jobs))
    for (g, n, mode), (src, used) in zip(jobs, res):
        write_if_changed(os.path.join(OUT, f'p4g_{g.name}.F90'), src)
    # gram matrices (tree projections of the structures)
    write_driver(comp, subs, verts, sg, tg)
    print('graphs:', len(jobs))

def iso_code(iso):
    t = []
    for k, v in iso.items():
        v = complex(v)
        c = f'({v.real:.17g}d0,{v.imag:.17g}d0)'
        t.append(c if k == '1' else f'{c}*{k}')
    return '+'.join(t) if t else '(0.d0,0.d0)'

def coup_code(g):
    t = []
    for k, v in g.coup.items():
        t.append(f'{k}**({v})')
    return '*'.join(t) if t else '1.d0'

def cplx(v):
    v = complex(v); return f'({v.real:.17g}d0,{v.imag:.17g}d0)'

def write_driver(comp, subs, verts, sg, tg):
    S = ['(-ep.ee)*gi_(1)', '(q.ep*qp.ee)*gi_(1)', '(-(g_(1,ep,ee)-ep.ee))', '(ep.ee*(g_(1,qp,q)-qp.q))',
         '(q.ep*(g_(1,qp,ee)-qp.ee)-qp.ee*(g_(1,q,ep)-q.ep))', '(q.ep*(g_(1,q,ee)-q.ee)-qp.ee*(g_(1,qp,ep)-qp.ep))']
    R = ['g_(1,ee)', 'P.ee*gi_(1)', 'q.ee*gi_(1)']
    KS = '#procedure KSUB\n#endprocedure\n'
    gramc = []
    for k in range(6):
        d = ft.run_form(f'gram{k+1}', KS + f'L GAM = {S[k]};\n.sort\n#call Proj\n#call Finish\n')
        gramc.append([ft.to_python(d[f'V{j}']) for j in range(1, 7)])
    vgr = []
    for k in range(3):
        d = ft.run_form(f'vgram{k+1}', KS + f'L GAM = {R[k]};\n.sort\n#call VProj\n#call VFinish\n')
        vgr.append([ft.to_python(d[f'V{j}']) for j in range(1, 4)])
    def tofort(e):
        e = e.replace('T[(0,0,0,0)]', '1.d0')
        e = re.sub(r'(\d+)\.0\b', r'\1.d0', e)
        e = re.sub(r'(?<![\w.])(\d+)\.(?![\dd])', r'\1.d0', e)
        return e
    L = []
    L.append('! CovOneLoopP4.F90 -- generated by ../n3lo/py/build_fortran.py; do not edit by hand.')
    L.append('! O(p^4) piN-loop Compton amplitude of covariant BChPT (kappa insertions, c1, c2, c3 graphs),')
    L.append('! on-shell renormalised (physical M, e, kappa), MSbar (Delta = 0, mu = M) for the c2 divergence.')
    L.append('! A(1:6) in the basis of hel_from_A (A7=A8=0); nucleon-mass units; alpha_em stripped.')
    L.append('MODULE p4mod\n  IMPLICIT NONE')
    L.append('  REAL*8 :: gA_p4, F_p4, mpi2_p4, kp, kn, c1, c2, c3')
    L.append('  COMPLEX*16 :: vc1(8), vc2(8), sL0, sL1')
    L.append(f'  COMPLEX*16 :: fac({len(comp)+len(subs)+len(tg)})')
    L.append('  LOGICAL :: p4_ready = .false.')
    L.append('  INTEGER :: p4_part = 0   ! 0 all, 1 kappa insertions, 2 c1, 3 c2, 4 c3 (diagnostics)')
    L.append('END MODULE p4mod\n')
    vlist = list(verts.items())
    vindex = {k: i + 1 for i, (k, _) in enumerate(vlist)}
    # setup
    L.append('SUBROUTINE p4_setup(gA_,F_,mpi2_,kp_,kn_,c1_,c2_,c3_)\n  USE p4mod\n  IMPLICIT NONE')
    L.append('  REAL*8, INTENT(in) :: gA_,F_,mpi2_,kp_,kn_,c1_,c2_,c3_')
    L.append('  COMPLEX*16 W(6), Wv(3), cA(3,3), Y(2), ab(2,5)')
    L.append('  REAL*8 G3(3,3), dls(3), sss(5), dlv, h, gA, F')
    L.append('  COMPLEX*16 Gc(3,3), Wc(3)')
    L.append('  INTEGER i, k, ipiv(3), info, vk')
    L.append('  REAL*8 iPq, PP')
    L.append('  INTEGER*4 oldkey')
    L.append('  INTEGER, EXTERNAL :: getversionkey')
    L.append('  gA_p4=gA_; F_p4=F_; mpi2_p4=mpi2_; kp=kp_; kn=kn_; c1=c1_; c2=c2_; c3=c3_; gA=gA_; F=F_')
    L.append('  oldkey=getversionkey(); CALL setversionkey(64)   ! Denner D0 (see CovOneLoopP4)')
    L.append('  dls(1)=1.d-4; dls(2)=2.d-4; dls(3)=3.d-4   ! quadratic extrapolation to q^2 = 0')
    for (key, v), i in zip(vlist, range(1, len(vlist) + 1)):
        fv = f'({cplx(v.pref)})*{coup_code(v)}*({iso_code(v.iso)})'
        L.append(f'  ! on-shell vertex constants: {v.name}')
        L.append('  DO k=1,3')
        L.append('    dlv=dls(k)')
        L.append(f'    CALL p4g_{v.name}(0.d0,0.d0,dlv,0.d0,mpi2_p4,(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),Wv)')
        L.append(f'    Wv=Wv*{fv}')
        L.append('    CALL vgram(dlv,G3)')
        L.append('    Gc=G3; Wc=Wv')
        L.append('    CALL zsolve3(Gc,Wc)')
        L.append('    cA(k,1)=Wc(1)+2.d0*Wc(2); cA(k,2)=Wc(2)')
        L.append('  ENDDO')
        L.append(f'  vc1({i})=3.d0*(cA(1,1)-cA(2,1))+cA(3,1); vc2({i})=3.d0*(cA(1,2)-cA(2,2))+cA(3,2)')
    fs = f'({cplx(sg.pref)})*{coup_code(sg)}*({iso_code(sg.iso)})'
    L.append('  ! self-energy constants L(M), L\'(M), L(p) = a p/ + b')
    L.append('  h=1.d-4; sss(1)=1.d0; sss(2)=1.d0+h; sss(3)=1.d0-h; sss(4)=1.d0+2.d0*h; sss(5)=1.d0-2.d0*h')
    L.append('  DO k=1,5')
    L.append(f'    CALL p4g_{sg.name}(0.d0,0.d0,0.d0,sss(k),mpi2_p4,(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),Y)')
    L.append(f'    Y=Y*{fs}')
    L.append('    ab(1,k)=Y(2)/sss(k); ab(2,k)=Y(1)')
    L.append('  ENDDO')
    L.append('  sL0=ab(1,1)+ab(2,1)')
    L.append('  sL1=ab(1,1)+2.d0*(8.d0*(ab(1,2)-ab(1,3)+ab(2,2)-ab(2,3))-(ab(1,4)-ab(1,5)+ab(2,4)-ab(2,5)))/(12.d0*h)')
    allc = comp + subs + tg
    for i, g in enumerate(allc, 1):
        if g in subs:
            if hasattr(g, 'vkey'): f = '-1.d0' if g.Ocpl == 'ch' else '-kp'
            else:
                c1_, c2_ = g.cpls
                f = '-' + '*'.join(['1.d0' if c == 'ch' else 'kp' for c in (c1_, c2_)])
            L.append(f'  fac({i}) = ({cplx(g.pref)})*({f})')
        else:
            L.append(f'  fac({i}) = ({cplx(g.pref)})*{coup_code(g)}*({iso_code(g.iso)})')
    L.append('  CALL setversionkey(oldkey)')
    L.append('  p4_ready=.true.\nEND SUBROUTINE p4_setup\n')
    # amplitude
    L.append('SUBROUTINE CovOneLoopP4(x,t,A)\n  USE p4mod\n  IMPLICIT NONE')
    L.append('  REAL*8, INTENT(in) :: x,t\n  COMPLEX*16, INTENT(out) :: A(8)')
    L.append('  COMPLEX*16 W(6), V(6), Gc(6,6), Z0\n  REAL*8 G6(6,6)\n  INTEGER*4 oldkey, ip\n  INTEGER, EXTERNAL :: getversionkey')
    L.append('  Z0=(0.d0,0.d0); V=Z0')
    L.append('  IF(.NOT.p4_ready) STOP \'CovOneLoopP4: call p4_setup first\'')
    L.append('  oldkey=getversionkey(); CALL setversionkey(64)   ! Denner D0: FF fails for D(0,M2,0,M2;u,s;m,m,M,M) below threshold')
    for i, g in enumerate(allc, 1):
        part = 1
        if g in tg: part = {'c1': 2, 'c2': 3, 'c3': 4}[[k for k in g.coup if k.startswith('c')][0]]
        args = 'Z0,Z0,Z0,Z0'
        if g in subs:
            if hasattr(g, 'vkey'): args = f'vc1({vindex[g.vkey]}),vc2({vindex[g.vkey]}),Z0,Z0'
            else: args = 'Z0,Z0,sL0,sL1'
        L.append(f'  IF(p4_part.EQ.0.OR.p4_part.EQ.{part}) THEN')
        L.append(f'    CALL p4g_{g.name}(x,t,0.d0,0.d0,mpi2_p4,{args},W); V=V+W*fac({i})')
        L.append('  ENDIF')
    L.append('  CALL setversionkey(oldkey)')
    L.append('  CALL cgram(x,t,G6)\n  Gc=G6\n  CALL zsolve6(Gc,V)\n  A=(0.d0,0.d0); A(1:6)=V')
    L.append('END SUBROUTINE CovOneLoopP4\n')
    # gram routines
    L.append('SUBROUTINE cgram(x,t,G)\n  IMPLICIT NONE\n  REAL*8, INTENT(in) :: x,t\n  REAL*8, INTENT(out) :: G(6,6)')
    L.append('  REAL*8 ixs,iys,iPq,PP\n  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t')
    for k in range(6):
        for j in range(6):
            L.append(f'  G({j+1},{k+1}) = {tofort(gramc[k][j])}')
    L.append('END SUBROUTINE cgram\n')
    L.append('SUBROUTINE vgram(dl,G)\n  IMPLICIT NONE\n  REAL*8, INTENT(in) :: dl\n  REAL*8, INTENT(out) :: G(3,3)\n  REAL*8 M\n  M=1.d0')
    for k in range(3):
        for j in range(3):
            L.append(f'  G({j+1},{k+1}) = {tofort(vgr[k][j])}')
    L.append('END SUBROUTINE vgram\n')
    for n in (3, 6):
        L.append(f'''SUBROUTINE zsolve{n}(A,b)
  ! Gaussian elimination with partial pivoting, b <- A^-1 b
  IMPLICIT NONE
  COMPLEX*16, INTENT(inout) :: A({n},{n}), b({n})
  COMPLEX*16 f, tmp({n})
  INTEGER i, j, k, p
  DO k=1,{n}
     p=k
     DO i=k+1,{n}
        IF(ABS(A(i,k)).GT.ABS(A(p,k))) p=i
     ENDDO
     IF(p.NE.k) THEN
        tmp=A(k,:); A(k,:)=A(p,:); A(p,:)=tmp
        f=b(k); b(k)=b(p); b(p)=f
     ENDIF
     DO i=k+1,{n}
        f=A(i,k)/A(k,k)
        A(i,k:{n})=A(i,k:{n})-f*A(k,k:{n})
        b(i)=b(i)-f*b(k)
     ENDDO
  ENDDO
  DO k={n},1,-1
     b(k)=(b(k)-SUM(A(k,k+1:{n})*b(k+1:{n})))/A(k,k)
  ENDDO
END SUBROUTINE zsolve{n}
''')
    write_if_changed(os.path.join(OUT, 'CovOneLoopP4.F90'), '\n'.join(L))

if __name__ == '__main__':
    main()
