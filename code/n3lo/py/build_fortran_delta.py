"""build_fortran_delta.py -- generate work/n3lo/fortran/dl/*.F90: Delta-loop Compton amplitudes beyond NNLO for code_v2.

Content (validated in the Python check path, tests/t_lpdelta.py, t_dA.py, t_dB.py):
  CovDeltaFull(x,t,A) : pi Delta loops with the complete minimal coupling (photon on the Delta charge, Delta-derivative
                        contacts, true isospin factors), on-shell renormalised, MS-bar (Delta=0, mu=M) otherwise.
                        Replaces CorrIsoCovOneDeltaLoop (s+u) when switched on.
  CovDeltaGND(x,t,A)  : gamma N Delta (gM, gE) vertex inside the pi N loops, second photon minimally coupled
                        everywhere, kappa shift subtracted (on-shell), MS-bar otherwise.
  dl_setup(gA,F,hA,MD,gM,gE,mpi2): per-graph factors and on-shell constants (call once).
A(1..6) in the basis of hel_from_A (A7=A8=0), nucleon-mass units, alpha_em stripped (as CovOneLoopP4).
The Gram routines cgram, vgram, zsolve3, zsolve6 are those of p4/CovOneLoopP4.F90.
"""
import os, re, sys, subprocess, hashlib, concurrent.futures as cf, numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import formtools as ft, fortgen, devalamp as dv, ddelta, delta
from build_fortran import LTIDX, fexpr, iso_code, coup_code, cplx, write_if_changed
from evalamp import coef_name
OUT = os.path.join(HERE, '..', 'fortran', 'dl'); os.makedirs(OUT, exist_ok=True)
FDIR = fortgen.FDIR

def tsymd(mask, kind, key, uv=False):
    return ('U' if uv else 'T') + str(mask) + coef_name(kind, key)

def _fortran_to_form(e):
    """FORM Fortran output (Format doublefortran) -> FORM input"""
    e = re.sub(r'(\d+)\.D0', r'\1', e)
    e = re.sub(r'(\d+)\.(?!\d)', r'\1', e)
    e = e.replace('**', '^')
    return e

def optimise_d(g, nproj, level='O2'):
    """FORM-optimised Fortran for a DGraph, starting from the cached projections (TCFT(mask,...) -> named symbols)"""
    _, d = dv.project(g)
    keys = set()
    for v in d.values():
        for m in re.findall(r'[TU]CFT\((\d+),(\d+),(\d+),(\d+),(\d+)\)', v):
            keys.add(tuple(int(z) for z in m))
    kinds = {m: g.topo_args(m)[0] for m in {k[0] for k in keys}}
    names = sorted({tsymd(k[0], kinds[k[0]], k[1:]) for k in keys} | {tsymd(k[0], kinds[k[0]], k[1:], True) for k in keys})
    decl = ('Symbols ' + ','.join(names) + ';\n') if names else ''
    ids = ''.join(f'id TCFT({k[0]},{k[1]},{k[2]},{k[3]},{k[4]}) = {tsymd(k[0], kinds[k[0]], k[1:])};\n'
                  f'id UCFT({k[0]},{k[1]},{k[2]},{k[3]},{k[4]}) = {tsymd(k[0], kinds[k[0]], k[1:], True)};\n' for k in keys)
    ids += 'id TCF(0,0,0,0) = 1;\n'
    src = '#include dcommon.h\n' + decl
    for j in range(1, nproj + 1):
        src += f'L V{j} = {_fortran_to_form(d.get(f"V{j}", "0"))};\n'
    src += ids + '.sort\n'
    src += f'Format {level},stats=off;\nFormat doublefortran;\nFormat 120;\n'
    out = os.path.join(FDIR, g.name + '.f90.inc')
    for j in range(1, nproj + 1):
        src += f'#optimize V{j}\n#write <{out}> "%O"\n#write <{out}> "      W({j}) = %e", V{j}\n#clearoptimize\n'
    src += '.end\n'
    h = hashlib.md5((src + open(os.path.join(ft.FORMDIR, 'dcommon.h')).read()).encode()).hexdigest()[:10]
    stamp = out + '.hash'
    if os.path.exists(out) and os.path.exists(stamp) and open(stamp).read() == h:
        return open(out).read(), keys, kinds
    if os.path.exists(out): open(out, 'w').close()
    fn = os.path.join(FDIR, g.name + '.frm'); open(fn, 'w').write(src)
    r = subprocess.run(['form', '-q', '-p', ft.FORMDIR, fn], capture_output=True, text=True, cwd=FDIR)
    if r.returncode != 0 or 'rror' in r.stdout: raise RuntimeError(g.name + ':\n' + r.stdout[-3000:])
    open(stamp, 'w').write(h)
    return open(out).read(), keys, kinds

_SPLIT = re.compile(r'(?<![\w.])(\d+)\.(?![\dDdEeQq_\w])')

def dbl(code):
    """FORM's optimised output writes rational constants as '4./3.': default-real (single precision) literals in
    Fortran.  Make them double precision (4.D0/3.D0); dyadic ones were exact anyway."""
    return _SPLIT.sub(r'\1.D0', code)

def graph_sub_d(g, nproj):
    code, keys, kinds = optimise_d(g, nproj)
    code = dbl(code)
    nz = max([int(z) for z in re.findall(r'\bZ(\d+)_', code)] + [0])
    used = sorted({tsymd(k[0], kinds[k[0]], k[1:]) for k in keys} | {tsymd(k[0], kinds[k[0]], k[1:], True) for k in keys})
    L = [f'SUBROUTINE dlg_{g.name}(x,t,dl,ss,mpi2,MD,F1s,F2s,L0s,L1s,W)', '  IMPLICIT NONE',
         '  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2,MD', '  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s',
         f'  COMPLEX*16, INTENT(out) :: W({nproj})',
         '  REAL*8 ixs,iys,iPq,PP,ixD,iyD,iMD,y,p1,p2,p3,p4,p12,p23,ma(4)']
    if nz: L.append(f'  COMPLEX*16 Z({nz})')
    if used: L.append('  COMPLEX*16 ' + ','.join(used))
    L += ['  COMPLEX*16, EXTERNAL :: B0i,C0i,D0i', '  COMPLEX*16 a0']
    if 'STUCK' in code: L += ['  REAL*8 STUCK', '  STUCK=1.d0          ! flags terms left at rank 4 (bubbles, dl_bhigh)']
    if g.mode == 'compton':
        L += ['  y=-x-0.5d0*t', '  ixs=0.5d0/x; iys=0.5d0/y; iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t',
              '  ixD=1.d0/(1.d0+2.d0*x-MD**2); iyD=1.d0/(1.d0+2.d0*y-MD**2); iMD=1.d0/(1.d0-MD**2)']
    for m in sorted(kinds):
        kind, inv, ms = g.topo_args(m)
        nm = {'A': [], 'B': ['p1'], 'C': ['p1', 'p2', 'p12'], 'D': ['p1', 'p2', 'p3', 'p4', 'p12', 'p23']}[kind]
        L.append(f'  ! topology {m} ({kind})')
        for n_, e in zip(nm, inv): L.append(f'  {n_} = {fexpr(e)}')
        for i, m2 in enumerate(ms):
            L.append(f'  ma({i+1}) = ' + {'mpi2': 'mpi2', '1': '1.d0', 'MD2': 'MD**2'}[m2])
        args = {'B': 'p1,ma(1),ma(2)', 'C': 'p1,p2,p12,ma(1),ma(2),ma(3)',
                'D': 'p1,p2,p3,p4,p12,p23,ma(1),ma(2),ma(3),ma(4)'}.get(kind)
        fn = {'B': 'B0i', 'C': 'C0i', 'D': 'D0i'}.get(kind)
        mk = sorted({k[1:] for k in keys if k[0] == m})
        if kind == 'A':
            L.append(f'  a0 = ma(1)*B0i({LTIDX["bb0"]},0.d0,0.d0,ma(1))')
        for key in mk:
            cn = coef_name(kind, key); tn = tsymd(m, kind, key); un = tsymd(m, kind, key, True)
            if kind == 'A':
                fin = {'aa0': 'a0', 'aa00': 'ma(1)*a0/4.d0+ma(1)**2/8.d0',
                       'aa0000': 'ma(1)**2*a0/24.d0+5.d0*ma(1)**3/144.d0'}[cn]
                uv = {'aa0': 'ma(1)', 'aa00': 'ma(1)**2/4.d0', 'aa0000': 'ma(1)**3/24.d0'}[cn]
                L.append(f'  {tn} = {fin}'); L.append(f'  {un} = {uv}')
            elif kind == 'B' and cn not in LTIDX:     # rank > 3 bubbles: Feynman-parameter form (bhigh.py)
                L.append(f'  CALL dl_bhigh({key[0] // 2},{key[1]},p1,ma(1),ma(2),{tn},{un})')
            else:
                L.append(f'  {tn} = {fn}({LTIDX[cn]},{args})'); L.append(f'  {un} = {fn}({LTIDX[cn]}+1,{args})')
    body = re.sub(r'\bZ(\d+)_', r'Z(\1)', code)
    L.append(body)
    L.append(f'END SUBROUTINE dlg_{g.name}\n')
    return '\n'.join(L).replace('\n     &', '&\n     &')

def graph_sub_old(h, nproj=6):
    """tree subtraction graphs (graphs.Graph, common.h): same routine shape"""
    code = dbl(fortgen.optimise(h, 'T', nproj))
    nz = max([int(z) for z in re.findall(r'\bZ(\d+)_', code)] + [0])
    L = [f'SUBROUTINE dlg_{h.name}(x,t,dl,ss,mpi2,MD,F1s,F2s,L0s,L1s,W)', '  IMPLICIT NONE',
         '  REAL*8, INTENT(in) :: x,t,dl,ss,mpi2,MD', '  COMPLEX*16, INTENT(in) :: F1s,F2s,L0s,L1s',
         f'  COMPLEX*16, INTENT(out) :: W({nproj})', '  REAL*8 ixs,iys,iPq,PP']
    if nz: L.append(f'  COMPLEX*16 Z({nz})')
    L.append('  ixs=0.5d0/x; iys=0.5d0/(-x-0.5d0*t); iPq=1.d0/(2.d0*x+0.5d0*t); PP=4.d0-t')
    L.append(re.sub(r'\bZ(\d+)_', r'Z(\1)', code))
    L.append(f'END SUBROUTINE dlg_{h.name}\n')
    return '\n'.join(L).replace('\n     &', '&\n     &')

def model_code(tag, m, routine):
    """setup fragment and amplitude routine for one DModel"""
    vlist = list(m.vert.items()); slist = list(m.selfe.items())
    vidx = {k: i + 1 for i, (k, _) in enumerate(vlist)}; sidx = {k: i + 1 for i, (k, _) in enumerate(slist)}
    nv, ns, ng = max(len(vlist), 1), max(len(slist), 1), len(m.graphs) + len(m.subs)
    decl = [f'  COMPLEX*16 :: {tag}_vc1({nv}), {tag}_vc2({nv}), {tag}_sL0({ns}), {tag}_sL1({ns}), {tag}_fac({ng})']
    S = [f'  ! ---- {routine}: on-shell vertex constants']
    for (key, v), i in zip(vlist, range(1, len(vlist) + 1)):
        fv = f'({cplx(v.pref)})*{coup_code(v)}*({iso_code(v.iso)})'
        S += ['  DO k=1,3', '    dlv=dls(k)',
              f'    CALL dlg_{v.name}(0.d0,0.d0,dlv,0.d0,mpi2_dl,MD_dl,(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),Wv)',
              f'    Wv=Wv*{fv}', '    CALL vgram(dlv,G3)', '    Gc=G3; Wc=Wv', '    CALL zsolve3(Gc,Wc)',
              '    cA(k,1)=Wc(1)+2.d0*Wc(2); cA(k,2)=Wc(2)', '  ENDDO',
              f'  {tag}_vc1({i})=3.d0*(cA(1,1)-cA(2,1))+cA(3,1); {tag}_vc2({i})=3.d0*(cA(1,2)-cA(2,2))+cA(3,2)']
    for (key, s), i in zip(slist, range(1, len(slist) + 1)):
        fs = f'({cplx(s.pref)})*{coup_code(s)}*({iso_code(s.iso)})'
        S += ['  h=1.d-4; sss(1)=1.d0; sss(2)=1.d0+h; sss(3)=1.d0-h; sss(4)=1.d0+2.d0*h; sss(5)=1.d0-2.d0*h', '  DO k=1,5',
              f'    CALL dlg_{s.name}(0.d0,0.d0,0.d0,sss(k),mpi2_dl,MD_dl,(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),Y)',
              f'    Y=Y*{fs}', '    ab(1,k)=Y(2)/sss(k); ab(2,k)=Y(1)', '  ENDDO',
              f'  {tag}_sL0({i})=ab(1,1)+ab(2,1)',
              f'  {tag}_sL1({i})=ab(1,1)+2.d0*(8.d0*(ab(1,2)-ab(1,3)+ab(2,2)-ab(2,3))-(ab(1,4)-ab(1,5)+ab(2,4)-ab(2,5)))/(12.d0*h)']
    allg = list(m.graphs) + list(m.subs)
    for i, g in enumerate(allg, 1):
        if g in m.subs: S.append(f'  {tag}_fac({i}) = ({cplx(g.pref)})*(-1.d0)')
        else: S.append(f'  {tag}_fac({i}) = ({cplx(g.pref)})*{coup_code(g)}*({iso_code(g.iso)})')
    A = [f'SUBROUTINE {routine}(x,t,A)', '  USE dlmod', '  IMPLICIT NONE', '  REAL*8, INTENT(in) :: x,t',
         '  COMPLEX*16, INTENT(out) :: A(8)', '  COMPLEX*16 W(6), V(6), Gc(6,6), Z0', '  REAL*8 G6(6,6)',
         '  INTEGER*4 oldkey', '  INTEGER, EXTERNAL :: getversionkey', '  Z0=(0.d0,0.d0); V=Z0',
         f"  IF(.NOT.dl_ready) STOP '{routine}: call dl_setup first'", '  oldkey=getversionkey(); CALL setversionkey(64)']
    for i, g in enumerate(allg, 1):
        args = 'Z0,Z0,Z0,Z0'
        if g in m.subs:
            args = (f'{tag}_vc1({vidx[g.vkey]}),{tag}_vc2({vidx[g.vkey]}),Z0,Z0' if hasattr(g, 'vkey')
                    else f'Z0,Z0,{tag}_sL0({sidx[g.skey]}),{tag}_sL1({sidx[g.skey]})')
        A.append(f'  CALL dlg_{g.name}(x,t,0.d0,0.d0,mpi2_dl,MD_dl,{args},W); V=V+W*{tag}_fac({i})')
    A += ['  CALL setversionkey(oldkey)', '  CALL cgram(x,t,G6)', '  Gc=G6', '  CALL zsolve6(Gc,V)',
          '  A=(0.d0,0.d0); A(1:6)=V', f'END SUBROUTINE {routine}\n']
    jobs = [(g, 6) for g in m.graphs] + [(v, 3) for v in m.vert.values()] + [(s, 2) for s in m.selfe.values()]
    return decl, S, A, jobs, list(m.subs)

def main(nproc=8, which=('af', 'bl')):
    cs = dv.dconsts()
    D = ddelta.Delta45(cs, nproc=nproc, which=tuple({'af': 'A', 'bl': 'B'}[w] for w in which))
    parts = {k: v for k, v in {'af': ('A_full', 'CovDeltaFull'), 'bl': ('Bloop', 'CovDeltaGND')}.items() if k in which}
    parts = {k: (D.m[m], r) for k, (m, r) in parts.items()}
    decls, setup, amps, jobs, subs = [], [], [], [], []
    for tag, (m, routine) in parts.items():
        d_, s_, a_, j_, sb_ = model_code(tag, m, routine)
        decls += d_; setup += s_; amps += a_; jobs += j_; subs += sb_
    with cf.ThreadPoolExecutor(nproc) as ex:
        srcs = list(ex.map(lambda j: graph_sub_d(*j), jobs))
    srcs += [graph_sub_old(h) for h in subs]
    names = [j[0].name for j in jobs] + [h.name for h in subs]
    for nm, src in zip(names, srcs):
        write_if_changed(os.path.join(OUT, f'dlg_{nm}.F90'), src)
    L = ['! CovDeltaLoops.F90 -- generated by ../n3lo/py/build_fortran_delta.py; do not edit by hand.',
         '! Delta-loop Compton amplitudes beyond NNLO: CovDeltaFull (complete minimal coupling of the pi Delta loops),',
         '! CovDeltaGND (gamma N Delta vertex inside the pi N loops).  See the generator for conventions.',
         'MODULE dlmod', '  IMPLICIT NONE', '  REAL*8 :: mpi2_dl, MD_dl', '  LOGICAL :: dl_ready = .false.'] + decls + \
        ['END MODULE dlmod\n',
         'SUBROUTINE dl_setup(gA,F,hA,MD,gM,gE,mpi2)', '  USE dlmod', '  IMPLICIT NONE',
         '  REAL*8, INTENT(in) :: gA,F,hA,MD,gM,gE,mpi2', '  REAL*8 iMp',
         '  COMPLEX*16 Wv(3), cA(3,3), Y(2), ab(2,5), Gc(3,3), Wc(3)', '  REAL*8 G3(3,3), dls(3), sss(5), dlv, h',
         '  INTEGER k', '  INTEGER*4 oldkey', '  INTEGER, EXTERNAL :: getversionkey',
         '  mpi2_dl=mpi2; MD_dl=MD; iMp=1.d0/(1.d0+MD)',
         '  oldkey=getversionkey(); CALL setversionkey(64)', '  dls(1)=1.d-4; dls(2)=2.d-4; dls(3)=3.d-4'] + setup + \
        ['  CALL setversionkey(oldkey)', '  dl_ready=.true.', 'END SUBROUTINE dl_setup\n'] + amps
    stub = {'af': 'CovDeltaFull', 'bl': 'CovDeltaGND'}
    for k, r in stub.items():                      # parts not generated in this run: stubs that stop
        if k not in which:
            L += [f'SUBROUTINE {r}(x,t,A)', '  IMPLICIT NONE', '  REAL*8, INTENT(in) :: x,t',
                  '  COMPLEX*16, INTENT(out) :: A(8)', f"  A=(0.d0,0.d0); STOP '{r}: not generated'", f'END SUBROUTINE {r}\n']
    L.append(BHIGH_F90)
    write_if_changed(os.path.join(OUT, 'CovDeltaLoops.F90'), '\n'.join(L))
    print('graph routines:', len(names))

BHIGH_F90 = '''
SUBROUTINE dl_bhigh(j,r,p2,m0,m1,fin,uv)
  ! two-point tensor coefficient B_{(00)^j 1^r}(p2,m0,m1) (LoopTools basis, mu=1, Delta_UV=0) from
  ! int_0^1 dx (-x)^r D^j/(2^j j!) (Delta_UV - ln D + H_j), D = x^2 p2 - x (p2+m0-m1) + m0 - i eps
  ! (../n3lo/py/bhigh.py; equal to LoopTools for ranks 0-3).  Split at the real zeros of D, cubic clustering.
  IMPLICIT NONE
  INTEGER, INTENT(in) :: j, r
  REAL*8, INTENT(in) :: p2, m0, m1
  COMPLEX*16, INTENT(out) :: fin, uv
  INTEGER, PARAMETER :: NG=48
  REAL*8, SAVE :: xg(NG), wg(NG)
  LOGICAL, SAVE :: init=.false.
  REAL*8 cuts(4), a, b, c, disc, rt, lo, hi, mid, s0, e0, u, w, x, D, jac, H, pre, fac
  COMPLEX*16 lnD
  INTEGER nc, i, k, m, sd
  IF(.NOT.init) THEN
     CALL dl_gauleg(NG,xg,wg); init=.true.
  ENDIF
  a=p2; b=-(p2+m0-m1); c=m0
  nc=2; cuts(1)=0.d0; cuts(2)=1.d0
  IF(ABS(a).GT.1.d-300) THEN
     disc=b*b-4.d0*a*c
     IF(disc.GT.0.d0) THEN
        DO k=-1,1,2
           rt=(-b+k*SQRT(disc))/(2.d0*a)
           IF(rt.GT.0.d0.AND.rt.LT.1.d0) THEN
              nc=nc+1; cuts(nc)=rt
           ENDIF
        ENDDO
     ENDIF
  ELSE IF(ABS(b).GT.0.d0) THEN
     rt=-c/b
     IF(rt.GT.0.d0.AND.rt.LT.1.d0) THEN
        nc=nc+1; cuts(nc)=rt
     ENDIF
  ENDIF
  DO i=1,nc-1                     ! sort
     DO k=i+1,nc
        IF(cuts(k).LT.cuts(i)) THEN
           x=cuts(i); cuts(i)=cuts(k); cuts(k)=x
        ENDIF
     ENDDO
  ENDDO
  H=0.d0
  DO k=1,j
     H=H+1.d0/k
  ENDDO
  fac=1.d0
  DO k=1,j
     fac=fac*2.d0*k
  ENDDO
  fin=(0.d0,0.d0); uv=(0.d0,0.d0)
  DO i=1,nc-1
     lo=cuts(i); hi=cuts(i+1); mid=0.5d0*(lo+hi)
     DO sd=1,2
        IF(sd.EQ.1) THEN
           s0=lo; e0=mid
        ELSE
           s0=mid; e0=hi
        ENDIF
        DO m=1,NG
           u=0.5d0*(xg(m)+1.d0); w=0.5d0*wg(m)
           IF(sd.EQ.1) THEN
              x=s0+(e0-s0)*u**3; jac=(e0-s0)*3.d0*u**2
           ELSE
              x=e0-(e0-s0)*(1.d0-u)**3; jac=(e0-s0)*3.d0*(1.d0-u)**2
           ENDIF
           D=x*x*p2-x*(p2+m0-m1)+m0
           lnD=LOG(ABS(D))
           IF(D.LT.0.d0) lnD=lnD-(0.d0,1.d0)*4.d0*ATAN(1.d0)
           pre=(-x)**r*D**j/fac
           fin=fin+w*jac*pre*(H-lnD)
           uv=uv+w*jac*pre
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE dl_bhigh

SUBROUTINE dl_gauleg(n,x,w)
  ! Gauss-Legendre nodes and weights on [-1,1]
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  REAL*8, INTENT(out) :: x(n), w(n)
  INTEGER i, k
  REAL*8 z, z1, p1, p2, p3, pp, pi
  pi=4.d0*ATAN(1.d0)
  DO i=1,(n+1)/2
     z=COS(pi*(i-0.25d0)/(n+0.5d0))
     DO
        p1=1.d0; p2=0.d0
        DO k=1,n
           p3=p2; p2=p1; p1=((2.d0*k-1.d0)*z*p2-(k-1.d0)*p3)/k
        ENDDO
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z; z=z1-p1/pp
        IF(ABS(z-z1).LT.1.d-15) EXIT
     ENDDO
     x(i)=-z; x(n+1-i)=z
     w(i)=2.d0/((1.d0-z*z)*pp*pp); w(n+1-i)=w(i)
  ENDDO
END SUBROUTINE dl_gauleg
'''

if __name__ == '__main__':
    main(which=tuple(sys.argv[1:]) or ('af', 'bl'))
