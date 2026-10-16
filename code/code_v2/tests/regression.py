#!/usr/bin/env python3
"""
Regression tests for ampl.exe v2 (called by run_regression.sh).

  a  static polarisabilities: amplitude extraction vs LMP2015 Tables 3-4 and vs the analytic
     App. C formulas; physical-mass and running options
  b  Thomson limit / LET: Xi- and mu-independence of A1(0,0) of the piDelta loops (bug B1),
     full NNLO vs Born+pi0+LEX (static NNLO polarisabilities) at 0.5-20 MeV
  c  unitarity: piN-loop multipoles vs tree-level photoproduction (isospin-symmetric and per
     channel with physical masses), forward Im A1 / Im A3 vs sigma_tot, sigma_1/2-sigma_3/2,
     gamma N Delta vertex Im parts (Fortran vs 60-digit PV2006, Watson at W = M_Delta)
  d  v1 vs v2 (driver grid below and above threshold); intentional differences listed
  e  gfortran vs ifort (v2)
  f  driver exactness (selftest, quadratic forms vs main program)
  g  pi0 gamma gamma transition form factor (&v2 pi0ff): pi0ff=f identical to the default, pi0ff=t leaves the
     forward points, the flags and the static polarisabilities of every mechanism unchanged, closed form vs
     LoopTools C0 (tools/check_pi0ff.F90)

usage: regression.py [a b c ...] (default: all)
Each check prints PASS / FAIL / INFO lines; the exit code is the number of FAILs.
"""
import os, sys, subprocess, tempfile, shutil, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
EXE = os.path.join(ROOT, 'execs', 'ampl.exe')
EXE_IF = os.path.join(ROOT, 'execs', 'ampl_ifort.exe')
EXE_V1 = os.path.join(HERE, 'ref_bin', 'ampl_v1.exe')
DRV = os.path.join(ROOT, 'driver', 'run_points.sh')
OUT = os.path.join(HERE, 'out')
os.makedirs(os.path.join(OUT, 'tmp'), exist_ok=True)
os.environ['V2SCRATCH'] = os.path.join(OUT, 'tmp')
sys.path.insert(0, os.path.join(ROOT, 'tools'))
sys.path.insert(0, os.path.join(ROOT, 'driver'))
NFAIL = 0


def report(ok, msg):
    global NFAIL
    tag = 'PASS' if ok is True else ('FAIL' if ok is False else 'INFO')
    if ok is False:
        NFAIL += 1
    print('  [%s] %s' % (tag, msg), flush=True)


def template(extra_v2='', fitdriver=None, onlylex='f', P=None):
    t = open(os.path.join(ROOT, 'driver', 'in_driver.template')).read()
    P = P or ['-999'] * 6
    t = (t.replace('@ALPHA@', P[0]).replace('@BETA@', P[1]).replace('@GEE@', P[2]).replace('@GMM@', P[3])
          .replace('@GEM@', P[4]).replace('@GME@', P[5]).replace('@ONLYLEX@', onlylex)
          .replace('@SELFTEST@', 'f').replace('@REFS@', '').replace('@ORDER@', '3.5').replace('@CONSTS@', 'LMP15')
          .replace('@DVERTEX@', '').replace('@DRUNNING@', '').replace('@V2X@', (', ' + extra_v2) if extra_v2 else ''))
    lines = [l for l in t.split('\n') if '&fitdriver' not in l]
    if fitdriver:
        lines.append(' &fitdriver ' + fitdriver + ' /')
    return '\n'.join(lines) + '\n'


def run_exe(indat, exe=EXE, files=()):
    d = tempfile.mkdtemp(prefix='reg_', dir=os.path.join(OUT, 'tmp'))
    shutil.copy(os.path.join(ROOT, 'driver', 'baryon.dat'), d)
    open(os.path.join(d, 'in.dat'), 'w').write(indat)
    for name, content in files:
        open(os.path.join(d, name), 'w').write(content)
    r = subprocess.run([exe], cwd=d, capture_output=True, text=True)
    return d, r.stdout


def driver(points, *keys, exe=EXE):
    d = tempfile.mkdtemp(prefix='drv_', dir=os.path.join(OUT, 'tmp'))
    p = os.path.join(d, 'pts.txt'); o = os.path.join(d, 'pts.out')
    open(p, 'w').write(points)
    subprocess.run([DRV, p, o, 'exe=' + exe] + list(keys), check=True, capture_output=True)
    res = {}
    for l in open(o):
        if l.startswith('#'):
            continue
        q = l.split()
        res[q[0]] = (float(q[5]), int(q[6]))
    shutil.rmtree(d, ignore_errors=True)
    return res


# ------------------------------------------------------------------------------------------- a
LMP_T = {'piN loops': [6.9, -1.8, -3.4, -0.1, 0.5, 0.9], 'Delta pole': [-0.1, 7.1, -0.4, 3.3, -0.4, 0.4],
         'piDelta loops': [4.4, -1.4, 0.4, -0.2, 0.1, -0.2], 'TOTAL': [11.2, 3.9, -3.3, 2.9, 0.2, 1.1]}
NAMES = ['aE1', 'bM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']


def test_a():
    print('(a) static polarisabilities')
    import static_pols
    d, log = run_exe(template(fitdriver='staticmode=t'))
    ana = {}
    for line in log.split('\n'):
        for nm, key in [('piN loops', '   piN loops '), ('Delta pole', '   Delta pole '), ('piDelta loops', '   piDelta loops'),
                        ('TOTAL', '   TOTAL ')]:
            if line.startswith(key) and 'numeric' not in line:
                ana[nm] = np.array([float(v) for v in line.split()[-6:]])
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res, tot = static_pols.main(os.path.join(d, 'static_As.out'))
    amp = {nm: res[nm][0] for nm in ['piN loops', 'Delta pole', 'piDelta loops']}
    amp['TOTAL'] = amp['piN loops'] + amp['Delta pole'] + amp['piDelta loops']
    for nm in ['piN loops', 'Delta pole', 'piDelta loops', 'TOTAL']:
        dl = np.max(np.abs(amp[nm] - np.array(LMP_T[nm])))
        da = np.max(np.abs(amp[nm] - ana[nm]))
        tol_a = 0.006 if nm in ('piN loops', 'TOTAL') else 0.002
        print('      %-14s amplitude %s' % (nm, ' '.join('%8.4f' % v for v in amp[nm])))
        report(bool(dl <= 0.06), '%-14s vs LMP2015 Tables 3-4 (rounded): max |diff| = %.3f (tol 0.06)' % (nm, dl))
        report(bool(da <= tol_a), '%-14s vs analytic App. C (current constants): max |diff| = %.4f (tol %.3f)' % (nm, da, tol_a))
    g0 = -amp['TOTAL'][2:].sum(); gpi = amp['TOTAL'][3] + amp['TOTAL'][5] - amp['TOTAL'][2] - amp['TOTAL'][4]
    report(None, 'total from amplitude: gamma0 = %.3f (LMP -0.9), gamma_pi = %.3f (LMP 7.2)' % (g0, gpi))
    # numerical extraction routine used for the physical-mass and running options
    for line in log.split('\n'):
        if 'piN numeric (both channels' in line:
            num = np.array([float(v) for v in line.split()[-6:]])
            report(bool(np.max(np.abs(num - ana['piN loops'])) < 0.006),
                   'piN_static_numeric vs analytic: max |diff| = %.4f' % np.max(np.abs(num - ana['piN loops'])))
        if 'Delta pole numeric' in line:
            num = np.array([float(v) for v in line.split()[-6:]])
            report(bool(np.max(np.abs(num - ana['Delta pole'])) < 1e-4),
                   'Delta-pole numeric extraction vs analytic: max |diff| = %.1e' % np.max(np.abs(num - ana['Delta pole'])))
    shutil.rmtree(d, ignore_errors=True)
    # option effects on the reference polarisabilities
    for opt in ['pimassphys=t', 'pimassphys=t,nmassphys=t', 'drunning=t']:
        d, log = run_exe(template(extra_v2=opt, fitdriver="pointsmode=t, pointsfile='p.txt', outfile='p.out'"),
                         files=[('p.txt', 'x 100 90 CS\n')])
        pr = [l for l in log.split('\n') if 'reference polarisabilities' in l][0]
        report(None, '%-26s Pref = %s' % (opt, pr.split(':')[-1].strip()))
        shutil.rmtree(d, ignore_errors=True)


# ------------------------------------------------------------------------------------------- b
def test_b():
    print('(b) Thomson limit / low-energy theorem')
    d, log = run_exe(template(fitdriver='staticmode=t'))
    c0 = [float(l.split('c0=')[1]) for l in log.split('\n') if 'diag_thomson: mudim' in l]
    spread = max(c0) - min(c0)
    report(bool(spread < 1e-9), 'piDelta loops A1(0,0) independent of Xi and mu (bug B1): spread %.1e over %d settings '
           '(v1: 4.9e-5); value %.1e (2-point extrapolation, LoopTools noise)' % (spread, len(c0), c0[0]))
    shutil.rmtree(d, ignore_errors=True)
    # full NNLO vs Born+pi0+contact terms carrying the analytic NNLO polarisabilities
    pts = ''.join('cs_%g_%g %g %g CS\n' % (E, th, E, th) for E in [0.5, 2, 5, 10, 12, 15, 17, 20] for th in [30, 60, 90, 150])
    full = driver(pts)
    Pref = [11.246404, 3.870657, -3.335106, 2.942725, 0.178734, 1.146936]
    lex = driver(pts, 'onlylex=t', 'alpha=%g' % Pref[0], 'beta=%g' % Pref[1], 'gEE=%g' % Pref[2], 'gMM=%g' % Pref[3],
                 'gEM=%g' % Pref[4], 'gME=%g' % Pref[5])
    for E in [0.5, 2, 5, 10, 12, 15, 17, 20]:
        dev = max(abs(full['cs_%g_%g' % (E, th)][0] / lex['cs_%g_%g' % (E, th)][0] - 1) for th in [30, 60, 90, 150])
        tol = 3e-4 * (E / 20.0) ** 4 + 1e-9
        report(bool(dev < tol), 'E=%5.1f MeV  |sigma(NNLO)/sigma(Born+pi0+LEX) - 1| max = %.1e  (tol %.1e = 3e-4 (E/20 MeV)^4: '
               'O(omega^4) terms of the Delta pole and loops)' % (E, dev, tol))
    # Thomson: E -> 0
    th = driver('t 0.05 90 CS\nt2 0.05 0 CS\n', 'order=2')
    # Thomson at 90 deg: (alpha/M)^2 (1+cos^2)/2 in nb/sr: r_p^2/2 = 0.5*(alpha*hbarc/M)^2
    rp = (1 / 137.036) * 0.19733 / 0.93827   # fm
    th90 = 0.5 * rp ** 2 * 1e7               # nb/sr  (1 fm^2 = 1e7 nb)
    report(bool(abs(th['t'][0] / th90 - 1) < 2e-4), 'Born (order 2) at 0.05 MeV, 90 deg: %.6f nb/sr vs Thomson %.6f' % (th['t'][0], th90))


# ------------------------------------------------------------------------------------------- c
def born_sigmas(W, mpi=0.13957, lmax=8):
    """sigma_tot and sigma_1/2 - sigma_3/2 (GeV^-2) of gamma p -> pi N (LO Born, both channels)"""
    from born_multipoles import calibrated
    B = calibrated(kp=0.0, kn=0.0, mpi=mpi)
    M = 0.93827
    s = W * W
    k = (s - M * M) / (2 * W)
    Eq = (s + mpi ** 2 - M * M) / (2 * W)
    q = np.sqrt(Eq * Eq - mpi * mpi)
    s12 = s32 = 0.0
    for ch in ('pi0p', 'pi+n'):
        m = B.multipoles(W, ch, lmax=lmax, ngauss=64)
        g = lambda key: m.get(key, 0.0)
        for l in range(0, lmax):
            Ep, Mp = g('E%d+' % l), g('M%d+' % l)
            Em, Mm = g('E%d-' % (l + 1)), g('M%d-' % (l + 1))
            A12p = -0.5 * ((l + 2) * Ep + l * Mp)
            A32p = 0.5 * np.sqrt(l * (l + 2)) * (Ep - Mp)
            A12m = -0.5 * (l * Em - (l + 2) * Mm)
            A32m = -0.5 * np.sqrt(l * (l + 2)) * (Em + Mm)
            s12 += (l + 1) * (abs(A12p) ** 2 + abs(A12m) ** 2)
            s32 += (l + 1) * (abs(A32p) ** 2 + abs(A32m) ** 2)
    f = 8 * np.pi * q / k
    return 0.5 * f * (s12 + s32), f * (s12 - s32)


def test_c():
    print('(c) unitarity')
    import check_unitarity_piN, check_unitarity_channels
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        w1 = check_unitarity_piN.main()
    report(bool(w1 < 1e-5), 'piN loops (single pion mass): Im f_XY^{l+-} = k sum|M_Born|^2 for 8 multipoles, '
           '160-350 MeV: max rel. dev. %.1e' % w1)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        w2 = check_unitarity_channels.main(Es=(147.0, 155.0, 200.0, 300.0))
    lines = [l for l in buf.getvalue().split('\n') if l.startswith(('pi0', 'pi+', 'sum'))]
    big = max(abs(float(x)) for l in lines for x in re.findall(r' [A-Z]{2}\d[+-] ([+-]\d\.\de[+-]\d+)', l))
    dom = max(abs(float(x)) for l in lines for x in re.findall(r' (?:MM1\+|EE1-|ME1\+) ([+-]\d\.\de[+-]\d+)', l))
    report(bool(dom < 5e-4), 'physical masses, per channel (pi0 p / pi+ n / sum) incl. 147 MeV (only pi0 p open): '
           'dominant multipoles max rel. dev. %.1e (all multipoles %.1e; small high-l ones near threshold)' % (dom, big))
    # forward Im A1, Im A3 at t = 0 (now evaluated exactly forward)
    M = 0.93827
    Es = [180.0, 250.0, 350.0]
    xt = ''.join('%.10f 0.0\n' % (E / 1000 / M) for E in Es)
    d, log = run_exe(template(fitdriver="ampsfile='xt.txt'"), files=[('xt.txt', xt)])
    D = np.loadtxt(os.path.join(d, 'amps.out'))
    shutil.rmtree(d, ignore_errors=True)
    aem = 1 / 137.036
    for row, E in zip(D, Es):
        x = row[0]
        ImA1 = row[2 + 16 * 2 + 8]; ImA3 = row[2 + 16 * 2 + 8 + 2]
        W = np.sqrt(M * M + 2 * M * E / 1000)
        st, ds = born_sigmas(W)
        st *= M * M; ds *= M * M                  # -> nucleon-mass units
        r1 = ImA1 / (x * st / (4 * np.pi * aem)); r3 = ImA3 / (x * ds / (8 * np.pi * aem))
        report(bool(abs(r1 - 1) < 2e-3 and abs(r3 - 1) < 2e-3),
               'E=%3.0f MeV, theta=0 exactly: Im A1/(nu sigma_tot/4pi alpha) = %.5f, Im A3/(nu dsigma/8pi alpha) = %.5f'
               % (E, r1, r3))
    # vertex corrections
    src = os.path.join(OUT, 'tmp', 'tv.f90')
    open(src, 'w').write('''program tv
  implicit none
  REAL*8 gA,hA,fpi,pi0mass,pionm
  COMMON/Consts/hA,gA,fpi,pi0mass,pionm
  INTEGER ivxmode
  COMMON /v2vx/ ivxmode
  REAL*8 ImgM, ImgE, Md, Ws(5)
  INTEGER i
  Ws=(/1.10d0,1.20d0,1.232d0,1.30d0,2.0d0/)
  hA=2.85d0; gA=1.27d0; fpi=0.09221d0/0.93827d0; pionm=0.13957d0/0.93827d0; pi0mass=0.13498d0/0.93827d0
  Md=1.232d0/0.93827d0; ivxmode=0
  do i=1,5
     call gnd_vertex_im(Ws(i)/0.93827d0,Md,ImgM,ImgE)
     write(*,'(3es24.15)') Ws(i), ImgM, ImgE
  enddo
end program
''')
    exe = os.path.join(OUT, 'tmp', 'tv')
    env = dict(os.environ, SDKROOT=subprocess.run(['xcrun', '--sdk', 'macosx', '--show-sdk-path'],
                                                  capture_output=True, text=True).stdout.strip())
    subprocess.run(['/opt/homebrew/bin/gfortran-16', '-O2', '-ffree-line-length-none', '-std=legacy', src,
                    os.path.join(ROOT, 'build-gfortran-native', 'As.o'), '-o', exe], check=True, env=env)
    F = np.loadtxt(subprocess.run([exe], capture_output=True, text=True).stdout.splitlines())
    import pv2006_vertex as PV
    dev = 0.0
    for W, gm, ge in F:
        a = PV.im_gME(W)
        dev = max(dev, abs(gm - a[0]) / abs(a[0]), abs(ge - a[1]) / abs(a[1]))
    report(bool(dev < 1e-9), 'gamma N Delta vertex Im gM, Im gE: Fortran (quad precision) vs 60-digit PV2006 at 5 W: '
           'max rel. dev. %.1e' % dev)
    # Watson at the peak
    from delta_pole_pw import F_quadratic, edem, MD, SMM
    from born_multipoles import calibrated
    B = calibrated(kp=0.0, kn=0.0)
    E = round((1.232 ** 2 - SMM ** 2) / (2 * SMM) * 1000, 4)
    R = F_quadratic([E])
    W = 1.232; m = 0.13957
    Eq = (W * W + m * m - SMM ** 2) / (2 * W); k = np.sqrt(Eq * Eq - m * m)
    m32, _, _ = B.mult32(W)
    MB = np.sqrt(2 / 3) * m32['M1+'].real; EB = np.sqrt(2 / 3) * m32['E1+'].real
    a = R[E]['MM'][0].real
    MG = MD * edem(W / SMM)
    kap = (W / SMM - 1) / (2 * (W / SMM + 1))
    ImG = MB * np.sqrt(k * MG / a); ImgE = EB / kap * np.sqrt(k * MG / a); ImgM = ImG + kap * ImgE
    pv = PV.im_gME(W)
    report(bool(abs(ImgM / pv[0] - 1) < 5e-3 and abs(ImgE / pv[1] - 1) < 5e-3),
           'Watson/unitarity construction at W=M_Delta (Born P33 multipoles x code Delta-pole partial wave): '
           'Im gM %.4f vs PV2006 %.4f, Im gE %.4f vs %.4f' % (ImgM, pv[0], ImgE, pv[1]))
    # Fortran Watson mode (vxmode=2, default) vs the Python construction off the peak
    src2 = os.path.join(OUT, 'tmp', 'tw.f90')
    open(src2, 'w').write(open(src).read().replace('ivxmode=0', 'ivxmode=2; aem_c=1.d0/137.036d0')
                          .replace('  COMMON /v2vx/ ivxmode\n', '  COMMON /v2vx/ ivxmode\n  REAL*8 aem_c\n  COMMON /v2aem/ aem_c\n'))
    exe2 = os.path.join(OUT, 'tmp', 'tw')
    subprocess.run(['/opt/homebrew/bin/gfortran-16', '-O2', '-ffree-line-length-none', '-std=legacy', src2,
                    os.path.join(ROOT, 'build-gfortran-native', 'As.o'), '-o', exe2], check=True, env=env)
    F2 = np.loadtxt(subprocess.run([exe2], capture_output=True, text=True).stdout.splitlines())
    Ws = [1.10, 1.20, 1.30]
    Es = [round((w * w - SMM ** 2) / (2 * SMM) * 1000, 4) for w in Ws]
    R = F_quadratic(Es)
    dev = 0.0
    for w, E in zip(Ws, Es):
        Eq = (w * w + m * m - SMM ** 2) / (2 * w); k = np.sqrt(Eq * Eq - m * m)
        m32, _, _ = B.mult32(w)
        MB = np.sqrt(2 / 3) * m32['M1+'].real; EB = np.sqrt(2 / 3) * m32['E1+'].real
        a = R[E]['MM'][0].real; MG = MD * edem(w / SMM); kap = (w / SMM - 1) / (2 * (w / SMM + 1))
        sq = np.sqrt(k * MG / a); ge = EB / kap * sq; gm = MB * sq + kap * ge
        row = [r for r in F2 if abs(r[0] - w) < 1e-9][0]
        dev = max(dev, abs(row[1] / gm - 1), abs(row[2] / ge - 1))
    report(bool(dev < 1e-4), 'vxmode=2 (Watson, default): Fortran (born_p33 + closed-form a(W)) vs Python construction from the '
           'code partial waves at W=1.1,1.2,1.3 GeV: max rel. dev. %.1e' % dev)
    # P33 unitarity of the complete amplitude with the vertex corrections (default settings)
    from pwa_run import run_pwa
    from delta_pole_pw import Dfun
    Es = [200.0, 250.0, 300.0, 330.0, 360.0]
    full = run_pwa(Es)
    piN = run_pwa(Es, oneD=False, anomaly=False, delta=(False, False))
    dev = 0.0
    for E in Es:
        w = np.sqrt(SMM ** 2 + 2 * SMM * E / 1000); Wm = w / SMM
        Eq = (w * w + m * m - SMM ** 2) / (2 * w); k = np.sqrt(Eq * Eq - m * m)
        m32, _, _ = B.mult32(w); MB = np.sqrt(2 / 3) * m32['M1+'].real
        om = (Wm * Wm - 1) / (2 * Wm)
        a = (1 / 137.036) * om ** 2 * (Wm + 1) ** 2 * (Wm + MD) / (8 * (1 + MD) ** 2) / SMM
        kap = (Wm - 1) / (2 * (Wm + 1)); MG = MD * edem(Wm)
        MDel = -np.sqrt(a * MG / k) * (2.97 + kap * 1.0 + 1j * MB * np.sqrt(k * MG / a)) / Dfun(Wm)
        pred = k * (abs(MB + MDel) ** 2 - abs(MB) ** 2)
        got = full[E]['MM1+'].imag - piN[E]['MM1+'].imag
        dev = max(dev, abs(got / pred - 1))
    report(bool(dev < 1e-5), 'complete amplitude (default, vertex corrections on): Im f_MM^{1+} = k(|M_B,1/2|^2 + '
           '|M_B,3/2 + M_Delta|^2) at 200-360 MeV (Watson theorem in P33): max rel. dev. %.1e' % dev)


# ------------------------------------------------------------------------------------------- d
def grid(Es, ths, obs=('CS', 'S3', 'S2X', 'S2Z'), frame='CM'):
    return ''.join('%s_%g_%g %g %g %s %s\n' % (o, E, th, E, th, o, frame) for E in Es for th in ths for o in obs)


def cmp(A, B, flags=None):
    out = []
    for k in A:
        a, b = A[k][0], B[k][0]
        if k.startswith('CS'):
            dv = abs(b / a - 1)
        else:
            dv = abs(b - a)
        out.append((dv, k, A[k][1]))
    return out


def test_d():
    print('(d) v1 vs v2 (CS: relative, asymmetries: absolute differences)')
    low = grid([20, 40, 60, 80, 100, 120, 140, 145], [0, 10, 30, 60, 90, 120, 150, 170, 180])
    high = grid([160, 200, 250, 300, 350], [10, 30, 60, 90, 120, 150, 170])
    v1l = driver(low, exe=EXE_V1); v1h = driver(high, exe=EXE_V1)
    # v2 with the v1 one-sided subtraction and raw loops: should equal v1 up to LoopTools noise
    v2l = driver(low, 'v2=iosmode=1,loopsafe=f'); v2h = driver(high, 'dvertex=f', 'v2=iosmode=1,loopsafe=f')
    v2d = driver(low)
    # regular points: outside the Gram-singular zones of both versions (v1 flag 0 and v2 flag 0)
    reg = [x for x in cmp(v1l, v2l) if v1l[x[1]][1] == 0 and v2d[x[1]][1] == 0]
    mx = max(reg)
    report(bool(mx[0] < 1e-6), 'below threshold, v2 with v1 subtraction (iosmode=1, loopsafe=f) vs v1, %d regular points: '
           'max %.1e (%s)' % (len(reg), mx[0], mx[1]))
    mxh = max(cmp(v1h, v2h))
    report(bool(mxh[0] < 1e-6), 'above threshold, dvertex=f, iosmode=1, loopsafe=f vs v1: max %.1e (%s)' % mxh[:2])
    # intentional changes
    c = cmp(v1l, v2d)
    regp = [x for x in c if v1l[x[1]][1] == 0 and v2d[x[1]][1] == 0]
    report(None, 'v2 default vs v1, below threshold, regular points (B1 on-shell subtraction fix): max %.1e (%s)' % max(regp)[:2])
    fb = [x for x in c if v2d[x[1]][1] != 0 or v1l[x[1]][1] != 0]
    report(None, 'near 0/180 deg (v1: observable extrapolation, v2: amplitude extrapolation): max %.1e (%s)' % max(fb)[:2])
    v2v = driver(high)
    report(None, 'above threshold, v2 default (gamma N Delta vertex corrections on) vs v1: max %.2e (%s)' % max(cmp(v1h, v2v))[:2])
    # order switch at 3.5 = default
    o35 = driver(low, 'order=3.5')
    report(bool(max(cmp(v2d, o35))[0] == 0.0), 'order=3.5 identical to the default')


# ------------------------------------------------------------------------------------------- e
def test_e():
    print('(e) gfortran vs ifort (v2)')
    if not os.path.exists(EXE_IF):
        report(False, 'ifort executable missing (make -f makefile.ifort)')
        return
    pts = grid([20, 40, 60, 100, 150, 250, 350], [0, 20, 60, 90, 120, 160, 180])
    A = driver(pts); B = driver(pts, exe=EXE_IF)
    for E in [20, 40, 60, 100, 150, 250, 350]:
        sub = [x for x in cmp(A, B) if '_%g_' % E in x[1]]
        r = max([x for x in sub if A[x[1]][1] == 0] or [(0.0, '-', 0)])
        z = max([x for x in sub if A[x[1]][1] != 0] or [(0.0, '-', 0)])
        tol = 1e-4 if E <= 20 else (1e-5 if E <= 40 else 1e-6)
        report(bool(r[0] < tol and z[0] < 10 * tol), 'E=%3d MeV: regular points max %.1e (%s) [tol %.0e]; '
               '0/180-deg zone (extrapolated) max %.1e (%s) [tol %.0e]' % (E, r[0], r[1], tol, z[0], z[1], 10 * tol))


# ------------------------------------------------------------------------------------------- f
def test_f():
    print('(f) fit driver exactness')
    d = tempfile.mkdtemp(prefix='self_', dir=os.path.join(OUT, 'tmp'))
    r = subprocess.run([DRV, os.path.join(ROOT, 'driver', 'test_points.txt'), os.path.join(d, 'o.out'), 'selftest=t'],
                       capture_output=True, text=True)
    v = [float(l.split('=')[-1]) for l in r.stdout.split('\n') if 'quadratic form vs direct' in l]
    report(bool(v and v[0] < 1e-6), 'selftest: quadratic form vs direct full-amplitude evaluation: %.1e' % (v[0] if v else -1))
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'driver', 'test_quadratic_exactness.py')], capture_output=True,
                       text=True)
    line = [l for l in r.stdout.split('\n') if l.startswith('RESULT')]
    report(bool(line and 'PASS' in line[0]), 'driver vs main program (random/huge shifts): %s' % (line[0] if line else r.stderr[-200:]))
    shutil.rmtree(d, ignore_errors=True)


# ------------------------------------------------------------------------------------------- g
def test_g():
    print('(g) pi0 gamma gamma transition form factor on the pi0-pole vertex (&v2 pi0ff)')
    import static_pols, io, contextlib
    pts = grid([20, 60, 100, 140, 150, 250, 300], [0, 30, 90, 150, 180])
    a = driver(pts); b = driver(pts, 'v2=pi0ff=f')
    report(bool(max(cmp(a, b))[0] == 0.0 and all(a[k][1] == b[k][1] for k in a)),
           'pi0ff=f identical to the default (values and flags, %d points)' % len(a))
    c = driver(pts, 'v2=pi0ff=t')
    fwd = max(x for x in cmp(a, c) if x[1].endswith('_0'))
    report(bool(fwd[0] == 0.0 and all(a[k][1] == c[k][1] for k in a)),
           'pi0ff=t: forward points (t = 0, F = 1) identical, loop-evaluation flags identical')
    bk = max(x for x in cmp(a, c) if x[1].startswith('CS') and x[1].endswith('_180'))
    report(None, 'pi0ff=t: largest relative change of dsigma at 180 deg: %.1e (%s)' % bk[:2])
    # static polarisabilities per mechanism (amplitude extraction) unchanged
    res = {}
    for ff in ('f', 't'):
        d, log = run_exe(template(extra_v2='pi0ff=' + ff, fitdriver='staticmode=t'))
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            r, _ = static_pols.main(os.path.join(d, 'static_As.out'))
        res[ff] = {nm: np.array(r[nm][0]) for nm in r}
        shutil.rmtree(d, ignore_errors=True)
    dmax = max(float(np.max(np.abs(res['t'][nm] - res['f'][nm]))) for nm in res['f'])
    an = res['t']['anomaly']
    report(bool(dmax < 1e-5), 'static polarisabilities of all six components unchanged with pi0ff=t: max |diff| = %.1e; '
           'pi0 pole: gE1E1 = -gM1M1 = gE1M2 = -gM1E2 = %.3f, gamma_pi(pi0) = %.2f' % (dmax, an[2], an[3] + an[5] - an[2] - an[4]))
    # closed form vs LoopTools C0
    LT = os.path.expanduser('~/.local/looptools-gfortran-native')
    exe = os.path.join(OUT, 'tmp', 'check_pi0ff')
    env = dict(os.environ, SDKROOT=subprocess.run(['xcrun', '--sdk', 'macosx', '--show-sdk-path'],
                                                  capture_output=True, text=True).stdout.strip())
    subprocess.run(['/opt/homebrew/bin/gfortran-16', '-O2', '-cpp', '-ffree-form', '-ffree-line-length-none', '-std=legacy',
                    '-I' + os.path.join(LT, 'include'), os.path.join(ROOT, 'tools', 'check_pi0ff.F90'),
                    os.path.join(LT, 'lib', 'libooptools.a'), '-o', exe], check=True, env=env)
    out = subprocess.run([exe], capture_output=True, text=True).stdout.splitlines()
    i0 = [i for i, l in enumerate(out) if l.startswith('#  t [GeV^2]')][0]
    i1 = [i for i, l in enumerate(out) if l.startswith('# backward')][0]
    d1 = [abs(float(l.split()[1]) - float(l.split()[2])) for l in out[i0 + 1:i1] if l.strip()]
    d2 = [abs(float(l.split()[2]) - float(l.split()[3])) for l in out[i1 + 2:] if l.strip() and not l.startswith('#')]
    f300 = [float(l.split()[4]) for l in out[i1 + 2:] if l.strip() and l.split()[0].startswith('300')]
    report(bool(d1 and d2 and max(d1 + d2) < 1e-12), 'closed-form F(t) vs -2 M^2 C0(0,0,t;M^2,M^2,M^2) (LoopTools) at %d t values '
           '(-2 ... +1 GeV^2 and t_min of 50-350 MeV): max |diff| = %.1e; F-1 = %.4f at 300 MeV/180 deg'
           % (len(d1) + len(d2), max(d1 + d2), f300[0] if f300 else float('nan')))


if __name__ == '__main__':
    which = sys.argv[1:] or list('abcdefg')
    for w in which:
        globals()['test_' + w]()
    print('REGRESSION SUMMARY: %d FAIL(s)' % NFAIL)
    sys.exit(NFAIL)
