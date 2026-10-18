"""Non-Born proton Compton multipoles of the v2 NNLO BChPT code (work/code_v2/execs/ampl.exe) in the
main-fit configuration (dvertex=t, pimassphys=t, nmassphys=t, consts='LMP15'), for given static
polarisabilities (contact-term LECs relative to the code's reference values).

Method: the main program's partial-wave projection (&multipoles pwa=t, SUBROUTINE pwa_gaga: 16-point
Gauss-Legendre rule in cos(theta_cm), J <= 5/2 with jmax=5) of the FULL amplitude, minus the same
projection of the Born + pi0-pole amplitude (run with loops and Delta switched off; identical
quadrature, so that the subtraction is exact).  The printed multipoles (f * 1e4 * m_pi/M, nucleon-mass
units) are converted to fm.  The projection conventions of pwa_gaga coincide with those of compton.py
(HMP16 Eqs. 3.19-3.22) -- checked in multipoles_figure.py (validation block).
Results are cached in work/pwa/fortran/runs/cache_v2 (key: hash of input files + executable mtime)."""
import os
import re
import hashlib
import subprocess
import tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
WORK = os.path.dirname(PWA)
EXE = os.path.join(WORK, 'code_v2', 'execs', 'ampl.exe')
CACHE = os.path.join(PWA, 'fortran', 'runs', 'cache_v2')
os.makedirs(CACHE, exist_ok=True)

SMM, HBARC, PIONM = 0.93827, 0.19733, 0.13957          # code constants (LMP15 set)
PRINT2FM = (HBARC / SMM) / (1e4 * PIONM / SMM)          # printed value -> fm
V2OPTS = "dvertex=t, consts='LMP15', pimassphys=t, nmassphys=t"
# reference (= NNLO) static polarisabilities of the v2 amplitude with physical pion/neutron masses, as
# printed by ampl.exe ("[v2] reference polarisabilities"); spin polarisabilities are passed at these values
# so that their LECs vanish identically
# 2026-10-18: values of the executable rebuilt with LoopTools-dq (code_v2, ampl.exe of 2026-10-18 00:57); the build before
# printed alphaE1=11.41543, betaM1=3.72913, gE1E1=-3.24588 (gM1M1, gE1M2, gM1E2 unchanged)
REF = dict(alphaE1=11.41553, betaM1=3.72902, gE1E1=-3.26075, gM1M1=2.93919, gE1M2=0.17008, gM1E2=1.13257)

BARYON = ("&hadr exch=t,t,f,nq=-1,1,1,rm=0.93827,wid=0.,cutf=2.2,gstr=13.8,gem=0./\n"
          "&hadr exch=f,f,f,nq=-1,1,1,rm=1.47,wid=.2,cutf=1.7,gstr=2.6,gem=-1./\n"
          "&hadr exch={dex},{dex},f,nq=1,3,3,rm=1.232,wid=0.12,cutf=2.6,gstr=0.36,gem=-1./\n")

INDAT = (" &energy Emin={Emin:.10f},Emax={Emax:.10f},Edel={Edel:.10f},kvar=3/\n"
         " &angles thmin=90.,thmax=90.,thdel=2.,cmsangle=t/\n"
         " &incl anomaly=t,oneNloop={loops},oneDloop={loops},relchpt=t,loopoints=640/\n"
         " &Delta gaugeinv=t,dng1=2.97,dng2=-1.,dg1=4.2,dg2=4.0,dz1=0.1,dz2=2.25/\n"
         " &constants gA=1.27d0,fpi=0.09221d0,pi0mass=0.13498d0,pionm=0.13957/\n"
         " &nucleon hbchpt=f,born=t,nuclamm=t,onlyLEX=f/\n"
         " &multipoles pwa={pwa},jmax=5/\n"
         " &waytocompute feynrule=f/\n"
         " &isospin proton=t/\n"
         " &scalarPs alphaE1fm={alphaE1:.6f},betaM1fm={betaM1:.6f},varyscalarPs=f/\n"
         " &spinPs gamEEfm={gE1E1:.6f},gamEMfm={gE1M2:.6f},gamMEfm={gM1E2:.6f},gamMMfm={gM1M1:.6f}/\n"
         "{driver}"
         " &v2 {v2} /\n")

PAT = re.compile(r'^\s*([-+0-9.Ee]+)\s+f\s*(EE|MM|EM|ME)\s+(\d)\s+([+-])\s+\(([^,]+),([^)]+)\)')


def _exe_tag():
    st = os.stat(EXE)
    return f'{EXE}:{st.st_size}:{st.st_mtime}'


def _run(files, outputs):
    key = hashlib.sha1((_exe_tag() + ''.join(files[k] for k in sorted(files))).encode()).hexdigest()[:20]
    cf = os.path.join(CACHE, key + '.npz')
    if os.path.exists(cf):
        z = np.load(cf, allow_pickle=True)
        return {k: str(z[k]) for k in outputs}
    with tempfile.TemporaryDirectory() as td:
        for k, v in files.items():
            open(os.path.join(td, k), 'w').write(v)
        r = subprocess.run([EXE], cwd=td, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(r.stdout[-3000:] + r.stderr[-3000:])
        out = {}
        for k in outputs:
            p = os.path.join(td, k)
            out[k] = open(p).read() if os.path.exists(p) else ''
        out['stdout'] = r.stdout
    np.savez(cf, **{k: np.array(v) for k, v in out.items()})
    return {k: out[k] for k in outputs}


def _pols(pols):
    p = dict(REF)
    if pols:
        p.update(pols)
    return p


def multipoles(E_MeV, pols=None, born_only=False, v2=V2OPTS):
    """Real and imaginary parts of the multipoles (fm) of the FULL amplitude (or of Born+pi0 only) at
    lab photon energies E_MeV (uniform grid).  Returns dict key -> complex array, keys ('EE',1,'-') etc."""
    E = np.atleast_1d(np.asarray(E_MeV, float))
    if len(E) > 2 and not np.allclose(np.diff(E), E[1] - E[0]):        # non-uniform: one run per energy
        parts = [multipoles([e], pols, born_only, v2) for e in E]
        res = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if not str(k).startswith('_')}
        res['_ref'], res['_stdout'] = parts[0]['_ref'], parts[0]['_stdout']
        return res
    step = E[1] - E[0] if len(E) > 1 else 1.0
    files = {'baryon.dat': BARYON.format(dex='f' if born_only else 't'),
             'in.dat': INDAT.format(Emin=E[0] / 1e3, Emax=(E[-1] + 0.25 * step) / 1e3, Edel=step / 1e3,
                                    loops='f' if born_only else 't', pwa='t', driver='', v2=v2, **_pols(pols))}
    out = _run(files, ['fort.14', 'fort.15', 'fort.16', 'stdout'])
    res = {}
    for fn in ('fort.14', 'fort.15', 'fort.16'):
        for line in out[fn].splitlines():
            m = PAT.match(line)
            if not m:
                continue
            Ep = float(m.group(1))
            i = int(np.argmin(abs(E - Ep)))
            if abs(E[i] - Ep) > 1e-6 * max(1.0, Ep):
                continue
            k = (m.group(2), int(m.group(3)), m.group(4))
            res.setdefault(k, np.full(len(E), np.nan, complex))[i] = complex(float(m.group(5)), float(m.group(6))) * PRINT2FM
    for k, v in res.items():
        assert np.all(np.isfinite(v)), f'missing energies for {k}'
    ref = re.search(r'reference polarisabilities \(aE1,bM1,gE1E1,gM1M1,gE1M2,gM1E2\):\s*([-0-9.\s]+)', out['stdout'])
    res['_ref'] = np.array(ref.group(1).split(), float) if ref else None
    res['_stdout'] = out['stdout']
    return res


def nonborn(E_MeV, pols=None, v2=V2OPTS):
    """non-Born multipoles (fm, complex) = full - (Born + pi0 pole), identical quadrature."""
    full = multipoles(E_MeV, pols, v2=v2)
    born = multipoles(E_MeV, born_only=True, v2=v2)
    out = {k: full[k] - born[k] for k in full if not str(k).startswith('_')}
    out['_ref'] = full['_ref']
    return out, born


def driver(points, pols=None, v2=V2OPTS):
    """fit-driver evaluation (direct full amplitude) of observables.
    points: list of (id, E_lab, theta, OBS, frame).  Returns dict id -> value at the given polarisabilities."""
    pts = ''.join(f'{i} {E:.6f} {t:.6f} {o} {fr}\n' for i, E, t, o, fr in points)
    files = {'baryon.dat': BARYON.format(dex='t'), 'points.txt': pts,
             'in.dat': INDAT.format(Emin=0.05, Emax=0.4, Edel=0.05, loops='t', pwa='f', v2=v2,
                                    driver=" &fitdriver pointsmode=t, pointsfile='points.txt', outfile='points.out', selftest=f /\n",
                                    **_pols(pols))}
    out = _run(files, ['points.out'])
    res = {}
    for line in out['points.out'].splitlines():
        if line.startswith('#') or not line.strip():
            continue
        p = line.split()
        res[p[0]] = float(p[5])
    return res
