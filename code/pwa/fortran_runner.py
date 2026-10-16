"""Run the (copied) Fortran BChPT code ampl.exe in helicity-dump mode and return amplitudes.
Configurations:
  'bchpt'  : complete NNLO BChPT (LEC shifts zero -> NNLO predictions for the static polarisabilities)
  'lex'    : onlyLEX=t  (Born + pi0 + contact polarisability terms with the given static values)
"""
import os, subprocess, tempfile, hashlib
import numpy as np
import compton as C

BASE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(BASE, 'fortran', 'execs', 'ampl.exe')
BARYON = os.path.join(BASE, 'fortran', 'execs', 'baryon.dat')
CACHE = os.path.join(BASE, 'fortran', 'runs', 'cache')
os.makedirs(CACHE, exist_ok=True)

NNLO = dict(alphaE1fm=11.2461, betaM1fm=3.87054, gamEEfm=-3.34375, gamEMfm=0.185504,
            gamMEfm=1.14328, gamMMfm=2.95934)

INDAT = """ &energy Emin=0.05,Emax=0.4,Edel=0.05,kvar=3/
 &angles thmin=20.,thmax=160.,thdel=2.,cmsangle=f/
 &incl anomaly=t,oneNloop={loops},oneDloop={loops},relchpt=t,loopoints=640/
 &Delta gaugeinv=t,dng1=2.97,dng2=-1.,dg1=4.2,dg2=4.0,dz1=0.1,dz2=2.25/
 &constants gA=1.27d0,fpi=0.09221d0,pi0mass=0.13498d0,pionm=0.13957/
 &nucleon hbchpt=f,born=t,nuclamm=t,onlyLEX={lex}/
 &multipoles pwa=f,jmax=3/
 &waytocompute feynrule=f/
 &isospin proton=t/
 &scalarPs alphaE1fm={alphaE1fm},betaM1fm={betaM1fm},varyscalarPs=f/
 &spinPs gamEEfm={gamEEfm},gamEMfm={gamEMfm},gamMEfm={gamMEfm},gamMMfm={gamMMfm}/
 &fitdriver helmode=t, pointsfile='pts.txt', outfile='hel.out', pointsmode=f /
"""


def run(E_MeV, theta_deg, frame='CM', config='bchpt', pols=None, use_cache=True):
    """Return Phi-normalised helicity arrays per component: M[N,6,i,j,hp,h] in fm.
    components: 0 Born, 1 pi0, 2 piN loops, 3 Delta pole, 4 piDelta loops, 5 contact/LEC."""
    E = np.atleast_1d(np.asarray(E_MeV, float))
    th = np.atleast_1d(np.asarray(theta_deg, float))
    E, th = np.broadcast_arrays(E, th)
    p = dict(NNLO)
    if pols:
        p.update(pols)
    lex = (config == 'lex')
    txt = INDAT.format(loops='f' if lex else 't', lex='t' if lex else 'f', **p)
    pts = ''.join(f'x{i} {e:.10f} {t:.10f} {frame}\n' for i, (e, t) in enumerate(zip(E.ravel(), th.ravel())))
    key = hashlib.sha1((txt + pts).encode()).hexdigest()[:16]
    cf = os.path.join(CACHE, key + '.npz')
    if use_cache and os.path.exists(cf):
        z = np.load(cf)
        return z['M'].reshape(E.shape + z['M'].shape[1:]), dict(z['k'].item()) if 'k' in z else None
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, 'in.dat'), 'w').write(txt)
        open(os.path.join(td, 'pts.txt'), 'w').write(pts)
        os.symlink(BARYON, os.path.join(td, 'baryon.dat'))
        r = subprocess.run([EXE], cwd=td, capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(os.path.join(td, 'hel.out')):
            raise RuntimeError(r.stdout[-2000:] + r.stderr[-2000:])
        ids, d = C.read_fortran_hel(os.path.join(td, 'hel.out'))
    M = C.fortran_to_phi_norm(d['M'], d['s'], d['EN'])
    np.savez(cf, M=M)
    return M.reshape(E.shape + M.shape[1:]), None
