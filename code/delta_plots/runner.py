"""
runner.py -- executable-agnostic runner of the BChPT Compton code in fit-driver mode, with every
mechanism of the amplitude switchable separately (needed for the Gram-matrix decomposition in
theory_dp.py).

Mechanisms (all at the order of the code, delta counting):
   Born   nucleon s+u pole graphs with Dirac+Pauli couplings      baryon.dat line 1  exch=t,t
   pi0    t-channel pi0 (WZW) pole                                 &incl anomaly
   piN    covariant O(p^3) piN loops                               &incl oneNloop
   Dpole  Delta(1232) s+u pole (gM, gE, energy-dependent width)    baryon.dat line 3  exch=t,t
   piD    covariant piDelta loops (MS-bar)                         &incl oneDloop
The driver evaluates each configuration with all contact terms switched off and adds the O(p^4)
contact terms (shifts d = P - Pref of the six dipole polarisabilities) as an exact quadratic form.

Results are cached in work/delta_plots/cache keyed by sha1(points, in.dat, baryon.dat, executable
path + size + mtime), so switching AMPL_EXE (e.g. to work/code_v2) never mixes results.
"""
import os, re, hashlib, subprocess, tempfile, shutil, json
import numpy as np
import config as C

MECHS = ('Born', 'pi0', 'piN', 'Dpole', 'piD')
PNAMES = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
PAIRS = [(a, b) for a in range(7) for b in range(a, 7)]
PIDX = {p: i for i, p in enumerate(PAIRS)}

# identical to work/fits/lib/bchpt.py DEFAULT_TEMPLATE + a trailing &v2 group (code_v2 options; executables without it
# ignore it).  @V2OPTS@ is replaced by config.v2opts(dvertex).
DEFAULT_TEMPLATE = """ &energy Emin=0.05,Emax=0.4,Edel=0.05,kvar=3/
 &angles thmin=20.,thmax=160.,thdel=2.,cmsangle=f/
 &incl anomaly=@ANOMALY@,oneNloop=@ONENLOOP@,oneDloop=@ONEDLOOP@,relchpt=t,loopoints=640/
 &Delta gaugeinv=t,dng1=@GM@,dng2=@GE@,dg1=4.2,dg2=4.0,dz1=0.1,dz2=2.25/
 &constants gA=1.270000d0,fpi=0.0922100d0,pi0mass=0.134980,pionm=0.139570/
 &nucleon hbchpt=f,born=t,nuclamm=t,onlyLEX=f/
 &multipoles pwa=f,jmax=3/
 &waytocompute feynrule=f/
 &isospin proton=t/
 &scalarPs alphaE1fm=-999.0,betaM1fm=-999.0,varyscalarPs=f/
 &spinPs gamEEfm=-999.0,gamEMfm=-999.0,gamMEfm=-999.0,gamMMfm=-999.0/
 &fitdriver pointsmode=t, pointsfile='points.txt', outfile='points.out', selftest=f/
 &v2 order=3.5, consts='LMP15', @V2OPTS@ /
"""

BARYON = """&hadr exch=@NEX@,@NEX@,f,nq=-1,1,1,rm=0.9383,wid=0.,cutf=2.2,gstr=13.8,gem=0./
&hadr exch=f,f,f,nq=-1,1,1,rm=1.47,wid=.2,cutf=1.7,gstr=2.6,gem=-1./
&hadr exch=@DEX@,@DEX@,f,nq=1,3,3,rm=1.232,wid=0.12,cutf=2.6,gstr=0.36,gem=-1./
&hadr exch=f,f,f,nq=1,1,1,rm=1.555,wid=0.12,cutf=2.0,gstr=0.6,gem=-0.3/
&hadr exch=f,f,f,nq=-1,3,1,rm=1.525,wid=.2,cutf=1.7,gstr=0.8,gem=-0.3/
&hadr exch=f,f,f,nq=1,3,1,rm=1.820,wid=0.12,cutf=1.9,gstr=0.80,gem=-0.25/
&hadr exch=f,f,f,nq=-1,5,1,rm=1.67,wid=.2,cutf=1.4,gstr=2.2,gem=3./
&hadr exch=f,f,f,nq=-1,3,3,rm=1.700,wid=.2,cutf=0.,gstr=0.5,gem=-.25/
"""


def tf(b):
    return 't' if b else 'f'


def make_inputs(mechs=MECHS, gM=C.GM0, template=None, dvertex=True, irun=None, gE=None):
    """irun: &v2 irun of code_v2 (running Delta parameters, 0-3); None = as DP_V2OPTS (default: no running);
    gE: gamma N Delta electric coupling dng2 (default config.GE0 = -1.0; the rendered in.dat for the default is unchanged)"""
    gE = C.GE0 if gE is None else float(gE)
    ind = template if template is not None else DEFAULT_TEMPLATE
    rep = {'@ANOMALY@': tf('pi0' in mechs), '@ONENLOOP@': tf('piN' in mechs), '@ONEDLOOP@': tf('piD' in mechs),
           '@GM@': repr(float(gM)), '@GE@': repr(float(gE)), '@V2OPTS@': C.v2opts(dvertex, irun)}
    for k, v in rep.items():
        ind = ind.replace(k, v)
    bar = BARYON.replace('@NEX@', tf('Born' in mechs)).replace('@DEX@', tf('Dpole' in mechs))
    return ind, bar


def points_text(points):
    """points: iterable of (id, E_lab, theta, OBS, frame)"""
    return ''.join('%-28s %12.6f %12.6f %-6s %s\n' % (pid, float(E), float(th), obs, fr) for pid, E, th, obs, fr in points)


def run(points, mechs=MECHS, gM=C.GM0, exe=None, template=None, name='run', verbose=True, dvertex=True, irun=None, gE=None):
    exe = exe or C.EXE
    if template is None and C.TEMPLATE and exe == C.EXE:
        template = open(C.TEMPLATE).read()
    ind, bar = make_inputs(mechs, gM, template, dvertex, irun, gE)
    ptxt = points_text(points)
    st = os.stat(exe)
    key = hashlib.sha1((ptxt + ind + bar + os.path.abspath(exe) + str(st.st_size) + str(int(st.st_mtime))).encode()).hexdigest()[:16]
    os.makedirs(C.CACHE, exist_ok=True)
    fout = os.path.join(C.CACHE, '%s_%s.out' % (name, key))
    if not os.path.exists(fout):
        wd = tempfile.mkdtemp(prefix='dplt_')
        try:
            open(os.path.join(wd, 'in.dat'), 'w').write(ind)
            open(os.path.join(wd, 'baryon.dat'), 'w').write(bar)
            open(os.path.join(wd, 'points.txt'), 'w').write(ptxt)
            with open(os.path.join(wd, 'run.log'), 'w') as log:
                r = subprocess.run([exe], cwd=wd, stdout=log, stderr=subprocess.STDOUT)
            if r.returncode != 0 or not os.path.exists(os.path.join(wd, 'points.out')):
                raise RuntimeError('ampl.exe failed in %s (see run.log)' % wd)
            shutil.copy(os.path.join(wd, 'points.out'), fout + '.tmp')
            os.replace(fout + '.tmp', fout)
            json.dump(dict(mechs=list(mechs), gM=gM, gE=(C.GE0 if gE is None else float(gE)), dvertex=dvertex, irun=irun, v2opts=C.v2opts(dvertex, irun),
                           exe=os.path.abspath(exe), npoints=len(points)), open(fout + '.json', 'w'), indent=1)
        except Exception:
            print('  [runner] run directory kept for inspection:', wd)
            raise
        else:
            shutil.rmtree(wd, ignore_errors=True)
        if verbose:
            print('  [runner] %-22s %5d points -> %s' % (name, len(points), os.path.basename(fout)))
    return DriverOutput(fout)


def ffloat(x):
    """Fortran real, incl. the E-less form of 3-digit exponents (e.g. '8.42-118')"""
    try:
        return float(x)
    except ValueError:
        m = re.match(r'^([+-]?[0-9.]+)([+-][0-9]{3})$', x)
        if m:
            return float(m.group(1) + 'E' + m.group(2))
        raise


class DriverOutput:
    """Reader of the fit-driver output (format of work/code/driver/README.md Sect. 4)."""

    def __init__(self, fn):
        self.fn = fn
        self.Pref = None; self.comp = {}
        ids, obs, E, thl, thc, val, flag, cn, cd = [], [], [], [], [], [], [], [], []
        for line in open(fn):
            if line.startswith('#'):
                if line.startswith('# Pref'):
                    self.Pref = np.array([float(v) for v in line.split('=')[1].split()])
                elif ':' in line and ('piN loops' in line or 'Delta pole' in line or 'piDelta loops' in line):
                    k, v = line[1:].split(':', 1)
                    try:
                        self.comp[k.strip()] = np.array([float(x) for x in v.split()])
                    except ValueError:
                        pass
                continue
            p = line.split()
            if not p:
                continue
            ids.append(p[0]); obs.append(p[1])
            E.append(float(p[2])); thl.append(float(p[3])); thc.append(float(p[4]))
            val.append(ffloat(p[5])); flag.append(int(p[6]))
            c = np.array([ffloat(v) for v in p[7:]])
            cn.append(c[:28]); cd.append(c[28:])
        self.ids = np.array(ids); self.obs = np.array(obs)
        self.E = np.array(E); self.theta_lab = np.array(thl); self.theta_cm = np.array(thc)
        self.value_file = np.array(val); self.flag = np.array(flag)
        self.CN = np.array(cn); self.CD = np.array(cd)

    @staticmethod
    def monomials(d):
        d = np.atleast_2d(d)
        dd = np.concatenate([np.ones((d.shape[0], 1)), d], axis=1)
        return np.stack([dd[:, a] * dd[:, b] for a, b in PAIRS], axis=1)

    def predict_array(self, P):
        """observable at absolute polarisabilities P (shape (6,) or (npoints, 6))"""
        P = np.asarray(P, float)
        if P.ndim == 1:
            m = self.monomials(P - self.Pref)[0]
            return (self.CN @ m) / (self.CD @ m)
        m = self.monomials(P - self.Pref)
        return np.sum(self.CN * m, axis=1) / np.sum(self.CD * m, axis=1)
