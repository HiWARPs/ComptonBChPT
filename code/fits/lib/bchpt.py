"""
bchpt.py -- thin, executable-agnostic interface to the NNLO BChPT Compton code (ampl.exe, fit-driver mode).

* Profiles (environment variable FIT_PROFILE, default 'v2'):
    v1 : work/code/execs/ampl.exe (audited v1), no &v2 namelist;
    v2 : work/code_v2/execs/ampl.exe with the MAIN v2 settings: gamma-N-Delta vertex corrections (dvertex=t),
         physical masses in the piN loops (pimassphys=t, nmassphys=t), LMP2015 constants (consts='LMP15').
  The executable can still be overridden with AMPL_EXE.
  An alternative in.dat template can be given with AMPL_TEMPLATE (must contain the same @KEY@
  placeholders as DEFAULT_TEMPLATE below); this allows re-running everything with work/code_v2.
* run(points, cfg) writes in.dat/baryon.dat/points.txt into a temporary directory, runs the
  executable and returns a DriverOutput.  Results are cached (release: <release>/amplitudes/) keyed by a
  hash of the inputs only (points, in.dat, baryon.dat); the executable used is recorded in the .out.json
  sidecar.  Without an executable, a missing table raises an error.
* Configurations (cfg dict, all optional):
     born, nuclamm, anomaly, oneNloop, oneDloop, delta_pole (bool)
     gM, gE, gA, fpi [GeV], pi0mass [GeV], pionm [GeV]
     refs  : dict of reference polarisabilities (refalphaE1=..)  [default: analytic NNLO of the code]
     Puser : 6 polarisabilities at which the 'value' column is evaluated (default -999 = Pref)
  Named presets: CONFIGS['NNLO'], ['NLO'], ['LO'], ['NNLO_noPiD'], ...
"""
import os, hashlib, subprocess, tempfile, shutil, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FITS = os.path.dirname(HERE)
WORK = os.path.dirname(FITS)
# [release] amplitude tables: <release>/amplitudes (next to code/) when it exists, else fits/cache; override
# with the environment variable BCHPT_CACHE
_AMPL = os.path.join(os.path.dirname(WORK), 'amplitudes')
CACHE = os.environ.get('BCHPT_CACHE') or (_AMPL if os.path.isdir(_AMPL) else os.path.join(FITS, 'cache'))
PROFILE = os.environ.get('FIT_PROFILE', 'v2')
EXES = {'v1': os.path.join(WORK, 'code', 'execs', 'ampl.exe'), 'v2': os.path.join(WORK, 'code_v2', 'execs', 'ampl.exe')}
DEFAULT_EXE = EXES[PROFILE]

PNAMES = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
PAIRS = [(a, b) for a in range(7) for b in range(a, 7)]

DEFAULT_TEMPLATE = """ &energy Emin=0.05,Emax=0.4,Edel=0.05,kvar=3/
 &angles thmin=20.,thmax=160.,thdel=2.,cmsangle=f/
 &incl anomaly=@ANOMALY@,oneNloop=@ONENLOOP@,oneDloop=@ONEDLOOP@,relchpt=t,loopoints=640/
 &Delta gaugeinv=t,dng1=@GM@,dng2=@GE@,dg1=4.2,dg2=4.0,dz1=0.1,dz2=2.25/
 &constants gA=@GA@,fpi=@FPI@,pi0mass=@PI0MASS@,pionm=@PIONM@/
 &nucleon hbchpt=f,born=@BORN@,nuclamm=@NUCLAMM@,onlyLEX=f/
 &multipoles pwa=f,jmax=3/
 &waytocompute feynrule=f/
 &isospin proton=t/
 &scalarPs alphaE1fm=@ALPHA@,betaM1fm=@BETA@,varyscalarPs=f/
 &spinPs gamEEfm=@GEE@,gamEMfm=@GEM@,gamMEfm=@GME@,gamMMfm=@GMM@/
 &fitdriver pointsmode=t, pointsfile='points.txt', outfile='points.out', selftest=f @REFS@/
"""

BARYON = """&hadr exch=t,t,f,nq=-1,1,1,rm=0.9383,wid=0.,cutf=2.2,gstr=13.8,gem=0./
&hadr exch=f,f,f,nq=-1,1,1,rm=1.47,wid=.2,cutf=1.7,gstr=2.6,gem=-1./
&hadr exch=@DEX@,@DEX@,f,nq=1,3,3,rm=1.232,wid=0.12,cutf=2.6,gstr=0.36,gem=-1./
&hadr exch=f,f,f,nq=1,1,1,rm=1.555,wid=0.12,cutf=2.0,gstr=0.6,gem=-0.3/
&hadr exch=f,f,f,nq=-1,3,1,rm=1.525,wid=.2,cutf=1.7,gstr=0.8,gem=-0.3/
&hadr exch=f,f,f,nq=1,3,1,rm=1.820,wid=0.12,cutf=1.9,gstr=0.80,gem=-0.25/
&hadr exch=f,f,f,nq=-1,5,1,rm=1.67,wid=.2,cutf=1.4,gstr=2.2,gem=3./
&hadr exch=f,f,f,nq=-1,3,3,rm=1.700,wid=.2,cutf=0.,gstr=0.5,gem=-.25/
"""

BASE = dict(born=True, nuclamm=True, anomaly=True, oneNloop=True, oneDloop=True, delta_pole=True,
            gM=2.97, gE=-1.0, gA=1.27, fpi=0.09221, pi0mass=0.13498, pionm=0.13957, refs=None, Puser=None,
            # v2 namelist &v2 (None = not written)
            order=None, consts=None, dvertex=None, drunning=None, hA=None, v2x=None)
PROFILE_BASE = {'v1': {}, 'v2': dict(dvertex=True, consts='LMP15', v2x='pimassphys=t,nmassphys=t')}

CONFIGS = {
    'NNLO': {},
    # order-by-order (delta counting): LO = Born without anomalous magnetic moment (Thomson/Dirac, O(p^2))
    'LO': dict(nuclamm=False, anomaly=False, oneNloop=False, oneDloop=False, delta_pole=False),
    # NLO = O(p^3): + kappa (Pauli) Born terms + pi0 pole + piN loops
    'NLO': dict(oneDloop=False, delta_pole=False),
    # Born (with kappa) + pi0 pole only: reference amplitude of the 'LEC-matched' truncation scheme
    'BornPi0': dict(oneNloop=False, oneDloop=False, delta_pole=False),
    # NNLO without piDelta loops / without Delta pole (diagnostics)
    'NNLO_noPiD': dict(oneDloop=False),
    'NNLO_noDpole': dict(delta_pole=False),
    # parameter variations
    'gM2.87': dict(gM=2.87), 'gM3.07': dict(gM=3.07), 'gM2.90': dict(gM=2.90),
    'gA1.2754': dict(gA=1.2754), 'fpi92.07': dict(fpi=0.09207),
    'PDGconst': dict(gA=1.2754, fpi=0.09207),
}
if PROFILE == 'v2':
    CONFIGS.update({
        'LO': dict(order=2.0),                     # v2 order switch (O(p^2) Born, kappa = 0)
        'NLO': dict(order=3.0),                    # O(p^3): + AMM, pi0 pole, piN loops (physical masses)
        'PDGconst': dict(consts='PDG24'), 'PDG24': dict(consts='PDG24'),
        'hA2.84': dict(hA=2.84), 'hA2.91': dict(hA=2.91),
        'novtx': dict(dvertex=False),              # without gamma-N-Delta vertex corrections
        'singlemass': dict(v2x=''),                # one (charged) pion mass in the loops (v1 physics)
        'pi0only': dict(v2x='pimassphys=t'),
        'v1equiv': dict(dvertex=False, v2x='iosmode=1,loopsafe=f'),
        'LO_flags': dict(nuclamm=False, anomaly=False, oneNloop=False, oneDloop=False, delta_pole=False),
        'NLO_flags': dict(oneDloop=False, delta_pole=False),
    })


def tf(b):
    return 't' if b else 'f'


def exe_path():
    return os.environ.get('AMPL_EXE') or DEFAULT_EXE


def template():
    p = os.environ.get('AMPL_TEMPLATE')
    if p:
        return open(p).read()
    return DEFAULT_TEMPLATE


def make_inputs(cfg):
    c = dict(BASE)
    c.update(PROFILE_BASE[PROFILE])
    c.update(cfg)
    P = c['Puser'] if c['Puser'] is not None else [-999.0] * 6
    refs = ''
    if c['refs']:
        refs = ', ' + ', '.join('%s=%.8f' % (k, v) for k, v in c['refs'].items())
    rep = {'@ANOMALY@': tf(c['anomaly']), '@ONENLOOP@': tf(c['oneNloop']), '@ONEDLOOP@': tf(c['oneDloop']),
           '@GM@': repr(float(c['gM'])), '@GE@': repr(float(c['gE'])), '@GA@': '%.6fd0' % c['gA'],
           '@FPI@': '%.7fd0' % c['fpi'], '@PI0MASS@': '%.6f' % c['pi0mass'], '@PIONM@': '%.6f' % c['pionm'],
           '@BORN@': tf(c['born']), '@NUCLAMM@': tf(c['nuclamm']),
           '@ALPHA@': repr(float(P[0])), '@BETA@': repr(float(P[1])), '@GEE@': repr(float(P[2])),
           '@GMM@': repr(float(P[3])), '@GEM@': repr(float(P[4])), '@GME@': repr(float(P[5])), '@REFS@': refs}
    ind = template()
    for k, v in rep.items():
        ind = ind.replace(k, v)
    v2 = []
    if c['order'] is not None:
        v2.append('order=%.1f' % c['order'])
    if c['consts'] is not None:
        v2.append("consts='%s'" % c['consts'])
    if c['dvertex'] is not None:
        v2.append('dvertex=%s' % tf(c['dvertex']))
    if c['drunning'] is not None:
        v2.append('drunning=%s' % tf(c['drunning']))
    if c['hA'] is not None:
        v2.append('hA=%.4f' % c['hA'])
    if c['v2x']:
        v2.append(c['v2x'])
    if v2:
        ind = ind.rstrip('\n') + '\n &v2 ' + ', '.join(v2) + '/\n'
    bar = BARYON.replace('@DEX@', tf(c['delta_pole']))
    return ind, bar


def points_text(points):
    """points: iterable of (id, E_lab, theta, OBS, frame)"""
    lines = []
    for p in points:
        pid, E, th, obs, fr = p
        lines.append('%-24s %12.6f %12.6f %-6s %s' % (pid, float(E), float(th), obs, fr))
    return '\n'.join(lines) + '\n'


def run(points, cfg=None, name=None, exe=None, verbose=True):
    cfg = cfg or {}
    if isinstance(cfg, str):
        name = name or cfg
        if cfg not in CONFIGS and cfg.startswith('gM'):
            CONFIGS[cfg] = dict(gM=float(cfg[2:]))      # e.g. 'gM2.85'
        cfg = CONFIGS[cfg]
    exe = exe or exe_path()
    ind, bar = make_inputs(cfg)
    ptxt = points_text(points)
    # [release] the table is keyed by its inputs only (points, in.dat, baryon.dat), not by the executable, so
    # that the tables shipped in amplitudes/ are found without the Fortran code.  To recompute after rebuilding
    # the code, delete the table(s) or the directory.  The .out.json sidecar records the executable and the
    # code version (environment variable BCHPT_CODE_ID) used.
    key = hashlib.sha1((ptxt + ind + bar).encode()).hexdigest()[:16]
    os.makedirs(CACHE, exist_ok=True)
    fout = os.path.join(CACHE, '%s_%s.out' % (name or 'run', key))
    if not os.path.exists(fout):
        if not os.path.isfile(exe):
            raise RuntimeError('amplitude table %s is not in %s and no executable is available (%s): build the '
                               'Fortran code (code/code_v2, see README.md) or set AMPL_EXE'
                               % (os.path.basename(fout), CACHE, exe))
        st = os.stat(exe)
        wd = tempfile.mkdtemp(prefix='bchpt_')
        try:
            open(os.path.join(wd, 'in.dat'), 'w').write(ind)
            open(os.path.join(wd, 'baryon.dat'), 'w').write(bar)
            open(os.path.join(wd, 'points.txt'), 'w').write(ptxt)
            with open(os.path.join(wd, 'run.log'), 'w') as log:
                r = subprocess.run([exe], cwd=wd, stdout=log, stderr=subprocess.STDOUT)
            if r.returncode != 0 or not os.path.exists(os.path.join(wd, 'points.out')):
                raise RuntimeError('ampl.exe failed in %s' % wd)
            shutil.copy(os.path.join(wd, 'points.out'), fout + '.tmp')
            os.replace(fout + '.tmp', fout)
            json.dump(dict(cfg=dict(BASE, **cfg), profile=PROFILE, name=name, npoints=len(points),
                           exe=os.path.abspath(exe), exe_size=st.st_size, exe_mtime=int(st.st_mtime),
                           code_id=os.environ.get('BCHPT_CODE_ID', ''), inputs=dict(in_dat=ind, baryon_dat=bar)),
                      open(fout + '.json', 'w'), indent=1, default=str)
        finally:
            shutil.rmtree(wd, ignore_errors=True)
        if verbose:
            print('  [bchpt] computed %d points (%s) -> %s' % (len(points), name, os.path.basename(fout)))
    return DriverOutput(fout)


class DriverOutput:
    """Reader of the fit-driver output (same format as work/code/driver/compton_qf.py)."""

    def __init__(self, fn):
        self.fn = fn
        self.Pref = None; self.Puser = None; self.comp = {}
        ids, obs, E, thl, thc, val, flag, cn, cd = [], [], [], [], [], [], [], [], []
        for line in open(fn):
            if line.startswith('#'):
                if line.startswith('# Pref'):
                    self.Pref = np.array([float(v) for v in line.split('=')[1].split()])
                elif line.startswith('# Puser'):
                    self.Puser = np.array([float(v) for v in line.split('=')[1].split()])
                elif ':' in line and ('piN loops' in line or 'Delta pole' in line or 'piDelta loops' in line):
                    k, v = line[1:].split(':', 1)
                    self.comp[k.strip()] = np.array([float(x) for x in v.split()])
                continue
            p = line.split()
            if not p:
                continue
            ids.append(p[0]); obs.append(p[1])
            E.append(float(p[2])); thl.append(float(p[3])); thc.append(float(p[4]))
            val.append(float(p[5])); flag.append(int(p[6]))
            c = np.array([float(v) for v in p[7:]])
            cn.append(c[:28]); cd.append(c[28:])
        self.ids = np.array(ids); self.obs = np.array(obs)
        self.E = np.array(E); self.theta_lab = np.array(thl); self.theta_cm = np.array(thc)
        self.value_file = np.array(val); self.flag = np.array(flag)
        self.CN = np.array(cn); self.CD = np.array(cd)
        self.index = {k: i for i, k in enumerate(self.ids)}

    def subset(self, idlist):
        idx = np.array([self.index[i] for i in idlist])
        return idx

    @staticmethod
    def monomials(d):
        d = np.atleast_2d(d)
        dd = np.concatenate([np.ones((d.shape[0], 1)), d], axis=1)
        return np.stack([dd[:, a] * dd[:, b] for a, b in PAIRS], axis=1)

    def predict_array(self, P, idx=None):
        P = np.asarray(P, dtype=float)
        single = (P.ndim == 1)
        m = self.monomials(np.atleast_2d(P) - self.Pref)
        CN = self.CN if idx is None else self.CN[idx]
        CD = self.CD if idx is None else self.CD[idx]
        out = (m @ CN.T) / (m @ CD.T)
        return out[0] if single else out

    def base_static(self):
        """static polarisabilities of the base amplitude actually contained in this run
        (sum of the components switched on).  For a complete NNLO run this equals Pref."""
        return self.Pref.copy()
