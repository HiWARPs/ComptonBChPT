#!/usr/bin/env python3
"""
Run ampl.exe in multipole mode (&multipoles pwa=t) and parse the Compton multipoles
f_EE^{l+-}, f_MM^{l+-}, f_EM^{1+}, f_ME^{1+} (c.m. frame) that pwa_gaga writes to fort.14/15/16.

The code prints f * 1e4 * m_pi/M (i.e. in units of 1e-4/m_pi); this module returns them in GeV^-1.
Usage (module):  run_pwa(energies_MeV, flags=dict(...), baryon=dict(...), v2=dict(...), exe=...)
returns {E_lab: {'EE1-':c, 'MM1-':c, 'EE2-':c, 'MM2-':c, 'EE1+':c, 'MM1+':c, 'EM1+':c, 'ME1+':c}}
"""
import os, re, subprocess, tempfile, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, '..', 'execs', 'ampl.exe')
SMM = 0.93827


def baryon_dat(nucl=(True, True), delta=(True, True), mdelta=1.232):
    t = lambda b: 't' if b else 'f'
    return ("&hadr exch=%s,%s,f,nq=-1,1,1,rm=0.93827,wid=0.,cutf=2.2,gstr=13.8,gem=0./\n"
            "&hadr exch=f,f,f,nq=-1,1,1,rm=1.47,wid=.2,cutf=1.7,gstr=2.6,gem=-1./\n"
            "&hadr exch=%s,%s,f,nq=1,3,3,rm=%.6f,wid=0.12,cutf=2.6,gstr=0.36,gem=-1./\n"
            % (t(nucl[0]), t(nucl[1]), t(delta[0]), t(delta[1]), mdelta))


def in_dat(Emin, Emax, Edel, anomaly=True, oneN=True, oneD=True, nuclamm=True, gM=2.97, gE=-1.0,
           v2=None, theta=90.0):
    t = lambda b: 't' if b else 'f'
    s = (" &energy Emin=%.8f,Emax=%.8f,Edel=%.8f,kvar=3/\n"
         " &angles thmin=%.3f,thmax=%.3f,thdel=2.,cmsangle=t/\n"
         " &incl anomaly=%s,oneNloop=%s,oneDloop=%s,relchpt=t,loopoints=640/\n"
         " &Delta gaugeinv=t,dng1=%.8f,dng2=%.8f,dg1=4.2,dg2=4.0,dz1=0.1,dz2=2.25/\n"
         " &constants gA=1.27d0,fpi=0.09221d0,pi0mass=0.13498d0,pionm=0.13957/\n"
         " &nucleon hbchpt=f,born=t,nuclamm=%s,onlyLEX=f/\n"
         " &multipoles pwa=t,jmax=3/\n &waytocompute feynrule=f/\n &isospin proton=t/\n"
         " &scalarPs alphaE1fm=0.,betaM1fm=0.,varyscalarPs=f/\n"
         " &spinPs gamEEfm=0.,gamEMfm=0.,gamMEfm=0.,gamMMfm=0./\n"
         % (Emin, Emax, Edel, theta, theta, t(anomaly), t(oneN), t(oneD), gM, gE, t(nuclamm)))
    if v2:
        s += " &v2 " + ", ".join("%s=%s" % kv for kv in v2.items()) + " /\n"
    return s


def parse(rundir):
    out = {}
    pat = re.compile(r'^\s*([-+0-9.Ee]+)\s+f\s*(EE|MM|EM|ME)\s+(\d)\s+([+-])\s+\(([^,]+),([^)]+)\)')
    for fn in ('fort.14', 'fort.15', 'fort.16'):
        p = os.path.join(rundir, fn)
        if not os.path.exists(p):
            continue
        for line in open(p):
            m = pat.match(line)
            if not m:
                continue
            E = round(float(m.group(1)), 6)
            key = '%s%s%s' % (m.group(2), m.group(3), m.group(4))
            val = complex(float(m.group(5)), float(m.group(6)))
            out.setdefault(E, {})[key] = val
    return out


def run_pwa(energies_MeV, exe=EXE, keep=False, **kw):
    """energies: list of lab photon energies [MeV]; non-uniform lists are run point by point."""
    E = sorted(energies_MeV)
    if len(E) > 2 and max(abs((E[i + 1] - E[i]) - (E[1] - E[0])) for i in range(len(E) - 1)) > 1e-9:
        out = {}
        for e1 in E:
            out.update(run_pwa([e1], exe=exe, keep=keep, **dict(kw)))
        return out
    if len(E) == 1:
        Emin = Emax = E[0] / 1000
        Edel = 0.01
    else:
        Edel = (E[1] - E[0]) / 1000
        Emin, Emax = E[0] / 1000, E[-1] / 1000
    bkw = {k: kw.pop(k) for k in ('nucl', 'delta', 'mdelta') if k in kw}
    d = tempfile.mkdtemp(prefix='pwa_', dir=os.environ.get('V2SCRATCH', None))
    open(os.path.join(d, 'baryon.dat'), 'w').write(baryon_dat(**bkw))
    open(os.path.join(d, 'in.dat'), 'w').write(in_dat(Emin, Emax, Edel, **kw))
    subprocess.run([os.path.abspath(exe)], cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    res = parse(d)
    if not keep:
        shutil.rmtree(d)
    # convert printed f*1e4*mpi/M  ->  GeV^-1 :   f[M^-1 units] = printed/(1e4*pionm), f[GeV^-1] = that/SMM
    pionm = 0.13957 / SMM
    conv = 1.0 / (1e4 * pionm) / SMM
    out = {E: {k: v * conv for k, v in d_.items()} for E, d_ in res.items()}
    # key by the nearest requested energy (the code accumulates E in floating point)
    keyed = {}
    for Er in E:
        kk = min(out.keys(), key=lambda x: abs(x - Er)) if out else None
        if kk is not None and abs(kk - Er) < 1e-3 * max(1.0, Er):
            keyed[float(Er)] = out[kk]
    return keyed


if __name__ == '__main__':
    r = run_pwa([250.0], oneD=False, anomaly=False, delta=(False, False))
    for E, d in r.items():
        print(E, d)
