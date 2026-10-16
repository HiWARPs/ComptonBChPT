#!/usr/bin/env python3
"""
End-to-end test of the fit driver against the ORIGINAL main-program path of ampl.exe.

1. run_points.sh computes the quadratic-form coefficients once (reference = the hard-coded
   *NNLOfm constants of ampl.f90, i.e. exactly the convention used by the main program).
2. For several random polarisability sets P_k the unmodified main loop of ampl.exe is run
   (no &fitdriver namelist), and its standard output files are read:
     xcs_Adistribution.output  (dsigma/dOmega_lab, nb/sr), xBa_Adistribution.output (Sigma_3 in %),
     xSigma2x_Adistribution.output, xSigma2z_Adistribution.output (in %).
3. The driver's quadratic forms evaluated at P_k are compared with these numbers.
Expected agreement: ~1e-7 relative (8 significant digits in the main-program output and its
single-precision acos for the lab->c.m. angle), independent of the size of the shifts.

usage: python3 test_quadratic_exactness.py [path/to/ampl.exe] [--hardrefs]
  --hardrefs : use the v1 hard-coded reference values (needed for v1 executables, whose main
               program used them); v2 computes the references analytically in both modes.
"""
import os, sys, subprocess, tempfile, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from compton_qf import DriverOutput

ARGS = [a for a in sys.argv[1:] if not a.startswith('--')]
HARDREFS = '--hardrefs' in sys.argv
EXE = ARGS[0] if ARGS else os.path.join(HERE, '..', 'execs', 'ampl.exe')
HARD = dict(refalphaE1='11.2461d0', refbetaM1='3.87054d0', refgamEE='-3.34375d0',
            refgamMM='2.95934d0', refgamEM='0.185504d0', refgamME='1.14328d0')
ENERGIES = [80.0, 140.0, 250.0, 300.0]
ANGLES = [40.0, 90.0, 140.0]
OBS = ['CS', 'S3', 'S2X', 'S2Z']

def main_program_run(P, E, workdir):
    tmpl = open(os.path.join(HERE, 'in_driver.template')).read()
    txt = (tmpl.replace('@ALPHA@', '%.10g' % P[0]).replace('@BETA@', '%.10g' % P[1])
              .replace('@GEE@', '%.10g' % P[2]).replace('@GMM@', '%.10g' % P[3])
              .replace('@GEM@', '%.10g' % P[4]).replace('@GME@', '%.10g' % P[5])
              .replace('@ONLYLEX@', 'f'))
    lines = [l for l in txt.split('\n') if '&fitdriver' not in l and '&v2' not in l]   # plain main program (v2 defaults)
    lines[0] = ' &energy Emin=%.6f,Emax=%.6f,Edel=1.,kvar=3/' % (E/1000, E/1000)
    lines[1] = ' &angles thmin=%.1f,thmax=%.1f,thdel=%.1f,cmsangle=f/' % (ANGLES[0], ANGLES[-1]+0.1, ANGLES[1]-ANGLES[0])
    open(os.path.join(workdir, 'in.dat'), 'w').write('\n'.join(lines) + '\n')
    subprocess.run([os.path.abspath(EXE)], cwd=workdir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    rd = lambda f: np.loadtxt(os.path.join(workdir, f), ndmin=2)
    cs = rd('xcs_Adistribution.output'); ba = rd('xBa_Adistribution.output')
    s2x = rd('xSigma2x_Adistribution.output'); s2z = rd('xSigma2z_Adistribution.output')
    out = {}
    for i, th in enumerate(ANGLES):
        out[('CS', E, th)] = cs[i, 1]; out[('S3', E, th)] = ba[i, 1] / 100
        out[('S2X', E, th)] = s2x[i, 1] / 100; out[('S2Z', E, th)] = s2z[i, 1] / 100
    return out

def main():
    work = tempfile.mkdtemp(prefix='cbchpt_test_')
    try:
        pts = os.path.join(work, 'pts.txt')
        with open(pts, 'w') as f:
            for E in ENERGIES:
                for th in ANGLES:
                    for o in OBS:
                        f.write('%s_%g_%g %g %g %s LAB\n' % (o, E, th, E, th, o))
        out = os.path.join(work, 'pts.out')
        cmd = [os.path.join(HERE, 'run_points.sh'), pts, out, 'exe=' + os.path.abspath(EXE)]
        if HARDREFS:
            cmd.append('ref=' + ','.join('%s=%s' % kv for kv in HARD.items()))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        T = DriverOutput(out)
        rng = np.random.default_rng(2026)
        Pref = T.Pref
        sets = [Pref.copy(), Pref + np.array([3.0, -2.0, 1.0, -1.5, 0.5, 0.8])]
        sets += [Pref + rng.normal(0, [2, 2, 2, 2, 1, 1]) for _ in range(2)]
        sets += [Pref + np.array([15.0, -12.0, 8.0, -9.0, 6.0, -7.0])]      # very large shifts
        worst = 0.0
        for P in sets:
            pred = T.predict_array(P)
            rel = []
            for E in ENERGIES:
                rd = os.path.join(work, 'main_%g' % E); os.makedirs(rd, exist_ok=True)
                shutil.copy(os.path.join(HERE, 'baryon.dat'), rd)
                direct = main_program_run(P, E, rd)
                for i, pid in enumerate(T.ids):
                    o, e, th = pid.split('_'); e = float(e); th = float(th)
                    if e != E: continue
                    d = direct[(o, e, th)]
                    # cross sections: relative deviation; asymmetries: absolute deviation
                    dev = abs(pred[i] - d) / abs(d) if o == 'CS' else abs(pred[i] - d)
                    rel.append(dev)
                    if dev > 1e-6: print('     %-14s driver=% .9e  main=% .9e  dev=%.1e' % (pid, pred[i], d, dev))
            worst = max(worst, max(rel))
            print('P = %s   max deviation (CS rel., asym. abs.) driver vs main program = %.2e' % (np.round(P, 3), max(rel)))
        print('RESULT: worst deviation %.2e (%s)' % (worst, 'PASS' if worst < 5e-6 else 'CHECK'))
    finally:
        shutil.rmtree(work, ignore_errors=True)

if __name__ == '__main__':
    main()
