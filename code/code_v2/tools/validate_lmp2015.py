#!/usr/bin/env python3
"""
Validation of the Delta-region physics of ampl.exe v2 against the LMP2015 figures (the draft
figure files LP_ms3/*.pdf contain the vector paths of the published curves/bands):
  * Sigma_2x, Sigma_3, Sigma_2z vs theta_cm at 224, 245, 265, 286, 310, 334 MeV (LMP2015 Figs. 15-17):
    max |v2 - band mid-line| with and without the gamma N Delta vertex corrections;
  * Im of the l=1 Compton multipoles (LMP2015 Fig. 10): red solid (with vertex corrections),
    red dashed (without), MAID.
usage: validate_lmp2015.py [--figdir ../../../LP_ms3]
"""
import os, sys, subprocess, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'figdigit'))
import bands
FIGDIR = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'LP_ms3'))
if '--figdir' in sys.argv:
    FIGDIR = sys.argv[sys.argv.index('--figdir') + 1]
DRV = os.path.join(HERE, '..', 'driver', 'run_points.sh')
ES = [224, 245, 265, 286, 310, 334]
THS = np.arange(10, 180, 10)


def run(keys):
    d = tempfile.mkdtemp(prefix='val_')
    p = os.path.join(d, 'p.txt'); o = os.path.join(d, 'p.out')
    with open(p, 'w') as f:
        for E in ES:
            for th in THS:
                for ob in ('S2X', 'S3', 'S2Z'):
                    f.write('%s_%d_%d %d %d %s CM\n' % (ob, E, th, E, th, ob))
    subprocess.run([DRV, p, o] + keys, check=True, capture_output=True)
    res = {}
    for l in open(o):
        if not l.startswith('#'):
            q = l.split(); res[q[0]] = float(q[5])
    return res


def main():
    A = run(['dvertex=f']); B = run(['dvertex=t'])
    figs = {'S2X': 'fig_theory_S2x_v2', 'S3': 'fig_test_legs_S3_v2', 'S2Z': 'fig_theory_S2z_v2'}
    print('max |v2 - LMP2015 band mid-line| over theta_cm = 10..170 deg   [without VC | with VC]  (band half-width)')
    for ob, fn in figs.items():
        P = bands.panels(os.path.join(FIGDIR, fn + '.pdf'), ES)
        for E in ES:
            th, lo, hi, _ = P[E]
            mid = np.interp(THS, th, (lo + hi) / 2); w = np.interp(THS, th, (hi - lo) / 2)
            a = np.array([A['%s_%d_%d' % (ob, E, t)] for t in THS]); b = np.array([B['%s_%d_%d' % (ob, E, t)] for t in THS])
            print('  %-4s E=%3d MeV:  %.4f | %.4f   inside band: %3.0f%% | %3.0f%%   (median half-width %.3f)' %
                  (ob, E, np.max(np.abs(a - mid)), np.max(np.abs(b - mid)), 100 * np.mean(np.abs(a - mid) <= w),
                   100 * np.mean(np.abs(b - mid) <= w), np.median(w)))
    import cmp_multipoles_lmp
    print('\nIm l=1 Compton multipoles vs LMP2015 Fig. 10 (draft fig_pion_production_multipoles_l1_LS_v2):')
    cmp_multipoles_lmp.main()


if __name__ == '__main__':
    main()
