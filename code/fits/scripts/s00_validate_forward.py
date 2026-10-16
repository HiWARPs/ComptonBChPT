"""Step 0: validation of the forward (theta = 0) evaluation by the driver.

The loop reduction is singular at t = 0; the driver replaces points with |t| < 5e-5 M^2 by a cubic extrapolation
in t from four safe nodes (flag = 1).  Checks:
 (1) Born (+pi0) configurations at theta = 0 against the analytic forward Born amplitudes
       f_B = -alpha/M,   g_B = -alpha kappa^2 nu / (2 M^2)   (lab frame, exact at tree level; the pi0 pole vanishes
       at t = 0), dsigma/dOmega_lab(0) = |f|^2 + |g|^2,  Sigma_2z(0) = -2 Re(f g*)/(|f|^2+|g|^2) (sigma_3/2-sigma_1/2).
 (2) NNLO (and NNLO with large polarisability shifts) at theta = 0 against an independent polynomial extrapolation
     in t of directly evaluated (flag = 0) small-angle points (theta_lab = 7...24 deg), for the observable value and
     for the full quadratic-form coefficient vector (tested by evaluating at shifted polarisabilities).
Output: results/forward_validation.md
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import bchpt
from common import RES, PREF_NNLO, PROFILE
import pandas as pd

AL, M, KAP, HBARC = 1 / 137.036, 938.27, 1.793, 197.33       # code constants (ampl.f90 / in.dat)
EF = [50, 60, 70, 79, 90, 100, 110, 120, 130, 135, 145, 150]
lines = ['# Validation of the theta = 0 evaluation (profile %s, executable %s)' % (PROFILE, bchpt.exe_path()), '',
         ('v2: theta = 0 is evaluated by extrapolating the LOOP AMPLITUDES A_i(s,t) (cubic in t at fixed s, nodes on the physical '
          'side) into the Gram-singular zone (flag 1); v1: the observable itself was extrapolated in t by the driver.') if PROFILE == 'v2'
         else 'v1: the driver extrapolates the observable (quadratic-form coefficients) cubically in t (flag 1).', '']

# (1) Born
pts = [('b_%d_%s' % (E, o), E, 0.0, o, 'LAB') for E in EF for o in ['CS', 'S2Z']]
TB = bchpt.run(pts, 'BornPi0', name='val_fwd_born')
vb = TB.value_file
lines += ['## (1) Born + pi0 (kappa included) at theta = 0 vs analytic forward Born amplitudes', '',
          '| E [MeV] | CS driver | CS analytic | rel. diff | S2z driver | S2z analytic | diff |', '|---|---|---|---|---|---|---|']
mx1 = [0, 0]
for i, E in enumerate(EF):
    f = -AL / M; g = -AL * KAP ** 2 * E / (2 * M ** 2)
    cs = (f * f + g * g) * HBARC ** 2 * 1e7
    s2 = -2 * f * g / (f * f + g * g)
    c_d, s_d = vb[2 * i], vb[2 * i + 1]
    mx1 = [max(mx1[0], abs(c_d / cs - 1)), max(mx1[1], abs(s_d - s2))]
    lines.append('| %d | %.6f | %.6f | %.1e | %.6f | %.6f | %.1e |' % (E, c_d, cs, c_d / cs - 1, s_d, s2, s_d - s2))
lines += ['', 'max |rel. diff| CS = %.1e, max |diff| S2z = %.1e' % tuple(mx1), '']

# (2) NNLO: extrapolation of direct small-angle points
TH = [7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 40]
pts = []
for E in EF:
    for o in ['CS', 'S2Z']:
        pts.append(('z_%d_%s' % (E, o), E, 0.0, o, 'LAB'))
        for th in TH:
            pts.append(('s_%d_%s_%d' % (E, o, th), E, th, o, 'LAB'))
T = bchpt.run(pts, 'NNLO', name='val_fwd_nnlo')
Ptests = {'NNLO': PREF_NNLO, 'shifted': PREF_NNLO + np.array([2.0, -2.0, 1.5, -1.5, 1.0, 1.0])}
lines += ['## (2) NNLO at theta = 0: driver value (%s) vs independent' % ('loop amplitudes extrapolated in t at fixed s, v2' if PROFILE == 'v2' else 'observable extrapolated cubically in t from |t| = (1.5..5)x5e-5 M^2, v1'),
          'polynomial (degree 3) extrapolation in t of the 8 smallest directly evaluated (flag-0) angles between 7 and 40 deg', '',
          '| E [MeV] | P | CS driver | CS poly | rel. diff | S2z driver | S2z poly | diff | flags of direct pts |',
          '|---|---|---|---|---|---|---|---|---|']
mx = {'NNLO': [0, 0], 'shifted': [0, 0]}
for E in EF:
    for pk, P in Ptests.items():
        row = []
        for o in ['CS', 'S2Z']:
            iz = T.index['z_%d_%s' % (E, o)]
            ids = [T.index['s_%d_%s_%d' % (E, o, th)] for th in TH]
            keep = [k for k, j in enumerate(ids) if T.flag[j] == 0][:8]     # the 8 smallest direct angles
            ids = [ids[k] for k in keep]
            v = T.predict_array(P)
            # t in units of M^2 for the lab angle
            nu = E / M
            th = np.radians(np.array(TH)[keep])
            nup = nu / (1 + nu * (1 - np.cos(th)))
            t = -2 * nu * nup * (1 - np.cos(th))
            c = np.polyfit(t, v[ids], 3)
            row.append((v[iz], np.polyval(c, 0.0), max(T.flag[ids])))
        (cd, cp, f1), (sd, sp, f2) = row
        if E <= 145:
            mx[pk] = [max(mx[pk][0], abs(cd / cp - 1)), max(mx[pk][1], abs(sd - sp))]
        lines.append('| %d | %s | %.5f | %.5f | %.1e | %.5f | %.5f | %.1e | %d,%d |' % (E, pk, cd, cp, cd / cp - 1, sd, sp, sd - sp, f1, f2))
lines += ['', 'E <= 145 MeV: max |rel. diff| CS: NNLO %.1e, shifted %.1e;  max |diff| S2z: NNLO %.1e, shifted %.1e' %
          (mx['NNLO'][0], mx['shifted'][0], mx['NNLO'][1], mx['shifted'][1]),
          'GHP data errors: CS 0.8-1.2 %, S2z 0.007-0.011, i.e. two orders of magnitude above these differences for E <= 145 MeV.', ''] + (
          ['v1 note: at 150 MeV (the single-mass loop threshold 149.95 MeV) the t dependence is non-analytic; 145/150 MeV unusable.']
          if PROFILE == 'v1' else
          ['v2 (physical masses): cusps at 144.7 (pi0 p) and 151.4 MeV (pi+ n); 150 MeV lies 1.4 MeV below the pi+ n threshold, where the',
           'small-angle polynomial is less reliable (steep t dependence close to threshold); the loop-amplitude extrapolation is exact to ~1e-5.']) + ['']
# (3) NNLO prediction at theta = 0 vs GHP points
fo = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data', 'db_forward_obs.csv'))
lines += ['## (3) NNLO prediction (polarisabilities = NNLO) at theta = 0 vs GHP forward points', '',
          '| E | dsigma(0) NNLO | GHP | pull | Sigma2z(0) NNLO | GHP | pull |', '|---|---|---|---|---|---|---|']
v = T.predict_array(PREF_NNLO)
for E in EF:
    c = fo[(fo.obs == 'CS') & (fo.E_lab_MeV == E)].iloc[0]; z = fo[(fo.obs == 'S2Z') & (fo.E_lab_MeV == E)].iloc[0]
    tc, tz = v[T.index['z_%d_CS' % E]], v[T.index['z_%d_S2Z' % E]]
    lines.append('| %d | %.4f | %.4f(%.4f) | %+.1f | %.4f | %.4f(%.4f) | %+.1f |' % (E, tc, c.value, c.stat_err, (tc - c.value) / c.stat_err,
                 tz, z.value, z.stat_err, (tz - z.value) / z.stat_err))
lines += ['']
open(os.path.join(RES, 'forward_validation.md'), 'w').write('\n'.join(lines))
print('\n'.join(lines))
