"""running2_observables.py -- effect of the running Delta parameters (code_v2 &v2 irun = 1, 2, 3; work/notes/delta_running2.md)
on the parameter-free NNLO prediction (main v2 settings of the fits: physical pion/neutron masses, vertex corrections):
dsigma/dOmega_lab at theta_lab = 60, 90, 133.5 deg, c.m. cross sections at 75, 90, 138 deg, Sigma_3 at 90 deg c.m.,
Sigma_2x, Sigma_2z at theta_lab = 90 deg, for E_lab = 100 ... 340 MeV; below threshold (100, 120, 140 MeV) the lab cross
sections at 60/90/135 deg.  Two evaluations of every running variant: (a) at its own static polarisabilities (driver
Pref of that configuration: irun=1 changes all six, irun=2 gamma_E1E1 and gamma_M1M1, irun=3 none); (b) at the static
polarisabilities of the default amplitude (contact terms restore the NNLO static values; isolates the energy dependence).
Run: .venv/bin/python running2_observables.py  -> running2_observables.{json,md}"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
FITS = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(FITS, 'lib'))
import numpy as np
import pandas as pd
import bchpt

ES = [100.0, 150.0, 200.0, 250.0, 280.0, 300.0, 320.0, 340.0]
PTS = []
for E in ES:
    for th in (60.0, 90.0, 133.5):
        PTS.append(('CS_%g_%g' % (E, th), E, th, 'CS', 'LAB'))
    for th in (75.0, 90.0, 138.0):
        PTS.append(('CSCM_%g_%g' % (E, th), E, th, 'CSCM', 'CM'))
    PTS.append(('S3_%g_90cm' % E, E, 90.0, 'S3', 'CM'))
    PTS.append(('S2X_%g_90' % E, E, 90.0, 'S2X', 'LAB'))
    PTS.append(('S2Z_%g_90' % E, E, 90.0, 'S2Z', 'LAB'))
for E in (100.0, 120.0, 140.0):
    for th in (60.0, 90.0, 135.0):
        PTS.append(('CSlow_%g_%g' % (E, th), E, th, 'CS', 'LAB'))

CF = ['NNLO', 'running1', 'running2', 'running3']
T = {c: bchpt.run(PTS, c, name='run2obs_' + c) for c in CF}
P0 = T['NNLO'].Pref.copy()
OUT = dict(Pref={c: T[c].Pref.tolist() for c in CF}, points=[list(p) for p in PTS], values={})
for c in CF:
    OUT['values'][c + '_ownPref'] = T[c].predict_array(T[c].Pref).tolist()
    OUT['values'][c + '_P0'] = T[c].predict_array(P0).tolist()
json.dump(OUT, open(os.path.join(HERE, 'running2_observables.json'), 'w'), indent=1)

ref = np.array(OUT['values']['NNLO_ownPref'])
ids = [p[0] for p in PTS]
L = ['# Running Delta parameters: observables of the parameter-free NNLO prediction', '',
     'Produced by `running2_observables.py` (fits library, main v2 settings: physical masses, vertex corrections).',
     'Static polarisabilities (driver Pref) per configuration (alpha, beta [1e-4 fm^3]; gammas [1e-4 fm^4]):', '']
for c in CF:
    L.append('* %-9s %s' % (c + ':', ' '.join('%8.4f' % v for v in T[c].Pref)))
L += ['', 'Columns: default value; then for irun = 1, 2, 3 the ratio to the default (cross sections) or the difference',
      '(asymmetries), evaluated (a) at the own static polarisabilities of the running amplitude [own] and (b) at the',
      'static polarisabilities of the default amplitude [P0] (contact terms restoring the NNLO values).', '']


def row_tab(sel, title, isratio):
    L.append('## ' + title); L.append('')
    hdr = '| E [MeV] | theta | default | r1 own | r1 P0 | r2 own | r2 P0 | r3 own | r3 P0 |'
    L.append(hdr); L.append('|---|---|---|---|---|---|---|---|---|')
    for i in sel:
        p = PTS[i]
        vals = []
        for c in CF[1:]:
            for k in ('_ownPref', '_P0'):
                v = OUT['values'][c + k][i]
                vals.append(v / ref[i] if isratio else v - ref[i])
        fmt = '%.4f' if isratio else '%+.4f'
        L.append('| %g | %g%s | %.4f | %s |' % (p[1], p[2], ' cm' if p[4] == 'CM' else '', ref[i], ' | '.join(fmt % v for v in vals)))
    L.append('')


for obs, title, isr in [('CS', 'dsigma/dOmega_lab [nb/sr], theta_lab = 60, 90, 133.5 deg (ratios)', True),
                        ('CSCM', 'dsigma/dOmega_cm [nb/sr], theta_cm = 75, 90, 138 deg (ratios)', True),
                        ('S3', 'Sigma_3 at 90 deg c.m. (differences)', False),
                        ('S2X', 'Sigma_2x at theta_lab = 90 deg (differences)', False),
                        ('S2Z', 'Sigma_2z at theta_lab = 90 deg (differences)', False),
                        ('CSlow', 'below threshold: dsigma/dOmega_lab at 60, 90, 135 deg, 100-140 MeV (ratios)', True)]:
    sel = [i for i, p in enumerate(PTS) if p[0].startswith(obs + '_')]
    row_tab(sel, title, isr)
open(os.path.join(HERE, 'running2_observables.md'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
