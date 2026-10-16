"""Per-beamtime A2 Compton data (March 2018 and July 2018 beamtimes) from
E. Mornacchi, PhD thesis, JGU Mainz (2021), Appendix A (unpolarized dsigma/dOmega_lab,
Tables A.1-A.10) and Appendix B (beam asymmetry Sigma_3, Tables B.1-B.6).
The PRL 128 (2022) 132503 supplement gives the combined (stat-weighted) values.
Usage: python3 -I build_db_mornacchi_thesis.py <thesis.txt> <outdir>
"""
import re, sys, csv, math, os
txt, out = sys.argv[1], sys.argv[2]
L = open(txt, encoding='utf-8', errors='replace').read().splitlines()
start = next(i for i, l in enumerate(L) if l.strip().startswith('Appendix A'))
Mp = 938.272
def th_cm(E, thl):
    # exact 2-body Compton kinematics: lab -> cm angle
    E, c = E/1000.0, math.cos(math.radians(thl)); M = Mp/1000.0
    Ep = E/(1 + E/M*(1 - c)); t = -2*E*Ep*(1 - c); s = M*M + 2*M*E
    w = (s - M*M)/(2*math.sqrt(s)); cc = 1 + t/(2*w*w)
    return math.degrees(math.acos(max(-1, min(1, cc))))
hdr = ['dataset','texkey','obs','E_lab_MeV','dE_MeV','theta_lab_deg','theta_cm_deg','frame_of_cs',
       'value','stat_err','sys_err_pt','norm_err_frac','units','source','notes']
rows, tab, beam, obs = [], None, None, None
num = re.compile(r'^\s*(\d+\.\d)\s+(\d+\.\d)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s*$')
for i in range(start, len(L)):
    l = L[i]
    m = re.match(r'\s*Table ([AB])\.(\d+):', l)
    if m:
        cap = l + ' ' + L[i+1]
        tab = f'{m.group(1)}.{m.group(2)}'
        beam = 'Mar2018' if 'March 2018' in cap else ('Jul2018' if 'July 2018' in cap else None)
        obs = 'CS' if m.group(1) == 'A' else 'S3'
        er = re.search(r'ωγ = ([\d.]+) − ([\d.]+)', cap)
        rng = (float(er.group(1)), float(er.group(2))) if er else (None, None)
        continue
    m = num.match(l)
    if m and tab and beam:
        E, thl, v, st, sy, co = map(float, m.groups())
        if obs == 'CS':
            nf = 0.03; units = 'nb/sr'; frame = 'lab'
            note = f'thesis Tab. {tab}; corr. 3% abs = {co}'
        else:
            nf = 0.05; units = '1'; frame = 'na'
            note = f'thesis Tab. {tab}; corr. (5/sqrt3)% abs = {co}; scale err uniform half-width 5% (rms 2.9%)'
        dE = f'{rng[1]-rng[0]:.1f}' if rng[0] else ''
        rows.append([f'A2_Mornacchi2021thesis_{obs}_{beam}', 'Mornacchi:2021thesis', obs, E, dE, thl,
                     f'{th_cm(E, thl):.3f}', frame, v, st, sy, nf, units,
                     'E. Mornacchi PhD thesis, JGU Mainz 2021, App. ' + tab[0], note])
with open(os.path.join(out, 'db_mornacchi_thesis.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(hdr); w.writerows(rows)
from collections import Counter
print(len(rows), Counter((r[0], r[3]) for r in rows))
