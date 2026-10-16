#!/usr/bin/env python3
"""Cross-check the OdL thesis transcription against the author's independent legacy files
(a0xx.leon01, Olmos_*.dat embedded in xmgrace files).
Usage: python3 -I check_odl.py transcription_weighted_mean.csv author_agr_embedded_data.txt"""
import sys
tr = [l.strip().split(',') for l in open(sys.argv[1]) if l[0].isdigit()]
tr = {(int(a), float(e)): (float(v), float(s)) for a, e, v, s in tr}
leg = open(sys.argv[2]).read().split('\n\n')
mp = {'a059.leon01': 59, 'a085.leon01': 85, 'a133.leon01': 133, 'Olmos_107.dat': 107, 'Olmos_155.dat': 155,
      'Olmos_59.dat': 59, 'Olmos_85.dat': 85, 'Olmos_135.dat': 133}
n = bad = 0; covered = set()
for b in leg:
    for k, ang in mp.items():
        if k in b.split('\n')[0]:
            for l in b.strip().split('\n')[1:]:
                E, v, s = map(float, l.split()[:3])
                n += 1; covered.add((ang, E))
                if (ang, E) not in tr or abs(tr[(ang, E)][0] - v) > 1e-9 or abs(tr[(ang, E)][1] - s) > 1e-9:
                    bad += 1; print('LEGACY-diff', k, E, v, s, tr.get((ang, E)))
print('legacy points compared:', n, 'mismatches:', bad, '; distinct transcription points covered:', len(covered), 'of', len(tr))
