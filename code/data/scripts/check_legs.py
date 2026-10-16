#!/usr/bin/env python3
"""Check LEGS transcription against pdftotext token stream (page 46) and author's legacy files."""
import sys, re
tr = [l.strip().split(',') for l in open(sys.argv[1]) if l[:2] in ('II', 'II')]
rows = {}
for t, E, th, v, e in tr:
    rows.setdefault((t, E), []).extend([v, e])
txt = open(sys.argv[2]).read()
iII = txt.find('TABLE II. Point cross'); iIII = txt.find('TABLE III. Point beam'); iIV = txt.find('TABLE IV. Point cross')
bad = 0; n = 0
for tab, seg in (('II', txt[iII:iIII]), ('III', txt[iIII:iIV])):
    for line in seg.split('\n'):
        tok = line.split()
        if tok and re.match(r'^\d{3}\.\d$', tok[0]):
            nums = re.findall(r'-?\d+\.\d+', line)[1:]
            mine = rows.get((tab, tok[0]))
            n += 1
            if mine != nums:
                bad += 1; print('DIFF', tab, tok[0], '\n  pdf :', nums, '\n  mine:', mine)
print('rows compared', n, 'diffs', bad)
# legacy
leg = open(sys.argv[3]).read().split('\n\n')
d = {(t, float(E), int(th)): (float(v), float(e)) for t, E, th, v, e in tr}
checks = {'a90.legs01': ('II', 'E', 90), 'asym90.legs01': ('III', 'E', 90), 'e_286.legs': ('II', 'th', 286.5), 'sigma286.legs01': ('III', 'th', 286.5)}
for b in leg:
    h = b.split('\n')[0]
    for k, (t, mode, fix) in checks.items():
        if k in h:
            nn = 0; nb = 0
            for l in b.strip().split('\n')[1:]:
                x, v, e = map(float, l.split()[:3])
                key = (t, x, fix) if mode == 'E' else (t, fix, int(x))
                nn += 1
                if key not in d or abs(d[key][0] - v) > 1e-9 or abs(d[key][1] - e) > 1e-9:
                    nb += 1; print('LEGACY-DIFF', k, x, v, e, d.get(key))
            print('legacy', k, 'n', nn, 'mismatch', nb)
