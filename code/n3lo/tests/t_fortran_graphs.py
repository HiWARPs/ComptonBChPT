"""t_fortran_graphs.py -- single generated Delta-graph routines (code_v2 dl objects) vs the Python projections."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import delta, devalamp as dv
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'fortran')
names = open(os.path.join(H, 'test_graphs.names')).read().split()
gs = {g.name: g for g in delta.all_graphs(maxgND=0)}
cs = dv.dconsts()
comp = dv.compile_graphs([gs[n] for n in names])
pts = [(0.05, -0.008), (0.1, -0.05), (0.15, -0.002)]
inp = '\n'.join(f'{x} {t} {k}' for k in range(1, len(names) + 1) for x, t in pts)
out = subprocess.run([os.path.join(H, 'test_graphs.exe')], input=inp, capture_output=True, text=True).stdout
rows = [l for l in out.split('\n') if len(l.split()) == 12]
i = 0
for n in names:
    for x, t in pts:
        v = np.array([float(z) for z in rows[i].split()]); W = v[0::2] + 1j * v[1::2]; i += 1
        R = comp[n].raw(dv.env(x, t, cs))
        print(f'{n:24s} x={x} t={t}: rel dev {np.max(abs(W - R)) / np.max(abs(R)):.2e}')
