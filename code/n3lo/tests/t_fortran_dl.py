"""t_fortran_dl.py -- generated Fortran (fortran/dl, test_dl.exe) vs the Python check path for the Delta loops beyond
NNLO: part 1 CovDeltaFull (complete minimal coupling), part 2 CovDeltaGND (gamma N Delta in the pi N loops)."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import ddelta, devalamp as dv
parts = [int(a) for a in sys.argv[1:]] or [1, 2]
pts = [(0.05, -0.008), (0.1, -0.05), (0.2, -0.1), (0.3, -0.3), (0.15, -0.002), (0.22, -0.18), (0.4, -0.5)]
inp = '\n'.join(f'{x} {t} {p}' for p in parts for x, t in pts)
exe = os.path.join(os.path.dirname(__file__), '..', 'fortran', 'test_dl.exe')
out = subprocess.run([exe], input=inp, capture_output=True, text=True).stdout
rows = [l for l in out.split('\n') if len(l.split()) == 12]
cs = dv.dconsts()
D = ddelta.Delta45(cs, which=tuple({1: 'A', 2: 'B'}[p] for p in parts))
k = 0
for p in parts:
    m = D.m['A_full' if p == 1 else 'Bloop']
    for x, t in pts:
        v = np.array([float(z) for z in rows[k].split()]); Af = v[0::2] + 1j * v[1::2]; k += 1
        A = m.A(x, t)
        print(f'part {p} x={x} t={t}: rel dev {np.max(abs(Af - A)) / np.max(abs(A)):.2e}')
