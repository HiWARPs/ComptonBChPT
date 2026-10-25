"""t_dA.py -- pi Delta loops with the complete minimal coupling (photon on the Delta charge): gauge invariance,
UV poles (contact form?), low-energy limit, static polarisabilities, and the correction dA = complete - NNLO code."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import delta, dassemble, devalamp as dv, statpol, ddelta
from fractions import Fraction as Fr

cs = dv.dconsts()
D = ddelta.Delta45(cs, which=('A',))
mf = D.m['A_full']
if 'nogi' not in sys.argv:
    mf.gi_compile()
    for x, t in [(0.1, -0.05), (0.25, -0.2)]:
        W, big = mf.GI(x, t)
        print(f'gauge invariance x={x} t={t}: max|sum| = {np.max(abs(W)):.2e}, largest graph {big:.2e}')
def lec(x, t, da, db):
    y = -x - t / 2
    return np.array([0.5 * (da + db) * (x ** 2 + y ** 2) + 0.5 * t * db, -db - (da + db) * t / 8, 0, 0, 0, 0])
rows, rhs = [], []
for x, t in [(0.1, -0.05), (0.2, -0.1), (0.3, -0.3), (0.05, -0.01)]:
    Ap = mf.A(x, t, pole=True); A = mf.A(x, t)
    print(f'x={x} t={t}: |A_pole|/|A| = {np.max(abs(Ap)) / np.max(abs(A)):.2e}  A_pole = {np.round(Ap.real, 6)}')
    rows += list(np.vstack([lec(x, t, 1, 0), lec(x, t, 0, 1)]).T); rhs += list(Ap.real)
c, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
print(f'UV pole as delta-alpha, delta-beta contact term: {c}, max residual {np.max(abs(np.array(rhs) - np.array(rows) @ c)):.2e}')
for x in (0.02, 0.01, 0.005):
    print(f'low energy x={x}: A = {np.round(mf.A(x, -0.5 * x * x).real, 7)}')
pp_, dpp = statpol.extract(lambda x, t: mf.A(x, t, pole=True))
print(f'A_full UV-pole coefficient as static pols {np.round(pp_, 5)} (spread {dpp.max():.1e})')
for name in ('A_full', 'A_lp'):
    p, dp = statpol.extract(lambda x, t: D.m[name].A(x, t))
    print(f'{name:7s} static pols {np.round(p, 4)}  g0,gpi {np.round(statpol.g0pi(p), 4)}  spread {dp.max():.1e}')
p, dp = statpol.extract(lambda x, t: D.parts(x, t)['dA'])
print(f'dA      static pols {np.round(p, 4)}  g0,gpi {np.round(statpol.g0pi(p), 4)}  spread {dp.max():.1e}')
