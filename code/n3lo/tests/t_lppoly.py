"""t_lppoly.py -- is (generator - code) for the NNLO pi Delta loops a polynomial in x, t (a scheme difference)?"""
import sys, os, subprocess, numpy as np, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import delta, dassemble, devalamp as dv, formtools as ft, evalamp as ev
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
S8 = [ft.compile_proj(ft.parse_fortran(open(os.path.join(HERE, '..', 'out', f'gram{k+1}.log')).read())) for k in range(8)]
cs = dv.dconsts()
gs = [g for g in delta.all_graphs(maxgND=0) if g.order == Fr(7, 2)]
FUDGE = {'1PI': 1.0, 'R2': -2.0, 'R1': 2 / 3, 'vertex': -2.0, 'self': 2 / 3}
m = dassemble.DModel(gs, cs, iso_override=FUDGE, tag='lp')
gafac = cs['hA'] ** 2 / (4 * cs['MD'] * np.pi * cs['F']) ** 2
pts = [(x, t) for x in (0.06, 0.1, 0.14, 0.18, 0.22, 0.26) for t in (-0.03, -0.07, -0.11, -0.17, -0.23)]
inp = '\n'.join(f'{x} {t} {cs["mpi2"]} {gafac}' for x, t in pts)
out = subprocess.run([os.path.join(HERE, 'covd.exe')], input=inp, capture_output=True, text=True).stdout
rows = [l for l in out.split('\n') if len(l.split()) == 16]
D, P = [], []
for (x, t), r in zip(pts, rows):
    v = np.array([float(z) for z in r.split()]); Ac = v[0::2] + 1j * v[1::2]
    en = dv.env(x, t, cs); en.update(T={(0, 0, 0, 0): 1.0}, U={})
    G8 = np.array([[eval(S8[k][j], en) for k in range(8)] for j in range(6)])
    G6 = ev.gram(en)
    V, _ = m.V(x, t); Vp, _ = m.V(x, t, pole=True)
    D.append((V - G8 @ Ac).real); P.append(Vp.real)
D, P = np.array(D), np.array(P)
print(len(pts), 'points')
for deg in (2, 3, 4, 5):
    mons = [(i, j) for i in range(deg + 1) for j in range(deg + 1) if i + j <= deg]
    X = np.array([[x ** i * t ** j for i, j in mons] for x, t in pts])
    for k in range(6):
        c, *_ = np.linalg.lstsq(X, D[:, k], rcond=None)
        res = D[:, k] - X @ c
        print(f'deg {deg} A{k+1}: max|diff| {np.max(abs(D[:, k])):.3e}  poly-fit residual {np.max(abs(res)):.2e}')
print('ratio diff / pole part:')
for k in range(6):
    print(f'A{k+1}:', np.round(D[:6, k] / np.where(abs(P[:6, k]) > 1e-12, P[:6, k], np.nan), 4))
