import sys, os, numpy as np, subprocess, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, kin, formtools as ft
cs = ev.consts()
t0 = time.time()
m = assemble.Model(3, cs)
print('compiled in', round(time.time() - t0, 1), 's;  L0, L1 =', m.L0, m.L1)
for k, v in m.vc.items(): print('  vertex', k, v)
# 8-structure Gram for the code amplitude
S8 = []
for k in range(8):
    S8.append(ft.compile_proj(ft.parse_fortran(open(os.path.join(os.path.dirname(__file__), '..', 'out', f'gram{k+1}.log')).read())))
gafac = cs['gA'] ** 2 / (4 * np.pi * cs['F']) ** 2
pts = [(0.05, -0.008), (0.1, -0.05), (0.2, -0.1), (0.3, -0.3), (0.15, -0.002)]
inp = '\n'.join(f'{x} {t} {cs["mpi2"]} {gafac}' for x, t in pts)
out = subprocess.run([os.path.join(os.path.dirname(__file__), 'covp3.exe')], input=inp, capture_output=True, text=True).stdout
rows = [l for l in out.split('\n') if l.strip() and not l.strip().startswith('=') and 'Z' not in l and len(l.split()) == 16]
for (x, t), r in zip(pts, rows):
    v = np.array([float(z) for z in r.split()]); Ac = v[0::2] + 1j * v[1::2]
    en = kin.env(x, t, cs['mpi2']); en.update(T={(0, 0, 0, 0): 1.0}, U={}, M=1.0)
    G8 = np.array([[eval(S8[k][j], en) for k in range(8)] for j in range(6)])
    Vc = G8 @ Ac
    A, V, G, parts = m.A(x, t, parts=True)
    print(f'x={x} t={t}:  rel dev of projections = {np.max(np.abs(V - Vc)) / np.max(np.abs(Vc)):.3e}')
    print('   mine  V:', np.round(V, 8))
    print('   code  V:', np.round(Vc, 8))

print('--- is the residual Born-like?  fit V_mine - V_code = a*V(charge Born) + b*V(magnetic Born)')
import importlib
born = ft.compile_proj(ft.parse_fortran(open(os.path.join(os.path.dirname(__file__), '..', 'out', 'born.log')).read()))
for (x, t), r in zip(pts, rows):
    v = np.array([float(z) for z in r.split()]); Ac = v[0::2] + 1j * v[1::2]
    en = kin.env(x, t, cs['mpi2']); en.update(T={(0, 0, 0, 0): 1.0}, U={}, M=1.0)
    G8 = np.array([[eval(S8[k][j], en) for k in range(8)] for j in range(6)])
    Vc = G8 @ Ac
    A, V, G, parts = m.A(x, t, parts=True)
    d = V - Vc
    e1 = dict(en, Q=1.0, kap=0.0); e2 = dict(en, Q=0.0, kap=1.0); e3 = dict(en, Q=1.0, kap=1.0)
    b1 = np.array([eval(c, e1) for c in born]); b3 = np.array([eval(c, e3) for c in born]); b2 = np.array([eval(c, e2) for c in born])
    bq = b1; bk = (b3 - b1 - b2)   # Q*kap cross term
    X = np.vstack([bq, bk]).T
    coef, res, *_ = np.linalg.lstsq(X, d, rcond=None)
    print(f'x={x} t={t}: |d|/|V| = {np.max(abs(d))/np.max(abs(V)):.2e}, fit a,b = {np.round(coef,8)}, residual {np.max(abs(d - X@coef))/np.max(abs(V)):.1e}')
