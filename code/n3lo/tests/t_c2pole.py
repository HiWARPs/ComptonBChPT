import sys, os, numpy as np, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import evalamp as ev, tchan, kin
cs = ev.consts(dict(ev.LMP15, c1=0.0, c2=1.0, c3=0.0))
gs = [g for g in tchan.graphs() if 'c2' in g.coup]
cmp = ev.compile_graphs(gs)
pts = [(x, t) for x in (0.08, 0.15, 0.22, 0.3) for t in (-0.01, -0.05, -0.1)]
A1, A2, X = [], [], []
for x, t in pts:
    en = kin.env(x, t, cs['mpi2']); en['M'] = 1.0
    Vp = sum(cmp[g.name].raw(en, pole=True) * cmp[g.name].coup(cs) * cmp[g.name].isof(cs) for g in gs)
    Ap = np.linalg.solve(ev.gram(en), Vp).real
    xi = (2 * x + t / 2) / 2
    A1.append(Ap[0]); A2.append(Ap[1]); X.append((xi, t))
    assert np.max(abs(Ap[2:])) < 1e-10
def fit(vals, mons):
    M = np.array([[xi**a * t**b for a, b in mons] for xi, t in X])
    c, *_ = np.linalg.lstsq(M, np.array(vals), rcond=None)
    return c, np.max(abs(np.array(vals) - M @ c))
mons = [(a, b) for a in (0, 2, 4) for b in (0, 1, 2) if a + 2 * b <= 4]
for nm, v in (('A1', A1), ('A2', A2)):
    c, r = fit(v, mons)
    print(nm, 'pole = ', ' + '.join(f'{ci:.6f} xi^{a} t^{b}' for ci, (a, b) in zip(c, mons) if abs(ci) > 1e-9), '  resid', f'{r:.1e}')
print('16 pi^2 F^2 =', 16 * np.pi**2 * cs['F']**2)
