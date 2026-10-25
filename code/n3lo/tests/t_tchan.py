import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import evalamp as ev, tchan, runform, formtools as ft, kin
cs = ev.consts(dict(ev.LMP15, c1=1.0, c2=1.0, c3=1.0))
gs = tchan.graphs()
cmp = ev.compile_graphs(gs)
gi = runform.project_all(gs, gi=True)
gic = {n: [compile(ft.to_python(d.get(f'V{j}', '0')), 'w', 'eval') for j in range(1, 9)] for n, d in gi.items()}
def lec(x, t, da, db):
    y = -x - t/2
    return np.array([0.5*(da+db)*(x**2+y**2)+0.5*t*db, -db-(da+db)*t/8, 0, 0, 0, 0])
for ci in ('c1', 'c2', 'c3'):
    rows, rhs = [], []
    for x, t in [(0.1, -0.05), (0.2, -0.1), (0.3, -0.3), (0.15, -0.02)]:
        en = kin.env(x, t, cs['mpi2']); en['M'] = 1.0
        V = np.zeros(6, complex); Vp = np.zeros(6, complex); W = np.zeros(8, complex); big = 0
        for g in gs:
            if g.coup.get(ci) is None: continue
            c = cmp[g.name]; f = c.coup(cs) * c.isof(cs)
            V += c.raw(en) * f; Vp += c.raw(en, pole=True) * f
            T, U = c.tensors(en); e = dict(en, T=T, U=U)
            w = np.array([eval(z, e) for z in gic[g.name]]) * f; W += w; big = max(big, np.max(abs(w)))
        G = ev.gram(en)
        A = np.linalg.solve(G, V); Ap = np.linalg.solve(G, Vp)
        print(f'{ci} x={x} t={t}: GI {np.max(abs(W))/big:.1e}  A={np.round(A.real,6)}  A_pole={np.round(Ap.real,6)}')
        rows += list(np.vstack([lec(x, t, 1, 0), lec(x, t, 0, 1)]).T); rhs += list(Ap.real)
    c, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    print(f'   pole = contact form? (da, db) = {np.round(c, 8)}, residual {np.max(abs(np.array(rhs) - np.array(rows) @ c)):.1e}')
