"""t_dB.py -- the gamma N Delta graphs of O(p^9/2): gauge invariance, UV poles, renormalisation of the Delta-pole
graphs with the one-loop gamma N -> Delta vertex (on-shell: Re residue at s = MD^2 unchanged), static polarisabilities
compared with the O(p^4) kappa insertions."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import delta, dassemble, devalamp as dv, statpol, ddelta, evalamp as ev

cs = dv.dconsts()
D = ddelta.Delta45(cs, which=('B',))
for k in ('Bloop', 'Bpole'):
    m = D.m[k]
    print(f'{k}: {len(m.graphs)} graphs, orders {sorted({str(g.order) for g in m.graphs})}')
    if 'nogi' not in sys.argv:
        m.gi_compile()
        for x, t in [(0.1, -0.05), (0.25, -0.2)]:
            W, big = m.GI(x, t)
            print(f'  gauge invariance x={x} t={t}: max|sum| = {np.max(abs(W)):.2e}, largest graph {big:.2e}')
    for x, t in [(0.1, -0.05), (0.2, -0.1), (0.05, -0.01)]:
        Ap = m.A(x, t, pole=True); A = m.A(x, t)
        print(f'  x={x} t={t}: A_pole = {np.round(Ap.real, 7)}   A = {np.round(A.real, 7)}')

# ---- tree Delta pole for the counterterms
tree = {}
for c1, c2 in [('gM', 'gM'), ('gM', 'gE'), ('gE', 'gM'), ('gE', 'gE')]:
    tree[(c1, c2)] = delta.dpole_graphs(c1, c2, tag='ct')
tc = dv.compile_graphs([g for v in tree.values() for g in v])
def tree_amp(x, t, gM, gE):
    cc = dict(cs, gM=gM, gE=gE)
    en = dv.env(x, t, cs)
    V = sum(tc[g.name].raw(en) * tc[g.name].coup(cc) * tc[g.name].isof(cc) for v in tree.values() for g in v)
    return np.linalg.solve(ev.gram(en), V)
def dtree(x, t):
    """dA/dgM, dA/dgE of the tree Delta pole at the physical couplings"""
    h = 1e-4
    dM = (tree_amp(x, t, cs['gM'] + h, cs['gE']) - tree_amp(x, t, cs['gM'] - h, cs['gE'])) / (2 * h)
    dE = (tree_amp(x, t, cs['gM'], cs['gE'] + h) - tree_amp(x, t, cs['gM'], cs['gE'] - h)) / (2 * h)
    return dM, dE

# ---- on-shell renormalisation of Bpole: Re residue at s = MD^2
mp = D.m['Bpole']
xD = (cs['MD'] ** 2 - 1) / 2
rows, rhs = [], []
for t in (-0.02, -0.08, -0.15):
    res = []
    for eta in (1e-3, -1e-3, 2e-3, -2e-3):
        x = xD * (1 + eta)
        res.append((2 * x - 2 * xD) * mp.A(x, t))
    R = (4 * (res[0] + res[1]) / 2 - (res[2] + res[3]) / 2) / 3      # Richardson in eta
    rt = []
    for eta in (1e-3, -1e-3):
        x = xD * (1 + eta)
        dM, dE = dtree(x, t)
        rt.append(((2 * x - 2 * xD) * dM, (2 * x - 2 * xD) * dE))
    RM = (rt[0][0] + rt[1][0]) / 2; RE = (rt[0][1] + rt[1][1]) / 2
    rows += list(np.vstack([RM.real, RE.real]).T); rhs += list(-R.real)
    print(f't={t}: Re residue Bpole {np.round(R.real, 6)}  Im {np.round(R.imag, 6)}')
cM, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
print(f'on-shell counterterms delta gM, delta gE = {cM}, max residual {np.max(abs(np.array(rhs) - np.array(rows) @ cM)):.2e}'
      f' (of max residue {np.max(abs(np.array(rhs))):.2e})')

def Bpole_os(x, t):
    dM, dE = dtree(x, t)
    return mp.A(x, t) + cM[0] * dM + cM[1] * dE

for k in ('Bloop', 'Bpole'):
    pp_, dpp = statpol.extract(lambda x, t: D.m[k].A(x, t, pole=True))
    print(f'{k:16s} UV-pole coefficient as static pols {np.round(pp_, 5)} (spread {dpp.max():.1e}); '
          f'nonzero gamma entries need O(p^5) spin-pol counterterms')
pk = {'kappa (O(p^4), for comparison)': np.array([-2.84, 3.50, -0.01, 0.48, -0.24, 0.55])}
out = {'Bloop': lambda x, t: D.m['Bloop'].A(x, t), 'Bpole MS-bar': lambda x, t: mp.A(x, t), 'Bpole on-shell': Bpole_os}
for k, f in out.items():
    p, dp = statpol.extract(f)
    print(f'{k:16s} static pols {np.round(p, 4)}  g0,gpi {np.round(statpol.g0pi(p), 4)}  spread {dp.max():.1e}')
for k, p in pk.items():
    print(f'{k:16s} {p}  g0,gpi {np.round(statpol.g0pi(p), 3)}')
