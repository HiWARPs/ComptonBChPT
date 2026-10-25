import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lt'))
import ltpy, evalamp as ev, tchan, statpol, kin
cs = ev.consts(dict(ev.LMP15, c1=0.0, c2=1.38, c3=0.0))
gs = [g for g in tchan.graphs() if 'c2' in g.coup]
cmp = ev.compile_graphs(gs)
def A(x, t):
    en = kin.env(x, t, cs['mpi2']); en['M'] = 1.0
    V = sum(cmp[g.name].raw(en) * cmp[g.name].coup(cs) * cmp[g.name].isof(cs) for g in gs)
    return np.linalg.solve(ev.gram(en), V)
for mu in (1.0, 1232/938.27):
    ltpy.setmudim(mu**2); ltpy.clear()
    p, _ = statpol.extract(A)
    print('mu/M=%.3f  alpha=%.4f beta=%.4f' % (mu, p[0], p[1]))
