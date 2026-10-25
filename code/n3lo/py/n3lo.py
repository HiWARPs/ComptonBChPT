"""n3lo.py -- the O(p^4) piN-loop amplitude (kappa insertions + c_i graphs), Python check path."""
import numpy as np, assemble, evalamp as ev, tchan, kin

CI_FIT2 = dict(c1=-0.99, c2=1.38, c3=-2.33)   # Yao et al., JHEP 05 (2016) 038, Fit II (N + LO Delta), GeV^-1
CI_FIT3 = dict(c1=-1.31, c2=0.78, c3=-2.55)   # same, Fit III (N + Delta)
CI_NODELTA = dict(c1=-1.22, c2=3.58, c3=-6.04)  # same, Fit I (Delta-less)

class P4:
    def __init__(self, ci=CI_FIT2, base=ev.LMP15, nproc=8):
        self.cs = ev.consts(dict(base, **ci))
        self.tg = tchan.graphs()
        self.m = assemble.Model(4, self.cs, nproc, extra_graphs=self.tg)
        self.p3 = None
    def parts(self, x, t, pole=False):
        """A (6) of: kappa insertions, c1, c2, c3 graphs"""
        A, V, G, by = self.m.A(x, t, parts=True, pole=pole)
        en = kin.env(x, t, self.cs['mpi2'])
        out = {}
        for key in ('c1', 'c2', 'c3'):
            Vk = sum(v for n, v in by.items() if n.startswith(f't_{key}_'))
            out[key] = np.linalg.solve(G, Vk)
        out['kappa'] = A - out['c1'] - out['c2'] - out['c3']
        out['total'] = A
        return out
