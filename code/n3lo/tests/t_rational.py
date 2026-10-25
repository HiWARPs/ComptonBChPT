import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, statpol
cs = ev.consts()
m = assemble.Model(4, cs)
orig = ev.Compiled.tensors
def no_rat(self, en):
    T, U = orig(self, en)
    class Z(dict):
        def __missing__(s, k): return 0.0
    return T, Z()
for lab in ('with rational terms', 'without rational terms'):
    if lab.startswith('without'): ev.Compiled.tensors = no_rat
    p, dp = statpol.extract(lambda x, t: m.A(x, t))
    print(f'{lab:24s}', np.round(p, 4), ' g0,gpi', np.round(statpol.g0pi(p), 4))
