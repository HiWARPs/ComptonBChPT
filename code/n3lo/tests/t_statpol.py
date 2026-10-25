import sys, os, numpy as np, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import statpol, assemble, evalamp as ev, kin
cs = ev.consts()
m3 = assemble.Model(3, cs)
p, dp = statpol.extract(lambda x, t: m3.A(x, t))
print('piN loops O(p^3), this code :', np.round(p, 4), 'g0,gpi', np.round(statpol.g0pi(p), 4), ' spread', dp.max())
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'code_v2', 'driver'))
from analytic_pols import pols
a = pols()['piN loops']; print('piN loops analytic LMP15 App C:', np.round(a, 4))
