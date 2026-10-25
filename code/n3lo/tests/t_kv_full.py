import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, statpol
base = ev.consts()
for lab, kp, kn in (('kv', 0.5, -0.5), ('ks', 0.5, 0.5), ('phys', base['kp'], base['kn'])):
    cs = dict(base, kp=kp, kn=kn)
    m = assemble.Model(4, cs)
    p, dp = statpol.extract(lambda x, t: m.A(x, t))
    print(lab, ' '.join('%8.4f' % v for v in list(p) + list(statpol.g0pi(p))), ' spread %.1e' % dp.max())
