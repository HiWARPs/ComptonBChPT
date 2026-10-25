import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, statpol, kin
cs = ev.consts()
m = assemble.Model(4, cs)
nu = 0.05; f = 2.0; t = -f*nu**2; x = nu - t/4
A, V, G, by = m.A(x, t, parts=True)
en = kin.env(x, t, cs['mpi2'])
for n, v in sorted(by.items(), key=lambda z: -np.max(abs(z[1])))[:6]:
    c = m.cmp[n]
    args = [fn(x, t, 0, 0) for fn in c.inv]
    print(n, f'{np.max(abs(v)):.3e}', c.kind, np.round(args, 6), c.ms)
print('x,t,s,u', x, t, 1+2*x, 2-(1+2*x)-t)
