import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, statpol
cs = ev.consts()
for order in (3, 4):
    m = assemble.Model(order, cs)
    print('order', order)
    for nu in (0.04, 0.02, 0.01):
        t = -nu**2; x = nu - t/4
        A, V, G, by = m.A(x, t, parts=True)
        L = statpol.lvov(A, nu, t).real
        print(f'  nu={nu}: A={np.array2string(A.real, precision=5)}   lvov*nu^0 = {np.array2string(L[:6], precision=3)}')
    if order == 4:
        nu = 0.02; t = -nu**2; x = nu - t/4
        A, V, G, by = m.A(x, t, parts=True)
        groups = {}
        for n, v in by.items():
            key = ('sub_' if n.endswith('_sub') else '') + m.cmp[n].g.__dict__.get('kind', 'sub')
            groups.setdefault(n.replace('_sub',''), np.zeros(6, complex))
            groups[n.replace('_sub','')] += v
        for n, v in sorted(groups.items(), key=lambda z: -np.max(abs(np.linalg.solve(G, z[1])))):
            a = np.linalg.solve(G, v)
            print(f'   {n:28s} {np.array2string(a.real, precision=4)}')
