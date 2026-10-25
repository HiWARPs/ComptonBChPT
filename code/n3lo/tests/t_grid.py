import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import assemble, evalamp as ev, statpol
cs = ev.consts()
m = assemble.Model(4, cs)
for f in (0.5, 1.0, 2.0):
    for nu in np.linspace(0.02, 0.12, 11):
        t = -f*nu**2; x = nu - t/4
        A = m.A(x, t); L = statpol.lvov(A, nu, t).real
        print(f'f={f} nu={nu:.3f}  L1..L6 = {np.array2string(L[:6], precision=3, max_line_width=200)}')
