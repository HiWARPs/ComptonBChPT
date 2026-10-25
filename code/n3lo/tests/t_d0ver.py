import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lt'))
import ltpy
mp = (139.57/938.27)**2
worst = 0
for nu in np.linspace(0.01, 0.5, 25):
    for f in np.linspace(0.05, 3.9, 30):
        t = -f*nu**2; x = nu - t/4; s = 1+2*x; u = 2-s-t
        if t < -4*x**2/(1+2*x): continue
        args = (0., 1., 0., 1., u, s, mp, mp, 1., 1.)
        ltpy.setversion(0); ltpy.clear(); a, _ = ltpy.D(*args)
        ltpy.setversion(64); ltpy.clear(); b, _ = ltpy.D(*args)
        d = max(abs(a[k]-b[k])/(abs(b[k])+1e-30) for k in ('dd0','dd1','dd22','dd00','dd123','dd00001'))
        if d > 1e-6: print(f'nu={nu:.3f} f={f:.2f}: rel diff a/b = {d:.2e}  dd0 a={a["dd0"]:.6g} b={b["dd0"]:.6g}')
        worst = max(worst, d if d < 1e-3 else 0)
print('max agreement elsewhere', worst)
