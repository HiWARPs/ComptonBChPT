import sys, os, numpy as np, subprocess, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import n3lo
pts = [(0.05, -0.008), (0.1, -0.05), (0.2, -0.1), (0.3, -0.3), (0.15, -0.002), (0.22, -0.18), (0.4, -0.5)]
inp = '\n'.join(f'{x} {t} {p}' for x, t in pts for p in (0, 1, 2, 3, 4))
t0 = time.time()
out = subprocess.run([os.path.join(os.path.dirname(__file__), '..', 'fortran', 'test_p4.exe')], input=inp, capture_output=True, text=True).stdout
tf = time.time() - t0
rows = [l for l in out.split('\n') if len(l.split()) == 12]
hdr = [l for l in out.split('\n') if l.startswith('#')]
print('\n'.join(hdr)); print(f'Fortran: {len(rows)} evaluations in {tf:.2f} s')
m = n3lo.P4()
print('python vc:', [np.round(v, 8) for v in m.m.vc.values()], 'L0,L1', m.m.L0, m.m.L1)
k = 0
for x, t in pts:
    r = m.parts(x, t)
    for p, key in zip((0, 1, 2, 3, 4), ('total', 'kappa', 'c1', 'c2', 'c3')):
        v = np.array([float(z) for z in rows[k].split()]); Af = v[0::2] + 1j*v[1::2]; k += 1
        d = np.max(abs(Af - r[key])) / max(np.max(abs(r[key])), 1e-30)
        if key == 'total' or d > 1e-8: print(f'x={x} t={t} {key}: rel dev {d:.2e}')
