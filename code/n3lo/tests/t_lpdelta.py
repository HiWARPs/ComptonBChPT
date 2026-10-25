"""t_lpdelta.py -- the NNLO pi Delta loops of the code (CorrIsoCovOneDeltaLoop: photons on the pions, gamma pi N Delta
contacts, 1PR graphs with the 'corrected isospin' factor 4/3 for all graphs) reproduced by the Delta generator."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import delta, dassemble, devalamp as dv, formtools as ft, kin, dgraphs
if len(sys.argv) > 1 and sys.argv[1] == '4d': dgraphs.RS4D = True
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
S8 = [ft.compile_proj(ft.parse_fortran(open(os.path.join(HERE, '..', 'out', f'gram{k+1}.log')).read())) for k in range(8)]
cs = dv.dconsts()
gs = [g for g in delta.all_graphs(maxgND=0) if g.order == Fr(7, 2)]
if dgraphs.RS4D:
    for g in gs: g.name += '_4d'
print(len(gs), 'graphs of O(p^7/2):', sorted({g.kind for g in gs}))
FUDGE = {'1PI': 1.0, 'R2': -2.0, 'R1': 2 / 3, 'vertex': -2.0, 'self': 2 / 3}   # facv, face of the code
m = dassemble.DModel(gs, cs, iso_override=FUDGE, tag='lp' + ('4d' if dgraphs.RS4D else ''))
gafac = cs['hA'] ** 2 / (4 * cs['MD'] * np.pi * cs['F']) ** 2
pts = [(0.05, -0.008), (0.1, -0.05), (0.2, -0.1), (0.3, -0.3), (0.15, -0.002)]
inp = '\n'.join(f'{x} {t} {cs["mpi2"]} {gafac}' for x, t in pts)
out = subprocess.run([os.path.join(HERE, 'covd.exe')], input=inp, capture_output=True, text=True).stdout
rows = [l for l in out.split('\n') if len(l.split()) == 16]
for (x, t), r in zip(pts, rows):
    v = np.array([float(z) for z in r.split()]); Ac = v[0::2] + 1j * v[1::2]
    en = dv.env(x, t, cs); en.update(T={(0, 0, 0, 0): 1.0}, U={})
    G8 = np.array([[eval(S8[k][j], en) for k in range(8)] for j in range(6)])
    Vc = G8 @ Ac
    V, by = m.V(x, t)
    print(f'x={x} t={t}: rel dev {np.max(abs(V - Vc)) / np.max(abs(Vc)):.2e}')
    print('   mine', np.round(V, 9)); print('   code', np.round(Vc, 9))

print('--- by class (LP sections, s + u halves)')
CL = {1: lambda g: g.kind == '1PI' and g.cfg['pi'] and len(g.cfg['pi']) == 2,
      2: lambda g: g.kind == '1PI' and g.cfg['V1'] and g.cfg['V2'],
      3: lambda g: g.kind == '1PI' and len(g.cfg['pi']) == 1,
      4: lambda g: g.kind == 'R2', 5: lambda g: g.kind == 'R1', 6: lambda g: g.cfg['sg']}
NM = {1: 'Box', 2: 'KR', 3: 'KP', 4: 'VR', 5: 'SE', 6: 'TP'}
for (x, t) in pts[1:4]:
    inp = '\n'.join(f'{x} {t} {cs["mpi2"]} {gafac} {k}' for k in range(1, 7))
    out = subprocess.run([os.path.join(HERE, 'covdsec.exe')], input=inp, capture_output=True, text=True).stdout
    rows = [l for l in out.split('\n') if len(l.split()) == 16]
    en = dv.env(x, t, cs); en.update(T={(0, 0, 0, 0): 1.0}, U={})
    G8 = np.array([[eval(S8[k][j], en) for k in range(8)] for j in range(6)])
    for k, r in zip(range(1, 7), rows):
        v = np.array([float(z) for z in r.split()]); Vc = G8 @ (v[0::2] + 1j * v[1::2])
        names = {g.name for g in gs if CL[k](g)}
        V, _ = m.V(x, t, names=names)
        print(f'x={x} t={t} {NM[k]:3s} ({len(names)} graphs): rel dev {np.max(abs(V - Vc)) / max(np.max(abs(Vc)), 1e-30):.2e}'
              f'   mine {np.round(V.real, 6)}   code {np.round(Vc.real, 6)}')
