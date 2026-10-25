"""t_dstatpol.py -- static polarisabilities of the NNLO pi Delta loops: generator (corrected isospin) vs LMP2015 App. C
(analytic expressions of the code) and vs the code's own amplitude (covd.exe) through the same L'vov extraction."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import delta, dassemble, devalamp as dv, formtools as ft, evalamp as ev, statpol
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'code_v2', 'driver'))
from analytic_pols import pols
HERE = os.path.dirname(os.path.abspath(__file__))
S8 = [ft.compile_proj(ft.parse_fortran(open(os.path.join(HERE, '..', 'out', f'gram{k+1}.log')).read())) for k in range(8)]
cs = dv.dconsts()
gs = [g for g in delta.all_graphs(maxgND=0) if g.order == Fr(7, 2)]
FUDGE = {'1PI': 1.0, 'R2': -2.0, 'R1': 2 / 3, 'vertex': -2.0, 'self': 2 / 3}
m = dassemble.DModel(gs, cs, iso_override=FUDGE, tag='lp')
gafac = cs['hA'] ** 2 / (4 * cs['MD'] * np.pi * cs['F']) ** 2
def code_A(x, t):
    out = subprocess.run([os.path.join(HERE, 'covd.exe')], input=f'{x} {t} {cs["mpi2"]} {gafac}\n', capture_output=True, text=True).stdout
    r = [l for l in out.split('\n') if len(l.split()) == 16][0]
    v = np.array([float(z) for z in r.split()]); Ac = v[0::2] + 1j * v[1::2]
    en = dv.env(x, t, cs); en.update(T={(0, 0, 0, 0): 1.0}, U={})
    G8 = np.array([[eval(S8[k][j], en) for k in range(8)] for j in range(6)])
    return np.linalg.solve(ev.gram(en), G8 @ Ac)
p, dp = statpol.extract(lambda x, t: m.A(x, t))
print('generator, corrected isospin :', np.round(p, 4), ' spread', np.round(dp.max(), 5))
pc, dpc = statpol.extract(code_A)
print('code (CorrIso), same route   :', np.round(pc, 4), ' spread', np.round(dpc.max(), 5))
print('LMP2015 App. C (analytic)    :', np.round(pols()['piDelta loops'], 4))
