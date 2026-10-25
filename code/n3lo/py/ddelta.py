"""ddelta.py -- Delta corrections beyond NNLO (Python check path).

Pieces (all renormalised as in dassemble.py, MS-bar at mu = M for the remaining UV poles):
  dA     : pi Delta loops with the complete minimal coupling (photon on the Delta charge, Delta-derivative contacts,
           true isospin factors) minus the NNLO code's pi Delta loops (photon on the pion only, 'corrected isospin');
           starts at O(p^4), gauge invariant only as a whole.
  Bloop  : gamma N Delta (gM, gE) vertex inside the pi N / pi Delta loop, the second photon minimally coupled
           everywhere (O(p^9/2); the partners of the O(p^4) kappa insertions).
  Bpole  : Delta-pole graphs with the one-loop gamma N -> Delta vertex (pi N loop, O(p^9/2)), MS-bar.
"""
import numpy as np, delta, dassemble, devalamp as dv
from fractions import Fraction as Fr

FUDGE = {'1PI': 1.0, 'R2': -2.0, 'R1': 2 / 3, 'vertex': -2.0, 'self': 2 / 3}   # facv, face of CorrIsoCovOneDeltaLoop

def gnd_on_leg(g):
    """True if the gamma N Delta vertex sits on an external leg (Delta tree propagator): the Delta-pole graphs"""
    return any(g.cpl.get((w, i)) in ('gM', 'gE') or str(g.cpl.get((w, i)))[:2] in ('gM', 'gE')
               for w in ('in', 'out') for i in range(len(g.cfg[w])))

class Delta45:
    def __init__(self, cs=None, nproc=8, which=('A', 'B')):
        self.cs = cs or dv.dconsts()
        self.m = {}
        if 'A' in which:
            allA = delta.all_graphs(maxgND=0)
            self.m['A_full'] = dassemble.DModel(allA, self.cs, nproc, tag='af')
            self.m['A_lp'] = dassemble.DModel([g for g in allA if g.order == Fr(7, 2)], self.cs, nproc,
                                              iso_override=FUDGE, tag='lp')
        if 'B' in which:
            allB = [g for g in delta.all_graphs(maxgND=1) if g.ngND == 1]
            self.m['Bloop'] = dassemble.DModel([g for g in allB if not gnd_on_leg(g)], self.cs, nproc, tag='bl')
            self.m['Bpole'] = dassemble.DModel([g for g in allB if gnd_on_leg(g)], self.cs, nproc, tag='bp')
    def parts(self, x, t, pole=False):
        out = {k: m.A(x, t, pole=pole) for k, m in self.m.items()}
        if 'A_full' in out: out['dA'] = out['A_full'] - out['A_lp']
        return out
