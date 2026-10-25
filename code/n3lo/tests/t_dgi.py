"""t_dgi.py -- gauge invariance (incoming polarisation -> q) of the Delta-loop graph sets."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'py'))
import delta, dassemble, devalamp as dv
from fractions import Fraction as Fr

cs = dv.dconsts()
sets = {}
allA = delta.all_graphs(maxgND=0)
sets['NNLO, corrected isospin (code)'] = ([g for g in allA if g.order == Fr(7, 2)],
                                          {'1PI': 1.0, 'R2': -2.0, 'R1': 2 / 3, 'vertex': -2.0, 'self': 2 / 3})
sets['NNLO graphs, true isospin'] = ([g for g in allA if g.order == Fr(7, 2)], None)
sets['all minimal couplings, true isospin'] = (allA, None)
which = sys.argv[1:] or list(sets)
for k in which:
    gs, fud = sets[k]
    m = dassemble.DModel(gs, cs, iso_override=fud, tag='gi' + str(abs(hash(k)) % 1000))
    m.gi_compile()
    for x, t in [(0.1, -0.05), (0.25, -0.2)]:
        W, big = m.GI(x, t)
        print(f'{k:40s} x={x} t={t}: max|sum| = {np.max(abs(W)):.2e}   largest single graph {big:.2e}')
