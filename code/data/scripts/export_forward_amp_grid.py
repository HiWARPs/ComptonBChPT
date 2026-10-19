#!/usr/bin/env python3
"""Forward amplitudes f(nu), g(nu) of the ForwardRCS 2026 update (github.com/HiWARPs/ForwardRCS; a local
clone given by the environment variable FORWARDRCS_DIR, only READ here) ABOVE the pion-production threshold, as data of the
main BChPT fit (author request 2026-10-19: forward sum rules up to 360 MeV, at the level of the amplitudes): Re f, Im f,
Re g, Im g at the lab energies E = 160, 165, ..., 360 MeV with
  - the complete exp+fit covariance between all energies and between Re and Im (Jacobians d f/d p of the unpolarized fit,
    d g/d p of the polarized fit, times their parameter covariances; f and g share no parameter), as in
    pcompton/tools/export_forward.py for Re f, Re g below threshold;
  - the separate systematics as signed modes on the same grid, read from the tables of the update (copies in
    work/data/raw/ForwardRCS2026/): data selection of f (unpol_selection.txt, dRe f, dIm f), PWA model of g below 0.2 GeV
    (pol_g_pwa_model.txt, dRe g, dIm g), and the one-sided elastic-channel shift of f (unpol_missing.txt, dRe f, dIm f).
The diagonal of the covariance is checked against the errors of unpol_f.txt / pol_g.txt.
Run with the interpreter of the ForwardRCS environment (its requirements.txt):
  FORWARDRCS_DIR=<clone of HiWARPs/ForwardRCS> <its python> -I work/data/scripts/export_forward_amp_grid.py
-> work/data/raw/ForwardRCS2026/forward_amp_grid.json (units mcb GeV = mub GeV)."""
import os, sys, io, json, contextlib, warnings
HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, '..', 'raw', 'ForwardRCS2026')
PC = os.environ['FORWARDRCS_DIR']
E_MEV = [160.0 + 5 * i for i in range(41)]            # 5-MeV grid; the fit uses every 10 MeV (variants: 5, 20)
os.chdir(PC)                                        # the code reads data/... with relative paths
sys.path.insert(0, PC)
warnings.filterwarnings('ignore')
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    import unpol2026, pol2026
    wt, wp = unpol2026.WholeRegion2026(), pol2026.WholeRegionPol2026('bw')

k = np.array(E_MEV) / 1000.
f = np.array([complex(wt.f(x)) for x in k])
g = np.array([complex(wp.g(x)) for x in k])
Jf = np.array([[complex(wt.df(i, x)) for i in range(len(wt.get_p()))] for x in k])
Jg = np.array([[complex(wp.dg(i, x)) for i in range(len(wp.get_p()))] for x in k])
Jf2, Jg2 = np.vstack([Jf.real, Jf.imag]), np.vstack([Jg.real, Jg.imag])     # rows: Re(E1..En), Im(E1..En)
Cf, Cg = Jf2 @ wt.get_pcov() @ Jf2.T, Jg2 @ wp.get_pcov() @ Jg2.T


def table(name, cols):
    d = {}
    for line in open(os.path.join(RAW, name)):
        p = line.split()
        if p and p[0][0].isdigit():
            d[round(1000 * float(p[0]), 3)] = [float(p[c]) for c in cols]
    return np.array([d[round(e, 3)] for e in E_MEV])


ef, eg = table('unpol_f.txt', (1, 2, 3, 4)), table('pol_g.txt', (1, 2, 3, 4))
for name, val, tab, C in (('f', f, ef, Cf), ('g', g, eg, Cg)):
    assert np.allclose(val.real, tab[:, 0], atol=1e-6) and np.allclose(val.imag, tab[:, 2], atol=1e-6), name + ' values'
    err = np.sqrt(np.diag(C))
    assert np.allclose(err, np.r_[tab[:, 1], tab[:, 3]], rtol=2e-3, atol=1e-7), (name + ' errors', err, tab[:, [1, 3]])
sel, pwa, miss = table('unpol_selection.txt', (1, 2)), table('pol_g_pwa_model.txt', (1, 2)), table('unpol_missing.txt', (1, 2))
out = dict(info=__doc__.split('\n')[0] + ' Rows of cov_f, cov_g and the modes: Re(E_1..E_n), Im(E_1..E_n).',
           E=E_MEV, re_f=f.real.tolist(), im_f=f.imag.tolist(), re_g=g.real.tolist(), im_g=g.imag.tolist(),
           cov_f=Cf.tolist(), cov_g=Cg.tolist(), sel_f=np.r_[sel[:, 0], sel[:, 1]].tolist(),
           pwa_g=np.r_[pwa[:, 0], pwa[:, 1]].tolist(), miss_f=np.r_[miss[:, 0], miss[:, 1]].tolist(),
           rank=dict(f=int(np.linalg.matrix_rank(Cf)), g=int(np.linalg.matrix_rank(Cg))))
json.dump(out, open(os.path.join(RAW, 'forward_amp_grid.json'), 'w'), indent=1)
print('E  Re f  Im f  Re g  Im g  (errors)')
for i, e in enumerate(E_MEV):
    print('%5.0f  %8.4f(%6.4f) %8.4f(%6.4f)  %8.4f(%6.4f) %8.4f(%6.4f)' % (e, f[i].real, np.sqrt(Cf[i, i]), f[i].imag,
          np.sqrt(Cf[i + len(E_MEV), i + len(E_MEV)]), g[i].real, np.sqrt(Cg[i, i]), g[i].imag, np.sqrt(Cg[i + len(E_MEV), i + len(E_MEV)])))
print('rank of cov_f, cov_g:', out['rank'], 'of', 2 * len(E_MEV))
