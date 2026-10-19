"""Figure: total Compton cross section sigma(gamma p -> gamma p) to 450 MeV, from sigma_tot.json (sigma_tot.py).
Since 2026-10-17 (user: Appendix F dropped, the prediction kept in the main text) the figure has the upper panel only;
the two-panel version with the cumulative Baldin integral is kept as work/backups/plot_sigma_tot_two_panels.py.  Single column, Computer Modern via the shared style.
Main curve: NNLO prediction (solid, dark grey) with the band of the O(p^4) LECs of natural size
(delta alpha_E1 = delta beta_M1 = +-0.7); NNLO fit dashed (blue); Born + pi0 pole; Thomson.  Legend entries and
colours as in Figs. cs and sigma3 of the paper (work/fits/scripts/s06_figures.py: TH_LAB, C_NNLO, C_F1, C_BAND).
The page size is pinned to that of the previous version (223.281 x 243.343 pt, tight-cropped 3.404 x 3.9 in):
the figure is included in the paper without scaling.
Run:  python3 plot_sigma_tot.py   (any cwd)
Output: paper/figs/fig_sigma_tot.pdf (+ fig_sigma_tot.png preview here)"""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import os, json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter, MultipleLocator
plt.style.use((_WORK + '/figstyle.mplstyle'))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = (_ROOT + '/paper/figs/fig_sigma_tot.pdf')
C_NNLO, C_F1, C_BAND, C_BORN, C_THOM, C_NB = '#444444', '#1f4e9c', '#9fb6dd', '#c0392b', '#e08e0b', '#2e8b57'
TH_LAB = ['NNLO B$\\chi$PT', 'B$\\chi$PT fit']      # band meaning in the caption (user, 2026-10-17)
R = json.load(open(os.path.join(HERE, 'sigma_tot.json'))); E = np.array(R['E'])
s, da, db, sf, sb, sth = (np.array(R[k]) for k in ('NNLO', 'NNLO_lec_da', 'NNLO_lec_db', 'NNLO_fit', 'BornPi0', 'Thomson'))
ds = np.hypot(da, db)


def cum(sig, Emin):
    out = np.full_like(E, np.nan); m = E >= Emin
    nu = E[m] / 197.327; f = sig[m] * 1e-4 / nu**2
    c = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(nu))]) / (2 * np.pi**2) * 1e4
    out[m] = c; return out


W, H = 223.281 / 72, 170.0 / 72                         # single panel [in]
fig, a1 = plt.subplots(1, 1, figsize=(W, H), layout='constrained')
fig.get_layout_engine().set(h_pad=0.01, w_pad=0.01)
# ---- upper panel: total cross section
a1.fill_between(E, s - ds, s + ds, color=C_BAND, alpha=0.6, lw=0, zorder=1)
a1.plot(E, sf, '--', color=C_F1, lw=1.0, zorder=2)
a1.plot(E, s, '-', color=C_NNLO, lw=1.2, zorder=3)
a1.plot(E, sb, ls='-.', color=C_BORN, lw=1.0, zorder=2)
a1.plot(E, sth, ls=':', color=C_THOM, lw=1.1, zorder=2)
a1.set_yscale('log'); a1.set_ylim(0.1, 5); a1.set_ylabel(r'$\sigma(\gamma p\to\gamma p)$ [$\mu$b]')
a1.set_xlabel(r'$E_\gamma$ [MeV]'); a1.set_xlim(0, 450)
a1.yaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1, 2, 5]))
a1.yaxis.set_major_formatter(FixedFormatter(['0.1', '0.2', '0.5', '1', '2', '5']))
a1.yaxis.set_minor_formatter(NullFormatter())
a1.axvline(144.7, color='0.5', lw=0.5, ls=':'); a1.axvline(309, color='0.5', lw=0.5, ls=':')
a1.text(147, 0.43, r'$\pi^0p$', fontsize=8, color='0.35'); a1.text(312, 0.115, r'$2\pi$', fontsize=8, color='0.35')
h1 = [(Patch(facecolor=C_BAND, alpha=0.6, lw=0), Line2D([], [], color=C_NNLO, lw=1.2)),
      Line2D([], [], color=C_F1, ls='--', lw=1.0),
      Line2D([], [], color=C_BORN, ls='-.', lw=1.0),
      Line2D([], [], color=C_THOM, ls=':', lw=1.1)]
fig.legend(h1, TH_LAB + [r'Born + $\pi^0$', 'Thomson'], loc='outside upper center', ncol=2,
           handletextpad=0.4, columnspacing=1.0)
with plt.rc_context({'savefig.bbox': 'standard'}):      # exact page size (no cropping)
    fig.savefig(OUT)
    fig.savefig(os.path.join(HERE, 'fig_sigma_tot.png'), dpi=300)
print('saved', OUT)
