"""Figure: total Compton cross section sigma(gamma p -> gamma p) to 450 MeV and its cumulative contribution
to the Baldin integral, from sigma_tot.json (sigma_tot.py).  Single column, Computer Modern via the shared style.
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
TH_LAB = ['NNLO B$\\chi$PT, $\\delta\\alpha_{E1},\\delta\\beta_{M1}=\\pm0.7$', 'NNLO fit']
R = json.load(open(os.path.join(HERE, 'sigma_tot.json'))); E = np.array(R['E'])
s, da, db, sf, sb, sth = (np.array(R[k]) for k in ('NNLO', 'NNLO_lec_da', 'NNLO_lec_db', 'NNLO_fit', 'BornPi0', 'Thomson'))
ds = np.hypot(da, db)


def cum(sig, Emin):
    out = np.full_like(E, np.nan); m = E >= Emin
    nu = E[m] / 197.327; f = sig[m] * 1e-4 / nu**2
    c = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(nu))]) / (2 * np.pi**2) * 1e4
    out[m] = c; return out


W, H = 223.281 / 72, 243.343 / 72                       # previous page size [in]
fig, (a1, a2) = plt.subplots(2, 1, figsize=(W, H), sharex=True, layout='constrained',
                             gridspec_kw=dict(height_ratios=[1.35, 1]))
fig.get_layout_engine().set(hspace=0.03, h_pad=0.01, w_pad=0.01)
# ---- upper panel: total cross section
a1.fill_between(E, s - ds, s + ds, color=C_BAND, alpha=0.6, lw=0, zorder=1)
a1.plot(E, sf, '--', color=C_F1, lw=1.0, zorder=2)
a1.plot(E, s, '-', color=C_NNLO, lw=1.2, zorder=3)
a1.plot(E, sb, ls='-.', color=C_BORN, lw=1.0, zorder=2)
a1.plot(E, sth, ls=':', color=C_THOM, lw=1.1, zorder=2)
a1.set_yscale('log'); a1.set_ylim(0.1, 5); a1.set_ylabel(r'$\sigma(\gamma p\to\gamma p)$ [$\mu$b]')
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
# ---- lower panel: cumulative Baldin integral
c150, c50 = cum(s, 150), cum(s, 50)
e150 = np.hypot(cum(da, 150), cum(db, 150)); e50 = np.hypot(cum(da, 50), cum(db, 50))
a2.fill_between(E, c150 - e150, c150 + e150, color=C_BAND, alpha=0.6, lw=0, zorder=1)
a2.fill_between(E, c50 - e50, c50 + e50, color=C_BAND, alpha=0.6, lw=0, zorder=1)
a2.plot(E, cum(sf, 150), '--', color=C_F1, lw=1.0, zorder=2, label='NNLO fit, from 150 MeV')
a2.plot(E, c150, '-', color=C_NNLO, lw=1.2, zorder=3, label='from 150 MeV')
a2.plot(E, c50, ls=':', color=C_NNLO, lw=1.2, zorder=3, label='from 50 MeV')
a2.plot(E, cum(s - sb, 150), ls='-.', color=C_NB, lw=1.0, zorder=2, label='non-Born, from 150 MeV')
a2.set_ylabel(r'$B_{\gamma p}(E_\gamma)$ [$10^{-4}$\,fm$^3$]')
a2.set_xlabel(r'$E_\gamma$ [MeV]'); a2.set_xlim(0, 450); a2.set_ylim(0, 0.085)
a2.yaxis.set_major_locator(MultipleLocator(0.02)); a2.yaxis.set_minor_locator(MultipleLocator(0.01))
a2.axvline(144.7, color='0.5', lw=0.5, ls=':'); a2.axvline(309, color='0.5', lw=0.5, ls=':')
hh, ll = a2.get_legend_handles_labels()
order = ['from 150 MeV', 'from 50 MeV', 'NNLO fit, from 150 MeV', 'non-Born, from 150 MeV']
hh = [hh[ll.index(k)] for k in order]
hh[0] = (Patch(facecolor=C_BAND, alpha=0.6, lw=0), hh[0])
a2.legend(hh, order, loc='upper left', handletextpad=0.4)
fig.align_ylabels([a1, a2])
with plt.rc_context({'savefig.bbox': 'standard'}):      # exact page size (no cropping)
    fig.savefig(OUT)
    fig.savefig(os.path.join(HERE, 'fig_sigma_tot.png'), dpi=300)
print('saved', OUT)
