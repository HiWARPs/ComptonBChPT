"""make_fig.py -- paper figure: O(alpha) QED correction factor for gamma p -> gamma p (point-like proton).
Writes paper/figs/fig_radcorr.pdf and work/radcorr/results/fig_radcorr.pdf (+ .png preview)."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
ROOT = os.path.dirname(WORK)
sys.path.insert(0, HERE)
import radcorr as rc

plt.style.use(os.path.join(WORK, 'figstyle.mplstyle'))
W1 = 3.404
COL = ['#1f4e9c', '#c0392b', '#2e8b57', '#e08e0b', '#6a3d9a']

om = np.linspace(20, 160, 29)
thetas = [0, 45, 90, 135, 180]
fig, axs = plt.subplots(2, 1, figsize=(W1, 4.3), sharex=True, layout='constrained')

# (a) Delta E = 3 MeV, several angles; dashed: virtual part only (lambda = M), shown for 180 deg
ax = axs[0]
for c, th in zip(COL, thetas):
    y = np.array([100 * rc.delta(w, th, 3.0) for w in om])
    ax.plot(om, y, color=c, label=r'$\theta_{\rm lab}=%d^\circ$' % th)
ax.axhline(0, color='0.6', lw=0.5)
ax.set_ylabel(r'$\delta\;[\%]$')
ax.legend(loc='upper left', ncol=2, fontsize=7)
ax.text(0.97, 0.06, r'$\Delta E = 3$ MeV', transform=ax.transAxes, ha='right', va='bottom')

# (b) angular dependence at 100 MeV for Delta E = 1, 3, 10 MeV and 3% of omega
ax = axs[1]
th = np.linspace(0, 180, 37)
for c, (dE, lab, ls) in zip(COL, [(1.0, r'$\Delta E=1$ MeV', '-'), (3.0, r'$3$ MeV', '-'), (10.0, r'$10$ MeV', '-')]):
    y = np.array([100 * rc.delta(100.0, t, dE) for t in th])
    ax.plot(th, y, color=c, ls=ls, label=lab)
yv = np.array([100 * rc.delta_virtual_only(100.0, t, lam=rc.M_P) for t in th])
ax.plot(th, yv, color='0.35', ls='--', label=r'virtual only ($\lambda=M$)')
ax.axhline(0, color='0.6', lw=0.5)
ax.set_xlabel(r'$\theta_{\rm lab}$ [deg]')
ax.set_ylabel(r'$\delta\;[\%]$')
ax.set_xlim(0, 180); ax.set_xticks([0, 45, 90, 135, 180])
ax.legend(loc='lower left', fontsize=7)
ax.text(0.97, 0.94, r'$\omega=100$ MeV', transform=ax.transAxes, ha='right', va='top')
axs[0].set_xlim(20, 160)
axs[0].set_xlabel(r'$\omega$ [MeV]')
axs[0].xaxis.set_label_position('bottom'); axs[0].tick_params(labelbottom=True)
for a in axs:
    a.set_ylim(-0.04, 0.045)
os.makedirs(os.path.join(ROOT, 'paper', 'figs'), exist_ok=True)
os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
fig.savefig(os.path.join(ROOT, 'paper', 'figs', 'fig_radcorr.pdf'))
fig.savefig(os.path.join(HERE, 'results', 'fig_radcorr.pdf'))
fig.savefig(os.path.join(HERE, 'results', 'fig_radcorr.png'), dpi=200)
print('written paper/figs/fig_radcorr.pdf')
