"""Figure: total Compton cross section sigma(gamma p -> gamma p) to 450 MeV and its cumulative contribution
to the Baldin integral.  Single column (246 pt), Computer Modern via the shared style."""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use((_WORK + '/figstyle.mplstyle'))
R = json.load(open('sigma_tot.json')); E = np.array(R['E'])
s, ds, sb, sth = (np.array(R[k]) for k in ('NNLO', 'NNLO_err', 'BornPi0', 'Thomson'))
def cum(sig, Emin):
    out = np.full_like(E, np.nan); m = E >= Emin
    nu = E[m] / 197.327; f = sig[m] * 1e-4 / nu**2
    c = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(nu))]) / (2 * np.pi**2) * 1e4
    out[m] = c; return out
fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.404, 3.9), sharex=True, gridspec_kw=dict(height_ratios=[1.35, 1], hspace=0.08))
a1.fill_between(E, s - ds, s + ds, color='#c0392b', alpha=0.25, lw=0)
a1.plot(E, s, color='#c0392b', label=r'B$\chi$PT fit')
a1.plot(E, sb, color='#1f4e9c', ls='--', label=r'Born + $\pi^0$')
a1.plot(E, sth, color='#444444', ls=':', label='Thomson')
a1.set_yscale('log'); a1.set_ylim(0.1, 5); a1.set_ylabel(r'$\sigma(\gamma p\to\gamma p)$ [$\mu$b]')
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
a1.yaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1, 2, 5]))
a1.yaxis.set_major_formatter(FixedFormatter(['0.1', '0.2', '0.5', '1', '2', '5']))
a1.yaxis.set_minor_formatter(NullFormatter())
a1.axvline(144.7, color='0.5', lw=0.5, ls=':'); a1.axvline(309, color='0.5', lw=0.5, ls=':')
a1.text(147, 0.43, r'$\pi^0p$', fontsize=8, color='0.35'); a1.text(312, 0.115, r'$2\pi$', fontsize=8, color='0.35')
a1.legend(loc='upper left')
for Emin, ls, lab in ((150, '-', r'from 150 MeV'), (50, '--', r'from 50 MeV')):
    a2.plot(E, cum(s, Emin), color='#c0392b', ls=ls, label=lab)
a2.plot(E, cum(s - sb, 150), color='#2e8b57', ls='-.', label=r'non-Born, from 150 MeV')
a2.set_ylabel(r'$B_{\gamma p}(E_\gamma)$ [$10^{-4}$\,fm$^3$]')
a2.set_xlabel(r'$E_\gamma$ [MeV]'); a2.set_xlim(0, 450); a2.set_ylim(0, 0.085)
a2.axvline(144.7, color='0.5', lw=0.5, ls=':'); a2.axvline(309, color='0.5', lw=0.5, ls=':')
a2.legend(loc='upper left')
fig.align_ylabels([a1, a2])
fig.savefig((_ROOT + '/paper/figs/fig_sigma_tot.pdf'))
print('saved')
