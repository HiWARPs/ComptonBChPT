"""Step 11: Fig. predS2abs -- the ABSOLUTE double-polarization asymmetries Sigma_2x, Sigma_2z of proton Compton scattering
below the pion-production threshold (Sec. IV G of the paper): the parameter-free NNLO BChPT prediction with its
uncertainty band, the NNLO fit, the partial-wave analysis (PWA) of this work and the empirical forward points.  Replaces
Table predS2 (Appendix) of the paper; the difference figure fig_pred_S2_diff.pdf (s10) is unchanged.

Outputs
  paper/figs/fig_pred_S2_abs.pdf         Fig. predS2abs (double column, 510 pt): Sigma_2x (upper row), Sigma_2z (lower row)
                                         versus theta_lab at E_lab = 100, 120, 140 MeV (copy: work/fits/figs/fig_pred_S2_abs.pdf)
  work/fits/figs/fig_pred_S2_abs.png     preview (200 dpi) for inspection
  work/fits/results/pred_S2_abs_checks.json   maximal differences PWA - NNLO, fit - NNLO, forward points, y ranges

Curves (per panel; all from files, nothing recomputed here)
  (1) NNLO prediction (solid '#444444') with the band (fill '#9fb6dd', alpha 0.6) = truncation (+) LEC uncertainty in
      quadrature, and (2) the NNLO fit (dashed '#1f4e9c'): columns NNLO, sig, F2 of work/fits/results/pred_S2_curves.csv,
      written by s10_fig_pred_S2.py (grid theta_lab = 0, 5, ..., 180 deg; the same numbers as Table predS2);
  (3) PWA fit (dotted olive '#6b8e23'): work/pwa/results/pwa_curves_S2.csv, written by work/pwa/scenarios/pwa_curves_S2.py
      (the main PWA fit of the paper, residuals of the bin nearest to the panel energy: 99, 117, 143 MeV; its statistical
      fit error is in the file but NOT drawn: 0.01-0.03 at 90 deg);
  (4) forward points (open black stars at theta = 0, Sigma_2z panels only): the empirical Sigma_2z(0) of Gryniuk-Hagelstein-
      Pascalutsa (database rows 'GHP_fwd_S2Z'), cubic spline through the 12 tabulated values, error interpolated linearly
      (exactly forward_S2Z() of s10).  Sigma_2x vanishes at theta = 0.
Conventions: A2 (driver codes S2X, S2Z, LAB frame): circularly polarized beam, Sigma_2x = (sigma_R - sigma_L)/(sigma_R +
sigma_L) with the target spin in the scattering plane (perpendicular to the beam), Sigma_2z with the target spin along the
beam; the PWA values use compton.S2_from_phis, validated against the driver (work/pwa/validation/validate_S2.py).
Run by hand after s10 and pwa_curves_S2:  .venv/bin/python scripts/s11_fig_pred_S2_abs.py   (~5 s)"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.interpolate import CubicSpline
from common import FITS, RES, FIGS, STYLE, load_db  # noqa

plt.style.use(STYLE)
ROOT = os.path.dirname(os.path.dirname(FITS))
PAPER_FIGS = os.path.join(ROOT, 'paper', 'figs')
PWA_CSV = os.path.join(ROOT, 'work', 'pwa', 'results', 'pwa_curves_S2.csv')
os.makedirs(FIGS, exist_ok=True)

C_NNLO, C_F2, C_BAND, C_PWA = '#444444', '#1f4e9c', '#9fb6dd', '#6b8e23'   # prediction, fit, band, PWA
W2 = 510 / 72.27                     # double-column width [in] (revtex4-2 prd), included without scaling
ENERGIES = [100.0, 120.0, 140.0]
OBS = [('S2X', '$\\Sigma_{2x}$'), ('S2Z', '$\\Sigma_{2z}$')]
YLIM = {'S2X': (-0.9, 0.1), 'S2Z': (-1.08, 1.08)}
YTICKS = {'S2X': [-0.8, -0.6, -0.4, -0.2, 0.0], 'S2Z': [-1.0, -0.5, 0.0, 0.5, 1.0]}

# ---------------------------------------------------------------- input curves
nn = pd.read_csv(os.path.join(RES, 'pred_S2_curves.csv'))
pw = pd.read_csv(PWA_CSV)
assert set(pw.obs) == {'S2X', 'S2Z'} and set(np.unique(pw.E)) == set(ENERGIES), 'PWA curves: unexpected content'


def nnlo(o, E):
    s = nn[(nn.dobs == o) & np.isclose(nn.E, E)].sort_values('theta')
    assert len(s) == 37 and np.isclose(s.theta.iloc[0], 0) and np.isclose(s.theta.iloc[-1], 180)
    return s.theta.values.astype(float), s.NNLO.values, s.sig.values, s.F2.values


def pwa(o, E):
    s = pw[(pw.obs == o) & np.isclose(pw.E, E)].sort_values('theta_lab')
    assert len(s) == 73 and np.isclose(s.theta_lab.iloc[0], 0) and np.isclose(s.theta_lab.iloc[-1], 180)
    return s.theta_lab.values, s.value.values, s.err.values


# ---------------------------------------------------------------- forward points Sigma_2z(0) of GHP (as s10)
df = load_db()
fw = df[df.dataset == 'GHP_fwd_S2Z'].sort_values('E')
assert len(fw) == 12 and np.all(fw.err_syspt == 0) and np.allclose(fw.err_pt, fw.err_stat)


def forward_S2Z(E):
    return float(CubicSpline(fw.E.values, fw.value.values)(E)), float(np.interp(E, fw.E.values, fw.err_pt.values))


# ---------------------------------------------------------------- figure
chk = {'max_abs_PWA_minus_NNLO': {}, 'max_abs_fit_minus_NNLO': {}, 'band_max': {}, 'pwa_err_max': {}, 'forward': {}}
fig, axs = plt.subplots(2, 3, figsize=(W2, 4.35), sharex=True, layout='constrained')
h_fwd = None
for r, (o, sym) in enumerate(OBS):
    for c, E in enumerate(ENERGIES):
        ax = axs[r, c]
        th, y, sig, f2 = nnlo(o, E)
        thp, yp, ep = pwa(o, E)
        ax.axhline(0, color='0.6', lw=0.4, zorder=0)
        ax.fill_between(th, y - sig, y + sig, color=C_BAND, alpha=0.6, lw=0, zorder=1)
        ax.plot(th, y, '-', color=C_NNLO, lw=1.2, zorder=2)
        ax.plot(th, f2, '--', color=C_F2, lw=1.0, zorder=3)
        ax.plot(thp, yp, ':', color=C_PWA, lw=1.2, zorder=4)
        if o == 'S2Z':
            v, e = forward_S2Z(E)
            h_fwd = ax.errorbar([0.0], [v], [e], fmt='*', color='k', mfc='white', ms=6.5, lw=0.7, zorder=6)
            chk['forward']['%g' % E] = dict(value=v, err=e, NNLO=float(y[0]), sig=float(sig[0]), PWA=float(yp[0]), pwa_err=float(ep[0]))
        if o == 'S2X':   # curves vanish at the ends and are high (near zero) at small angles: the lower left corner is free
            ax.text(0.04, 0.06, '%g MeV' % E, transform=ax.transAxes, va='bottom', ha='left')
        else:            # Sigma_2z is negative (-0.45 ... -0.9) at small angles: the upper left corner is free
            ax.text(0.04, 0.93, '%g MeV' % E, transform=ax.transAxes, va='top', ha='left')
        ax.set_xlim(-6, 180); ax.set_xticks([0, 45, 90, 135, 180])
        ax.set_ylim(*YLIM[o]); ax.set_yticks(YTICKS[o])
        if c == 0:
            ax.set_ylabel(sym)
        if r == 1:
            ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')
        # bookkeeping: the numbers quoted in the text
        m = np.isin(thp, th)
        d = yp[m] - y
        i = int(np.argmax(np.abs(d)))
        chk['max_abs_PWA_minus_NNLO']['%s_%g' % (o, E)] = dict(value=float(abs(d[i])), theta=float(th[i]), PWA=float(yp[m][i]),
                                                               NNLO=float(y[i]), band=float(sig[i]), pwa_err=float(ep[m][i]))
        j = int(np.argmax(np.abs(f2 - y)))
        chk['max_abs_fit_minus_NNLO']['%s_%g' % (o, E)] = dict(value=float(abs(f2[j] - y[j])), theta=float(th[j]))
        chk['band_max']['%s_%g' % (o, E)] = float(sig.max())
        chk['pwa_err_max']['%s_%g' % (o, E)] = float(ep.max())
        assert (y - sig).min() > YLIM[o][0] and (y + sig).max() < YLIM[o][1] and yp.min() > YLIM[o][0] and yp.max() < YLIM[o][1]
handles = [(Patch(facecolor=C_BAND, alpha=0.6, lw=0), Line2D([], [], color=C_NNLO, lw=1.2)),
           Line2D([], [], color=C_F2, ls='--', lw=1.0), Line2D([], [], color=C_PWA, ls=':', lw=1.2), h_fwd]
labels = ['NNLO B$\\chi$PT', 'NNLO fit', 'PWA fit', 'forward point']
fig.legend(handles, labels, loc='outside lower center', ncol=4, columnspacing=1.6, handlelength=2.2)
with plt.rc_context({'savefig.bbox': 'standard'}):
    for f in (os.path.join(PAPER_FIGS, 'fig_pred_S2_abs.pdf'), os.path.join(FIGS, 'fig_pred_S2_abs.pdf')):
        fig.savefig(f); print('wrote', f)
    fig.savefig(os.path.join(FIGS, 'fig_pred_S2_abs.png'), dpi=200); print('wrote', os.path.join(FIGS, 'fig_pred_S2_abs.png'))
plt.close(fig)
chk['ylim'] = YLIM
json.dump(chk, open(os.path.join(RES, 'pred_S2_abs_checks.json'), 'w'), indent=1)
print('max |PWA - NNLO| (at theta), |fit - NNLO|, band max, PWA fit error max:')
for o, _ in OBS:
    for E in ENERGIES:
        k = '%s_%g' % (o, E); a = chk['max_abs_PWA_minus_NNLO'][k]; b = chk['max_abs_fit_minus_NNLO'][k]
        print('  %s %3.0f MeV: PWA-NNLO %.4f at %3.0f deg (PWA %+.4f, NNLO %+.4f, band %.4f, PWA err %.4f) | fit-NNLO %.4f at %3.0f deg | band max %.4f | PWA err max %.4f' % (
            o, E, a['value'], a['theta'], a['PWA'], a['NNLO'], a['band'], a['pwa_err'], b['value'], b['theta'], chk['band_max'][k], chk['pwa_err_max'][k]))
print('forward points:', json.dumps(chk['forward'], indent=1))
