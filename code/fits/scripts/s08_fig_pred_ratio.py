"""Figure: the parameter-free NNLO prediction versus the data near the new-data energies (paper Fig. "predratio").

Ratio data/prediction versus the lab angle at the energies of Fig. "cs" (s06_figures.py, fig1a_cs_A2_HIGS): the HIgS
energy (unpolarised cross section at 81.3 MeV, polarised dsigma_par, dsigma_perp at 83.4 MeV), the five A2 2022
energies and 149.87 MeV; all cross-section data within +-2.5 MeV of the panel energy, Pugh 1957 omitted (as in Fig. cs).
Since 2026-10-27 (author): the figure is about the prediction alone; each data set is divided by the normalisation
the prediction prefers (scenarios/nnlo_prediction_table.json[set]['n_th'], Table pred: normalisation alone fitted with
the truncation covariance); the theory curves carry no normalisation.  The BChPT fit is not drawn (it is in Figs. cs and
pulls).  The prediction is the NNLO amplitude at P = PREF_NNLO (no free polarisability), evaluated per data point with
the driver ('db' run).
Drawn in every panel: the line at 1; the truncation band 1 +- sigma_th/prediction (theory.truncation via
truncation_for on a fine angular grid, same settings as the fit covariance: scheme B, sqrtQ, Lambda_b = 750 MeV,
cbar = rms); the PWA fit of this work (bare curve, without its normalisations, from work/pwa/results/pwa_curves.csv as in Fig. cs) divided by the prediction (dotted olive); and the
forward SR at theta = 0 (open star; the 'GHP_fwd_CS' rows interpolated to the panel energy as in Fig. cs, no
normalisation).  The curves of the HIgS panel are those of the unpolarised cross section at 81.3 MeV (the polarised
points are divided by their own predictions; at theta = 0 all three cross sections coincide).  The three HIgS
measurements and repeated points at the same angle (beam times) are displaced by a few degrees for visibility.
Writes work/fits/figs/fig_pred_ratio.pdf and paper/figs/fig_pred_ratio.pdf (final printed size, double column).
Not part of the make pipeline; run by hand:  .venv/bin/python scripts/s08_fig_pred_ratio.py
"""
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
from common import *  # noqa
import database as DB

plt.style.use(STYLE)
os.makedirs(FIGS, exist_ok=True)
PAPER_FIGS = os.path.join(os.path.dirname(os.path.dirname(FITS)), 'paper', 'figs')
W1, W2 = 246 / 72.27, 510 / 72.27          # final printed widths [in] (as in s06_figures.py)
C_BAND, C_PWA = '#9fb6dd', '#6b8e23'
YLIM = (0.75, 1.25)                         # common y-range (narrowed 2026-10-17); points outside are marked and reported
OWN_YLIM = {149.87: (0.75, 1.50)}           # own range: just above the pi0 threshold the forward SR is 1.33 x prediction
PANELS = [81.3, 92.2, 103.2, 113.6, 124.1, 135.3, 149.87]   # energies of Fig. cs
WINDOW, OMIT = 2.5, ('Pugh1957',)
HIGS = ['HIGS_2022_CIRC', 'HIGS_2022_LIN']

df = load_db()
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
FIT = FL['F2_d_all_thB']
NORM = {k: v['n_th'] for k, v in json.load(open(os.path.join(FITS, 'scenarios', 'nnlo_prediction_table.json'))).items()
        if isinstance(v, dict) and 'n_th' in v}
PWA = pd.read_csv(os.path.join(FITS, '..', 'pwa', 'results', 'pwa_curves.csv'))

# ---- data points: ratio to the prediction, each set at the normalisation the prediction prefers
T = outputs(df, 'NNLO')                     # cached 'db' run of all database rows
near = np.array([np.abs(df.E.values - E) <= WINDOW for E in PANELS])
sel = df[near.any(0) & df.dobs.isin(['CS', 'CSCM', 'SPAR', 'SPERP']) & ~df.dataset.isin(FWD + list(OMIT))
         & (df.theta_lab > 0)].copy()
sel['panel'] = [PANELS[int(np.argmin(np.abs(np.array(PANELS) - E)))] for E in sel.E]
assert set(sel.dataset) <= set(NORM), set(sel.dataset) - set(NORM)
sel['n'] = [NORM[ds] for ds in sel.dataset]
idx = T.subset(sel.pid.values)
sel['pred'] = T.predict_array(PREF_NNLO, idx)
assert np.allclose(sel['pred'].values, T.value_file[idx])
sel['ratio'] = sel.value / (sel.n * sel.pred)
sel['eratio'] = sel.err_pt / (sel.n * sel.pred)
# truncation sigma at the data points: identical to FIT['sig_th'] at the pids that are in the fit below threshold
_, tr = truncation_for(df, sel, tag='db')
sel['sig_th'] = tr['sigma']
fit_sig = dict(zip(FIT['pids'], FIT['sig_th']))
_infit = sel.pid.isin(fit_sig).values
assert not sel[_infit].dataset.isin(HIGS).any()
assert np.allclose(sel.sig_th.values[_infit], [fit_sig[p] for p in sel.pid[_infit]])

# ---- curves on a fine angular grid: CS at the panel energies
th = np.arange(0.0, 180.01, 2.5)


def grid_rows(spec):
    """as in s06_figures.py: list of (tag, E, theta, dobs, frame) -> DataFrame for outputs()/truncation_for()"""
    rows = []
    for tag, E, t, o, fr in spec:
        thl = t if fr == 'LAB' else float(DB.th_cm2lab(E, t))
        thc = t if fr == 'CM' else float(DB.th_lab2cm(E, t))
        rows.append(dict(pid='%s_%s_%.3f_%.3f_%s' % (tag, o, E, t, fr), dataset='grid', E=E, theta=t, theta_lab=thl,
                         theta_cm=thc, dobs=o, frame=fr, value=0.0, err_pt=1.0, norm=0.0, flags=''))
    return pd.DataFrame(rows)


def forward_point(E):
    """as in s06_figures.py: empirical forward cross section ('GHP_fwd_CS' rows) at energy E, cubic spline through
    the tabulated values, uncertainty linearly interpolated"""
    d = df[df.dataset == 'GHP_fwd_CS'].sort_values('E')
    return float(CubicSpline(d.E.values, d.value.values)(E)), float(np.interp(E, d.E.values, d.err_pt.values))


g = grid_rows([('pr', E, t, 'CS', 'LAB') for E in PANELS for t in th])
Tg = outputs(g, 'NNLO', tag='fig_pred_ratio')
# flag 1 = loop amplitudes extrapolated in t inside the Gram-singular zones near 0/180 deg (driver README, Sect. 6);
# benign (the forward grid of s06 is entirely flag 1).  flag 2 (E < 15 MeV) must not occur here.
assert (Tg.flag <= 1).all(), 'driver flags on the grid: %s' % np.unique(Tg.flag)
print('grid: %d of %d points with the loops extrapolated near 0/180 deg (flag 1)' % ((Tg.flag == 1).sum(), len(Tg.flag)))
pred_g = Tg.predict_array(PREF_NNLO)
_, trg = truncation_for(g, g, tag='fig_pred_ratio')
sig_g = trg['sigma']

# ---- figure
fig, axs = plt.subplots(2, 4, figsize=(W2, 4.6), sharex=True, layout='constrained')
axs = axs.ravel()
# marker, legend label, x-offset [deg] of the three HIgS measurements (the offsets separate them at the common angles)
HSTY = {'CS': ('s', 'HI$\\gamma$S 22 $d\\sigma/d\\Omega$ (81.3 MeV)', -2.5), 'SPAR': ('^', 'HI$\\gamma$S 22 $d\\sigma_\\parallel$ (83.4 MeV)', 2.5),
        'SPERP': ('v', 'HI$\\gamma$S 22 $d\\sigma_\\perp$ (83.4 MeV)', 0.0)}
leg = {}                        # data legend: label -> handle
outside = []                    # points outside YLIM (reported)
curves = {}                     # report: per panel (fit, PWA, forward) ratios


def jitter(d, step=1.5):
    """x-offsets [deg]: per HIgS observable (HSTY) plus, for repeated measurements of the same (dataset, dobs) at the
    same angle (HIgS beam times), -step, 0, +step ... (symmetric); at most 4 deg in total"""
    off = np.array([HSTY[o][2] if ds in HIGS else 0.0 for ds, o in zip(d.dataset, d.dobs)])
    for (ds, o, t), grp in d.groupby(['dataset', 'dobs', 'theta_lab']):
        k = len(grp)
        off[d.index.get_indexer(grp.index)] += (np.arange(k) - (k - 1) / 2) * step
    return off


for i, E in enumerate(PANELS):
    ax = axs[i]
    m = np.isclose(g.E.values, E)
    ax.fill_between(th, 1 - sig_g[m] / pred_g[m], 1 + sig_g[m] / pred_g[m], color=C_BAND, alpha=0.6, lw=0, zorder=1)
    ax.axhline(1.0, color='k', lw=0.6, zorder=2)
    pw = PWA[(PWA.obs == 'CS') & np.isclose(PWA.E, E)].sort_values('theta_lab')
    assert len(pw), 'no PWA curve at %.2f MeV' % E
    pw_r = pw.value.values / np.interp(pw.theta_lab.values, th, pred_g[m])
    ax.plot(pw.theta_lab, pw_r, ':', color=C_PWA, lw=1.1, zorder=3)
    yl = OWN_YLIM.get(E, YLIM)
    v, e = forward_point(E)
    hf = ax.errorbar([0.0], [v / pred_g[m][0]], [e / pred_g[m][0]], fmt='*', color='k', mfc='white', ms=6.5, lw=0.7,
                     zorder=6)
    leg['forward SR'] = hf
    curves[E] = (pw_r, v / pred_g[m][0], e / pred_g[m][0])
    assert yl[0] < pw_r.min() and max(pw_r.max(), curves[E][1] + curves[E][2]) < yl[1], E
    d = sel[np.isclose(sel.panel, E)].copy()
    d['x'] = d.theta_lab + jitter(d)
    for (ds, o), dd in d.groupby(['dataset', 'dobs'], sort=False):
        if ds in HIGS:
            mk, lab, _ = HSTY[o]; c = MARK[ds][1]; mfc = c
        else:
            (mk, c), lab = MARK[ds], SHORT[ds]
            mfc = 'white' if ds.startswith('Hallin') else c
        vis = (dd.ratio >= yl[0]) & (dd.ratio <= yl[1])
        h = ax.errorbar(dd.x[vis], dd.ratio[vis], dd.eratio[vis], fmt=mk, color=c, mfc=mfc, ms=3.5, lw=0.7, zorder=5)
        leg.setdefault(lab, h)
        for _, r in dd[~vis].iterrows():
            outside.append(r)
            assert r.ratio > yl[1], r.pid                 # arrows are drawn for points above the frame only
            ax.annotate('', xy=(r.x, yl[1] - 0.004), xytext=(r.x, yl[1] - 0.06),
                        arrowprops=dict(arrowstyle='->', color=c, lw=0.8, shrinkA=0, shrinkB=0))
    ax.set_title('%.1f MeV' % E, pad=2)
    ax.set_xlim(-6, 180); ax.set_ylim(*yl)
    ax.yaxis.set_tick_params(labelleft=i % 4 == 0 or E in OWN_YLIM)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    if i % 4 == 0:
        ax.set_ylabel('data / NNLO prediction')
    if i >= 3:
        ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')
        ax.xaxis.set_tick_params(labelbottom=True)

# legends: theory above the panels, data in the empty eighth panel (as in Fig. cs)
axs[7].axis('off')
h_th = [(Patch(facecolor=C_BAND, alpha=0.6, lw=0), Line2D([], [], color='k', lw=0.6)),
        Line2D([], [], color=C_PWA, ls=':', lw=1.1)]
l_th = ['NNLO prediction $\\pm\\,\\sigma_{\\rm th}$', 'PWA fit / prediction']
fig.legend(h_th, l_th, loc='outside upper center', ncol=2, columnspacing=2.0, handletextpad=0.5)
order = [k for k in ['A2 22'] + [v[1] for v in HSTY.values()] + ['A2 01', 'SAL 95', 'SAL 93', 'Moscow 74', 'forward SR']
         if k in leg] + [k for k in leg if k not in ['A2 22', 'A2 01', 'SAL 95', 'SAL 93', 'Moscow 74', 'forward SR']
                         and k not in [v[1] for v in HSTY.values()]]
axs[7].legend([leg[k] for k in order], order, loc='center', ncol=1, handletextpad=0.3, frameon=False)

with plt.rc_context({'savefig.bbox': 'standard'}):      # exact size, no cropping (included without scaling)
    for dd in (FIGS, PAPER_FIGS):
        fig.savefig(os.path.join(dd, 'fig_pred_ratio.pdf'))
plt.close(fig)

# ---- report
print('fig_pred_ratio.pdf written to %s and %s' % (FIGS, PAPER_FIGS))
print('normalisations applied to the data (preferred by the prediction, n_th):')
for ds, d in sel.groupby('dataset'):
    print('  %-24s n=%.3f  N=%2d  ratio %.3f-%.3f  sigma_th/pred %.3f-%.3f' % (
        ds, d.n.iloc[0], len(d), d.ratio.min(), d.ratio.max(), (d.sig_th / d.pred).min(), (d.sig_th / d.pred).max()))
if outside:
    print('points above the y-range (arrows; values to be quoted in the caption):')
    for r in outside:
        print('  %s  E=%.1f theta=%.0f %s  value %.2f(%.2f)  n %.3f  prediction %.3f  ratio %.2f(%.2f)' % (
            r.pid, r.E, r.theta_lab, r.dobs, r.value, r.err_pt, r.n, r.pred, r.ratio, r.eratio))
for E in PANELS:
    p, v, e = curves[E]
    m = np.isclose(g.E.values, E)
    print('  E=%6.2f: PWA/pred %.3f-%.3f (0 deg %.3f); forward SR/pred %.3f(%.3f); sigma_th/pred %.3f-%.3f' % (
        E, p.min(), p.max(), p[0], v, e, (sig_g[m] / pred_g[m]).min(), (sig_g[m] / pred_g[m]).max()))
