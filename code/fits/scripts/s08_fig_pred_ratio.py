"""Figure: the parameter-free NNLO prediction versus the new data (paper Fig. "predratio").

Ratio data/prediction at the published normalisation (n = 1) versus the lab angle, for the HIgS data of Li et al.
2022 (panel 1: unpolarised cross section at 81.3 MeV and the polarised cross sections dsigma_par, dsigma_perp at
83.4 MeV) and the A2 2022 unpolarised cross sections at their five energies (panels 2-6).  The prediction is the
NNLO amplitude at P = PREF_NNLO (no free polarisability), evaluated per data point with the driver ('db' run).
Drawn in every panel: the line at 1; the truncation band 1 +- sigma_th/prediction (theory.truncation via
truncation_for on a fine angular grid, same settings as the fit covariance: scheme B, sqrtQ, Lambda_b = 650 MeV,
cbar = rms); the NNLO fit (fits_L 'F2_d_all_thB') divided by the prediction, both computed with predict_array on
the grid (solid red); and the normalisation n_th the prediction prefers for the set (scenarios/
nnlo_prediction_table.json, 'n_th': fit of the normalisation alone with the truncation covariance) as a dashed
line.  The band and the fit/prediction curve of the HIgS panel are those of the unpolarised cross section at
81.3 MeV (reference; the polarised points are divided by their own predictions).  Points of the same set at the
same angle (HIgS beam times) are displaced by +-2 deg for visibility.
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
from common import *  # noqa
import database as DB

plt.style.use(STYLE)
os.makedirs(FIGS, exist_ok=True)
PAPER_FIGS = os.path.join(os.path.dirname(os.path.dirname(FITS)), 'paper', 'figs')
W1, W2 = 246 / 72.27, 510 / 72.27          # final printed widths [in] (as in s06_figures.py)
C_NNLO, C_F2, C_BAND = '#444444', '#c0392b', '#9fb6dd'
C_A2, C_HIGS = '#1f4e9c', '#d62e8a'         # MARK colours of s06: A2 22 blue circles, HIgS 22 magenta squares
YLIM = (0.7, 1.3)                           # common y-range; points outside are marked and reported

df = load_db()
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
FIT = FL['F2_d_all_thB']
P_FIT = np.array(FIT['P'])
NTH = json.load(open(os.path.join(FITS, 'scenarios', 'nnlo_prediction_table.json')))
n_th = {k: NTH[k]['n_th'] for k in ['A2_2022_CS', 'HIGS_2022_CIRC', 'HIGS_2022_LIN']}

# ---- data points: ratio to the prediction at n = 1 (prediction = driver value at P = Pref, i.e. value_file)
T = outputs(df, 'NNLO')                     # cached 'db' run of all database rows
sel = df[df.dataset.isin(['A2_2022_CS', 'HIGS_2022_CIRC', 'HIGS_2022_LIN'])].copy()
idx = T.subset(sel.pid.values)
sel['pred'] = T.predict_array(PREF_NNLO, idx)
assert np.allclose(sel['pred'].values, T.value_file[idx])
sel['ratio'] = sel.value / sel.pred
sel['eratio'] = sel.err_pt / sel.pred
# truncation sigma at the data points (identical to FIT['sig_th'] at the same pids; checked below)
_, tr = truncation_for(df, sel, tag='db')
sel['sig_th'] = tr['sigma']
fit_sig = dict(zip(FIT['pids'], FIT['sig_th']))
assert np.allclose(sel.sig_th.values, [fit_sig[p] for p in sel.pid])

# ---- curves on a fine angular grid: CS at the panel energies
E_A2 = sorted(df[df.dataset == 'A2_2022_CS'].E.unique())          # 92.2, 103.2, 113.6, 124.1, 135.3
E_H = 81.3
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


g = grid_rows([('pr', E, t, 'CS', 'LAB') for E in [E_H] + list(E_A2) for t in th])
Tg = outputs(g, 'NNLO', tag='fig_pred_ratio')
# flag 1 = loop amplitudes extrapolated in t inside the Gram-singular zones near 0/180 deg (driver README, Sect. 6);
# benign (the forward grid of s06 is entirely flag 1).  flag 2 (E < 15 MeV) must not occur here.
assert (Tg.flag <= 1).all(), 'driver flags on the grid: %s' % np.unique(Tg.flag)
print('grid: %d of %d points with the loops extrapolated near 0/180 deg (flag 1)' % ((Tg.flag == 1).sum(), len(Tg.flag)))
pred_g = Tg.predict_array(PREF_NNLO)
fit_g = Tg.predict_array(P_FIT)
_, trg = truncation_for(g, g, tag='fig_pred_ratio')
sig_g = trg['sigma']

# ---- figure
PANELS = [('HIGS', E_H)] + [('A2', E) for E in E_A2]
fig, axs = plt.subplots(2, 3, figsize=(W2, 4.3), sharex=True, sharey=True, layout='constrained')
axs = axs.ravel()
# marker, legend label, x-offset [deg] of the three HIgS measurements (the offsets separate them at the common angles)
HSTY = {'CS': ('s', 'HI$\\gamma$S 22 $d\\sigma/d\\Omega$ (81.3 MeV)', -2.5), 'SPAR': ('^', 'HI$\\gamma$S 22 $d\\sigma_\\parallel$ (83.4 MeV)', 2.5),
        'SPERP': ('v', 'HI$\\gamma$S 22 $d\\sigma_\\perp$ (83.4 MeV)', 0.0)}
hd, ld = [], []                 # data legend (lower row)
outside = []                    # points outside YLIM (reported)


def jitter(d, step=1.5):
    """x-offsets [deg]: per observable (HSTY) plus, for repeated measurements of the same (dataset, dobs) at the
    same angle (HIgS beam times), -step, 0, +step ... (symmetric); at most 4 deg in total"""
    off = np.array([HSTY[o][2] for o in d.dobs])
    for (ds, o, t), grp in d.groupby(['dataset', 'dobs', 'theta_lab']):
        k = len(grp)
        off[d.index.get_indexer(grp.index)] += (np.arange(k) - (k - 1) / 2) * step
    return off


for i, (kind, E) in enumerate(PANELS):
    ax = axs[i]
    m = np.isclose(g.E.values, E)
    ax.fill_between(th, 1 - sig_g[m] / pred_g[m], 1 + sig_g[m] / pred_g[m], color=C_BAND, alpha=0.6, lw=0, zorder=1)
    ax.axhline(1.0, color='k', lw=0.6, zorder=2)
    ax.plot(th, fit_g[m] / pred_g[m], '-', color=C_F2, lw=1.2, zorder=3)
    if kind == 'A2':
        ax.axhline(n_th['A2_2022_CS'], color=C_A2, ls='--', lw=0.9, zorder=2)
        d = sel[(sel.dataset == 'A2_2022_CS') & np.isclose(sel.E, E)]
        h = ax.errorbar(d.theta_lab, d.ratio, d.eratio, fmt='o', color=C_A2, ms=3.5, lw=0.7, zorder=5)
        if 'A2 22 $d\\sigma/d\\Omega$' not in ld:
            hd.append(h); ld.append('A2 22 $d\\sigma/d\\Omega$')
        ax.text(0.04, 0.93, 'A2, %.1f MeV' % E, transform=ax.transAxes, va='top')
    else:
        ax.axhline(n_th['HIGS_2022_CIRC'], color=C_HIGS, ls='--', lw=0.9, zorder=2)
        ax.axhline(n_th['HIGS_2022_LIN'], color=C_HIGS, ls=':', lw=0.9, zorder=2)
        d = sel[sel.dataset.isin(['HIGS_2022_CIRC', 'HIGS_2022_LIN'])].copy()
        d['x'] = d.theta_lab + jitter(d)
        for o in ['CS', 'SPAR', 'SPERP']:
            dd = d[d.dobs == o]
            mk, lab, _ = HSTY[o]
            vis = (dd.ratio >= YLIM[0]) & (dd.ratio <= YLIM[1])
            h = ax.errorbar(dd.x[vis], dd.ratio[vis], dd.eratio[vis], fmt=mk, color=C_HIGS, ms=3.5, lw=0.7, zorder=5)
            hd.append(h); ld.append(lab)
            for _, r in dd[~vis].iterrows():
                outside.append(r)
        ax.text(0.04, 0.93, 'HI$\\gamma$S', transform=ax.transAxes, va='top')
    # points of this panel outside the common range (also from the A2 sets, if any)
    if kind == 'A2':
        for _, r in d[(d.ratio < YLIM[0]) | (d.ratio > YLIM[1])].iterrows():
            outside.append(r)
    ax.set_xlim(0, 180); ax.set_ylim(*YLIM)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(45))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    if i % 3 == 0:
        ax.set_ylabel('data / NNLO prediction')
    if i >= 3:
        ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')

# mark the points above the range (arrow at their angle, values in the panel)
if outside:
    ax = axs[0]
    above = [r for r in outside if r.ratio > YLIM[1]]
    if above:
        xs = sorted(set(float(r.theta_lab) for r in above))
        for x in xs:
            ax.annotate('', xy=(x, YLIM[1] - 0.004), xytext=(x, YLIM[1] - 0.06),
                        arrowprops=dict(arrowstyle='->', color=C_HIGS, lw=0.8, shrinkA=0, shrinkB=0))
        txt = ', '.join('%.2f(%.0f)' % (r.ratio, 100 * r.eratio) for r in sorted(above, key=lambda r: r.ratio))
        assert all(r.dobs == 'SPAR' and r.theta_lab == 90.0 for r in above), [r.pid for r in above]
        ax.text(0.97, 0.93, '$d\\sigma_\\parallel(90^\\circ)$:\n%s' % txt, transform=ax.transAxes, ha='right', va='top',
                linespacing=1.1)

# legends: theory (upper row), data (lower row); both outside the panels
h_th = [Line2D([], [], color=C_F2, lw=1.2), Patch(facecolor=C_BAND, alpha=0.6, lw=0),
        Line2D([], [], color=C_A2, ls='--', lw=0.9), Line2D([], [], color=C_HIGS, ls='--', lw=0.9),
        Line2D([], [], color=C_HIGS, ls=':', lw=0.9)]
l_th = ['NNLO fit / prediction', 'prediction $\\pm\\,\\sigma_{\\rm th}$',
        '$n_{\\rm th}=%.3f$ (A2)' % n_th['A2_2022_CS'], '$n_{\\rm th}=%.3f$ (HI$\\gamma$S unpol.)' % n_th['HIGS_2022_CIRC'],
        '$n_{\\rm th}=%.3f$ (HI$\\gamma$S pol.)' % n_th['HIGS_2022_LIN']]
fig.legend(h_th, l_th, loc='outside upper center', ncol=3, columnspacing=2.0, handletextpad=0.5)
fig.legend(hd, ld, loc='outside lower center', ncol=4, columnspacing=1.5, handletextpad=0.3)

with plt.rc_context({'savefig.bbox': 'standard'}):      # exact size, no cropping (included without scaling)
    for dd in (FIGS, PAPER_FIGS):
        fig.savefig(os.path.join(dd, 'fig_pred_ratio.pdf'))
plt.close(fig)

# ---- report
print('fig_pred_ratio.pdf written to %s and %s' % (FIGS, PAPER_FIGS))
print('y-range %s; n_th: %s' % (YLIM, {k: round(v, 4) for k, v in n_th.items()}))
print('ratio ranges (data/prediction, n = 1):')
for ds, d in sel.groupby('dataset'):
    print('  %-16s N=%2d  ratio %.3f-%.3f  (with errors %.3f-%.3f)  sigma_th/pred %.3f-%.3f' % (
        ds, len(d), d.ratio.min(), d.ratio.max(), (d.ratio - d.eratio).min(), (d.ratio + d.eratio).max(),
        (d.sig_th / d.pred).min(), (d.sig_th / d.pred).max()))
if outside:
    print('points outside the y-range (marked by an arrow and quoted in the panel):')
    for r in outside:
        print('  %s  E=%.1f theta=%.0f %s  value %.2f(%.2f)  prediction %.3f  ratio %.3f(%.3f)' % (
            r.pid, r.E, r.theta_lab, r.dobs, r.value, r.err_pt, r.pred, r.ratio, r.eratio))
else:
    print('all points inside the y-range')
for E in [E_H] + list(E_A2):
    m = np.isclose(g.E.values, E)
    print('  E=%.1f: fit/prediction %.3f-%.3f; sigma_th/prediction %.3f-%.3f' % (
        E, (fit_g[m] / pred_g[m]).min(), (fit_g[m] / pred_g[m]).max(), (sig_g[m] / pred_g[m]).min(), (sig_g[m] / pred_g[m]).max()))
