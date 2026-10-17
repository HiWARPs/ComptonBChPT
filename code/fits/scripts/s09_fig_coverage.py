"""Figure: coverage check of the truncation uncertainty (paper Fig. "coverage").

Uses the NNLO fit 'F2_d_all_thB' of results/fits_L.json (main database = pre-2015 + A2 2022 minus the 3-sigma outlier,
HIgS 2022 excluded since 2026-10-17; truncation covariance, scheme B).  For each ANGULAR data point of the fit (forward
points of FWD excluded; 245 points): the residual r = y - n t (data minus model with the fitted normalisation,
'resid'), the experimental error sigma_exp = err_pt of the database row (statistical (+) point-to-point systematic)
and the truncation sigma_th = 'sig_th' (sqrt of the diagonal of the truncation covariance).  In 10-MeV bins of the
photon lab energy from 50 to 150 MeV (the last bin includes 150 MeV; a bin with fewer than MIN_N points is merged
with its neighbour):
  upper panel: the pulls r/sqrt(sigma_exp^2 + sigma_th^2) of all points (gray dots) and the bin-wise rms of the
               pull with sigma_exp only (open markers, dashed) and with sigma_exp (+) sigma_th (filled, solid);
  lower panel: for the cross sections only (dobs CS, CSCM; SPAR, SPERP would be included if polarised cross sections
               were in the fit -- the HIgS ones are not; the asymmetries Sigma_3 are excluded) the bin-wise rms of
               r/t, sigma_exp/t and sigma_th/t in per cent
               of the model t = n * theory (residual solid black, experimental dotted, truncation dashed red).
Prints the table of bin values and writes it to results/s09_coverage_bins.txt.
Writes work/fits/figs/fig_coverage.pdf and paper/figs/fig_coverage.pdf (final printed size, single column).
Variant:  'median' as the first argument replaces the rms of the lower panel by the median of the absolute values
(robust against the few points with large relative errors); written as fig_coverage_median.pdf and
s09_coverage_bins_median.txt in work/fits only.  The rms of the pulls (upper panel) is unchanged.
Not part of the make pipeline; run by hand:  .venv/bin/python scripts/s09_fig_coverage.py [median]
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *  # noqa

plt.style.use(STYLE)
os.makedirs(FIGS, exist_ok=True)
PAPER_FIGS = os.path.join(os.path.dirname(os.path.dirname(FITS)), 'paper', 'figs')
W1, W2 = 246 / 72.27, 510 / 72.27          # final printed widths [in] (as in s06_figures.py)
C_NNLO, C_F1, C_F2, C_BAND = '#444444', '#1f4e9c', '#c0392b', '#9fb6dd'
EDGES = np.arange(50.0, 150.01, 10.0)
MIN_N = 5
CS_LIKE = ['CS', 'CSCM', 'SPAR', 'SPERP']
STAT = sys.argv[1] if len(sys.argv) > 1 else 'rms'
assert STAT in ('rms', 'median'), STAT
SUFFIX = '' if STAT == 'rms' else '_' + STAT

df = load_db()
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
FIT = FL['F2_d_all_thB']
s = df.set_index('pid').loc[FIT['pids']].reset_index()
s['resid'] = FIT['resid']; s['model'] = FIT['model']; s['sig_th'] = FIT['sig_th']
ang = s[~s.dataset.isin(FWD)].copy()
# main database since 2026-10-17: pre-2015 + A2 2022 minus the 3-sigma outlier = 245 angular points (261 with HIgS until then)
assert len(ang) == 245, len(ang)
assert not ang.dataset.isin(HIGS).any(), 'HIgS data are not part of the main database'
nlow = int((ang.E < EDGES[0]).sum())
ang = ang[(ang.E >= EDGES[0]) & (ang.E <= EDGES[-1])].copy()
ang['pull_exp'] = ang.resid / ang.err_pt
ang['pull_tot'] = ang.resid / np.sqrt(ang.err_pt ** 2 + ang.sig_th ** 2)
ang['cs'] = ang.dobs.isin(CS_LIKE)

# ---- bins (merge bins with fewer than MIN_N points into their neighbour)
bins = [(EDGES[i], EDGES[i + 1]) for i in range(len(EDGES) - 1)]


def inbin(lo, hi):
    return (ang.E >= lo) & ((ang.E < hi) if hi < EDGES[-1] else (ang.E <= hi))


merged = []
while len(bins) > 1:
    counts = [int(inbin(lo, hi).sum()) for lo, hi in bins]
    small = [i for i, n in enumerate(counts) if n < MIN_N]
    if not small:
        break
    i = small[0]                                   # merge with the neighbour that has fewer points
    if i == 0:
        j = 1
    elif i == len(bins) - 1:
        j = i - 1
    else:
        j = i - 1 if counts[i - 1] <= counts[i + 1] else i + 1
    a, b = min(i, j), max(i, j)
    new = (bins[a][0], bins[b][1])
    merged.append('%.0f-%.0f MeV (%d + %d points)' % (new[0], new[1], counts[a], counts[b]))
    bins = bins[:a] + [new] + bins[b + 1:]


def rms(x):
    x = np.asarray(x, float)
    return np.sqrt(np.mean(x ** 2)) if len(x) else np.nan


def medabs(x):
    x = np.asarray(x, float)
    return np.median(np.abs(x)) if len(x) else np.nan


agg = rms if STAT == 'rms' else medabs            # statistic of the lower panel (the pulls always use the rms)
rows = []
for lo, hi in bins:
    m = inbin(lo, hi); d = ang[m]; c = d[d.cs]
    rows.append(dict(lo=lo, hi=hi, Ec=0.5 * (lo + hi), N=int(m.sum()), Ncs=int(len(c)),
                     rms_pull_exp=rms(d.pull_exp), rms_pull_tot=rms(d.pull_tot),
                     rms_res=100 * agg(c.resid / c.model), rms_exp=100 * agg(c.err_pt / c.model),
                     rms_th=100 * agg(c.sig_th / c.model)))
tab = pd.DataFrame(rows)

# ---- figure
fig, axs = plt.subplots(2, 1, figsize=(W1, 4.2), sharex=True, layout='constrained')
ax = axs[0]
h_dots, = ax.plot(ang.E, ang.pull_tot, 'o', color='0.65', ms=1.6, mew=0, alpha=0.8, zorder=1)
ax.axhline(1.0, color='k', lw=0.6, zorder=2)
h_pe, = ax.plot(tab.Ec, tab.rms_pull_exp, 's--', color=C_F1, mfc='white', ms=4, lw=1.0, zorder=4)
h_pt, = ax.plot(tab.Ec, tab.rms_pull_tot, 'o-', color=C_F1, ms=4, lw=1.2, zorder=5)
ax.set_ylabel('pull')
ax.set_ylim(-3.6, 3.6)
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
ax = axs[1]
h_r, = ax.plot(tab.Ec, tab.rms_res, 'o-', color='k', ms=3, lw=1.2)
h_e, = ax.plot(tab.Ec, tab.rms_exp, 's:', color=C_F1, ms=3, lw=1.2)
h_t, = ax.plot(tab.Ec, tab.rms_th, 'D--', color=C_F2, ms=3, lw=1.1)
ax.set_ylabel('%s [\\%%]' % STAT)
ax.set_ylim(0, max(12.0, 1.1 * tab[['rms_res', 'rms_exp', 'rms_th']].max().max()))
ax.set_xlabel('$E_\\gamma$ [MeV]'); ax.set_xlim(50, 150)
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
# one legend for both panels (column-wise: upper-panel entries left, lower-panel entries right)
S = STAT
fig.legend([h_pe, h_pt, h_dots, h_r, h_e, h_t],
           ['rms pull, $\\sigma_{\\rm exp}$ only', 'rms pull, $\\sigma_{\\rm exp}\\oplus\\sigma_{\\rm th}$',
            'pulls, $\\sigma_{\\rm exp}\\oplus\\sigma_{\\rm th}$',
            '%s residual $|y-n\\,t|/t$' % S, '%s $\\sigma_{\\rm exp}/t$' % S, '%s $\\sigma_{\\rm th}/t$' % S],
           loc='outside upper center', ncol=2, columnspacing=1.5, handletextpad=0.5, handlelength=2.2)

with plt.rc_context({'savefig.bbox': 'standard'}):      # exact size, no cropping (included without scaling)
    for dd in ((FIGS, PAPER_FIGS) if STAT == 'rms' else (FIGS,)):
        fig.savefig(os.path.join(dd, 'fig_coverage%s.pdf' % SUFFIX))
plt.close(fig)

# ---- report
lines = ['coverage check, NNLO fit F2_d_all_thB: %d angular points, %d of them in 50-150 MeV (%d below 50 MeV omitted)'
         % (len(ang) + nlow, len(ang), nlow),
         'bins: %s' % ('%d-MeV bins' % (EDGES[1] - EDGES[0]) + ('; merged: %s' % ', '.join(merged) if merged else '; no bin merged')),
         'cross-section-like points (CS, CSCM, SPAR, SPERP) in the lower panel: %d; asymmetries excluded there: %d'
         % (int(ang.cs.sum()), int((~ang.cs).sum())),
         'lower-panel statistic: %s of the relative quantities (in per cent of the model t = n * theory)' % S,
         '%-9s %4s %4s %12s %12s %11s %11s %11s' % ('E [MeV]', 'N', 'Ncs', 'rms pull', 'rms pull', S + ' |r|/t', S + ' s_exp/t', S + ' s_th/t'),
         '%-9s %4s %4s %12s %12s %11s %11s %11s' % ('', '', '', '(s_exp)', '(s_exp+s_th)', '[%]', '[%]', '[%]')]
for r in tab.itertuples():
    lines.append('%3.0f-%3.0f  %4d %4d %12.2f %12.2f %11.2f %11.2f %11.2f' % (
        r.lo, r.hi, r.N, r.Ncs, r.rms_pull_exp, r.rms_pull_tot, r.rms_res, r.rms_exp, r.rms_th))
d = ang; c = ang[ang.cs]
lines.append('%-9s %4d %4d %12.2f %12.2f %11.2f %11.2f %11.2f' % (
    'all', len(d), len(c), rms(d.pull_exp), rms(d.pull_tot), 100 * agg(c.resid / c.model), 100 * agg(c.err_pt / c.model),
    100 * agg(c.sig_th / c.model)))
txt = '\n'.join(lines) + '\n'
print(txt)
open(os.path.join(RES, 's09_coverage_bins%s.txt' % SUFFIX), 'w').write(txt)
print('fig_coverage%s.pdf written to %s' % (SUFFIX, ' and '.join((FIGS, PAPER_FIGS) if STAT == 'rms' else (FIGS,))))
