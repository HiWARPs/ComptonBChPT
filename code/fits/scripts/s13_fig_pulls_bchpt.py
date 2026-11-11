"""Fig. pulls of the paper: chi^2 contributions of the individual data points in the BChPT fits up to 360 MeV (the
analogue of the PWA figure pwapulls, lower panel).  Upper panel: the main fit, which leaves out the LEGS data; the LEGS
points are drawn as open symbols with their pulls evaluated at the main-fit parameters (normalisations 1), as the
excluded data.  Lower panel: the same fit with the LEGS data included (variant combined_withLEGS).  Pulls are
UNWEIGHTED: residual / sqrt(err_pt^2 + truncation variance (below threshold) or NNLO band (Delta region)); the
forward-SR amplitude block above threshold is not a set of data points and is left out.  Prints, per data set, how
many points exceed the 3 sigma level (chi^2_i > 9), the criterion that removed the single A2 2001 point from the main
database (scripts/s03_fits_L.py).  Output: work/fits/figs/fig_pulls_bchpt.pdf (-> paper/figs/fit_pulls.pdf via
paper/figs/collect_figs.sh) and results/pulls_bchpt.json.
Run from work/fits with the project venv:  ../../.venv/bin/python scripts/s13_fig_pulls_bchpt.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lib')); sys.path.insert(0, os.path.join(HERE, '..', 'scenarios'))
os.chdir(os.path.join(HERE, '..'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from common import STYLE, FIGS, FWD, MARK, SHORT, save_json
import combined_fit as CF

plt.style.use(STYLE)
W1 = 246 / 72.27                        # printed single-column width [in]
NSIG2 = 9.0                             # (3 sigma)^2


def legs_at(res):
    """pulls of the LEGS points at the parameters of a fit that excludes them (normalisations 1, band in quadrature)"""
    sD, sth = CF.rows_delta(res['Emax'], legs=True)
    m = sD.dataset.isin(CF.LEGS).values
    sD, sth = sD[m], sth[m]
    idx = CF.T0.subset(sD.pid.values)
    O = np.array([CF.G[g].predict_array(np.array(res['P']), idx) for g in CF.GE_GRID])
    t = CubicSpline(CF.grid, O, axis=0)(res['gE'][0])
    return pd.DataFrame(dict(dataset=list(sD.dataset), E=sD.E.values, theta_lab=sD.theta_lab.values,
                             pull=(sD.value.values - t) / np.sqrt(sD.err_pt.values ** 2 + sth ** 2)))


fits = {'main': CF.fit('combined'), 'legs': CF.fit('combined_withLEGS', legs=True)}
fig, axs = plt.subplots(2, 1, figsize=(W1, 5.3), sharex=True, layout='constrained')
reg, counts = {}, {}
for ax, key, title in ((axs[0], 'main', 'Main fit (LEGS data excluded)'), (axs[1], 'legs', 'Fit with the LEGS data')):
    d = pd.DataFrame({k: fits[key]['points'][k] for k in ('dataset', 'E', 'theta_lab', 'pull')})
    d['open'] = False
    if key == 'main':
        ex = legs_at(fits['main']); ex['open'] = True
        d = pd.concat([d, ex], ignore_index=True)
    c2 = d.pull ** 2
    low = (d.E <= 150) & ~d.dataset.isin(FWD)
    h, = ax.plot(d.E[low], c2[low], ls='none', marker='o', ms=1.8, color='0.6', mew=0)
    reg.setdefault('below threshold', h)
    fw = d.dataset.isin(FWD)
    h, = ax.plot(d.E[fw], c2[fw], ls='none', marker='*', ms=4.5, color='k', mfc='white', mew=0.6)
    reg.setdefault('forward SRs', h)
    cnt = {'below threshold': (int(low.sum()), int((c2[low] > NSIG2).sum()))}
    for ds, g in d[(d.E > 150)].groupby('dataset', sort=False):
        mk, c = MARK.get(ds, ('o', 'k'))
        lab = SHORT.get(ds, ds) + (r' $\Sigma_3$' if ds.endswith('_S3') else '')
        for o, gg in g.groupby('open'):
            h, = ax.plot(gg.E, gg.pull ** 2, ls='none', marker=mk, ms=3.0, mew=0.6, color=c, mfc='none' if o else c)
            reg.setdefault(lab + (' (excluded)' if o else ''), h)
        cnt[ds] = (len(g), int((g.pull ** 2 > NSIG2).sum()), float((g.pull ** 2).max()))
    counts[key] = cnt
    ax.axhline(NSIG2, ls='--', color='0.3', lw=0.8)
    ax.text(362, NSIG2 + 0.4, r'$3\sigma$', ha='right', va='bottom')
    ax.set_ylabel(r'$\chi^2_i$')
    ax.set_ylim(0, max(21.0, 1.1 * float(c2.max())))
    ax.set_title(title, loc='left', fontsize=8, pad=2.5)
    print(key, {k: v for k, v in cnt.items()})
    for _, r in d[c2 > NSIG2].sort_values('pull', key=abs, ascending=False).iterrows():
        print('   %-20s E=%6.1f th=%5.1f chi2_i=%5.1f%s' % (r.dataset, r.E, r.theta_lab, r.pull ** 2, ' (excluded)' if r.open else ''))
axs[1].set_xlabel(r'$E_\gamma$ [MeV]')
axs[1].set_xlim(45, 365)
labs = [l for l in reg if l not in ('below threshold', 'forward SRs')] + ['below threshold', 'forward SRs']
fig.legend([reg[l] for l in labs], labs, loc='outside upper center', ncol=3, handlelength=1.0, handletextpad=0.4,
           columnspacing=1.0, fontsize=6.5)
with plt.rc_context({'savefig.bbox': 'standard'}):       # keep the full 246 pt width
    fig.savefig(os.path.join(FIGS, 'fig_pulls_bchpt.pdf'))
save_json(dict(note='per set: (N, N with chi2_i > 9, max chi2_i); pulls unweighted, with band', counts=counts,
               P=dict(main=fits['main']['P'], legs=fits['legs']['P'])), 'pulls_bchpt.json')
print('written', os.path.join(FIGS, 'fig_pulls_bchpt.pdf'))
