"""Publication figures for the updated PWA (work/pwa/figs/*.pdf), shared style work/figstyle.mplstyle.
All fits are on the MAIN DATABASE (update_fits.select drops UF.MAIN_DROP = TAPS 108.42/133 and, since 2026-10-17, the
HIGS 2022 sets by default; 'new' = A2 + HIGS is kept as a test with with_higs=True); the pulls figure draws the
excluded TAPS point as an open symbol with its chi2 at the fitted parameters (the HIGS points are not drawn there: they
are not part of the database of the fit).  The HIGS panels of fig_cs show the PWA fit as a prediction.  The paper's
configuration (legacy normalisations floated) is obtained by running this file through
scenarios/make_figures_legacynorm.py."""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import compton as C
import fortran_runner as FR
import pwa_model as PM
import plb_inputs as PI
import reproduce_plb as RP
import update_fits as UF
from bchpt_multipoles import BChPT

plt.style.use((_WORK + '/figstyle.mplstyle'))
COL = dict(main='#1f4e9c', pruned='#c0392b', new='#2e8b57', plb1='#7f7f7f', plb1pp='#e08e0b', bchpt='k', born='#6a3d9a')
HB = 197.327
W1, W2 = 246 / 72.27, 510 / 72.27      # final printed widths [in] (single / double column); included without scaling


def save_exact(fig, path):
    """save at the exact figure size (no cropping)"""
    with plt.rc_context({'savefig.bbox': 'standard'}):
        fig.savefig(path)
    plt.close(fig)


# experiment labels and symbols of the pull figure: one legend for both panels, the same label and symbol for
# a data set in the KLP18 and in the updated database (labels and symbols as in the Delta-region figures)
_TEAL, _GREY, _MAG, _PURPLE, _BROWN = '#138d90', '#3d3d3d', '#c2185b', '#6a3d9a', '#8c564b'
PULL_LABEL = {'TAPS2001': 'TAPS 01', 'OlmosDeLeon': 'TAPS 01', 'Federspiel1991': 'Illinois 91', 'Federspiel': 'Illinois 91',
              'Zieger1992': 'MAMI 92', 'Zieger': 'MAMI 92', 'Baranov1974_th90': 'Moscow 74', 'Baranov1974_th150': 'Moscow 74',
              'Baranov': 'Moscow 74', 'MacGibbon1995_tagged': 'SAL 95', 'MacGibbon1995_untagged': 'SAL 95', 'McGibbon': 'SAL 95',
              'Hallin1993': 'SAL 93', 'Hallin': 'SAL 93', 'Oxley': 'Chicago 58', 'Goldansky': 'Moscow 60',
              'Bernardini': 'Illinois 60', 'Pugh': 'MIT 57', 'Hyman': 'MIT 59', 'A2_2022_CS': r'A2 22 $d\sigma$',
              'A2_2022_S3': r'A2 22 $\Sigma_3$', 'HIGS_2022_circ': r'HI$\gamma$S 22 circ.', 'HIGS_2022_lin': r'HI$\gamma$S 22 lin.'}
PULL_STYLE = {'TAPS 01': ('h', _TEAL, True), 'Illinois 91': ('^', _PURPLE, False), 'MAMI 92': ('*', _TEAL, True),
              'Moscow 74': ('v', _BROWN, False), 'SAL 95': ('o', _GREY, False), 'SAL 93': ('s', _GREY, False),
              'Chicago 58': ('>', _BROWN, False), 'Moscow 60': ('<', _BROWN, False), 'Illinois 60': ('X', _PURPLE, True),
              'MIT 57': ('d', _PURPLE, False), 'MIT 59': ('d', _PURPLE, True), r'A2 22 $d\sigma$': ('s', _TEAL, True),
              r'A2 22 $\Sigma_3$': ('s', _TEAL, False), r'HI$\gamma$S 22 circ.': ('p', _MAG, True),
              r'HI$\gamma$S 22 lin.': ('p', _MAG, False)}


# ------------------------------------------------------------------ fits
def fits():
    from main_fits import MAIN
    d = UF.load()
    F = {}
    m = UF.build_model(UF.select(d), **MAIN); F['main'] = (m, m.fit())
    m = UF.build_model(UF.select(d, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points']), **MAIN)
    F['pruned'] = (m, m.fit())
    m = UF.build_model(UF.select(d, new_only=True, with_higs=True), **MAIN); F['new'] = (m, m.fit())   # A2 + HIGS (test)
    d0 = PI.load_klp2018()
    m = RP.make_model(RP.select(d0, '1')); F['plb1'] = (m, m.fit())
    m = RP.make_model(RP.select(d0, "1''")); F['plb1pp'] = (m, m.fit())
    return F


def nearest_bin(m, E):
    return int(np.argmin(abs(m.bins - E)))


# ------------------------------------------------------------------ (i) pulls
def place_labels(fig, ax, xy, txt, avoid=(), hlines=()):
    """annotate the points xy with txt next to them, choosing for each label (largest y first) the first
    position, among a few around the point, that stays inside the axes and overlaps neither a data marker,
    another label, the artists `avoid`, nor the horizontal lines `hlines` (layout only, no data change)"""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pts = np.vstack([ax.transData.transform(np.c_[l.get_xdata(), l.get_ydata()]) for l in ax.lines
                     if l.get_marker() not in (None, 'None', '') and len(l.get_xdata())])
    taken = [a.get_window_extent(r) for a in avoid]
    hl = [ax.transData.transform((0, y))[1] for y in hlines]
    axbb = ax.get_window_extent(r)
    cands = [(3, 2, 'left'), (3, -8, 'left'), (-3, 2, 'right'), (-3, -8, 'right'), (3, 6, 'left'), (-3, 6, 'right'),
             (3, -12, 'left'), (-3, -12, 'right'), (5, 10, 'left'), (-5, 10, 'right')]
    order = np.argsort([-y for _, y in xy])
    for i in order:
        best = None
        for dx, dy, ha in cands:
            a = ax.annotate(txt[i], xy[i], xytext=(dx, dy), textcoords='offset points', ha=ha)
            bb = a.get_window_extent(r).expanded(1.0, 1.0)
            pad = 1.5
            nbad = (bb.x0 < axbb.x0 or bb.x1 > axbb.x1 or bb.y0 < axbb.y0 or bb.y1 > axbb.y1) * 10
            nbad += sum(bb.overlaps(t) for t in taken) * 5
            nbad += sum(bb.y0 - pad < y < bb.y1 + pad for y in hl) * 3
            nbad += int(np.sum((pts[:, 0] > bb.x0 - pad) & (pts[:, 0] < bb.x1 + pad) &
                               (pts[:, 1] > bb.y0 - pad) & (pts[:, 1] < bb.y1 + pad)))
            if best is None or nbad < best[0]:
                if best is not None:
                    best[1].remove()
                best = (nbad, a, bb)
            else:
                a.remove()
            if nbad == 0:
                break
        taken.append(best[2])


def fig_pulls(F):
    fig, axs = plt.subplots(2, 1, figsize=(W1, 5.3), sharex=True, layout='constrained')
    reg = {}
    for ax, key, title in ((axs[0], 'plb1', 'KLP18 database, Fit 1 (reproduced)'), (axs[1], 'main', 'Updated database, PWA fit')):
        m, r = F[key]
        pulls, _ = m.pulls(r.p)
        d = m.d
        sets = list(dict.fromkeys(d.dataset))
        jitter = lambda i: 0.6 * (i - len(sets) / 2) / len(sets) * 4
        for i, s in enumerate(sets):
            sel = (d.dataset == s).values
            E = d.E.values[sel] + jitter(i)
            lab = PULL_LABEL.get(s, s.replace('_', ' '))
            mk, c, filled = PULL_STYLE.get(lab, ('o', 'k', True))
            h, = ax.plot(E, pulls[sel]**2, ls='none', marker=mk, ms=2.9 * (1.3 if mk in '*X' else 1.0), mew=0.6, color=c,
                         mfc=c if filled else 'none')
            reg.setdefault(lab, h)
        # points excluded from the main database (UF.MAIN_DROP): open symbol, chi2 contribution evaluated at the
        # fitted parameters (same parameter layout: same normalisation groups and bins), labelled "excluded"
        xy_ex, txt_ex, chi2_ex = [], [], 0.0
        if key == 'main' and UF.MAIN_DROP:
            from main_fits import MAIN
            ma = UF.build_model(UF.select(UF.load(), keep_outlier=True), **MAIN)
            assert ma.names == m.names, 'parameter layout of the all-points model differs from the main fit'
            pa, _ = ma.pulls(r.p)
            da = ma.d
            for ds, E0, th0 in UF.MAIN_DROP:
                j = int(np.where((da.dataset == ds) & (abs(da.E - E0) < 0.05) & (abs(da.theta_lab - th0) < 0.5))[0][0])
                lab = PULL_LABEL.get(ds, ds.replace('_', ' '))
                mk, c, _ = PULL_STYLE.get(lab, ('o', 'k', True))
                x = da.E.values[j] + jitter(sets.index(ds) if ds in sets else 0)
                ax.plot([x], [pa[j]**2], ls='none', marker=mk, ms=2.9 * (1.3 if mk in '*X' else 1.0), mew=0.6, color=c, mfc='none')
                xy_ex.append((da.E.values[j], pa[j]**2))
                txt_ex.append(f'{lab} ({da.E.values[j]:.0f}, {da.theta_lab.values[j]:.0f}' + r'$^\circ$),' + '\nexcluded')   # two lines: fits the panel
                chi2_ex = max(chi2_ex, float(pa[j]**2))
                print(f'fig_pulls: excluded point {ds} {E0} MeV {th0} deg: pull {pa[j]:.2f}, chi2 {pa[j]**2:.1f} at the main-fit parameters')
        ax.axhline(9, ls='--', color='0.3', lw=0.8)
        t3 = ax.text(154, 9.3, r'$3\sigma$', ha='right', va='bottom')
        ax.set_ylabel(r'$\chi^2_i$')
        ax.set_ylim(0, max(21.0, 1.12 * float(np.max(pulls**2)), 1.12 * chi2_ex))
        ax.set_xlim(50, 155)
        ax.set_title(title, loc='left', fontsize=8, pad=2.5)
        out = np.where(np.abs(pulls) > 3)[0]
        txt = [f'{PULL_LABEL.get(d.dataset.values[j], d.dataset.values[j])} ({d.E.values[j]:.0f}, '
               f'{d.theta_lab.values[j]:.0f}' + r'$^\circ$)' for j in out]
        place_labels(fig, ax, [(d.E.values[j], pulls[j]**2) for j in out] + xy_ex, txt + txt_ex, avoid=[t3], hlines=[9.0])
    labs = list(PULL_STYLE)
    labs = [l for l in labs if l in reg] + [l for l in reg if l not in labs]
    fig.legend([reg[l] for l in labs], labs, loc='outside upper center', ncol=3, handlelength=1.0, handletextpad=0.4,
               columnspacing=1.0, markerscale=1.1)
    axs[1].set_xlabel(r'$E_\gamma$ [MeV]')
    axs[1].set_xlim(50, 155)
    save_exact(fig, 'figs/fig_pulls.pdf')


# ------------------------------------------------------------------ (ii) cross sections
def bchpt_obs(E, th, obs):
    M, _ = FR.run(E, th, 'LAB', 'bchpt')
    k = C.kin(E, th, 'lab')
    return C.obs_from_M(M.sum(axis=1), k)[obs]


def born_obs(E, th, obs):
    k = C.kin(E, th, 'lab')
    return C.obs_from_M(C.born_amplitudes(E, th, 'lab'), k)[obs]


def band(ax, m, r, E, th, obs, col, label, alpha=0.3):
    v, e = m.curve_err(E, th, nearest_bin(m, E), r.p, r.cov, obs)
    ax.fill_between(th, v - e, v + e, color=col, alpha=alpha, lw=0)
    ax.plot(th, v, color=col, lw=1.0, label=label)


def fig_cs(F):
    d = UF.load()
    th = np.linspace(20, 160, 71)
    panels = [(81.3, 'HIGS_2022_circ', 'cs_lab', r'$81.3$ MeV'), (83.4, 'HIGS_2022_lin', None, r'$83.4$ MeV (lin.)'),
              (92.2, 'A2_2022_CS', 'cs_lab', r'$92.2$ MeV'), (103.2, 'A2_2022_CS', 'cs_lab', r'$103.2$ MeV'),
              (113.6, 'A2_2022_CS', 'cs_lab', r'$113.6$ MeV'), (124.1, 'A2_2022_CS', 'cs_lab', r'$124.1$ MeV'),
              (135.3, 'A2_2022_CS', 'cs_lab', r'$135.3$ MeV')]
    taps = d[d.dataset == 'TAPS2001']
    fig, axs = plt.subplots(2, 4, figsize=(7.0, 3.9), sharex=True, sharey=True)
    axs = axs.ravel()
    for ax, (E, ds, obs, lab) in zip(axs, panels):
        if obs is None:
            for o, ls, mkr in (('cs_par_lab', '-', 'o'), ('cs_perp_lab', '-', 's')):
                band(ax, *F['main'], E, th, o, COL['main'], None)
                band(ax, *F['pruned'], E, th, o, COL['pruned'], None, alpha=0.2)
                ax.plot(th, bchpt_obs(E, th, o), color=COL['bchpt'], ls='--', lw=0.9)
                x = d[(d.dataset == ds) & (d.obs == ('CS_PAR' if o == 'cs_par_lab' else 'CS_PERP'))]
                ax.errorbar(x.theta_lab, x.value, x.err, ls='none', marker=mkr, color='k', mfc='w' if 'par' in o else 'k', ms=3, lw=0.7)
            ax.text(0.55, 0.42, r'HI$\gamma$S: $\parallel$ open,' + '\n' + r'$\perp$ filled', transform=ax.transAxes, fontsize=6)
        else:
            band(ax, *F['main'], E, th, obs, COL['main'], 'PWA, MAIN')
            band(ax, *F['pruned'], E, th, obs, COL['pruned'], 'PWA, pruned', alpha=0.2)
            m1, r1 = F['plb1']
            ax.plot(th, m1.curve(E, th, nearest_bin(m1, E), r1.p), color=COL['plb1'], ls=':', lw=1.0, label='KLP18 Fit 1')
            ax.plot(th, bchpt_obs(E, th, obs), color=COL['bchpt'], ls='--', lw=0.9, label=r'B$\chi$PT')
            ax.plot(th, born_obs(E, th, obs), color=COL['born'], ls='-.', lw=0.8, label=r'Born+$\pi^0$')
            x = d[(d.dataset == ds) & (abs(d.E - E) < 0.05) & (d.obs == 'CS')]
            ax.errorbar(x.theta_lab, x.value, x.err, ls='none', marker='o' if 'A2' in ds else 's', color='k', ms=3, lw=0.7,
                        label='A2 2022' if 'A2' in ds else r'HI$\gamma$S 2022')
            t = taps[abs(taps.E - E) < 4.5]
            if len(t):
                ax.errorbar(t.theta_lab, t.value, t.err, ls='none', marker='D', color='#e08e0b', mfc='w', ms=3, lw=0.7,
                            label=r'TAPS 2001 ($|\Delta E|<5$ MeV)')
                ax.text(0.05, 0.06, f'TAPS: {t.E.iloc[0]:.1f}', transform=ax.transAxes, fontsize=6, color='#e08e0b')
        ax.text(0.95, 0.06 if obs else 0.9, lab, transform=ax.transAxes, ha='right', fontsize=7.5)
        ax.set_ylim(0, 26)
        ax.set_xlim(15, 165)
    axs[7].axis('off')
    h, l = axs[2].get_legend_handles_labels()
    h2, l2 = axs[0].get_legend_handles_labels()
    hh = dict(zip(l, h)); hh.update({k: v for k, v in zip(l2, h2) if k not in hh})
    t = axs[3].get_legend_handles_labels()
    for k, v in zip(*t[::-1]):
        hh.setdefault(v, k)
    axs[7].legend(hh.values(), hh.keys(), loc='center', fontsize=7)
    for ax in axs[4:7]:
        ax.set_xlabel(r'$\theta_{\rm lab}$ [deg]')
    axs[3].set_xlabel(r'$\theta_{\rm lab}$ [deg]')
    axs[3].xaxis.set_tick_params(labelbottom=True)
    for ax in (axs[0], axs[4]):
        ax.set_ylabel(r'$d\sigma/d\Omega_{\rm lab}$ [nb/sr]')
    fig.savefig('figs/fig_cs_A2_HIGS.pdf')
    plt.close(fig)


# ------------------------------------------------------------------ (iii) Sigma3
def fig_s3(F):
    d = UF.load()
    th = np.linspace(20, 160, 71)
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.4), sharey=True)
    for ax, E in zip(axs, (92.2, 108.5, 129.5)):
        band(ax, *F['main'], E, th, 'S3', COL['main'], 'PWA, MAIN')
        band(ax, *F['pruned'], E, th, 'S3', COL['pruned'], 'PWA, pruned', alpha=0.2)
        for key, ls, lab in (('plb1', ':', 'KLP18 Fit 1'), ('plb1pp', (0, (4, 2)), "KLP18 Fit 1''")):
            m, r = F[key]
            ax.plot(th, m.curve(E, th, nearest_bin(m, E), r.p, obs='S3'), color=COL[key], ls=ls, lw=1.0, label=lab)
        ax.plot(th, bchpt_obs(E, th, 'S3'), color=COL['bchpt'], ls='--', lw=0.9, label=r'B$\chi$PT')
        ax.plot(th, born_obs(E, th, 'S3'), color=COL['born'], ls='-.', lw=0.8, label=r'Born+$\pi^0$')
        x = d[(d.dataset == 'A2_2022_S3') & (abs(d.E - E) < 0.05)]
        ax.errorbar(x.theta_lab, x.value, x.err, ls='none', marker='o', color='k', ms=3, lw=0.7, label='A2 2022')
        ax.text(0.95, 0.92, f'{E} MeV', transform=ax.transAxes, ha='right', fontsize=7.5)
        ax.set_xlabel(r'$\theta_{\rm lab}$ [deg]')
        ax.set_xlim(15, 165); ax.set_ylim(-1.05, 0.25)
    axs[0].set_ylabel(r'$\Sigma_3$')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=7, fontsize=6.5, bbox_to_anchor=(0.5, -0.08))
    fig.savefig('figs/fig_S3_A2.pdf')
    plt.close(fig)


# ------------------------------------------------------------------ (iv) polarisability plane
def ellipse(ax, x, y, sx, sy, rho, col, label, ls='-'):
    cov = np.array([[sx**2, rho * sx * sy], [rho * sx * sy, sy**2]])
    w, v = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    e = Ellipse((x, y), 2 * np.sqrt(w[1]), 2 * np.sqrt(w[0]), angle=ang, fill=False, color=col, lw=1.0, ls=ls)
    ax.add_patch(e)
    ax.plot(x, y, marker='o', color=col, ms=3, label=label, ls='none')


def fig_plane(F, R):
    fig, ax = plt.subplots(figsize=(3.4, 3.9))
    for key, lab, col, ls in (('plb1', 'KLP18 Fit 1 (repr.)', COL['plb1'], '-'), ('plb1pp', "KLP18 Fit 1'' (repr.)", COL['plb1pp'], '-'),
                              ('main', 'This work: MAIN (all data)', COL['main'], '-'), ('pruned', 'This work: pruned', COL['pruned'], '-'),
                              ('new', 'This work: A2+HIGS only', COL['new'], '--')):
        m, r = F[key]
        ellipse(ax, r.der['alpha_m_beta'], r.der['gamma_pi'], r.der_err['alpha_m_beta'], r.der_err['gamma_pi'],
                r.corr('alpha_m_beta', 'gamma_pi'), col, lab, ls)
    # published KLP18 values (Table III)
    for (x, y, sx, sy, lab, mk) in ((10.5, 7.2, 0.4, 1.0, 'KLP18 Fit 1 (publ.)', 's'), (6.4, 3.5, 0.6, 1.2, "KLP18 Fit 1'' (publ.)", 'D')):
        ax.errorbar(x, y, sy, sx, color='0.4', marker=mk, ms=3, mfc='w', lw=0.7, ls='none', label=lab)
    ax.errorbar(7.3, 7.2, 1.7, 1.0, color='k', marker='*', ms=6, lw=0.8, ls='none', label=r'B$\chi$PT (LMP15)')
    ax.errorbar(10.3, 9.9, 2.0, 1.2, color='#17becf', marker='^', ms=4, lw=0.8, ls='none', label='DR bootstrap (MRPP22)')
    ax.plot(10.5, 7.8, marker='v', color='#8c564b', ms=4, ls='none', label='DR (OdL01/BGLMN98)')
    ax.axvspan(9.2 - 0.5, 9.2 + 0.5, color='0.85', alpha=0.6, lw=0, label='PDG 2026')
    ax.set_xlabel(r'$\alpha_{E1}-\beta_{M1}$ [$10^{-4}$ fm$^3$]')
    ax.set_ylabel(r'$\gamma_\pi$ [$10^{-4}$ fm$^4$]')
    ax.set_xlim(4.5, 13); ax.set_ylim(1, 12.5)
    ax.legend(fontsize=5.6, loc='upper center', ncol=2, handletextpad=0.3, columnspacing=0.8, bbox_to_anchor=(0.5, -0.17))
    fig.savefig('figs/fig_polplane.pdf')
    plt.close(fig)


# ------------------------------------------------------------------ (v) dynamical polarisabilities
DYN = [('alphaE1', r'$\alpha_{E1}(\omega)$', [('EE', 1, '+', 2), ('EE', 1, '-', 1)], 2, 3),
       ('betaM1', r'$\beta_{M1}(\omega)$', [('MM', 1, '+', 2), ('MM', 1, '-', 1)], 2, 3),
       ('gE1E1', r'$\gamma_{E1E1}(\omega)$', [('EE', 1, '+', 1), ('EE', 1, '-', -1)], 3, 4),
       ('gM1M1', r'$\gamma_{M1M1}(\omega)$', [('MM', 1, '+', 1), ('MM', 1, '-', -1)], 3, 4),
       ('gE1M2', r'$\gamma_{E1M2}(\omega)$', [('EM', 1, '+', 6)], 3, 4),
       ('gM1E2', r'$\gamma_{M1E2}(\omega)$', [('ME', 1, '+', 6)], 3, 4)]


def dyn_from_mult(f, E, comb, n):
    w = (E / HB) / np.sqrt(1 + 2 * E / (C.MN_GEV * 1e3))
    return sum(c * f[(a, b, s)] for a, b, s, c in comb) / w**n * 1e4


def dyn_pols(m, r):
    out = {k[0]: [] for k in DYN}
    for ib, E in enumerate(m.bins):
        A, a = m._affine_for(np.array([float(E)]), ib)
        A, a = A[0], a[0]
        npm = A.shape[1]
        for key, lab, comb, n, _ in DYN:
            w = (E / HB) / np.sqrt(1 + 2 * E / (C.MN_GEV * 1e3))
            c = np.zeros(len(PM.ALLKEYS))
            for a_, b_, s_, cc in comb:
                c[PM.ALLKEYS.index((a_, b_, s_))] += cc
            g = c @ A / w**n * 1e4
            val = (c @ (A @ r.p[:npm] + a)) / w**n * 1e4
            err = np.sqrt(g @ r.cov[:npm, :npm] @ g)
            out[key].append((E, val, err))
    return {k: np.array(v) for k, v in out.items()}


def fig_dyn(F, save=True):
    B = BChPT()
    Eg = np.arange(30, 151, 2.0)
    stB = dict(alphaE1=11.25, betaM1=3.87, gE1E1=-3.34, gM1M1=2.96, gE1M2=0.19, gM1E2=1.14)
    ylim = dict(alphaE1=(8, 20), betaM1=(0, 12), gE1E1=(-16, 2), gM1M1=(-1, 7), gE1M2=(-1.5, 2.5), gM1E2=(-1, 6))
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.0), sharex=True)
    res = {}
    for key, col, lab, off in (('main', COL['main'], 'This work: MAIN', -1.0), ('pruned', COL['pruned'], 'This work: pruned', 1.0),
                               ('plb1', COL['plb1'], 'KLP18 Fit 1 (repr.)', -2.0), ('plb1pp', COL['plb1pp'], "KLP18 Fit 1'' (repr.)", 2.0)):
        m, r = F[key]
        dp = dyn_pols(m, r)
        res[key] = {k: v.tolist() for k, v in dp.items()}
        st = r.der
        stat = dict(alphaE1=(st['alphaE1'], r.der_err['alphaE1']), betaM1=(st['betaM1'], r.der_err['betaM1']),
                    gE1E1=(st['gE1E1'], r.der_err['gE1E1']), gM1M1=(st['gM1M1'], r.der_err['gM1M1']),
                    gE1M2=(st['gE1M2'], r.der_err['gE1M2']), gM1E2=(st['gM1E2'], r.der_err['gM1E2']))
        for ax, (k, labk, comb, n, _) in zip(axs.ravel(), DYN):
            v = dp[k]
            ax.errorbar(v[:, 0] + off, v[:, 1], v[:, 2], ls='none', marker='o', ms=2.5, color=col, lw=0.7,
                        label=lab, mfc='w' if key.startswith('plb') else col)
            ax.errorbar(0 + off, stat[k][0], stat[k][1], ls='none', marker='s', ms=2.5, color=col, lw=0.7, mfc='w' if key.startswith('plb') else col)
    for ax, (k, labk, comb, n, _) in zip(axs.ravel(), DYN):
        f = {kk: B(kk, Eg) for kk in PM.L1KEYS}
        ax.plot(Eg, dyn_from_mult(f, Eg, comb, n), color='k', ls='--', lw=0.9, label=r'B$\chi$PT (NNLO)')
        ax.plot(0, stB[k], marker='*', color='k', ms=6, ls='none')
        ax.text(0.05, 0.9, labk, transform=ax.transAxes, fontsize=8)
        ax.set_xlim(-5, 155)
        ax.set_ylim(*ylim[k])
    for ax in axs[1]:
        ax.set_xlabel(r'$E_\gamma$ [MeV]')
    axs[0, 0].set_ylabel(r'[$10^{-4}$ fm$^3$]'); axs[1, 0].set_ylabel(r'[$10^{-4}$ fm$^4$]')
    axs[0, 2].legend(fontsize=5.5, loc='lower left')
    axs[0, 2].set_ylim(-16, 2)
    if save:
        fig.savefig('figs/fig_dynpols.pdf')
    plt.close(fig)
    return res


def forward_model(m, r, ib):
    """model non-Born forward amplitudes fbar, gbar (fm) at bin ib, with 1-sigma errors and covariance."""
    E = np.array([float(m.bins[ib])])
    npm = m.ng + m.nr * m.nb
    p0 = r.p[:npm]
    f0, g0 = m.forward_amps(E, ib, p0)
    Jf = np.zeros(npm); Jg = np.zeros(npm)
    for j in list(range(m.ng)) + [m.ng + ib * m.nr + k for k in range(m.nr)]:
        q = p0.copy(); q[j] += 1.0
        fj, gj = m.forward_amps(E, ib, q)
        Jf[j] = fj[0] - f0[0]; Jg[j] = gj[0] - g0[0]
    C_ = r.cov[:npm, :npm]
    return f0[0], g0[0], Jf, Jg, C_


def fig_forward(F):
    """forward observables dsigma/dOmega_lab(0) = |f|^2+|g|^2 and Sigma_2z(0) = -2 f g/(f^2+g^2) vs energy:
    GHP (table points and smooth representation), GHP shifted by the fitted nuisances, PWA fits at the bin centres."""
    from sumrules import ForwardAmps, MP_GHP, ALPHA, KAPPA, HBARC_MEV
    FA = ForwardAmps()
    Mfm = MP_GHP / HBARC_MEV
    def born(E):
        return -ALPHA / Mfm, -ALPHA * KAPPA**2 * (E / HBARC_MEV) / (2 * Mfm**2)
    def obs(f, g):
        return (f**2 + g**2) * 1e7, -2 * f * g / (f**2 + g**2)
    fo = pd.read_csv((_WORK + '/data/db_forward_obs.csv'))
    m, r = F['main']
    nf, ng = r.p[m.names.index('nf')], r.p[m.names.index('ng')]
    E = np.linspace(50, 150, 201)
    fB, gB = born(E)
    ef, eg = FA.err(E)
    cs0, s0 = obs(FA.fbar(E) + fB, FA.gbar(E) + gB)
    cs1, s1 = obs(FA.fbar(E) + nf * ef + fB, FA.gbar(E) + ng * eg + gB)
    fig, axs = plt.subplots(2, 2, figsize=(W2, 3.9), sharex=True, layout='constrained')
    (a, b), (c, d_) = axs
    a.plot(E, cs0, color='k', lw=0.9, label='GHP')
    a.plot(E, cs1, color=COL['main'], ls='--', lw=0.9, label=rf'GHP shifted ($n_f={nf:.2f}$, $n_g={ng:.2f}$)')
    x = fo[fo.obs == 'CS']; a.errorbar(x.E_lab_MeV, x.value, x.stat_err, ls='none', marker='o', ms=3, color='k', mfc='w', lw=0.7, label='forward points')
    b.plot(E, cs1 / cs0 - 1, color=COL['main'], ls='--', lw=0.9)
    b.axhline(0, color='k', lw=0.6)
    xc = x.value.values / np.interp(x.E_lab_MeV.values, E, cs0) - 1
    b.errorbar(x.E_lab_MeV, xc, x.stat_err / x.value, ls='none', marker='o', ms=3, color='k', mfc='w', lw=0.7)
    c.plot(E, s0, color='k', lw=0.9); c.plot(E, s1, color=COL['main'], ls='--', lw=0.9)
    y = fo[fo.obs == 'S2Z']; c.errorbar(y.E_lab_MeV, y.value, y.stat_err, ls='none', marker='o', ms=3, color='k', mfc='w', lw=0.7)
    d_.plot(E, s1 - s0, color=COL['main'], ls='--', lw=0.9); d_.axhline(0, color='k', lw=0.6)
    d_.errorbar(y.E_lab_MeV, y.value - np.interp(y.E_lab_MeV.values, E, s0), y.stat_err, ls='none', marker='o', ms=3, color='k', mfc='w', lw=0.7)
    for key, col, off, lab in (('main', COL['main'], -0.8, 'PWA fit'), ('pruned', COL['pruned'], 0.8, 'PWA, pruned fit')):
        mm, rr = F[key]
        for ib, Eb in enumerate(mm.bins):
            fb_, gb_, Jf, Jg, C_ = forward_model(mm, rr, ib)
            fBb, gBb = born(Eb)
            f, g = fb_ + fBb, gb_ + gBb
            cs, s2 = obs(f, g)
            dcs = np.array([2 * f * 1e7, 2 * g * 1e7]); den = f**2 + g**2
            ds = np.array([-2 * g / den + 2 * f * g * 2 * f / den**2, -2 * f / den + 2 * f * g * 2 * g / den**2])
            J = np.vstack([Jf, Jg])
            Cfg = J @ C_ @ J.T
            ecs = np.sqrt(dcs @ Cfg @ dcs); es = np.sqrt(ds @ Cfg @ ds)
            ref_cs = np.interp(Eb, E, cs0); ref_s = np.interp(Eb, E, s0)
            kw = dict(ls='none', marker='s', ms=2.8, color=col, lw=0.7, label=lab if ib == 0 else None)
            a.errorbar(Eb + off, cs, ecs, **kw); b.errorbar(Eb + off, cs / ref_cs - 1, ecs / ref_cs, **kw)
            c.errorbar(Eb + off, s2, es, **kw); d_.errorbar(Eb + off, s2 - ref_s, es, **kw)
    a.set_ylabel(r'$d\sigma/d\Omega_{\rm lab}(0^\circ)$ [nb/sr]')
    b.set_ylabel(r'relative to GHP')
    c.set_ylabel(r'$\Sigma_{2z}(0^\circ)$')
    d_.set_ylabel(r'$\Sigma_{2z}-\Sigma_{2z}^{\rm GHP}$')
    b.set_ylim(-0.06, 0.06); d_.set_ylim(-0.05, 0.05)
    for ax in (c, d_):
        ax.set_xlabel(r'$E_\gamma$ [MeV]')
    h, l = a.get_legend_handles_labels()
    fig.legend(h, l, loc='outside upper center', ncol=5, columnspacing=2.0)
    save_exact(fig, 'figs/fig_forward.pdf')


if __name__ == '__main__':
    F = fits()
    fig_forward(F)
    R = json.load(open('results/update_fits.json'))
    fig_pulls(F)
    fig_cs(F)
    fig_s3(F)
    fig_plane(F, R)
    dyn = fig_dyn(F)
    json.dump(dyn, open('results/dynamical_pols.json', 'w'), indent=1, default=float)
    for k, (m, r) in F.items():
        print(k, r.summary(k))


def fig_repro(F):
    """KLP18 Fig. 4 analogue: published curves (digitised from the PDF vectors) vs this reproduction."""
    import digitize_plb as Dg
    cv = Dg.curves()
    d0 = PI.load_klp2018()
    th = np.linspace(0, 180, 91)
    fig, axs = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
    for ax, E in zip(axs, (59, 109)):
        for key, nm, col in (('plb1', 'fit1', COL['plb1']), ('plb1pp', 'fit1pp', COL['plb1pp'])):
            t, v = cv[(E, nm)]
            ax.plot(t, v, color=col, lw=2.2, alpha=0.45, label=f'KLP18 {"Fit 1" if nm == "fit1" else "Fit 1$^{\\prime\\prime}$"} (published)')
            m, r = F[key]
            ax.plot(th, m.curve(float(E), th, PI.BINS.index(E), r.p), color=col, lw=0.9, ls='--',
                    label=f'{"Fit 1" if nm == "fit1" else "Fit 1$^{\\prime\\prime}$"} (this work)')
        t, v = cv[(E, 'BChPT')]
        ax.plot(t, v, color='k', ls=':', lw=1.0, label=r'B$\chi$PT (KLP18 Fig.)')
        thf = np.linspace(4, 176, 87)   # NNLO loop amplitudes numerically unstable within ~3 deg of 0/180 deg
        ax.plot(thf, bchpt_obs(float(E), thf, 'cs_lab'), color='k', lw=0.6, label=r'B$\chi$PT (Fortran, this work)')
        x = d0[d0.bin == E]
        ax.errorbar(x.theta_lab, x.value, x.err, ls='none', marker='o', ms=2.5, color='0.2', lw=0.6, label='KLP18 data')
        ax.text(0.04, 0.06, f'{E} MeV', transform=ax.transAxes)
        ax.set_ylabel(r'$d\sigma/d\Omega_{\rm lab}$ [nb/sr]')
        ax.set_ylim(5, 28)
    axs[1].set_xlabel(r'$\theta_{\rm lab}$ [deg]')
    axs[0].legend(fontsize=5.3, loc='upper center', ncol=2)
    fig.savefig('figs/fig_klp18_reproduction.pdf')
    plt.close(fig)


if __name__ == '__main__':
    fig_repro(F)
