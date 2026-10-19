"""
make_figs.py -- Delta-region figures and chi^2 table (see work/notes/delta_plots.md).

usage:  [AMPL_EXE=...] [DP_ALPHA=.. DP_BETA=.. DP_COV_AB=.. DP_SPINFIT=.. DP_V2OPTS=..] python make_figs.py
        (or ./make_all.sh; see config.py for the environment variables)
Writes  paper/figs/delta/*.pdf, work/delta_plots/png/*.png (renders for inspection),
        work/delta_plots/results/{chi2_table.md,chi2_table.csv,chi2_table.tex,data_with_theory.csv,
        validation.json,curves_*.csv}
"""
import os, sys, json, unicodedata
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import config as C
import data_dp as D
import theory_dp as TD

plt.style.use(C.STYLE)
for d in (C.FIGDIR, C.PNG, C.RESULTS):
    os.makedirs(d, exist_ok=True)

# colours (2026-10-17 evening, user: one presentation for all figures): the same as the below-threshold figures of work/fits
# (s06_figures: prediction line '#444444' with the band '#9fb6dd' at alpha 0.6, NNLO fit dashed '#1f4e9c'); the tree-level
# curve therefore moves from blue to green; the coupling band (R_EM, g_M) is light red, distinct from the blue NNLO band.
C_MAIN, C_BAND, C_GM = '#444444', '#c5d3eb', '#f9c9c3'     # '#c5d3eb' = s06's '#9fb6dd' at alpha 0.6 on white, drawn opaque here
BAND_ALPHA = 1.0                                               # so that the overlap with the light red coupling band stays clean
C_GM_EDGE = '#e0877d'
C_TREE, C_NOVC = '#2e8b57', '#4d4d4d'
C_FIT = '#1f4e9c'                  # the NNLO-fit curve: dashed dark blue, as in the below-threshold figures
# legend labels (short; the captions of the paper give the full description).  Since 2026-10-16 the main curve is
# the parameter-free NNLO prediction (config.MAIN = 'pred'); in the legacy mode 'fit' the labels are swapped.
LAB_MAIN = r'NNLO B$\chi$PT' if C.MAIN == 'pred' else r'NNLO fit'
LAB_FIT = r'NNLO fit' if C.MAIN == 'pred' else r'NNLO prediction'
LAB_GM = r'$\gamma N\Delta$ couplings'      # light band = (g_M, g_E) grid of config.BAND_GRID; the ranges are stated in the captions (user, 2026-10-17)
LAB_TREE = r'no loops'                  # tree level: Born, pi0 pole, Delta pole (user wording, 2026-10-17)
LAB_NOVC = r'no vertex corrections'
LAB_RUN = r'running $\Delta$ parameters'
# extra dotted curve (config.EXTRA, 2026-10-17): 'run3' = the prediction with the running Delta parameters (irun=3), dotted
# in the colour of the main curve, legend entry right after 'NNLO fit'; 'novc' = without the vertex corrections (grey
# dotted, legend entry after the tree-level curve, the figures until 2026-10-17); 'none' = no extra curve.
EXTRA_KEY = {'novc': 'noVC', 'none': None}.get(C.EXTRA, C.EXTRA)      # key of the curve in theory_dp.Theory.curves()
C_RUN = C_MAIN
EXTRA_COLOR = C_NOVC if C.EXTRA == 'novc' else C_RUN
EXTRA_LABEL = LAB_NOVC if C.EXTRA == 'novc' else LAB_RUN
EXTRA_AFTER_FIT = C.EXTRA != 'novc'
EXTRA_LW = 1.0
FIT_DASHES = (3.0, 1.6)
FIT_LW = 0.9


def has_extra(cv):
    return EXTRA_KEY is not None and EXTRA_KEY in cv
MS = 3.2
W1, W2 = 246 / 72.27, 510 / 72.27      # final printed widths [in] (single / double column); included without scaling
E_2PI = 2 * 134.9768 * (1 + 134.9768 / C.MP)   # pi0 pi0 p threshold, lab photon energy [MeV] (308.8)
SHOWN = {}       # figure -> {dataset: number of points shown}
SHOWN_PIDS = set()


def mark_shown(fig, sel):
    SHOWN_PIDS.update(sel.pid)
    d = SHOWN.setdefault(fig, {})
    for k, v in sel.groupby('dataset').size().items():
        d[k] = d.get(k, 0) + int(v)

# ============================================================== kinematic grids
E_GRID = np.r_[np.arange(10.0, 140.0, 5.0), np.arange(140.0, 156.0, 0.25), np.arange(156.0, 800.01, 2.0)]   # Fig. 1: to 800 MeV (450 until 2026-10-17)
TH_GRID = np.r_[0.0, 1.0, np.arange(2.0, 178.01, 2.0), 179.0, 180.0]
FIG1 = [('LAB', 60.0), ('LAB', 90.0), ('LAB', 133.5), ('CM', 75.0), ('CM', 90.0), ('CM', 125.0)]   # 125 cm (LARA, LEGS, SAL) replaces 138 cm since 2026-10-17
FIG2A = [160.0, 182.0, 213.1, 224.3, 234.7, 244.8, 254.1, 265.1]             # rows 1-2 of the angular distributions (part (a) until 2026-10-18)
FIG2B = [275.6, 286.5, 298.2, 310.1, 323.8, 333.8, 350.5, 370.5, 389.0, 444.4, 454.9]   # rows 3-5 (part (b) until 2026-10-18); 444.4/454.9 = LARA angular distributions (2026-10-18)
FIG2 = FIG2A + FIG2B
LEGS_E = [213.1, 224.3, 234.7, 244.8, 254.1, 265.1, 275.6, 286.5, 298.2, 310.1, 323.8, 333.8]
E3B = np.r_[np.arange(10.0, 140.0, 5.0), np.arange(140.0, 156.0, 0.25), np.arange(156.0, 800.01, 2.0)]   # Fig. 3b: to 800 MeV (350 until 2026-10-17)
FIG4 = [('S2X', 288.0, 'A2 15 (273-303 MeV)'), ('S2X', 330.595, 'Martel thesis (315-346 MeV)'),
        ('S2Z', 275.0, 'A2 20 (265-285 MeV)'), ('S2Z', 295.0, 'A2 20 (285-305 MeV)')]
E5 = np.r_[np.arange(50.0, 140.0, 5.0), np.arange(140.0, 156.0, 0.25), np.arange(156.0, 400.01, 2.0)]
FIG5 = [('S2X', 90.0), ('S2X', 130.0), ('S2Z', 90.0), ('S2Z', 130.0)]
E6 = np.arange(150.0, 400.01, 2.0)
FIG6 = [(0, 80, 65.0), (80, 100, 90.0), (100, 125, 112.5), (125, 181, 140.0)]


def gid(*a):
    return 'g_' + '_'.join(str(x) for x in a)


def build_points(df):
    pts = []
    for fr, th in FIG1:
        for E in E_GRID:
            pts.append((gid('f1', fr, th, '%.2f' % E), E, th, 'CS' if fr == 'LAB' else 'CSCM', fr))
    for E in FIG2:
        for th in TH_GRID:
            pts.append((gid('f2', E, th), E, th, 'CSCM', 'CM'))
    for E in LEGS_E:
        for th in TH_GRID:
            pts.append((gid('f3', E, th), E, th, 'S3', 'CM'))
    for E in E3B:
        pts.append((gid('f3b', '%.2f' % E), E, 90.0, 'S3', 'CM'))
    for o, E, _ in FIG4:
        for th in TH_GRID:
            pts.append((gid('f4', o, E, th), E, th, o, 'LAB'))
        # energy grid across the bin (bin-averaged curve)
    for o, E, _ in FIG4:
        lo, hi = BIN[E]
        for Eb in np.linspace(lo, hi, 7):
            for th in TH_GRID:
                pts.append((gid('f4b', o, E, '%.3f' % Eb, th), Eb, th, o, 'LAB'))
    for o, th in FIG5:
        for E in E5:
            pts.append((gid('f5', o, th, '%.2f' % E), E, th, o, 'LAB'))
    for a, b, th in FIG6:
        for E in E6:
            pts.append((gid('f6', th, '%.1f' % E), E, th, 'CSCM', 'CM'))
    pts += D.driver_points(df[df.obs.isin(['CS', 'CSCM', 'S3', 'S2X', 'S2Z'])])
    return pts


BIN = {288.0: (272.73, 303.32), 330.595: (315.25, 345.94), 275.0: (265.0, 285.0), 295.0: (285.0, 305.0)}


# ============================================================== forward points above the pion threshold
# work/data/db_forward_obs_above.csv (GHP forward amplitudes above threshold, when available): columns E_lab_MeV, obs (CS: lab
# dsigma/dOmega(0) in nb/sr; S2Z: Sigma_2z(0)), value, stat_err.  Drawn as open black stars at theta = 0 (as the forward points
# of the below-threshold figures) in Fig. 2 (c.m. cross section: converted with the 0-deg Jacobian) and Fig. 4 (Sigma_2z).
FWD_ABOVE = os.path.join(C.DATA, 'db_forward_obs_above.csv')
FWD_DF = pd.read_csv(FWD_ABOVE) if os.path.exists(FWD_ABOVE) else None
FWD_LABEL = 'forward point'


def forward_above(obs, E0, window):
    """forward points of observable obs within +-window MeV of E0: list of (E, value, err)"""
    if FWD_DF is None:
        return []
    q = FWD_DF[(FWD_DF.obs == obs) & (np.abs(FWD_DF.E_lab_MeV - E0) <= window + 1e-9)]
    return [(float(r.E_lab_MeV), float(r.value), float(r.stat_err)) for r in q.itertuples()]


def draw_forward(ax, pts, reg, cm=False, zorder=6):
    for E, v, e in pts:
        if cm:                                   # lab -> c.m. at theta = 0: dsigma_cm = dsigma_lab / J(0)
            J = D.jac_cm_to_lab(E, 0.0); v, e = v / J, e / J
        ax.errorbar([0.0], [v], [e], fmt='*', color='k', mfc='white', ms=6.5, mew=0.7, lw=0.7, zorder=zorder)
    if pts and FWD_LABEL not in reg:
        reg[FWD_LABEL] = Line2D([], [], ls='none', marker='*', color='k', mfc='white', ms=6.5, mew=0.7)


# ============================================================== drawing helpers
E_BAND_MAX = 450.0     # the bands of the vs-E figures end here (user decision 2026-10-17: upper limit of applicability)


def draw_theory(ax, x, cv, idx, tree=True, fit=True, gm=True, extra=False, band_xmax=None):
    """band_xmax: draw the two bands only for x <= band_xmax (vs-E figures: E_BAND_MAX); the curves continue"""
    x = np.asarray(x, float)
    m = np.ones(len(x), bool) if band_xmax is None else (x <= band_xmax + 1e-9)
    if gm:
        ax.fill_between(x[m], cv['cpl_lo'][idx][m], cv['cpl_hi'][idx][m], color=C_GM, lw=0, zorder=1)
        # edges above the opaque blue band: where the coupling band lies inside it (e.g. above 440 MeV) it stays visible
        ax.plot(x[m], cv['cpl_lo'][idx][m], '-', color=C_GM_EDGE, lw=0.35, zorder=2.5)
        ax.plot(x[m], cv['cpl_hi'][idx][m], '-', color=C_GM_EDGE, lw=0.35, zorder=2.5)
    ax.fill_between(x[m], cv['lo'][idx][m], cv['hi'][idx][m], color=C_BAND, lw=0, alpha=BAND_ALPHA, zorder=2)
    if tree:
        ax.plot(x, cv['tree'][idx], '--', color=C_TREE, lw=0.9, zorder=3, dashes=(4, 2))
    if extra and has_extra(cv):
        ax.plot(x, cv[EXTRA_KEY][idx], ':', color=EXTRA_COLOR, lw=EXTRA_LW, zorder=3.5)
    if fit:
        ax.plot(x, cv['fit'][idx], '--', color=C_FIT, lw=FIT_LW, zorder=3.8, dashes=FIT_DASHES)
    ax.plot(x, cv['main'][idx], '-', color=C_MAIN, lw=1.15, zorder=4)


def theory_handles(fit=True, tree=True, gm=True, extra=False):
    """legend order (user decisions 2026-10-16/17): main (band) / NNLO fit / running Delta parameters / tree / gM band;
    with the legacy extra curve 'novc': main / NNLO fit / tree / no VC / gM band"""
    h, l = [], []
    extra = extra and EXTRA_KEY is not None        # DP_EXTRA=none: never an extra entry
    hx = Line2D([], [], color=EXTRA_COLOR, ls=':', lw=EXTRA_LW)
    h.append((Patch(facecolor=C_BAND, alpha=BAND_ALPHA, edgecolor='none'), Line2D([], [], color=C_MAIN, lw=1.15))); l.append(LAB_MAIN)
    if fit:
        h.append(Line2D([], [], color=C_FIT, ls='--', dashes=FIT_DASHES, lw=FIT_LW)); l.append(LAB_FIT)
    if extra and EXTRA_AFTER_FIT:
        h.append(hx); l.append(EXTRA_LABEL)
    if tree:
        h.append(Line2D([], [], color=C_TREE, ls='--', dashes=(4, 2), lw=0.9)); l.append(LAB_TREE)
    if extra and not EXTRA_AFTER_FIT:
        h.append(hx); l.append(EXTRA_LABEL)
    if gm:
        h.append(Patch(facecolor=C_GM, edgecolor=C_GM_EDGE, lw=0.4)); l.append(LAB_GM)
    return h, l


def draw_points(ax, label, x, y, e, reg, xerr=None, zorder=5, ms=MS):
    m, c, filled = D.STYLE[label]
    msz = ms * (1.35 if m in ('*', 'P', 'X') else 1.0)
    h = ax.errorbar(x, y, yerr=e, xerr=xerr, fmt=m, color=c, mfc=c if filled else 'white', mec=c, ms=msz, mew=0.7,
                    elinewidth=0.6, capsize=0, zorder=zorder, label=label)
    if label not in reg:
        reg[label] = Line2D([], [], ls='none', marker=m, color=c, mfc=c if filled else 'white', mec=c, ms=msz, mew=0.7)
    return h


def data_handles(reg):
    labs = [l for l in D.ORDER if l in reg] + [l for l in reg if l not in D.ORDER]
    return [reg[l] for l in labs], labs


def save(fig, name, exact=False):
    """exact=True: keep the exact figure size (no cropping) -- the paper includes these PDFs without scaling"""
    with plt.rc_context({'savefig.bbox': 'standard'} if exact else {}):
        fig.savefig(os.path.join(C.FIGDIR, name + '.pdf'))
        fig.savefig(os.path.join(C.PNG, name + '.png'), dpi=200)
    plt.close(fig)
    print('  wrote', name)


def legend_rows(fig, th, tl, dh, dl, ncol_t=None, ncol_d=None):
    """theory legend in a row above the panels, data legend in a row (or rows) below (constrained layout)"""
    fig.legend(th, tl, loc='outside upper center', ncol=ncol_t or len(tl), handlelength=2.0, columnspacing=1.0)
    if dh:
        fig.legend(dh, dl, loc='outside lower center', ncol=ncol_d or len(dl), handlelength=1.0, handletextpad=0.35,
                   columnspacing=1.1)


def legend_below(fig, th, tl, dh, dl, ncol_t=3, ncol_d=6, y_t=0.0, y_d=None, fs=7.0):
    lt = fig.legend(th, tl, loc='upper center', ncol=ncol_t, bbox_to_anchor=(0.5, y_t), fontsize=fs,
                    handlelength=2.2, columnspacing=1.2, borderaxespad=0.)
    fig.add_artist(lt)
    if dh:
        fig.canvas.draw()
        bb = lt.get_window_extent().transformed(fig.transFigure.inverted())
        fig.legend(dh, dl, loc='upper center', ncol=ncol_d, bbox_to_anchor=(0.5, bb.y0 - 0.012), fontsize=fs,
                   handlelength=1.0, columnspacing=1.0, handletextpad=0.3, borderaxespad=0.)


# ============================================================== figures
def fig1(df, cv, ix):
    fig, axs = plt.subplots(2, 3, figsize=(W2, 4.7), sharex=True, layout='constrained')
    reg = {}
    cs = df[df.is_cs & (df.dataset != 'Hallin1993_TabI')]
    for k, (fr, th0) in enumerate(FIG1):
        ax = axs.ravel()[k]
        idx = np.array([ix[gid('f1', fr, th0, '%.2f' % E)] for E in E_GRID])
        draw_theory(ax, E_GRID, cv, idx, extra=True, band_xmax=E_BAND_MAX)
        col = 'theta_lab' if fr == 'LAB' else 'theta_cm'
        sel = cs[(cs.frame == fr) & (np.abs(cs[col] - th0) <= C.WIN_THETA + 1e-6)]
        mark_shown('Fig. 1 (dsigma vs E)', sel)
        for lab in [l for l in D.ORDER if l in set(sel.label)]:
            s = sel[sel.label == lab]
            draw_points(ax, lab, s.E.values, s.value.values, s.err.values, reg)
        frl = 'lab' if fr == 'LAB' else 'cm'
        ax.text(0.05, 0.93, r'$%g^\circ$ %s' % (th0, 'lab' if fr == 'LAB' else 'c.m.'), transform=ax.transAxes,
                va='top', ha='left')
        ax.axvline(E_2PI, color='0.45', lw=0.5, ls=':', zorder=0.5)       # two-pion (pi0 pi0 p) threshold
        mb = E_GRID <= E_BAND_MAX
        ymax = max(np.max(cv['cpl_hi'][idx][mb]), np.max(cv['main'][idx]), np.max(cv['fit'][idx]), np.max(cv['tree'][idx]),
                   (sel.value + sel.err).max() if len(sel) else 0)
        ax.set_ylim(0, 1.08 * ymax)
        ax.set_xlim(0, 800)
        if k % 3 == 0:
            ax.set_ylabel(r'$d\sigma/d\Omega_{\rm %s}$ [nb/sr]' % frl)
        if k >= 3:
            ax.set_xlabel(r'$E_\gamma$ [MeV]')
        ax.xaxis.set_major_locator(plt.MultipleLocator(200))
    th, tl = theory_handles(extra=has_extra(cv))
    dh, dl = data_handles(reg)
    legend_rows(fig, th, tl, dh, dl, ncol_d=9)
    save(fig, 'fig_delta_cs_vsE', exact=True)
    return reg


def fig2(df, cv, ix, energies, name):
    """all angular distributions in one grid (parts (a) and (b) until 2026-10-18), equally spaced rows; the legend fills
    the free slot(s) of the last row, under the x labels of the panel above (user, 2026-10-18)"""
    nc = 4
    nr = int(np.ceil(len(energies) / nc))
    assert len(energies) < nr * nc, 'fig2: the legend needs a free slot in the last row'
    fig, axs = plt.subplots(nr, nc, figsize=(W2, 1.22 * nr + 0.3), sharex=True, layout='constrained')
    # rect: room on the right for the half of the '180' tick label of the panel above the free slot (labelled after the layout)
    fig.get_layout_engine().set(h_pad=0.03, w_pad=0.03, hspace=0.02, wspace=0.02, rect=(0, 0, 0.988, 1))
    axs = axs.ravel()
    reg = {}
    cs = df[df.is_cs & (df.dataset != 'Hallin1993_TabI')].copy()
    # each point assigned to the nearest panel energy of ALL Fig. 2 panels (no point shown twice)
    Ep = np.array(FIG2)
    near = Ep[np.argmin(np.abs(cs.E.values[:, None] - Ep[None, :]), axis=1)]
    cs['Epanel'] = near
    cs = cs[np.abs(cs.E - cs.Epanel) <= 3.0 + 1e-6]
    conv_labels = set()
    for k, E0 in enumerate(energies):
        ax = axs[k]
        idx = np.array([ix[gid('f2', E0, th)] for th in TH_GRID])
        draw_theory(ax, TH_GRID, cv, idx, extra=C.EXTRA_THETA)     # the extra curve also here (default for 'run3')
        sel = cs[np.isclose(cs.Epanel, E0)]
        mark_shown('Fig. 2 (dsigma vs theta)', sel)
        th, v, e, conv = D.as_frame(sel, 'CM')
        sel = sel.assign(thp=th, vp=v, ep=e, conv=conv)
        conv_labels |= set(sel[sel.conv].label)
        for lab in [l for l in D.ORDER if l in set(sel.label)]:
            s = sel[sel.label == lab]
            draw_points(ax, lab, s.thp.values, s.vp.values, s.ep.values, reg)
        draw_forward(ax, forward_above('CS', E0, 3.0), reg, cm=True)
        ax.text(0.05, 0.93, r'%g MeV' % E0, transform=ax.transAxes, va='top', ha='left')
        fwd_pts = forward_above('CS', E0, 3.0)
        ymax = max(np.max(cv['cpl_hi'][idx]), np.max(cv['hi'][idx]), np.max(cv['tree'][idx]),
                   (sel.vp + sel.ep).max() if len(sel) else 0,
                   max([(v + e) / D.jac_cm_to_lab(E, 0.0) for E, v, e in fwd_pts], default=0))
        ax.set_ylim(0, 1.17 * ymax)
        ax.set_xlim(-6 if FWD_DF is not None else 0, 180)
        ax.xaxis.set_major_locator(plt.MultipleLocator(60))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(15))
        if k % nc == 0:
            ax.set_ylabel(r'$d\sigma/d\Omega_{\rm cm}$ [nb/sr]')
    for k in range(len(energies), len(axs)):
        axs[k].axis('off')
    bottom = [k for k in range(len(energies)) if k >= len(energies) - nc]      # lowest panel of each column
    above = [k for k in bottom if k + nc < len(axs)]                            # ... with a free slot below it
    for k in bottom:
        if k not in above:
            axs[k].set_xlabel(r'$\theta_{\rm cm}$ [deg]')
            axs[k].xaxis.set_tick_params(labelbottom=True)
    # freeze the layout before labelling the panels above the free slots, so that all rows stay equally spaced;
    # their x labels then occupy the top of the free slot, and the legend goes below them
    fig.get_layout_engine().execute(fig)
    fig.set_layout_engine('none')
    for k in above:
        axs[k].set_xlabel(r'$\theta_{\rm cm}$ [deg]')
        axs[k].xaxis.set_tick_params(labelbottom=True)
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    r = fig.canvas.get_renderer()
    y_top = min(inv.transform(axs[k].get_tightbbox(r).get_points())[0, 1] for k in above)
    x_left = inv.transform(axs[len(energies) - 1].get_tightbbox(r).get_points())[1, 0]   # right end of the last panel
    th, tl = theory_handles(extra=C.EXTRA_THETA and has_extra(cv))
    dh, dl = data_handles(reg)
    # theory in one column, the data sets in two columns below it (one slot is too narrow for 2 x (curve + label))
    kw = dict(fontsize=7.0, labelspacing=0.25, borderaxespad=0.0, loc='upper left')
    lt = fig.legend(th, tl, bbox_to_anchor=(x_left, y_top - 0.012), handlelength=1.8, handletextpad=0.4, **kw)
    fig.canvas.draw()
    bt = inv.transform(lt.get_window_extent(r).get_points())
    ld = fig.legend(dh, dl, bbox_to_anchor=(x_left, bt[0, 1] - 0.006), ncol=2, handlelength=1.0, handletextpad=0.3,
                    columnspacing=0.8, **kw)
    fig.canvas.draw()
    bd = inv.transform(ld.get_window_extent(r).get_points())
    # common left edge: the pair centred under the panel above, but inside the free region
    w = max(bt[1, 0] - bt[0, 0], bd[1, 0] - bd[0, 0])
    pa = axs[above[-1]].get_position()
    x0 = min(max((pa.x0 + pa.x1 - w) / 2, x_left), 1.0 - w)
    lt.set_bbox_to_anchor((x0, y_top - 0.012), transform=fig.transFigure)
    ld.set_bbox_to_anchor((x0, bt[0, 1] - 0.006), transform=fig.transFigure)
    fig.canvas.draw()
    lb = np.r_[inv.transform(lt.get_window_extent(r).get_points()), inv.transform(ld.get_window_extent(r).get_points())]
    print('  %s: legend box x %.3f-%.3f (free %.3f-1), y %.3f-%.3f (free 0-%.3f)' % (
        name, lb[:, 0].min(), lb[:, 0].max(), x_left, lb[:, 1].min(), lb[:, 1].max(), y_top))
    assert lb[:, 0].min() >= x_left - 1e-3 and lb[:, 0].max() <= 1.0 + 1e-3 and lb[:, 1].min() >= 0.0, 'fig2: legend overflows'
    save(fig, name, exact=True)
    return conv_labels


def fig3a(df, cv, ix):
    fig, axs = plt.subplots(3, 4, figsize=(W2, 4.45), sharex=True, layout='constrained')
    fig.get_layout_engine().set(h_pad=0.03, w_pad=0.03, hspace=0.02, wspace=0.02)
    axs = axs.ravel()
    reg = {}
    s3 = df[df.obs == 'S3']
    for k, E0 in enumerate(LEGS_E):
        ax = axs[k]
        idx = np.array([ix[gid('f3', E0, th)] for th in TH_GRID])
        draw_theory(ax, TH_GRID, cv, idx, extra=True)
        sel = s3[np.abs(s3.E - E0) <= 0.05]
        mark_shown('Fig. 3a (Sigma3 vs theta)', sel)
        for lab in [l for l in D.ORDER if l in set(sel.label)]:
            s = sel[sel.label == lab]
            draw_points(ax, lab, s.theta_cm.values, s.value.values, s.err.values, reg)
        ax.axhline(0, color='0.6', lw=0.4, zorder=0)
        ax.text(0.05, 0.94, r'%g' % E0, transform=ax.transAxes, va='top', ha='left')    # energy in MeV (units in the caption)
        ax.set_xlim(0, 180)
        ax.xaxis.set_major_locator(plt.MultipleLocator(60))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(15))
        if k % 4 == 0:
            ax.set_ylabel(r'$\Sigma_3$')
        if k >= 8:
            ax.set_xlabel(r'$\theta_{\rm cm}$ [deg]')
    # common y range per row, with headroom for the energy labels
    for r, (lo, hi) in enumerate([(-0.40, 0.46), (-0.25, 0.42), (-0.05, 0.60)]):    # third row 0.80 -> 0.60 (2026-10-17)
        for a in axs[4 * r:4 * r + 4]:
            a.set_ylim(min(lo, a.get_ylim()[0]), hi)
    th, tl = theory_handles(extra=has_extra(cv))
    dh, dl = data_handles(reg)
    hh, ll = th + dh, tl + dl                       # 3 columns, filled column-wise -> rows: (entries 0, 1, 2) / (3, 4, 5),
    o = [0, 3, 1, 4, 2, 5] if len(hh) == 6 else list(range(len(hh)))   # i.e. (main, NNLO fit, running) / (tree, g_M, LEGS)
    fig.legend([hh[i] for i in o], [ll[i] for i in o], loc='outside upper center', ncol=3, handlelength=2.0,
               columnspacing=2.0)
    save(fig, 'fig_delta_S3_vs_theta', exact=True)


def fig3b(df, cv, ix):
    fig, ax = plt.subplots(figsize=(W1, 2.45), layout='constrained')
    reg = {}
    idx = np.array([ix[gid('f3b', '%.2f' % E)] for E in E3B])
    draw_theory(ax, E3B, cv, idx, extra=True, band_xmax=E_BAND_MAX)
    s3 = df[(df.obs == 'S3') & (np.abs(df.theta_cm - 90.0) <= C.WIN_THETA + 1e-6)]
    mark_shown('Fig. 3b (Sigma3(90) vs E)', s3)
    for lab in [l for l in D.ORDER if l in set(s3.label)]:
        s = s3[s3.label == lab]
        xerr = None
        if lab == 'Frascati 68':
            xerr = np.array([[s.E.values[0] - s.Elo.values[0]], [s.Ehi.values[0] - s.E.values[0]]])
        draw_points(ax, lab, s.E.values, s.value.values, s.err.values, reg, xerr=xerr)
    ax.axhline(0, color='0.6', lw=0.4, zorder=0)
    ax.set_xlim(0, 800); ax.set_ylim(-1.05, 0.75)
    ax.xaxis.set_major_locator(plt.MultipleLocator(200))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel(r'$E_\gamma$ [MeV]'); ax.set_ylabel(r'$\Sigma_3(\theta_{\rm cm}=90^\circ)$')
    th, tl = theory_handles(extra=has_extra(cv))
    dh, dl = data_handles(reg)
    # one legend (curves and data) inside the axes, lower centre-right: free of the curves (user, 2026-10-17)
    ax.legend(th + dh, tl + dl, loc='lower center', bbox_to_anchor=(0.62, 0.02), ncol=2, handlelength=1.8, handletextpad=0.5,
              columnspacing=1.0, borderaxespad=0.3, fontsize=7.5)
    save(fig, 'fig_delta_S3_90_vsE', exact=True)
    return s3


def fig4(df, cv, ix):
    # single-column figure (2x2 panels, rows share the y axis)
    fig, axs = plt.subplots(2, 2, figsize=(W1, 4.05), sharex=True, sharey='row', layout='constrained')
    axs = axs.ravel()
    reg = {}
    for k, (o, E0, title) in enumerate(FIG4):
        ax = axs[k]
        idx = np.array([ix[gid('f4', o, E0, th)] for th in TH_GRID])
        draw_theory(ax, TH_GRID, cv, idx, extra=True)
        # bin-averaged NNLO (bremsstrahlung weight 1/E): <N>/<D>; only reported (differs by <= 0.01)
        lo, hi = BIN[E0]
        Eb = np.linspace(lo, hi, 7)
        w = 1.0 / Eb; w /= w.sum()
        avg = np.array([np.sum([w[j] * cv['num'][ix[gid('f4b', o, E0, '%.3f' % Eb[j], th)]] for j in range(7)]) /
                        np.sum([w[j] * cv['den'][ix[gid('f4b', o, E0, '%.3f' % Eb[j], th)]] for j in range(7)])
                        for th in TH_GRID])
        cv.setdefault('binavg', {})[(o, E0)] = avg
        sel = df[(df.obs == o) & np.isclose(df.E, E0)]
        mark_shown('Fig. 4 (Sigma2x,2z vs theta)', sel)
        for lab in [l for l in D.ORDER if l in set(sel.label)]:
            s = sel[sel.label == lab]
            draw_points(ax, lab, s.theta_lab.values, s.value.values, s.err.values, reg, ms=3.8)
        if o == 'S2Z':
            lo_, hi_ = BIN[E0]
            draw_forward(ax, forward_above('S2Z', 0.5 * (lo_ + hi_), 0.5 * (hi_ - lo_)), reg)
        if o == 'S2Z':
            # Paudyal 10.4% correlated scale uncertainty: grey bars at the bottom
            yb = -0.09
            for t, v in zip(sel.theta_lab, sel.value):
                ax.add_patch(Rectangle((t - 2.5, yb), 5.0, 0.104 * abs(v), color='0.6', lw=0, zorder=2))
        ax.axhline(0, color='0.6', lw=0.4, zorder=0)
        ax.set_xlim(-6 if FWD_DF is not None else 0, 180)
        ax.xaxis.set_major_locator(plt.MultipleLocator(60))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(15))
        sym = r'\Sigma_{2x}' if o == 'S2X' else r'\Sigma_{2z}'
        if k % 2 == 0:
            ax.set_ylabel(r'$%s$' % sym)
        if o == 'S2X':
            ax.text(0.04, 0.07, r'%s MeV' % ('%g' % round(E0, 1)), transform=ax.transAxes, va='bottom', ha='left')
        else:
            ax.text(0.5, 0.93, r'%s MeV' % ('%g' % round(E0, 1)), transform=ax.transAxes, va='top', ha='center')
        if k >= 2:
            ax.set_xlabel(r'$\theta_{\rm lab}$ [deg]')
    axs[0].set_ylim(-0.75, 0.75); axs[1].set_ylim(-0.75, 0.75)
    axs[2].set_ylim(-0.1, 1.05); axs[3].set_ylim(-0.1, 1.05)
    th, tl = theory_handles(extra=has_extra(cv))
    dh, dl = data_handles(reg)
    dh.append(Patch(facecolor='0.6', edgecolor='none')); dl.append(r'A2 20 scale unc.\ (10.4\%)')
    legend_rows(fig, th, tl, dh, dl, ncol_t=2, ncol_d=2)
    save(fig, 'fig_delta_S2x_S2z', exact=True)


def fig5(df, cv, ix):
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.4), sharex=True)
    axs = axs.ravel()
    reg = {}
    for k, (o, th0) in enumerate(FIG5):
        ax = axs[k]
        idx = np.array([ix[gid('f5', o, th0, '%.2f' % E)] for E in E5])
        draw_theory(ax, E5, cv, idx, extra=True)
        sel = df[(df.obs == o) & (np.abs(df.theta_lab - th0) <= C.WIN_THETA + 1e-6)]
        mark_shown('Fig. 5 (Sigma2x,2z vs E)', sel)
        for lab in [l for l in D.ORDER if l in set(sel.label)]:
            s = sel[sel.label == lab]
            xerr = np.vstack([s.E.values - s.Elo.values, s.Ehi.values - s.E.values])
            draw_points(ax, lab, s.E.values, s.value.values, s.err.values, reg, xerr=xerr, ms=3.8)
        ax.axhline(0, color='0.6', lw=0.4, zorder=0)
        ax.set_xlim(50, 400)
        sym = r'\Sigma_{2x}' if o == 'S2X' else r'\Sigma_{2z}'
        ax.set_ylabel(r'$%s(\theta_{\rm lab}=%g^\circ)$' % (sym, th0))
        if k >= 2:
            ax.set_xlabel(r'$E_\gamma$ [MeV]')
    fig.tight_layout(rect=(0, 0.1, 1, 1), h_pad=0.5, w_pad=1.0)
    th, tl = theory_handles(extra=has_extra(cv))
    dh, dl = data_handles(reg)
    fig.legend(th + dh, tl + dl, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.092), fontsize=7.0,
               handlelength=2.2, columnspacing=1.2, borderaxespad=0.)
    save(fig, 'fig_delta_S2_vsE')


def fig6(df, cv, ix):
    """data / NNLO at the exact kinematics of every Delta-region cross-section point (E > 150 MeV)"""
    cs = df[df.is_cs & (df.E > C.E_DELTA_MIN) & (df.dataset != 'Hallin1993_TabI')].copy()
    ii = np.array([ix[p] for p in cs.pid])
    th = cv['main'][ii]
    cs['r'] = cs.value / th; cs['re'] = cs.err / th
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.7), sharex=True, sharey=True)
    axs = axs.ravel()
    reg = {}
    for k, (a, b, thr) in enumerate(FIG6):
        ax = axs[k]
        s = cs[(cs.theta_cm >= a) & (cs.theta_cm < b)].sort_values('E')
        mark_shown('Fig. 6 (dsigma ratio)', s)
        idx = np.array([ix[gid('f6', thr, '%.1f' % E)] for E in E6])
        m0 = cv['main'][idx]
        ax.fill_between(E6, cv['cpl_lo'][idx] / m0, cv['cpl_hi'][idx] / m0, color=C_GM, lw=0, zorder=1)
        ax.plot(E6, cv['cpl_lo'][idx] / m0, '-', color=C_GM_EDGE, lw=0.35, zorder=1.5)
        ax.plot(E6, cv['cpl_hi'][idx] / m0, '-', color=C_GM_EDGE, lw=0.35, zorder=1.5)
        ax.fill_between(E6, cv['lo'][idx] / m0, cv['hi'][idx] / m0, color=C_BAND, lw=0, alpha=0.85, zorder=2)
        ax.axhline(1.0, color=C_MAIN, lw=1.0, zorder=3)
        for lab in [l for l in D.ORDER if l in set(s.label)]:
            q = s[s.label == lab]
            draw_points(ax, lab, q.E.values, q.r.values, q.re.values, reg, ms=2.8)
        lo_ = r'%d^\circ\leq' % a if a > 0 else ''
        hi_ = r'<%d^\circ' % b if b < 181 else ''
        ax.set_title(r'$%s\theta_{\rm cm}%s$  (bands: $\theta_{\rm cm}=%g^\circ$)' % (lo_, hi_, thr), fontsize=8, pad=3)
        ax.set_xlim(150, 400); ax.set_ylim(0.45, 1.75)
        if k % 2 == 0:
            ax.set_ylabel(r'data / NNLO')
        if k >= 2:
            ax.set_xlabel(r'$E_\gamma$ [MeV]')
    fig.tight_layout(rect=(0, 0.13, 1, 1), h_pad=0.6, w_pad=0.6)
    th = [(Patch(facecolor=C_BAND, edgecolor='none'), Line2D([], [], color=C_MAIN, lw=1.0)),
          Patch(facecolor=C_GM, edgecolor=C_GM_EDGE, lw=0.4)]
    tl = [(r'NNLO B$\chi$PT prediction (=1) and its uncertainty' if C.MAIN == 'pred' else
           r'NNLO B$\chi$PT, fitted $\alpha_{E1},\beta_{M1}$ (=1) and its uncertainty'), LAB_GM]
    dh, dl = data_handles(reg)
    legend_below(fig, th, tl, dh, dl, ncol_t=2, ncol_d=7, y_t=0.115, y_d=0.075)
    save(fig, 'fig_delta_cs_ratio')
    return cs


# ============================================================== chi^2 table
def chi2_rows(df, cv, ix, Emin, Emax, curve='main'):
    s = df[(df.E > Emin) & (df.E <= Emax) & df.obs.isin(['CS', 'CSCM', 'S3', 'S2X', 'S2Z'])].copy()
    s['set'] = s.dataset.replace({'Hallin1993_TabII': 'Hallin1993_TabII(D)', 'Hallin1993_TabII_D': 'Hallin1993_TabII(D)'})
    rows = []
    for ds in dict.fromkeys(s.set):
        q = s[s.set == ds]
        ii = np.array([ix[p] for p in q.pid])
        t = cv[curve][ii]; y = q.value.values; e = q.err.values
        sth = 0.5 * (cv['hi'][ii] - cv['lo'][ii])
        dn = float(q.norm.iloc[0])
        c0 = float(np.sum(((y - t) / e) ** 2))

        def nfit(err):
            if dn <= 0:
                return 1.0, float(np.sum(((y - t) / err) ** 2))
            a = np.sum(t * t / err ** 2) + 1 / dn ** 2
            b = np.sum(t * y / err ** 2) + 1 / dn ** 2
            n = b / a
            return n, float(np.sum(((y - n * t) / err) ** 2) + ((n - 1) / dn) ** 2)
        n1, c1 = nfit(e)
        n2, c2 = nfit(np.sqrt(e ** 2 + sth ** 2))
        rows.append(dict(set=ds, label=q.label.iloc[0], obs='/'.join(sorted(set(q.obs))), N=len(q),
                         Erange=('%.0f' % q.E.min()) if q.E.max() - q.E.min() < 0.5 else '%.0f-%.0f' % (q.E.min(), q.E.max()), norm=dn, in_fitdb=bool(q.in_fitdb.iloc[0]),
                         published=bool(q.published.iloc[0]), chi2_n1=c0, chi2N_n1=c0 / len(q),
                         n=n1, chi2_nfit=c1, chi2N_nfit=c1 / len(q), chi2N_th=c2 / len(q), n_th=n2))
    return pd.DataFrame(rows)


SETNAME = {
    'A2_Martel2015_S2x': 'A2 15 $\\Sigma_{2x}$ (Martel), 288 MeV', 'Martel_thesis_S2x': 'Martel thesis $\\Sigma_{2x}$, 331 MeV (unpubl.)',
    'A2_Paudyal2020_S2z': 'A2 20 $\\Sigma_{2z}$ (Paudyal)', 'LEGS2001_S3': 'LEGS 01 $\\Sigma_3$', 'LEGS2001_CS': 'LEGS 01 $d\\sigma$',
    'Hunger1997_75': 'MAMI 97 75° c.m. (Hünger)', 'Hunger1997_90': 'MAMI 97 90° c.m. (Hünger)', 'Wissmann1999': 'MAMI 99 (Wissmann)',
    'Camen2002': 'MAMI 02 (Camen/Kossert)', 'Hallin1993_TabII(D)': 'SAL 93 (Hallin) Tab. II', 'Hallin1993_TabI': 'SAL 93 (Hallin) Tab. I [not indep.]',
    'OlmosdeLeon2001': 'A2 01 (Olmos de León)', 'Genzel1972': 'Bonn 72 (Genzel)', 'Gray1967': 'Illinois 67 (Gray)',
    'Baranov1966': 'Moscow 66 (Baranov)', 'Baranov1966_refined': 'Moscow 66 refined, fwd pts', 'Kolb2000': 'SAL 00 (Kolb)',
    'Nagashima1964': 'Tokyo 64 (Nagashima)', 'Barbiellini1968': 'Frascati 68 $\\Sigma_3$ (sign flipped)',
    'MacGibbon1995_untagged': 'SAL 95 (MacGibbon)'}
SETORDER = ['A2_Martel2015_S2x', 'Martel_thesis_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_S3', 'Barbiellini1968', 'LEGS2001_CS',
            'Hunger1997_75', 'Hunger1997_90', 'Wissmann1999', 'Camen2002', 'Kolb2000', 'Hallin1993_TabII(D)', 'Hallin1993_TabI',
            'OlmosdeLeon2001', 'Genzel1972', 'Gray1967', 'Baranov1966', 'Baranov1966_refined', 'Nagashima1964']


# rows of Table 'delta' of the paper (Appendix G), in its order and with its labels (\cite keys as in main.tex)
PAPERROWS = [('A2_Martel2015_S2x', r'A2 $\Si_{2x}$~\cite{A2:2014iky}                  '),
             ('A2_Paudyal2020_S2z', r'A2 $\Si_{2z}$~\cite{A2:2019bqm}'),
             ('LEGS2001_S3', r'LEGS $\Si_3$~\cite{Blanpied:2001ae}'),
             ('LEGS2001_CS', r'LEGS $\dd\si/\dd\W$~\cite{Blanpied:2001ae}'),
             ('Hunger1997_75', r'MAMI $75^\circ$~\cite{Hunger:1997pt}'),
             ('Hunger1997_90', r'MAMI $90^\circ$~\cite{Hunger:1997pt}'),
             ('Wissmann1999', r'MAMI $131^\circ$~\cite{Wissmann:1999vi}'),
             ('Camen2002', r'MAMI $136^\circ$~\cite{Camen:2001st}'),
             ('Hallin1993_TabII(D)', r'SAL~\cite{Hallin:1993ft}'),
             ('OlmosdeLeon2001', r'A2 (2001), 157--164\,MeV~\cite{deLeon:2001dnx}')]


def _pf(v):
    """number format of Table delta in main.tex: two decimals, one decimal from 10 on"""
    return ('%.1f' % v) if v >= 9.995 else ('%.2f' % v)


def chi2_table(df, cv, ix):
    """chi^2/N per data set.  Columns: the main curve (config.MAIN; its band for 'incl. th.'), the other of the two
    (prediction/fit) with n = 1 and n fitted, tree level, without vertex corrections."""
    tabs = {}
    for key, (a, b) in {'350': (C.E_DELTA_MIN, C.E_DELTA_MAX), '400': (C.E_DELTA_MIN, 400.0)}.items():
        t = chi2_rows(df, cv, ix, a, b)
        tf = chi2_rows(df, cv, ix, a, b, curve='fit')
        t['chi2N_fit_n1'] = tf.chi2N_n1.values
        t['chi2N_fit_nfit'] = tf.chi2N_nfit.values
        t['n_fit'] = tf.n.values
        t['chi2N_fit_th'] = tf.chi2N_th.values         # NNLO-fit curve with the band of the MAIN curve (reference only)
        tt = chi2_rows(df, cv, ix, a, b, curve='tree')
        t['chi2N_tree_nfit'] = tt.chi2N_nfit.values
        if 'noVC' in cv:
            t['chi2N_noVC_nfit'] = chi2_rows(df, cv, ix, a, b, curve='noVC').chi2N_nfit.values
        if has_extra(cv) and EXTRA_KEY != 'noVC':
            # running Delta parameters: n = 1 and n fitted, point errors only (the paper's '+run.' column, with the band
            # recomputed for the running amplitude, is results/running/chi2_table_irun3.csv of chi2_running.py)
            tr = chi2_rows(df, cv, ix, a, b, curve=EXTRA_KEY)
            t['chi2N_%s_n1' % EXTRA_KEY] = tr.chi2N_n1.values
            t['chi2N_%s_nfit' % EXTRA_KEY] = tr.chi2N_nfit.values
            t['n_%s' % EXTRA_KEY] = tr.n.values
        t['order'] = [SETORDER.index(x) if x in SETORDER else 99 for x in t.set]
        t = t.sort_values('order').drop(columns='order')
        tabs[key] = t
        t.to_csv(os.path.join(C.RESULTS, 'chi2_table_E%s.csv' % key), index=False)
    mainlab = 'NNLO prediction' if C.MAIN == 'pred' else 'NNLO fit'
    fitlab = 'NNLO fit' if C.MAIN == 'pred' else 'NNLO prediction'
    # markdown
    L = ['main curve: %s (band: %s)' % (mainlab, 'LMP15-type envelope, nothing fitted' if C.MAIN == 'pred' else
                                           'envelope without static alpha, beta (+) fit error'),
         'other curve: %s (alpha, beta = %.3f, %.3f%s)' % (fitlab, C.ALPHA_FIT, C.BETA_FIT,
                                                           (', from ' + os.path.relpath(C.FITS_JSON, C.ROOT)) if C.FITS_JSON else '')]
    for key, t in tabs.items():
        L.append('\n### %s < E_lab <= %s MeV\n' % (int(C.E_DELTA_MIN), key))
        hdr = ('| data set | obs | N | E range | scale δ | %s: χ²/N (n=1) | χ²/N (n fitted) | n | χ²/N incl. th. band | '
               '%s: χ²/N (n=1) | χ²/N (n fitted) | n | χ²/N Born+π⁰+Δ |' % (mainlab, fitlab))
        if 'chi2N_noVC_nfit' in t:
            hdr += ' χ²/N without VC |'
        runcol = 'chi2N_%s_n1' % EXTRA_KEY if EXTRA_KEY and ('chi2N_%s_n1' % EXTRA_KEY) in t else None
        if runcol:
            hdr += ' running Δ (irun=%d): χ²/N (n=1) | χ²/N (n fitted) | n |' % C.EXTRA_IRUN
        hdr += ' |'
        L.append(hdr)
        L.append('|' + '---|' * (hdr.count('|') - 1))
        for r in t.itertuples():
            nm = SETNAME.get(r.set, r.set).replace('$\\Sigma_{2x}$', 'Σ2x').replace('$\\Sigma_{2z}$', 'Σ2z').replace('$\\Sigma_3$', 'Σ3').replace('$d\\sigma$', 'dσ')
            sc = ('%.1f%%' % (100 * r.norm)) if r.norm > 0 else '–'
            nn = ('%.3f' % r.n) if r.norm > 0 else '1'
            c1 = ('%.2f' % r.chi2N_nfit) if r.norm > 0 else '–'
            nf = ('%.3f' % r.n_fit) if r.norm > 0 else '1'
            cf = ('%.2f' % r.chi2N_fit_nfit) if r.norm > 0 else '–'
            line = '| %s | %s | %d | %s | %s | %.2f | %s | %s | %.2f | %.2f | %s | %s | %.2f |' % (
                nm, r.obs, r.N, r.Erange, sc, r.chi2N_n1, c1, nn, r.chi2N_th, r.chi2N_fit_n1, cf, nf, r.chi2N_tree_nfit)
            if 'chi2N_noVC_nfit' in t:
                line += ' %.2f |' % r.chi2N_noVC_nfit
            if runcol:
                rr = getattr(r, 'chi2N_%s_nfit' % EXTRA_KEY); rn = getattr(r, 'n_%s' % EXTRA_KEY)
                line += ' %.2f | %s | %s |' % (getattr(r, runcol), ('%.2f' % rr) if r.norm > 0 else '–', ('%.3f' % rn) if r.norm > 0 else '1')
            fr = ''
            line += ' %s |' % fr
            L.append(line)
    open(os.path.join(C.RESULTS, 'chi2_table.md'), 'w').write('\n'.join(L) + '\n')
    # LaTeX (main range): all sets, main curve
    t = tabs['350']
    T = ['%% main curve: %s' % mainlab, '\\begin{tabular}{lrrcccc}', '\\hline',
         'Data set & $N$ & $E_\\gamma$ [MeV] & $\\chi^2/N$ ($n=1$) & $\\chi^2/N$ ($n$ fitted) & $n$ & $\\chi^2/N$ (incl.\\ th.)\\\\', '\\hline']
    for r in t.itertuples():
        nm = unicodedata.normalize('NFC', SETNAME.get(r.set, r.set)).replace('°', '$^\\circ$').replace('ü', '\\"u').replace('é', "\\'e").replace('ó', "\\'o")
        c1 = ('%.2f' % r.chi2N_nfit) if r.norm > 0 else '--'
        nn = ('%.2f' % r.n) if r.norm > 0 else '--'
        T.append('%s & %d & %s & %.2f & %s & %s & %.2f\\\\' % (nm, r.N, r.Erange.replace('-', '--'), r.chi2N_n1, c1, nn, r.chi2N_th))
    T += ['\\hline', '\\end{tabular}']
    open(os.path.join(C.RESULTS, 'chi2_table.tex'), 'w').write('\n'.join(T) + '\n')
    # LaTeX rows in the format of Table 'delta' of the paper (Appendix G), main curve and the other curve
    tt = t.set_index('set')
    P = ['%% rows of Table delta (paper Appendix G) in the format of main.tex; generated by work/delta_plots/make_figs.py',
         '%% columns: Data set & N & chi^2/N (n = 1) & chi^2/N (n fitted) & chi^2/N (n fitted, theory band in quadrature)',
         '%% --- main curve: %s (band: %s) ---' % (mainlab, 'LMP15-type envelope' if C.MAIN == 'pred' else 'envelope (+) fit error')]
    for s, lab in PAPERROWS:
        r = tt.loc[s]
        P.append('%s & %d & %s & %s & %s\\\\' % (lab, r.N, _pf(r.chi2N_n1), _pf(r.chi2N_nfit) if r.norm > 0 else '--', _pf(r.chi2N_th)))
    P.append('%% --- other curve: %s, n = 1 and n fitted (no own band; last column: with the band of the main curve) ---' % fitlab)
    for s, lab in PAPERROWS:
        r = tt.loc[s]
        P.append('%s & %d & %s & %s & %s\\\\' % (lab, r.N, _pf(r.chi2N_fit_n1), _pf(r.chi2N_fit_nfit) if r.norm > 0 else '--', _pf(r.chi2N_fit_th)))
    open(os.path.join(C.RESULTS, 'chi2_table_paper.tex'), 'w').write('\n'.join(P) + '\n')
    return tabs


def provenance(df):
    figs = list(SHOWN)
    L = ['| data set | label | obs (frame of angle) | N | E range [MeV] | angles | point error | scale | in fitdb | ' + ' | '.join(
        f.split(' (')[0] for f in figs) + ' |', '|' + '---|' * (9 + len(figs))]
    for ds, q in df[df.obs.isin(['CS', 'CSCM', 'S3', 'S2X', 'S2Z'])].groupby('dataset', sort=False):
        th = q.theta_lab if q.frame.iloc[0] == 'LAB' else q.theta_cm
        ang = sorted(set(np.round(th, 1)))
        angs = ('%g' % ang[0]) if len(ang) == 1 else ('%g-%g (%d)' % (ang[0], ang[-1], len(ang)))
        er = ('%.0f' % q.E.min()) if q.E.max() - q.E.min() < 0.5 else '%.0f-%.0f' % (q.E.min(), q.E.max())
        cols = [str(SHOWN[f].get(ds, '')) for f in figs]
        L.append('| %s | %s | %s (%s) | %d | %s | %s | %s | %s | %s | %s |' % (
            ds, q.label.iloc[0].replace('$\\gamma$', 'γ'), q.obs.iloc[0], q.frame.iloc[0], len(q), er, angs,
            q.note.iloc[0].split(';')[0][:60], ('%.1f%%' % (100 * q.norm.iloc[0])) if q.norm.iloc[0] > 0 else '–',
            'yes' if q.in_fitdb.iloc[0] else 'no', ' | '.join(cols)))
    open(os.path.join(C.RESULTS, 'provenance.md'), 'w').write('\n'.join(L) + '\n')


CITE = {'LEGS 01': 'Blanpied:2001ae', 'MAMI 97': 'Hunger:1997pt', 'MAMI 99': 'Wissmann:1999vi',
        'MAMI 02': 'Camen:2001st,Kossert:2002ws', 'SAL 93': 'Hallin:1993ft', 'SAL 95': 'MacGibbon:1995in',
        'SAL 00': 'Kolb:2000ix', 'Bonn 72': 'Genzel:1972ba', 'Illinois 67': 'Gray:1967zzd',
        'Moscow 66': 'Baranov:1966aaa,Baranov:1966oqa', 'Tokyo 64': 'Nagashima:1964zz', 'Frascati 68': 'Barbiellini:1968rlx',
        'A2 01': 'deLeon:2001dnx', 'A2 22': 'A2:2021wmm', 'HI$\\gamma$S 22': 'Li:2022vnz', 'Illinois 91': 'Federspiel:1991yd',
        'Moscow 74': 'Baranov:1974ec,Baranov:1975ju', 'MIT 57': 'Pugh:1957zz', 'MIT 59': 'Hyman:1959zz',
        'Moscow 60': 'Goldansky:1960zz', 'Chicago 58': 'Oxley:1958zz', 'Illinois 60': 'Bernardini:1960wya',
        'A2 15': 'A2:2014iky', 'Martel thesis (unpubl.)': 'Martel:2013phd', 'A2 20': 'A2:2019bqm'}
OBSTEX = {'CS': r'$d\sigma/d\Omega_{\rm lab}$', 'CSCM': r'$d\sigma/d\Omega_{\rm cm}$', 'S3': r'$\Si_3$',
          'S2X': r'$\Si_{2x}$', 'S2Z': r'$\Si_{2z}$'}


def data_table(df):
    """LaTeX table of the data sets shown in the Delta-region figures (label = legend entry)."""
    shown = set()
    for d in SHOWN.values():
        shown |= set(d)
    q = df[df.pid.isin(SHOWN_PIDS)]
    T = [r'\begin{tabular}{llllrr}', r'\hline',
         r'Label & Ref. & Observable & $\theta$ frame & $E_\gamma$ [MeV] & $N$\\', r'\hline']
    for lab in [l for l in D.ORDER if l in set(q.label)]:
        z = q[q.label == lab]
        obs = ', '.join(OBSTEX[o] for o in dict.fromkeys(z.obs))
        fr = '/'.join(sorted(set(z.frame.str.lower().str.replace('cm', 'c.m.'))))
        er = ('%.0f' % z.E.min()) if z.E.max() - z.E.min() < 0.5 else '%.0f--%.0f' % (z.E.min(), z.E.max())
        T.append(r'%s & \cite{%s} & %s & %s & %s & %d\\' % (lab, CITE[lab], obs, fr, er, len(z)))
    T += [r'\hline', r'\end{tabular}']
    open(os.path.join(C.RESULTS, 'data_table.tex'), 'w').write('\n'.join(T) + '\n')


# ============================================================== main
def main():
    df = D.build()
    df.to_csv(os.path.join(C.RESULTS, 'data_display.csv'), index=False)
    pts = build_points(df)
    # unique ids
    seen = set(); P = []
    for p in pts:
        if p[0] not in seen:
            seen.add(p[0]); P.append(p)
    print('executable :', C.EXE)
    print('&v2 options:', C.v2opts(True))
    print('main curve :', 'NNLO prediction (P = Pref, nothing fitted)' if C.MAIN == 'pred' else 'NNLO fit')
    print('NNLO fit   : alpha, beta =', C.ALPHA_FIT, C.BETA_FIT, ' cov', C.COV_AB.ravel(), ' from', C.FITS_JSON)
    print('band       :', 'LMP15-type envelope (static alpha, beta kept, no fit error)' if C.MAIN == 'pred' else
          'envelope without static alpha, beta (+) fit error in quadrature')
    print('spin-fit set (validation only, not drawn):', C.P_SPINFIT)
    print('extra curve:', {'run3': 'running Delta parameters, irun=3 (values and slopes fixed at s0; static pols unchanged)',
                           'run2': 'running Delta parameters, irun=2, at the static pols of the main curve',
                           'run1': 'running Delta parameters, irun=1, at the static pols of the main curve',
                           'novc': 'without the gammaNDelta vertex corrections', 'none': 'none'}[C.EXTRA],
          '; also in the angular distributions' if C.EXTRA_THETA else '')
    print('points     :', len(P))
    T = TD.Theory(P, tag='dp')
    val = T.validate()
    json.dump(val, open(os.path.join(C.RESULTS, 'validation.json'), 'w'), indent=1)
    for k, v in val.items():
        print('  validation %-28s %s' % (k, v))
    cv = T.curves()
    json.dump(dict(main=C.MAIN, Pref=list(T.Pref), P_main=list(T.P_main()), P_fit=list(T.P_fit()), fits_json=C.FITS_JSON,
                   fit_key=C.FIT_KEY, static=T.static.tolist(), exe=C.EXE, v2opts=C.v2opts(True),
                   P_spinfit_validation_only=list(C.P_SPINFIT), cov_ab_fit=C.COV_AB.tolist(),
                   remove_static_ab=C.REMOVE_STATIC_AB, add_ab_error=C.ADD_AB_ERROR,
                   Lambda_chi=C.LAMBDA_CHI, hasVC=T.hasVC, Pref_band={k: list(v) for k, v in T.coupling_runs().items() if str(k).startswith('Pref')}, band_grid=[list(x) for x in C.BAND_GRID],
                   extra=C.EXTRA, extra_irun=C.EXTRA_IRUN, extra_in_theta=C.EXTRA_THETA,
                   extra_v2opts=(C.v2opts(True, C.EXTRA_IRUN) if C.EXTRA_IRUN is not None else None),
                   Pref_extra=(list(T.full_run.Pref) if T.full_run is not None else None)),
              open(os.path.join(C.RESULTS, 'run_info.json'), 'w'), indent=1)
    ix = {p[0]: i for i, p in enumerate(P)}
    nflag = int(np.sum(cv['flag'])); print('  extrapolated (flag=1) points:', nflag)
    extra_cols = tuple(dict.fromkeys(k for k in ('noVC', EXTRA_KEY) if k and k in cv))     # 'noVC' always, plus 'run3'
    # data with theory
    dd = df[df.obs.isin(['CS', 'CSCM', 'S3', 'S2X', 'S2Z'])].copy()
    ii = np.array([ix[p] for p in dd.pid])
    for k in ('main', 'lo', 'hi', 'cpl_lo', 'cpl_hi', 'fit', 'tree') + extra_cols:
        dd['th_' + k] = cv[k][ii]
    dd.to_csv(os.path.join(C.RESULTS, 'data_with_theory.csv'), index=False)
    # curves for reproducibility
    cur = pd.DataFrame({'id': [p[0] for p in P], 'E': [p[1] for p in P], 'theta': [p[2] for p in P],
                        'obs': [p[3] for p in P], 'frame': [p[4] for p in P]})
    for k in ('main', 'lo', 'hi', 'trunc_lo', 'trunc_hi', 'ab', 'cpl_lo', 'cpl_hi', 'fit', 'tree', 'flag') + extra_cols:
        cur[k] = cv[k]
    cur[cur.id.str.startswith('g_')].to_csv(os.path.join(C.RESULTS, 'curves_grid.csv'), index=False, float_format='%.6g')
    ext = D.build_ext450()
    print('  Fig. 1: plotting-only points 400 < E <= 450 MeV:', ext.groupby('dataset').E.agg(['size', 'min', 'max']).to_dict('index'))
    lara = D.build_lara()
    print('  Fig. 1: plotting-only LARA (Bonn 01) points, digitized:', lara.groupby('theta_cm').E.agg(['size', 'min', 'max']).to_dict('index'))
    fig1(pd.concat([df, ext, lara], ignore_index=True), cv, ix)
    lara_ang = D.build_lara_ang()
    print('  Fig. 2: plotting-only LARA angular distributions (Bonn 01), digitized:', lara_ang.groupby('E').size().to_dict())
    df2 = pd.concat([df, ext, lara_ang], ignore_index=True)       # + the plotting-only points above 400 MeV
    c2 = fig2(df2, cv, ix, FIG2, 'fig_delta_cs_vstheta')
    print('  Fig. 2: lab-frame sets converted to c.m.:', sorted(c2))
    fig3a(df, cv, ix)
    fig3b(df, cv, ix)
    fig4(df, cv, ix)
    fig5(df, cv, ix)
    fig6(df, cv, ix)
    tabs = chi2_table(df, cv, ix)
    pd.set_option('display.width', 250)
    print(tabs['350'][['set', 'N', 'Erange', 'chi2N_n1', 'chi2N_nfit', 'n', 'chi2N_th', 'chi2N_fit_n1', 'chi2N_fit_nfit', 'n_fit']].round(3).to_string(index=False))
    provenance(df)
    data_table(df)
    # bin-averaged vs nominal asymmetries (Fig. 4), for the notes
    rep = {}
    for (o, E0), avg in cv.get('binavg', {}).items():
        idx = np.array([ix[gid('f4', o, E0, th)] for th in TH_GRID])
        rep['%s_%g' % (o, E0)] = float(np.max(np.abs(avg - cv['main'][idx])))
    json.dump(rep, open(os.path.join(C.RESULTS, 'binavg_maxdiff.json'), 'w'), indent=1)
    print('  bin-averaged vs nominal-energy asymmetry, max |diff|:', {k: round(v, 3) for k, v in rep.items()})


if __name__ == '__main__':
    main()
