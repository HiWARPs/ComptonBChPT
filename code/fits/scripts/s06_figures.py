"""Step 6: figures (PDF, shared style) in work/fits/figs/."""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from common import *  # noqa
import bchpt, database as DB

plt.style.use(STYLE)
os.makedirs(FIGS, exist_ok=True)
df = load_db()
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
FD = json.load(open(os.path.join(RES, 'delta_region.json')))
P_F1 = np.array(FL['F1_d_all_thB']['P']); P_F2 = np.array(FL['F2_d_all_thB']['P'])
P_MAIN = P_F2
P_F1n = np.array(FL['F1_d_all']['P'])
MP = 938.272

MARK = {'OlmosdeLeon2001': ('D', '#c0392b'), 'Federspiel1991': ('^', '#2e8b57'), 'Zieger1992': ('*', '#6a3d9a'),
        'Baranov1974_75_th90': ('v', '#8c564b'), 'Baranov1974_75_th150': ('v', '#8c564b'),
        'Hallin1993_TabII': ('s', '#444444'), 'Hallin1993_TabII_D': ('s', '#444444'),
        'MacGibbon1995_tagged': ('o', '#e08e0b'), 'MacGibbon1995_untagged': ('o', '#e08e0b'),
        'Pugh1957': ('p', '#17becf'), 'Hyman1959': ('h', '#9a9a2e'), 'Goldansky1960': ('<', '#2e8b57'),
        'Oxley1958': ('>', '#6a3d9a'), 'Bernardini1960': ('X', '#8c564b'),
        'A2_2022_CS': ('o', '#1f4e9c'), 'A2_2022_S3': ('o', '#1f4e9c'), 'HIGS_2022_CIRC': ('s', '#d62e8a'),
        'HIGS_2022_LIN': ('s', '#d62e8a'), 'LEGS2001_CS': ('d', '#e08e0b'), 'LEGS2001_S3': ('d', '#e08e0b'),
        'Hunger1997_75': ('o', '#1f4e9c'), 'Hunger1997_90': ('o', '#1f4e9c'), 'Wissmann1999': ('^', '#1f4e9c'),
        'Camen2002': ('v', '#1f4e9c'), 'Genzel1972': ('h', '#2e8b57'), 'Gray1967': ('<', '#6a3d9a'),
        'Baranov1966': ('>', '#8c564b'), 'Kolb2000': ('*', '#c0392b'), 'Nagashima1964': ('p', '#17becf'),
        'A2_Martel2015_S2x': ('o', '#1f4e9c'), 'A2_Paudyal2020_S2z': ('o', '#1f4e9c')}
C_NNLO, C_F1, C_F2, C_BAND = '#444444', '#1f4e9c', '#c0392b', '#9fb6dd'
from matplotlib.lines import Line2D   # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
W1, W2 = 246 / 72.27, 510 / 72.27     # final printed widths [in]: single / double column (revtex4-2 prd)
# short experiment labels for the legends (as in the Delta-region figures, work/delta_plots/data_dp.LABEL)
SHORT = {'OlmosdeLeon2001': 'TAPS 01', 'Federspiel1991': 'Illinois 91', 'Zieger1992': 'MAMI 92',
         'Baranov1974_75_th90': 'Moscow 74', 'Baranov1974_75_th150': 'Moscow 74', 'Hallin1993_TabII': 'SAL 93',
         'Hallin1993_TabII_D': 'SAL 93', 'MacGibbon1995_tagged': 'SAL 95', 'MacGibbon1995_untagged': 'SAL 95',
         'Pugh1957': 'MIT 57', 'Hyman1959': 'MIT 59', 'Goldansky1960': 'Moscow 60', 'Oxley1958': 'Chicago 58',
         'Bernardini1960': 'Illinois 60', 'A2_2022_CS': 'A2 22', 'A2_2022_S3': 'A2 22', 'HIGS_2022_CIRC': 'HI$\\gamma$S 22',
         'HIGS_2022_LIN': 'HI$\\gamma$S 22', 'LEGS2001_CS': 'LEGS 01', 'LEGS2001_S3': 'LEGS 01', 'Hunger1997_75': 'MAMI 97',
         'Hunger1997_90': 'MAMI 97', 'Wissmann1999': 'MAMI 99', 'Camen2002': 'MAMI 02', 'Kolb2000': 'SAL 00',
         'Genzel1972': 'Bonn 72', 'Gray1967': 'Illinois 67', 'Baranov1966': 'Moscow 66', 'Nagashima1964': 'Tokyo 64'}
# theory legend entries shared by the angular-distribution figures (Figs. cs, sigma3 and higs of the paper); the
# forward and Emax figures use the first two / the first entry only (no PWA curve there)
TH_LAB = ['NNLO B$\\chi$PT', 'NNLO fit', 'PWA fit']      # band = LEC uncertainty, stated in the captions (user, 2026-10-17)
DLEC = 0.7   # natural size of the O(p^4) LECs (LP 2010, Eq. 25): band of the prediction in Figs. cs and sigma3
# PWA fit of this work (work/pwa/scenarios/main_fits_legacynorm.py, fit MAIN): bare curves without normalisation
# factors, evaluated at the panel energies with the residual functions of the nearest PWA bin by
# work/pwa/scenarios/pwa_curves.py (columns obs, E, theta_lab, value; obs = CS, S3, SPAR, SPERP as in the fit database)
C_PWA = '#6b8e23'   # olive green
PWA = pd.read_csv(os.path.join(FITS, '..', 'pwa', 'results', 'pwa_curves.csv'))


def pwa_curve(obs, E, thmin=0.0):
    """(theta_lab, value) of the PWA-fit curve for observable obs at lab energy E, or None if the CSV has none"""
    s = PWA[(PWA.obs == obs) & np.isclose(PWA.E, E) & (PWA.theta_lab >= thmin - 1e-9)].sort_values('theta_lab')
    return (s.theta_lab.values, s.value.values) if len(s) else None
DORDER = ['A2 22', 'HI$\\gamma$S 22', 'TAPS 01', 'SAL 95', 'SAL 93', 'Illinois 91', 'MAMI 92', 'Moscow 74', 'MIT 57', 'MIT 59',
          'Moscow 60', 'Chicago 58', 'Illinois 60']   # legend order of the data sets


def tex(s):
    """LaTeX-safe (usetex) version of a database label"""
    for a, b in (('\u2225,\u22a5', '$_{\\parallel,\\perp}$'), ('\u03b3', '$\\gamma$'), ('\u03c3', '$\\sigma$'),
                 ('\u03a33', '$\\Sigma_3$'), ('\u03a32x', '$\\Sigma_{2x}$'), ('\u03a32z', '$\\Sigma_{2z}$'),
                 ('\u00b0', '$^\\circ$'), ('\u00e9', "\\'e"), ('\u00fc', '\\"u'), ('%', '\\%')):
        s = s.replace(a, b)
    return s


def th_handles():
    return [(Patch(facecolor=C_BAND, alpha=0.6, lw=0), Line2D([], [], color=C_NNLO, lw=1.2)),
            Line2D([], [], color=C_F1, ls='--', lw=1.0),
            Line2D([], [], color=C_PWA, ls=':', lw=1.1)]


def save_exact(fig, name):
    """save at the exact figure size (no cropping): the PDF is included in the paper without scaling"""
    with plt.rc_context({'savefig.bbox': 'standard'}):
        fig.savefig(os.path.join(FIGS, name + '.pdf'))
    plt.close(fig)


def grid_rows(spec):
    """spec: list of (tag, E, theta, dobs, frame) -> DataFrame usable by outputs()/truncation_for()"""
    rows = []
    for tag, E, th, o, fr in spec:
        thl = th if fr == 'LAB' else float(DB.th_cm2lab(E, th))
        thc = th if fr == 'CM' else float(DB.th_lab2cm(E, th))
        rows.append(dict(pid='%s_%s_%.3f_%.3f_%s' % (tag, o, E, th, fr), dataset='grid', E=E, theta=th, theta_lab=thl,
                         theta_cm=thc, dobs=o, frame=fr, value=0.0, err_pt=1.0, norm=0.0, flags=''))
    return pd.DataFrame(rows)


def curves(g, tag, P_list, band=True, extra_cfgs=()):
    T = outputs(g, 'NNLO', tag=tag)
    out = {'NNLO': T.value_file.copy()}
    for k, P in P_list.items():
        out[k] = T.predict_array(P)
    if band:
        _, tr = truncation_for(g, g, tag=tag)
        out['sig'] = tr['sigma']
    ea, eb = np.eye(6)[0], np.eye(6)[1]
    da = (T.predict_array(PREF_NNLO + DLEC * ea) - T.predict_array(PREF_NNLO - DLEC * ea)) / 2
    db = (T.predict_array(PREF_NNLO + DLEC * eb) - T.predict_array(PREF_NNLO - DLEC * eb)) / 2
    out['lec'] = np.sqrt(da ** 2 + db ** 2)      # uncertainty of the prediction from LECs of natural size
    for c in extra_cfgs:
        Tc = outputs(g, c, tag=tag)
        out[c] = Tc.predict_array(np.r_[P_list['F1'][:2], Tc.Pref[2:]])
    return out


def to_lab(r):
    """data row -> (theta_lab, value_lab, err_lab) for unpolarised cross sections"""
    if r['dobs'] == 'CSCM':
        b = r['E'] / (r['E'] + MP)
        c = np.cos(np.radians(r['theta_lab']))
        J = (1 - b * b) / (1 - b * c) ** 2
        return r['theta_lab'], r['value'] * J, r['err_pt'] * J
    return r['theta_lab'], r['value'], r['err_pt']


def to_cm(r):
    if r['dobs'] == 'CS':
        b = r['E'] / (r['E'] + MP)
        c = np.cos(np.radians(r['theta_lab']))
        J = (1 - b * b) / (1 - b * c) ** 2
        return r['theta_cm'], r['value'] / J, r['err_pt'] / J
    return r['theta_cm'], r['value'], r['err_pt']


def vsE(conv):
    def f(r):
        _, v, e = conv(r)
        return r['E'], v, e
    return f


def plot_data(ax, rows, conv=to_lab, seen=None, ms=3.5):
    seen = set() if seen is None else seen
    for ds, d in rows.groupby('dataset'):
        m, c = MARK.get(ds, ('o', 'k'))
        xs, ys, es = zip(*[conv(r) for _, r in d.iterrows()])
        ex = SHORT.get(ds, tex(d['expt'].iloc[0]))
        lab = ex if ex not in seen else None
        seen.add(ex)
        ax.errorbar(xs, ys, es, fmt=m, color=c, mfc='white' if ds.startswith('Hallin') or ds in ('Pugh1957', 'Hyman1959') else c,
                    ms=ms, lw=0.7, label=lab, zorder=5)
    return seen


# ---------------------------------------------------------------- (i) angular distributions dsigma/dOmega_lab
def forward_point(E):
    """empirical forward cross section of GHP (rows 'GHP_fwd_CS' of the database) interpolated to energy E:
    cubic spline through the 12 tabulated values, uncertainty linearly interpolated (the quoted errors)."""
    from scipy.interpolate import CubicSpline
    d = df[df.dataset == 'GHP_fwd_CS'].sort_values('E')
    return float(CubicSpline(d.E.values, d.value.values)(E)), float(np.interp(E, d.E.values, d.err_pt.values))


def fig_cs(energies, fname, window, title_fmt='%.1f MeV', ylim=None, omit=(), fwd=False, own_range=()):
    """angular distributions; ylim: common y-range of all panels (extended minimally if anything would be clipped);
    omit: data sets not drawn (they stay in the fit); fwd: add the GHP forward point at theta = 0 (open star);
    own_range: energies whose panel gets its own y-range instead of the common one."""
    th = np.arange(0.0 if fwd else 20.0, 175.1, 2.5)
    spec = [('cs', E, t, 'CS', 'LAB') for E in energies for t in th]
    g = grid_rows(spec)
    cv = curves(g, 'fig_cs_' + fname, {'F1': P_F1, 'F2': P_F2})
    n = len(energies); nc = 4 if n > 6 else (3 if n > 4 else 2)
    nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(W2, 2.0 * nr + 0.55), sharex=True, layout='constrained')
    axs = np.atleast_1d(axs).ravel()
    seen = set()
    ymin, ymax = np.inf, -np.inf          # extremes of everything drawn in the common-range panels (clipping check)
    pext = []                              # per-panel extremes
    for i, E in enumerate(energies):
        ax = axs[i]
        m = np.isclose(g.E.values, E)
        ymin0, ymax0 = ymin, ymax
        ax.fill_between(th, cv['NNLO'][m] - cv['lec'][m], cv['NNLO'][m] + cv['lec'][m], color=C_BAND, alpha=0.6, lw=0)
        ax.plot(th, cv['NNLO'][m], '-', color=C_NNLO, lw=1.2)
        ax.plot(th, cv['F2'][m], '--', color=C_F1, lw=1.0)
        pw = pwa_curve('CS', E, thmin=th[0])
        if pw is not None:
            ax.plot(*pw, ':', color=C_PWA, lw=1.1)
        else:
            print('%s: no PWA curve at %.2f MeV (not in pwa_curves.csv)' % (fname, E))
        sel = df[(df.dobs.isin(['CS', 'CSCM'])) & (np.abs(df.E - E) <= window) & ~df.dataset.isin(FWD) & (df.theta_lab > 0)
                 & ~df.dataset.isin(omit)]
        seen = plot_data(ax, sel, seen=seen)
        ymin = min(ymin, cv['F2'][m].min(), cv['NNLO'][m].min())      # curves and data only (the band edge may be clipped)
        ymax = max(ymax, cv['F2'][m].max(), cv['NNLO'][m].max())
        if pw is not None:
            ymin = min(ymin, pw[1].min()); ymax = max(ymax, pw[1].max())
        if len(sel):
            vals = np.array([to_lab(r)[1:] for _, r in sel.iterrows()])
            ymin = min(ymin, (vals[:, 0] - vals[:, 1]).min()); ymax = max(ymax, (vals[:, 0] + vals[:, 1]).max())
        if fwd:
            v, e = forward_point(E)
            ax.errorbar([0.0], [v], [e], fmt='*', color='k', mfc='white', ms=6.5, lw=0.7, zorder=6,
                        label='forward point' if 'forward point' not in seen else None)
            seen.add('forward point')
            ymin = min(ymin, v - e); ymax = max(ymax, v + e)
        pext.append((ymin, ymax))
        if E in own_range:                 # this panel does not constrain the common range
            ymin, ymax = ymin0, ymax0
        ax.text(0.04, 0.92, title_fmt % E, transform=ax.transAxes, va='top')
        ax.set_xlim(-6 if fwd else 15, 180)
        if ylim is None:
            lo = min(cv['F2'][m].min(), sel.value.min() if len(sel) else 99) * 0.75
            hi = max(cv['F2'][m].max(), (sel.value + sel.err_pt).max() if len(sel) else 0) * 1.12
            ax.set_ylim(max(0, lo), hi)
        if i % nc == 0:
            ax.set_ylabel('$d\\sigma/d\\Omega_{\\rm lab}$ [nb/sr]')
        if i >= n - nc:
            ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')
    if ylim is not None:                 # common range; extend minimally (to 0.5 nb/sr) if anything would be clipped
        lo, hi = ylim
        if ymin < lo:
            lo = np.floor(ymin * 2) / 2
        if ymax > hi:
            hi = np.ceil(ymax * 2) / 2
        print('%s: drawn range %.2f-%.2f nb/sr; y-axis %s -> (%.1f, %.1f)%s' % (
            fname, ymin, ymax, ylim, lo, hi, '' if (lo, hi) == tuple(ylim) else '  [EXTENDED to avoid clipping]'))
        for i, ax in enumerate(axs[:n]):
            if energies[i] in own_range:
                pmin, pmax = pext[i]
                ax.set_ylim(np.floor(pmin * 2) / 2, np.ceil(pmax * 2) / 2)
                print('%s: panel %s MeV with its own range (%.1f, %.1f)' % (fname, energies[i], np.floor(pmin * 2) / 2, np.ceil(pmax * 2) / 2))
            else:
                ax.set_ylim(lo, hi)
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
    for j in range(n, len(axs)):
        axs[j].axis('off')
    h, l = [], []
    for ax in axs:
        hh, ll = ax.get_legend_handles_labels(); h += hh; l += ll
    o = [k for k in DORDER if k in l] + [k for k in l if k not in DORDER]
    h, l = [h[l.index(k)] for k in o], o
    fig.legend(th_handles(), TH_LAB, loc='outside upper center', ncol=3, columnspacing=2.0)
    if len(axs) > n:          # the data legend goes into the empty panel slot
        axs[n].legend(h, l, loc='center', ncol=1, handletextpad=0.3, frameon=False)
    else:
        fig.legend(h, l, loc='outside lower center', ncol=min(len(l), 8), handletextpad=0.3, columnspacing=1.6)
    save_exact(fig, fname)


# paper Fig. "cs": common y-range, Pugh 1957 (MIT 57) omitted for clarity, GHP forward point at theta = 0
fig_cs([81.3, 92.2, 103.2, 113.6, 124.1, 135.3, 149.87], 'fig1a_cs_A2_HIGS', window=2.5, ylim=(4.5, 24.5), omit=('Pugh1957',), fwd=True,
       own_range=(149.87,))
fig_cs([58.94, 79.19, 108.42, 149.87], 'fig1b_cs_older', window=5.0)


# ---------------------------------------------------------------- (ii) Sigma_3 (A2 + HIGS-derived)
def higs_sigma3():
    s = df[df.dataset == 'HIGS_2022_LIN']
    out = []
    for th, d in s.groupby('theta_lab'):
        par = d[d.dobs == 'SPAR']; per = d[d.dobs == 'SPERP']
        wp = 1 / par.err_pt ** 2
        P, dP = np.sum(wp * par.value) / wp.sum(), 1 / np.sqrt(wp.sum())
        Q, dQ = per.value.iloc[0], per.err_pt.iloc[0]
        S = (P - Q) / (P + Q)
        dS = 2 / (P + Q) ** 2 * np.sqrt((Q * dP) ** 2 + (P * dQ) ** 2)
        out.append((th, S, dS))
    return np.array(out)


S3H = higs_sigma3()
th = np.arange(20, 175.1, 2.5)
Es = [83.4, 92.2, 108.5, 129.5]
g = grid_rows([('s3', E, t, 'S3', 'LAB') for E in Es for t in th])
cv = curves(g, 'fig_s3', {'F1': P_F1, 'F2': P_F2, 'F1n': P_F1n})
fig, axs = plt.subplots(2, 2, figsize=(W1, 3.9), sharex=True, sharey=True, layout='constrained')   # 3.9 in: 3-row legend (was 3.75 with 2 rows)
axs = axs.ravel()
for i, E in enumerate(Es):
    ax = axs[i]; m = np.isclose(g.E.values, E)
    ax.fill_between(th, cv['NNLO'][m] - cv['lec'][m], cv['NNLO'][m] + cv['lec'][m], color=C_BAND, alpha=0.6, lw=0)
    ax.plot(th, cv['NNLO'][m], '-', color=C_NNLO, lw=1.2, label=TH_LAB[0])
    ax.plot(th, cv['F2'][m], '--', color=C_F1, lw=1.0, label='NNLO fit')
    pw = pwa_curve('S3', E, thmin=th[0])
    if pw is not None:
        ax.plot(*pw, ':', color=C_PWA, lw=1.1, label='PWA fit')
    else:
        print('fig2_sigma3: no PWA curve at %.2f MeV (not in pwa_curves.csv)' % E)
    if E == 83.4:
        hH = ax.errorbar(S3H[:, 0], S3H[:, 1], S3H[:, 2], fmt='s', color='#d62e8a', ms=3.5, label='HI$\\gamma$S 22 (derived)')
    else:
        d = df[(df.dataset == 'A2_2022_S3') & np.isclose(df.E, E)]
        hA = ax.errorbar(d.theta_lab, d.value, d.err_pt, fmt='o', color='#1f4e9c', ms=3.5, label='A2 22')
    ax.text(0.97, 0.05, '%.1f MeV' % E, transform=ax.transAxes, ha='right', va='bottom')   # bottom right: free of curves and data at all four energies
    ax.set_xlim(15, 180); ax.set_ylim(-1.08, 0.08)                                           # slightly below -1 to slightly above 0 (user, 2026-10-17)
    if i >= 2:
        ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')
    if i % 2 == 0:
        ax.set_ylabel('$\\Sigma_3$')
# left column: theory (3 entries), right column: data (2); columnspacing 1.0 and handlelength 1.5 keep the legend within W1
# (3.35 in; with the defaults 1.5 / 1.8 its box measures 3.47 in > W1 = 3.40 in)
leg = fig.legend(th_handles() + [hA, hH], TH_LAB + ['A2 22', 'HI$\\gamma$S 22 (derived)'], loc='outside upper center', ncol=2,
                 columnspacing=1.0, handlelength=1.5)
fig.canvas.draw()
_wl = leg.get_window_extent(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).width
print('fig2_sigma3: legend width %.2f in within the figure width %.2f in%s' % (_wl, W1, '' if _wl <= W1 else '  [OVERFLOW]'))
save_exact(fig, 'fig2_sigma3')

# ---------------------------------------------------------------- (iii) HIGS linearly polarised cross sections
# paper Fig. "higs": left dsigma_par, dsigma_perp at 83.4 MeV; right unpolarised (circular beam) at 81.3 MeV.  Curves as
# in Fig. cs: NNLO prediction (solid, with the natural-LEC band), NNLO fit (dashed), plus the fit to the HIGS data alone
# (dash-dotted, with its fitted normalisations).
# The GHP forward point (open star, as in Fig. cs) is interpolated to the panel energy; at theta = 0 the polarised
# cross sections both equal the unpolarised one, so the single star applies to both curves of the left panel.
th = np.arange(0.0, 175.1, 2.5)
g = grid_rows([('hl', 83.4, t, o, 'LAB') for o in ['SPAR', 'SPERP'] for t in th] + [('hc', 81.3, t, 'CS', 'LAB') for t in th])
cv = curves(g, 'fig_higs', {'F2': P_F2, 'H': np.r_[np.array(FL['F1_HIGS_free']['P'])[:2], PREF_NNLO[2:]]}, band=False)
nH = {k: FL['F1_HIGS_free']['perset'][k]['n'] for k in HIGS}
cv['H'] = np.where((g.dobs == 'CS').values, nH['HIGS_2022_CIRC'], nH['HIGS_2022_LIN']) * cv['H']
C_HFIT = '#2e8b57'
# distinct colour and symbol for each of the three HIGS measurements (circular beam as in Fig. cs)
HSTY = {'SPAR': ('o', '#e08e0b', 'HI$\\gamma$S 22 $d\\sigma_\\parallel$'), 'SPERP': ('^', '#6a3d9a', 'HI$\\gamma$S 22 $d\\sigma_\\perp$'),
        'CS': ('s', '#d62e8a', 'HI$\\gamma$S 22 (circ.)')}
fig, axs = plt.subplots(1, 2, figsize=(W2, 2.75), layout='constrained')
hd, ld = [], []                      # data legend (lower row): the three HIGS measurements and the forward point
for ax, E, obs in [(axs[0], 83.4, ['SPAR', 'SPERP']), (axs[1], 81.3, ['CS'])]:
    for o in obs:
        m = (g.dobs == o).values
        ax.fill_between(th, cv['NNLO'][m] - cv['lec'][m], cv['NNLO'][m] + cv['lec'][m], color=C_BAND, alpha=0.6, lw=0)
        ax.plot(th, cv['NNLO'][m], '-', color=C_NNLO, lw=1.2)
        ax.plot(th, cv['F2'][m], '--', color=C_F1, lw=1.0)
        ax.plot(th, cv['H'][m], '-.', color=C_HFIT, lw=1.0)
        pw = pwa_curve(o, E)
        if pw is not None:
            ax.plot(*pw, ':', color=C_PWA, lw=1.1)
        else:
            print('fig3_higs_linear: no PWA curve for %s at %.2f MeV (not in pwa_curves.csv)' % (o, E))
        print('fig3_higs_linear %-5s %.1f MeV: LEC band max %.2f nb/sr, max %.1f%% of the prediction (at %.0f deg)' % (
            o, E, cv['lec'][m].max(), 100 * (cv['lec'][m] / cv['NNLO'][m]).max(), th[np.argmax(cv['lec'][m] / cv['NNLO'][m])]))
        mk, c, lab = HSTY[o]
        d = df[(df.dataset == ('HIGS_2022_CIRC' if o == 'CS' else 'HIGS_2022_LIN')) & (df.dobs == o)]
        hd.append(ax.errorbar(d.theta_lab, d.value, d.err_pt, fmt=mk, color=c, ms=3.5, lw=0.7, zorder=5)); ld.append(lab)
    v, e = forward_point(E)
    hf = ax.errorbar([0.0], [v], [e], fmt='*', color='k', mfc='white', ms=6.5, lw=0.7, zorder=6)
    ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]'); ax.set_xlim(-6, 180)
hd.append(hf); ld.append('forward point')
axs[0].set_ylabel('$d\\sigma/d\\Omega_{\\rm lab}$ [nb/sr]')
axs[0].text(0.5, 0.95, '$E_\\gamma=83.4$ MeV, linear polarisation', transform=axs[0].transAxes, ha='center', va='top')
axs[0].set_ylim(0, 26)
axs[1].text(0.04, 0.95, '$E_\\gamma=81.3$ MeV, unpolarised', transform=axs[1].transAxes, va='top')
axs[1].set_ylim(4, 21); axs[1].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
fig.legend(th_handles() + [Line2D([], [], color=C_HFIT, ls='-.', lw=1.0)], TH_LAB + ['fit to HI$\\gamma$S only'],
           loc='outside upper center', ncol=4, columnspacing=2.0)
fig.legend(hd, ld, loc='outside lower center', ncol=4, handletextpad=0.3, columnspacing=1.6)
save_exact(fig, 'fig3_higs_linear')


# ---------------------------------------------------------------- (iv) alpha-beta plane
def ellipse(ax, P, cov, **kw):
    c = np.array(cov)[:2, :2]
    w, v = np.linalg.eigh(c)
    k = np.sqrt(2.30)     # 68% joint CL in 2D
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    e = Ellipse((P[0], P[1]), 2 * k * np.sqrt(w[1]), 2 * k * np.sqrt(w[0]), angle=ang, fill=kw.pop('fill', False), **kw)
    ax.add_patch(e)
    return e


fig, ax = plt.subplots(figsize=(W1, 3.85), layout='constrained')
a = np.linspace(8, 16, 50)
h_bal = ax.fill_between(a, 14.0 - 0.2 - a, 14.0 + 0.2 - a, color='#dddddd', lw=0, label='Baldin SR, average 14.0(2)')
# the two photoabsorption fits of Gryniuk et al. 2015 (fit I: MAID below the 2pi threshold; fit II: SAID), hatched
with plt.rc_context({'hatch.linewidth': 0.5}):
    h_b1 = ax.fill_between(a, 14.29 - 0.27 - a, 14.29 + 0.27 - a, facecolor='none', edgecolor='#8c8c8c', hatch='////', lw=0,
                           label='Baldin SR, fit I 14.29(27)')
    h_b2 = ax.fill_between(a, 13.85 - 0.22 - a, 13.85 + 0.22 - a, facecolor='none', edgecolor='#8c8c8c', hatch='\\\\\\\\', lw=0,
                           label='Baldin SR, fit II 13.85(22)')
# the NNLO fit of this work (68% joint CL)
r = FL['F2_d_all_thB']
h_fit = ellipse(ax, r['P'], r['covP'], ec=C_F1, ls='-', lw=2.0, label='B$\\chi$PT fit, this work')
ax.plot(*r['P'][:2], 'o', color=C_F1, ms=2)
# NNLO prediction: the amplitude used here (physical masses, second row of Table NNLOpols) with the LMP15 uncertainty
h_nnlo = ax.errorbar(PREF_NNLO[0], PREF_NNLO[1], 0.7, 0.7, fmt='P', color='k', ms=5, lw=0.7, label='NNLO B$\\chi$PT')
# determinations that involve the new data (A2 2022, HIGS 2022), and the PWA of this work
others = [((11.59, 0.80), (2.58, 0.77), '#c0392b', 'X', 'PWA, this work'),
          ((10.99, np.hypot(0.16, np.hypot(0.47, np.hypot(0.17, 0.34)))), (3.14, np.hypot(0.21, np.hypot(0.24, np.hypot(0.20, 0.35)))), '#1f4e9c', 'v', 'A2 22'),
          ((13.8, 1.25), (0.2, 1.25), '#d62e8a', '^', 'HI$\\gamma$S 22'),
          ((12.7, 0.8), (2.4, 0.6), '#e08e0b', '>', 'MRPP 22 (DR)'),
          ((11.5, 0.4), (2.31, 0.29), '#6a3d9a', 's', 'PDG 2026')]
hs = []
for (av, ae), (bv, be), c, m, lab in others:
    hs.append(ax.errorbar(av, bv, be, ae, fmt=m, color=c, mfc=c, ms=4.5, lw=0.7, label=lab, clip_on=True))
ax.set_xlim(9.5, 15.2); ax.set_ylim(-0.5, 5.0)
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1)); ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
ax.set_xlabel('$\\alpha_{E1}$ [$10^{-4}$ fm$^3$]'); ax.set_ylabel('$\\beta_{M1}$ [$10^{-4}$ fm$^3$]')
handles = [h_bal, h_b1, h_b2, h_nnlo, h_fit] + hs
fig.legend(handles, [h.get_label() for h in handles], loc='outside lower center', ncol=2, handlelength=2.0,
           handletextpad=0.5, columnspacing=1.0)
save_exact(fig, 'fig4_alpha_beta')

# ---------------------------------------------------------------- (v) Delta region: Sigma_2x, Sigma_2z vs angle
th = np.arange(20, 175.1, 2.5)
spec = [('s2x', 288.0, t, 'S2X', 'LAB') for t in th] + [('s2z', E, t, 'S2Z', 'LAB') for E in [275.0, 295.0] for t in th]
g = grid_rows(spec)
SP = {'S_A2_LEGS_gpi': ('-', C_F1, 'spin fit: A2 $\\Sigma_{2x,2z}$+LEGS $\\Sigma_3$+$\\gamma_\\pi$'),
      'S_A2_2free': ('-.', C_F2, 'spin fit: A2 $\\Sigma_{2x,2z}$ ($\\gamma_{E1E1},\\gamma_{M1M1}$)')}
T = outputs(g, 'NNLO', tag='fig_delta')
vnv = outputs(g, 'novtx', tag='fig_delta').predict_array(np.r_[P_MAIN[:2], PREF_NNLO[2:]]) if PROFILE == 'v2' else None
Tg1, Tg2 = outputs(g, 'gM2.85', tag='fig_delta'), outputs(g, 'gM3.10', tag='fig_delta')
vals = {k: T.predict_array(np.array(FD[k]['P'])) for k in SP}
vn = T.predict_array(np.r_[P_MAIN[:2], PREF_NNLO[2:]])
v1 = Tg1.predict_array(np.r_[P_MAIN[:2], Tg1.Pref[2:]]); v2 = Tg2.predict_array(np.r_[P_MAIN[:2], Tg2.Pref[2:]])
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.4))
for ax, (o, E, ds) in zip(axs, [('S2X', 288.0, 'A2_Martel2015_S2x'), ('S2Z', 275.0, 'A2_Paudyal2020_S2z'),
                                ('S2Z', 295.0, 'A2_Paudyal2020_S2z')]):
    m = ((g.dobs == o) & np.isclose(g.E, E)).values
    ax.fill_between(th, np.minimum(v1[m], v2[m]), np.maximum(v1[m], v2[m]), color='#dddddd', lw=0, label='$g_M=2.85$-$3.10$')
    ax.plot(th, vn[m], '--', color=C_NNLO, lw=1.0, label='NNLO ($\\alpha,\\beta$ of NNLO fit)')
    if vnv is not None:
        ax.plot(th, vnv[m], ':', color='#8c564b', lw=1.0, label='same, no $\\gamma N\\Delta$ vertex corr.')
    for k, (ls, c, lab) in SP.items():
        ax.plot(th, vals[k][m], ls, color=c, lw=1.1, label=lab)
    d = df[(df.dataset == ds) & np.isclose(df.E, E)]
    ax.errorbar(d.theta_lab, d.value, d.err_pt, fmt='o', color='k', ms=3.5, label='A2 data')
    ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]'); ax.set_xlim(15, 180)
    ax.set_ylabel('$\\Sigma_{2x}$' if o == 'S2X' else '$\\Sigma_{2z}$')
    ax.text(0.05, 0.92, '$E_\\gamma=%.0f$ MeV' % E, transform=ax.transAxes, va='top')
axs[0].set_ylim(-0.75, 0.35); axs[1].set_ylim(-0.6, 1.05); axs[2].set_ylim(-0.6, 1.05)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(rect=(0, 0.12, 1, 1)); fig.savefig(os.path.join(FIGS, 'fig5a_delta_S2x_S2z.pdf')); plt.close(fig)

# LEGS Sigma_3 vs angle at three energies
Es = [254.1, 286.5, 310.1]
th = np.arange(30, 165.1, 2.5)
g = grid_rows([('ls3', E, t, 'S3', 'CM') for E in Es for t in th])
T = outputs(g, 'NNLO', tag='fig_legs3')
vnv = outputs(g, 'novtx', tag='fig_legs3').predict_array(np.r_[P_MAIN[:2], PREF_NNLO[2:]]) if PROFILE == 'v2' else None
Tg1, Tg2 = outputs(g, 'gM2.85', tag='fig_legs3'), outputs(g, 'gM3.10', tag='fig_legs3')
vn = T.predict_array(np.r_[P_MAIN[:2], PREF_NNLO[2:]])
v1 = Tg1.predict_array(np.r_[P_MAIN[:2], Tg1.Pref[2:]]); v2 = Tg2.predict_array(np.r_[P_MAIN[:2], Tg2.Pref[2:]])
vs = T.predict_array(np.array(FD['S_A2_LEGS_gpi']['P']))
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.2), sharey=True)
for ax, E in zip(axs, Es):
    m = np.isclose(g.E.values, E)
    ax.fill_between(th, np.minimum(v1[m], v2[m]), np.maximum(v1[m], v2[m]), color='#dddddd', lw=0, label='$g_M=2.85$-$3.10$')
    ax.plot(th, vn[m], '--', color=C_NNLO, lw=1.0, label='NNLO ($\\alpha,\\beta$ of NNLO fit)')
    if vnv is not None:
        ax.plot(th, vnv[m], ':', color='#8c564b', lw=1.0, label='same, no vertex corr.')
    ax.plot(th, vs[m], '-', color=C_F1, label='spin fit (A2+LEGS+$\\gamma_\\pi$)')
    d = df[(df.dataset == 'LEGS2001_S3') & np.isclose(df.E, E)]
    ax.errorbar(d.theta_cm, d.value, d.err_pt, fmt='d', color='#e08e0b', ms=3.5, label='LEGS 01')
    ax.set_xlabel('$\\theta_{\\rm cm}$ [deg]'); ax.text(0.05, 0.92, '$E_\\gamma=%.1f$ MeV' % E, transform=ax.transAxes, va='top')
axs[0].set_ylabel('$\\Sigma_3$')
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=5, fontsize=6.5, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(os.path.join(FIGS, 'fig5b_delta_LEGS_S3.pdf')); plt.close(fig)

# dsigma vs energy at fixed angles up to 350 MeV
E = np.arange(140, 352, 2.0)
panels = [('CM', 75.0, 'CSCM', '$\\theta_{\\rm cm}=75^\\circ$'), ('CM', 90.0, 'CSCM', '$\\theta_{\\rm cm}=90^\\circ$'),
          ('CM', 115.0, 'CSCM', '$\\theta_{\\rm cm}=115^\\circ$'), ('LAB', 133.0, 'CS', '$\\theta_{\\rm lab}\\approx 131$-$136^\\circ$')]
g = grid_rows([('dcs%d' % i, e, a, o, fr) for i, (fr, a, o, _) in enumerate(panels) for e in E])
T = outputs(g, 'NNLO', tag='fig_dcs')
vnv = outputs(g, 'novtx', tag='fig_dcs').predict_array(np.r_[P_MAIN[:2], PREF_NNLO[2:]]) if PROFILE == 'v2' else None
Tg1, Tg2 = outputs(g, 'gM2.85', tag='fig_dcs'), outputs(g, 'gM3.10', tag='fig_dcs')
Tm = outputs(g, 'gM2.90', tag='fig_dcs')
vn = T.predict_array(np.r_[P_MAIN[:2], PREF_NNLO[2:]])
v0 = T.predict_array(PREF_NNLO)
v1 = Tg1.predict_array(np.r_[P_MAIN[:2], Tg1.Pref[2:]]); v2 = Tg2.predict_array(np.r_[P_MAIN[:2], Tg2.Pref[2:]])
vm = Tm.predict_array(np.r_[P_MAIN[:2], Tm.Pref[2:]])
fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6))
axs = axs.ravel(); seen = set()
for i, (fr, a, o, lab) in enumerate(panels):
    ax = axs[i]; m = g.pid.str.startswith('dcs%d_' % i).values
    ax.fill_between(E, np.minimum(v1[m], v2[m]), np.maximum(v1[m], v2[m]), color='#dddddd', lw=0, label='$g_M=2.85$-$3.10$' if i == 0 else None)
    ax.plot(E, v0[m], '--', color=C_NNLO, lw=1.0, label='NNLO prediction' if i == 0 else None)
    ax.plot(E, vn[m], '-', color=C_F1, label='NNLO, $\\alpha,\\beta$ of NNLO fit ($g_M=2.97$)' if i == 0 else None)
    ax.plot(E, vm[m], '-.', color='#2e8b57', lw=1.0, label='same, $g_M=2.90$' if i == 0 else None)
    if vnv is not None:
        ax.plot(E, vnv[m], ':', color='#8c564b', lw=1.0, label='$g_M=2.97$, no vertex corr.' if i == 0 else None)
    if fr == 'CM':
        sel = df[df.dobs.isin(['CS', 'CSCM']) & (np.abs(df.theta_cm - a) <= 3) & (df.E > 140) & (df.E <= 352) & ~df.dataset.isin(FWD)]
        seen = plot_data(ax, sel, conv=vsE(to_cm), seen=seen, ms=3)
        ax.set_ylabel('$d\\sigma/d\\Omega_{\\rm cm}$ [nb/sr]')
    else:
        sel = df[df.dobs.isin(['CS', 'CSCM']) & (np.abs(df.theta_lab - a) <= 3.5) & (df.E > 140) & (df.E <= 352) & ~df.dataset.isin(FWD)]
        seen = plot_data(ax, sel, conv=vsE(to_lab), seen=seen, ms=3)
        ax.set_ylabel('$d\\sigma/d\\Omega_{\\rm lab}$ [nb/sr]')
    ax.text(0.04, 0.92, lab, transform=ax.transAxes, va='top')
    ax.set_xlim(140, 352)
    if i >= 2:
        ax.set_xlabel('$E_\\gamma$ [MeV]')
h, l = [], []
for ax in axs:
    hh, ll = ax.get_legend_handles_labels(); h += hh; l += ll
fig.legend(h, l, loc='lower center', ncol=4, fontsize=6.5, bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(os.path.join(FIGS, 'fig5c_delta_cs_vs_E.pdf')); plt.close(fig)

# ---------------------------------------------------------------- (vi) pulls / chi^2 per data set (main F1)
for key, fname in [('F2_d_all_thB', 'fig6a_pulls_main_F2'), ('F1_d_all_thB', 'fig6b_pulls_F1')]:
    r = FL[key]
    s = df.set_index('pid').loc[r['pids']]
    pull = np.array(r['resid']) / s.err_pt.values
    sets = list(dict.fromkeys(s.dataset))
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw=dict(width_ratios=[2.6, 1]))
    ax = axs[0]
    for j, ds in enumerate(sets):
        m = (s.dataset == ds).values
        mk, c = MARK.get(ds, ('o', 'k'))
        x = j + 0.25 * (np.random.default_rng(j).random(m.sum()) - 0.5)
        ax.plot(x, pull[m], mk, color=c, ms=2.5, alpha=0.8)
        ax.text(j, 3.9, '%.1f/%d' % (r['perset'][ds]['chi2'], r['perset'][ds]['N']), rotation=90, fontsize=5.5, ha='center', va='bottom')
    ax.axhline(0, color='k', lw=0.5); ax.axhline(3, color='0.6', lw=0.5, ls=':'); ax.axhline(-3, color='0.6', lw=0.5, ls=':')
    ax.set_xticks(range(len(sets))); ax.set_xticklabels([tex(DB.EXPT.get(k, k)) for k in sets], rotation=70, ha='right', fontsize=5.5)
    ax.set_ylabel('pull $(y-n\\,t)/\\sigma_{\\rm pt}$'); ax.set_ylim(-5.2, 6.2); ax.set_xlim(-0.7, len(sets) - 0.3)
    ax = axs[1]
    nonf = ~s.dataset.isin(FWD).values
    ax.hist(pull[nonf], bins=np.arange(-4, 4.01, 0.5), color=C_BAND, edgecolor='white', lw=0.5, label='angular data')
    xx = np.linspace(-4, 4, 200)
    ax.plot(xx, nonf.sum() * 0.5 * np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), 'k-', lw=0.8, label='N(0,1)')
    ax.set_xlabel('pull'); ax.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(os.path.join(FIGS, fname + '.pdf')); plt.close(fig)

# ---------------------------------------------------------------- (vii) forward points (paper Fig. "fwd")
# Left: dsigma(0) and Sigma_2z(0) vs energy; right: the same relative to the empirical (GHP) values, which are
# represented by cubic splines through the tabulated forward amplitudes Re f, Re g.  Curves as in Fig. cs (NNLO
# prediction solid with the natural-LEC band, NNLO fit dashed) plus the NNLO fit without the truncation covariance
# (dotted); in the relative panels the band is transformed like the curves (divided by the empirical cross section,
# shifted by the empirical Sigma_2z).  The thin black dashed line is the empirical input shifted by the fitted nuisances n_f, n_g.
import forward as FWDLIB                      # noqa: E402
from scipy.interpolate import CubicSpline     # noqa: E402
Ef = np.arange(40, 151, 1.0)
g = grid_rows([('fw', e, 0.0, o, 'LAB') for o in ['CS', 'S2Z'] for e in Ef])
cv = curves(g, 'fig_fwd', {'F2': P_F2, 'F2n': np.array(FL['F2_d_all']['P'])}, band=False)
Eg, fv, dfv, gv, dgv = FWDLIB.amplitudes()
emp = dict(zip(['CS', 'S2Z'], FWDLIB.obs(CubicSpline(Eg, fv)(Ef), CubicSpline(Eg, gv)(Ef))))
# GHP points shifted by the fitted nuisances of the main fit: y - lam_f s_f - lam_g s_g (model = n t + lam.s)
rm = FL['F2_d_all_thB']
dg = df[df.dataset.isin(FWD)].sort_values(['dataset', 'E']).reset_index(drop=True)
_, sh, _ = FWDLIB.model(dg, mode='nuisance')
dg['shifted'] = dg.value.values - rm['lam_shift'][0] * sh[0] - rm['lam_shift'][1] * sh[1]
fig, axs = plt.subplots(2, 2, figsize=(W2, 4.3), sharex=True, layout='constrained')
for row, (o, yl, yr) in enumerate([('CS', '$d\\sigma/d\\Omega_{\\rm lab}(0^\\circ)$ [nb/sr]', 'relative to empirical'),
                                    ('S2Z', '$\\Sigma_{2z}(0^\\circ)$', '$\\Sigma_{2z}-\\Sigma_{2z}^{\\rm emp}$')]):
    m = (g.dobs == o).values
    d = dg[dg.dataset == ('GHP_fwd_CS' if o == 'CS' else 'GHP_fwd_S2Z')]
    used = (d.E <= (145 if PROFILE == 'v2' else 135)).values
    e_emp = emp[o]
    ref = np.interp(d.E.values, Ef, e_emp)
    if o == 'CS':      # relative deviation for the cross section, difference for the asymmetry
        rel = lambda y, at=None: y / (e_emp if at is None else at) - 1
        yd, ed = d.value.values / ref - 1, d.err_pt.values / ref
        ys = d.shifted.values / ref - 1
    else:
        rel = lambda y, at=None: y - (e_emp if at is None else at)
        yd, ed = d.value.values - ref, d.err_pt.values
        ys = d.shifted.values - ref
    if o == 'CS':
        print('fig7_forward dsigma(0): LEC band max %.3f nb/sr = %.2f%% of the prediction (at %.0f MeV); |pred-emp| max %.1f%%'
              % (cv['lec'][m].max(), 100 * (cv['lec'][m] / cv['NNLO'][m]).max(), Ef[np.argmax(cv['lec'][m] / cv['NNLO'][m])],
                 100 * np.abs(cv['NNLO'][m] / e_emp - 1).max()))
    else:
        print('fig7_forward Sigma_2z(0): LEC band max %.4f (at %.0f MeV); |pred-emp| max %.4f'
              % (cv['lec'][m].max(), Ef[np.argmax(cv['lec'][m])], np.abs(cv['NNLO'][m] - e_emp).max()))
    for ax, f in [(axs[row, 0], lambda y: y), (axs[row, 1], rel)]:
        ax.fill_between(Ef, f(cv['NNLO'][m] - cv['lec'][m]), f(cv['NNLO'][m] + cv['lec'][m]), color=C_BAND, alpha=0.6, lw=0)
        ax.plot(Ef, f(cv['NNLO'][m]), '-', color=C_NNLO, lw=1.2)
        ax.plot(Ef, f(cv['F2'][m]), '--', color=C_F1, lw=1.0)
        ax.plot(Ef, f(cv['F2n'][m]), ':', color=C_F1, lw=1.1, label='NNLO fit, no th.\\ cov.' if (row, f is not rel) == (0, True) else None)
        ax.axvline(144.7, color='0.6', lw=0.5, ls=':')
    ax = axs[row, 0]
    ax.errorbar(d.E[used], d.value[used], d.err_pt[used], fmt='o', color='k', ms=3, lw=0.7, label='forward points (used)' if row == 0 else None)
    ax.errorbar(d.E[~used], d.value[~used], d.err_pt[~used], fmt='o', mfc='white', color='k', ms=3, lw=0.7,
                label='forward point (not used)' if row == 0 else None)
    ax.set_ylabel(yl)
    ax = axs[row, 1]
    ax.axhline(0, color='k', lw=0.6)
    ax.plot(d.E.values, ys, '--', color='k', lw=0.7, label='empirical, shifted by fitted $n_f$, $n_g$' if row == 0 else None)
    ax.errorbar(d.E[used], yd[used], ed[used], fmt='o', color='k', ms=3, lw=0.7)
    ax.errorbar(d.E[~used], yd[~used], ed[~used], fmt='o', mfc='white', color='k', ms=3, lw=0.7)
    ax.set_ylabel(yr)
    ax.set_ylim(*((-0.065, 0.065) if o == 'CS' else (-0.06, 0.06)))
for ax in axs[1]:
    ax.set_xlabel('$E_\\gamma$ [MeV]'); ax.set_xlim(40, 152)
h, l = [], []
for ax in (axs[0, 0], axs[0, 1]):
    hh, ll = ax.get_legend_handles_labels(); h += hh; l += ll
order = ['NNLO fit, no th.\\ cov.', 'empirical, shifted by fitted $n_f$, $n_g$', 'forward points (used)', 'forward point (not used)']
h, l = th_handles()[:2] + [h[l.index(k)] for k in order], TH_LAB[:2] + order    # column-wise pairs: theory, theory, data (no PWA curve here)
fig.legend(h, l, loc='outside upper center', ncol=3, columnspacing=1.5, handletextpad=0.5)
save_exact(fig, 'fig7_forward')

# ---------------------------------------------------------------- (viii) energy-cut scan
# NNLO fit (forward points as data) with and without the theoretical covariance; the Baldin-sum-rule fits (F1) are no
# longer shown (dropped 2026-10-16; their scan results stay in fits_L.json, keys scan_F1_*).
fig, axs = plt.subplots(2, 1, figsize=(W1, 4.0), sharex=True, layout='constrained')   # 4.0 in: 3-row legend (was 4.3 in with 5 rows)
for th_, ls, mk in [('noth', ':', 'o'), ('thB', '-', 's')]:
    Ems = [90, 100, 110, 120, 130, 140, 150, 160, 170]
    rr = [FL['scan_F2_%s_E%d' % (th_, e)] for e in Ems]
    off = -0.6 if th_ == 'noth' else 0.6
    for ax, i in zip(axs, [0, 1]):
        ax.errorbar(np.array(Ems) + off, [r['P'][i] for r in rr], [r['dP'][i] for r in rr], fmt=mk + ls, color=C_F1, ms=3,
                    mfc='white' if th_ == 'noth' else C_F1, lw=0.8,
                    label='NNLO fit, %s' % ('no th.\\ cov.' if th_ == 'noth' else 'with th.\\ cov.'))
axs[0].set_ylabel('$\\alpha_{E1}$ [$10^{-4}$ fm$^3$]'); axs[1].set_ylabel('$\\beta_{M1}$ [$10^{-4}$ fm$^3$]')
for ax, p0 in zip(axs, PREF_NNLO[:2]):       # NNLO prediction with the natural-LEC band, as in Figs. cs, sigma3
    ax.axhspan(p0 - DLEC, p0 + DLEC, color=C_BAND, alpha=0.6, lw=0, zorder=0)
    ax.axhline(p0, color=C_NNLO, lw=1.2, zorder=1)
    ax.axvline(144.7, color='0.6', lw=0.5, ls=':')
axs[1].set_xlabel('$E_{\\rm max}$ [MeV]')
h, l = axs[1].get_legend_handles_labels()
# three entries in one column: two columns measure 3.9 in (> W1 = 3.4 in) even with columnspacing 1.0; the extent is checked below
leg = fig.legend(th_handles()[:1] + h, TH_LAB[:1] + l, loc='outside upper center', ncol=1, handlelength=2.6)
fig.canvas.draw()
_wl = leg.get_window_extent(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).width
print('fig8_Emax_scan: legend width %.2f in within the figure width %.2f in%s' % (_wl, W1, '' if _wl <= W1 else '  [OVERFLOW]'))
save_exact(fig, 'fig8_Emax_scan')
print('figures written to', FIGS)
