"""Forward Compton scattering off the proton, 50-450 MeV (paper figure pwa_forward, redesign of 2026-10-17; since 2026-10-18 with
both empirical evaluations): the 2026 update of the empirical forward amplitudes (public repository github.com/HiWARPs/ForwardRCS,
tables in work/data/raw/ForwardRCS2026: unpol_f.txt, pol_g.txt, obs.txt with exp+fit errors; their systematics unpol_selection.txt
(data selection, f), pol_g_pwa_model.txt (PWA model of g below 0.2 GeV, propagated to the observables), obs_selection.txt; band =
exp+fit (+) systematics in quadrature; green), and
the empirical forward amplitudes of Gryniuk, Hagelstein, Pascalutsa (GHP: PRD 92 (2015) 074031 and PRD 94 (2016) 034043,
Supplemental Material in work/data/raw/GHP2016_supplement, photoabsorption fits I and II combined as in
work/data/scripts/build_db_forward_above.py: mean, error = mean quoted error (+) half the difference in quadrature; orange),
the forward points used in the fits (work/data/db_forward_obs.csv, 50-150 MeV; since 2026-10-18 from the 2026 update), the PWA main fit of the paper at its bin
centres (legacy normalisations floated, configuration main_fits.MAIN, as in multipoles_figure.py; non-Born f, g of
make_figures.forward_model + GHP Born terms) and the parameter-free NNLO BChPT prediction (work/pwa/results/
nnlo_forward_amplitudes.csv, written by nnlo_forward.py; exact forward amplitudes of the v2 code at the NNLO static
polarisabilities).
Panels: left dsigma/dOmega_lab(0) (log) and Sigma_2z(0); right f(nu) and g(nu), Re solid, Im dashed.  GHP conventions.
Outputs: paper/figs/pwa_forward.pdf (exact double-column width W2 = 510 pt), work/pwa/figs/pwa_forward.png (preview),
         work/pwa/results/forward_figure.json (the plotted PWA points and the NNLO-empirical comparison numbers, GHP and 2026),
         a 'Figure' section appended to work/notes/nnlo_forward_amplitudes.md.
Run:  python3 work/pwa/scenarios/forward_figure.py
"""
import os
import sys
import json

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
WORK = os.path.dirname(PWA)
ROOT = os.path.dirname(WORK)
for p in (PWA, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
os.chdir(PWA)

import numpy as np                                    # noqa: E402
import pandas as pd                                   # noqa: E402
import matplotlib                                     # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402
from matplotlib.patches import Patch                  # noqa: E402
from matplotlib.ticker import MultipleLocator         # noqa: E402
import make_figures as MF                             # noqa: E402  (style, W2, forward_model, COL)
import update_fits as UF                              # noqa: E402
from main_fits import MAIN                            # noqa: E402
from sumrules import MP_GHP, ALPHA, KAPPA, HBARC_MEV, UBGEV_PER_FM   # noqa: E402
import nnlo_forward as NF                             # noqa: E402
from nnlo_forward import fmt_rows                     # noqa: E402

RAW = os.path.join(WORK, 'data', 'raw', 'GHP2016_supplement')
RAW26 = os.path.join(WORK, 'data', 'raw', 'ForwardRCS2026')
FWD = os.path.join(WORK, 'data', 'db_forward_obs.csv')
PDF = os.path.join(ROOT, 'paper', 'figs', 'pwa_forward.pdf')
PNG = os.path.join(PWA, 'figs', 'pwa_forward.png')
OUTJSON = os.path.join(PWA, 'results', 'forward_figure.json')
EMIN, EMAX, ETH = 50.0, 450.0, 150.0
COL = dict(ghp='#d9800b', band='#d9800b', nnlo='#444444', pwa=MF.COL['main'],    # GHP in orange, distinct from the grey NNLO curve (2026-10-17)
           r26='#1b9e77', band26='#1b9e77')                                      # 2026 update in green (2026-10-18)
LEG = {'Oxley': 0.10, 'Pugh': 0.10, 'Hyman': 0.10, 'Bernardini': 0.10, 'Goldansky': 0.08}


# ---------------------------------------------------------------- GHP supplement
def read_sup(name):
    """'nu [GeV]  columns...' -> (E [MeV], array of the remaining columns)"""
    rows = []
    for line in open(os.path.join(RAW, name)):
        p = line.split()
        if not p:
            continue
        try:
            rows.append([float(x) for x in p])
        except ValueError:
            continue
    a = np.array(rows)
    return np.round(a[:, 0] * 1e3, 3), a[:, 1:]


def combine(y1, e1, y2, e2):
    """fits I and II: mean, error = mean quoted error (+) half-difference in quadrature"""
    return 0.5 * (y1 + y2), np.sqrt((0.5 * (e1 + e2))**2 + (0.5 * (y1 - y2))**2)


def ghp(Emin=EMIN, Emax=EMAX):
    E, f1 = read_sup('f_1.txt'); E2, f2 = read_sup('f_2.txt'); Eg, g = read_sup('g.txt')
    Ed1, d1 = read_sup('ds_1.txt'); Ed2, d2 = read_sup('ds_2.txt'); Es1, s1 = read_sup('Sigma2z_ds_1.txt'); Es2, s2 = read_sup('Sigma2z_ds_2.txt')
    for x in (E2, Eg, Ed1, Ed2, Es1, Es2):
        assert np.allclose(x, E)
    sel = (E >= Emin - 1e-6) & (E <= Emax + 1e-6)
    out = dict(E=E[sel])
    out['Ref'], out['dRef'] = combine(f1[sel, 0], f1[sel, 1], f2[sel, 0], f2[sel, 1])
    out['Imf'], out['dImf'] = combine(f1[sel, 2], f1[sel, 3], f2[sel, 2], f2[sel, 3])
    out['Reg'], out['dReg'], out['Img'], out['dImg'] = g[sel, 0], g[sel, 1], g[sel, 2], g[sel, 3]       # one table only
    out['cs'], out['dcs'] = combine(d1[sel, 0], d1[sel, 1], d2[sel, 0], d2[sel, 1])
    r1, r2 = s1[sel, 0] / d1[sel, 0], s2[sel, 0] / d2[sel, 0]                   # Sigma_2z = (-2 Re f g*)/dsigma per fit
    e1 = np.sqrt((s1[sel, 1] / d1[sel, 0])**2 + (r1 * d1[sel, 1] / d1[sel, 0])**2)
    e2 = np.sqrt((s2[sel, 1] / d2[sel, 0])**2 + (r2 * d2[sel, 1] / d2[sel, 0])**2)
    out['s2z'], out['ds2z'] = combine(r1, e1, r2, e2)
    return out


def rcs2026(Emin=EMIN, Emax=EMAX):
    """the 2026 update (same keys as ghp()): values of unpol_f.txt, pol_g.txt, obs.txt; errors = exp+fit (+) data selection of f
    (+) PWA model of g (signed modes, propagated linearly to dsigma(0) and Sigma_2z(0)), in quadrature"""
    def rd(fn):
        rows = [[float(x) for x in l.split()] for l in open(os.path.join(RAW26, fn)) if l.strip() and l.split()[0][0].isdigit()]
        a = np.array(rows)
        return {round(1e3 * x, 3): r for x, r in zip(a[:, 0], a[:, 1:])}
    f, fs, g, gp, o, os_ = (rd(n) for n in ('unpol_f.txt', 'unpol_selection.txt', 'pol_g.txt', 'pol_g_pwa_model.txt', 'obs.txt',
                                            'obs_selection.txt'))
    E = np.array([e for e in sorted(o) if Emin - 1e-6 <= e <= Emax + 1e-6 and all(e in t for t in (f, fs, g, gp, os_))])
    A = lambda t, j: np.array([t[e][j] for e in E])
    K = 1e3 / 389.379
    Rf, dRf, If, dIf, Rg, dRg, Ig, dIg = A(f, 0), A(f, 1), A(f, 2), A(f, 3), A(g, 0), A(g, 1), A(g, 2), A(g, 3)
    sRf, sIf, pRg, pIg = A(fs, 0), A(fs, 1), A(gp, 0), A(gp, 1)
    cs, dcs, s2, ds2 = A(o, 0), A(o, 1), A(o, 4), A(o, 5)
    pcs = 2 * K * (Rg * pRg + Ig * pIg)                              # PWA-model mode of g -> dsigma(0), -2 Re f g*, Sigma_2z
    pn = -2 * K * (Rf * pRg + If * pIg)
    ps2 = pn / cs - s2 * pcs / cs
    assert np.allclose(K * (Rf**2 + If**2 + Rg**2 + Ig**2), cs, rtol=2e-6)
    return dict(E=E, Ref=Rf, dRef=np.hypot(dRf, sRf), Imf=If, dImf=np.hypot(dIf, sIf), Reg=Rg, dReg=np.hypot(dRg, pRg),
                Img=Ig, dImg=np.hypot(dIg, pIg), cs=cs, dcs=np.sqrt(dcs**2 + A(os_, 0)**2 + pcs**2),
                s2z=s2, ds2z=np.sqrt(ds2**2 + A(os_, 2)**2 + ps2**2))


# ---------------------------------------------------------------- PWA main fit at the bin centres
def load_leg(per_beamtime=False):
    d = UF.load(per_beamtime=per_beamtime)
    for s, e in LEG.items():
        m = d.dataset == s
        d.loc[m, 'norm_group'] = s
        d.loc[m, 'norm_err'] = e
    return d


def born_ghp(E):
    """GHP Born terms (fm) with the constants used for the Born subtraction of the fitted forward amplitudes"""
    Mfm = MP_GHP / HBARC_MEV
    return -ALPHA / Mfm, -ALPHA * KAPPA**2 * (np.asarray(E, float) / HBARC_MEV) / (2 * Mfm**2)


def pwa_points():
    m = UF.build_model(UF.select(load_leg()), **MAIN)
    r = m.fit()
    _ref = json.load(open(os.path.join(PWA, 'results', 'main_fits_legacynorm.json')))['MAIN']
    assert abs(float(r.der['alpha_m_beta']) - _ref['val']['alpha_m_beta']) < 1e-3 and int(r.npts) == _ref['npts'], (r.der, _ref['val'])
    rows = []
    for ib, Eb in enumerate(m.bins):
        fb, gb, Jf, Jg, C_ = MF.forward_model(m, r, ib)
        fB, gB = born_ghp(Eb)
        f, g = fb + fB, gb + gB
        J = np.vstack([Jf, Jg]); Cfg = J @ C_ @ J.T
        den = f**2 + g**2
        cs, s2 = den * 1e7, -2 * f * g / den
        dcs = np.array([2 * f, 2 * g]) * 1e7
        ds = np.array([-2 * g / den + 4 * f**2 * g / den**2, -2 * f / den + 4 * f * g**2 / den**2])
        rows.append(dict(E=float(Eb), f=f * UBGEV_PER_FM, df=np.sqrt(Cfg[0, 0]) * UBGEV_PER_FM, g=g * UBGEV_PER_FM, dg=np.sqrt(Cfg[1, 1]) * UBGEV_PER_FM,
                         cs=cs, dcs=float(np.sqrt(dcs @ Cfg @ dcs)), s2z=s2, ds2z=float(np.sqrt(ds @ Cfg @ ds))))
    return pd.DataFrame(rows), dict(chi2=float(r.chi2), npts=int(r.npts), nf=float(r.p[m.names.index('nf')]), ng=float(r.p[m.names.index('ng')]))


# ---------------------------------------------------------------- figure
def save_both(fig):
    with plt.rc_context({'savefig.bbox': 'standard'}):
        fig.savefig(PDF)
        fig.savefig(PNG, dpi=300)
    plt.close(fig)


def thr(ax):
    ax.axvline(ETH, color='k', ls=':', lw=0.6, zorder=1)


def figure(G, T, P, R=None):
    E = T.E_MeV.values
    fo = pd.read_csv(FWD)
    fig, axs = plt.subplots(2, 2, figsize=(MF.W2, 4.5), sharex=True, layout='constrained')
    (a, c), (b, d) = axs
    band = dict(color=COL['band'], alpha=0.25, lw=0, zorder=2)
    ghpl = dict(color=COL['ghp'], lw=0.9, zorder=3)
    band26 = dict(color=COL['band26'], alpha=0.35, lw=0, zorder=2.5)
    r26l = dict(color=COL['r26'], lw=0.9, zorder=3.5)
    nnl = dict(color=COL['nnlo'], lw=1.1, zorder=4)
    pts = dict(ls='none', marker='o', ms=3.2, color='k', mfc='w', mew=0.7, elinewidth=0.7, zorder=5)
    pwa = dict(ls='none', marker='s', ms=2.8, color=COL['pwa'], mfc=COL['pwa'], elinewidth=0.7, zorder=6)
    # --- (a) forward cross section
    a.fill_between(G['E'], G['cs'] - G['dcs'], G['cs'] + G['dcs'], **band)
    a.plot(G['E'], G['cs'], **ghpl)
    a.fill_between(R['E'], R['cs'] - R['dcs'], R['cs'] + R['dcs'], **band26)
    a.plot(R['E'], R['cs'], **r26l)
    a.plot(E, T.dsigma0_nb, **nnl)
    x = fo[fo.obs == 'CS']; h_pts = a.errorbar(x.E_lab_MeV, x.value, x.stat_err, **pts)
    h_pwa = a.errorbar(P.E, P.cs, P.dcs, **pwa)
    a.set_yscale('log'); a.set_ylim(0.5, 1500)
    a.set_ylabel(r'$d\sigma/d\Omega_{\rm lab}(0^\circ)$ [nb/sr]')
    # --- (b) Sigma_2z
    b.fill_between(G['E'], np.clip(G['s2z'] - G['ds2z'], -1, 1), np.clip(G['s2z'] + G['ds2z'], -1, 1), **band)   # |Sigma_2z| <= 1
    b.plot(G['E'], G['s2z'], **ghpl)
    b.fill_between(R['E'], np.clip(R['s2z'] - R['ds2z'], -1, 1), np.clip(R['s2z'] + R['ds2z'], -1, 1), **band26)
    b.plot(R['E'], R['s2z'], **r26l)
    b.plot(E, T.Sigma2z0, **nnl)
    y = fo[fo.obs == 'S2Z']; b.errorbar(y.E_lab_MeV, y.value, y.stat_err, **pts)
    b.errorbar(P.E, P.s2z, P.ds2z, **pwa)
    b.axhline(0, color='k', lw=0.4, zorder=1)
    b.set_ylim(-1.1, 1.0); b.yaxis.set_major_locator(MultipleLocator(0.5))
    b.set_ylabel(r'$\Sigma_{2z}(0^\circ)$')
    # --- (c) f(nu), (d) g(nu): Re solid, Im dashed
    for ax, re, dre, im, dim, nre, nim, pv, pe, lab, k in (
            (c, G['Ref'], G['dRef'], G['Imf'], G['dImf'], T.Re_f, T.Im_f, P.f, P.df, 'f', 'f'),
            (d, G['Reg'], G['dReg'], G['Img'], G['dImg'], T.Re_g, T.Im_g, P.g, P.dg, 'g', 'g')):
        ax.fill_between(G['E'], re - dre, re + dre, **band)
        ax.fill_between(G['E'], im - dim, im + dim, **band)
        ax.plot(G['E'], re, **ghpl)
        ax.plot(G['E'], im, ls='--', **ghpl)
        ax.fill_between(R['E'], R['Re' + k] - R['dRe' + k], R['Re' + k] + R['dRe' + k], **band26)
        ax.fill_between(R['E'], R['Im' + k] - R['dIm' + k], R['Im' + k] + R['dIm' + k], **band26)
        ax.plot(R['E'], R['Re' + k], **r26l)
        ax.plot(R['E'], R['Im' + k], ls='--', **r26l)
        ax.plot(E, nre, **nnl)
        ax.plot(E, nim, ls='--', **nnl)
        ax.errorbar(P.E, pv, pe, **pwa)
        ax.axhline(0, color='k', lw=0.4, zorder=1)
        ax.set_ylabel(rf'${lab}(\nu)$ [$\mu$b$\,$GeV]')
    c.set_ylim(-9, 16); c.yaxis.set_major_locator(MultipleLocator(5))
    d.set_ylim(-7, 4); d.yaxis.set_major_locator(MultipleLocator(2))
    c.text(395, 12.6, r'Im$\,f$', ha='center', va='center')
    c.text(300, -6.6, r'Re$\,f$', ha='center', va='center')
    d.text(195, 2.3, r'Im$\,g$', ha='center', va='center')
    d.text(210, -4.8, r'Re$\,g$', ha='center', va='center')
    for ax in axs.ravel():
        thr(ax)
        ax.set_xlim(EMIN, EMAX)
        ax.xaxis.set_major_locator(MultipleLocator(100)); ax.xaxis.set_minor_locator(MultipleLocator(25))
    for ax in (b, d):
        ax.set_xlabel(r'$E_\gamma$ [MeV]')
    h = [(Patch(facecolor=COL['band'], alpha=0.25, lw=0), Line2D([], [], color=COL['ghp'], lw=0.9)),
         (Patch(facecolor=COL['band26'], alpha=0.35, lw=0), Line2D([], [], color=COL['r26'], lw=0.9)),
         h_pts, h_pwa, Line2D([], [], color=COL['nnlo'], lw=1.1)]
    l = ['empirical, GHP 2015/16', 'empirical, 2026 update', 'forward points', 'PWA fit', r'NNLO B$\chi$PT']
    leg = fig.legend(h, l, loc='outside upper center', ncol=5, columnspacing=1.6, handlelength=2.2)
    fig.canvas.draw()
    wl = leg.get_window_extent(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).width
    print('legend width %.2f in (figure %.2f in)%s' % (wl, MF.W2, '' if wl <= MF.W2 else '  [OVERFLOW]'))
    save_both(fig)


def comparison(G, T):
    """maximal |NNLO - GHP| of f and g (Re, Im, modulus of the complex difference) below and above threshold, in mub GeV,
    and the maximal deviations in units of the GHP error; also dsigma(0) and Sigma_2z(0)."""
    Ei = T.E_MeV.values
    gi = {k: np.interp(Ei, G['E'], G[k]) for k in G if k != 'E'}
    f = T.Re_f.values + 1j * T.Im_f.values; g = T.Re_g.values + 1j * T.Im_g.values
    fe = gi['Ref'] + 1j * gi['Imf']; ge = gi['Reg'] + 1j * gi['Img']
    out = {}
    for tag, sel in (('below', Ei <= 150), ('above', Ei > 150)):
        d = {}
        for name, a, b, e in (('Re_f', f.real, fe.real, gi['dRef']), ('Im_f', f.imag, fe.imag, gi['dImf']),
                              ('Re_g', g.real, ge.real, gi['dReg']), ('Im_g', g.imag, ge.imag, gi['dImg'])):
            i = np.argmax(abs(a - b)[sel]); j = np.where(sel)[0][i]
            d[name] = dict(max_abs=float(abs(a - b)[j]), E=float(Ei[j]), nnlo=float(a[j]), ghp=float(b[j]), ghp_err=float(e[j]),
                           max_pull=float(np.max(abs(a - b)[sel] / np.maximum(e[sel], 1e-9))))
        for name, a, b in (('abs_f', f, fe), ('abs_g', g, ge)):
            i = np.argmax(abs(a - b)[sel]); j = np.where(sel)[0][i]
            d[name] = dict(max_abs=float(abs(a - b)[j]), E=float(Ei[j]))
        for name, a, b, e in (('cs_rel', T.dsigma0_nb.values, gi['cs'], gi['dcs']), ('s2z', T.Sigma2z0.values, gi['s2z'], gi['ds2z'])):
            dev = (a / b - 1) if name == 'cs_rel' else (a - b)
            i = np.argmax(abs(dev)[sel]); j = np.where(sel)[0][i]
            d[name] = dict(max_abs=float(abs(dev[j])), E=float(Ei[j]), nnlo=float(a[j]), ghp=float(b[j]), ghp_err=float(e[j]))
        out[tag] = d
    # peak positions and heights of Im f
    out['Imf_peak'] = dict(nnlo=dict(E=float(Ei[np.argmax(f.imag)]), value=float(np.max(f.imag))),
                           ghp=dict(E=float(G['E'][np.argmax(G['Imf'])]), value=float(np.max(G['Imf']))))
    out['Img_min'] = dict(nnlo=dict(E=float(Ei[np.argmin(g.imag)]), value=float(np.min(g.imag))),
                          ghp=dict(E=float(G['E'][np.argmin(G['Img'])]), value=float(np.min(G['Img']))))
    return out


def note_section(G, T, cmp_, info):
    """append (or replace) the figure section of work/notes/nnlo_forward_amplitudes.md"""
    E = T.E_MeV.values
    gi = {k: np.interp(E, G['E'], G[k]) for k in G if k != 'E'}
    rows = []
    for e in (50, 100, 130, 150, 200, 250, 300, 325, 350, 400, 450):
        i = int(np.argmin(abs(E - e)))
        rows.append([f'{e}', f'{T.Re_f[i]:+.3f} / {gi["Ref"][i]:+.3f} ({gi["dRef"][i]:.3f})', f'{T.Im_f[i]:+.2f} / {gi["Imf"][i]:+.2f} ({gi["dImf"][i]:.2f})',
                     f'{T.Re_g[i]:+.3f} / {gi["Reg"][i]:+.3f} ({gi["dReg"][i]:.3f})', f'{T.Im_g[i]:+.2f} / {gi["Img"][i]:+.2f} ({gi["dImg"][i]:.2f})',
                     f'{T.dsigma0_nb[i]:.2f} / {gi["cs"][i]:.2f} ({100 * (T.dsigma0_nb[i] / gi["cs"][i] - 1):+.0f} %)', f'{T.Sigma2z0[i]:+.3f} / {gi["s2z"][i]:+.3f} ({gi["ds2z"][i]:.3f})'])
    b, a = cmp_['below'], cmp_['above']
    txt = f"""
## Figure `paper/figs/pwa_forward.pdf` (script `work/pwa/scenarios/forward_figure.py`)

2 x 2 panels, double-column width, 50-450 MeV: left dsigma/dOmega_lab(0) (log) and Sigma_2z(0); right f(nu) and g(nu) (Re solid,
Im dashed).  Green: the 2026 update (ForwardRCS; band = exp+fit (+) data selection (+) PWA model of g); the NNLO-vs-2026 numbers are in
work/pwa/results/forward_figure.json ('nnlo_vs_2026').  Orange: GHP empirical amplitudes with the band of the combined fits I/II (mean, error = mean quoted (+) half-difference;
Sigma_2z formed per fit from the supplement tables -2 Re f g* and dsigma, its band clipped to |Sigma_2z| <= 1; g has a single table).
Open circles: the 12 forward points of the fits (db_forward_obs.csv; 2026 update since 2026-10-18).  Blue squares: PWA main fit of the paper (legacy normalisations
floated, chi^2 = {info['chi2']:.1f} for N = {info['npts']}, n_f = {info['nf']:.2f}, n_g = {info['ng']:.2f}) at its bin centres, non-Born f, g from
make_figures.forward_model + GHP Born terms, errors from the fit covariance.  Dark grey: NNLO BChPT (primary columns of the CSV).
Dotted line: pion-production threshold (150 MeV).  PNG preview work/pwa/figs/pwa_forward.png; plotted points and comparison numbers
in work/pwa/results/forward_figure.json.

### NNLO vs GHP (mub GeV; NNLO / GHP (GHP error))

{fmt_rows(rows, ['E [MeV]', 'Re f', 'Im f', 'Re g', 'Im g', 'dsigma(0) [nb/sr] (dev.)', 'Sigma_2z(0)'])}

Maximal deviations |NNLO - GHP|:
- below threshold (50-150 MeV): |f| {b['abs_f']['max_abs']:.3f} at {b['abs_f']['E']:.0f} MeV (Re f {b['Re_f']['nnlo']:+.3f} vs {b['Re_f']['ghp']:+.3f} +- {b['Re_f']['ghp_err']:.3f});
  |g| {b['abs_g']['max_abs']:.3f} at {b['abs_g']['E']:.0f} MeV (Re g {b['Re_g']['nnlo']:+.3f} vs {b['Re_g']['ghp']:+.3f} +- {b['Re_g']['ghp_err']:.3f});
  dsigma(0) {100 * b['cs_rel']['max_abs']:.0f} % at {b['cs_rel']['E']:.0f} MeV (the minimum of the cross section; 9 % at 135 MeV, 1 % at 50 MeV);
  Sigma_2z {b['s2z']['max_abs']:.3f} at {b['s2z']['E']:.0f} MeV.
- above threshold (150-450 MeV): |f| {a['abs_f']['max_abs']:.2f} at {a['abs_f']['E']:.0f} MeV (Re f {a['Re_f']['nnlo']:+.2f} vs {a['Re_f']['ghp']:+.2f} +- {a['Re_f']['ghp_err']:.2f} at {a['Re_f']['E']:.0f} MeV;
  Im f {a['Im_f']['nnlo']:+.2f} vs {a['Im_f']['ghp']:+.2f} +- {a['Im_f']['ghp_err']:.2f} at {a['Im_f']['E']:.0f} MeV);
  |g| {a['abs_g']['max_abs']:.2f} at {a['abs_g']['E']:.0f} MeV (Re g {a['Re_g']['nnlo']:+.2f} vs {a['Re_g']['ghp']:+.2f} +- {a['Re_g']['ghp_err']:.2f} at {a['Re_g']['E']:.0f} MeV;
  Im g {a['Im_g']['nnlo']:+.2f} vs {a['Im_g']['ghp']:+.2f} +- {a['Im_g']['ghp_err']:.2f} at {a['Im_g']['E']:.0f} MeV);
  dsigma(0) {100 * a['cs_rel']['max_abs']:.0f} % at {a['cs_rel']['E']:.0f} MeV; Sigma_2z {a['s2z']['max_abs']:.3f} at {a['s2z']['E']:.0f} MeV.
- Im f peak: NNLO {cmp_['Imf_peak']['nnlo']['value']:.2f} at {cmp_['Imf_peak']['nnlo']['E']:.0f} MeV, GHP {cmp_['Imf_peak']['ghp']['value']:.2f} at {cmp_['Imf_peak']['ghp']['E']:.0f} MeV;
  Im g minimum: NNLO {cmp_['Img_min']['nnlo']['value']:.2f} at {cmp_['Img_min']['nnlo']['E']:.0f} MeV, GHP {cmp_['Img_min']['ghp']['value']:.2f} at {cmp_['Img_min']['ghp']['E']:.0f} MeV.
(For scale: GHP's fits I and II differ by up to 0.9 (Re f) and 1.4 (Im f) mub GeV above threshold.)
"""
    marker = '\n## Figure `paper/figs/pwa_forward.pdf`'
    old = open(NF.NOTE).read() if os.path.exists(NF.NOTE) else ''
    if marker in old:
        old = old[:old.index(marker)]
    open(NF.NOTE, 'w').write(old.rstrip('\n') + '\n' + txt)


def main():
    G = ghp()
    R = rcs2026()
    T = NF.table()
    P, info = pwa_points()
    figure(G, T, P, R)
    cmp_ = comparison(G, T)
    cmp26 = comparison(R, T)
    json.dump(dict(pwa_points=P.round(6).to_dict('records'), pwa_fit=info, nnlo_vs_ghp=cmp_, nnlo_vs_2026=cmp26,
                   note=__doc__.split('Outputs:')[0]), open(OUTJSON, 'w'), indent=1)
    note_section(G, T, cmp_, info)
    print('NNLO vs 2026 update:')
    for tag in ('below', 'above'):
        print(tag, 'threshold:')
        for k, v in cmp26[tag].items():
            print('  ', k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items()})
    print('Im f peak:', cmp26['Imf_peak'], ' Im g min:', cmp26['Img_min'])
    print('PWA fit:', info)
    print(P.round(4).to_string(index=False))
    for tag in ('below', 'above'):
        print(tag, 'threshold:')
        for k, v in cmp_[tag].items():
            print('  ', k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items()})
    print('Im f peak:', cmp_['Imf_peak'], ' Im g min:', cmp_['Img_min'])
    print('written:', PDF, PNG, OUTJSON)


if __name__ == '__main__':
    main()
