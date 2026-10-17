"""Step 10: predictions of the double-polarization asymmetries Sigma_2x, Sigma_2z BELOW the pion-production
threshold (Sec. IV G, Fig. predS2, and Appendix I, Table predS2, of the paper: planned MAMI measurement with an
active polarized target, Ref. Mornacchi:PhD).

Outputs
  paper/figs/fig_pred_S2_diff.pdf            Fig. predS2 (double column): Sigma_2x (upper row), Sigma_2z (lower row)
                                             versus theta_lab at E_lab = 100, 120, 140 MeV, RELATIVE TO THE NNLO PREDICTION
                                             (copy: work/fits/figs/fig10b_pred_S2_diff.pdf)
  work/fits/figs/fig_pred_S2_abs_unused.pdf  the same as absolute asymmetries (not used in the paper: the curves coincide)
  work/fits/results/pred_S2_table.tex        LaTeX tabular (Table predS2): NNLO-prediction values at theta_lab = 60, 90,
                                             120, 150 deg with the uncertainty (truncation (+) LECs) and the shift for
                                             gamma_pi + 2;  .md: the same in markdown;  .csv: all columns
  work/fits/results/pred_S2_table_fit.*      the previous table relative to the NNLO fit (uncertainty = truncation (+)
                                             fit, shift of the fit for gamma_pi + 2), regenerated for reference
  work/fits/results/pred_S2_curves.csv       all curves and uncertainties on the grid (reproducibility)
  work/fits/results/pred_S2_sensitivity.csv  d Sigma / d P_k at the NNLO prediction (theta_lab = 90, 105, 120 deg)
  work/fits/results/pred_S2_checks.json      the checks and the numbers quoted in Sec. IV G (printed at the end)

Curves (per panel)
  (1) NNLO prediction : P = PREF_NNLO (no free parameter; physical masses, vertex corrections, LMP15 constants) -- the
                        REFERENCE (zero line of the difference figure); band = quadrature of
                        - the truncation uncertainty (scheme B, default settings of theory.truncation, reference value 1
                          for the asymmetries; orders y0, y2, y3 at the NNLO static polarizabilities) and
                        - the uncertainty from the O(p^4) LECs of natural size, delta alpha_E1 = delta beta_M1 = +-0.7,
                          propagated by central differences (exactly as curves(...)['lec'] of s06_figures.py does for
                          the cross sections of Figs. cs and sigma3);
  (2) NNLO fit        : P from fits_L.json 'F2_d_all_thB' (alpha_E1, beta_M1 fitted, spin polarizabilities at their NNLO
                        values); its own uncertainty (truncation (+) fit, covP) is only used for the reference table;
  (3) NNLO prediction with gamma_pi + 2: (1) with gM1M1 + 1 and gE1E1 - 1 (gamma_0 unchanged);
  (4) forward points  : the empirical Sigma_2z(0) of Gryniuk-Hagelstein-Pascalutsa (database rows 'GHP_fwd_S2Z',
                        err_pt = quadrature propagation of the quoted uncertainties of Re f, Re g, strongly correlated
                        between energies), interpolated to the panel energy with a cubic spline through the 12 tabulated
                        values and the error linearly (as forward_point() of s06_figures.py does for the forward cross
                        section); drawn at theta = 0 in the Sigma_2z panels (Sigma_2x vanishes there).
Observable codes 'S2X', 'S2Z' of the driver, LAB frame (angle = lab scattering angle), i.e. the conventions of the
A2 Delta-region data (work/code_v2/driver/README.md Sect. 3):
  Sigma_2x = (sigma_R - sigma_L)/(sigma_R + sigma_L), target spin along +x (in the scattering plane, perpendicular to
             the beam, towards the scattered photon);  Sigma_2z: target spin along +z (beam direction);
  R/L = photon helicity +1/-1 (Lensky-McGovern-Pascalutsa 2015, Babusci et al. 1998; A2: Martel 2015, Paudyal 2020).
Not part of the `make` pipeline; run by hand:  .venv/bin/python scripts/s10_fig_pred_S2.py
History: until 2026-10-15 the reference was the NNLO fit (band = truncation (+) fit) and a curve with the gamma's of the
Delta-region fit was drawn; that version is kept in work/backups/*_fitref.*.
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
ROOT = os.path.dirname(os.path.dirname(FITS))
PAPER_FIGS = os.path.join(ROOT, 'paper', 'figs')
os.makedirs(RES, exist_ok=True)

C_NNLO, C_F2, C_BAND, C_GPI = '#444444', '#1f4e9c', '#9fb6dd', '#2e8b57'   # prediction, fit, band, gamma_pi + 2
DLEC = 0.7                           # natural size of the O(p^4) LECs (as s06_figures.DLEC: band of Figs. cs, sigma3)
W2 = 510 / 72.27                     # double-column width [in] (revtex4-2 prd), included without scaling
TAG = 'pred_S2'

# ---------------------------------------------------------------- parameter sets
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
fit = FL['F2_d_all_thB']
P_F2 = np.array(fit['P']); COV_F2 = np.array(fit['covP'])
assert fit['free'] == ['alphaE1', 'betaM1'], fit['free']
D_GPI = np.array([0.0, 0.0, -1.0, +1.0, 0.0, 0.0])       # gE1E1 - 1, gM1M1 + 1: gamma_pi + 2, gamma_0 unchanged
P_GPI = PREF_NNLO + D_GPI                                # the prediction with gamma_pi + 2 (dotted curve, table)
P_F2_GPI = P_F2 + D_GPI                                  # the fit with gamma_pi + 2 (reference table only)
for P0, P1 in ((PREF_NNLO, P_GPI), (P_F2, P_F2_GPI)):
    assert abs(gammapi(P1) - gammapi(P0) - 2.0) < 1e-12 and abs(gamma0(P1) - gamma0(P0)) < 1e-12

# ---------------------------------------------------------------- grid
ENERGIES = [100.0, 120.0, 140.0]              # panels
E_ALL = ENERGIES + [130.0]                    # 130 MeV only for the comparison with work/delta_plots (fig5)
TH = np.arange(0.0, 180.01, 5.0)              # 37 lab angles
OBS = ['S2X', 'S2Z']
TH_TAB = [60.0, 90.0, 120.0, 150.0]
TH_SENS = [90.0, 105.0, 120.0]


def grid_rows(spec):
    """spec: list of (tag, E, theta, dobs, frame) -> DataFrame usable by outputs()/truncation_for() (as s06_figures)"""
    rows = []
    for tag, E, th, o, fr in spec:
        thl = th if fr == 'LAB' else float(DB.th_cm2lab(E, th))
        thc = th if fr == 'CM' else float(DB.th_lab2cm(E, th))
        rows.append(dict(pid='%s_%s_%.3f_%.3f_%s' % (tag, o, E, th, fr), dataset='grid', E=E, theta=th, theta_lab=thl,
                         theta_cm=thc, dobs=o, frame=fr, value=0.0, err_pt=1.0, norm=0.0, flags=''))
    return pd.DataFrame(rows)


g = grid_rows([(TAG, E, t, o, 'LAB') for o in OBS for E in E_ALL for t in TH])
T = outputs(g, 'NNLO', tag=TAG)
assert np.allclose(T.Pref, PREF_NNLO), (T.Pref, PREF_NNLO)
g['flag'] = T.flag

# ---------------------------------------------------------------- curves
cv = {'NNLO': T.value_file.copy(), 'F2': T.predict_array(P_F2), 'gpi': T.predict_array(P_GPI), 'F2gpi': T.predict_array(P_F2_GPI)}
assert np.allclose(cv['NNLO'], T.predict_array(PREF_NNLO), atol=1e-9)
# truncation uncertainty: orders y0 (Born + pi0 + contact), y2 (NLO + contact), y3 (NNLO), all at the NNLO static
# polarizabilities, exactly as order_values() does for the database rows; reference value 1 for the asymmetries
_, tr = truncation_for(g, g, tag=TAG)
assert np.all(tr['yref'] == 1.0)
cv['sig_tr'] = tr['sigma']
# LEC uncertainty of the prediction: delta alpha_E1 = delta beta_M1 = +-DLEC, central differences, quadrature
ea, eb = np.eye(6)[0], np.eye(6)[1]
cv['d_alpha'] = (T.predict_array(PREF_NNLO + DLEC * ea) - T.predict_array(PREF_NNLO - DLEC * ea)) / 2
cv['d_beta'] = (T.predict_array(PREF_NNLO + DLEC * eb) - T.predict_array(PREF_NNLO - DLEC * eb)) / 2
cv['sig_lec'] = np.sqrt(cv['d_alpha'] ** 2 + cv['d_beta'] ** 2)
cv['sig'] = np.sqrt(cv['sig_tr'] ** 2 + cv['sig_lec'] ** 2)          # band of the prediction
# fit uncertainty of (alpha_E1, beta_M1): linear propagation with the fit covariance (reference table only)
eps = 0.01
grad = []
for k in range(2):
    dp = np.zeros(6); dp[k] = eps
    grad.append((T.predict_array(P_F2 + dp) - T.predict_array(P_F2 - dp)) / (2 * eps))
grad = np.array(grad)
cv['sig_ab'] = np.sqrt(np.einsum('ip,ij,jp->p', grad, COV_F2[:2, :2], grad))
cv['sig_fit'] = np.sqrt(cv['sig_tr'] ** 2 + cv['sig_ab'] ** 2)
cv['shift_gpi'] = cv['gpi'] - cv['NNLO']            # prediction: gamma_pi + 2 minus prediction
cv['shift_F2gpi'] = cv['F2gpi'] - cv['F2']          # fit: the same for the fit (reference table)
for k, v in cv.items():
    g[k] = v
g['y0'], g['y2'], g['y3'] = tr['y0'], tr['y2'], tr['y3']
g.to_csv(os.path.join(RES, 'pred_S2_curves.csv'), index=False, float_format='%.6g')


def sel(o, E, th=None):
    m = (g.dobs == o).values & np.isclose(g.E.values, E)
    if th is not None:
        m &= np.isclose(g.theta.values, th)
    return m


# ---------------------------------------------------------------- forward points Sigma_2z(0) of GHP
df = load_db()
fw = df[df.dataset == 'GHP_fwd_S2Z'].sort_values('E')
assert len(fw) == 12 and np.all(fw.err_syspt == 0) and np.allclose(fw.err_pt, fw.err_stat)


def forward_S2Z(E):
    """empirical Sigma_2z(0) of GHP interpolated to energy E: cubic spline through the 12 tabulated values, uncertainty
    linearly interpolated (the quoted errors) -- as forward_point() of s06_figures.py for the forward cross section."""
    return float(CubicSpline(fw.E.values, fw.value.values)(E)), float(np.interp(E, fw.E.values, fw.err_pt.values))


FWD = {}
for E in E_ALL:
    v, e = forward_S2Z(E)
    i0 = np.flatnonzero(sel('S2Z', E, 0.0))[0]
    tab = fw[np.isclose(fw.E, E)]
    FWD[E] = dict(value=v, err=e, tabulated=bool(len(tab)), linear=float(np.interp(E, fw.E.values, fw.value.values)),
                  pred=float(cv['NNLO'][i0]), sig_pred=float(cv['sig'][i0]), fit=float(cv['F2'][i0]),
                  diff=float(v - cv['NNLO'][i0]), pull_exp=float((v - cv['NNLO'][i0]) / e),
                  pull_tot=float((v - cv['NNLO'][i0]) / np.hypot(e, cv['sig'][i0])),
                  diff_fit=float(v - cv['F2'][i0]), pull_fit_exp=float((v - cv['F2'][i0]) / e))
    if len(tab):
        assert abs(v - float(tab.value.iloc[0])) < 1e-12 and abs(e - float(tab.err_pt.iloc[0])) < 1e-12


# ---------------------------------------------------------------- figure (Fig. predS2): relative to the NNLO prediction
def draw(axs, relative):
    """relative=True: Delta Sigma = Sigma - Sigma^pred (the paper figure); False: absolute asymmetries"""
    ref = cv['NNLO'] if relative else np.zeros_like(cv['NNLO'])
    h_fwd = None
    for r, (o, sym) in enumerate([('S2X', '$\\Sigma_{2x}$'), ('S2Z', '$\\Sigma_{2z}$')]):
        ymax = 0.0
        for c, E in enumerate(ENERGIES):
            ax = axs[r, c]; m = sel(o, E)
            if relative:
                ax.axhline(0, color=C_NNLO, lw=1.2, zorder=2)
            else:
                ax.axhline(0, color='0.6', lw=0.4, zorder=0)
                ax.plot(TH, cv['NNLO'][m], '-', color=C_NNLO, lw=1.2, zorder=2)
            ax.fill_between(TH, cv['NNLO'][m] - ref[m] - cv['sig'][m], cv['NNLO'][m] - ref[m] + cv['sig'][m],
                            color=C_BAND, alpha=0.6, lw=0, zorder=1)
            ax.plot(TH, cv['F2'][m] - ref[m], '--', color=C_F2, lw=1.0, zorder=3)
            ax.plot(TH, cv['gpi'][m] - ref[m], ':', color=C_GPI, lw=1.2, zorder=3)
            ymax = max(ymax, np.abs(cv['NNLO'][m] - ref[m] + cv['sig'][m]).max(), np.abs(cv['NNLO'][m] - ref[m] - cv['sig'][m]).max(),
                       np.abs(cv['F2'][m] - ref[m]).max(), np.abs(cv['gpi'][m] - ref[m]).max())
            if o == 'S2Z':                       # forward point (Sigma_2x vanishes at theta = 0)
                i0 = np.flatnonzero(sel(o, E, 0.0))[0]
                v, e = FWD[E]['value'] - ref[i0], FWD[E]['err']
                h_fwd = ax.errorbar([0.0], [v], [e], fmt='*', color='k', mfc='white', ms=6.5, lw=0.7, zorder=6)
                ymax = max(ymax, abs(v) + e)
            if relative:   # top-right corner: free in all panels (band and curves vanish at 180 deg); the top-left one is
                ax.text(0.96, 0.93, '%g MeV' % E, transform=ax.transAxes, va='top', ha='right')   # covered by the band at 140 MeV
            else:
                ax.text(0.04, 0.93, '%g MeV' % E, transform=ax.transAxes, va='top')
            ax.set_xlim(-6, 180); ax.set_xticks([0, 45, 90, 135, 180])
            if c == 0:
                ax.set_ylabel(('$\\Delta$' if relative else '') + sym)
            if r == 1:
                ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')
        yl = (np.ceil((ymax + 0.004) * 100) / 100) if relative else (np.ceil((ymax + 0.04) * 10) / 10)   # symmetric, common per row
        for c in range(3):
            axs[r, c].set_ylim(-yl, yl)
    handles = [(Patch(facecolor=C_BAND, alpha=0.6, lw=0), Line2D([], [], color=C_NNLO, lw=1.2)),
               Line2D([], [], color=C_F2, ls='--', lw=1.0), Line2D([], [], color=C_GPI, ls=':', lw=1.2), h_fwd]
    labels = ['NNLO B$\\chi$PT', 'NNLO fit', 'NNLO B$\\chi$PT, $\\gamma_\\pi+2$',
              'forward point']
    return handles, labels


for relative, files in ((True, [os.path.join(PAPER_FIGS, 'fig_pred_S2_diff.pdf'), os.path.join(FIGS, 'fig10b_pred_S2_diff.pdf')]),
                        (False, [os.path.join(FIGS, 'fig_pred_S2_abs_unused.pdf')])):
    fig, axs = plt.subplots(2, 3, figsize=(W2, 4.35), sharex=True, layout='constrained')
    handles, labels = draw(axs, relative)
    fig.legend(handles, labels, loc='outside lower center', ncol=4, columnspacing=1.6, handlelength=2.2)
    with plt.rc_context({'savefig.bbox': 'standard'}):
        for f in files:
            fig.savefig(f)
            print('wrote', f)
    plt.close(fig)

# ---------------------------------------------------------------- tables
rows = []
for E in ENERGIES:
    for th in TH_TAB:
        d = dict(E=E, theta=th)
        for o in OBS:
            m = sel(o, E, th); assert m.sum() == 1
            i = np.flatnonzero(m)[0]
            d[o] = cv['NNLO'][i]; d[o + '_sig'] = cv['sig'][i]; d[o + '_tr'] = cv['sig_tr'][i]; d[o + '_lec'] = cv['sig_lec'][i]
            d[o + '_dgpi'] = cv['shift_gpi'][i]
            d[o + '_fit'] = cv['F2'][i]; d[o + '_fit_sig'] = cv['sig_fit'][i]; d[o + '_fit_ab'] = cv['sig_ab'][i]
            d[o + '_fit_dgpi'] = cv['shift_F2gpi'][i]
        rows.append(d)
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(RES, 'pred_S2_table.csv'), index=False, float_format='%.4f')


def pm(v, e):
    return '$%+.3f \\pm %.3f$' % (v, e)


def write_table(name, head_tex, head_md, col):
    """col: (value, uncertainty, shift) column names of tab, without the observable prefix"""
    tex = head_tex + ['\\begin{tabular}{ccrrrr}', '\\hline',
                      '$E_\\gamma$ [MeV] & $\\theta_{\\rm lab}$ & $\\Si_{2x}$ & $\\De\\Si_{2x}(\\ga_\\pi+2)$ & $\\Si_{2z}$ & $\\De\\Si_{2z}(\\ga_\\pi+2)$ \\\\',
                      '\\hline']
    md = [head_md, '|---|---|---|---|---|---|']
    for _, r in tab.iterrows():
        v = [r['%s%s' % (o, c)] for o in OBS for c in col]
        tex.append('%g & $%g^\\circ$ & %s & $%+.3f$ & %s & $%+.3f$ \\\\' % (r.E, r.theta, pm(v[0], v[1]), v[2], pm(v[3], v[4]), v[5]))
        md.append('| %g | %g° | %+.3f ± %.3f | %+.3f | %+.3f ± %.3f | %+.3f |' % (r.E, r.theta, v[0], v[1], v[2], v[3], v[4], v[5]))
    tex += ['\\hline', '\\end{tabular}']
    open(os.path.join(RES, name + '.tex'), 'w').write('\n'.join(tex) + '\n')
    open(os.path.join(RES, name + '.md'), 'w').write('\n'.join(md) + '\n')
    return md


md = write_table('pred_S2_table',
                 ['% generated by scripts/s10_fig_pred_S2.py: NNLO prediction (P = PREF_NNLO, no free parameter), uncertainty =',
                  '% truncation (+) LECs (delta alpha_E1 = delta beta_M1 = 0.7) in quadrature; Delta = shift for gamma_pi + 2',
                  '% (gM1M1 + 1, gE1E1 - 1, gamma_0 fixed); lab angle; conventions of the A2 Sigma_2x, Sigma_2z data'],
                 '| E [MeV] | θ_lab | Σ2x (NNLO pred. ± tot.) | ΔΣ2x (γπ+2) | Σ2z (NNLO pred. ± tot.) | ΔΣ2z (γπ+2) |',
                 ('', '_sig', '_dgpi'))
print('\n'.join(md))
md_fit = write_table('pred_S2_table_fit',
                     ['% generated by scripts/s10_fig_pred_S2.py: NNLO fit (F2_d_all_thB), uncertainty = truncation (+) fit (quadrature);',
                      '% Delta = shift for gamma_pi + 2 (gM1M1 + 1, gE1E1 - 1); lab angle; conventions of the A2 Sigma_2x, Sigma_2z data'],
                     '| E [MeV] | θ_lab | Σ2x (NNLO fit ± tot.) | ΔΣ2x (γπ+2) | Σ2z (NNLO fit ± tot.) | ΔΣ2z (γπ+2) |',
                     ('_fit', '_fit_sig', '_fit_dgpi'))
tab.to_csv(os.path.join(RES, 'pred_S2_table_fit.csv'), index=False, float_format='%.4f')

# ---------------------------------------------------------------- sensitivity: d Sigma / d P_k at the NNLO prediction
# (central differences with step 0.5 units; the observable is a ratio of quadratic forms, so this is nearly exact)
sens = []
for E in ENERGIES:
    for th in TH_SENS:
        for o in OBS:
            i = np.flatnonzero(sel(o, E, th))[0]
            d = dict(obs=o, E=E, theta=th, value=float(cv['NNLO'][i]))
            for k, nm in enumerate(bchpt.PNAMES):
                dp = np.zeros(6); dp[k] = 0.5
                d['d_' + nm] = float((T.predict_array(PREF_NNLO + dp, [i]) - T.predict_array(PREF_NNLO - dp, [i]))[0])
            # gamma_pi at fixed gamma_0 (gM1M1 + 1/2, gE1E1 - 1/2 per unit) and gamma_0 at fixed gamma_pi (gM1M1 - 1/2, gE1E1 - 1/2)
            d['d_gpi_g0fixed'] = 0.5 * (d['d_gM1M1'] - d['d_gE1E1'])
            d['d_g0_gpifixed'] = -0.5 * (d['d_gM1M1'] + d['d_gE1E1'])
            d['gpi_units_per_0.01'] = 0.01 / abs(d['d_gpi_g0fixed'])     # precision 0.01 on Sigma -> units of gamma_pi
            sens.append(d)
sens = pd.DataFrame(sens)
sens.to_csv(os.path.join(RES, 'pred_S2_sensitivity.csv'), index=False, float_format='%.4f')
pd.set_option('display.width', 250)
print('\nsensitivity d Sigma / d P_k (per unit of 1e-4 fm^3 or fm^4) at the NNLO prediction:')
print(sens.round(4).to_string(index=False))

# ---------------------------------------------------------------- checks and the numbers of Sec. IV G
chk = {}
# (a) Sigma_2x must vanish at 0 and 180 deg; Sigma_2z(0) is the forward asymmetry (compare with the empirical GHP values)
chk['S2X_at_0_180'] = {('%g' % E): [float(cv['NNLO'][sel('S2X', E, 0.0)][0]), float(cv['NNLO'][sel('S2X', E, 180.0)][0])] for E in E_ALL}
chk['S2Z_at_180'] = {('%g' % E): float(cv['NNLO'][sel('S2Z', E, 180.0)][0]) for E in E_ALL}
chk['forward_S2Z'] = {('%g' % E): FWD[E] for E in E_ALL}
chk['flag_points'] = int(g.flag.sum())
chk['flag_angles'] = sorted(set(g.theta[g.flag > 0].tolist()))
# (b) numbers quoted in Sec. IV G: maximal differences from the prediction, band sizes, sensitivities
secG = {}
for o in OBS:
    for E in ENERGIES:
        m = sel(o, E)
        th_m = TH[np.argmax(np.abs(cv['F2'][m] - cv['NNLO'][m]))]
        th_g = TH[np.argmax(np.abs(cv['gpi'][m] - cv['NNLO'][m]))]
        mid = m & (g.theta.values >= 90.0) & (g.theta.values <= 120.0)
        secG['%s_%g' % (o, E)] = dict(
            max_abs_fit_minus_pred=float(np.abs(cv['F2'][m] - cv['NNLO'][m]).max()), at_theta_fit=float(th_m),
            max_abs_gpi_minus_pred=float(np.abs(cv['gpi'][m] - cv['NNLO'][m]).max()), at_theta_gpi=float(th_g),
            band_max=float(cv['sig'][m].max()), band_min=float(cv['sig'][m].min()),
            band_tr_max=float(cv['sig_tr'][m].max()), band_lec_max=float(cv['sig_lec'][m].max()),
            band_at_90=float(cv['sig'][sel(o, E, 90.0)][0]), band_tr_at_90=float(cv['sig_tr'][sel(o, E, 90.0)][0]),
            band_lec_at_90=float(cv['sig_lec'][sel(o, E, 90.0)][0]),
            gpi_shift_90_120=[float(cv['shift_gpi'][mid].min()), float(cv['shift_gpi'][mid].max())],
            fit_minus_pred_90_120=[float((cv['F2'] - cv['NNLO'])[mid].min()), float((cv['F2'] - cv['NNLO'])[mid].max())],
            c2=float(tr['c2'][sel(o, E, 90.0)][0]), c3=float(tr['c3'][sel(o, E, 90.0)][0]))
s2x = sens[sens.obs == 'S2X']; s2z = sens[(sens.obs == 'S2Z') & (sens.theta == 90.0)]
secG['S2X_dgpi_90_120'] = {('%g' % E): [float(s2x[s2x.E == E].d_gpi_g0fixed.min()), float(s2x[s2x.E == E].d_gpi_g0fixed.max())] for E in ENERGIES}
secG['S2X_gpi_units_per_0.01_90_120'] = {('%g' % E): [float(s2x[s2x.E == E]['gpi_units_per_0.01'].min()), float(s2x[s2x.E == E]['gpi_units_per_0.01'].max())] for E in ENERGIES}
secG['S2Z_90_dgpi'] = {('%g' % E): float(s2z[s2z.E == E].d_gpi_g0fixed.iloc[0]) for E in ENERGIES}
secG['S2Z_90_dg0'] = {('%g' % E): float(s2z[s2z.E == E].d_g0_gpifixed.iloc[0]) for E in ENERGIES}
secG['S2Z_90_dbeta'] = {('%g' % E): float(s2z[s2z.E == E].d_betaM1.iloc[0]) for E in ENERGIES}
secG['S2Z_90_dalpha'] = {('%g' % E): float(s2z[s2z.E == E].d_alphaE1.iloc[0]) for E in ENERGIES}
chk['secIVG'] = secG
# (c) comparison with work/delta_plots (fig5 curves: Sigma_2x,2z vs E at theta_lab = 90, 130 deg; 'main' = NNLO fit there)
cg = os.path.join(ROOT, 'work', 'delta_plots', 'results', 'curves_grid.csv')
if os.path.exists(cg):
    dp = pd.read_csv(cg)
    dp = dp[dp.id.str.startswith('g_f5_')]
    comp = []
    for o in OBS:
        for E in E_ALL:
            for th in (90.0, 130.0):
                a = dp[(dp.obs == o) & np.isclose(dp.E, E) & np.isclose(dp.theta, th)]
                if not len(a):
                    continue
                i = np.flatnonzero(sel(o, E, th))[0]
                a = a.iloc[0]
                comp.append(dict(obs=o, E=E, theta=th, F2=float(cv['F2'][i]), dp_main=float(a['main']), dF2=float(cv['F2'][i] - a['main']),
                                 NNLO=float(cv['NNLO'][i]), sig=float(cv['sig'][i]), sig_tr=float(cv['sig_tr'][i]), sig_lec=float(cv['sig_lec'][i]),
                                 dp_half=float(0.5 * (a['hi'] - a['lo'])), dp_trunc_half=float(0.5 * (a['trunc_hi'] - a['trunc_lo'])), dp_ab=float(a['ab'])))
    comp = pd.DataFrame(comp)
    comp.to_csv(os.path.join(RES, 'pred_S2_vs_delta_plots.csv'), index=False, float_format='%.5f')
    chk['delta_plots'] = dict(max_abs_dev_F2=float(comp.dF2.abs().max()), n=int(len(comp)), file=cg)
    print('\ncomparison with work/delta_plots fig5 curves (same executable, P_main = F2_d_all_thB):')
    print(comp.round(4).to_string(index=False))
json.dump(chk, open(os.path.join(RES, 'pred_S2_checks.json'), 'w'), indent=1)
print('\nchecks:', json.dumps(chk, indent=1))
print('PREF  =', PREF_NNLO, ' gamma0 = %.3f gamma_pi = %.3f' % (gamma0(PREF_NNLO), gammapi(PREF_NNLO)))
print('P_F2  =', P_F2, ' gamma0 = %.3f gamma_pi = %.3f' % (gamma0(P_F2), gammapi(P_F2)))
print('band at 90 deg (S2X/S2Z): trunc, LEC (d_alpha, d_beta), total | fit: ab, total:')
for o in OBS:
    for E in ENERGIES:
        i = np.flatnonzero(sel(o, E, 90.0))[0]
        print('  %s %5.0f MeV: %.4f %.4f (%+.4f %+.4f) %.4f | %.4f %.4f   c2=%.3f c3=%.3f' % (
            o, E, cv['sig_tr'][i], cv['sig_lec'][i], cv['d_alpha'][i], cv['d_beta'][i], cv['sig'][i], cv['sig_ab'][i], cv['sig_fit'][i], tr['c2'][i], tr['c3'][i]))
