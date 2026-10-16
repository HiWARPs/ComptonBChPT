"""Step 10: predictions of the double-polarization asymmetries Sigma_2x, Sigma_2z BELOW the pion-production
threshold (outlook of the paper: planned MAMI measurement with an active polarized target, Ref. Mornacchi:PhD).

Outputs
  paper/figs/fig_pred_S2.pdf              figure (double column): Sigma_2x (upper row), Sigma_2z (lower row) versus
                                          theta_lab at E_lab = 100, 120, 140 MeV
  work/fits/results/pred_S2_table.tex     LaTeX tabular: NNLO-fit values at theta_lab = 60, 90, 120, 150 deg with the
                                          total uncertainty (truncation (+) fit) and the shift for gamma_pi + 2
  work/fits/results/pred_S2_table.md      the same in markdown
  work/fits/results/pred_S2_curves.csv    all curves and uncertainties on the grid (reproducibility)
  work/fits/results/pred_S2_checks.json   the checks printed at the end (0/180 deg, forward point, delta_plots)

Curves (per panel)
  (1) NNLO fit        : P from fits_L.json 'F2_d_all_thB' (alpha_E1, beta_M1 fitted, spin polarizabilities at their
                        NNLO values); band = quadrature of the truncation uncertainty (scheme B, default settings of
                        theory.truncation, reference value 1 for asymmetries) and the (alpha, beta) fit uncertainty
                        propagated linearly with covP;
  (2) NNLO prediction : P = PREF_NNLO (no free parameter);
  (3) gamma's from the Delta-region fit: the complete parameter set of delta_region.json 'S_A2_LEGS_gpi'
                        (= Eq. spinfit of the paper: gE1E1 = -3.1, gM1M1 = 2.7, gE1M2 = -1.3, gM1E2 = 2.6, with
                        alpha_E1 = 11.01, beta_M1 = 3.10);
  (4) NNLO fit, gamma_pi + 2: (1) with gM1M1 + 1 and gE1E1 - 1 (gamma_0 unchanged).
Observable codes 'S2X', 'S2Z' of the driver, LAB frame (angle = lab scattering angle), i.e. the conventions of the
A2 Delta-region data (work/code_v2/driver/README.md Sect. 3):
  Sigma_2x = (sigma_R - sigma_L)/(sigma_R + sigma_L), target spin along +x (in the scattering plane, perpendicular to
             the beam, towards the scattered photon);  Sigma_2z: target spin along +z (beam direction);
  R/L = photon helicity +1/-1 (Lensky-McGovern-Pascalutsa 2015, Babusci et al. 1998; A2: Martel 2015, Paudyal 2020).
Not part of the `make` pipeline; run by hand:  .venv/bin/python scripts/s10_fig_pred_S2.py
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
ROOT = os.path.dirname(os.path.dirname(FITS))
PAPER_FIGS = os.path.join(ROOT, 'paper', 'figs')
os.makedirs(RES, exist_ok=True)

C_NNLO, C_F2, C_BAND, C_SPIN, C_GPI = '#444444', '#c0392b', '#9fb6dd', '#1f4e9c', '#2e8b57'
W2 = 510 / 72.27                     # double-column width [in] (revtex4-2 prd), included without scaling
TAG = 'pred_S2'

# ---------------------------------------------------------------- parameter sets
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
FD = json.load(open(os.path.join(RES, 'delta_region.json')))
fit = FL['F2_d_all_thB']
P_F2 = np.array(fit['P']); COV_F2 = np.array(fit['covP'])
assert fit['free'] == ['alphaE1', 'betaM1'], fit['free']
P_SPIN = np.array(FD['S_A2_LEGS_gpi']['P'])
# Eq. (spinfit) of the paper: gE1E1 = -3.1(6), gM1M1 = 2.7(4), gE1M2 = -1.3(8), gM1E2 = 2.6(4); alpha = 11.01, beta = 3.10
assert np.allclose(np.round(P_SPIN, 2), [11.01, 3.10, -3.09, 2.69, -1.26, 2.60], atol=0.006), P_SPIN
D_GPI = np.array([0.0, 0.0, -1.0, +1.0, 0.0, 0.0])       # gE1E1 - 1, gM1M1 + 1: gamma_pi + 2, gamma_0 unchanged
P_GPI = P_F2 + D_GPI
assert abs(gammapi(P_GPI) - gammapi(P_F2) - 2.0) < 1e-12 and abs(gamma0(P_GPI) - gamma0(P_F2)) < 1e-12

# ---------------------------------------------------------------- grid
ENERGIES = [100.0, 120.0, 140.0]              # panels
E_ALL = ENERGIES + [130.0]                    # 130 MeV only for the comparison with work/delta_plots (fig5)
TH = np.arange(0.0, 180.01, 5.0)              # 37 lab angles
OBS = ['S2X', 'S2Z']
TH_TAB = [60.0, 90.0, 120.0, 150.0]


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
cv = {'NNLO': T.value_file.copy(), 'F2': T.predict_array(P_F2), 'spin': T.predict_array(P_SPIN), 'gpi': T.predict_array(P_GPI)}
assert np.allclose(cv['NNLO'], T.predict_array(PREF_NNLO), atol=1e-9)
# truncation uncertainty: orders y0 (Born + pi0 + contact), y2 (NLO + contact), y3 (NNLO), all at the NNLO static
# polarizabilities, exactly as order_values() does for the database rows; reference value 1 for the asymmetries
_, tr = truncation_for(g, g, tag=TAG)
assert np.all(tr['yref'] == 1.0)
cv['sig_tr'] = tr['sigma']
# fit uncertainty of (alpha_E1, beta_M1): linear propagation with the fit covariance (central differences)
eps = 0.01
grad = []
for k in range(2):
    dp = np.zeros(6); dp[k] = eps
    grad.append((T.predict_array(P_F2 + dp) - T.predict_array(P_F2 - dp)) / (2 * eps))
grad = np.array(grad)
cv['sig_ab'] = np.sqrt(np.einsum('ip,ij,jp->p', grad, COV_F2[:2, :2], grad))
cv['sig'] = np.sqrt(cv['sig_tr'] ** 2 + cv['sig_ab'] ** 2)
cv['shift_gpi'] = cv['gpi'] - cv['F2']
for k, v in cv.items():
    g[k] = v
g['y0'], g['y2'], g['y3'] = tr['y0'], tr['y2'], tr['y3']
g.to_csv(os.path.join(RES, 'pred_S2_curves.csv'), index=False, float_format='%.6g')


def sel(o, E, th=None):
    m = (g.dobs == o).values & np.isclose(g.E.values, E)
    if th is not None:
        m &= np.isclose(g.theta.values, th)
    return m


# ---------------------------------------------------------------- figure
fig, axs = plt.subplots(2, 3, figsize=(W2, 4.35), sharex=True, layout='constrained')
for r, (o, sym) in enumerate([('S2X', '$\\Sigma_{2x}$'), ('S2Z', '$\\Sigma_{2z}$')]):
    ymax = 0.0
    for c, E in enumerate(ENERGIES):
        ax = axs[r, c]; m = sel(o, E)
        ax.axhline(0, color='0.6', lw=0.4, zorder=0)
        ax.fill_between(TH, cv['F2'][m] - cv['sig'][m], cv['F2'][m] + cv['sig'][m], color=C_BAND, alpha=0.6, lw=0, zorder=1)
        ax.plot(TH, cv['NNLO'][m], '--', color=C_NNLO, lw=1.0, zorder=3)
        ax.plot(TH, cv['spin'][m], '-.', color=C_SPIN, lw=1.0, zorder=3)
        ax.plot(TH, cv['gpi'][m], ':', color=C_GPI, lw=1.2, zorder=3)
        ax.plot(TH, cv['F2'][m], '-', color=C_F2, lw=1.2, zorder=4)
        ymax = max(ymax, np.abs(cv['F2'][m] + cv['sig'][m]).max(), np.abs(cv['F2'][m] - cv['sig'][m]).max(),
                   np.abs(cv['NNLO'][m]).max(), np.abs(cv['spin'][m]).max(), np.abs(cv['gpi'][m]).max())
        ax.text(0.04, 0.93, '%g MeV' % E, transform=ax.transAxes, va='top')
        ax.set_xlim(0, 180); ax.set_xticks([0, 45, 90, 135, 180])
        if c == 0:
            ax.set_ylabel(sym)
        if r == 1:
            ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')
    yl = np.ceil((ymax + 0.04) * 10) / 10         # symmetric, common per row
    for c in range(3):
        axs[r, c].set_ylim(-yl, yl)
handles = [(Patch(facecolor=C_BAND, alpha=0.6, lw=0), Line2D([], [], color=C_F2, lw=1.2)),
           Line2D([], [], color=C_NNLO, ls='--', lw=1.0), Line2D([], [], color=C_SPIN, ls='-.', lw=1.0),
           Line2D([], [], color=C_GPI, ls=':', lw=1.2)]
labels = ['NNLO fit $\\pm$ (truncation $\\oplus$ fit)', 'NNLO prediction', "$\\gamma$'s from the $\\Delta$-region fit",
          'NNLO fit, $\\gamma_\\pi+2$']
fig.legend(handles, labels, loc='outside lower center', ncol=4, columnspacing=1.6, handlelength=2.2)
with plt.rc_context({'savefig.bbox': 'standard'}):
    fig.savefig(os.path.join(PAPER_FIGS, 'fig_pred_S2.pdf'))
    fig.savefig(os.path.join(FIGS, 'fig10_pred_S2.pdf'))
plt.close(fig)
print('wrote', os.path.join(PAPER_FIGS, 'fig_pred_S2.pdf'))

# companion figure: the same curves relative to the NNLO fit (the asymmetries below threshold are Born-dominated and
# the four curves of fig_pred_S2 nearly coincide; here the sensitivity to the spin polarizabilities is visible)
fig, axs = plt.subplots(2, 3, figsize=(W2, 4.35), sharex=True, layout='constrained')
for r, (o, sym) in enumerate([('S2X', '$\\Sigma_{2x}$'), ('S2Z', '$\\Sigma_{2z}$')]):
    ymax = 0.0
    for c, E in enumerate(ENERGIES):
        ax = axs[r, c]; m = sel(o, E)
        ax.axhline(0, color=C_F2, lw=1.2, zorder=4)
        ax.fill_between(TH, -cv['sig'][m], cv['sig'][m], color=C_BAND, alpha=0.6, lw=0, zorder=1)
        ax.plot(TH, cv['NNLO'][m] - cv['F2'][m], '--', color=C_NNLO, lw=1.0, zorder=3)
        ax.plot(TH, cv['spin'][m] - cv['F2'][m], '-.', color=C_SPIN, lw=1.0, zorder=3)
        ax.plot(TH, cv['gpi'][m] - cv['F2'][m], ':', color=C_GPI, lw=1.2, zorder=3)
        ymax = max(ymax, cv['sig'][m].max(), np.abs(cv['NNLO'][m] - cv['F2'][m]).max(), np.abs(cv['spin'][m] - cv['F2'][m]).max(),
                   np.abs(cv['gpi'][m] - cv['F2'][m]).max())
        ax.text(0.04, 0.93, '%g MeV' % E, transform=ax.transAxes, va='top')
        ax.set_xlim(0, 180); ax.set_xticks([0, 45, 90, 135, 180])
        if c == 0:
            ax.set_ylabel('$\\Delta$' + sym)
        if r == 1:
            ax.set_xlabel('$\\theta_{\\rm lab}$ [deg]')
    yl = np.ceil((ymax + 0.004) * 100) / 100
    for c in range(3):
        axs[r, c].set_ylim(-yl, yl)
labels_d = ['NNLO fit $\\pm$ (truncation $\\oplus$ fit)', 'NNLO prediction', "$\\gamma$'s from the $\\Delta$-region fit",
            'NNLO fit, $\\gamma_\\pi+2$']
fig.legend(handles, labels_d, loc='outside lower center', ncol=4, columnspacing=1.6, handlelength=2.2)
with plt.rc_context({'savefig.bbox': 'standard'}):
    fig.savefig(os.path.join(PAPER_FIGS, 'fig_pred_S2_diff.pdf'))
    fig.savefig(os.path.join(FIGS, 'fig10b_pred_S2_diff.pdf'))
plt.close(fig)
print('wrote', os.path.join(PAPER_FIGS, 'fig_pred_S2_diff.pdf'))

# ---------------------------------------------------------------- table
rows = []
for E in ENERGIES:
    for th in TH_TAB:
        d = dict(E=E, theta=th)
        for o in OBS:
            m = sel(o, E, th); assert m.sum() == 1
            i = np.flatnonzero(m)[0]
            d[o] = cv['F2'][i]; d[o + '_sig'] = cv['sig'][i]; d[o + '_tr'] = cv['sig_tr'][i]; d[o + '_ab'] = cv['sig_ab'][i]
            d[o + '_dgpi'] = cv['shift_gpi'][i]; d[o + '_NNLO'] = cv['NNLO'][i]; d[o + '_spin'] = cv['spin'][i]
        rows.append(d)
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(RES, 'pred_S2_table.csv'), index=False, float_format='%.4f')


def pm(v, e):
    return '$%+.3f \\pm %.3f$' % (v, e)


tex = ['% generated by scripts/s10_fig_pred_S2.py: NNLO fit (F2_d_all_thB), uncertainty = truncation (+) fit (quadrature);',
       '% Delta = shift for gamma_pi + 2 (gM1M1 + 1, gE1E1 - 1); lab angle; conventions of the A2 Sigma_2x, Sigma_2z data',
       '\\begin{tabular}{ccrrrr}', '\\hline',
       '$E_\\gamma$ [MeV] & $\\theta_{\\rm lab}$ & $\\Si_{2x}$ & $\\De\\Si_{2x}(\\ga_\\pi+2)$ & $\\Si_{2z}$ & $\\De\\Si_{2z}(\\ga_\\pi+2)$ \\\\',
       '\\hline']
md = ['| E [MeV] | θ_lab | Σ2x (NNLO fit ± tot.) | ΔΣ2x (γπ+2) | Σ2z (NNLO fit ± tot.) | ΔΣ2z (γπ+2) |',
      '|---|---|---|---|---|---|']
for _, r in tab.iterrows():
    tex.append('%g & $%g^\\circ$ & %s & $%+.3f$ & %s & $%+.3f$ \\\\' % (
        r.E, r.theta, pm(r.S2X, r.S2X_sig), r.S2X_dgpi, pm(r.S2Z, r.S2Z_sig), r.S2Z_dgpi))
    md.append('| %g | %g° | %+.3f ± %.3f | %+.3f | %+.3f ± %.3f | %+.3f |' % (
        r.E, r.theta, r.S2X, r.S2X_sig, r.S2X_dgpi, r.S2Z, r.S2Z_sig, r.S2Z_dgpi))
tex += ['\\hline', '\\end{tabular}']
open(os.path.join(RES, 'pred_S2_table.tex'), 'w').write('\n'.join(tex) + '\n')
open(os.path.join(RES, 'pred_S2_table.md'), 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))

# ---------------------------------------------------------------- sensitivity: d Sigma / d P_k at the NNLO fit
# (central differences with step 0.5 units; the observable is a ratio of quadratic forms, so this is nearly exact)
sens = []
for E in ENERGIES:
    for th in (90.0, 120.0):
        for o in OBS:
            i = np.flatnonzero(sel(o, E, th))[0]
            d = dict(obs=o, E=E, theta=th, value=float(cv['F2'][i]))
            for k, nm in enumerate(bchpt.PNAMES):
                dp = np.zeros(6); dp[k] = 0.5
                d['d_' + nm] = float((T.predict_array(P_F2 + dp, [i]) - T.predict_array(P_F2 - dp, [i]))[0])
            # gamma_pi at fixed gamma_0 (gM1M1 + 1/2, gE1E1 - 1/2 per unit) and gamma_0 at fixed gamma_pi (gM1M1 - 1/2, gE1E1 - 1/2)
            d['d_gpi_g0fixed'] = 0.5 * (d['d_gM1M1'] - d['d_gE1E1'])
            d['d_g0_gpifixed'] = -0.5 * (d['d_gM1M1'] + d['d_gE1E1'])
            sens.append(d)
sens = pd.DataFrame(sens)
sens.to_csv(os.path.join(RES, 'pred_S2_sensitivity.csv'), index=False, float_format='%.4f')
print('\nsensitivity d Sigma / d P_k (per unit of 1e-4 fm^3 or fm^4) at the NNLO fit:')
print(sens.round(4).to_string(index=False))

# ---------------------------------------------------------------- checks
chk = {}
# (a) Sigma_2x must vanish at 0 and 180 deg; Sigma_2z(0) is the forward asymmetry (compare with the empirical GHP values)
chk['S2X_at_0_180'] = {('%g' % E): [float(cv['F2'][sel('S2X', E, 0.0)][0]), float(cv['F2'][sel('S2X', E, 180.0)][0])] for E in E_ALL}
chk['S2Z_at_180'] = {('%g' % E): float(cv['F2'][sel('S2Z', E, 180.0)][0]) for E in E_ALL}
df = load_db()
fw = df[df.dataset == 'GHP_fwd_S2Z'].sort_values('E')
chk['S2Z_at_0'] = {}
for E in E_ALL:
    i = np.flatnonzero(sel('S2Z', E, 0.0))[0]
    emp = fw[np.isclose(fw.E, E)]
    chk['S2Z_at_0']['%g' % E] = dict(F2=float(cv['F2'][i]), sig=float(cv['sig'][i]), NNLO=float(cv['NNLO'][i]),
                                     GHP=(float(emp.value.iloc[0]) if len(emp) else None), GHP_err=(float(emp.err_pt.iloc[0]) if len(emp) else None))
chk['flag_points'] = int(g.flag.sum())
chk['flag_angles'] = sorted(set(g.theta[g.flag > 0].tolist()))
# (b) comparison with work/delta_plots (fig5 curves: Sigma_2x,2z vs E at theta_lab = 90, 130 deg)
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
                                 spin=float(cv['spin'][i]), dp_spin=float(a['spin']), dspin=float(cv['spin'][i] - a['spin']),
                                 sig=float(cv['sig'][i]), sig_tr=float(cv['sig_tr'][i]), sig_ab=float(cv['sig_ab'][i]),
                                 dp_half=float(0.5 * (a['hi'] - a['lo'])), dp_trunc_half=float(0.5 * (a['trunc_hi'] - a['trunc_lo'])), dp_ab=float(a['ab'])))
    comp = pd.DataFrame(comp)
    comp.to_csv(os.path.join(RES, 'pred_S2_vs_delta_plots.csv'), index=False, float_format='%.5f')
    chk['delta_plots'] = dict(max_abs_dev_F2=float(comp.dF2.abs().max()), max_abs_dev_spin=float(comp.dspin.abs().max()),
                              n=int(len(comp)), file=cg)
    pd.set_option('display.width', 200)
    print('\ncomparison with work/delta_plots fig5 curves (same executable, P_main = F2_d_all_thB, spin fit S_A2_LEGS_gpi):')
    print(comp.round(4).to_string(index=False))
json.dump(chk, open(os.path.join(RES, 'pred_S2_checks.json'), 'w'), indent=1)
print('\nchecks:', json.dumps(chk, indent=1))
print('P_F2  =', P_F2, ' gamma0 = %.3f gamma_pi = %.3f' % (gamma0(P_F2), gammapi(P_F2)))
print('P_SPIN=', P_SPIN, ' gamma0 = %.3f gamma_pi = %.3f' % (gamma0(P_SPIN), gammapi(P_SPIN)))
print('PREF  =', PREF_NNLO)
print('band at 90 deg (S2X/S2Z): trunc, ab, total:')
for o in OBS:
    for E in ENERGIES:
        i = np.flatnonzero(sel(o, E, 90.0))[0]
        print('  %s %5.0f MeV: %.4f %.4f %.4f   c2=%.3f c3=%.3f' % (o, E, cv['sig_tr'][i], cv['sig_ab'][i], cv['sig'][i], tr['c2'][i], tr['c3'][i]))
