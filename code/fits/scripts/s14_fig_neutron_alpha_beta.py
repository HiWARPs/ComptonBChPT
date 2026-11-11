"""Figure: neutron scalar polarizabilities in the alpha-beta plane (paper Fig. "neutronab", Sec. IV D since 2026-11-10).

Shown:
  * neutron Baldin sum rule 15.2, uncertainty 0.3-0.5 (Levchuk & L'vov, NPA 674 (2000) 449, Eq. (3.11b) and text;
    the deuteron analyses use +-0.5), gray band of +-0.5;
  * NNLO BChPT prediction for the neutron (paper Table "neutron": 13.71, 4.63; re-examined uncertainties 3.16, 1.17);
  * "BChPT, proton LECs": the main proton fit plus the NNLO neutron-minus-proton difference (only the piN loops differ
    at this order, Delta graphs isospin symmetric), i.e. the O(p^4) counterterms of the proton fit taken as isoscalar.
    Covariance = main-fit covariance (combined_fit.json 'combined') + its systematic (Eq. mainresult: 0.45, 0.37,
    anticorrelated, as for alpha-beta in Tab. summary: the forward SRs fix the sum in every variant) + the truncation
    of the n-p difference by the rule of Sec. II C applied to the sum and the difference directly, delta^2 |n - p|
    (delta = 0.48; referee 2026-11-09 E2; before: syst uncorrelated, rule on alpha and beta separately);
  * empirical: PDG 2026; Kossert et al. 2003 (quasi-free, EPJA 16, 259: 12.5 +-1.8 +1.1-0.6 +-1.1); Myers et al. 2014
    (COMPTON@MAX-lab, world elastic deuteron data, as quoted in Godagama et al. Eq. (3): 11.55 +-1.25 +-0.2 +-0.8);
    Godagama et al. 2026 (HIgS, arXiv:2611.06609, Eq. (9): 12.25 +-1.1 +-0.3 +-0.8).  Errors added in quadrature;
    the Baldin-constrained values are anticorrelated (drawn, as in fig4_alpha_beta, as crosses).
  * Alarcon et al. 2020 (PRD 102, 014006, arXiv:2005.09518; texkey Alarcon:2020wjg), Eqs. (16b), (17): the O(p^4)
    counterterm of beta fitted to the neutron Baldin SR 15.2(5) with alpha at its NNLO value -> beta_n = 1.5(5)
    (no theory error), alpha_n = 15.2 - 1.5 = 13.7; the error is in beta only (added 2026-11-10, author request).
Check printed: the two routes to the "proton LECs" point (proton fit + (n-p) vs neutron NNLO + delta alpha, delta beta
of Eq. LECvalues) agree to 0.03 (they differ only by the g_E dependence of the Delta pole in alpha).
Writes work/fits/figs/fig_neutron_alpha_beta.pdf and paper/figs/fig_neutron_alpha_beta.pdf (single column).
Run (from work/fits):  ../../.venv/bin/python scripts/s14_fig_neutron_alpha_beta.py
"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from common import *  # noqa

plt.style.use(STYLE)
PAPER_FIGS = os.path.join(os.path.dirname(os.path.dirname(FITS)), 'paper', 'figs')
W1 = 246 / 72.27
C_NNLO, C_FIT = '#444444', '#1f4e9c'                    # as in s06 (prediction black, fit blue)
DELTA = (1232.0 - 938.272) / 775.26  # delta = Delta/m_rho = 0.379 (author, 2026-11-11; LMP15 used m_pi/Delta = 0.48)

# NNLO BChPT, paper Table "neutron" (rows "NNLO BChPT, neutron", "uncertainty, this work", "Neutron - proton")
N_NNLO, N_ERR = np.array([13.71, 4.63]), np.array([1.97, 0.82])       # uncertainties: delta = Delta/m_rho (2026-11-11)
N_MINUS_P = np.array([2.30, 0.90])
BALDIN_N = (15.2, 0.5)

MF = json.load(open(os.path.join(RES, 'combined_fit.json')))['combined']
P_p = np.array(MF['P'][:2])
SYST = np.array([0.4501, -0.369])                       # anticorrelated systematic (Eq. mainresult)
A = np.array([[0.5, 0.5], [0.5, -0.5]])                 # (sum, difference) -> (alpha, beta)
TRUNC = A @ np.diag((DELTA ** 2 * np.abs([N_MINUS_P.sum(), N_MINUS_P[0] - N_MINUS_P[1]])) ** 2) @ A.T
cov = np.array(MF['covP'])[:2, :2] + np.outer(SYST, SYST) + TRUNC
P_n = P_p + N_MINUS_P
alt = N_NNLO + np.array([-0.649, -0.470])               # Eq. LECvalues (delta = Delta/m_rho refit, 2026-11-11)
s_sum, s_dif = np.sqrt(cov.sum()), np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1])
print('proton LECs: alpha_n = %.2f(%.2f), beta_n = %.2f(%.2f), rho = %.2f, sum = %.2f(%.2f) [%.2f sigma above BSR],'
      ' diff = %.2f(%.2f); via NNLO_n + LECs: %.2f, %.2f'
      % (P_n[0], np.sqrt(cov[0, 0]), P_n[1], np.sqrt(cov[1, 1]), cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]),
         P_n.sum(), s_sum, (P_n.sum() - BALDIN_N[0]) / np.hypot(s_sum, BALDIN_N[1]), P_n[0] - P_n[1], s_dif, *alt))
assert np.all(np.abs(P_n - alt) < 0.05), 'the two routes to the proton-LEC point disagree'

# Baldin + chiral alpha-beta (2026-11-11): the isovector O(p^4) correction lowers the sum to the neutron Baldin SR and
# leaves the difference, as the proton fit does with the isoscalar one; sum and difference are independent inputs.
D_n, var_D = P_n[0] - P_n[1], s_dif ** 2
P_B = A @ np.array([BALDIN_N[0], D_n])
cov_B = A @ np.diag([BALDIN_N[1] ** 2, var_D]) @ A.T
x_iv = (alt.sum() - BALDIN_N[0]) / 2                    # check: NNLO_n + LECs, isovector correction split equally
assert np.allclose(P_B.sum(), BALDIN_N[0]) and np.isclose(P_B[0] - P_B[1], D_n)
assert np.all(np.abs(P_B - (alt - x_iv)) < 0.05), 'the two routes to the Baldin point disagree'
dn_mp = abs(N_MINUS_P[0] - N_MINUS_P[1])                # delta instead of delta^2 for the isovector part of alpha-beta
var_D_d = var_D - (DELTA ** 2 * dn_mp) ** 2 + (DELTA * dn_mp) ** 2
print('Baldin + alpha-beta: alpha_n = %.2f(%.2f), beta_n = %.2f(%.2f), rho = %.2f; isovector shift of the sum %.2f '
      '(each of alpha, beta %.2f); via NNLO_n + LECs: %.2f, %.2f; with delta: sigma = %.2f'
      % (P_B[0], np.sqrt(cov_B[0, 0]), P_B[1], np.sqrt(cov_B[1, 1]), cov_B[0, 1] / cov_B[0, 0], BALDIN_N[0] - P_n.sum(),
         -x_iv, *(alt - x_iv), np.sqrt(BALDIN_N[1] ** 2 + var_D_d) / 2))


def ellipse(ax, P, c, **kw):                            # 68% joint CL, as s06
    w, v = np.linalg.eigh(c)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    k = np.sqrt(2.30)
    e = Ellipse(P, 2 * k * np.sqrt(w[1]), 2 * k * np.sqrt(w[0]), angle=ang, fill=False, **kw)
    ax.add_patch(e)
    return e


fig, ax = plt.subplots(figsize=(W1, 3.6), layout='constrained')
a = np.linspace(5, 20, 50)
h_bal = ax.fill_between(a, BALDIN_N[0] - BALDIN_N[1] - a, BALDIN_N[0] + BALDIN_N[1] - a, color='#dddddd', lw=0,
                        label='Baldin SR, 15.2(5)')
h_nnlo = ax.errorbar(*N_NNLO, N_ERR[1], N_ERR[0], fmt='P', color='k', ms=5, lw=0.7, label='NNLO B$\\chi$PT')
h_fit = ellipse(ax, P_n, cov, ec=C_FIT, ls='-', lw=2.0, label='B$\\chi$PT, proton LECs')
ax.plot(*P_n, 'o', color=C_FIT, ms=2)
h_bfit = ellipse(ax, P_B, cov_B, ec=C_FIT, ls='--', lw=1.4, label='B$\\chi$PT, LECs + Baldin SR')
ax.plot(*P_B, 'o', color=C_FIT, ms=2)
q = lambda *x: float(np.sqrt(np.sum(np.square(x))))
emp = [  # (alpha, (lo, hi)), (beta, (lo, hi)), colour, marker, label
    ((11.8, (1.1, 1.1)), (3.7, (1.2, 1.2)), '#6a3d9a', 's', 'PDG 2026'),
    ((12.5, (q(1.8, 0.6, 1.1), q(1.8, 1.1, 1.1))), (2.7, (q(1.8, 1.1, 1.1), q(1.8, 0.6, 1.1))), '#e08e0b', '>',
     'Kossert 03 (quasi-free)'),
    ((11.55, (q(1.25, 0.2, 0.8),) * 2), (3.65, (q(1.25, 0.2, 0.8),) * 2), '#6b8e23', 'v', 'Myers 14'),
    ((12.25, (q(1.1, 0.3, 0.8),) * 2), (2.95, (q(1.1, 0.3, 0.8),) * 2), '#d62e8a', '^', 'HI$\\gamma$S 26'),
]
hs = []
for (av, (al, ah)), (bv, (bl, bh)), c, m, lab in emp:
    hs.append(ax.errorbar(av, bv, [[bl], [bh]], [[al], [ah]], fmt=m, color=c, mfc=c, ms=4.5, lw=0.7, label=lab))
AHLP = (15.2 - 1.5, 1.5, 0.5)                           # Alarcon:2020wjg Eqs. (16b), (17): beta fitted to the Baldin SR
assert abs(AHLP[0] - N_NNLO[0]) < 0.05, 'Alarcon:2020wjg keeps alpha_n at its NNLO value'
h_ahlp = ax.errorbar(AHLP[0], AHLP[1], AHLP[2], fmt='D', color=C_FIT, mfc='white', ms=4, lw=0.7,
                     label='B$\\chi$PT, Baldin fit [Alarc\\\'on 20]')
ax.set_xlim(8.5, 17.5); ax.set_ylim(0, 6.5)
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1)); ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
ax.set_xlabel('$\\alpha_{E1}^{(n)}$ [$10^{-4}$ fm$^3$]'); ax.set_ylabel('$\\beta_{M1}^{(n)}$ [$10^{-4}$ fm$^3$]')
handles = [h_bal, h_nnlo, h_fit, h_bfit, h_ahlp] + hs
fig.legend(handles, [h.get_label() for h in handles], loc='outside lower center', ncol=2, handlelength=2.0,
           handletextpad=0.5, columnspacing=1.0)
for d in (FIGS, PAPER_FIGS):
    with plt.rc_context({'savefig.bbox': 'standard'}):
        fig.savefig(os.path.join(d, 'fig_neutron_alpha_beta.pdf'))
print('wrote fig_neutron_alpha_beta.pdf to', FIGS, 'and', PAPER_FIGS)
