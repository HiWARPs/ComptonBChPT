"""NNLO BChPT prediction for the static dipole polarisabilities of the NEUTRON, with the uncertainty re-examined as for
the proton (scenarios/spin_pol_uncertainties.py; author request 2026-10-20: "include a subsection on predictions for
neutron polarizabilities, scalar and spin with the revised uncertainty estimate analogously to the proton").
Mechanisms (Lensky, McGovern, Pascalutsa, EPJC 75 (2015) 604 = LMP15, Sect. 4.3 and App. C):
  piN loops   O(p^3): the neutron expressions of LMP15 App. C.1(b) (alpha, beta also in the author's notebooks
              ~/MySci/Math/ComptonScattering/{EM_Polarizabilities_Results_Neutron, Neutron_ScalarPolarizabilities_NNLO_BChPT}.nb,
              identical; the 2015 notebook evaluated with wolframscript agrees to 1e-9, check (iv)).  Only the pi^- p loop contributes (pi^0 n: no charges at this order), so the physical-mass
              evaluation is the formula with the charged-pion mass; nucleon mass M_n (M_p in the loop changes the values by
              < 0.01, printed).
  piDelta loops O(p^{7/2}) and Delta pole: equal to the proton ones at this order (LMP15 Sect. 4.3); evaluated with
              code_v2/driver/analytic_pols.py at M_N = M_n (at M_p they coincide with the driver's proton values, checked).
Uncertainty (as for the proton): LMP15 rule Max[d^2|piN|, d|Delta graphs|, d^2|total|], d = 0.48 (reproduces LMP15's
neutron uncertainties, checked); re-examined: the Delta pole carries the uncertainty of the couplings (g_M = 2.97 +- 0.10,
R_EM = default -+ 0.5 %, exact relation lib/rem.py), in quadrature with Max[d^2|piN|, d|piDelta|, d^2|total|].
Writes results/neutron_pols.json."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(HERE, '..', 'lib'), os.path.join(HERE, '..', '..', 'code_v2', 'driver')]
import warnings
import numpy as np
import analytic_pols as AP
import rem as RM
from common import save_json, RES

D = 0.48
MP, MNEU, MPI = 0.93827, 0.939565, 0.13957
LAB = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2', 'gamma0', 'gammapi']


def piN_neutron(gA=1.27, F=0.09221, MN=MNEU, mpi=MPI, aem=1/137.036, hbarc=0.19733):
    """O(p^3) piN-loop contribution to the neutron (LMP15 Eqs. (C9)-(C14)); alpha, beta [1e-4 fm^3], gammas [1e-4 fm^4]"""
    e2 = 4 * np.pi * aem; mu = mpi / MN; s = np.sqrt(4 - mu**2); A = np.arccos(mu**2 / 2 - 1); L = np.log(mu**2)
    ca = e2 * gA**2 / (192 * np.pi**3 * F**2 * MN)
    a = ca / (mu * s**5) * (mu * s * (-mu**2 + 3 * (mu**2 - 4)**2 * L + 16) + (-3 * mu**6 + 30 * mu**4 - 98 * mu**2 + 80) * A)
    b = -ca / (mu * s**3) * (mu * s * (3 * (mu**2 - 4) * L - 11) - (3 * mu**4 - 18 * mu**2 + 2) * A)
    k = e2 * gA**2 / (np.pi**3 * F**2 * MN**2 * mu**2 * s**5)
    g1 = k / 96 * (s * (2 * (2 * mu**4 - 7 * mu**2 + 8) - mu**2 * (mu**2 - 4)**2 * L) + mu * (mu**6 - 10 * mu**4 + 24 * mu**2 - 12) * A)
    g2 = k / 192 * (s * (6 * mu**4 - 52 * mu**2 - 3 * (mu**2 - 4)**2 * mu**2 * L + 16) + mu * (3 * mu**6 - 30 * mu**4 + 94 * mu**2 - 40) * A)
    g3 = -k / 384 * (s * (mu**2 * (mu**2 - 4)**2 * L + 2 * (mu**4 + 10 * mu**2 - 8)) - mu * (mu**6 - 10 * mu**4 + 46 * mu**2 - 40) * A)
    g4 = k / 384 * (s * (-6 * mu**4 + 76 * mu**2 + 5 * (mu**2 - 4)**2 * mu**2 * L - 16) - mu * (5 * mu**6 - 50 * mu**4 + 166 * mu**2 - 88) * A)
    f3, f4 = hbarc**3 * 1e4, hbarc**4 * 1e4
    return np.array([a * f3, b * f3, (-g1 - g3) * f4, g4 * f4, g3 * f4, (g2 + g4) * f4])


def ext(v):
    return np.r_[v, -(v[2] + v[3] + v[4] + v[5]), -v[2] + v[3] - v[4] + v[5]]


def lmp_rule(piN, delta_graphs, tot):
    return np.maximum.reduce([D**2 * np.abs(piN), D * np.abs(delta_graphs), D**2 * np.abs(tot)])


with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    # ---- checks ----
    # (i) LMP15 Tables 5, 6 (O(p^3) piN rows, their parameter set = the current constants with isospin-symmetric masses)
    LMP_piN = np.array([9.4, -1.1, -4.7, -0.2, 0.6, 1.3, 3.0, 5.3])
    assert np.all(np.abs(ext(piN_neutron(MN=MP)) - LMP_piN) < 0.051), ext(piN_neutron(MN=MP))
    # (ii) heavy-baryon limit: neutron and proton piN loops coincide for M_N -> infinity (at fixed m_pi)
    assert np.allclose(piN_neutron(MN=1e4), AP.pols(MN=1e4)['piN loops'], rtol=2e-3)
    # (iv) the author's notebook Neutron_ScalarPolarizabilities_NNLO_BChPT.nb (2015), evaluated with wolframscript
    #      (alpha_em = 1/137.04, M_N = 938.27, m_pi = 139.57): piN, Delta pole and piDelta alpha, beta agree to 1e-9
    assert np.allclose(piN_neutron(MN=MP, aem=1/137.04)[:2], [9.38065790174812, -1.1012159293598598], rtol=1e-9)
    assert np.allclose(AP.pols(aem=1/137.04)['Delta pole'][:2], [-0.10970355263633684, 7.149884931959332], rtol=1e-9)
    assert np.allclose(AP.pols(aem=1/137.04)['piDelta loops'][:2], [4.4252053617252844, -1.4465975671009133], rtol=1e-7)
    # (iii) LMP15 neutron totals and uncertainties (their Tables 5, 6) from their rounded mechanism rows
    LMP_tot = np.array([13.7, 4.6, -4.7, 2.9, 0.2, 1.6, 0.03, 9.0]); LMP_err = np.array([3.1, 2.7, 1.1, 1.5, 0.2, 0.4, 1.4, 2.0])
    lp = dict(piN=LMP_piN, piD=np.array([4.4, -1.4, 0.4, -0.2, 0.1, -0.2, -0.1, -0.9]), pole=np.array([-0.1, 7.1, -0.4, 3.3, -0.4, 0.4, -2.8, 4.5]))
    assert np.all(np.abs(lp['piN'] + lp['piD'] + lp['pole'] - LMP_tot) < 0.11)
    rule_lmp = lmp_rule(lp['piN'], lp['piD'] + lp['pole'], LMP_tot)
    print('LMP15 rule on the LMP15 neutron rows:', np.round(rule_lmp, 2), ' published:', LMP_err)
    assert np.all(np.abs(rule_lmp - LMP_err)[[0, 1, 2, 3, 5, 6, 7]] < 0.1)   # gE1E2: 0.14 from the rounded rows vs 0.2
    # ---- this work: physical masses, M_N = M_n ----
    piN = ext(piN_neutron())
    aD = AP.pols(MN=MNEU)
    piD, pole = ext(aD['piDelta loops']), ext(aD['Delta pole'])
    tot = piN + piD + pole
    print('M_p instead of M_n in the piN loops changes them by', np.round(ext(piN_neutron(MN=MP)) - piN, 3))
    gE_lo, gE_hi = RM.ge_of_rem(RM.rem(-1.0) - 0.5), RM.ge_of_rem(RM.rem(-1.0) + 0.5)
    s_gM = 0.5 * np.abs(ext(AP.pols(MN=MNEU, gM=3.07)['Delta pole']) - ext(AP.pols(MN=MNEU, gM=2.87)['Delta pole']))
    s_gE = 0.5 * np.abs(ext(AP.pols(MN=MNEU, gE=gE_hi)['Delta pole']) - ext(AP.pols(MN=MNEU, gE=gE_lo)['Delta pole']))
    lmp = lmp_rule(piN, piD + pole, tot)
    trunc = np.maximum.reduce([D**2 * np.abs(piN), D * np.abs(piD), D**2 * np.abs(tot)])
    reex = np.sqrt(trunc**2 + s_gM**2 + s_gE**2)
    # proton (driver, physical masses) for the comparison and the p-n differences
    U = json.load(open(os.path.join(RES, 'spin_pol_uncertainties.json')))
    prot = np.array(U['NNLO'])

rows = dict(piN=piN, piDelta=piD, Delta_pole=pole, NNLO=tot, gM_pm0p1=s_gM, REM_pm0p5=s_gE, LMP15_rule=lmp,
            trunc_without_Dpole=trunc, reexamined=reex, LMP15_published=LMP_tot, LMP15_published_err=LMP_err,
            proton_NNLO=prot, n_minus_p=tot - prot)
print('%-22s' % '' + ''.join('%9s' % l[:8] for l in LAB))
for k, v in rows.items():
    print('%-22s' % k + ''.join('%9.3f' % x for x in v))
save_json(dict(labels=LAB, delta=D, M_n=MNEU, m_pi=MPI, gE_REM0p5=(gE_lo, gE_hi), **rows), 'neutron_pols.json')
