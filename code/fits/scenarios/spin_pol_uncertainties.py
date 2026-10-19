"""Re-examination of the uncertainties of the NNLO static polarisabilities (author question, 2026-10-19: "are the
uncertainties of the NNLO prediction for the spin polarisabilities adequate, gamma_0 seems to be too large").
Decomposition of the static values of the amplitude used in the fits (physical masses, Table NNLOpols second row)
into piN loops (O(p^3)), Delta pole and piDelta loops (O(p^{7/2})), and their sensitivity to the gammaNDelta couplings
(g_M = 2.97 +- 0.10, the Delta -> N gamma width; g_E such that R_EM = -2.23 +- 0.5 %, exact relation of lib/rem.py).
Rules:
  LMP15  (Lensky, McGovern, Pascalutsa 2015, Sect. 4.4): Max[d^2 |piN|, d |Delta graphs|, d^2 |total|], d = m_pi/Delta
         = 0.48; reproduces the published uncertainties (check below).
  REEX   (this work): the Delta-pole term is fixed in size by the gammaNDelta couplings determined from the Delta
         region, so its uncertainty is that of g_M and g_E, not d times its size (precedent: LMP15 did not apply the
         rule to alpha, beta either -- it would give 3.1, 2.7 -- and kept +-0.7, for their own reasons); the truncation part is
         Max[d^2 |piN|, d |piDelta|, d^2 |total|], combined in quadrature with the coupling variations.
The rules are applied to gamma_0 and gamma_pi directly (as LMP15 did), and for comparison the four spin
uncertainties are propagated in quadrature (independent).  Writes results/spin_pol_uncertainties.json."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import bchpt
import rem as RM
from fitter import C_GAMMA0, C_GAMMAPI
from common import save_json

D = 0.48
GE_PM = (RM.ge_of_rem(RM.rem(-1.0) - 0.5), RM.ge_of_rem(RM.rem(-1.0) + 0.5))      # g_E for R_EM = default -+ 0.5 %
LAB = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2', 'gamma0', 'gammapi']
PT = [('pref', 100.0, 90.0, 'CS', 'LAB')]


def pref(cfg, tag):
    return bchpt.run(PT, cfg, name='pref_' + tag, verbose=False).Pref.copy()


def ext(P):
    return np.r_[P, C_GAMMA0 @ P, C_GAMMAPI @ P]


R = {c: pref(c, c) for c in ['NNLO', 'NLO', 'BornPi0', 'NNLO_noPiD', 'NNLO_noDpole', 'gM2.87', 'gM3.07']}
for g in GE_PM:
    R['gE%.4f' % g] = pref(dict(gE=g), 'gE%.4f' % g)
piN = ext(R['NLO'] - R['BornPi0']); dpole = ext(R['NNLO_noPiD'] - R['NLO']); piD = ext(R['NNLO_noDpole'] - R['NLO'])
tot = ext(R['NNLO'] - R['BornPi0'])
assert np.allclose(piN + dpole + piD, tot, atol=2e-3), 'mechanisms do not add up'
s_gM = 0.5 * np.abs(ext(R['gM3.07'] - R['gM2.87']))
s_gE = 0.5 * np.abs(ext(R['gE%.4f' % GE_PM[1]] - R['gE%.4f' % GE_PM[0]]))
lmp = np.maximum.reduce([D ** 2 * np.abs(piN), D * np.abs(dpole + piD), D ** 2 * np.abs(tot)])
PUB = np.array([0.70, 0.70, 0.80, 1.50, 0.20, 0.30, 1.40, 1.70])            # Table NNLOpols, first row
assert np.all(np.abs(lmp[2:] - PUB[2:]) < 0.06), ('LMP15 rule not reproduced', lmp)
trunc = np.maximum.reduce([D ** 2 * np.abs(piN), D * np.abs(piD), D ** 2 * np.abs(tot)])
reex = np.sqrt(trunc ** 2 + s_gM ** 2 + s_gE ** 2)
quad = lambda s, c: float(np.sqrt(np.sum((np.abs(c[2:6]) * s[2:6]) ** 2)))     # independent spin pols -> g0, gpi
rows = dict(NNLO=tot, piN=piN, Delta_pole=dpole, piDelta=piD, gM_pm0p1=s_gM, REM_pm0p5=s_gE, LMP15_rule=lmp,
            LMP15_published=PUB, trunc_without_Dpole=trunc, reexamined=reex)
print('%-20s' % '' + ''.join('%9s' % l[:8] for l in LAB))
for k, v in rows.items():
    print('%-20s' % k + ''.join('%9.3f' % x for x in v))
for nm, s in (('LMP15', lmp), ('reexamined', reex)):
    print('%s, four spin pols independent: sigma(gamma_0) = %.2f, sigma(gamma_pi) = %.2f'
          % (nm, quad(s, np.r_[0, 0, -1, -1, -1, -1.]), quad(s, np.r_[0, 0, -1, 1, -1, 1.])))
print('empirical gamma_0 (ForwardRCS 2026) = -1.02(10); NNLO %.2f' % tot[6])
save_json(dict(labels=LAB, delta=D, gE_REM0p5=GE_PM, **rows,
               quad_g0=dict(LMP15=quad(lmp, -C_GAMMA0), reexamined=quad(reex, -C_GAMMA0)),
               quad_gpi=dict(LMP15=quad(lmp, C_GAMMAPI), reexamined=quad(reex, C_GAMMAPI))), 'spin_pol_uncertainties.json')
