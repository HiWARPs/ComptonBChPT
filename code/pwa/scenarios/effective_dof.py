"""Effective number of parameters of the PWA fit and of the BChPT fit below threshold (author question 2026-10-20: "is
chi^2 per point an adequate comparison between the BChPT and PWA fit, given that the latter has many more free parameters,
and they fit over different regions").  Both fits are rebuilt at their minima: the PWA fit of the paper
(scenarios/main_fits_legacynorm.py MAIN; 3 global parameters, 6 residual functions in each of 11 bins, 15 normalisations,
2 forward nuisances, priors) and the BChPT Baldin fit without truncation covariance (fits_L.json 'F1_d_all'; alpha, beta,
15 normalisations), i.e. the pair of the comparison "Which of the two frameworks describes the data better?" (same error
model, normalisations floated).  The effective number of parameters determined by a set of data rows B is the sum of
their leverages, p_eff(B) = tr(J_B H^-1 J_B^T), H = J^T J of all residual rows (data, priors, constraints, forward
SRs); summed over all rows it equals the number of parameters.  Compared on the common points: the 181 points of the
modern experiments and the 237 points of the PWA database (the BChPT database has 8 more, Federspiel et al. at 32 MeV).
Writes results/effective_dof.{json,md}.  Run from work/pwa/."""
import json
import os
import sys
import numpy as np
from scipy.stats import chi2 as chi2dist

PWA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.dirname(PWA)
sys.path[:0] = [PWA, os.path.join(PWA, 'scenarios')]
os.chdir(PWA)
import update_fits as UF
import main_fits_legacynorm as L

MODERN = {'pwa': ['A2_2022_CS', 'A2_2022_S3', 'TAPS2001', 'Hallin1993', 'MacGibbon1995_tagged', 'MacGibbon1995_untagged'],
          'bchpt': ['A2_2022_CS', 'A2_2022_S3', 'OlmosdeLeon2001', 'Hallin1993_TabII', 'MacGibbon1995_tagged',
                    'MacGibbon1995_untagged']}


def leverages(J):
    return np.einsum('ij,jk,ik->i', J, np.linalg.inv(J.T @ J), J)


# ---------------- PWA ----------------
m = UF.build_model(UF.select(L.load_leg()), UF.BINS_A, **L.MAIN)
r = m.fit(with_prior=True)
rp, Jp = m.residuals(r.p, True, jac=True)
Jp = Jp[:, r.free]
nd = int(m.mask.sum())
lp = leverages(Jp)[:nd]
dsp, Ep = m.d.dataset.values[np.where(m.mask)[0]], m.d.E.values[np.where(m.mask)[0]]
chip = rp[:nd] ** 2                        # whitened rows; block sums are exact (correlated blocks lie within one set)
ref = json.load(open('results/main_fits_legacynorm.json'))['MAIN']['chi2']
assert abs(chip.sum() - ref) < 1e-6 * ref, (chip.sum(), ref)

# ---------------- BChPT (fits_L.json F1_d_all) ----------------
os.chdir(os.path.join(WORK, 'fits'))
sys.path[:0] = [os.path.join(WORK, 'fits', 'scripts'), os.path.join(WORK, 'fits', 'lib')]
import s03_fits_L as S
from common import build_fit, PREF_NNLO
from fitter import C_BALDIN
s = S.rows_for(S.df, S.ALL, 'F1')
F = build_fit(s, S.T, PREF_NNLO, ['alphaE1', 'betaM1'], constraints=[(C_BALDIN, S.BALDIN[0], S.BALDIN[1], 'Baldin')],
              fwd_mode='nuisance', thcov=None)
q, _, _ = F.fit()
Jb = F._jac(F.residuals, q)
rb = F.residuals(q)
lb = leverages(Jb)[:F.N]
chib = rb[:F.N] ** 2                       # diagonal errors: one row per point
dsb, Eb = s.dataset.values, s.E.values
ref = json.load(open(os.path.join(WORK, 'fits', 'results', 'fits_L.json')))['F1_d_all']['chi2_data']
assert abs(chib.sum() - ref) < 1e-6 * ref, (chib.sum(), ref)
os.chdir(PWA)

# ---------------- comparison on common points ----------------
common_b = ~((dsb == 'Federspiel1991') & (Eb < 55.0))         # the 8 points at 32 MeV are outside the PWA bins
assert common_b.sum() == nd, (common_b.sum(), nd)
out = dict(npar=dict(pwa=int(Jp.shape[1]), bchpt=int(Jb.shape[1])),
           p_eff_all_rows=dict(pwa=float(leverages(Jp).sum()), bchpt=float(leverages(Jb).sum())))
md = ['### Effective number of parameters: PWA fit vs BChPT `F1_d_all` (Baldin sum rule, no truncation covariance)\n',
      'p_eff = sum of the leverages of the data rows (diagonal of the hat matrix); nu = N - p_eff.\n',
      '| points | N | chi2 PWA | p_eff PWA | chi2/nu PWA | chi2 BChPT | p_eff BChPT | chi2/nu BChPT | Delta chi2 | Delta p_eff | p | Delta AIC |',
      '|---|---|---|---|---|---|---|---|---|---|---|---|']
for lab, mp, mb in (('modern (A2 2022, A2 2001, Saskatoon)', np.isin(dsp, MODERN['pwa']), np.isin(dsb, MODERN['bchpt'])),
                    ('PWA database (237)', np.ones(nd, bool), common_b)):
    assert mp.sum() == mb.sum()
    N = int(mp.sum())
    cp, pp, cb, pb = chip[mp].sum(), lp[mp].sum(), chib[mb].sum(), lb[mb].sum()
    dc, dp = cb - cp, pp - pb
    pv = float(chi2dist.sf(dc, dp))       # chance of a chi^2 gain >= dc from dp extra effective parameters
    daic = (cp + 2 * pp) - (cb + 2 * pb)
    out[lab] = dict(N=N, chi2_pwa=float(cp), p_eff_pwa=float(pp), chi2_bchpt=float(cb), p_eff_bchpt=float(pb),
                    chi2_per_nu_pwa=float(cp / (N - pp)), chi2_per_nu_bchpt=float(cb / (N - pb)), dchi2=float(dc),
                    dp_eff=float(dp), p_gain=pv, dAIC_pwa_minus_bchpt=float(daic))
    md.append('| %s | %d | %.1f | %.1f | %.3f | %.1f | %.1f | %.3f | %.1f | %.1f | %.2f | %+.1f |'
              % (lab, N, cp, pp, cp / (N - pp), cb, pb, cb / (N - pb), dc, dp, pv, daic))
# PWA: where its data-determined dof sit (share of each parameter group, diag of H^-1 J_d^T J_d)
Jd = Jp[:nd]
share = np.diag(np.linalg.inv(Jp.T @ Jp) @ (Jd.T @ Jd))
names = np.array(m.names)[r.free]
groups = {'global (alpha-beta, gE1E1, gM1E2)': ['D', 'gEE', 'gME'], 'forward nuisances': ['nf', 'ng']}
for k in ('fR1', 'fR2', 'fR3', 'fR4', 'fR5', 'fR6'):
    groups['residual functions ' + k] = [n for n in names if n.startswith(k + '@')]
groups['normalisations'] = [n for n in names if n.startswith('norm:')]
out['pwa_share'] = {g: float(sum(share[list(names).index(n)] for n in ns)) for g, ns in groups.items()}
md += ['', 'PWA: data-determined dof by parameter group (sum %.2f):' % share.sum()] + \
      ['- %s: %.2f' % (g, v) for g, v in out['pwa_share'].items()]
json.dump(out, open('results/effective_dof.json', 'w'), indent=1)
open('results/effective_dof.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))
