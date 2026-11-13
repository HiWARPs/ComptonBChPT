"""Rows of the paper's Table 'variants' (and the main-fit numbers quoted in the text) from results/combined_fit.json.
Written 2026-11-11 for the delta = Delta/m_rho refit, so that the table is regenerated, not retyped.
Run from work/fits: ../../.venv/bin/python scripts/s15_variants_rows.py"""
import json, os, numpy as np
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
J = json.load(open(os.path.join(RES, 'combined_fit.json')))
CG0, CGP = np.array([0, 0, -1, -1, -1, -1.]), np.array([0, 0, -1, 1, -1, 1.])


def pe(v, e, nd):                                       # 11.00(50), -2.2(1.0), 8.2(8)
    if e is None or e == 0:
        return '$%.*f$' % (nd, v)
    if nd == 2:
        return '$%.2f(%d)$' % (v, round(e * 100))
    return '$%.1f(%.1f)$' % (v, e) if e >= 0.95 else '$%.1f(%d)$' % (v, round(e * 10))


def row(label, key, gpi_err=True):
    c = J[key]; P = np.array(c['P']); C = np.array(c['covP']); e = np.sqrt(np.diag(C))
    rem = c['REM'] if isinstance(c['REM'], list) else [c['REM'], 0.0]
    gp = CGP @ P; dgp = np.sqrt(CGP @ C @ CGP)
    return '%s & %s & %s & %s & %s\\\\' % (label, pe(P[0], e[0], 2), pe(P[1], e[1], 2), pe(rem[0], rem[1], 1),
                                          pe(gp, dgp if (gpi_err and dgp > 1e-6) else None, 1))


# 2026-11-12 (author): no variants in the weight w in the paper (rows w = 1, 0.1, 0.03, 0.01, 'w = 0.01, 1x',
# 'SRs > nu0, w = 1', 'With LEGS, w = 1' removed; the fits remain in combined_fit.json)
ROWS = [('Main fit', 'combined'), None,
        ('Spin pol.\\ NNLO', 'combined_spinNNLO'), ('Priors $1\\times$', 'combined_spin1sig'),
        ('Priors $2\\times$', 'combined_spin2sig'), ('$3\\times$ Ref.~\\cite{Lensky:2015awa}', 'combined_spinLMP15'), None,
        ('$R_{EM}$ fixed', 'combined_REMfixed'), ('$\\EL\\leq 350$\\,MeV', 'combined_E350'),
        ('No SRs $>\\nu_0$', 'combined_noFSR'), ('5-MeV grid', 'combined_FSR5'), ('20-MeV grid', 'combined_FSR20'), None,
        ('Without LARA', 'combined_noLARA'), ('With LEGS', 'combined_withLEGS')]
for r in ROWS:
    print('\\hline' if r is None else row(*r))
c = J['combined']; P = np.array(c['P']); C = np.array(c['covP'])
print('% main fit: P', np.round(P, 3), 'err', np.round(np.sqrt(np.diag(C)), 3), 'REM', np.round(c['REM'], 3),
      'g0 %.3f(%.3f) gpi %.3f(%.3f) chi2 %.1f/%d' % (CG0 @ P, np.sqrt(CG0 @ C @ CG0), CGP @ P, np.sqrt(CGP @ C @ CGP), c['chi2'], c['dof']))
