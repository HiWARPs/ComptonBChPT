"""Main-fit numbers quoted in the text of the paper (Sects. IV B, IV F, V, App. D), collected from the results files, so
that a refit is propagated by comparison, not by retyping (written 2026-11-12 for the LARA refit).
Run from work/fits after scenarios/combined_fit.py, combined_fit_parts.py, weight_scan.py d038 ... and delta_plots."""
import json, os
import numpy as np
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
DP = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'delta_plots', 'results')
J = json.load(open(os.path.join(RES, 'combined_fit.json')))
FL = json.load(open(os.path.join(RES, 'fits_L.json')))['F2_d_all_thB']
CG0, CGP = np.array([0, 0, -1, -1, -1, -1.]), np.array([0, 0, -1, 1, -1, 1.])
NM = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']


def fit_line(k):
    c = J[k]; P = np.array(c['P']); C = np.array(c['covP']); e = np.sqrt(np.diag(C))
    s = '%-26s ' % k + ' '.join('%s=%.2f(%.2f)' % (n, p, d) for n, p, d in zip(NM, P, e))
    s += ' a+b=%.2f(%.2f) a-b=%.2f(%.2f)' % (*c['der']['apb'], *c['der']['amb'])
    s += ' g0=%.2f(%.2f) gpi=%.2f(%.2f)' % (CG0 @ P, np.sqrt(CG0 @ C @ CG0), CGP @ P, np.sqrt(CGP @ C @ CGP))
    s += ' gE=%.3f(%.3f) REM=%.2f(%.2f) rho=%.2f chi2=%.1f/%d N=%d (%d+%d)' % (*c['gE'], *c['REM'], c['rho_ab'], c['chi2'],
                                                                                c['dof'], c['N'], c['N_L'], c['N_D'])
    return s


print('fit below threshold: a=%.2f b=%.2f' % tuple(FL['P'][:2]))
for k in [k for k in J if k.startswith('combined')]:
    print(fit_line(k))
m = J['combined']
print('main - below threshold: da=%.2f db=%.2f' % (m['P'][0] - FL['P'][0], m['P'][1] - FL['P'][1]))
print('chi2 of the sets below threshold at the main fit (weighted fit, own normalizations): %.1f' % m['chi2_below150_sets'])
for k in ('delta_only', 'delta_only_withLEGS', 'delta_only_noLARA'):
    if k in J:
        print('%-22s REM=%.2f(%.2f)%% gE=%.3f(%.3f) N=%d chi2=%.1f' % (k, *J[k]['REM'], *J[k]['gE'], J[k]['N'], J[k]['chi2']))
print('per set (main):', {d: (round(v['chi2'], 1), v.get('N')) for d, v in m['perset'].items()})
# half-ranges over w and over the prior widths (Sect. IV B)
ws = ['combined_w0.01', 'combined_w0.03', 'combined_w0.1', 'combined', 'combined_w1']
pr = ['combined_spin1sig', 'combined_spin2sig', 'combined']
for lab, ks in (('w 0.01..1', ws), ('priors 1..3', pr)):
    A = np.array([J[k]['P'][:2] for k in ks])
    print('half-ranges %s: alpha %.2f beta %.2f' % (lab, *(0.5 * (A.max(0) - A.min(0)))))
print('with LEGS shift: da=%.2f db=%.2f d(a-b)=%.2f' % (J['combined_withLEGS']['P'][0] - m['P'][0],
                                                     J['combined_withLEGS']['P'][1] - m['P'][1],
                                                     J['combined_withLEGS']['der']['amb'][0] - m['der']['amb'][0]))
# unweighted chi^2 parts (combined_fit_parts.json)
fp = os.path.join(RES, 'combined_fit_parts.json')
if os.path.exists(fp):
    print('parts:', json.dumps(json.load(open(fp)), indent=None)[:3000])
fw = os.path.join(RES, 'weight_scan_d038.json')
if os.path.exists(fw):
    W = json.load(open(fw))
    for k, r in W.items():
        p = r['parts']
        print('weight scan %-18s a=%.2f b=%.2f REM=%.2f parts L=%.1f D=%.1f FSR=%.1f total=%.1f' % (
            k, r['P'][0], r['P'][1], r['REM'][0], p['L'], p['D'], p['FSR'], p['L'] + p['D'] + p['FSR']))
for f in ('chi2_table_paper.tex',):
    fn = os.path.join(DP, f)
    if os.path.exists(fn):
        print(open(fn).read())

# ---- checklist pulls (Table checklist): |found - predicted| / sqrt(sig_found^2 + sig_pred^2); predicted = NNLO with
#      physical masses (Tab. NNLOpols, second row) and the re-examined uncertainties (third row); alpha, beta found with
#      the systematic uncertainty of Eq. mainresult (0.45, 0.37) in quadrature; R_EM: default -2.23(50)%, PDG -2.5(3)%
PRED = np.array([11.42, 3.73, -3.26, 2.94, 0.17, 1.13]); SP = np.array([0.70, 0.70, 0.47, 0.48, 0.08, 0.19])
SYST = np.array([0.45, 0.37, 0, 0, 0, 0])
P = np.array(m['P']); C = np.array(m['covP']); e = np.sqrt(np.diag(C) + SYST ** 2)
for n, f, de, p, sp in zip(NM, P, e, PRED, SP):
    print('pull %-8s found %.2f(%.2f) pred %.2f(%.2f): %.2f sigma' % (n, f, de, p, sp, abs(f - p) / np.hypot(de, sp)))
gp, dgp = CGP @ P, np.sqrt(CGP @ C @ CGP)
print('pull gpi      found %.2f(%.2f) pred 7.16(1.08): %.2f sigma' % (gp, dgp, abs(gp - 7.16) / np.hypot(dgp, 1.08)))
r, dr = m['REM']
print('pull R_EM     found %.2f(%.2f): default %.2f sigma, PDG %.2f sigma' % (r, dr, abs(r + 2.23) / np.hypot(dr, 0.5),
                                                                          abs(r + 2.5) / np.hypot(dr, 0.3)))
print('alpha, beta total: %.2f +- %.2f (fit) +- 0.45 = %.2f; %.2f +- %.2f +- 0.37 = %.2f; a-b total (syst linear, anticorr.) %.2f(%.2f)' % (
    P[0], np.sqrt(C[0, 0]), np.hypot(np.sqrt(C[0, 0]), 0.45), P[1], np.sqrt(C[1, 1]), np.hypot(np.sqrt(C[1, 1]), 0.37),
    m['der']['amb'][0], np.hypot(m['der']['amb'][1], 0.45 + 0.37)))
