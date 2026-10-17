"""gE_scan.py -- sensitivity of the Delta-region observables to the gamma-N-Delta electric coupling g_E (E2/M1 ratio):
complete NNLO prediction (driver Pref, physical masses) recomputed for g_E = 0, -0.5, -1.0 (paper), -1.5, -2.0 at
g_M = 2.97, on the figure grids and at the data points of the Delta-region chi^2 table (150 < E <= 350 MeV, n = 1 and
with the normalisation fitted within its quoted uncertainty).  Two evaluations per g_E: 'full' = the prediction at its
own static polarisabilities (as the g_M band of the figures), 'P0' = at the static polarisabilities of the paper
(contact terms absorb the g_E dependence of the static values; only the energy dependence changes).
R_EM at the Delta pole: R_EM = (M_D - M_N)/(2 (M_D + M_N)) * g_E/g_M = 0.0677 g_E/g_M (PV2006 Eq. for G_E*/G_M*),
i.e. g_E = -1.0 <-> R_EM = -2.3 %.   Results -> results/gE_scan/ (csv + md).   usage: python gE_scan.py"""
import os, sys, json
os.environ.setdefault('DP_EXTRA', 'none'); os.environ.setdefault('DP_MAIN', 'pred')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import config as C, data_dp as D, make_figs as MF, runner

GE_LIST = [0.0, -0.5, -1.0, -1.5, -2.0]
REM = lambda gE, gM=C.GM0: (1232.0 - C.MP) / (2 * (1232.0 + C.MP)) * gE / gM
OUT = os.path.join(C.RESULTS, 'gE_scan'); os.makedirs(OUT, exist_ok=True)

df = D.build(); pts = MF.build_points(df)
seen = set(); P = []
for p in pts:
    if p[0] not in seen:
        seen.add(p[0]); P.append(p)
ix = {p[0]: i for i, p in enumerate(P)}
runs = {}
for gE in GE_LIST:
    tmpl = runner.DEFAULT_TEMPLATE.replace('dng2=-1.0', 'dng2=%s' % repr(float(gE)))
    assert gE == -1.0 or tmpl != runner.DEFAULT_TEMPLATE
    runs[gE] = runner.run(P, template=tmpl, name='dp_gEscan_%+.1f' % gE)
P0 = runs[-1.0].Pref.copy()
print('static polarisabilities (Pref) vs g_E  [alpha, beta, gE1E1, gM1M1, gE1M2, gM1E2]:')
for gE in GE_LIST:
    print('  gE=%5.1f  R_EM=%5.2f%%  Pref=%s' % (gE, 100 * REM(gE), np.round(runs[gE].Pref, 3)))
cur = {gE: dict(full=runs[gE].predict_array(runs[gE].Pref), P0=runs[gE].predict_array(P0)) for gE in GE_LIST}

# ---- curves: relative / absolute changes w.r.t. gE = -1.0 on the figure grids
rows = []
def grid(prefix, E=None):
    return [p for p in P if p[0].startswith(prefix)]
for gE in GE_LIST:
    for mode in ('full', 'P0'):
        c = cur[gE][mode]; c0 = cur[-1.0][mode]
        # dsigma vs E at the c.m. angles (Fig. 1 lower row), 250-350 MeV
        for fr, th in MF.FIG1:
            sel = [p for p in P if p[0].startswith('g_f1_%s_%s_' % (fr, th)) and 250 <= p[1] <= 350]
            ii = np.array([ix[p[0]] for p in sel]); r = c[ii] / c0[ii] - 1
            rows.append(dict(gE=gE, mode=mode, obs='dsigma %s %g' % (fr, th), stat='max |rel. change| 250-350 MeV', value=float(np.max(np.abs(r)))))
            sel = [p for p in P if p[0].startswith('g_f1_%s_%s_' % (fr, th)) and 300 <= p[1] <= 340]
            ii = np.array([ix[p[0]] for p in sel]); r = c[ii] / c0[ii] - 1
            rows.append(dict(gE=gE, mode=mode, obs='dsigma %s %g' % (fr, th), stat='mean rel. change at the peak (300-340)', value=float(np.mean(r))))
        # Sigma_3 at the LEGS energies (Fig. 3a) and Sigma_3(90) vs E (Fig. 3b)
        sel = [p for p in P if p[0].startswith('g_f3_')]; ii = np.array([ix[p[0]] for p in sel])
        rows.append(dict(gE=gE, mode=mode, obs='Sigma3 (LEGS energies, all angles)', stat='max |abs. change|', value=float(np.max(np.abs(c[ii] - c0[ii])))))
        for E0 in (213.1, 254.1, 286.5, 323.8, 333.8):
            sel = [p for p in P if p[0].startswith('g_f3_%s_' % E0) and 60 <= p[2] <= 140]; ii = np.array([ix[p[0]] for p in sel])
            rows.append(dict(gE=gE, mode=mode, obs='Sigma3 at %g MeV, 60-140 deg' % E0, stat='max |abs. change|', value=float(np.max(np.abs(c[ii] - c0[ii])))))
        sel = [p for p in P if p[0].startswith('g_f3b_') and 250 <= p[1] <= 350]; ii = np.array([ix[p[0]] for p in sel])
        rows.append(dict(gE=gE, mode=mode, obs='Sigma3(90) 250-350 MeV', stat='max |abs. change|', value=float(np.max(np.abs(c[ii] - c0[ii])))))
        for o, E0, lab in MF.FIG4:
            sel = [p for p in P if p[0].startswith('g_f4_%s_%s_' % (o, E0)) and 40 <= p[2] <= 160]; ii = np.array([ix[p[0]] for p in sel])
            rows.append(dict(gE=gE, mode=mode, obs='%s %s' % (o, lab), stat='max |abs. change|', value=float(np.max(np.abs(c[ii] - c0[ii])))))
curves_tab = pd.DataFrame(rows); curves_tab.to_csv(os.path.join(OUT, 'gE_curve_changes.csv'), index=False)

# ---- chi^2 per data set vs gE (n = 1 and with the normalisation fitted), parabola -> preferred gE
s = df[(df.E > C.E_DELTA_MIN) & (df.E <= C.E_DELTA_MAX) & df.obs.isin(['CS', 'CSCM', 'S3', 'S2X', 'S2Z'])].copy()
s['set'] = s.dataset.replace({'Hallin1993_TabII': 'Hallin1993_TabII(D)', 'Hallin1993_TabII_D': 'Hallin1993_TabII(D)'})
def chi2_set(q, t):
    y, e = q.value.values, q.err.values; sn = float(q.norm.iloc[0])
    c1 = float(np.sum(((y - t) / e) ** 2))
    if sn > 0:
        n = (np.sum(y * t / e ** 2) + 1 / sn ** 2) / (np.sum(t * t / e ** 2) + 1 / sn ** 2)
        cn = float(np.sum(((y - n * t) / e) ** 2) + ((n - 1) / sn) ** 2)
    else:
        n, cn = 1.0, c1
    return c1, cn, float(n)
crow = []
for ds, lab in MF.PAPERROWS:
    q = s[s.set == ds]
    if not len(q):
        continue
    ii = np.array([ix[p] for p in q.pid])
    for mode in ('full', 'P0'):
        c1s, cns, ns = [], [], []
        for gE in GE_LIST:
            c1, cn, n = chi2_set(q, cur[gE][mode][ii]); c1s.append(c1); cns.append(cn); ns.append(n)
        g = np.array(GE_LIST)
        rec = dict(set=ds, mode=mode, N=len(q))
        for k, gE in enumerate(GE_LIST):
            rec['chi2_n1_%+.1f' % gE] = round(c1s[k], 2); rec['chi2_nfit_%+.1f' % gE] = round(cns[k], 2); rec['n_%+.1f' % gE] = round(ns[k], 3)
        for key, arr in (('n1', c1s), ('nfit', cns)):
            a, b, c = np.polyfit(g, arr, 2)
            gmin = -b / (2 * a) if a > 0 else np.nan
            rec['gE_min_' + key] = round(float(gmin), 2) if np.isfinite(gmin) else np.nan
            rec['gE_err_' + key] = round(float(1 / np.sqrt(a)), 2) if a > 0 else np.nan
            rec['REM_min_' + key] = round(100 * REM(gmin), 2) if np.isfinite(gmin) else np.nan
            rec['REM_err_' + key] = round(100 * abs(REM(1 / np.sqrt(a))), 2) if a > 0 else np.nan
            rec['dchi2_' + key] = round(float(np.polyval([a, b, c], gmin) - arr[2]), 2) if a > 0 else np.nan   # chi2(min) - chi2(-1.0)
        crow.append(rec)
chi_tab = pd.DataFrame(crow); chi_tab.to_csv(os.path.join(OUT, 'gE_chi2_sets.csv'), index=False)

# ---- report
with open(os.path.join(OUT, 'gE_scan.md'), 'w') as f:
    f.write('# Sensitivity of the Delta-region observables to g_E (E2/M1)\n\nR_EM = 0.0677 g_E/g_M; g_E = 0, -0.5, -1.0, -1.5, -2.0 <-> R_EM = %s %% (g_M = 2.97)\n\n' % ', '.join('%.1f' % (100 * REM(g)) for g in GE_LIST))
    f.write('## Static polarisabilities (Pref) vs g_E\n\n')
    for gE in GE_LIST:
        f.write('- g_E = %.1f: %s\n' % (gE, np.round(runs[gE].Pref, 3).tolist()))
    f.write('\n## Changes of the curves relative to g_E = -1.0\n\n')
    f.write(curves_tab.pivot_table(index=['obs', 'stat', 'mode'], columns='gE', values='value').round(3).to_string() + '\n')
    f.write('\n## chi^2 per data set (150 < E <= 350 MeV) vs g_E\n\n')
    pd.set_option('display.width', 250)
    cols = ['set', 'mode', 'N'] + ['chi2_nfit_%+.1f' % g for g in GE_LIST] + ['gE_min_nfit', 'gE_err_nfit', 'REM_min_nfit', 'REM_err_nfit', 'dchi2_nfit', 'gE_min_n1', 'gE_err_n1', 'dchi2_n1']
    f.write(chi_tab[cols].to_string(index=False) + '\n')
print(open(os.path.join(OUT, 'gE_scan.md')).read())
