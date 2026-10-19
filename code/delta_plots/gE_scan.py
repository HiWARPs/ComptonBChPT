"""gE_scan.py -- sensitivity of the Delta-region observables to the gamma-N-Delta electric coupling g_E (E2/M1 ratio):
complete NNLO prediction (driver Pref, physical masses) recomputed for g_E = 0, -0.5, -1.0 (paper), -1.5, -2.0 at
g_M = 2.97, on the figure grids and at the data points of the Delta-region chi^2 table (150 < E <= 350 MeV, n = 1 and
with the normalisation fitted within its quoted uncertainty).  Two evaluations per g_E: 'full' = the prediction at its
own static polarisabilities (as the g_M band of the figures), 'P0' = at the static polarisabilities of the paper
(contact terms absorb the g_E dependence of the static values; only the energy dependence changes).
R_EM at the Delta pole: R_EM = -G_E*/G_M* (exact, fits/lib/rem.py, PV2006; until 2026-10-20 the leading-order form 0.0677 g_E/g_M),
i.e. g_E = -1.0 <-> R_EM = -2.23 %.   Results -> results/gE_scan/ (csv + md).   usage: python gE_scan.py"""
import os, sys, json
os.environ.setdefault('DP_EXTRA', 'none'); os.environ.setdefault('DP_MAIN', 'pred')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import config as C, data_dp as D, make_figs as MF, runner

GE_LIST = [0.0, -0.5, -1.0, -1.5, -2.0]
REM, dREM = C.RM.rem, C.RM.drem_dge                         # [%]
gradR = lambda gE, gM: np.array([-REM(gE, gM) / (gM - C.RM.K * gE), dREM(gE, gM)])   # d R_EM/d(g_M, g_E) [%]
OUT = os.path.join(C.RESULTS, 'gE_scan'); os.makedirs(OUT, exist_ok=True)

df = D.build(); pts = MF.build_points(df)
seen = set(); P = []
for p in pts:
    if p[0] not in seen:
        seen.add(p[0]); P.append(p)
ix = {p[0]: i for i, p in enumerate(P)}
runs = {}
for gE in GE_LIST:
    runs[gE] = runner.run(P, gE=gE, name='dp_gEscan_%+.1f' % gE)
P0 = runs[-1.0].Pref.copy()
print('static polarisabilities (Pref) vs g_E  [alpha, beta, gE1E1, gM1M1, gE1M2, gM1E2]:')
for gE in GE_LIST:
    print('  gE=%5.1f  R_EM=%5.2f%%  Pref=%s' % (gE, REM(gE), np.round(runs[gE].Pref, 3)))
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
            rec['REM_min_' + key] = round(REM(gmin), 3) if np.isfinite(gmin) else np.nan
            rec['REM_err_' + key] = round(dREM(gmin) / np.sqrt(a), 3) if a > 0 and np.isfinite(gmin) else np.nan
            rec['dchi2_' + key] = round(float(np.polyval([a, b, c], gmin) - arr[2]), 2) if a > 0 else np.nan   # chi2(min) - chi2(-1.0)
        crow.append(rec)
chi_tab = pd.DataFrame(crow); chi_tab.to_csv(os.path.join(OUT, 'gE_chi2_sets.csv'), index=False)

# ---- joint (g_M, g_E) scan: 3 x 3 grid, quadratic form chi^2(gM, gE) per set -> joint minimum, errors, correlation
GM_J = [C.GM_LO, C.GM0, C.GM_HI]; GE_J = [-1.5, -1.0, -0.5]
jr = {}
for gM in GM_J:
    for gE in GE_J:
        if gM == C.GM0:
            jr[(gM, gE)] = runs[gE]
        elif gE == -1.0:
            jr[(gM, gE)] = runner.run(P, gM=gM, name='all_NNLO_VC_gM%.2f' % gM)        # the cached g_M-band runs of theory_dp
        else:
            jr[(gM, gE)] = runner.run(P, gM=gM, gE=gE, name='dp_gEscan_gM%.2f_gE%+.1f' % (gM, gE))
jc = {k: v.predict_array(v.Pref) for k, v in jr.items()}          # complete prediction at its own static polarisabilities
jrow = []
for ds, lab in MF.PAPERROWS:
    q = s[s.set == ds]
    if not len(q):
        continue
    ii = np.array([ix[p] for p in q.pid])
    for key in ('n1', 'nfit'):
        X, Y, Z = [], [], []
        for (gM, gE), c in jc.items():
            c1, cn, n = chi2_set(q, c[ii]); X.append(gM - C.GM0); Y.append(gE - C.GE0); Z.append(c1 if key == 'n1' else cn)
        X, Y, Z = map(np.array, (X, Y, Z))
        A = np.stack([np.ones_like(X), X, Y, X * X, X * Y, Y * Y], axis=1)
        coef, *_ = np.linalg.lstsq(A, Z, rcond=None)
        c0, bx, by, axx, axy, ayy = coef
        H = np.array([[2 * axx, axy], [axy, 2 * ayy]])                 # Hessian of chi^2
        rec = dict(set=ds, key=key, N=len(q), chi2_ref=round(float(Z[(X == 0) & (Y == 0)][0]), 2), fit_resid=round(float(np.max(np.abs(A @ coef - Z))), 2))
        try:
            ev = np.linalg.eigvalsh(H)
            if np.all(ev > 0):
                m = np.linalg.solve(H, -np.array([bx, by]))
                cov = 2 * np.linalg.inv(H)                            # Delta chi^2 = 1 errors
                rec.update(gM_min=round(float(m[0] + C.GM0), 3), gE_min=round(float(m[1] + C.GE0), 3),
                           gM_err=round(float(np.sqrt(cov[0, 0])), 3), gE_err=round(float(np.sqrt(cov[1, 1])), 3),
                           corr=round(float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])), 2),
                           REM_min=round(REM(m[1] + C.GE0, m[0] + C.GM0), 2),
                           REM_err=round(float(np.sqrt(gradR(m[1] + C.GE0, m[0] + C.GM0) @ cov @ gradR(m[1] + C.GE0, m[0] + C.GM0))), 2),
                           chi2_min=round(float(c0 + bx * m[0] + by * m[1] + axx * m[0] ** 2 + axy * m[0] * m[1] + ayy * m[1] ** 2), 2))
                # g_E at fixed g_M = 2.97 (conditional) for comparison with the 1D scan
                rec['gE_cond'] = round(float(-by / (2 * ayy) + C.GE0), 3) if ayy > 0 else np.nan
                rec['gE_cond_err'] = round(float(1 / np.sqrt(ayy)), 3) if ayy > 0 else np.nan
            else:
                rec.update(note='chi^2 not convex in (gM, gE): eigenvalues %s' % np.round(ev, 2).tolist())
        except np.linalg.LinAlgError:
            rec.update(note='singular')
        jrow.append(rec)
joint_tab = pd.DataFrame(jrow); joint_tab.to_csv(os.path.join(OUT, 'gE_gM_joint.csv'), index=False)

# ---- robust correlation measure: g_E minimum at each fixed g_M (parabola through g_E = -1.5, -1.0, -0.5), n fitted;
#      slope d(g_E min)/d(g_M) from g_M = 2.87 and 3.07 -> shift of the preferred R_EM per +0.1 in g_M
crow2 = []
for ds, lab in MF.PAPERROWS:
    q = s[s.set == ds]
    if not len(q):
        continue
    ii = np.array([ix[p] for p in q.pid])
    rec = dict(set=ds, N=len(q))
    mins = {}
    for gM in GM_J:
        g = np.array(GE_J); c2 = np.array([chi2_set(q, jc[(gM, gE)][ii])[1] for gE in GE_J])
        a_, b_, c_ = np.polyfit(g, c2, 2)
        mins[gM] = (-b_ / (2 * a_)) if a_ > 0 else np.nan
        rec['gE_min_gM%.2f' % gM] = round(float(mins[gM]), 2) if np.isfinite(mins[gM]) else np.nan
        rec['gE_err_gM%.2f' % gM] = round(float(1 / np.sqrt(a_)), 2) if a_ > 0 else np.nan
        rec['chi2min_gM%.2f' % gM] = round(float(np.polyval([a_, b_, c_], mins[gM])), 1) if a_ > 0 else np.nan
    if all(np.isfinite(v) for v in mins.values()):
        slope = (mins[C.GM_HI] - mins[C.GM_LO]) / (C.GM_HI - C.GM_LO)
        rec['dgE_per_0.1gM'] = round(float(0.1 * slope), 2)
        rec['dREM_per_0.1gM'] = round(float(0.1 * (REM(mins[C.GM_HI], C.GM_HI) - REM(mins[C.GM_LO], C.GM_LO)) / (C.GM_HI - C.GM_LO)), 3)
        for gM in GM_J:                                   # preferred R_EM at each fixed g_M (with that g_M in R_EM)
            rec['REM_min_gM%.2f' % gM] = round(float(REM(mins[gM], gM)), 2)
            rec['REM_err_gM%.2f' % gM] = round(float(dREM(mins[gM], gM) * rec['gE_err_gM%.2f' % gM]), 2)
    crow2.append(rec)
cond_tab = pd.DataFrame(crow2); cond_tab.to_csv(os.path.join(OUT, 'gE_cond_min_vs_gM.csv'), index=False)

# ---- report
with open(os.path.join(OUT, 'gE_scan.md'), 'w') as f:
    f.write('# Sensitivity of the Delta-region observables to g_E (E2/M1)\n\nR_EM = -G_E*/G_M* (exact); g_E = 0, -0.5, -1.0, -1.5, -2.0 <-> R_EM = %s %% (g_M = 2.97)\n\n' % ', '.join('%.2f' % REM(g) for g in GE_LIST))
    f.write('## Static polarisabilities (Pref) vs g_E\n\n')
    for gE in GE_LIST:
        f.write('- g_E = %.1f: %s\n' % (gE, np.round(runs[gE].Pref, 3).tolist()))
    f.write('\n## Changes of the curves relative to g_E = -1.0\n\n')
    f.write(curves_tab.pivot_table(index=['obs', 'stat', 'mode'], columns='gE', values='value').round(3).to_string() + '\n')
    f.write('\n## chi^2 per data set (150 < E <= 350 MeV) vs g_E\n\n')
    pd.set_option('display.width', 250)
    cols = ['set', 'mode', 'N'] + ['chi2_nfit_%+.1f' % g for g in GE_LIST] + ['gE_min_nfit', 'gE_err_nfit', 'REM_min_nfit', 'REM_err_nfit', 'dchi2_nfit', 'gE_min_n1', 'gE_err_n1', 'dchi2_n1']
    f.write(chi_tab[cols].to_string(index=False) + '\n')
    f.write('\n## Joint (g_M, g_E) scan (3 x 3 grid, quadratic form per set; errors = Delta chi^2 = 1; corr = correlation of g_M and g_E)\n\n')
    f.write(joint_tab.to_string(index=False) + '\n')
    f.write('\n## g_E minimum at fixed g_M = 2.87 / 2.97 / 3.07 (parabola in g_E, n fitted) and the shift of the preferred value per +0.1 in g_M\n\n')
    f.write(cond_tab.to_string(index=False) + '\n')
    f.write('\n## Sensitivity per 1 %% of R_EM (mean of the g_E = -1.5, -0.5 changes, each per its own change of R_EM)\n\n')
    t = curves_tab[curves_tab['mode'] == 'full'].pivot_table(index=['obs', 'stat'], columns='gE', values='value')
    t['per_1pct'] = 0.5 * (t[-0.5].abs() / abs(REM(-0.5) - REM(-1.0)) + t[-1.5].abs() / abs(REM(-1.5) - REM(-1.0)))
    f.write(t[['per_1pct']].round(3).to_string() + '\n')
print(open(os.path.join(OUT, 'gE_scan.md')).read())
