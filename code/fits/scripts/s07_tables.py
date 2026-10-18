"""Step 7: result tables (Markdown + LaTeX) and error budgets -> results/tables.md, results/tables.tex,
results/error_budget.json.  Main configuration = F2 (forward points as data, correlated nuisances) + EFT
truncation covariance (scheme B); F1 (Baldin penalty) is the comparison variant."""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
from common import RES, save_json, PROFILE, MAIN_DROP, BAL_TXT, G0_TXT
V2 = (PROFILE == 'v2')
FWDE = 145 if V2 else 135
import database as DB

FL = json.load(open(os.path.join(RES, 'fits_L.json')))
FD = json.load(open(os.path.join(RES, 'delta_region.json')))
MD, TX = [], []


import re
TEXMAP = [('**', ''), ('γM1M1', '$\\gamma_{M1M1}$'), ('γE1E1', '$\\gamma_{E1E1}$'), ('γE1M2', '$\\gamma_{E1M2}$'),
          ('γM1E2', '$\\gamma_{M1E2}$'), ('γπ', '$\\gamma_\\pi$'), ('γ0', '$\\gamma_0$'), ('γNΔ', '$\\gamma N\\Delta$'),
          ('HIγS', 'HI$\\gamma$S'), ('Σ2z(0)', '$\\Sigma_{2z}(0^\\circ)$'), ('Σ2x', '$\\Sigma_{2x}$'), ('Σ2z', '$\\Sigma_{2z}$'),
          ('Σ3', '$\\Sigma_3$'), ('dσ(0)', 'd$\\sigma(0^\\circ)$'), ('dσ', 'd$\\sigma$'), ('δα,δβ', '$\\delta\\alpha,\\delta\\beta$'),
          ('λ_f', '$\\lambda_f$'), ('λ_g', '$\\lambda_g$'), ('Λ_b', '$\\Lambda_b$'), ('Λb', '$\\Lambda_b$'), ('c̄', '$\\bar c$'),
          ('ℓE', '$\\ell_E$'), ('ℓθ', '$\\ell_\\theta$'), ('ℓ', '$\\ell$'), ('χ²', '$\\chi^2$'), ('mπ', '$m_\\pi$'),
          ('fπ', '$f_\\pi$'), ('√', '$\\surd$'), ('πΔ', '$\\pi\\Delta$'), ('Δ', '$\\Delta$'), ('δ', '$\\delta$'),
          ('±', '$\\pm$'), ('≤', '$\\le$'), ('→', '$\\to$'), ('−', '$-$'), ('½', '1/2'), ('°', '$^\\circ$'), ('π', '$\\pi$'),
          ('α', '$\\alpha$'), ('β', '$\\beta$'), ('γ', '$\\gamma$'), ('σ', '$\\sigma$'), ('|', '$|$'), ('_', '\\_'),
          ('%', '\\%'), ('ü', '\\"u'), ('é', "\\'e")]


def texify(x):
    out = ''
    i = 0
    while i < len(x):
        for a, b in TEXMAP:
            if x.startswith(a, i):
                out += b; i += len(a); break
        else:
            out += x[i]; i += 1
    return out.replace('$$', '')


def md(*l):
    MD.extend(l)


def tx(*l):
    TX.extend(l)


def pe(v, e, nd=2):
    return ('%.' + str(nd) + 'f(%d)') % (v, round(e * 10 ** nd)) if e < 1 else ('%.' + str(nd) + 'f(%.' + str(nd) + 'f)') % (v, e)


def complete(r):
    C = np.array(r['covP'])
    for k, J in [('apb', [1, 1, 0, 0, 0, 0]), ('amb', [1, -1, 0, 0, 0, 0]), ('g0', [0, 0, -1, -1, -1, -1]), ('gpi', [0, 0, -1, 1, -1, 1])]:
        J = np.array(J, float)
        r[k] = float(J @ np.array(r['P']))
        r['d' + k] = float(np.sqrt(max(J @ C @ J, 0)))
    return r


def row(r):
    r = complete(r)
    P, d = r['P'], r['dP']
    return dict(N=r['N'], a=pe(P[0], d[0]), b=pe(P[1], d[1]), apb=pe(r.get('apb', P[0] + P[1]), r.get('dapb', 0)),
                amb=pe(r.get('amb', P[0] - P[1]), r.get('damb', 0)), chi=('%.1f/%d' % (r['chi2'], r['dof'])),
                chin='%.2f' % (r['chi2'] / r['dof']), p='%.2f' % r['p'])


def table(title, names, labels=None, cols=('N', 'a', 'b', 'apb', 'amb', 'chi', 'p'), ref=None, extra=None, src=FL):
    hdr = {'N': 'N', 'a': 'αE1', 'b': 'βM1', 'apb': 'α+β', 'amb': 'α−β', 'chi': 'χ²/dof', 'p': 'p', 'da': 'Δα', 'db': 'Δβ'}
    thdr = {'gMM': '$\\gamma_{M1M1}$', 'gEE': '$\\gamma_{E1E1}$', 'g0': '$\\gamma_0$', 'gpi': '$\\gamma_\\pi$', 'N': '$N$', 'a': '$\\alpha_{E1}$', 'b': '$\\beta_{M1}$', 'apb': '$\\alpha+\\beta$', 'amb': '$\\alpha-\\beta$',
            'chi': '$\\chi^2/{\\rm dof}$', 'p': '$p$', 'da': '$\\Delta\\alpha$', 'db': '$\\Delta\\beta$'}
    cols = list(cols) + (['da', 'db'] if ref else []) + (list(extra[0]) if extra else [])
    md('', '### ' + title, '', '| fit | ' + ' | '.join(hdr.get(c, c) for c in cols) + ' |', '|---|' + '---|' * len(cols))
    tx('', '% ' + title, '\\begin{tabular}{l' + 'c' * len(cols) + '}', '\\hline',
       'Fit & ' + ' & '.join(thdr.get(c, c) for c in cols) + ' \\\\', '\\hline')
    for i, n in enumerate(names):
        if n not in src:
            continue
        r = src[n]
        d = row(r)
        if ref:
            d['da'] = '%+.2f' % (r['P'][0] - src[ref]['P'][0]); d['db'] = '%+.2f' % (r['P'][1] - src[ref]['P'][1])
        if extra:
            d.update(extra[1](r))
        lab = labels[i] if labels else n
        md('| %s | ' % lab + ' | '.join(str(d[c]) for c in cols) + ' |')
        tx(texify(lab) + ' & ' + ' & '.join(str(d[c]).replace('*', '$^*$') for c in cols) + ' \\\\')
    tx('\\hline', '\\end{tabular}')


# ------------------------------------------------------------------ main results
spin = (['gMM', 'gEE', 'g0', 'gpi'], lambda r: dict(gMM=pe(r['P'][3], r['dP'][3]) if r['dP'][3] > 0 else '%.2f*' % r['P'][3],
                                                    gEE=pe(r['P'][2], r['dP'][2]) if r['dP'][2] > 0 else '%.2f*' % r['P'][2],
                                                    g0=pe(r['g0'], r.get('dg0', 0)) if r.get('dg0', 0) > 0 else '%.2f*' % r['g0'],
                                                    gpi=pe(r['gpi'], r.get('dgpi', 0)) if r.get('dgpi', 0) > 0 else '%.2f*' % r['gpi']))
md('# Fit results (generated by scripts/s07_tables.py; code profile %s)' % PROFILE, '',
   'Units: α, β in 1e-4 fm³; γ in 1e-4 fm⁴.  Errors are the full fit errors (statistical ⊕ floating normalisations ⊕',
   'forward-point nuisances/Baldin ⊕ truncation covariance where included).  "*" = fixed input (NNLO value).',
   'Main configuration: F2 = MAIN database (pre-2015 world data + A2 2022, E ≤ 150 MeV, minus the 3σ outlier(s) MAIN_DROP = %s;' % (list(MAIN_DROP),),
   'HIγS 2022 EXCLUDED since 2026-10-17, author decision) + forward points (E ≤ %d MeV, correlated nuisances), no Baldin penalty,' % FWDE,
   'spin polarisabilities at NNLO, EFT truncation covariance (scheme B, LEC-matched).  "d_all" = main database; "d_allpts" = all',
   'points incl. the outlier; "j_withHIGS" = main database + HIγS 2022 (the former main database).')
table('Main fits (region L, E ≤ 150 MeV)',
      ['F2_d_all_thB', 'F2_d_allpts_thB', 'F2_d_all', 'F2_d_fwd_uncorr_thB', 'F2_d_fwd_uncorr', 'F2pi_d_all_thB', 'F3_d_all_thB', 'F3_d_all',
       'F3g0_d_all_thB', 'F1_d_all_thB', 'F1_d_all', 'F3b_d_all_thB', 'F3c_d_all_thB'],
      ['**MAIN**: F2 + trunc. cov., main database (no HIγS, 3σ outlier removed)', 'F2 + trunc. cov., all points incl. the 3σ outlier', 'F2, no trunc. cov.',
       'F2 fwd. uncorrelated + trunc. cov.', 'F2 fwd. uncorrelated, no trunc. cov.',
       'F2 + γπ free (fixed γ0) + trunc.', 'F3: F2 + γM1M1, γE1E1 free + trunc.', 'F3, no trunc. cov.',
       'F3g0: fwd CS(0) + γ0 = %s, γM1M1, γE1E1 free + trunc.' % G0_TXT, 'F1 (Baldin %s) + trunc. cov.' % BAL_TXT,
       'F1, no trunc. cov.', 'F3b: F1 + γπ free (fixed γ0) + trunc.', 'F3c: F1 + γM1M1 free + trunc.'],
      extra=spin)

# ------------------------------------------------------------------ database variants
VN = [('a_pre2015', '(a) pre-2015 only'), ('b_A2only', '(b) A2 2022 only'), ('c_HIGSonly', '(c) HIγS 2022 only (test; not in the main database)'),
      ('d_all', '(d) MAIN database (pre-2015 + A2 2022, 3σ outlier removed, no HIγS)'), ('d_allpts', "(d') main database incl. the 3σ outlier"),
      ('e_noTAPS', '(e) main − TAPS'), ('f_GMPFpruned', '(f) GMPF/McGovern pruning'),
      ('g_A2perbeam', '(g) A2 per beamtime'), ('h_3sigma', '(h) 3σ rule on the main database (check)'),
      ('j_withHIGS', '(j) main + HIγS 2022 (the former main database)')]
for kind, lab in [('F2', 'main configuration (F2)'), ('F1', 'F1 (Baldin)')]:
    names = []; labs = []
    for v, l in VN:
        nm = '%s_%s_thB' % (kind, v)
        if v == 'h_3sigma':
            nm = '%s_h_3sigma_thBmain' % kind
        names.append(nm); labs.append(l)
    names += ['%s_F4_all_thB' % kind, '%s_F4_pruned_thB' % kind, '%s_F4_all_Hallin170_thB' % kind]
    labs += ['F4: E ≤ 170 MeV (Hallin ≤ 150)', 'F4 + GMPF pruning', 'F4 with Hallin ≤ 170']
    table('Database variants, %s, with truncation covariance' % lab, names, labs, ref='%s_d_all_thB' % kind)
    names = ['%s_%s' % (kind, v) for v, _ in VN] + ['%s_F4_all' % kind, '%s_F4_pruned' % kind, '%s_F4_all_Hallin170' % kind]
    table('Database variants, %s, without truncation covariance' % lab, names, labs, ref='%s_d_all' % kind)

# ------------------------------------------------------------------ normalisations, per-set chi2 of main fit
for key in ['F2_d_all_thB', 'F1_d_all_thB']:
    r = FL[key]
    md('', '### Per-data-set χ² and normalisations: %s' % key, '', '| set | N | χ² | χ²/N | n_k | δ_k (prior) | pull (n−1)/δ |', '|---|---|---|---|---|---|---|')
    tx('', '% per-set ' + key, '\\begin{tabular}{lrrrcc}', '\\hline', 'Data set & $N$ & $\\chi^2$ & $\\chi^2/N$ & $n_k$ & $\\delta_k$ \\\\', '\\hline')
    db = {}
    import pandas as pd
    dfdb = pd.read_csv(os.path.join(os.path.dirname(RES), 'db', 'fitdb.csv'), keep_default_na=False)
    for ds, v in r['perset'].items():
        dl = float(dfdb[dfdb.dataset == ds]['norm'].iloc[0])
        pull = (v['n'] - 1) / dl if dl > 0 else 0
        md('| %s | %d | %.1f | %.2f | %s | %s | %s |' % (DB.EXPT.get(ds, ds), v['N'], v['chi2'], v['chi2'] / v['N'],
           pe(v['n'], v['dn'], 3) if dl > 0 else '1 (fixed)', '%.1f%%' % (100 * dl) if dl > 0 else '–', '%+.2f' % pull if dl > 0 else '–'))
        tx('%s & %d & %.1f & %.2f & %s & %s \\\\' % (ds.replace('_', '\\_'), v['N'], v['chi2'], v['chi2'] / v['N'],
           pe(v['n'], v['dn'], 3) if dl > 0 else '1', ('%.1f\\%%' % (100 * dl)) if dl > 0 else '--'))
    tx('\\hline', '\\end{tabular}')
    if len(r.get('lam_shift', [])):
        md('', 'Forward-amplitude nuisance pulls: λ_f = %.2f ± %.2f, λ_g = %.2f ± %.2f (prior N(0,1)); χ²_data = %.1f, χ²_norm = %.1f, χ²_nuis = %.1f.'
           % (r['lam_shift'][0], r['lam_shift_err'][0], r['lam_shift'][1], r['lam_shift_err'][1], r['chi2_data'], r['chi2_norm'], r['chi2_shift']))

# ------------------------------------------------------------------ error budgets
def budget(kind):
    m = FL['%s_d_all_thB' % kind]; n = FL['%s_d_all' % kind]
    out = {}
    tot = np.array(m['dP'][:2])
    stat = np.array(n['dP_stat'][:2])                         # experimental statistics (+p2p syst), no theory cov.
    statth = np.array(m['dP_stat'][:2])                       # stat (+) truncation (nuisances fixed)
    trunc = np.sqrt(np.clip(statth ** 2 - stat ** 2, 0, None))
    norm = np.array(m['dP_norm'][:2])
    con = np.array(m['dP_fwd'][:2]) if kind == 'F2' else np.array(m['dP_baldin'][:2])
    # truncation-model variants (spread = rms deviation from the main fit)
    tv = ['L600', 'L900', 'fixed048', 'max', 'l50', 'l12', 'l100']
    dv = np.array([np.array(FL['%s_d_all_thB_%s' % (kind, v)]['P'][:2]) - np.array(m['P'][:2]) for v in tv])
    tmodel = np.sqrt(np.mean(dv ** 2, axis=0))
    # database selection: half range over d (main), d' (all points incl. the outlier), e, f, g, h, F4
    dbv = ['%s_d_all_thB', '%s_d_allpts_thB', '%s_e_noTAPS_thB', '%s_f_GMPFpruned_thB', '%s_g_A2perbeam_thB', '%s_h_3sigma_thBmain',
           '%s_F4_all_thB']
    Pv = np.array([FL[x % kind]['P'][:2] for x in dbv])
    dbs = 0.5 * (Pv.max(0) - Pv.min(0))
    # spin polarisabilities: freeing gamma_pi at fixed gamma_0 (data-driven, gamma_pi = 7.3(1.0))
    sp = 'F2pi_d_all_thB' if kind == 'F2' else 'F3b_d_all_thB'
    spinp = np.sqrt(np.clip(np.array(FL[sp]['dP'][:2]) ** 2 - tot ** 2, 0, None))
    # gM +- 0.1
    gm = 0.5 * np.abs(np.array(FL['%s_d_all_gM3.07_thB' % kind]['P'][:2]) - np.array(FL['%s_d_all_gM2.87_thB' % kind]['P'][:2]))
    const = np.abs(np.array(FL['%s_d_all_PDGconst_thB' % kind]['P'][:2]) - np.array(m['P'][:2]))   # v2: consts=PDG24
    # hA: piDelta loops scale as hA^2; hA uncertainty ~1% (Gamma_Delta 114-120 MeV) -> 2% of the full non-static piDelta effect
    if ('%s_d_all_hA2.84_thB' % kind) in FL:     # v2: hA = 2.84-2.91 (Gamma_Delta = 114-120 MeV)
        hA = 0.5 * np.abs(np.array(FL['%s_d_all_hA2.91_thB' % kind]['P'][:2]) - np.array(FL['%s_d_all_hA2.84_thB' % kind]['P'][:2]))
    else:
        hA = 0.02 * np.abs(np.array(FL['%s_d_all_NNLO_noPiD_thB' % kind]['P'][:2]) - np.array(m['P'][:2]))
    fwdmodel = np.abs(np.array(FL['F2_d_fwd_uncorr_thB']['P'][:2]) - np.array(m['P'][:2])) / 2 if kind == 'F2' else np.zeros(2)
    comps = [('stat. (exp., point-to-point)', stat), ('normalisations', norm),
             ('forward points (nuisances λ_f, λ_g)' if kind == 'F2' else 'Baldin sum rule ' + BAL_TXT, con),
             ('EFT truncation (covariance in χ²)', trunc), ('= total fit error', tot),
             ('truncation-model choice (Λ_b, δ, c̄, ℓ)', tmodel), ("database selection (half range d,d',e,f,g,h,F4)", dbs),
             ('spin polarisabilities (γπ free: ' + ('%.1f(%.1f)' % (FL[sp]['gpi'], FL[sp]['dgpi'])) + ')', spinp),
             ('γNΔ coupling gM ± 0.1', gm),
             ('constants → PDG24 (gA, fπ, masses, hA refit)' if V2 else 'constants gA, fπ → PDG', const),
             ('hA ± (Γ_Δ = 114–120 MeV: hA = 2.84–2.91)' if V2 else 'hA (Δ width), est.', hA)]
    if kind == 'F2':
        comps.append(('forward-error model (½|corr−uncorr|)', fwdmodel))
    extra = [c for c in comps[5:]]          # everything beyond the total fit error
    grand = np.sqrt(tot ** 2 + sum(v ** 2 for _, v in extra))
    comps.append(('**TOTAL**', grand))
    # alpha-beta and alpha+beta
    out = dict(P=m['P'][:2], comps=[(k, list(v)) for k, v in comps], total=list(grand), apb=m['apb'], amb=m['amb'])
    return out


EB = {}
for kind in ['F2', 'F1']:
    b = budget(kind); EB[kind] = b
    title = 'main fit (F2: forward points)' if kind == 'F2' else 'variant F1 (Baldin)'
    md('', '### Error budget: %s' % title, '', 'Central values: αE1 = %.2f, βM1 = %.2f' % tuple(b['P']), '',
       '| source | δαE1 | δβM1 |', '|---|---|---|')
    tx('', '% error budget ' + kind, '\\begin{tabular}{lcc}', '\\hline', 'Source & $\\delta\\alpha_{E1}$ & $\\delta\\beta_{M1}$ \\\\', '\\hline')
    for k, v in b['comps']:
        md('| %s | %.2f | %.2f |' % (k, v[0], v[1]))
        tx('%s & %.2f & %.2f \\\\' % (texify(k), v[0], v[1]))
    tx('\\hline', '\\end{tabular}')
save_json(EB, 'error_budget.json')

# ------------------------------------------------------------------ HIGS test
HT = FL['HIGS_test']
md('', '### HIγS 2022 test (HIγS = Li:2022vnz, CIRC + LIN sets, 16 points; NOT part of the main database since 2026-10-17)', '',
   'main = main database (no HIγS); +HIγS = variant (j) main + HIγS (the former main fit); Δ = shift of the result when the HIγS data are ADDED;',
   'χ²_PG = χ²(+HIγS) − χ²(HIγS alone) − χ²(main), dof 2; χ²(HIγS pred.) = the 16 HIγS points predicted from the MAIN fit, at n = 1 and',
   'with the two normalisations fitted within their 2.24% prior (penalty included).', '',
   '| configuration | α (main) | β (main) | α−β (main) | α (+HIγS) | β (+HIγS) | Δα | Δβ | χ²/dof (main) | χ²/dof (+HIγS) | χ²(HIγS) in +HIγS fit | χ²_PG (dof 2) | p_PG | χ²(HIγS pred.) n=1 / n fitted |',
   '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for k in ['F2_thB', 'F2', 'F1_thB', 'F1', 'F3_thB', 'F3']:
    h = HT[k]
    md('| %s | %s | %s | %.2f | %s | %s | %+.2f | %+.2f | %.1f/%d | %.1f/%d | %.1f | %.2f | %.3f | %.1f / %.1f |' % (
        k.replace('_thB', ' + trunc. cov.'), pe(h['P_main'][0], h['dP_main'][0]), pe(h['P_main'][1], h['dP_main'][1]),
        h['P_main'][0] - h['P_main'][1], pe(h['P_withHIGS'][0], h['dP_withHIGS'][0]), pe(h['P_withHIGS'][1], h['dP_withHIGS'][1]),
        h['shift'][0], h['shift'][1], h['chi2_main'], h['dof_main'], h['chi2_withHIGS'], h['dof_withHIGS'],
        sum(h['chi2_HIGS_in_withHIGS'].values()), h['PG'], h['PG_p'], h['chi2_HIGS_pred_n1'], h['chi2_HIGS_pred_nfit']))
for k in ['F1_fixsum', 'F1_fixsum_thB']:
    h = HT[k]
    md('', ('%s (α+β = ' + BAL_TXT.split('(')[0] + ' fixed in all three fits; one shared parameter α−β): χ²_PG = %.2f (dof 1, p = %.3f); α−β: +HIγS %.2f, main %.2f, HIγS alone %.2f')
       % (k, h['PG'], h['PG_p'], h['amb_withHIGS'], h['amb_main'], h['amb_HIGS']))
md('', 'HIγS alone (α, β free, no sum constraint, + trunc. cov.): α = %s, β = %s, χ² = %.1f/%d (cf. Li et al. HBχPT: 15.4(1.8), 2.1(2.0)); without trunc. cov.: α = %s, β = %s, χ² = %.1f/%d.'
   % (pe(HT['F2_thB']['P_HIGSfree'][0], HT['F2_thB']['dP_HIGSfree'][0]), pe(HT['F2_thB']['P_HIGSfree'][1], HT['F2_thB']['dP_HIGSfree'][1]),
      HT['F2_thB']['chi2_HIGSfree'], HT['F2_thB']['dof_HIGSfree'],
      pe(HT['F2']['P_HIGSfree'][0], HT['F2']['dP_HIGSfree'][0]), pe(HT['F2']['P_HIGSfree'][1], HT['F2']['dP_HIGSfree'][1]),
      HT['F2']['chi2_HIGSfree'], HT['F2']['dof_HIGSfree']))
h = HT['F2_thB']
md('', 'HIγS normalisations in the +HIγS fit (F2 + trunc. cov.): ' + ', '.join('%s %s' % (k, pe(v[0], v[1], 3)) for k, v in h['n_HIGS_in_withHIGS'].items()))
md('', '#### Pulls of the 16 HIγS points predicted by the MAIN fit (F2 + trunc. cov., no HIγS in the fit)', '',
   'n_HIγS fitted (with 2.24%% prior, P fixed): CIRC %.3f, LIN %.3f' % tuple(h['n_HIGS_fit'].values()), '',
   '| point | obs | E | θ_lab | data | prediction (n=1) | pull (n=1) | pull (n fitted) |', '|---|---|---|---|---|---|---|---|')
for i, p in enumerate(h['HIGS_pids']):
    md('| %s | %s | %.1f | %.0f | %.2f ± %.2f | %.2f | %+.2f | %+.2f |' % (p, h['HIGS_obs'][i], h['HIGS_E'][i], h['HIGS_th'][i], h['HIGS_y'][i],
       h['HIGS_err'][i], h['HIGS_pred'][i], h['pull_n1'][i], h['pull_nfit'][i]))

# ------------------------------------------------------------------ theory-error variants
table('Truncation-error model variants (main configuration F2)',
      ['F2_d_all_thB', 'F2_d_all_thB_L600', 'F2_d_all_thB_L900', 'F2_d_all_thB_fixed048', 'F2_d_all_thB_max', 'F2_d_all_thB_l50',
       'F2_d_all_thB_l12', 'F2_d_all_thB_l100', 'F2_d_all_thB_diag', 'F2_d_all_thA', 'F2_d_all_thA_diag', 'F2_d_all'],
      ['default: scheme B, Λb = 750 MeV, p = max(mπ,E)/Λb, c̄ rms, ℓE = 25 MeV, ℓθ = 40°', 'Λb = 600 MeV', 'Λb = 900 MeV',
       'δ = mπ/Δ = 0.48 fixed', 'c̄ = max(|c2|,|c3|)', 'ℓE = 50 MeV, ℓθ = 50°', 'ℓE = 12 MeV, ℓθ = 25°', 'ℓE = 100 MeV, ℓθ = 100°',
       'uncorrelated (diagonal)', 'scheme A (naive, double counts LECs)', 'scheme A, diagonal', 'no truncation covariance'],
      ref='F2_d_all_thB')
table('Truncation-error model variants (F1)',
      ['F1_d_all_thB', 'F1_d_all_thB_L600', 'F1_d_all_thB_L900', 'F1_d_all_thB_fixed048', 'F1_d_all_thB_max', 'F1_d_all_thB_l50',
       'F1_d_all_thB_l12', 'F1_d_all_thB_l100', 'F1_d_all_thB_diag', 'F1_d_all_thA', 'F1_d_all_thA_diag', 'F1_d_all'],
      ['default', 'Λb = 600 MeV', 'Λb = 900 MeV', 'δ = 0.48 fixed', 'c̄ = max', 'ℓ = 50/50', 'ℓ = 12/25', 'ℓ = 100/100',
       'diagonal', 'scheme A', 'scheme A diagonal', 'no truncation covariance'], ref='F1_d_all_thB')

# ------------------------------------------------------------------ order-by-order, parameters, spin inputs
for kind in ['F2', 'F1']:
    names = ['%s_d_all_thB', '%s_d_all_NLO_thB', '%s_d_all_NLO_spinNNLO_thB', '%s_d_all_NNLO_noPiD_thB', '%s_d_all_gM2.87_thB',
             '%s_d_all_gM3.07_thB', '%s_d_all_gA1.2754_thB', '%s_d_all_fpi92.07_thB', '%s_d_all_PDGconst_thB',
             '%s_d_all_hA2.84_thB', '%s_d_all_hA2.91_thB', '%s_d_all_singlemass_thB',
             '%s_d_all_spinPaudyal_thB', '%s_d_all_spinMRPP22_thB', '%s_d_all', '%s_d_all_NLO', '%s_d_all_NLO_spinNNLO',
             '%s_d_all_NNLO_noPiD', '%s_d_all_gM2.87', '%s_d_all_gM3.07', '%s_d_all_PDGconst', '%s_d_all_hA2.84', '%s_d_all_hA2.91',
             '%s_d_all_singlemass', '%s_d_all_spinPaudyal', '%s_d_all_spinMRPP22']
    PDGL = 'PDG24 constants (consts=PDG24, hA = 2.873)' if V2 else 'gA, fπ PDG'
    labs = ['NNLO (reference) + trunc.', 'NLO amplitude + δα,δβ + trunc.', 'NLO + δα,δβ, spin pols matched to NNLO + trunc.',
            'NNLO without πΔ loops (static pols matched) + trunc.', 'gM = 2.87 + trunc.', 'gM = 3.07 + trunc.', 'gA = 1.2754 + trunc.',
            'fπ = 92.07 MeV + trunc.', PDGL + ' + trunc.', 'hA = 2.84 + trunc.', 'hA = 2.91 + trunc.',
            'single (charged) pion mass in loops (as in the 2015 code) + trunc.', 'spin pols = Paudyal 2020 + trunc.', 'spin pols = MRPP 2022 + trunc.',
            'NNLO (no trunc.)', 'NLO + δα,δβ', 'NLO, spin matched', 'NNLO w/o πΔ loops', 'gM = 2.87', 'gM = 3.07', PDGL,
            'hA = 2.84', 'hA = 2.91', 'single pion mass', 'spin = Paudyal 2020', 'spin = MRPP 2022']
    table('Order-by-order fits, parameter and spin-polarisability variations (%s)' % kind, [n % kind for n in names], labs,
          ref='%s_d_all_thB' % kind)

# ------------------------------------------------------------------ forward model, benchmarks, single observables
table('Forward-point treatment (main configuration)',
      ['F2_d_all_thB', 'F2_d_fwd_uncorr_thB', 'F2_d_fwdonly_CS_thB', 'F2_d_fwd135_thB', 'F2_d_fwd135_uncorr_thB',
       'F2_d_fwd150_thB', 'F2_d_fwd150_uncorr_thB', 'F2_d_all', 'F2_d_fwd_uncorr', 'F2_d_fwdonly_CS', 'F2_d_fwd135', 'F2_d_fwd_150'],
      ['nuisance model, E ≤ %d MeV, + trunc. (MAIN)' % FWDE, 'uncorrelated, E ≤ %d, + trunc.' % FWDE, 'dσ(0) only (no Σ2z(0)), + trunc.',
       'nuisance, E ≤ 135 MeV (without 145, 150), + trunc.', 'uncorrelated, E ≤ 135, + trunc.',
       'nuisance, E ≤ 150 MeV, + trunc.', 'uncorrelated, E ≤ 150, + trunc.', 'nuisance, no trunc.',
       'uncorrelated, no trunc.', 'dσ(0) only, no trunc.', 'nuisance, E ≤ 135, no trunc.', 'nuisance, E ≤ 150, no trunc.'], ref='F2_d_all_thB')
table('Benchmarks and single-observable / subset fits',
      ['bench_Mornacchi', 'bench_Mornacchi_NNLOspin', 'LM14_type', 'LM14_type_B14', 'F1_A2_CSonly', 'F1_A2_S3only', 'F1_TAPSonly',
       'F1_pre2015+HIGS', 'F1_pre2015+A2', 'F1_A2+HIGS', 'F1_d_free', 'F1_A2_CSonly_thB', 'F1_A2_S3only_thB', 'F1_TAPSonly_thB',
       'F1_pre2015+HIGS_thB', 'F1_pre2015+A2_thB', 'F1_A2+HIGS_thB', 'F1_d_free_thB'],
      ['Mornacchi benchmark: A2 per beamtime, Baldin 13.8(4), spin = Paudyal [thesis: 10.66(15)(41), 3.40(20)(22), 154.6/187]',
       'same, spin pols at NNLO', 'LM14-type: pre-2015, GMPF pruning, E ≤ 170, Baldin 13.8(4) [LM14: 10.6(5), 3.2(5)]',
       'same with Baldin ' + BAL_TXT, 'A2 dσ only (Baldin)', 'A2 Σ3 only (Baldin)', 'TAPS only (Baldin)', 'pre-2015 + HIγS (Baldin)',
       'pre-2015 + A2 (Baldin)', 'A2 + HIγS (Baldin)', 'all, no sum constraint', 'A2 dσ only + trunc.', 'A2 Σ3 only + trunc.',
       'TAPS only + trunc.', 'pre-2015 + HIγS + trunc.', 'pre-2015 + A2 + trunc.', 'A2 + HIγS + trunc.', 'all, no sum constraint + trunc.'])

# ------------------------------------------------------------------ outliers and energy scan
md('', '### Outlier analysis (3σ rule, two passes)', '')
if 'MAIN_DROP' in FL:
    r = FL['MAIN_DROP']
    md('- MAIN_DROP (held fixed; removed from the main database): %s; re-derived from the all-points fit %s: %s.  Points with |pull| > 2.5 (pass 1): %s'
       % (r['pids'], r['determined_from'], r['removed'], ', '.join('%s (%+.2f)' % (a, b) for a, b in r['pulls_gt2_5'][0].items()) or 'none'))
for k in ['F2_h_3sigma_thBmain', 'F1_h_3sigma_thBmain', 'F2_h_3sigma', 'F1_h_3sigma']:
    r = FL[k]
    rem = r.get('removed') or []
    md('- %s (search on the main database): removed %s; %s α = %s, β = %s, χ²/dof = %.1f/%d.%s' % (
       k, rem if rem else 'nothing', 'then' if rem else 'fit unchanged (= d_all):', pe(r['P'][0], r['dP'][0]), pe(r['P'][1], r['dP'][1]),
       r['chi2'], r['dof'], ('  Points with |pull| > 2.5 (pass 1): ' + (', '.join('%s (%+.2f)' % (a, b) for a, b in r['pulls_gt2.5'][0].items()) or 'none'))
       if r.get('pulls_gt2.5') else ''))
md('', '### Upper-energy-cut scan (main database, no HIγS; Hallin restricted to ≤ 150 MeV)', '', '| E_max | F2 + trunc.: α | β | F2 no trunc.: α | β | F1 + trunc.: α | β | F1 no trunc.: α | β |', '|---|---|---|---|---|---|---|---|---|')
for E in [90, 100, 110, 120, 130, 140, 150, 160, 170]:
    cells = []
    for k in ['F2_thB', 'F2_noth', 'F1_thB', 'F1_noth']:
        r = FL['scan_%s_%s_E%d' % (k.split('_')[0], k.split('_')[1], E)]
        cells += [pe(r['P'][0], r['dP'][0]), pe(r['P'][1], r['dP'][1])]
    md('| %d | ' % E + ' | '.join(cells) + ' |')

# ------------------------------------------------------------------ Delta region
for lab in ['NNLO', 'NNLO_novtx', 'F2', 'F2_novtx', 'F1', 'F2_noth']:
    if 'chi2_' + lab not in FD:
        continue
    md('', '### Δ region (150 < E ≤ 350 MeV): χ² per data set, %s' % ({'NNLO': 'parameter-free NNLO prediction',
       'NNLO_novtx': 'parameter-free NNLO prediction WITHOUT γNΔ vertex corrections',
       'F2_novtx': 'α, β of the main fit, WITHOUT vertex corrections'}.get(lab, 'NNLO with α, β of %s (spin pols NNLO)' % lab)), '',
       '| set | obs | N | E range | χ² (n=1) | χ² (n fitted) | n | χ²/N |', '|---|---|---|---|---|---|---|---|')
    for r in FD['chi2_' + lab]:
        md('| %s | %s | %d | %s | %.1f | %.1f | %.3f | %.2f |' % (DB.EXPT.get(r['set'], r['set']), r['obs'], r['N'], r['Erange'],
           r['chi2_n1'], r['chi2_nfit'], r['n'], r['chi2N']))
md('', '### Δ region: exploratory spin-polarisability fits (α−β from main fit, α+β = %s, γ0 = %s)' % (BAL_TXT, G0_TXT), '',
   '| fit | data | N | γE1E1 | γM1M1 | γE1M2 | γM1E2 | γ0 | γπ | gM | χ²/dof |', '|---|---|---|---|---|---|---|---|---|---|---|')
DESC = {'S_A2_2free': ('Σ2x+Σ2z', 'γE1E1, γM1M1 free'), 'S_A2': ('Σ2x+Σ2z', '4 free (under-determined)'),
        'S_A2_gpi': ('Σ2x+Σ2z, γπ = 8.0(1.8)', '4 free'), 'S_A2_LEGS': ('Σ2x+Σ2z+LEGS Σ3 (< 2π thr.)', '4 free'),
        'S_A2_LEGS_gpi': ('Σ2x+Σ2z+LEGS Σ3, γπ = 8.0(1.8)', '4 free'), 'S_A2_LEGS_2free': ('Σ2x+Σ2z+LEGS Σ3', 'γE1E1, γM1M1 free')}
for base in ['S_A2_2free', 'S_A2', 'S_A2_gpi', 'S_A2_LEGS_2free', 'S_A2_LEGS', 'S_A2_LEGS_gpi']:
    for sfx, sl in [('', ''), ('_vtx', ' + vertex th. err.'), ('_gM', ' + gM float (prior 2.97(10))'), ('_gMwide', ' + gM float (prior 2.97(30))'),
                    ('_novtx', ' WITHOUT vertex corr.'), ('_gMwide_novtx', ' + gM float (2.97(30)), WITHOUT vertex corr.')]:
        k = base + sfx
        if k not in FD:
            continue
        r = FD[k]
        P, d = r['P'], r['dP']
        f = lambda i: pe(P[i], d[i]) if d[i] > 0 else '%.2f*' % P[i]
        gm = pe(*r['extra']['gM'], 3) if r['extra'] else '2.97*'
        md('| %s%s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %.1f/%d |' % (k, sl, DESC[base][0] + '; ' + DESC[base][1], r['N'], f(2), f(3), f(4), f(5),
           pe(r['g0'], r['dg0']), pe(r['gpi'], r['dgpi']), gm, r['chi2'], r['dof']))
md('', '| γNΔ coupling from Δ-region data alone (spin pols NNLO; α, β as above; prior 2.97(50)) | gM | χ²/dof | gM without vertex corr. | χ²/dof |',
   '|---|---|---|---|---|')
for k, l in [('gM_only_A2asym', 'A2 Σ2x + Σ2z'), ('gM_only_LEGS_S3', 'LEGS Σ3 (< 2π thr.)'), ('gM_only_MAMIcs', 'MAMI dσ (Hünger, Wissmann, Camen)'),
             ('gM_only_LEGScs', 'LEGS dσ'), ('gM_only_Hallin', 'Hallin dσ (> 170 MeV)')]:
    r = FD[k]; rn = FD.get(k + '_novtx')
    md('| %s | %s | %.1f/%d | %s | %s |' % (l, pe(*r['extra']['gM'], 3), r['chi2'], r['dof'],
       pe(*rn['extra']['gM'], 3) if rn else '–', ('%.1f/%d' % (rn['chi2'], rn['dof'])) if rn else '–'))
md('', 'Literature: Martel 2015 (DR): γE1E1 = −3.5(1.2), γM1M1 = 3.16(85), γE1M2 = −0.7(1.2), γM1E2 = 1.99(29); Martel BχPT: −2.6(8), 2.7(5);',
   'Paudyal 2020 HDPV: −3.18(52), 2.98(43), −0.44(67), 1.58(43); Paudyal BχPT: −2.65(43), 2.43(42), −1.32(72), 2.47(42);',
   'MRPP 2022: −3.0(6)(4), 3.7(5)(1), −1.2(1.0)(3), 2.0(7)(4); NNLO prediction (LMP15): −3.3(8), 2.9(1.5), 0.2(2), 1.1(3).')

open(os.path.join(RES, 'tables.md'), 'w').write('\n'.join(MD) + '\n')
open(os.path.join(RES, 'tables.tex'), 'w').write('% generated by scripts/s07_tables.py\n' + '\n'.join(TX) + '\n')
print('\n'.join(MD[:60]))
