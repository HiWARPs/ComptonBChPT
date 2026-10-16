"""Markdown + LaTeX tables from results/plb_reproduction.json, results/update_fits.json, results/main_fits.json and
results/main_fits_legacynorm.json.

Two configurations are tabulated, both on the MAIN DATABASE (all data at E <= 150 MeV minus the 3-sigma outlier
TAPS 108.42/133 = update_fits.MAIN_DROP; 'db_allpts' keeps it):
  P-*  the PWA fit of the paper: all normalisations floated, incl. the pre-1961 sets (main_fits_legacynorm.json);
       its method uncertainty goes to results/main_total_errors.json;
  D/M/S/V/I/F/C  the fixed-legacy-normalisation configuration (main_fits.json; Sections 5-10 of pwa_report.md); in the
       paper only as the row "normalizations of the oldest data fixed"; errors in results/main_total_errors_fixedlegacy.json.
Writes results/tables.md, results/tables.tex."""
import json
import numpy as np
import pandas as pd
import reproduce_plb as RP
import update_fits as UF

P = json.load(open('results/plb_reproduction.json'))
R = json.load(open('results/update_fits.json'))
MF = json.load(open('results/main_fits.json'))
ML = json.load(open('results/main_fits_legacynorm.json'))
md, tex = [], []
DROP = ', '.join(f'{ds} {E:g}/{th:g}' for ds, E, th in UF.MAIN_DROP)
METHOD = ['w30', 'bchptctr', 'l2zero', 'gEMfree', 'binsB', 'binmode', 'A2beamtimes', 'fwd_uncorr']
QK = ('alpha_m_beta', 'gamma_pi', 'betaM1', 'alphaE1')


def f(v, e=None, nd=1):
    return f'{v:.{nd}f}' if e is None else f'{v:.{nd}f}({int(round(e * 10**nd))})'


def ftex(v, e=None, nd=1):
    return f'${v:.{nd}f}$' if e is None or e == 0 else f'${v:.{nd}f}\\pm{e:.{nd}f}$'


md.append(f'Main database: all data at E <= 150 MeV minus the 3-sigma outlier ({DROP}; BChPT pid OlmosdeLeon2001_029). '
          'All fits below (except the KLP18 reproduction and the rows marked "all points") are on the main database.\n')

# ---------------- T1: reproduction of KLP18 ----------------
md.append('### Table R1. Reproduction of KLP18 (PLB 782 (2018) 34), Tables II and III\n')
md.append('MAP = joint fit with the KLP18 Gaussian widths as priors (errors: stage-2 = residuals fixed, as in KLP18 / full = marginalised); '
          'MC = KLP18 Monte-Carlo + stage-2 procedure, mean and spread over seeds; PLB = published.\n')
md.append('| Fit | quantity | MAP | err stage-2 | err full | MC mean (spread) | PLB |')
md.append('|---|---|---|---|---|---|---|')
tex.append('% Table R1: reproduction of KLP18\n\\begin{tabular}{l|ccccccc|c}\n\\hline')
tex.append('Fit & $\\alpha_{E1}$ & $\\beta_{M1}$ & $\\gamma_{E1E1}$ & $\\gamma_{M1M1}$ & $\\gamma_{M1E2}$ & $\\alpha_{E1}-\\beta_{M1}$ & $\\gamma_\\pi$ & $\\chi^2/N$\\\\\\hline')
for name in ['Fit 0', 'Fit 1', 'Fit 1_3s(PLB 4 pts)', "Fit 1'", "Fit 1''", 'Fit 2', "Fit 2''"]:
    r = P[name]
    ref = RP.PLB.get(name, RP.PLB.get('Fit 1_3s'))
    for k in RP.KEYS:
        if k == 'gE1M2' and name != 'Fit 0':
            continue
        mc = ''
        if 'mc' in r:
            v = np.array([x['val'][k] for x in r['mc']])
            mc = f'{v.mean():.2f} ({v.std():.2f})'
        pr = ref.get(k)
        prs = f'{pr[0]} ± {pr[1]}' if pr else ''
        md.append(f'| {name} | {k} | {r["val"][k]:.2f} | {r["err_stage2"][k]:.2f} | {r["err_full"][k]:.2f} | {mc} | {prs} |')
    md.append(f'| {name} | chi2/N (N) | {r["chi2pt"]:.3f} ({r["npts"]}); chi2/138 = {r["chi2"]/138:.3f} | | | | {ref["chi2pt"]} |')
    v, e = r['val'], r['err_stage2']
    row = ' & '.join(ftex(v[k], e[k]) for k in ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gM1E2', 'alpha_m_beta', 'gamma_pi'])
    nm = name.replace('_3s(PLB 4 pts)', '$_{3\\sigma}$').replace("''", "$''$").replace("1'", "1$'$")
    tex.append(f'{nm} (this repr.) & {row} & ${r["chi2pt"]:.2f}$\\\\')
    prow = ' & '.join(ftex(*ref[k]) if k in ref else '--' for k in ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gM1E2', 'alpha_m_beta', 'gamma_pi'])
    tex.append(f'{nm} (KLP18) & {prow} & ${ref["chi2pt"]}$\\\\\\hline')
tex.append('\\end{tabular}\n')

d = pd.read_csv('results/db_update.csv')
MV = [('MAIN', 'MAIN: all data of the main database, floating norms, forward points as chi2 data (correlated nuisances)'),
      ('db_allpts', f'all points incl. the 3-sigma outlier ({DROP})'),
      ('db_3sigma', '3-sigma rule rerun on the main database (check; = MAIN when nothing is flagged)'),
      ('db_pruned', 'McGovern pruning (no Oxley, Bernardini, Baranov 150, TAPS 108.4/133)'),
      ('db_modern', 'modern sets only (TAPS, MacGibbon, Hallin, Federspiel, A2, HIGS)'),
      ('db_new_only', 'A2 2022 + HIGS 2022 only'),
      ('db_A2_only', 'A2 2022 only'),
      ('db_noHIGS', 'without HIGS 2022'),
      ('db_noTAPS', 'without TAPS 2001'),
      ('db_noA2', 'without A2 2022'),
      ('db_noLegacy', 'without pre-1961 legacy sets (Oxley, Goldansky, Bernardini, Pugh, Hyman)'),
      ('db_pruned_noHIGS', 'pruned and without HIGS'),
      ('db_old_only', 'pre-2015 data only'),
      ('db_old_only_pruned', 'pre-2015 data only, pruned')]
ITEMS = [('statistical (incl. norms, forward-point nuisances n_f, n_g, residual priors)', None),
         ('residual priors x3 wider', 'w30'),
         ('residual priors centred on BChPT', 'bchptctr'),
         ('l=2 multipoles: BChPT -> 0', 'l2zero'),
         ('gamma_E1M2 free (vs 0)', 'gEMfree'),
         ('energy binning re-centred on A2/HIGS energies', 'binsB'),
         ('energy mode: data at bin energy (KLP18) vs actual energy', 'binmode'),
         ('A2 per beam time (4 scales) vs combined', 'A2beamtimes'),
         ('forward-point errors uncorrelated (vs correlated nuisances)', 'fwd_uncorr')]
assert [t for _, t in ITEMS if t] == METHOD


def config_tables(M, pre, title, errfile):
    """Tables D (per-dataset chi2), M (main + database variants), S (budget), V (method variants), I (influence),
    F (forward points), C (PG tests) of one MAIN configuration M; table ids prefixed by `pre`."""
    main = M['MAIN']
    md.append(f'\n## {title}\n')
    # ---- D: database and per-dataset chi2
    md.append(f'\n### Table {pre}D. Database (E <= 150 MeV) and chi^2 per dataset in the MAIN fit ({title})\n')
    md.append('N = points in the main fit (the excluded point is not counted).\n')
    md.append('| dataset | texkey | obs | N | E [MeV] | theta_lab [deg] | norm [%] | chi2 / N (main) | norm factor (main) |')
    md.append('|---|---|---|---|---|---|---|---|---|')
    tex.append(f'% Table {pre}D: database ({title})\n\\begin{{tabular}}{{llrrrrr}}\n\\hline')
    tex.append('Data set & obs. & $N$ & $E_\\gamma$ [MeV] & $\\theta_{\\rm lab}$ [deg] & scale [\\%] & $\\chi^2$ (main fit)\\\\\\hline')
    dm = UF.select(d)
    for (ds, obs), g in d.groupby(['dataset', 'obs'], sort=False):
        gm = dm[(dm.dataset == ds) & (dm.obs == obs)]
        pdx = main['per_dataset'].get(ds, {})
        nf = main['norms'].get(g.norm_group.iloc[0], {})
        nfs = f'{nf["N"]:.3f}({int(round(nf["err"]*1000))})' if nf else '--'
        md.append(f'| {ds} | {g.texkey.iloc[0]} | {obs} | {len(gm)} | {g.E.min():.1f}-{g.E.max():.1f} | {g.theta_lab.min():.0f}-{g.theta_lab.max():.0f} | '
                  f'{100*g.norm_err.iloc[0]:.1f} | {pdx.get("chi2", np.nan):.1f}/{pdx.get("N", 0)} | {nfs} |')
        tex.append(f'{ds.replace("_", " ")} \\cite{{{g.texkey.iloc[0]}}} & {obs.replace("_", " ")} & {len(gm)} & {g.E.min():.0f}--{g.E.max():.0f} & '
                   f'{g.theta_lab.min():.0f}--{g.theta_lab.max():.0f} & {100*g.norm_err.iloc[0]:.1f} & {pdx.get("chi2", np.nan):.1f}\\\\')
    tex.append('\\hline\n\\end{tabular}\n')
    # ---- M: MAIN and database variants
    md.append(f'\n### Table {pre}M. MAIN fit and database variants ({title}; errors full = stat incl. norms, nuisances, residual priors)\n')
    md.append('| fit | description | N | chi2/N | alpha_E1 | beta_M1 | alpha-beta | gamma_pi | gamma_E1E1 | gamma_M1M1 | gamma_M1E2 | alpha+beta | gamma_0 |')
    md.append('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    tex.append(f'% Table {pre}M: MAIN and database variants ({title})\n\\begin{{tabular}}{{l r c cccc ccc cc}}\n\\hline')
    tex.append('Fit & $N$ & $\\chi^2/N$ & $\\alpha_{E1}$ & $\\beta_{M1}$ & $\\alpha_{E1}-\\beta_{M1}$ & $\\gamma_\\pi$ & $\\gamma_{E1E1}$ & $\\gamma_{M1M1}$ & $\\gamma_{M1E2}$ & $\\alpha_{E1}+\\beta_{M1}$ & $\\gamma_0$\\\\\\hline')
    for t, desc in MV:
        r = M[t]; v, e = r['val'], r['err']
        md.append(f'| {t} | {desc} | {r["npts"]} | {r["chi2pt"]:.2f} | {v["alphaE1"]:.2f}({e["alphaE1"]:.2f}) | {v["betaM1"]:.2f}({e["betaM1"]:.2f}) | '
                  f'{v["alpha_m_beta"]:.2f}({e["alpha_m_beta"]:.2f}) | {v["gamma_pi"]:.2f}({e["gamma_pi"]:.2f}) | {v["gE1E1"]:.2f}({e["gE1E1"]:.2f}) | '
                  f'{v["gM1M1"]:.2f}({e["gM1M1"]:.2f}) | {v["gM1E2"]:.2f}({e["gM1E2"]:.2f}) | {v["alpha_p_beta"]:.2f}({e["alpha_p_beta"]:.2f}) | '
                  f'{v["gamma0"]:.3f}({e["gamma0"]:.3f}) |')
        tex.append(f'{desc.replace("_", " ")} & {r["npts"]} & {r["chi2pt"]:.2f} & ' +
                   ' & '.join(ftex(v[k], e[k], 2 if k == 'gamma0' else 1) for k in ['alphaE1', 'betaM1', 'alpha_m_beta', 'gamma_pi', 'gE1E1', 'gM1M1', 'gM1E2', 'alpha_p_beta', 'gamma0']) + '\\\\')
    tex.append('\\hline\n\\end{tabular}\n')
    md.append(f'\n3-sigma rule rerun on the main database removes: {M["db_3sigma"].get("removed") or "nothing"}; '
              f'on all points it flags: {[list(x) for x in UF.MAIN_DROP]} (= MAIN_DROP).')
    # ---- S: uncertainty budget
    base = main['val']
    md.append(f'\n### Table {pre}S. Uncertainty budget of the MAIN fit ({title}; shift = variant - MAIN)\n')
    md.append('| source | d(alpha-beta) | d(gamma_pi) | d(beta_M1) | d(alpha_E1) |')
    md.append('|---|---|---|---|---|')
    tot = {k: 0.0 for k in QK}
    meth = {k: 0.0 for k in QK}
    for lab, t in ITEMS:
        if t is None:
            dv = {k: main['err'][k] for k in QK}
        else:
            dv = {k: M[t]['val'][k] - base[k] for k in QK}
            for k in QK:
                meth[k] += dv[k]**2
        for k in QK:
            tot[k] += dv[k]**2
        md.append(f'| {lab} | {dv["alpha_m_beta"]:+.2f} | {dv["gamma_pi"]:+.2f} | {dv["betaM1"]:+.2f} | {dv["alphaE1"]:+.2f} |')
    md.append(f'| method/theory total (quadrature) | {np.sqrt(meth["alpha_m_beta"]):.2f} | {np.sqrt(meth["gamma_pi"]):.2f} | {np.sqrt(meth["betaM1"]):.2f} | {np.sqrt(meth["alphaE1"]):.2f} |')
    md.append(f'| **total** | **{np.sqrt(tot["alpha_m_beta"]):.2f}** | **{np.sqrt(tot["gamma_pi"]):.2f}** | **{np.sqrt(tot["betaM1"]):.2f}** | **{np.sqrt(tot["alphaE1"]):.2f}** |')
    for lab, t in (('variant (not added): exact sum-rule constraints (previous main)', 'exact_constraints'),
                   ('variant (not added): no floating normalisations', 'no_float_norm'),
                   (f'database (not added): all points incl. the outlier ({DROP})', 'db_allpts'),
                   ('database (not added): McGovern pruning', 'db_pruned'),
                   ('database (not added): modern sets only', 'db_modern'),
                   ('database (not added): without HIGS', 'db_noHIGS'),
                   ('database (not added): without legacy sets', 'db_noLegacy')):
        md.append(f'| {lab} | {M[t]["val"]["alpha_m_beta"]-base["alpha_m_beta"]:+.2f} | {M[t]["val"]["gamma_pi"]-base["gamma_pi"]:+.2f} | '
                  f'{M[t]["val"]["betaM1"]-base["betaM1"]:+.2f} | {M[t]["val"]["alphaE1"]-base["alphaE1"]:+.2f} |')
    md.append(f'\nPosteriors of the sum-rule quantities (controlled by the forward-point nuisances): alpha+beta = '
              f'{base["alpha_p_beta"]:.2f} +- {main["err"]["alpha_p_beta"]:.2f}, gamma_0 = {base["gamma0"]:.3f} +- {main["err"]["gamma0"]:.3f}; '
              f'n_f = {main["nf"]:.2f}, n_g = {main["ng"]:.2f}; rho(alpha-beta, gamma_pi) = {main["corr_amb_gpi"]:.2f}.')
    SYS = {k: dict(central=float(base[k]), stat=float(main['err'][k]), method=float(np.sqrt(meth[k])), total=float(np.sqrt(tot[k])))
           for k in QK}
    SYS['alpha_p_beta'] = dict(central=float(base['alpha_p_beta']), stat=float(main['err']['alpha_p_beta']))
    SYS['gamma0'] = dict(central=float(base['gamma0']), stat=float(main['err']['gamma0']))
    SYS['_config'] = title
    SYS['_main_database'] = f'all data at E <= 150 MeV minus {DROP}; N = {main["npts"]}'
    SYS['_method_variants'] = METHOD
    json.dump(SYS, open(errfile, 'w'), indent=1)
    # ---- V: method variants
    md.append(f'\n### Table {pre}V. Method variants around MAIN ({title}; full values)\n')
    md.append('| variant | N | chi2/N | alpha-beta | beta_M1 | alpha_E1 | gamma_pi | gamma_M1M1 | alpha+beta | gamma_0 |')
    md.append('|---|---|---|---|---|---|---|---|---|---|')
    for t in ('MAIN', 'w30', 'bchptctr', 'bchptctr_w30', 'l2zero', 'gEMfree', 'binsB', 'binmode', 'A2beamtimes', 'fwd_uncorr',
              'exact_constraints', 'no_float_norm'):
        r = M[t]; v, e = r['val'], r['err']
        md.append(f'| {t} | {r["npts"]} | {r["chi2pt"]:.3f} | {v["alpha_m_beta"]:.2f}({e["alpha_m_beta"]:.2f}) | {v["betaM1"]:.2f}({e["betaM1"]:.2f}) | '
                  f'{v["alphaE1"]:.2f}({e["alphaE1"]:.2f}) | {v["gamma_pi"]:.2f}({e["gamma_pi"]:.2f}) | {v["gM1M1"]:.2f}({e["gM1M1"]:.2f}) | '
                  f'{v["alpha_p_beta"]:.2f} | {v["gamma0"]:.3f} |')
    # ---- I: influence
    md.append(f'\n### Table {pre}I. Influence of single sets/points ({title}; shift = variant - MAIN)\n')
    md.append('| removed / added | N | alpha-beta | beta_M1 | gamma_pi | d(alpha-beta) | d(beta_M1) |')
    md.append('|---|---|---|---|---|---|---|')
    md.append(f'| (none: MAIN) | {main["npts"]} | {base["alpha_m_beta"]:.2f} | {base["betaM1"]:.2f} | {base["gamma_pi"]:.2f} | | |')
    for k, v in M['influence'].items():
        md.append(f'| {k} | {v["N"]} | {v["val"]["alpha_m_beta"]:.2f} | {v["val"]["betaM1"]:.2f} | {v["val"]["gamma_pi"]:.2f} | '
                  f'{v["val"]["alpha_m_beta"]-base["alpha_m_beta"]:+.2f} | {v["val"]["betaM1"]-base["betaM1"]:+.2f} |')
    # ---- F: forward points
    md.append(f'\n### Table {pre}F. Impact of the forward (sum-rule) amplitudes of GHP 2015/16 ({title})\n')
    md.append('| treatment | alpha-beta | gamma_pi | beta_M1 | alpha+beta | gamma_0 | n_f | n_g | chi2/N |')
    md.append('|---|---|---|---|---|---|---|---|---|')
    for t, lab in (('MAIN', '**MAIN (floating norms, correlated nuisances)**'), ('fwd_uncorr', 'MAIN with uncorrelated forward errors'),
                   ('exact_constraints', 'MAIN with exact constraints (previous main)')):
        r = M[t]; v, e = r['val'], r['err']
        md.append(f'| {lab} | {v["alpha_m_beta"]:.2f}({e["alpha_m_beta"]:.2f}) | {v["gamma_pi"]:.2f}({e["gamma_pi"]:.2f}) | '
                  f'{v["betaM1"]:.2f}({e["betaM1"]:.2f}) | {v["alpha_p_beta"]:.2f} | {v["gamma0"]:.3f} | '
                  f'{r.get("nf", float("nan")):.2f} | {r.get("ng", float("nan")):.2f} | {r["chi2pt"]:.3f} |')
    for t, lab in (('all_std', 'no norms: exact constraints (KLP18)'), ('all_std_fwdchi2', 'no norms: chi2 points, correlated nuisances'),
                   ('all_std_fwdchi2_uncorr', 'no norms: chi2 points, uncorrelated errors'),
                   ('all_std_nofwd', 'no forward points (only alpha+beta, gamma_0 priors), wide fR2,fR3'),
                   ('all_std_nofwd_w10', 'no forward points, KLP18-size priors on fR2,fR3')):
        r = R[t]; v, e = r['val'], r['err']
        md.append(f'| {lab} | {v["alpha_m_beta"]:.2f}({e["alpha_m_beta"]:.2f}) | {v["gamma_pi"]:.2f}({e["gamma_pi"]:.2f}) | '
                  f'{v["betaM1"]:.2f}({e["betaM1"]:.2f}) | {v["alpha_p_beta"]:.2f} | {v["gamma0"]:.3f} | {r.get("nf", float("nan")):.2f} | '
                  f'{r.get("ng", float("nan")):.2f} | {r["chi2pt"]:.3f} |')
    # ---- C: PG tests
    md.append(f'\n### Table {pre}C. Compatibility (parameter-goodness-of-fit) tests, data chi2 only ({title})\n')
    md.append('| subset | N sub / rest | chi2(all) | chi2(subset) | chi2(rest) | Delta chi2 | alpha-beta subset / rest | beta_M1 subset / rest | gamma_pi subset / rest |')
    md.append('|---|---|---|---|---|---|---|---|---|')
    for k, v in M['PG'].items():
        md.append(f'| {k} | {v["N_sub"]} / {v["N_rest"]} | {v["chi2_all"]:.1f} | {v["chi2_sub"]:.1f} | {v["chi2_rest"]:.1f} | {v["dchi2"]:.1f} | '
                  f'{v["amb_sub"]:.2f}({v["e_amb_sub"]:.2f}) / {v["amb_rest"]:.2f}({v["e_amb_rest"]:.2f}) | {v["b_sub"]:.2f} / {v["b_rest"]:.2f} | '
                  f'{v["gpi_sub"]:.2f} / {v["gpi_rest"]:.2f} |')
    return SYS


SYS_P = config_tables(ML, 'P-', 'paper configuration: all normalisations floated, incl. the pre-1961 sets',
                      'results/main_total_errors.json')
SYS_F = config_tables(MF, '', 'fixed-legacy-normalisation configuration (pre-1961 sets not floated; Sections 5-10 of pwa_report.md)',
                      'results/main_total_errors_fixedlegacy.json')

# ---------------- exact-constraint variants (update_fits.py, no legacy norms) ----------------
VAR = [('old_raw', 'pre-2015 data, TAPS stat. errors only (cf. KLP18 Fit 1)'),
       ('old_std', 'pre-2015 data, stat(+)syst'),
       ('old_pruned', "pre-2015 data, pruned (cf. KLP18 Fit 1'')"),
       ('TAPS_only', 'TAPS 2001 only'),
       ('new_only', 'A2 2022 + HIGS 2022 only'),
       ('new_only_A2', 'A2 2022 (dsigma + Sigma3) only'),
       ('all_std', 'all data'),
       ('all_std_norm', 'all data, floating norms (exact sum-rule constraints)'),
       ('all_3sigma_norm', 'all data, floating norms, 3 sigma outliers removed'),
       ('all_pruned_norm', 'all data, floating norms, McGovern pruning'),
       ('all_noTAPS', 'all data without TAPS 2001'),
       ('all_noA2', 'all data without A2 2022'),
       ('all_std_A2beamtimes_norm', 'A2 per beam time (separate norms)'),
       ('all_std_norm_fwdchi2', 'forward amplitudes as chi2 points (nuisance)'),
       ('all_pruned_norm_fwdchi2', 'pruned, forward amplitudes as chi2 points'),
       ('all_std_norm_w30', 'residual priors x3 wider'),
       ('all_std_norm_bchptctr', 'residual priors centred on BChPT'),
       ('all_std_gEMfree', 'gamma_E1M2 free')]
md.append('\n## Exact-constraint configuration (update_fits.py; sum rules as exact constraints, legacy norms fixed)\n')
md.append('\n### Table U-old. Variant: sum rules as EXACT constraints (KLP18 treatment; previous main = all_std_norm). gamma_E1M2 = 0 unless stated; errors full; main database\n')
md.append('| fit | description | N | chi2/N | alpha_E1 | beta_M1 | alpha-beta | gamma_pi | gamma_E1E1 | gamma_M1M1 | gamma_M1E2 | rho(a-b,gpi) |')
md.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
tex.append('% Table U-old: exact sum-rule constraints variant\n\\begin{tabular}{l r c cccc ccc}\n\\hline')
tex.append('Fit & $N$ & $\\chi^2/N$ & $\\alpha_{E1}$ & $\\beta_{M1}$ & $\\alpha_{E1}-\\beta_{M1}$ & $\\gamma_\\pi$ & $\\gamma_{E1E1}$ & $\\gamma_{M1M1}$ & $\\gamma_{M1E2}$\\\\\\hline')
for t, desc in VAR:
    r = R[t]; v, e = r['val'], r['err']
    md.append(f'| {t} | {desc} | {r["npts"]} | {r["chi2pt"]:.2f} | {v["alphaE1"]:.2f}({e["alphaE1"]:.2f}) | {v["betaM1"]:.2f}({e["betaM1"]:.2f}) | '
              f'{v["alpha_m_beta"]:.2f}({e["alpha_m_beta"]:.2f}) | {v["gamma_pi"]:.2f}({e["gamma_pi"]:.2f}) | {v["gE1E1"]:.2f}({e["gE1E1"]:.2f}) | '
              f'{v["gM1M1"]:.2f}({e["gM1M1"]:.2f}) | {v["gM1E2"]:.2f}({e["gM1E2"]:.2f}) | {r["corr_amb_gpi"]:.2f} |')
    tex.append(f'{desc.replace("_", " ")} & {r["npts"]} & {r["chi2pt"]:.2f} & ' +
               ' & '.join(ftex(v[k], e[k]) for k in ['alphaE1', 'betaM1', 'alpha_m_beta', 'gamma_pi', 'gE1E1', 'gM1M1', 'gM1E2']) + '\\\\')
tex.append('\\hline\n\\end{tabular}\n')

open('results/tables.md', 'w').write('\n'.join(md) + '\n')
open('results/tables.tex', 'w').write('\n'.join(tex) + '\n')
print('\n'.join(md))
print('paper configuration:', json.dumps(SYS_P, indent=1))
