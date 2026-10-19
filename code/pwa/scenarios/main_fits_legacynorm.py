"""MAIN configuration of the updated PWA (user decision 2026-10-14) = the PWA fit of the paper:
all data of the MAIN DATABASE, floating normalisations of ALL sets (also the pre-1961 sets: 10 %, Goldansky 8 %),
forward amplitudes Re f, Re g (sumrules.SOURCE) as chi^2 DATA at the bin centres with correlated nuisances n_f, n_g
(alpha+beta = AB_SR + SIG_AB n_f, gamma0 = G0_SR + SIG_G0 n_g move with them; sumrules.py: 2026 14.009(208), -1.025(102); ghp2016 14.0(2), -0.929(105)), f^R_2, f^R_3 free per bin (wide prior),
KLP18 widths for f^R_1,4,5,6.

MAIN DATABASE (author decisions 2026-10-14 and 2026-10-17; same as the BChPT fits): pre-2015 world data + A2 2022 at
E <= 150 MeV MINUS the single 3-sigma outlier TAPS 108.42 MeV / 133 deg (update_fits.MAIN_DROP = BChPT pid
OlmosdeLeon2001_029).  The HIGS 2022 sets (Li et al., 16 points) are EXCLUDED from the main database since 2026-10-17
(author decision; until then the main database was all data minus MAIN_DROP, N = 253).  update_fits.select drops
MAIN_DROP and the HIGS sets by default, so MAIN, every method variant, every database variant (subsets containing
TAPS), the influence and the PG tests are all relative to the main database (N = 237); the GMPF-pruned selections
never contained the point.
  'db_allpts'      : main database incl. the outlier (N = 238, keep_outlier=True);
  'MAIN_DROP_check': the iterated 3-sigma search on db_allpts, which must re-derive exactly MAIN_DROP (else: error);
  'db_3sigma'      : the iterated 3-sigma search rerun on the main database (check: nothing else flagged -> = MAIN);
  'db_withHIGS'    : main database + HIGS 2022 (N = 253; the former main fit; with_higs=True); 'db_pruned_withHIGS';
  'db_new_only'    : A2 2022 + HIGS 2022 (test; HIGS not in the main database); 'db_A2_only' = the new data of the main database;
  'PG'             : parameter-goodness-of-fit; for the HIGS-containing subsets 'all' is the with-HIGS fit;
  'gEM_profile'    : MAIN with gamma_E1M2 fixed at -3 ... 5 (chi2 and parameters);
  'method'         : the "method" uncertainty of the paper = quadrature sum of the shifts of the variants
                     w30, bchptctr, l2zero, gEMfree, binsB, binmode, A2beamtimes, fwd_uncorr; 'total' = stat (+) method.
Runs: uncertainty budget around MAIN, database variants, influence of single sets, parameter-goodness-of-fit.
Writes results/main_fits_legacynorm.json (and results/update/{pulls,p,cov}_L_<tag>.*).  Run from work/pwa/ (4-5 s)."""
import sys as _sys; _sys.path.insert(0, '.')
import update_fits as _UF
_LEG = _UF.LEGACY_NORM   # Baranov et al. 2001 delta (was 10%, Goldansky 8%)
def load_leg(per_beamtime=False):
    d = _UF.load(per_beamtime=per_beamtime)
    for s, e in _LEG.items():
        m = d.dataset == s
        d.loc[m, 'norm_group'] = s; d.loc[m, 'norm_err'] = e
    return d

import json
import os
import numpy as np
import pandas as pd
import update_fits as UF

MAIN = dict(sr_mode='chi2', fwd_err='nuisance', float_norm=True, prior_w=(10., 1000., 1000., 10., 1., 1.))
KEYS = UF.KEYS
METHOD_VARIANTS = ['w30', 'bchptctr', 'l2zero', 'gEMfree', 'binsB', 'binmode', 'A2beamtimes', 'fwd_uncorr']
LEGACY = ['Oxley', 'Goldansky', 'Bernardini', 'Pugh', 'Hyman']
HIGS = list(UF.HIGS)                                   # not in the main database (2026-10-17); test only
A2 = ['A2_2022_CS', 'A2_2022_S3']
NEW = A2 + HIGS                                        # all new data (A2 + HIGS): test, not a subset of the main database
MODERN = ['TAPS2001', 'MacGibbon1995_tagged', 'MacGibbon1995_untagged', 'Hallin1993', 'Federspiel1991'] + A2


def cfg(**kw):
    c = dict(MAIN)
    c.update(kw)
    return c


def line(r):
    v, e = r['val'], r['err']
    return (f'{r["tag"]:34s} N={r["npts"]:3d} chi2/N={r["chi2pt"]:.3f} a-b={v["alpha_m_beta"]:5.2f}({e["alpha_m_beta"]:.2f}) '
            f'b={v["betaM1"]:4.2f}({e["betaM1"]:.2f}) gpi={v["gamma_pi"]:5.2f}({e["gamma_pi"]:.2f}) '
            f'a+b={v["alpha_p_beta"]:.2f}({e["alpha_p_beta"]:.2f}) g0={v["gamma0"]:.3f}({e["gamma0"]:.3f}) '
            f'gEE={v["gE1E1"]:.2f} gMM={v["gM1M1"]:.2f}({e["gM1M1"]:.2f}) gME={v["gM1E2"]:.2f} gEM={v["gE1M2"]:.2f}')


def outlier_search(d, start_drop=(), keep_outlier=False, tag='L_tmp', **kw):
    """iterated 3-sigma rule (as in KLP18): refit, drop every point with |pull| > 3, repeat until none is left.
    Returns the list of removed (dataset, E, theta_lab).  keep_outlier=True starts from the all-points database."""
    allrm = list(start_drop)
    for it in range(5):
        res, m, r = UF.fit_variant(tag, UF.select(d, drop_points=allrm, keep_outlier=keep_outlier), **kw)
        pp = pd.read_csv(f'results/update/pulls_{tag}.csv')
        o = pp[abs(pp.pull) > 3]
        if len(o) == 0:
            break
        allrm += [(a, b, c) for a, b, c in zip(o.dataset, o.E, o.theta_lab)]
    return allrm


def main():
    d = load_leg()
    out = {}

    def go(tag, dd, bins=UF.BINS_A, **kw):
        res, m, r = UF.fit_variant('L_' + tag, dd, bins, **kw)
        res['tag'] = tag
        out[tag] = res
        print(line(res), flush=True)
        return res, m, r

    # ---------------- MAIN (main database) ----------------
    out['MAIN_DROP'] = [list(x) for x in UF.MAIN_DROP]
    go('MAIN', UF.select(d), **MAIN)
    # ---------------- method / theory variants (budget), all on the main database ----------------
    go('w30', UF.select(d), **cfg(prior_w=(30., 1000., 1000., 30., 3., 3.)))
    go('bchptctr', UF.select(d), **cfg(prior_center='bchpt'))
    go('bchptctr_w30', UF.select(d), **cfg(prior_center='bchpt', prior_w=(30., 1000., 1000., 30., 3., 3.)))
    go('l2zero', UF.select(d), **cfg(l2='zero'))
    go('gEMfree', UF.select(d), **cfg(fix_gEM=None))
    # profile in the fixed gamma_E1M2 (paper: "fixing gamma_E1M2 anywhere between -3 and 5 ...")
    out['gEM_profile'] = []
    for g in (-3., -2., -1., 0., 1., 2., 3., 5.):
        r, _, _ = UF.fit_variant('L_tmp', UF.select(d), **cfg(fix_gEM=g))
        out['gEM_profile'].append(dict(gE1M2=g, chi2=r['chi2'], val=r['val'], err=r['err']))
    go('binsB', UF.select(d), bins=UF.BINS_B, **MAIN)
    go('binmode', UF.select(d), **cfg(energy_mode='bin'))
    go('A2beamtimes', UF.select(load_leg(per_beamtime=True)), **MAIN)
    go('fwd_uncorr', UF.select(d), **cfg(fwd_err='uncorr'))
    go('exact_constraints', UF.select(d), sr_mode='constraint', float_norm=True)          # old main
    go('no_float_norm', UF.select(d), **cfg(float_norm=False))
    # method uncertainty of the paper: quadrature sum of the shifts of METHOD_VARIANTS; total = stat (+) method
    base = out['MAIN']
    meth = {k: float(np.sqrt(sum((out[t]['val'][k] - base['val'][k])**2 for t in METHOD_VARIANTS))) for k in KEYS}
    out['method'] = dict(variants=METHOD_VARIANTS,
                         shifts={t: {k: out[t]['val'][k] - base['val'][k] for k in KEYS} for t in METHOD_VARIANTS},
                         method=meth, stat={k: base['err'][k] for k in KEYS},
                         total={k: float(np.hypot(base['err'][k], meth[k])) for k in KEYS})
    print('   method (quadrature of %s):' % ', '.join(METHOD_VARIANTS),
          {k: round(meth[k], 2) for k in ('alpha_m_beta', 'betaM1', 'alphaE1', 'gamma_pi')},
          ' total:', {k: round(out['method']['total'][k], 2) for k in ('alpha_m_beta', 'betaM1', 'alphaE1', 'gamma_pi')})
    # ---------------- database variants ----------------
    # all points incl. the outlier, and the verification that the iterated 3-sigma rule on it flags exactly MAIN_DROP
    go('db_allpts', UF.select(d, keep_outlier=True), **MAIN)
    pulls = pd.read_csv('results/update/pulls_L_db_allpts.csv')
    o3 = pulls[abs(pulls.pull) > 3]
    out['db_allpts_outliers'] = o3[['dataset', 'obs', 'E', 'theta_lab', 'value', 'err', 'pred', 'pull']].to_dict('records')
    derived = outlier_search(d, keep_outlier=True, **MAIN)
    out['MAIN_DROP_check'] = dict(derived=[list(x) for x in derived], ok=sorted(derived) == sorted(UF.MAIN_DROP))
    print('   3-sigma search on all points flags:', derived, '-> MAIN_DROP verified' if out['MAIN_DROP_check']['ok'] else '!!! MISMATCH')
    if not out['MAIN_DROP_check']['ok']:
        raise RuntimeError(f'iterated 3-sigma search on the all-points database gives {derived}, not MAIN_DROP={UF.MAIN_DROP}')
    # 3-sigma rule rerun on the main database (check): expected to flag nothing, i.e. db_3sigma = MAIN
    pulls = pd.read_csv('results/update/pulls_L_MAIN.csv')
    o3 = pulls[abs(pulls.pull) > 3]
    out['MAIN_outliers'] = o3[['dataset', 'obs', 'E', 'theta_lab', 'value', 'err', 'pred', 'pull']].to_dict('records')
    allrm = outlier_search(d, **MAIN)
    go('db_3sigma', UF.select(d, drop_points=allrm), **MAIN)
    out['db_3sigma']['removed'] = [list(x) for x in allrm]
    print('   3-sigma search on the main database flags:', allrm or 'nothing (db_3sigma = MAIN)')
    go('db_pruned', UF.select(d, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points']), **MAIN)
    go('db_modern', UF.select(d, sets=MODERN), **MAIN)
    go('db_new_only', UF.select(d, sets=NEW), **MAIN)                        # A2 + HIGS (test; HIGS not in the main database)
    go('db_A2_only', UF.select(d, sets=A2), **MAIN)
    go('db_HIGS_only', UF.select(d, sets=HIGS), **MAIN)                      # HIGS 2022 alone (test; Table pwa)
    go('db_withHIGS', UF.select(d, with_higs=True), **MAIN)                  # main database + HIGS 2022 (the former main fit)
    go('db_noTAPS', UF.select(d, exclude_sets=['TAPS2001']), **MAIN)
    go('db_noA2', UF.select(d, exclude_sets=A2), **MAIN)
    go('db_noLegacy', UF.select(d, exclude_sets=LEGACY), **MAIN)
    go('db_old_only', UF.select(d, old_only=True), **MAIN)
    go('db_old_only_pruned', UF.select(d, old_only=True, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points']), **MAIN)
    go('db_pruned_withHIGS', UF.select(d, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points'], with_higs=True), **MAIN)
    # ---------------- influence of single sets / points (shift = variant - MAIN) ----------------
    infl = {}
    for nm, kw in (('drop Oxley', dict(exclude_sets=['Oxley'])),
                   ('drop Oxley 60/150', dict(drop_points=[('Oxley', 60.0, 150.0)])),
                   ('keep TAPS 108.4/133 (all points)', dict(keep_outlier=True)),
                   ('add HIGS 2022 (former main database)', dict(with_higs=True)),
                   ('drop Baranov 150deg', dict(exclude_sets=['Baranov1974_th150'])),
                   ('drop Bernardini', dict(exclude_sets=['Bernardini'])),
                   ('drop Goldansky', dict(exclude_sets=['Goldansky'])),
                   ('drop Oxley+Goldansky', dict(exclude_sets=['Oxley', 'Goldansky'])),
                   ('drop A2 Sigma3', dict(exclude_sets=['A2_2022_S3']))):
        res, m, r = UF.fit_variant('L_tmp', UF.select(d, **kw), **MAIN)
        infl[nm] = dict(N=res['npts'], chi2=res['chi2'], val=res['val'], err=res['err'],
                        shift={k: res['val'][k] - base['val'][k] for k in KEYS})
        print(f'   influence {nm:34s} a-b={res["val"]["alpha_m_beta"]:.2f} ({infl[nm]["shift"]["alpha_m_beta"]:+.2f}) '
              f'b={res["val"]["betaM1"]:.2f} ({infl[nm]["shift"]["betaM1"]:+.2f}) gpi={res["val"]["gamma_pi"]:.2f}')
    out['influence'] = infl
    # ---------------- parameter goodness of fit (subsets of the main database) ----------------
    pg = {}
    # 'all' = subset + rest: MAIN for subsets of the main database, db_withHIGS for the HIGS-containing subsets
    for name, sub, allkey in (('TAPS', ['TAPS2001'], 'MAIN'), ('A2', A2, 'MAIN'), ('HIGS', HIGS, 'db_withHIGS'),
                              ('new vs old', NEW, 'db_withHIGS')):
        rA, _, _ = UF.fit_variant('L_tmp', UF.select(d, sets=sub), **MAIN)
        rB, _, _ = UF.fit_variant('L_tmp', UF.select(d, exclude_sets=sub), **MAIN)
        pg[name] = dict(all_key=allkey, chi2_all=out[allkey]['chi2'], chi2_sub=rA['chi2'], chi2_rest=rB['chi2'],
                        dchi2=out[allkey]['chi2'] - rA['chi2'] - rB['chi2'], N_sub=rA['npts'], N_rest=rB['npts'],
                        amb_sub=rA['val']['alpha_m_beta'], amb_rest=rB['val']['alpha_m_beta'],
                        e_amb_sub=rA['err']['alpha_m_beta'], e_amb_rest=rB['err']['alpha_m_beta'],
                        gpi_sub=rA['val']['gamma_pi'], gpi_rest=rB['val']['gamma_pi'],
                        b_sub=rA['val']['betaM1'], b_rest=rB['val']['betaM1'])
        print('   PG', name, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in pg[name].items()}, flush=True)
    out['PG'] = pg
    for f in ('pulls_L_tmp.csv', 'p_L_tmp.npy', 'cov_L_tmp.npy'):
        if os.path.exists('results/update/' + f):
            os.remove('results/update/' + f)
    json.dump(out, open('results/main_fits_legacynorm.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
