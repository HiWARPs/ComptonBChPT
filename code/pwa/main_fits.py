"""MAIN configuration of the updated PWA (user decision 2026-10-14):
all data, floating normalisations, GHP forward amplitudes Re f, Re g as chi^2 DATA at the bin centres with
correlated nuisances n_f, n_g (alpha+beta = AB_SR + SIG_AB n_f, gamma0 = G0_SR + SIG_G0 n_g move with them; sumrules.py: 2026 14.009(208), -1.025(102); ghp2016 14.0(2), -0.929(105)),
f^R_2, f^R_3 free per bin (wide prior), KLP18 widths for f^R_1,4,5,6.
Runs: uncertainty budget around MAIN, database variants, influence of single sets, parameter-goodness-of-fit.
Writes results/main_fits.json (and results/update/{pulls,p,cov}_M_<tag>.*).
NOTE (2026-10-14, later): this is the FIXED-legacy-normalisation configuration (Oxley, Goldansky, Bernardini, Pugh,
Hyman not floated); the paper's PWA fit floats them too (scenarios/main_fits_legacynorm.py).  It appears in the paper
only as the row "normalizations of the oldest data fixed".  All selections are on the MAIN DATABASE (update_fits.select
drops MAIN_DROP = TAPS 108.42/133 and, since 2026-10-17, the HIGS 2022 sets by default); 'db_allpts' keeps the point,
'db_3sigma' = 3-sigma rule rerun on the main database (check), 'db_withHIGS' = main database + HIGS 2022 (the former
main), 'db_new_only' = A2 + HIGS (test), 'db_A2_only' = the new data of the main database."""
import json
import numpy as np
import pandas as pd
import update_fits as UF

MAIN = dict(sr_mode='chi2', fwd_err='nuisance', float_norm=True, prior_w=(10., 1000., 1000., 10., 1., 1.))
KEYS = UF.KEYS
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


def main():
    d = UF.load()
    out = {}

    def go(tag, dd, bins=UF.BINS_A, **kw):
        res, m, r = UF.fit_variant('M_' + tag, dd, bins, **kw)
        res['tag'] = tag
        out[tag] = res
        print(line(res), flush=True)
        return res, m, r

    # ---------------- MAIN ----------------
    go('MAIN', UF.select(d), **MAIN)
    # ---------------- method / theory variants (budget) ----------------
    go('w30', UF.select(d), **cfg(prior_w=(30., 1000., 1000., 30., 3., 3.)))
    go('bchptctr', UF.select(d), **cfg(prior_center='bchpt'))
    go('bchptctr_w30', UF.select(d), **cfg(prior_center='bchpt', prior_w=(30., 1000., 1000., 30., 3., 3.)))
    go('l2zero', UF.select(d), **cfg(l2='zero'))
    go('gEMfree', UF.select(d), **cfg(fix_gEM=None))
    go('binsB', UF.select(d), bins=UF.BINS_B, **MAIN)
    go('binmode', UF.select(d), **cfg(energy_mode='bin'))
    go('A2beamtimes', UF.select(UF.load(per_beamtime=True)), **MAIN)
    go('fwd_uncorr', UF.select(d), **cfg(fwd_err='uncorr'))
    go('exact_constraints', UF.select(d), sr_mode='constraint', float_norm=True)          # old main
    go('no_float_norm', UF.select(d), **cfg(float_norm=False))
    # ---------------- database variants ----------------
    out['MAIN_DROP'] = [list(x) for x in UF.MAIN_DROP]
    go('db_allpts', UF.select(d, keep_outlier=True), **MAIN)                            # all points incl. the outlier
    pulls = pd.read_csv('results/update/pulls_M_MAIN.csv')
    o3 = pulls[abs(pulls.pull) > 3]
    rm = [(a, b, c) for a, b, c in zip(o3.dataset, o3.E, o3.theta_lab)]
    allrm = list(rm)
    for it in range(5):                       # iterate the 3 sigma rule until no point exceeds 3 sigma
        res, m, r = UF.fit_variant('M_tmp', UF.select(d, drop_points=allrm), **MAIN)
        pp = pd.read_csv('results/update/pulls_M_tmp.csv')
        o = pp[abs(pp.pull) > 3]
        if len(o) == 0:
            break
        allrm += [(a, b, c) for a, b, c in zip(o.dataset, o.E, o.theta_lab)]
    go('db_3sigma', UF.select(d, drop_points=allrm), **MAIN)
    out['db_3sigma']['removed'] = [list(x) for x in allrm]
    out['MAIN_outliers'] = o3[['dataset', 'obs', 'E', 'theta_lab', 'value', 'err', 'pred', 'pull']].to_dict('records')
    go('db_pruned', UF.select(d, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points']), **MAIN)
    go('db_modern', UF.select(d, sets=MODERN), **MAIN)
    go('db_new_only', UF.select(d, sets=NEW), **MAIN)                        # A2 + HIGS (test; HIGS not in the main database)
    go('db_A2_only', UF.select(d, sets=A2), **MAIN)
    go('db_withHIGS', UF.select(d, with_higs=True), **MAIN)                  # main database + HIGS 2022 (the former main)
    go('db_noTAPS', UF.select(d, exclude_sets=['TAPS2001']), **MAIN)
    go('db_noA2', UF.select(d, exclude_sets=A2), **MAIN)
    go('db_noLegacy', UF.select(d, exclude_sets=LEGACY), **MAIN)
    go('db_old_only', UF.select(d, old_only=True), **MAIN)
    go('db_old_only_pruned', UF.select(d, old_only=True, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points']), **MAIN)
    go('db_pruned_withHIGS', UF.select(d, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points'], with_higs=True), **MAIN)
    # ---------------- influence of single sets / points ----------------
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
        res, m, r = UF.fit_variant('M_tmp', UF.select(d, **kw), **MAIN)
        infl[nm] = dict(N=res['npts'], chi2=res['chi2'], val=res['val'], err=res['err'],
                        shift={k: res['val'][k] - out['MAIN']['val'][k] for k in KEYS})
        print(f'   influence {nm:34s} a-b={res["val"]["alpha_m_beta"]:.2f} b={res["val"]["betaM1"]:.2f} gpi={res["val"]["gamma_pi"]:.2f}')
    out['influence'] = infl
    # ---------------- parameter goodness of fit ----------------
    pg = {}
    # 'all' = subset + rest: MAIN for subsets of the main database, db_withHIGS for the HIGS-containing subsets
    for name, sub, allkey in (('TAPS', ['TAPS2001'], 'MAIN'), ('A2', A2, 'MAIN'), ('HIGS', HIGS, 'db_withHIGS'),
                              ('new vs old', NEW, 'db_withHIGS')):
        rA, _, _ = UF.fit_variant('M_tmp', UF.select(d, sets=sub), **MAIN)
        rB, _, _ = UF.fit_variant('M_tmp', UF.select(d, exclude_sets=sub), **MAIN)
        pg[name] = dict(all_key=allkey, chi2_all=out[allkey]['chi2'], chi2_sub=rA['chi2'], chi2_rest=rB['chi2'],
                        dchi2=out[allkey]['chi2'] - rA['chi2'] - rB['chi2'], N_sub=rA['npts'], N_rest=rB['npts'],
                        amb_sub=rA['val']['alpha_m_beta'], amb_rest=rB['val']['alpha_m_beta'],
                        e_amb_sub=rA['err']['alpha_m_beta'], e_amb_rest=rB['err']['alpha_m_beta'],
                        gpi_sub=rA['val']['gamma_pi'], gpi_rest=rB['val']['gamma_pi'],
                        b_sub=rA['val']['betaM1'], b_rest=rB['val']['betaM1'])
        print('   PG', name, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in pg[name].items()}, flush=True)
    out['PG'] = pg
    for f in ('pulls_M_tmp.csv', 'p_M_tmp.npy', 'cov_M_tmp.npy'):
        import os
        if os.path.exists('results/update/' + f):
            os.remove('results/update/' + f)
    json.dump(out, open('results/main_fits.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
