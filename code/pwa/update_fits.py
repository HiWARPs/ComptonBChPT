"""Updated PWA fits (A2 2022 + pre-2015 world data; HIGS 2022 as a test), variants for database choices,
normalisations, outliers, sum-rule treatment, binning and regularisation.
Writes results/update_fits.json, results/update/pulls_<tag>.csv and results/update/p_<tag>.npy

MAIN DATABASE (author decisions 2026-10-14 and 2026-10-17, same as the BChPT fits, work/fits/lib/common.py):
pre-2015 world data + A2 2022 at E <= 150 MeV MINUS the single 3-sigma outlier, the TAPS (Olmos de Leon 2001) cross
section at 108.42 MeV, 133 deg (MAIN_DROP; BChPT pid OlmosdeLeon2001_029).  The HIGS 2022 sets (Li et al.,
HIGS_2022_circ and HIGS_2022_lin, 16 points) are EXCLUDED from the main database since 2026-10-17 (author decision);
until then the main database was all data minus MAIN_DROP.  `select()` drops MAIN_DROP unless keep_outlier=True and
drops the HIGS sets unless with_higs=True or they are named explicitly in `sets`, so EVERY selection built through it
(main fit, method variants, subsets containing TAPS, influence and PG tests, figures) is a subset of the main
database; the GMPF-pruned selections never contained the point.  keep_outlier=True gives the all-points database
(variant 'db_allpts'); with_higs=True gives the former main database (variants 'all_withHIGS' here, 'db_withHIGS' in
main_fits*.py); sets=HIGS gives the HIGS-alone test ('new_only_HIGS')."""
import json
import sys
import numpy as np
import pandas as pd
import pwa_model as PM
from sumrules import ForwardAmps
import compton as C

BINS_A = [59, 69, 79, 89, 99, 109, 117, 127, 135, 143, 150]          # KLP18 bins
BINS_B = [59, 69, 81, 92, 103, 113, 124, 135, 143, 150]              # re-centred on the HIGS/A2 energies
KEYS = ['alphaE1', 'betaM1', 'alpha_m_beta', 'gamma_pi', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2', 'alpha_p_beta', 'gamma0']
MAIN_DROP = [('TAPS2001', 108.42, 133.0)]      # (dataset, E [MeV], theta_lab [deg]) of the point(s) excluded from the main database
HIGS = ['HIGS_2022_circ', 'HIGS_2022_lin']     # HIGS 2022 (Li et al.): excluded from the main database since 2026-10-17 (test only)
PRUNE = dict(drop_sets=['Oxley', 'Bernardini', 'Baranov1974_th150'], drop_points=[('TAPS2001', 108.42, 133.0)])


def load(per_beamtime=False):
    return pd.read_csv('results/db_update_A2beamtimes.csv' if per_beamtime else 'results/db_update.csv')


def select(d, sets=None, exclude_sets=(), drop_points=(), odl_stat_only=False, new_only=False, old_only=False,
           keep_outlier=False, with_higs=False):
    """data selection: MAIN_DROP is removed from every selection unless keep_outlier=True, and the HIGS 2022 sets are
    removed unless with_higs=True or they are named explicitly in `sets` (see module docstring)"""
    d = d.copy()
    keep = np.ones(len(d), bool)
    if sets is not None:
        keep &= d.dataset.isin(sets).values
    keep &= ~d.dataset.isin(list(exclude_sets)).values
    if not with_higs and not (sets is not None and any(s in HIGS for s in sets)):
        keep &= ~d.dataset.isin(HIGS).values            # HIGS 2022 not in the main database (2026-10-17)
    for ds, E, th in list(drop_points) + ([] if keep_outlier else list(MAIN_DROP)):
        keep &= ~((d.dataset == ds) & (abs(d.E - E) < 0.05) & (abs(d.theta_lab - th) < 0.5)).values
    if new_only:
        keep &= (d.group == 'new').values
    if old_only:
        keep &= (d.group == 'old').values
    d = d[keep].copy()
    if odl_stat_only:
        m = (d.dataset == 'TAPS2001').values
        d.loc[m, 'err'] = d.loc[m, 'stat']
    return d


def build_model(d, bins=BINS_A, **kw):
    opts = dict(l2='bchpt', energy_mode='actual', fix_gEM=0.0, sr_mode='constraint', fwd=ForwardAmps())
    opts.update(kw)
    d = d.copy()
    b = np.array(bins, float)
    d['bin'] = b[np.argmin(abs(d.E.values[:, None] - b[None, :]), axis=1)]
    return PM.Model(d, bins, **opts)


def summarize(tag, m, r, extra=None):
    pulls, pred = m.pulls(r.p)
    d = m.d.copy(); d['pred'] = pred; d['pull'] = pulls
    d.to_csv(f'results/update/pulls_{tag}.csv', index=False)
    np.save(f'results/update/p_{tag}.npy', r.p)
    np.save(f'results/update/cov_{tag}.npy', r.cov)
    r2 = m.fit_globals_only(r.p)
    per = d.groupby('dataset').apply(lambda g: pd.Series(dict(N=len(g), chi2=float(np.sum(g.pull**2)))),
                                     include_groups=False)
    norms = {}
    if m.ngroups:
        npm = m.ng + m.nr * m.nb
        for j, g in enumerate(m.ngroups):
            norms[g] = dict(N=float(1 + m.gsig[j] * r.p[npm + j]), err=float(m.gsig[j] * np.sqrt(r.cov[npm + j, npm + j])))
    res = dict(tag=tag, npts=int(len(d)), chi2=float(r.chi2), chi2pt=float(r.chi2 / len(d)), npar=int(r.free.sum()),
               val={k: float(r.der[k]) for k in KEYS}, err={k: float(r.der_err[k]) for k in KEYS},
               err_stage2={k: float(r2.der_err[k]) for k in KEYS},
               corr_amb_gpi=r.corr('alpha_m_beta', 'gamma_pi'),
               per_dataset={k: dict(N=int(v.N), chi2=float(v.chi2)) for k, v in per.iterrows()},
               outliers3=d[abs(d.pull) > 3][['dataset', 'obs', 'E', 'theta_lab', 'value', 'err', 'pred', 'pull']].to_dict('records'),
               norms=norms)
    if 'nf' in m.names:
        res['nf'] = float(r.p[m.names.index('nf')]); res['ng'] = float(r.p[m.names.index('ng')])
    if extra:
        res.update(extra)
    return res


def fit_variant(tag, d, bins=BINS_A, **kw):
    m = build_model(d, bins, **kw)
    r = m.fit(with_prior=True)
    return summarize(tag, m, r), m, r


def line(res):
    v, e = res['val'], res['err']
    return (f'{res["tag"]:28s} N={res["npts"]:3d} chi2/N={res["chi2pt"]:.3f}  a={v["alphaE1"]:5.2f}({e["alphaE1"]:.2f}) '
            f'b={v["betaM1"]:5.2f}({e["betaM1"]:.2f}) a-b={v["alpha_m_beta"]:5.2f}({e["alpha_m_beta"]:.2f}) '
            f'gpi={v["gamma_pi"]:5.2f}({e["gamma_pi"]:.2f}) gEE={v["gE1E1"]:5.2f}({e["gE1E1"]:.2f}) '
            f'gMM={v["gM1M1"]:5.2f}({e["gM1M1"]:.2f}) gEM={v["gE1M2"]:5.2f} gME={v["gM1E2"]:5.2f}({e["gM1E2"]:.2f})')


def main():
    d = load()
    out = {}

    def go(tag, dd, bins=BINS_A, **kw):
        res, m, r = fit_variant(tag, dd, bins, **kw)
        out[tag] = res
        print(line(res), flush=True)
        return res, m, r

    # ---------- baselines ----------
    go('old_raw', select(d, old_only=True, odl_stat_only=True))            # analogue of KLP18 Fit 1
    go('old_std', select(d, old_only=True))
    go('old_pruned', select(d, old_only=True, exclude_sets=PRUNE['drop_sets'], drop_points=PRUNE['drop_points']))
    go('new_only', select(d, new_only=True, with_higs=True))               # A2 2022 + HIGS 2022 (test: HIGS not in the main database)
    go('new_only_A2', select(d, sets=['A2_2022_CS', 'A2_2022_S3']))        # the new data of the main database
    go('new_only_A2CS', select(d, sets=['A2_2022_CS']))
    go('new_only_HIGS', select(d, sets=['HIGS_2022_circ', 'HIGS_2022_lin']))
    go('TAPS_only', select(d, sets=['TAPS2001']))
    # ---------- main updated fits ----------
    res, m, r = go('all_raw', select(d, odl_stat_only=True))
    res, m, r = go('all_std', select(d))
    go('all_std_norm', select(d), float_norm=True)
    go('all_pruned', select(d, exclude_sets=PRUNE['drop_sets'], drop_points=PRUNE['drop_points']))
    go('all_pruned_norm', select(d, exclude_sets=PRUNE['drop_sets'], drop_points=PRUNE['drop_points']), float_norm=True)
    # 3-sigma outliers of all_std (one pass) and iterated
    pulls = pd.read_csv('results/update/pulls_all_std.csv')
    o3 = pulls[abs(pulls.pull) > 3]
    dp = [(a, b, c) for a, b, c in zip(o3.dataset, o3.E, o3.theta_lab)]
    go('all_3sigma', select(d, drop_points=dp))
    pulls = pd.read_csv('results/update/pulls_all_std_norm.csv')
    o3 = pulls[abs(pulls.pull) > 3]
    dpn = [(a, b, c) for a, b, c in zip(o3.dataset, o3.E, o3.theta_lab)]
    go('all_3sigma_norm', select(d, drop_points=dpn), float_norm=True)
    out['all_3sigma_norm']['removed'] = dpn
    dd = select(d)
    allrm = []
    for it in range(5):
        res, m, r = fit_variant(f'tmp', dd)
        pp = pd.read_csv('results/update/pulls_tmp.csv')
        o = pp[abs(pp.pull) > 3]
        if len(o) == 0:
            break
        allrm += [(a, b, c) for a, b, c in zip(o.dataset, o.E, o.theta_lab)]
        dd = select(d, drop_points=allrm)
    go('all_3sigma_iter', dd)
    out['all_3sigma_iter']['removed'] = allrm
    go('all_noTAPS', select(d, exclude_sets=['TAPS2001']))
    go('all_noA2', select(d, exclude_sets=['A2_2022_CS', 'A2_2022_S3']))
    go('all_noA2S3', select(d, exclude_sets=['A2_2022_S3']))
    go('all_withHIGS', select(d, with_higs=True))                          # main database + HIGS 2022 (the former main)
    go('all_withHIGS_norm', select(d, with_higs=True), float_norm=True)    # (former all_std_norm)
    go('all_noLegacy', select(d, exclude_sets=['Oxley', 'Goldansky', 'Bernardini', 'Pugh', 'Hyman']))
    go('modern', select(d, sets=['TAPS2001', 'MacGibbon1995_tagged', 'MacGibbon1995_untagged', 'Hallin1993', 'Federspiel1991',
                                 'A2_2022_CS', 'A2_2022_S3']))
    # ---------- per-beamtime A2 ----------
    dbt = load(per_beamtime=True)
    go('all_std_A2beamtimes_norm', select(dbt), float_norm=True)
    go('all_std_A2beamtimes', select(dbt))
    # ---------- theory/method variants on all_std ----------
    go('all_std_gEMfree', select(d), fix_gEM=None)
    go('all_std_l2zero', select(d), l2='zero')
    go('all_std_binsB', select(d), bins=BINS_B)
    go('all_std_binmode', select(d), energy_mode='bin')
    go('all_std_prior_x0.5', select(d), prior_w=(5., 5., 5., 5., .5, .5))
    go('all_std_prior_x2', select(d), prior_w=(20., 20., 20., 20., 2., 2.))
    go('all_std_noresid', select(d), prior_w=(1e-4,) * 6)
    go('all_std_ab13.8', select(d), ab=13.8)
    go('all_std_ab14.2', select(d), ab=14.2)
    go('all_std_g0-0.83', select(d), g0=-0.829)
    go('all_std_g0-1.03', select(d), g0=-1.029)
    # forward amplitudes as chi2 data (nuisance-correlated / uncorrelated); fR2, fR3 then free (wide prior)
    PW_CHI2 = (10., 1000., 1000., 10., 1., 1.)
    go('all_std_fwdchi2', select(d), sr_mode='chi2', fwd_err='nuisance', prior_w=PW_CHI2)
    go('all_std_fwdchi2_uncorr', select(d), sr_mode='chi2', fwd_err='uncorr', prior_w=PW_CHI2)
    go('all_std_norm_fwdchi2', select(d), sr_mode='chi2', float_norm=True, prior_w=PW_CHI2)
    go('all_pruned_norm_fwdchi2', select(d, exclude_sets=PRUNE['drop_sets'], drop_points=PRUNE['drop_points']),
       sr_mode='chi2', float_norm=True, prior_w=PW_CHI2)
    go('new_only_fwdchi2', select(d, new_only=True), sr_mode='chi2', prior_w=PW_CHI2)
    go('all_std_nofwd', select(d), sr_mode='chi2', fwd_err='uncorr', fwd_scale=1e4, prior_w=PW_CHI2)
    go('all_std_norm_nofwd', select(d), sr_mode='chi2', fwd_err='uncorr', fwd_scale=1e4, prior_w=PW_CHI2, float_norm=True)
    go('all_std_nofwd_w10', select(d), sr_mode='chi2', fwd_err='uncorr', fwd_scale=1e4)
    # ---------- regularisation of the residual functions ----------
    REG = {'w30': dict(prior_w=(30., 30., 30., 30., 3., 3.)),
           'bchptctr': dict(prior_center='bchpt'),
           'bchptctr_w30': dict(prior_center='bchpt', prior_w=(30., 30., 30., 30., 3., 3.))}
    for rn, kw in REG.items():
        go(f'all_std_{rn}', select(d), **kw)
        go(f'all_std_norm_{rn}', select(d), float_norm=True, **kw)
        go(f'all_pruned_norm_{rn}', select(d, exclude_sets=PRUNE['drop_sets'], drop_points=PRUNE['drop_points']),
           float_norm=True, **kw)
        go(f'new_only_{rn}', select(d, new_only=True), **kw)
        go(f'old_raw_{rn}', select(d, old_only=True, odl_stat_only=True), **kw)
        go(f'old_pruned_{rn}', select(d, old_only=True, exclude_sets=PRUNE['drop_sets'], drop_points=PRUNE['drop_points']), **kw)
    go('old_noTAPS', select(d, old_only=True, exclude_sets=['TAPS2001']))
    # ---------- parameter-goodness-of-fit (compatibility) tests, data chi2 only ----------
    pg = {}
    # 'all' = the fit to subset + rest: the main database (all_std) for subsets of it, the with-HIGS fit for the
    # HIGS-containing subsets (HIGS 2022 is not part of the main database since 2026-10-17)
    for name, sub, allkey in (('TAPS', ['TAPS2001'], 'all_std'), ('A2', ['A2_2022_CS', 'A2_2022_S3'], 'all_std'),
                              ('HIGS', HIGS, 'all_withHIGS'), ('new', ['A2_2022_CS', 'A2_2022_S3'] + HIGS, 'all_withHIGS')):
        rA, _, _ = fit_variant('pgA', select(d, sets=sub))
        rB, _, _ = fit_variant('pgB', select(d, exclude_sets=sub))
        pg[name] = dict(all_key=allkey, chi2_all=out[allkey]['chi2'], chi2_sub=rA['chi2'], chi2_rest=rB['chi2'],
                        dchi2=out[allkey]['chi2'] - rA['chi2'] - rB['chi2'],
                        amb_sub=rA['val']['alpha_m_beta'], amb_rest=rB['val']['alpha_m_beta'],
                        gpi_sub=rA['val']['gamma_pi'], gpi_rest=rB['val']['gamma_pi'])
        print('PG', name, pg[name], flush=True)
    out['PG'] = pg
    json.dump(out, open('results/update_fits.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
