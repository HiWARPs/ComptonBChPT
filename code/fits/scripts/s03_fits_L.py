"""Step 3: fits in the low-energy region L (E_lab <= 150 MeV; variant <= 170 MeV).

F1 : alpha, beta free, Baldin alpha+beta = BALDIN (lib/common.py: 14.01(21) ForwardRCS 2026; 14.0(2) GHP15 with FWD_SOURCE=ghp2016) as penalty (spin pols at NNLO)
F2 : alpha, beta free, forward points dsigma(0), Sigma_2z(0) (E <= FWD_EMAX) with correlated nuisances, no Baldin
F3 : F2 + gamma_M1M1 and gamma_E1E1 free (gamma_0 constrained by the forward Sigma_2z points)
F3b: F1 + gamma_pi free at fixed gamma_0 (direction dgM1M1 = -dgE1E1 = 1/2)
F3c: F1 + gamma_M1M1 free alone (MPG-2013 style)
F4 : F1/F2 with E <= 170 MeV (OdL to 164 MeV, Hallin <= 150 MeV as in MPG Set III; variant Hallin <= 170)
Database variants a-j (see report).  Also: Mornacchi benchmark, theory-covariance fits, order-by-order fits,
parameter variations.

MAIN DATABASE (since 2026-10-17): 'd_all' denotes the MAIN database = pre-2015 world data (OLD_L) + A2 2022 (A2),
E_lab <= 150 MeV, MINUS the 3-sigma outlier(s) MAIN_DROP (lib/common.py; the single TAPS point OlmosdeLeon2001_029 at
108.42 MeV, 133 deg), plus the GHP forward points (F2/F3 kinds).  The HIgS 2022 sets (Li et al., HIGS_2022_CIRC and
HIGS_2022_LIN, 16 points) are EXCLUDED from the main database since 2026-10-17 (author decision); until then the
main database was OLD_L + A2 + HIGS (minus MAIN_DROP).  The KEY NAMES of the 'd_all' family are UNCHANGED, so that
every consumer of 'F2_d_all_thB' (s05-s11, scenarios, delta_plots, sigma_tot, pwa/scenarios/multipoles_figure.py,
radcorr) picks up the new main fit automatically.  The former main database is kept as the variant 'j_withHIGS'
(= OLD_L + A2 + HIGS, F1/F2 with and without the truncation covariance, F3 in the HIGS test); the HIgS data remain
as the test 'c_HIGSonly' / '*_HIGS_free*' (HIgS alone) and in the 'HIGS test' block (shift of alpha, beta when HIgS
is added, chi2 of the 16 HIgS points predicted from the main fit at n = 1 and with fitted normalisations,
parameter-goodness-of-fit of HIgS vs the main database).  The former variant 'i_noHIGS' coincides with 'd_all' and
is dropped.  'g_A2perbeam' = OLD_L + A2PB (no HIgS); 'e_noTAPS', 'f_GMPFpruned', the F4, forward-point,
truncation-model, order, parameter, spin and E_max-scan variants and the error budget are all built from the new
main database (ALL).

MAIN_DROP (since 2026-10-14) was determined once from the all-points fit (F2 + truncation covariance, two-pass
3-sigma rule) and is held fixed; this script re-derives it from the all-points fit at every run ('F2_d_allpts_thB',
now on the main database without HIgS) and stops if the list differs.  The drop is applied inside rows_for(), i.e.
to EVERY fit built through do_fit() -- the whole 'd_all' family (d_all, d_all_thB, the truncation-model, order,
parameter, spin, forward-point, F4 and E_max-scan variants, F2pi/F3*, fixsum, HIGS test) with UNCHANGED key names,
and also to every subset that contains the TAPS set (a_pre2015, g_A2perbeam, j_withHIGS, F1_TAPSonly,
pre2015+HIGS, pre2015+A2, LM14_type*; the GMPF-pruned sets never contained the point, which carries the flag
GMPF_out).  Subsets without TAPS (b_A2only, c_HIGSonly, e_noTAPS, A2+HIGS, bench_*) are unaffected.  The all-points
database is kept as the variant 'd_allpts' (keep_outlier=True), with and without the truncation covariance, for F1
and F2.  The 'h_3sigma' / 'h_3sigma_thBmain' variants run the 3-sigma search on the main database (check that
nothing else is flagged; 'removed' is then empty and the entry equals the d_all fit).  Scripts that select rows
directly (not through rows_for) must drop MAIN_DROP and the HIgS sets themselves.
"""
import os, sys, copy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
from common import *  # noqa
import theory as TH
import bchpt

df = load_db()
dfp = load_db(per_beamtime=True)
T = outputs(df, 'NNLO')
Tp = outputs(dfp, 'NNLO', tag='dbpb')
FWD_EMAX = 145.0 if PROFILE == 'v2' else 135.0   # v2: physical cusps -> 145 MeV usable
RESULTS = {}


def thcov_for(s, dfx, tag, **kw):
    return truncation_for(dfx, s, tag=tag, **kw)


def rows_for(dfx, sets, kind, Emax=150.0, excl=(), fwdEmax=FWD_EMAX, extra_mask=None, keep_outlier=False):
    """rows of the MAIN database for the fit: selection minus MAIN_DROP (unless keep_outlier), + forward points."""
    s = select(dfx, sets, Emax=Emax, exclude_flags=excl, extra_mask=extra_mask)
    if not keep_outlier and len(MAIN_DROP):
        s = s[~s['pid'].isin(MAIN_DROP)]
    if kind in ('F2', 'F3', 'F2pi'):
        s = pd.concat([s, select(dfx, FWD, Emax=fwdEmax)])
    return s


def do_fit(name, sets, kind='F1', Emax=150.0, excl=(), P0=PREF_NNLO, dfx=None, TT=None, thcov=False,
           fwd_mode='nuisance', fwdEmax=FWD_EMAX, baldin=BALDIN, extra_mask=None, Toffset=None, store=True,
           decompose=True, drop_pids=(), thopts=None, keep_outlier=False, extra_cons=(), norm_free=()):
    dfx = df if dfx is None else dfx
    TT = T if TT is None else TT
    s = rows_for(dfx, sets, kind, Emax, excl, fwdEmax, extra_mask, keep_outlier=keep_outlier)
    if len(drop_pids):
        s = s[~s['pid'].isin(drop_pids)]
    if len(norm_free):                                   # normalisations of these sets free (prior width 1000 %)
        s = s.copy(); s.loc[s['dataset'].isin(norm_free), 'norm'] = 10.0
    free = ['alphaE1', 'betaM1']
    cons = []
    if kind in ('F1', 'F3b', 'F3c') and baldin is not None:
        cons.append((C_BALDIN, baldin[0], baldin[1], 'Baldin'))
    if kind == 'F3':
        free += ['gM1M1', 'gE1E1']
    if kind in ('F3b', 'F2pi'):
        free += ['gpi@g0']
    if kind == 'F3c':
        free += ['gM1M1']
    cons += list(extra_cons)
    C_th = None
    if thcov:
        C_th, _ = thcov_for(s, dfx, 'db' if dfx is df else 'dbpb', **(thopts or {}))
    F = build_fit(s, TT, P0, free, constraints=cons, fwd_mode=fwd_mode, thcov=C_th, Toffset=Toffset)
    R = F.run(decompose=decompose)
    out = dict(name=name, kind=kind, N=len(s), Emax=Emax, sets=list(dict.fromkeys(s['dataset'])), free=free,
               P=R.P, dP=R.dP, covP=R.covP, chi2=R.chi2, dof=R.dof, p=R.pvalue, chi2_data=R.chi2_data,
               chi2_norm=R.chi2_norm, chi2_con=R.chi2_con, chi2_shift=R.chi2_shift,
               perset=R.perset, g0=gamma0(R.P), gpi=gammapi(R.P), pids=list(s['pid']),
               resid=R.resid, model=R.model, contrib=R.contrib, lam_shift=R.lam_shift,
               lam_shift_err=np.sqrt(np.clip(np.diag(R.cov)[-len(R.lam_shift):], 0, None)) if len(R.lam_shift) else [],
               sig_th=np.sqrt(np.diag(C_th)) if C_th is not None else np.zeros(len(s)))
    if decompose:
        out['dP_stat'] = np.sqrt(np.clip(np.diag(R.covP_stat), 0, None))
        out['dP_statnorm'] = np.sqrt(np.clip(np.diag(R.covP_statnorm), 0, None))
        if kind in ('F1', 'F3b', 'F3c') and baldin is not None:
            # Baldin part: refit with alpha+beta fixed (sigma -> 1e-4)
            F2_ = build_fit(s, TT, P0, free, constraints=[(C_BALDIN, baldin[0], 1e-4, 'Baldin')], fwd_mode=fwd_mode,
                            thcov=C_th, Toffset=Toffset)
            R2 = F2_.run(decompose=False)
            out['dP_fixsum'] = R2.dP
            out['dP_baldin'] = np.sqrt(np.clip(R.dP ** 2 - R2.dP ** 2, 0, None))
        else:
            out['dP_fwd'] = np.sqrt(np.clip(R.dP ** 2 - out['dP_statnorm'] ** 2, 0, None))
        out['dP_norm'] = np.sqrt(np.clip(out['dP_statnorm'] ** 2 - out['dP_stat'] ** 2, 0, None))
    out['rho_ab'] = R.covP[0, 1] / np.sqrt(R.covP[0, 0] * R.covP[1, 1])
    out['apb'] = R.P[0] + R.P[1]; out['dapb'] = float(np.sqrt(R.covP[0, 0] + R.covP[1, 1] + 2 * R.covP[0, 1]))
    out['amb'] = R.P[0] - R.P[1]; out['damb'] = float(np.sqrt(R.covP[0, 0] + R.covP[1, 1] - 2 * R.covP[0, 1]))
    J0 = np.array([0, 0, -1, -1, -1, -1.]); Jp = np.array([0, 0, -1, 1, -1, 1.])
    out['dg0'] = float(np.sqrt(J0 @ R.covP @ J0)); out['dgpi'] = float(np.sqrt(Jp @ R.covP @ Jp))
    if store:
        RESULTS[name] = out
    print('%-28s N=%3d  a=%6.2f(%.2f) b=%5.2f(%.2f)%s  chi2/dof=%6.1f/%d=%.2f p=%.2f' % (
        name, len(s), R.P[0], R.dP[0], R.P[1], R.dP[1],
        ('  gMM=%.2f(%.2f) gEE=%.2f(%.2f) g0=%.2f gpi=%.2f' % (R.P[3], R.dP[3], R.P[2], R.dP[2], gamma0(R.P), gammapi(R.P))
         if kind.startswith('F3') else ''), R.chi2, R.dof, R.chi2 / max(R.dof, 1), R.pvalue))
    return out


# 'd_all' = MAIN database = pre-2015 world data + A2 2022 (minus MAIN_DROP, + forward points for F2/F3);
# HIgS 2022 excluded since 2026-10-17 (author decision).  ALLH = the former main database (with HIgS) = 'j_withHIGS'.
ALL = list(MAIN_L_SETS)                      # = OLD_L + A2 (lib/common.py)
ALLH = list(MAIN_L_SETS_WITH_HIGS)           # = OLD_L + A2 + HIGS (former main, variant j)
assert ALL == OLD_L + A2 and ALLH == OLD_L + A2 + HIGS
NOTAPS = [k for k in ALL if k != 'OlmosdeLeon2001']
LEGACY = ['Pugh1957', 'Oxley1958', 'Hyman1959', 'Goldansky1960', 'Bernardini1960']   # bremsstrahlung data of 1957-1960
NOLEGACY = [k for k in ALL if k not in LEGACY]
VARIANTS = {
    'a_pre2015': dict(sets=OLD_L),
    'b_A2only': dict(sets=A2),
    'c_HIGSonly': dict(sets=HIGS),                            # test: HIgS alone (not part of the main database)
    'd_all': dict(sets=ALL),                                  # MAIN database = ALL minus MAIN_DROP (rows_for)
    'd_allpts': dict(sets=ALL, keep_outlier=True),            # (d') all points incl. the 3-sigma outlier
    'e_noTAPS': dict(sets=NOTAPS),
    'f_GMPFpruned': dict(sets=ALL, excl=('GMPF_out',)),
    'g_A2perbeam': dict(sets=OLD_L + A2PB, dfx='pb'),         # main database with the A2 data per beam time
    'j_withHIGS': dict(sets=ALLH),                            # (j) main database + HIgS 2022 (the former main)
    'k_noLegacy': dict(sets=NOLEGACY),                        # (k) main database without the data prior to 1961
}


def run_variant(kind, vname, v, **kw):
    dfx = dfp if v.get('dfx') == 'pb' else df
    TT = Tp if v.get('dfx') == 'pb' else T
    return do_fit('%s_%s' % (kind, vname), v['sets'], kind, excl=v.get('excl', ()), dfx=dfx, TT=TT,
                  keep_outlier=v.get('keep_outlier', False), **kw)


def flagged(res, nsig=3.0, with_th=False):
    """pids of the non-forward points with |pull| > nsig w.r.t. the fit `res`; pull = r / err_pt (with_th=False) or
    r / sqrt(err_pt^2 + sigma_th^2) (with_th=True).  Also returns the |pull| > 2.5 dictionary."""
    s = df.set_index('pid').loc[res['pids']]
    isf = s['dataset'].isin(FWD).values
    den = np.sqrt(s['err_pt'].values ** 2 + (np.array(res['sig_th']) ** 2 if with_th else 0.0))
    pull = np.array(res['resid']) / den
    bad = [p for p, x, f in zip(res['pids'], pull, isf) if abs(x) > nsig and not f]
    big = {p: float(x) for p, x, f in zip(res['pids'], pull, isf) if abs(x) > 2.5 and not f}
    return bad, big


def outlier_search(base, kind, nsig=3.0, with_th=False, keep_outlier=False, name=None):
    """two-pass 3-sigma rule starting from the fit `base`: flag, refit without the flagged points, flag again.
    Returns (last fit (a copy; = base if nothing was flagged), removed list, [|pull|>2.5 dict per pass])."""
    removed, hist, cur = [], [], base
    for it in range(2):
        bad, big = flagged(cur, nsig, with_th)
        hist.append(big)
        new = [b for b in bad if b not in removed]
        if not new:
            break
        removed += new
        cur = do_fit(name or '%s_outlier_search' % kind, ALL, kind, drop_pids=removed, thcov=with_th,
                     keep_outlier=keep_outlier, store=name is not None)
    return dict(cur), removed, hist


def outlier_variant(kind):
    """(h): 3-sigma rule (pull = r/err_pt) applied to the MAIN database, starting from the (d) fit; two passes.
    With MAIN_DROP already removed this is a check that nothing else is flagged ('removed' is then empty and the
    entry is a copy of the d_all fit)."""
    cur, removed, hist = outlier_search(RESULTS['%s_d_all' % kind], kind, with_th=False, name='%s_h_3sigma' % kind)
    cur['name'] = '%s_h_3sigma' % kind; cur['removed'] = removed; cur['pulls_gt2.5'] = hist
    RESULTS['%s_h_3sigma' % kind] = cur
    return removed


def outlier_variant_th(kind, nsig=3.0):
    """(h) with truncation covariance: pull_i = r_i / sqrt(err_pt^2 + sigma_th^2), two passes, on the MAIN database."""
    cur, removed, hist = outlier_search(RESULTS['%s_d_all_thB' % kind], kind, nsig, with_th=True,
                                        name='%s_h_3sigma_thBmain' % kind)
    cur['name'] = '%s_h_3sigma_thBmain' % kind; cur['removed'] = removed; cur['pulls_gt2.5'] = hist
    RESULTS['%s_h_3sigma_thBmain' % kind] = cur
    print('   outliers (%s, with trunc. cov., main database):' % kind, removed or 'none', ' |pull|>2.5 in pass 1:', hist[0])
    return removed


def verify_main_drop(nsig=3.0):
    """Determine the outlier list from the ALL-POINTS main-configuration fit (F2 + trunc. cov., two-pass 3-sigma
    rule, exactly the procedure of outlier_variant_th) and check that it equals the fixed MAIN_DROP."""
    base = RESULTS['F2_d_allpts_thB']
    cur, removed, hist = outlier_search(base, 'F2', nsig, with_th=True, keep_outlier=True)
    print('=== MAIN_DROP: 3-sigma outliers of the all-points fit F2_d_allpts_thB (N=%d): %s' % (base['N'], removed))
    print('    |pull| > 2.5 in pass 1: %s' % hist[0])
    for p in removed:
        r = df.set_index('pid').loc[p]
        print('    %s: %s  E = %.2f MeV  theta_lab = %.1f deg  value = %.3f(%.3f)' % (p, r['dataset'], r['E'], r['theta_lab'], r['value'], r['err_pt']))
    print('    MAIN_DROP (held fixed, lib/common.py) = %s' % (MAIN_DROP,))
    if sorted(removed) != sorted(MAIN_DROP):
        raise RuntimeError('MAIN_DROP %s differs from the outliers of the all-points fit %s: update MAIN_DROP in '
                           'lib/common.py (and the notes) deliberately, then rerun.' % (MAIN_DROP, removed))
    RESULTS['MAIN_DROP'] = dict(pids=list(MAIN_DROP), determined_from='F2_d_allpts_thB', removed=removed, pulls_gt2_5=hist)


if __name__ == '__main__':
    print('=== Mornacchi benchmark (A2 per beamtime, Baldin 13.8(4), spin pols = Paudyal 2020) ===')
    P0 = PREF_NNLO.copy(); P0[2:] = PAUDYAL
    do_fit('bench_Mornacchi', A2PB, 'F1', P0=P0, dfx=dfp, TT=Tp, baldin=(13.8, 0.4))
    do_fit('bench_Mornacchi_NNLOspin', A2PB, 'F1', dfx=dfp, TT=Tp, baldin=(13.8, 0.4))
    print('=== Lensky-McGovern 2014 type: pre-2015, GMPF pruning, E<=170 (Hallin<=150), Baldin 13.8(4) ===')
    do_fit('LM14_type', OLD_L, 'F1', Emax=170, excl=('GMPF_out',), baldin=(13.8, 0.4))
    do_fit('LM14_type_B14', OLD_L, 'F1', Emax=170, excl=('GMPF_out',))

    print('=== MAIN database = pre-2015 world data + A2 2022 (HIgS 2022 excluded since 2026-10-17) minus MAIN_DROP =', MAIN_DROP,
          '(all-points variant: d_allpts; with HIgS: j_withHIGS) ===')
    print('    angular sets:', ALL)
    for kind in ['F1', 'F2']:
        print('=== %s database variants ===' % kind)
        for vn, v in VARIANTS.items():
            run_variant(kind, vn, v)
        rem = outlier_variant(kind)
        print('   3-sigma outliers removed (%s, main database):' % kind, rem or 'none')
        do_fit('%s_h_3sigma_thB' % kind, ALL, kind, drop_pids=rem, thcov=True)

    print('=== F3 variants (spin polarisabilities) ===')
    for vn in ['d_all', 'f_GMPFpruned', 'a_pre2015', 'b_A2only']:
        run_variant('F3', vn, VARIANTS[vn])
    # F3g0: forward CS(0) only + gamma_0 = GHP_G0 penalty (2026: -1.025(102)), gM1M1 and gE1E1 free
    for vn in ['d_all', 'f_GMPFpruned']:
        for th in [False, True]:
            s_ = select(df, VARIANTS[vn]['sets'], exclude_flags=VARIANTS[vn].get('excl', ()))
            s_ = s_[~s_['pid'].isin(MAIN_DROP)]                      # main database (direct selection, not rows_for)
            s_ = pd.concat([s_, select(df, ['GHP_fwd_CS'], Emax=FWD_EMAX)])
            C_th = truncation_for(df, s_)[0] if th else None
            F = build_fit(s_, T, PREF_NNLO, ['alphaE1', 'betaM1', 'gM1M1', 'gE1E1'],
                          constraints=[(C_GAMMA0, GHP_G0[0], GHP_G0[1], 'gamma0')], thcov=C_th)
            R = F.run()
            nm = 'F3g0_%s%s' % (vn, '_thB' if th else '')
            RESULTS[nm] = dict(name=nm, kind='F3g0', N=len(s_), P=R.P, dP=R.dP, covP=R.covP, chi2=R.chi2, dof=R.dof,
                               p=R.pvalue, perset=R.perset, g0=gamma0(R.P), gpi=gammapi(R.P), pids=list(s_['pid']),
                               dP_stat=np.sqrt(np.diag(R.covP_stat)), dP_statnorm=np.sqrt(np.diag(R.covP_statnorm)),
                               resid=R.resid, model=R.model, contrib=R.contrib, chi2_data=R.chi2_data,
                               rho_ab=R.covP[0, 1] / np.sqrt(R.covP[0, 0] * R.covP[1, 1]))
            print('%-28s N=%3d  a=%6.2f(%.2f) b=%5.2f(%.2f)  gMM=%.2f(%.2f) gEE=%.2f(%.2f) g0=%.2f gpi=%.2f chi2/dof=%.1f/%d p=%.2f'
                  % (nm, len(s_), R.P[0], R.dP[0], R.P[1], R.dP[1], R.P[3], R.dP[3], R.P[2], R.dP[2], gamma0(R.P),
                     gammapi(R.P), R.chi2, R.dof, R.pvalue))
    for kind in ['F3b', 'F3c']:
        for vn in ['d_all', 'f_GMPFpruned']:
            run_variant(kind, vn, VARIANTS[vn])
            run_variant(kind, vn + '_thB', VARIANTS[vn], thcov=True)

    print('=== F4: E <= 170 MeV ===')
    for kind in ['F1', 'F2']:
        do_fit('%s_F4_all' % kind, ALL, kind, Emax=170, excl=(), extra_mask=~((df.dataset == 'Hallin1993_TabII') & (df.E > 150)))
        do_fit('%s_F4_all_Hallin170' % kind, ALL, kind, Emax=170)
        do_fit('%s_F4_pruned' % kind, ALL, kind, Emax=170, excl=('GMPF_out',))

    print('=== forward-point model variants (F2, d) ===')
    do_fit('F2_d_fwd_uncorr', ALL, 'F2', fwd_mode='uncorr')
    do_fit('F2_d_fwd_150', ALL, 'F2', fwdEmax=150.0)
    do_fit('F2_d_fwdonly_CS', ALL, 'F2', extra_mask=None, drop_pids=[p for p in df.pid if p.startswith('GHP_fwd_S2Z')])

    print('=== theory (truncation) covariance included: scheme B (LEC-matched, default) and A (naive) ===')
    for kind in ['F1', 'F2']:
        do_fit('%s_d_all_thB' % kind, ALL, kind, thcov=True)
        do_fit('%s_d_allpts_thB' % kind, ALL, kind, thcov=True, keep_outlier=True)   # (d') all points incl. outlier
        do_fit('%s_f_GMPFpruned_thB' % kind, ALL, kind, excl=('GMPF_out',), thcov=True)
        do_fit('%s_e_noTAPS_thB' % kind, NOTAPS, kind, thcov=True)
        do_fit('%s_g_A2perbeam_thB' % kind, OLD_L + A2PB, kind, dfx=dfp, TT=Tp, thcov=True)
        do_fit('%s_j_withHIGS_thB' % kind, ALLH, kind, thcov=True)          # (j) with HIgS 2022 = the former main fit
        do_fit('%s_k_noLegacy_thB' % kind, NOLEGACY, kind, thcov=True)      # (k) without the bremsstrahlung data of 1957-1960
        do_fit('%s_d_all_thB_diag' % kind, ALL, kind, thcov=True, thopts=dict(corr=False))
        do_fit('%s_d_all_thB_L600' % kind, ALL, kind, thcov=True, thopts=dict(Lam=600.0))
        do_fit('%s_d_all_thB_L900' % kind, ALL, kind, thcov=True, thopts=dict(Lam=900.0))
        do_fit('%s_d_all_thB_fixed048' % kind, ALL, kind, thcov=True, thopts=dict(mode='fixed'))
        do_fit('%s_d_all_thB_max' % kind, ALL, kind, thcov=True, thopts=dict(cbar='max'))
        do_fit('%s_d_all_thB_l50' % kind, ALL, kind, thcov=True, thopts=dict(lE=50.0, lth=50.0))
        do_fit('%s_d_all_thB_l12' % kind, ALL, kind, thcov=True, thopts=dict(lE=12.0, lth=25.0))
        do_fit('%s_d_all_thB_l100' % kind, ALL, kind, thcov=True, thopts=dict(lE=100.0, lth=100.0))
        do_fit('%s_d_all_thA' % kind, ALL, kind, thcov=True, thopts=dict(scheme='A'))
        do_fit('%s_d_all_thA_diag' % kind, ALL, kind, thcov=True, thopts=dict(scheme='A', corr=False))
    for vn in ['d_all', 'f_GMPFpruned']:
        run_variant('F3', vn + '_thB', VARIANTS[vn], thcov=True)
    do_fit('F2_d_fwd_uncorr_thB', ALL, 'F2', fwd_mode='uncorr', thcov=True)
    do_fit('F2_d_fwdonly_CS_thB', ALL, 'F2', thcov=True, drop_pids=[p for p in df.pid if p.startswith('GHP_fwd_S2Z')])
    do_fit('F1_F4_all_thB', ALL, 'F1', Emax=170, thcov=True, extra_mask=~((df.dataset == 'Hallin1993_TabII') & (df.E > 150)))
    do_fit('F2_F4_all_thB', ALL, 'F2', Emax=170, thcov=True, extra_mask=~((df.dataset == 'Hallin1993_TabII') & (df.E > 150)))
    for vn in ['a_pre2015', 'b_A2only', 'c_HIGSonly']:
        for kind in ['F1', 'F2']:
            run_variant(kind, vn + '_thB', VARIANTS[vn], thcov=True)

    print('=== order-by-order fits (base amplitude NLO + O(p^4) LECs) and parameter variations ===')
    VAR = ['NLO', 'NNLO_noPiD', 'gM2.87', 'gM3.07'] + (['PDGconst', 'hA2.84', 'hA2.91', 'singlemass'] if PROFILE == 'v2'
                                                      else ['gA1.2754', 'fpi92.07', 'PDGconst'])
    for cfg in VAR:
        TT = outputs(df, cfg)
        # total static polarisabilities P = static(base of cfg) + d ; driver expects Pref + d
        off = base_static(TT, cfg) - TT.Pref
        P0 = base_static(TT, cfg).copy()
        P0[2:] = PREF_NNLO[2:] if cfg not in ('NLO',) else P0[2:]
        if cfg.startswith('gM') or cfg in ('gA1.2754', 'fpi92.07', 'PDGconst', 'hA2.84', 'hA2.91', 'singlemass'):
            P0 = TT.Pref.copy()       # complete NNLO with modified parameter: spin pols = its own NNLO values
            off = np.zeros(6)
        for kind in ['F1', 'F2']:
            do_fit('%s_d_all_%s' % (kind, cfg), ALL, kind, TT=TT, P0=P0, Toffset=off, decompose=False)
            if cfg == 'NLO':   # NLO amplitude, but all four spin pols matched to NNLO via contact terms
                P0m = P0.copy(); P0m[2:] = PREF_NNLO[2:]
                do_fit('%s_d_all_NLO_spinNNLO' % kind, ALL, kind, TT=TT, P0=P0m, Toffset=off, decompose=False)
                do_fit('%s_d_all_NLO_spinNNLO_thB' % kind, ALL, kind, TT=TT, P0=P0m, Toffset=off, decompose=False,
                       thcov=True)
            do_fit('%s_d_all_%s_thB' % (kind, cfg), ALL, kind, TT=TT, P0=P0, Toffset=off, decompose=False, thcov=True)
    # spin polarisabilities at Paudyal / Mornacchi-2022 bootstrap values instead of NNLO
    for kind in ['F1', 'F2']:
        P0 = PREF_NNLO.copy(); P0[2:] = PAUDYAL
        do_fit('%s_d_all_spinPaudyal' % kind, ALL, kind, P0=P0, decompose=False)
        P0 = PREF_NNLO.copy(); P0[2:] = [-3.0, 3.7, -1.2, 2.0]
        do_fit('%s_d_all_spinMRPP22' % kind, ALL, kind, P0=P0, decompose=False)
    do_fit('F1_d_all_Baldin13.8', ALL, 'F1', baldin=(13.8, 0.4), decompose=False)
    print('=== MAIN configuration: F2 (forward points, nuisance model) + truncation covariance (scheme B) ===')
    verify_main_drop()                       # MAIN_DROP re-derived from the all-points fit; stops if it differs
    for kind in ['F2', 'F1']:
        outlier_variant_th(kind)             # 3-sigma search on the MAIN database: check that nothing else is flagged
    do_fit('F2_d_fwd150_thB', ALL, 'F2', fwdEmax=150.0, thcov=True)
    do_fit('F2_d_fwd150_uncorr_thB', ALL, 'F2', fwdEmax=150.0, thcov=True, fwd_mode='uncorr')
    do_fit('F2_d_fwd135_thB', ALL, 'F2', fwdEmax=135.0, thcov=True)
    do_fit('F2_d_fwd135_uncorr_thB', ALL, 'F2', fwdEmax=135.0, thcov=True, fwd_mode='uncorr')
    do_fit('F2_d_fwd135', ALL, 'F2', fwdEmax=135.0)
    do_fit('F2_d_fwd150', ALL, 'F2', fwdEmax=150.0)
    do_fit('F2pi_d_all_thB', ALL, 'F2pi', thcov=True)
    do_fit('F2pi_d_all', ALL, 'F2pi')
    do_fit('F2pi_f_GMPFpruned_thB', ALL, 'F2pi', excl=('GMPF_out',), thcov=True)
    for kind in ['F2', 'F1']:
        P0 = PREF_NNLO.copy(); P0[2:] = PAUDYAL
        do_fit('%s_d_all_spinPaudyal_thB' % kind, ALL, kind, P0=P0, thcov=True, decompose=False)
        P0 = PREF_NNLO.copy(); P0[2:] = [-3.0, 3.7, -1.2, 2.0]
        do_fit('%s_d_all_spinMRPP22_thB' % kind, ALL, kind, P0=P0, thcov=True, decompose=False)
    do_fit('F2_F4_pruned_thB', ALL, 'F2', Emax=170, excl=('GMPF_out',), thcov=True)
    do_fit('F2_F4_all_Hallin170_thB', ALL, 'F2', Emax=170, thcov=True)

    print('=== HIGS test: MAIN database (no HIgS) vs main + HIgS 2022 (Li:2022vnz, variant j), HIgS alone, PG test, HIgS pulls ===')
    # Since 2026-10-17 the HIgS data are NOT part of the main database: 'main' = d_all (OLD_L + A2 minus MAIN_DROP,
    # + forward points), 'withHIGS' = j_withHIGS (the former main).  shift = P(withHIGS) - P(main) is the shift of
    # the result when the HIgS data are ADDED; the 16 HIgS points are predicted from the main fit (n = 1 and n fitted
    # within the 2.24 % normalisation prior); PG = chi2(withHIGS) - chi2(HIgS alone) - chi2(main), dof 2 (alpha, beta).
    HT = {}
    from scipy.stats import chi2 as _c2, norm
    for th in [False, True]:
        sfx = '_thB' if th else ''
        for kind in ['F1', 'F2', 'F3']:
            r_main = RESULTS['%s_d_all%s' % (kind, sfx)]
            if '%s_j_withHIGS%s' % (kind, sfx) in RESULTS:
                r_with = RESULTS['%s_j_withHIGS%s' % (kind, sfx)]
            else:                                                    # F3 (with and without truncation covariance)
                r_with = run_variant(kind, 'j_withHIGS' + sfx, VARIANTS['j_withHIGS'], thcov=th)
            # HIGS alone, alpha and beta free, no sum constraint (the Baldin/forward information is attributed to 'rest');
            # the same fit for every kind (kept under the three names for the consumers, e.g. s06: 'F1_HIGS_free')
            r_h = do_fit('%s_HIGS_free%s' % (kind, sfx), HIGS, 'F1', baldin=(14.0, 1e3), thcov=th, decompose=False)
            pg = r_with['chi2'] - r_h['chi2'] - r_main['chi2']
            # the HIGS points predicted from the MAIN fit (no HIgS in the fit)
            sH = select(df, HIGS)
            idxH = T.subset(sH.pid.values)
            tH = T.predict_array(np.array(r_main['P']), idxH)
            pull1 = (sH.value.values - tH) / sH.err_pt.values
            nfit = {}
            pulln = np.empty(len(sH))
            for ds in HIGS:
                m = (sH.dataset == ds).values
                y, e, t, dn = sH.value.values[m], sH.err_pt.values[m], tH[m], sH[m].norm.iloc[0]
                n = (np.sum(t * y / e ** 2) + 1 / dn ** 2) / (np.sum(t * t / e ** 2) + 1 / dn ** 2)
                nfit[ds] = n
                pulln[m] = (y - n * t) / e
            chiH1 = float(np.sum(pull1 ** 2)); chiHn = float(np.sum(pulln ** 2) + sum(((nfit[k] - 1) / 0.0224) ** 2 for k in HIGS))
            HT[kind + sfx] = dict(P_main=r_main['P'], dP_main=r_main['dP'], P_withHIGS=r_with['P'], dP_withHIGS=r_with['dP'],
                                  shift=np.array(r_with['P']) - np.array(r_main['P']),           # shift when HIgS is ADDED
                                  chi2_main=r_main['chi2'], dof_main=r_main['dof'], chi2_withHIGS=r_with['chi2'], dof_withHIGS=r_with['dof'],
                                  chi2_HIGSfree=r_h['chi2'], dof_HIGSfree=r_h['dof'], P_HIGSfree=r_h['P'], dP_HIGSfree=r_h['dP'],
                                  PG=pg, PG_dof=2, PG_p=float(_c2.sf(pg, 2)),
                                  HIGS_pids=list(sH.pid), HIGS_E=list(sH.E), HIGS_th=list(sH.theta_lab), HIGS_obs=list(sH.dobs),
                                  HIGS_y=list(sH.value), HIGS_err=list(sH.err_pt), HIGS_pred=tH, pull_n1=pull1, pull_nfit=pulln,
                                  n_HIGS_fit=nfit, chi2_HIGS_pred_n1=chiH1, chi2_HIGS_pred_nfit=chiHn,
                                  chi2_HIGS_in_withHIGS={k: r_with['perset'][k]['chi2'] for k in HIGS if k in r_with['perset']},
                                  n_HIGS_in_withHIGS={k: (r_with['perset'][k]['n'], r_with['perset'][k]['dn']) for k in HIGS if k in r_with['perset']},
                                  norms_main={k: (v['n'], v['dn']) for k, v in r_main['perset'].items()},
                                  norms_withHIGS={k: (v['n'], v['dn']) for k, v in r_with['perset'].items()})
            print('   %s%s: main a=%.2f b=%.2f; +HIgS a=%.2f b=%.2f (shift %+.2f, %+.2f); PG chi2 = %.2f - %.2f - %.2f = %.2f (dof 2, p = %.3f); '
                  'HIgS predicted from the main fit: chi2 = %.1f/16 (n=1), %.1f (n fitted %s)'
                  % (kind, sfx, r_main['P'][0], r_main['P'][1], r_with['P'][0], r_with['P'][1], r_with['P'][0] - r_main['P'][0],
                     r_with['P'][1] - r_main['P'][1], r_with['chi2'], r_h['chi2'], r_main['chi2'], pg, _c2.sf(pg, 2), chiH1, chiHn,
                     ', '.join('%.3f' % v for v in nfit.values())))
        # fixed-sum PG (alpha+beta = BALDIN[0] exactly in all three fits; one parameter alpha-beta; dof 1)
        rw = do_fit('F1_j_withHIGS_fixsum' + sfx, ALLH, 'F1', baldin=(BALDIN[0], 1e-4), thcov=th, decompose=False)
        rm = do_fit('F1_d_all_fixsum' + sfx, ALL, 'F1', baldin=(BALDIN[0], 1e-4), thcov=th, decompose=False)
        rh = do_fit('F1_HIGS_fixsum' + sfx, HIGS, 'F1', baldin=(BALDIN[0], 1e-4), thcov=th, decompose=False)
        pg = rw['chi2'] - rm['chi2'] - rh['chi2']
        HT['F1_fixsum' + sfx] = dict(PG=pg, PG_dof=1, PG_p=float(_c2.sf(pg, 1)), chi2_withHIGS=rw['chi2'], chi2_main=rm['chi2'],
                                     chi2_HIGS=rh['chi2'], amb_withHIGS=rw['P'][0] - rw['P'][1], amb_main=rm['P'][0] - rm['P'][1],
                                     amb_HIGS=rh['P'][0] - rh['P'][1])
        print('   F1 fixed sum%s: PG = %.2f (dof 1, p = %.3f); alpha-beta: with HIgS %.2f, main %.2f, HIgS alone %.2f'
              % (sfx, pg, _c2.sf(pg, 1), rw['P'][0] - rw['P'][1], rm['P'][0] - rm['P'][1], rh['P'][0] - rh['P'][1]))
    RESULTS['HIGS_test'] = HT

    print('=== single-observable / subset variants (with and without theory covariance) ===')
    for th in [False, True]:
        sfx = '_thB' if th else ''
        do_fit('F1_A2_CSonly' + sfx, ['A2_2022_CS'], 'F1', thcov=th)
        do_fit('F1_A2_S3only' + sfx, ['A2_2022_S3'], 'F1', thcov=th)
        do_fit('F1_TAPSonly' + sfx, ['OlmosdeLeon2001'], 'F1', thcov=th)
        do_fit('F1_pre2015+HIGS' + sfx, OLD_L + HIGS, 'F1', thcov=th)
        do_fit('F1_pre2015+A2' + sfx, OLD_L + A2, 'F1', thcov=th)
        do_fit('F1_A2+HIGS' + sfx, A2 + HIGS, 'F1', thcov=th)
        do_fit('F1_d_free' + sfx, ALL, 'F1', baldin=None, thcov=th, decompose=False)   # no sum constraint (no dummy constraint: dof = N - 2 - norms)
    # HIgS normalisations free (Sec. IV B): the HIgS data alone + forward points, and the main database + HIgS
    do_fit('F2_c_HIGSonly_freenorm_thB', HIGS, 'F2', thcov=True, norm_free=HIGS, decompose=False)
    do_fit('F2_j_withHIGS_freenorm_thB', ALLH, 'F2', thcov=True, norm_free=HIGS, decompose=False)
    # the NNLO fit with and without the forward points: the angular data alone (main database, no forward points, no sum
    # constraint), the forward points alone (they fix alpha+beta only; alpha-beta held at the main fit, its value does not
    # matter), and the parameter-goodness-of-fit test of the two (one common direction, alpha+beta: dof 1)
    rm = RESULTS['F2_d_all_thB']
    ra = do_fit('F2_d_nofwd_thB', ALL, 'F2', thcov=True, fwdEmax=0.0, decompose=False)
    rf = do_fit('F2_fwdonly_thB', [], 'F2', thcov=True, decompose=False,
                extra_cons=[(np.array([1, -1, 0, 0, 0, 0.]), rm['amb'], 1e-3, 'alpha-beta fixed')])
    pg, pgdof = rm['chi2'] - ra['chi2'] - rf['chi2'], rm['dof'] - ra['dof'] - rf['dof']
    RESULTS['FWD_test'] = dict(PG=pg, PG_dof=pgdof, PG_p=float(_c2.sf(pg, pgdof)), PG_sigma=float(norm.isf(_c2.sf(pg, pgdof) / 2)),
                               apb_angular=(ra['apb'], ra['dapb']), apb_forward=(rf['apb'], rf['dapb']), apb_joint=(rm['apb'], rm['dapb']),
                               chi2_joint=rm['chi2'], chi2_angular=ra['chi2'], chi2_forward=rf['chi2'])
    print('   angular data vs forward points: alpha+beta %.2f(%.2f) vs %.2f(%.2f), joint %.2f(%.2f); PG = %.2f - %.2f - %.2f = %.2f '
          '(dof %d, p = %.3f, %.1f sigma)' % (ra['apb'], ra['dapb'], rf['apb'], rf['dapb'], rm['apb'], rm['dapb'], rm['chi2'], ra['chi2'],
                                               rf['chi2'], pg, pgdof, _c2.sf(pg, pgdof), RESULTS['FWD_test']['PG_sigma']))
    print('=== energy-cut scans ===')
    for th in [False, True]:
        for kind in ['F1', 'F2']:
            for Emax in [90, 100, 110, 120, 130, 140, 150, 160, 170]:
                em = ~((df.dataset == 'Hallin1993_TabII') & (df.E > 150))
                do_fit('scan_%s_%s_E%d' % (kind, 'thB' if th else 'noth', Emax), ALL, kind, Emax=Emax, thcov=th,
                       extra_mask=em, decompose=False)
    # alpha+beta: combination of Baldin and forward points is NOT done (double counting) -- by construction.
    save_json(RESULTS, 'fits_L.json')
    print('saved', len(RESULTS), 'fits')
