"""pi0ff_fits.py -- impact of the pi0 gamma gamma transition form factor of the proton triangle loop on the pi0-pole
vertex (work/code_v2, &v2 pi0ff=t, CHANGES.md Sect. F; derivation in work/notes/pi0_formfactor.md) on the fits:
 (1) NNLO fit of the paper's main configuration F2_d_all_thB (main database, forward points as data, truncation
     covariance) and the Baldin fit F1_d_all_thB, with and without pi0ff.  The truncation covariance (scheme B,
     LEC-matched, built from the order-by-order outputs Born+pi0 / NLO / NNLO) is recomputed with pi0ff=t, since the
     pi0 pole is part of the Born+pi0 reference; a variant keeps the pi0ff=f covariance to isolate that step.
 (2) Delta-region spin-polarisability fit of the paper (s05_delta.py key S_A2_LEGS_gpi: A2 Sigma_2x, Sigma_2z and
     LEGS Sigma_3 below the 2pi threshold, Gaussian priors on alpha+beta, alpha-beta [from the main fit], gamma_0,
     gamma_pi) with pi0ff=t (alpha-beta prior from the pi0ff main fit; variant with the paper's prior).
 (3) chi^2 per point of the Delta-region data sets (paper Table 'delta', columns n = 1 and n fitted) with pi0ff=t.
 (5) relative change of dsigma/dOmega_lab and change of Sigma_3 at 100, 140, 150, 250, 300 MeV.
Machinery: work/fits/lib (bchpt.run / outputs / build_fit) and scripts/s03_fits_L.do_fit (imported; its thcov_for is
temporarily replaced by the pi0ff version).  Cache: new entries work/fits/cache/pi0ff_* (in.dat carries pi0ff=t, so the
hash differs from every pi0ff=f entry); nothing of the pi0ff=f pipeline (results/, cache) is overwritten.
Run:  python3 pi0ff_fits.py   (~2 min)
Outputs (this directory): pi0ff_fits.json, pi0ff_fits.md"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
FITS = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(FITS, 'scripts')); sys.path.insert(0, os.path.join(FITS, 'lib'))
os.chdir(FITS)
import numpy as np
import pandas as pd
import bchpt
import theory as TH
from common import (load_db, outputs, select, build_fit, FWD, BALDIN, GHP_G0, C_BALDIN, C_GAMMA0, gamma0, gammapi,
                    PREF_NNLO, MAIN_DROP, RES)
import s03_fits_L as S          # import: runs/caches outputs(df,'NNLO') with the current executable

FF = 'pimassphys=t,nmassphys=t,pi0ff=t'          # main v2 settings + form factor
bchpt.CONFIGS.update({'NNLO_pi0ff': dict(v2x=FF), 'NLO_pi0ff': dict(order=3.0, v2x=FF),
                      'BornPi0_pi0ff': dict(oneNloop=False, oneDloop=False, delta_pole=False, v2x=FF)})
TAG = 'pi0ff'
df = S.df
T0 = S.T                                         # pi0ff=f (pipeline amplitude)
Tff = outputs(df, 'NNLO_pi0ff', tag=TAG)         # pi0ff=t
assert np.allclose(Tff.Pref, T0.Pref, atol=1e-6), (Tff.Pref, T0.Pref)
assert list(Tff.ids) == list(T0.ids)
FL = json.load(open(os.path.join(RES, 'fits_L.json')))
DR = json.load(open(os.path.join(RES, 'delta_region.json')))
OUT = {}


def conv(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))


# ---------------------------------------------------------------- truncation covariance with pi0ff=t (copy of common.order_values/truncation_for)
def order_values_ff(dfx, s):
    T3 = outputs(dfx, 'NNLO_pi0ff', tag=TAG)
    idx = T3.subset(s['pid'].values)
    y3 = T3.value_file[idx]

    def matched(cfg):                           # LEC-matched: all orders with the NNLO static polarisabilities
        TT = outputs(dfx, cfg, tag=TAG)
        off = TT.Pref - TT.Pref                 # base_static(TT) = TT.Pref in the v2 profile
        return TT.predict_array(T3.Pref - off, idx)
    return matched('BornPi0_pi0ff'), matched('NLO_pi0ff'), y3


def truncation_for_ff(dfx, s, mode='sqrtQ', Lam=650.0, cbar='rms', corr=True, lE=25.0, lth=40.0):
    y0, y2, y3 = order_values_ff(dfx, s)
    tr = TH.truncation(y0, y2, y3, s['dobs'].values, s['E'].values, mode=mode, Lam=Lam, cbar=cbar)
    if corr:
        C = TH.covariance(tr['sigma'], s['dobs'].values, s['E'].values, s['theta_cm'].values, lE=lE, lth=lth)
    else:
        C = np.diag(tr['sigma'] ** 2)
    tr['y0'], tr['y2'], tr['y3'] = y0, y2, y3
    return C, tr


thcov_std = S.thcov_for


def thcov_ff(s, dfx, tag, **kw):
    assert dfx is df
    return truncation_for_ff(dfx, s, **kw)


def strip(r):
    keep = ['name', 'kind', 'N', 'free', 'P', 'dP', 'covP', 'chi2', 'dof', 'p', 'chi2_data', 'chi2_norm', 'chi2_con',
            'chi2_shift', 'perset', 'g0', 'gpi', 'apb', 'dapb', 'amb', 'damb', 'rho_ab', 'dP_stat', 'dP_statnorm',
            'dP_fwd', 'dP_baldin', 'dP_norm', 'lam_shift']
    return {k: r[k] for k in keep if k in r}


# ---------------------------------------------------------------- (1) NNLO fits
print('=== (1) NNLO fits, main configuration (F2 + truncation covariance) and Baldin fit (F1 + truncation covariance) ===')
fits = {}
for kind in ['F2', 'F1']:
    r0 = S.do_fit('%s_d_all_thB' % kind, S.ALL, kind, thcov=True, store=False)                 # pi0ff=f (= fits_L.json)
    ref = FL['%s_d_all_thB' % kind]
    print('    reproduction of results/fits_L.json: max |dP| = %.1e, dchi2 = %.1e' % (
        np.max(np.abs(np.array(r0['P']) - np.array(ref['P']))), r0['chi2'] - ref['chi2']))
    S.thcov_for = thcov_ff
    r1 = S.do_fit('%s_d_all_thB_pi0ff' % kind, S.ALL, kind, thcov=True, TT=Tff, store=False)   # pi0ff=t, covariance recomputed
    S.thcov_for = thcov_std
    r2 = S.do_fit('%s_d_all_thB_pi0ff_covF' % kind, S.ALL, kind, thcov=True, TT=Tff, store=False)   # pi0ff=t amplitude, pi0ff=f covariance
    r3 = S.do_fit('%s_d_all' % kind, S.ALL, kind, store=False)                                 # no truncation covariance
    r4 = S.do_fit('%s_d_all_pi0ff' % kind, S.ALL, kind, TT=Tff, store=False)
    for r in (r0, r1, r2, r3, r4):
        fits[r['name']] = strip(r)
    # truncation-covariance change itself
    s_ = S.rows_for(df, S.ALL, kind)
    C0, tr0 = thcov_std(s_, df, 'db'); C1, tr1 = thcov_ff(s_, df, 'db')
    sig0, sig1 = np.sqrt(np.diag(C0)), np.sqrt(np.diag(C1))
    fits['%s_trunc_sigma_change' % kind] = dict(max_rel=float(np.max(np.abs(sig1 / sig0 - 1))),
                                                mean_rel=float(np.mean(sig1 / sig0 - 1)),
                                                max_abs_CS_rel=float(np.max(np.abs(sig1 - sig0)[s_.dobs.values == 'CS'] / np.abs(tr0['y3'])[s_.dobs.values == 'CS'])))
    print('    truncation sigma_th with pi0ff: max rel. change %.2e, mean %.2e' % (
        fits['%s_trunc_sigma_change' % kind]['max_rel'], fits['%s_trunc_sigma_change' % kind]['mean_rel']))
OUT['fits'] = fits
OUT['fits_L_reference'] = {k: strip(FL[k]) for k in ['F2_d_all_thB', 'F1_d_all_thB', 'F2_d_all', 'F1_d_all']}

# ---------------------------------------------------------------- (2) Delta-region spin fit (copy of s05_delta.spin_fit, S_A2_LEGS_gpi)
E2PI = 320.7
GPI = (8.0, 1.8)
AL = ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_S3']


def spin_fit(name, Tx, F2, sets=AL, gpi=GPI, legs_Emax=E2PI):
    m = df.dataset.isin(sets) & (df.E > 150) & (df.E <= 350)
    m &= ~((df.dataset == 'LEGS2001_S3') & (df.E >= legs_Emax))
    s = df[m].copy()
    P0 = np.array(F2['P']).copy(); cP = np.array(F2['covP'])
    amb, damb = P0[0] - P0[1], float(np.sqrt(cP[0, 0] + cP[1, 1] - 2 * cP[0, 1]))
    cons = [(C_BALDIN, BALDIN[0], BALDIN[1], 'Baldin'), (np.array([1, -1, 0, 0, 0, 0.]), amb, damb, 'a-b(F2)'),
            (C_GAMMA0, GHP_G0[0], GHP_G0[1], 'gamma0'), (np.array([0, 0, -1, 1, -1, 1.]), gpi[0], gpi[1], 'gammapi')]
    free = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']
    idx = Tx.subset(s.pid.values)
    F = build_fit(s, Tx, P0, free, constraints=cons)
    R = F.run(decompose=False)
    J = np.array([0, 0, -1, 1, -1, 1.]); J0 = np.array([0, 0, -1, -1, -1, -1.])
    out = dict(name=name, sets=sets, N=len(s), P=R.P, dP=R.dP, covP=R.covP, chi2=R.chi2, dof=R.dof, p=R.pvalue,
               g0=gamma0(R.P), gpi=gammapi(R.P), dgpi=float(np.sqrt(J @ R.covP @ J)), dg0=float(np.sqrt(J0 @ R.covP @ J0)),
               perset=R.perset, chi2_data=R.chi2_data, chi2_con=R.chi2_con, prior_amb=(amb, damb))
    print('%-32s N=%3d gEE=%5.2f(%.2f) gMM=%5.2f(%.2f) gEM=%5.2f(%.2f) gME=%5.2f(%.2f) g0=%5.2f gpi=%5.2f(%.2f) a=%.2f b=%.2f chi2=%.2f/%d'
          % (name, len(s), R.P[2], R.dP[2], R.P[3], R.dP[3], R.P[4], R.dP[4], R.P[5], R.dP[5], gamma0(R.P), gammapi(R.P),
             out['dgpi'], R.P[0], R.P[1], R.chi2, R.dof))
    return out


print('=== (2) Delta-region spin-polarisability fit S_A2_LEGS_gpi ===')
sp = {}
sp['S_A2_LEGS_gpi'] = spin_fit('S_A2_LEGS_gpi', T0, fits['F2_d_all_thB'])                 # = delta_region.json
sp['S_A2_LEGS_gpi_pi0ff'] = spin_fit('S_A2_LEGS_gpi_pi0ff', Tff, fits['F2_d_all_thB_pi0ff'])
sp['S_A2_LEGS_gpi_pi0ff_priorF'] = spin_fit('S_A2_LEGS_gpi_pi0ff_priorF', Tff, fits['F2_d_all_thB'])
ref = DR['S_A2_LEGS_gpi']
print('    reproduction of results/delta_region.json: max |dP| = %.1e, dchi2 = %.1e' % (
    np.max(np.abs(np.array(sp['S_A2_LEGS_gpi']['P']) - np.array(ref['P']))), sp['S_A2_LEGS_gpi']['chi2'] - ref['chi2']))
OUT['spin'] = sp
OUT['spin_reference'] = {k: ref[k] for k in ['P', 'dP', 'chi2', 'dof', 'gpi', 'g0']}

# ---------------------------------------------------------------- (3) chi^2 per data set in the Delta region (copy of s05_delta.chi2_table)
D_ALL = ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_S3', 'LEGS2001_CS', 'Hunger1997_75', 'Hunger1997_90',
         'Wissmann1999', 'Camen2002', 'Hallin1993_TabII', 'Hallin1993_TabII_D', 'OlmosdeLeon2001', 'Gray1967',
         'Genzel1972', 'Baranov1966', 'Kolb2000', 'Nagashima1964']


def chi2_table(P, Tx):
    s = df[df.dataset.isin(D_ALL) & (df.E > 150) & (df.E <= 350)].copy()
    s.loc[s.dataset == 'Hallin1993_TabII', 'dataset'] = 'Hallin1993_TabII_D'
    idx = Tx.subset(s.pid.values)
    th = Tx.predict_array(np.asarray(P, float), idx)
    rows = []
    for ds in dict.fromkeys(s.dataset):
        m = (s.dataset == ds).values
        y, e, t = s.value.values[m], s.err_pt.values[m], th[m]
        c0 = float(np.sum(((y - t) / e) ** 2))
        dn = float(s[m].norm.iloc[0])
        if dn > 0:
            a = np.sum(t * t / e ** 2) + 1 / dn ** 2
            b = np.sum(t * y / e ** 2) + 1 / dn ** 2
            n = b / a
            c1 = float(np.sum(((y - n * t) / e) ** 2) + ((n - 1) / dn) ** 2)
        else:
            n, c1 = 1.0, c0
        rows.append(dict(set=ds, obs='/'.join(sorted(set(s[m].dobs))), N=int(m.sum()),
                         Erange='%.0f-%.0f' % (s[m].E.min(), s[m].E.max()), chi2_n1=c0, chi2_nfit=c1, n=n,
                         chi2N_n1=c0 / m.sum(), chi2N_nfit=c1 / m.sum(), floated=bool(dn > 0)))
        if ds == 'LEGS2001_S3':
            mb = m & (s.E.values < E2PI)
            c = float(np.sum(((s.value.values[mb] - th[mb]) / s.err_pt.values[mb]) ** 2))
            rows.append(dict(set='LEGS2001_S3 (E<2pi thr.)', obs='S3', N=int(mb.sum()), Erange='<%.0f' % E2PI,
                             chi2_n1=c, chi2_nfit=c, n=1.0, chi2N_n1=c / mb.sum(), chi2N_nfit=c / mb.sum(), floated=False))
    return pd.DataFrame(rows)


print('=== (3) chi^2 per data set in the Delta region (150-350 MeV), alpha, beta from the NNLO fit, spin pols at NNLO ===')
P_F2 = np.array(fits['F2_d_all_thB']['P']); P_F2ff = np.array(fits['F2_d_all_thB_pi0ff']['P'])
tab0 = chi2_table(P_F2, T0)                 # paper (= delta_region.json chi2_F2)
tabf = chi2_table(P_F2ff, Tff)              # pi0ff=t, own fit
tabfP = chi2_table(P_F2, Tff)               # pi0ff=t at the paper's polarisabilities
ref = pd.DataFrame(DR['chi2_F2'])
print('    reproduction of delta_region.json chi2_F2: max |dchi2| = %.1e' % np.max(np.abs(tab0.chi2_n1.values - ref.chi2_n1.values)))
comb = tab0[['set', 'obs', 'N', 'Erange', 'floated']].copy()
comb['n1_f'] = tab0.chi2N_n1; comb['n1_t'] = tabf.chi2N_n1; comb['nfit_f'] = tab0.chi2N_nfit; comb['nfit_t'] = tabf.chi2N_nfit
comb['n_f'] = tab0.n; comb['n_t'] = tabf.n
comb['n1_t_paperP'] = tabfP.chi2N_n1; comb['nfit_t_paperP'] = tabfP.chi2N_nfit
comb['dchi2_n1'] = tabf.chi2_n1 - tab0.chi2_n1; comb['dchi2_nfit'] = tabf.chi2_nfit - tab0.chi2_nfit
pd.set_option('display.width', 250)
print(comb.round(4).to_string(index=False))
OUT['delta_chi2'] = comb.to_dict('records')

# ---------------------------------------------------------------- (5) observables
print('=== (5) relative change of dsigma/dOmega_lab and change of Sigma_3 (main-fit polarisabilities) ===')
ES = [100.0, 140.0, 150.0, 250.0, 300.0]
THS = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 170.0, 180.0]
pts = [('%s_%g_%g' % (o, E, th), E, th, o, 'LAB') for E in ES for th in THS for o in ('CS', 'S3', 'S2X', 'S2Z')]
A = bchpt.run(pts, 'NNLO', name='pi0ff_obs_f'); B = bchpt.run(pts, 'NNLO_pi0ff', name='pi0ff_obs_t')
va = A.predict_array(P_F2); vb = B.predict_array(P_F2)
obs = []
for i, p in enumerate(pts):
    o = p[3]
    obs.append(dict(id=p[0], obs=o, E=p[1], theta_lab=p[2], theta_cm=float(A.theta_cm[i]), val_f=float(va[i]), val_t=float(vb[i]),
                    rel=float(vb[i] / va[i] - 1) if o == 'CS' else None, diff=float(vb[i] - va[i])))
OUT['observables'] = obs
ot = pd.DataFrame(obs)
for o in ('CS', 'S3', 'S2X', 'S2Z'):
    piv = ot[ot.obs == o].pivot(index='E', columns='theta_lab', values='rel' if o == 'CS' else 'diff')
    print('--- %s: %s' % (o, 'dsigma(pi0ff)/dsigma - 1' if o == 'CS' else 'Sigma(pi0ff) - Sigma'))
    print(piv.map(lambda v: '%+.2e' % v).to_string())

json.dump(OUT, open(os.path.join(HERE, 'pi0ff_fits.json'), 'w'), indent=1, default=conv)

# ---------------------------------------------------------------- markdown
L = ['# pi0 gamma gamma form factor (pi0ff): fit impact', '',
     'Produced by `pi0ff_fits.py` (results in `pi0ff_fits.json`).  Units: alpha, beta 1e-4 fm^3; gamma 1e-4 fm^4.',
     'F2_d_all_thB = main configuration of the paper (main database, forward points as data, truncation covariance);',
     'F1_d_all_thB = Baldin fit; `_pi0ff` = amplitude and truncation covariance with the form factor; `_pi0ff_covF` = amplitude',
     'with the form factor but the truncation covariance of the pi0ff=f pipeline; `F?_d_all` = without truncation covariance.', '',
     '## (1) NNLO fits', '',
     '| fit | N | alpha_E1 | beta_M1 | alpha+beta | alpha-beta | chi^2/dof | d alpha | d beta | d chi^2 |',
     '|---|---|---|---|---|---|---|---|---|---|']
for kind in ['F2', 'F1']:
    base = fits['%s_d_all_thB' % kind]
    for sfx in ['_d_all_thB', '_d_all_thB_pi0ff', '_d_all_thB_pi0ff_covF', '_d_all', '_d_all_pi0ff']:
        r = fits[kind + sfx]
        b = base if 'thB' in sfx else fits['%s_d_all' % kind]
        L.append('| %s | %d | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.2f/%d | %+.4f | %+.4f | %+.3f |' % (
            kind + sfx, r['N'], r['P'][0], r['dP'][0], r['P'][1], r['dP'][1], r['apb'], r['dapb'], r['amb'], r['damb'],
            r['chi2'], r['dof'], r['P'][0] - b['P'][0], r['P'][1] - b['P'][1], r['chi2'] - b['chi2']))
L.append('')
for kind in ['F2', 'F1']:
    t = fits['%s_trunc_sigma_change' % kind]
    L.append('%s: truncation sigma_th recomputed with pi0ff=t: max relative change %.1e (mean %.1e).' % (kind, t['max_rel'], t['mean_rel']))
L.append('')
L += ['## (2) Delta-region spin-polarisability fit (S_A2_LEGS_gpi: A2 Sigma_2x, Sigma_2z, LEGS Sigma_3 < 2pi threshold; priors alpha+beta = 14.0(2), alpha-beta from the main fit, gamma_0 = -0.929(105), gamma_pi = 8.0(1.8))', '',
      '| fit | prior alpha-beta | gE1E1 | gM1M1 | gE1M2 | gM1E2 | gamma_0 | gamma_pi | alpha | beta | chi^2/dof |',
      '|---|---|---|---|---|---|---|---|---|---|---|']
for k, r in sp.items():
    L.append('| %s | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f(%.3f) | %.3f | %.3f | %.2f/%d |' % (
        k, r['prior_amb'][0], r['prior_amb'][1], r['P'][2], r['dP'][2], r['P'][3], r['dP'][3], r['P'][4], r['dP'][4], r['P'][5], r['dP'][5],
        r['g0'], r['dg0'], r['gpi'], r['dgpi'], r['P'][0], r['P'][1], r['chi2'], r['dof']))
r0_, r1_ = sp['S_A2_LEGS_gpi'], sp['S_A2_LEGS_gpi_pi0ff']
L += ['', 'Shifts (pi0ff - paper): dgE1E1 = %+.4f, dgM1M1 = %+.4f, dgE1M2 = %+.4f, dgM1E2 = %+.4f, dgamma_pi = %+.4f, dchi^2 = %+.3f' % tuple(
    [r1_['P'][i] - r0_['P'][i] for i in range(2, 6)] + [r1_['gpi'] - r0_['gpi'], r1_['chi2'] - r0_['chi2']]), '']
L += ['## (3) chi^2 per point in the Delta region (150 < E <= 350 MeV; alpha, beta from the respective NNLO fit, spin pols at NNLO)', '',
      '| data set | obs | N | E range | n=1 (f) | n=1 (pi0ff) | n fitted (f) | n fitted (pi0ff) | n (f) | n (pi0ff) | d chi^2 (n=1) | d chi^2 (n fit) |',
      '|---|---|---|---|---|---|---|---|---|---|---|---|']
for r in OUT['delta_chi2']:
    L.append('| %s | %s | %d | %s | %.3f | %.3f | %s | %s | %.3f | %.3f | %+.3f | %+.3f |' % (
        r['set'], r['obs'], r['N'], r['Erange'], r['n1_f'], r['n1_t'], ('%.3f' % r['nfit_f']) if r['floated'] else '--',
        ('%.3f' % r['nfit_t']) if r['floated'] else '--', r['n_f'], r['n_t'], r['dchi2_n1'], r['dchi2_nfit']))
L += ['', '## (5) Observables: dsigma(pi0ff)/dsigma - 1 (CS, lab) and Sigma(pi0ff) - Sigma, main-fit polarisabilities', '']
for o in ('CS', 'S3', 'S2X', 'S2Z'):
    piv = ot[ot.obs == o].pivot(index='E', columns='theta_lab', values='rel' if o == 'CS' else 'diff')
    L.append('%s (rows: E_lab [MeV]; columns: theta_lab [deg])' % o)
    L.append('')
    L.append('| E | ' + ' | '.join('%g' % c for c in piv.columns) + ' |')
    L.append('|---|' + '---|' * len(piv.columns))
    for E, row in piv.iterrows():
        L.append('| %g | ' % E + ' | '.join('%+.2e' % v for v in row.values) + ' |')
    L.append('')
open(os.path.join(HERE, 'pi0ff_fits.md'), 'w').write('\n'.join(L) + '\n')
print('written', os.path.join(HERE, 'pi0ff_fits.json'), os.path.join(HERE, 'pi0ff_fits.md'))
