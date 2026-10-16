"""run_radcorr_fit.py -- impact of the O(alpha) QED corrections on the polarisability fits.

(1) reproduces the main fit (F2_d_all_thB: all data + GHP forward points, truncation covariance) and the
    'without truncation covariance' (F2_d_all) and Baldin (F1_d_all_thB, F1_d_all) fits of work/fits/results/fits_L.json
    with the same machinery;
(2) re-runs them with the theory prediction of every cross-section row (obs CS, CS_PAR, CS_PERP; not the GHP forward
    points, not Sigma_3 / Sigma_2z) multiplied by 1 + delta_i(omega_i, theta_lab,i; Delta E), delta from radcorr.py;
(3) reports the shifts of alpha_E1, beta_M1, alpha+beta, alpha-beta and chi^2 for several Delta E choices.

Outputs: work/radcorr/results/fit_shifts.{json,md}.  Nothing in work/fits/results is touched.
"""
import os, sys, json, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(WORK, 'fits', 'lib'))
sys.path.insert(0, os.path.join(WORK, 'fits', 'scripts'))
sys.path.insert(0, HERE)
import radcorr as rc
import s03_fits_L as s03                       # module-level: loads db and the cached NNLO driver outputs
from common import build_fit, FWD, PREF_NNLO, C_BALDIN, BALDIN, gamma0, gammapi

OUT = os.path.join(HERE, 'results')
os.makedirs(OUT, exist_ok=True)
REF = json.load(open(os.path.join(WORK, 'fits', 'results', 'fits_L.json')))
CS_OBS = {'CS': None, 'CS_PAR': 'par', 'CS_PERP': 'perp'}

# --- soft-photon cut per data set ('experimental' choice; see notes/radiative_corrections.md) -------------------
# A2 2022: missing-mass window 920-955 MeV (Mornacchi thesis Sect. 7.2) -> an undetected photon of energy E raises
#          M_miss by ~E, so E < ~17 MeV survives the cut; we take 15 MeV.  All others: window not quantified
#          -> default 3% of omega.
DE_EXPT = {'A2_2022_CS': 15.0}
DE_DEFAULT_FRAC = 0.03

_cache = {}


def delta_row(E, thl, dE, pol):
    key = (round(float(E), 6), round(float(thl), 6), round(float(dE), 6), pol)
    if key not in _cache:
        _cache[key] = rc.delta(E, thl, dE, pol=pol)
    return _cache[key]


def delta_vector(s, spec):
    """spec: float (uniform Delta E in MeV), 'frac:x' (Delta E = x*omega), or 'expt' (per data set)."""
    d = np.zeros(len(s))
    for i, r in enumerate(s.itertuples()):
        if r.dataset in FWD or r.obs not in CS_OBS:
            continue
        pol = CS_OBS[r.obs]
        if spec in ('expt', 'expt_unpol'):
            dE = DE_EXPT.get(r.dataset, DE_DEFAULT_FRAC * r.E)
            if spec == 'expt_unpol':
                pol = None                      # unpolarised delta also for the HIGS dsigma_par/perp rows
        elif isinstance(spec, str) and spec.startswith('frac:'):
            dE = float(spec[5:]) * r.E
        else:
            dE = float(spec)
        d[i] = delta_row(r.E, r.theta_lab, dE, pol)
    return d


def do_fit_rc(kind, thcov, delta=None):
    """replica of s03.do_fit(name, ALL, kind, thcov=thcov) with the theory multiplied by (1 + delta)."""
    df, T = s03.df, s03.T
    s = s03.rows_for(df, s03.ALL, kind, 150.0, (), s03.FWD_EMAX, None)
    free = ['alphaE1', 'betaM1']
    cons = [(C_BALDIN, BALDIN[0], BALDIN[1], 'Baldin')] if kind == 'F1' else []
    C_th = s03.thcov_for(s, df, 'db')[0] if thcov else None
    idx = T.subset(s['pid'].values)
    fac = np.ones(len(s)) if delta is None else 1.0 + delta

    def theory_fn(P, ex):
        return T.predict_array(P, idx) * fac
    F = build_fit(s, T, PREF_NNLO, free, constraints=cons, fwd_mode='nuisance', thcov=C_th, theory_fn=theory_fn)
    R = F.run(decompose=False)
    out = dict(kind=kind, thcov=thcov, N=len(s), P=R.P[:2].tolist(), dP=R.dP[:2].tolist(), chi2=float(R.chi2), dof=int(R.dof),
               apb=float(R.P[0] + R.P[1]), amb=float(R.P[0] - R.P[1]),
               dapb=float(np.sqrt(R.covP[0, 0] + R.covP[1, 1] + 2 * R.covP[0, 1])),
               damb=float(np.sqrt(R.covP[0, 0] + R.covP[1, 1] - 2 * R.covP[0, 1])),
               norms={k: float(v['n']) for k, v in R.perset.items() if v['dn'] > 0})
    return s, out


FITS = [('F2', True, 'F2_d_all_thB', 'Main (forward points + truncation cov.)'),
        ('F2', False, 'F2_d_all', 'Main without truncation covariance'),
        ('F1', True, 'F1_d_all_thB', 'Baldin (alpha+beta = 14.0(2)) + truncation cov.'),
        ('F1', False, 'F1_d_all', 'Baldin without truncation covariance')]
SPECS = [('none', None), ('1 MeV', 1.0), ('3 MeV', 3.0), ('10 MeV', 10.0), ('3% of omega', 'frac:0.03'),
         ('1% of omega', 'frac:0.01'), ('10% of omega', 'frac:0.10'), ('expt (A2: 15 MeV; else 3%)', 'expt'), ('expt, unpol. delta for HIGS par/perp', 'expt_unpol')]

if __name__ == '__main__':
    results = {}
    lines = ['# Impact of the O(alpha) QED corrections on the fits (run_radcorr_fit.py)', '',
             'Theory of every cross-section row multiplied by 1 + delta_i (point-like proton, virtual + soft, lab cut Delta E);',
             'Sigma_3, Sigma_2z and the GHP forward points unchanged.  Units 1e-4 fm^3.  Reference values from fits_L.json.', '']
    t0 = time.time()
    for kind, th, key, label in FITS:
        ref = REF[key]
        lines += ['## %s  [`%s`]' % (label, key), '',
                  'reference (fits_L.json): alpha = %.4f(%.4f), beta = %.4f(%.4f), a+b = %.4f(%.4f), a-b = %.4f(%.4f), chi2/dof = %.3f/%d' % (
                      ref['P'][0], ref['dP'][0], ref['P'][1], ref['dP'][1], ref['apb'], ref['dapb'], ref['amb'], ref['damb'], ref['chi2'], ref['dof']), '',
                  '| Delta E | <delta> CS rows [%] | max delta [%] | alpha | beta | a+b | a-b | chi2 | d alpha | d beta | d(a+b) | d(a-b) | d chi2 |',
                  '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
        base = None
        results[key] = dict(label=label, reference=dict(P=ref['P'][:2], dP=ref['dP'][:2], chi2=ref['chi2'], dof=ref['dof'],
                                                        apb=ref['apb'], amb=ref['amb']), variants={})
        for sname, spec in SPECS:
            s, _ = do_fit_rc(kind, th, None) if spec is None and base is None else (None, None)
            if spec is None:
                s, out = do_fit_rc(kind, th, None)
                base = out
                dvec = np.zeros(out['N'])
                # exact reproduction check
                ok = (abs(out['P'][0] - ref['P'][0]) < 1e-6 and abs(out['P'][1] - ref['P'][1]) < 1e-6
                      and abs(out['chi2'] - ref['chi2']) < 1e-5 and out['dof'] == ref['dof'])
                out['reproduces_reference'] = bool(ok)
                print('%-14s %-28s a=%.5f b=%.5f chi2=%.4f/%d  reproduces fits_L.json: %s' % (key, sname, out['P'][0], out['P'][1], out['chi2'], out['dof'], ok))
            else:
                s = s03.rows_for(s03.df, s03.ALL, kind, 150.0, (), s03.FWD_EMAX, None)
                dvec = delta_vector(s, spec)
                _, out = do_fit_rc(kind, th, dvec)
                print('%-14s %-28s a=%.5f b=%.5f chi2=%.4f  da=%+.5f db=%+.5f  <delta>=%.2e' % (
                    key, sname, out['P'][0], out['P'][1], out['chi2'], out['P'][0] - base['P'][0], out['P'][1] - base['P'][1], dvec[dvec != 0].mean()))
            iscs = dvec != 0
            out['delta_mean'] = float(dvec[iscs].mean()) if iscs.any() else 0.0
            out['delta_max'] = float(np.abs(dvec).max())
            out['delta_by_set'] = {ds: float(dvec[(s['dataset'] == ds).values].mean()) for ds in s['dataset'].unique()
                                   if ((s['dataset'] == ds).values & iscs).any()}
            out['shift'] = dict(alpha=out['P'][0] - base['P'][0], beta=out['P'][1] - base['P'][1], apb=out['apb'] - base['apb'],
                                amb=out['amb'] - base['amb'], chi2=out['chi2'] - base['chi2'])
            results[key]['variants'][sname] = out
            sh = out['shift']
            lines.append('| %s | %+.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %+.4f | %+.4f | %+.4f | %+.4f | %+.3f |' % (
                sname, 100 * out['delta_mean'], 100 * out['delta_max'], out['P'][0], out['P'][1], out['apb'], out['amb'], out['chi2'],
                sh['alpha'], sh['beta'], sh['apb'], sh['amb'], sh['chi2']))
        lines.append('')
    # mean delta per data set for the 'expt' and '3 MeV' choices (main fit rows)
    s = s03.rows_for(s03.df, s03.ALL, 'F2', 150.0, (), s03.FWD_EMAX, None)
    lines += ['## Mean delta [%] per data set (cross-section rows)', '', '| data set | N_CS | <E> | Delta E = 1 MeV | 3 MeV | 10 MeV | 3% | expt |', '|---|---|---|---|---|---|---|---|']
    vecs = {sp: delta_vector(s, sp) for sp in [1.0, 3.0, 10.0, 'frac:0.03', 'expt']}
    for ds in s['dataset'].unique():
        m = (s['dataset'] == ds).values & (vecs[3.0] != 0)
        if not m.any():
            continue
        lines.append('| %s | %d | %.1f | ' % (ds, m.sum(), s[m]['E'].mean()) + ' | '.join('%+.4f' % (100 * vecs[sp][m].mean()) for sp in [1.0, 3.0, 10.0, 'frac:0.03', 'expt']) + ' |')
    lines.append('')
    lines.append('run time %.0f s' % (time.time() - t0))
    json.dump(results, open(os.path.join(OUT, 'fit_shifts.json'), 'w'), indent=1)
    open(os.path.join(OUT, 'fit_shifts.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
