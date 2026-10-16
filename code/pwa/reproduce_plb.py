"""Reproduction of the fits of Krupina, Lensky, Pascalutsa, PLB 782 (2018) 34 (Table II and III).
Database: work/data/db_klp2018.csv (= author notebook `datafull`, 138 points).
Kinematics as in the original notebook: all data of an energy bin evaluated at the bin energy.
Sum rules: Re f, Re g of GHP as used in the notebook (plb_inputs.PLBForward), alpha+beta=14.0, gamma0=-0.929.
Fitting methods:
  MAP  : joint minimisation of chi^2 + sum_bins sum_i (fR_i/w_i)^2, w = (10,10,1,1) [KLP18 footnote 2 widths]
  MC   : KLP18 procedure: Monte-Carlo sampling of all parameters from the footnote-2 Gaussians, best chi^2
         retained, then minimisation over the global parameters with residuals fixed (several seeds)
  errors 'stage2' = global-parameter errors with residuals fixed (KLP18 convention); 'full' = marginalised
         over the residual functions with the Gaussian priors.
Writes results/plb_reproduction.json and prints tables."""
import json
import numpy as np
import pandas as pd
import pwa_model as PM
import plb_inputs as PI

PLB = {  # Table II / III of KLP18
    'Fit 0': dict(alphaE1=(12.2, 0.3), betaM1=(1.8, 0.3), gE1E1=(-1.6, 2.6), gM1M1=(1.8, 1.1), gE1M2=(-1.3, 3.7), gM1E2=(2.0, 0.7), chi2pt=1.35),
    'Fit 1': dict(alphaE1=(12.2, 0.3), betaM1=(1.8, 0.3), gE1E1=(-3.1, 0.7), gM1M1=(1.6, 0.3), gE1M2=(0.0, 0), gM1E2=(2.5, 0.7), chi2pt=1.35,
                  alpha_m_beta=(10.5, 0.4), gamma_pi=(7.2, 1.0)),
    'Fit 1_3s': dict(alphaE1=(11.8, 0.3), betaM1=(2.2, 0.3), gE1E1=(-2.7, 0.6), gM1M1=(1.5, 0.3), gE1M2=(0.0, 0), gM1E2=(2.2, 0.7), chi2pt=0.97),
    "Fit 1'": dict(alphaE1=(10.6, 0.3), betaM1=(3.4, 0.3), gE1E1=(-1.0, 0.8), gM1M1=(1.0, 0.3), gE1M2=(0.0, 0), gM1E2=(1.0, 0.7), chi2pt=0.99,
                   alpha_m_beta=(7.2, 0.6), gamma_pi=(3.0, 1.1)),
    "Fit 1''": dict(alphaE1=(10.2, 0.4), betaM1=(3.8, 0.4), gE1E1=(-1.2, 0.8), gM1M1=(0.6, 0.3), gE1M2=(0.0, 0), gM1E2=(1.6, 0.8), chi2pt=0.62,
                    alpha_m_beta=(6.4, 0.6), gamma_pi=(3.5, 1.2)),
    'Fit 2': dict(alphaE1=(11.7, 0.3), betaM1=(2.3, 0.3), gE1E1=(-2.6, 0.6), gM1M1=(1.1, 0.3), gE1M2=(0.0, 0), gM1E2=(2.4, 0.7), chi2pt=1.35),
    "Fit 2''": dict(alphaE1=(10.8, 0.4), betaM1=(3.2, 0.4), gE1E1=(-1.9, 0.8), gM1M1=(0.7, 0.3), gE1M2=(0.0, 0), gM1E2=(2.2, 0.8), chi2pt=0.69),
}
KEYS = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2', 'alpha_m_beta', 'gamma_pi']


def select(d, which):
    keep = np.ones(len(d), bool)
    err = d.err.values.copy()
    if which in ("1'",):
        drop = [('OlmosDeLeon', 89, 133), ('OlmosDeLeon', 109, 133), ('Oxley', 60, 120), ('Oxley', 60, 150)]
        for ex, E, th in drop:
            keep &= ~((d.dataset == ex) & (abs(d.E - E) < 0.6) & (abs(d.theta_lab - th) < 0.6)).values
    if which in ("1''", "2''"):
        keep &= ~(d.dataset == 'Oxley').values
        keep &= ~(d.dataset == 'Bernardini').values
        keep &= ~((d.dataset == 'Baranov') & (abs(d.theta_lab - 150) < 0.6)).values
        keep &= ~((d.dataset == 'OlmosDeLeon') & (abs(d.E - 109) < 0.6) & (abs(d.theta_lab - 133) < 0.6)).values
        odl = (d.dataset == 'OlmosDeLeon').values
        err[odl] = np.hypot(err[odl], 0.05 * d.value.values[odl])
    out = d[keep].copy()
    out['err'] = err[keep]
    return out


def make_model(d, l2='bchpt', gEM_free=False):
    return PM.Model(d, PI.BINS, l2=l2, energy_mode='bin', fix_gEM=None if gEM_free else 0.0,
                    sr_mode='constraint', fwd=PI.PLBForward())


def mc_fit(m, seed, ndraw=400000, batch=20000, means=None, widths=None):
    """KLP18 Monte-Carlo sweep + stage-2 minimisation."""
    rng = np.random.default_rng(seed)
    npm = m.ng + m.nr * m.nb
    mu = np.zeros(npm); sg = np.zeros(npm)
    gm = {'D': (9.5, 2.0), 'gEE': (-3.3, 2.0), 'gME': (1.1, 2.0), 'gEM': (0.2, 2.0)}
    for j, n in enumerate(m.gnames):
        mu[j], sg[j] = gm[n]
    w = {'fR1': 10., 'fR4': 10., 'fR5': 1., 'fR6': 1.}
    for b in range(m.nb):
        for r, nm in enumerate(m.rnames):
            sg[m.ng + b * m.nr + r] = w[nm]
    best = (np.inf, None)
    idx = np.where(m.mask)[0]
    u, V, lf, y, e = m.u[idx], m.V[idx], m.lf[idx] * 1e7, m.y[idx], m.e[idx]
    for _ in range(ndraw // batch):
        P = mu + sg * rng.standard_normal((batch, npm))
        Phi = u[:, :, None] + np.einsum('ikj,nj->ikn', V, P)
        cs = np.einsum('k,ikn->in', PM.WCS, Phi**2) * lf[:, None]
        chi2 = np.sum(((cs - y[:, None]) / e[:, None])**2, axis=0)
        k = np.argmin(chi2)
        if chi2[k] < best[0]:
            best = (chi2[k], P[k].copy())
    p = np.zeros(m.npar); p[:npm] = best[1]
    r = m.fit_globals_only(p)
    return best[0], r


def run(which, l2='bchpt', gEM_free=False, outlier_from=None, mc_seeds=0):
    d0 = PI.load_klp2018()
    d = select(d0, which.replace('Fit ', '').replace('_3s', ''))
    if which.endswith('_3s'):
        d = d[~outlier_from].copy()
    m = make_model(d, l2=l2, gEM_free=gEM_free)
    rmap = m.fit(with_prior=True)
    r2 = m.fit_globals_only(rmap.p)                    # stage-2 errors around the MAP point
    res = dict(npts=len(d), chi2=rmap.chi2, chi2pt=rmap.chi2 / len(d),
               val={k: rmap.der[k] for k in KEYS}, err_full={k: rmap.der_err[k] for k in KEYS},
               err_stage2={k: r2.der_err[k] for k in KEYS})
    if mc_seeds:
        mcs = []
        for s in range(mc_seeds):
            c2mc, rmc = mc_fit(m, s)
            mcs.append(dict(seed=s, chi2_mc=float(c2mc), chi2=float(rmc.chi2), val={k: rmc.der[k] for k in KEYS},
                            err={k: rmc.der_err[k] for k in KEYS}))
        res['mc'] = mcs
    return res, m, rmap


if __name__ == '__main__':
    out = {}
    r1, m1, f1 = run('Fit 1', mc_seeds=6)
    np.save('results/plb_fit1_p.npy', f1.p); np.save('results/plb_fit1_cov.npy', f1.cov)
    out['Fit 1'] = r1
    pulls, _ = m1.pulls(f1.p)
    d1 = m1.d.copy(); d1['pull'] = pulls
    d1.to_csv('results/plb_fit1_pulls.csv', index=False)
    outl = np.abs(pulls) > 3
    print('Fit 1 3-sigma outliers:')
    print(d1[outl][['dataset', 'E', 'bin', 'theta_lab', 'value', 'err', 'pull']].to_string())
    # map the outliers back onto the full database rows
    d0 = PI.load_klp2018()
    key0 = list(zip(d0.dataset, d0.E, d0.theta_lab))
    okeys = set(zip(d1[outl].dataset, d1[outl].E, d1[outl].theta_lab))
    omask = np.array([k in okeys for k in key0])
    out['Fit 0'] = run('Fit 0', gEM_free=True, mc_seeds=3)[0]
    out['Fit 1_3s'] = run('Fit 1_3s', outlier_from=omask, mc_seeds=3)[0]
    out["Fit 1'"] = run("Fit 1'", mc_seeds=3)[0]
    rpp, mpp, fpp = run("Fit 1''", mc_seeds=6)
    out["Fit 1''"] = rpp
    np.save('results/plb_fit1pp_p.npy', fpp.p); np.save('results/plb_fit1pp_cov.npy', fpp.cov)
    # compare angular distributions with the published curves (Fig. 5)
    import digitize_plb as Dg
    cv = Dg.curves()
    ths = np.arange(20., 161., 10.)
    print('\nmax |rel. dev.| of reproduced curves from published KLP18 Fig.5 curves (theta_lab 20-160):')
    for b, E in enumerate(PI.BINS):
        key = 149 if E == 150 else E
        row = []
        for nm, (mm, ff) in (('fit1', (m1, f1)), ('fit1pp', (mpp, fpp))):
            t, v = cv[(key, nm)]
            mine = mm.curve(float(E), ths, b, ff.p)
            row.append(f'{nm}: {np.max(abs(mine/np.interp(ths, t, v)-1)):.3f}')
        print(f'   {E:4d} MeV  ' + '  '.join(row))
    # 3-sigma fit with exactly the 4 outliers named in KLP18 (Fig. 2 caption)
    four = [('OlmosDeLeon', 89, 155), ('OlmosDeLeon', 109, 133), ('Oxley', 60, 150), ('Goldansky', 55, 150)]
    m4 = np.zeros(len(d0), bool)
    for ex, E, th in four:
        m4 |= ((d0.dataset == ex) & (abs(d0.E - E) < 0.6) & (abs(d0.theta_lab - th) < 0.6)).values
    out['Fit 1_3s(PLB 4 pts)'] = run('Fit 1_3s', outlier_from=m4, mc_seeds=3)[0]
    PLB['Fit 1_3s(PLB 4 pts)'] = PLB['Fit 1_3s']
    out['Fit 2'] = run('Fit 2', l2='zero', mc_seeds=3)[0]
    out["Fit 2''"] = run("Fit 2''", l2='zero', mc_seeds=3)[0]
    out['outliers_fit1'] = d1[outl][['dataset', 'E', 'theta_lab', 'pull']].to_dict('records')
    json.dump(out, open('results/plb_reproduction.json', 'w'), indent=1, default=float)
    # print comparison
    for name, r in out.items():
        if not name.startswith('Fit'):
            continue
        print(f'\n{name}: N={r["npts"]} chi2/pt = {r["chi2pt"]:.3f}  chi2/138 = {r["chi2"]/138:.3f}  (PLB {PLB[name]["chi2pt"]})')
        for k in KEYS:
            ref = PLB[name].get(k)
            s = f'   {k:13s} MAP {r["val"][k]:6.2f} +- {r["err_stage2"][k]:.2f} (stage2) +- {r["err_full"][k]:.2f} (full)'
            if ref:
                s += f'   PLB {ref[0]:5.1f} +- {ref[1]:.1f}'
            if 'mc' in r:
                v = np.array([x['val'][k] for x in r['mc']])
                s += f'   | MC+stage2 seeds: {v.mean():6.2f} (spread {v.std():.2f}, range {v.min():.2f}..{v.max():.2f})'
            print(s)
        if 'mc' in r:
            print('   MC chi2/pt per seed:', [round(x['chi2'] / r['npts'], 3) for x in r['mc']])
