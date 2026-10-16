"""Scenario: exclude HIGS 2022 -- forward points as chi2 data (correlated nuisances), floating norms."""
import sys, json
sys.path.insert(0, '.')
import update_fits as U
d = U.load(); H = ['HIGS_2022_circ', 'HIGS_2022_lin']
PW = (10., 1000., 1000., 10., 1., 1.)
out = json.load(open('scenarios/no_higs.json'))
for k in [k for k in out if 'fwdchi2' in k]: out.pop(k)
for tag, dd in [('all_fwdchi2_norm', U.select(d)), ('noHIGS_fwdchi2_norm', U.select(d, exclude_sets=H)),
                ('pruned_fwdchi2_norm', U.select(d, exclude_sets=U.PRUNE['drop_sets'], drop_points=U.PRUNE['drop_points'])),
                ('pruned_noHIGS_fwdchi2_norm', U.select(d, exclude_sets=U.PRUNE['drop_sets'] + H, drop_points=U.PRUNE['drop_points']))]:
    res, m, r = U.fit_variant('scen_' + tag, dd, sr_mode='chi2', fwd_err='nuisance', float_norm=True, prior_w=PW)
    v, e = res['val'], res['err']
    out[tag] = dict(N=res['npts'], chi2=res['chi2'], val={k: v[k] for k in U.KEYS}, err={k: e[k] for k in U.KEYS}, norms=res.get('norms'))
    print(f"{tag:28s} N={res['npts']:3d} chi2/N={res['chi2']/res['npts']:.3f}  a={v['alphaE1']:.2f}({e['alphaE1']:.2f}) "
          f"b={v['betaM1']:.2f}({e['betaM1']:.2f}) a-b={v['alpha_m_beta']:.2f}({e['alpha_m_beta']:.2f}) "
          f"a+b={v['alpha_p_beta']:.2f} gpi={v['gamma_pi']:.2f}({e['gamma_pi']:.2f})")
json.dump(out, open('scenarios/no_higs.json', 'w'), indent=1)
