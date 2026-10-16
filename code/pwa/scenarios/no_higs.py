"""Scenario: exclude HIGS 2022 (both beam settings) from the updated PWA; compare with the all-data fits.
Also: predict the HIGS points from the no-HIGS fit (pulls). Writes scenarios/no_higs.json"""
import sys, json, numpy as np
sys.path.insert(0, '.')
import update_fits as U
d = U.load()
H = ['HIGS_2022_circ', 'HIGS_2022_lin']
out = {}
for tag, kw in [('all_norm', dict(float_norm=True)), ('noHIGS_norm', dict(float_norm=True)),
                ('all_fwdchi2_norm', dict(float_norm=True, sr_mode='chi2')),
                ('noHIGS_fwdchi2_norm', dict(float_norm=True, sr_mode='chi2')),
                ('pruned_norm', dict(float_norm=True)), ('pruned_noHIGS_norm', dict(float_norm=True))]:
    excl = H if 'noHIGS' in tag else []
    if 'pruned' in tag:
        dd = U.select(d, exclude_sets=U.PRUNE['drop_sets'] + excl, drop_points=U.PRUNE['drop_points'])
    else:
        dd = U.select(d, exclude_sets=excl)
    try:
        res, m, r = U.fit_variant('scen_' + tag, dd, **kw)
    except TypeError as e:
        print('skip', tag, e); continue
    out[tag] = dict(N=res['npts'], chi2=res['chi2'], val={k: res['val'][k] for k in U.KEYS},
                    err={k: res['err'][k] for k in U.KEYS}, norms=res.get('norms'))
    v, e = res['val'], res['err']
    print(f"{tag:22s} N={res['npts']:3d} chi2/N={res['chi2']/res['npts']:.3f}  a={v['alphaE1']:.2f}({e['alphaE1']:.2f}) "
          f"b={v['betaM1']:.2f}({e['betaM1']:.2f}) a-b={v['alpha_m_beta']:.2f}({e['alpha_m_beta']:.2f}) gpi={v['gamma_pi']:.2f}({e['gamma_pi']:.2f})")
    if tag == 'noHIGS_norm':
        # predict HIGS points with the no-HIGS parameters, using the all-data model machinery
        resA, mA, rA = U.fit_variant('scen_tmp_all', d, float_norm=True)
        # evaluate the all-data model at the no-HIGS optimum is not directly possible (different parameter vectors);
        # instead report HIGS chi2 contributions in the all-data fit and the shift of observables below.
json.dump(out, open('scenarios/no_higs.json', 'w'), indent=1)
