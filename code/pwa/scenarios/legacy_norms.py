"""Referee request: PWA MAIN with floating normalisations also for the legacy sets (Oxley, Pugh, Hyman,
Bernardini: 10%; Goldansky: 8%; Zieger: as BChPT, not floated), as in the BChPT fits. Also influence of
dropping Oxley / Goldansky in this configuration. Writes scenarios/legacy_norms.json"""
import sys, json
sys.path.insert(0, '.')
import update_fits as UF
import main_fits as MF
d = UF.load()
leg = {'Oxley': 0.10, 'Pugh': 0.10, 'Hyman': 0.10, 'Bernardini': 0.10, 'Goldansky': 0.08}
dl = d.copy()
for s, e in leg.items():
    m = dl.dataset == s
    dl.loc[m, 'norm_group'] = s; dl.loc[m, 'norm_err'] = e
out = {}
for tag, dd in [('MAIN', UF.select(d)), ('MAIN_legacynorm', UF.select(dl)),
                ('legacynorm_noOxley', UF.select(dl, exclude_sets=['Oxley'])),
                ('legacynorm_noGoldansky', UF.select(dl, exclude_sets=['Goldansky'])),
                ('legacynorm_noLegacy', UF.select(dl, exclude_sets=list(leg))),
                ('legacynorm_pruned', UF.select(dl, exclude_sets=UF.PRUNE['drop_sets'], drop_points=UF.PRUNE['drop_points']))]:
    res, m, r = UF.fit_variant('scen_' + tag, dd, **MF.MAIN)
    v, e = res['val'], res['err']
    out[tag] = dict(N=res['npts'], chi2=res['chi2'], val={k: v[k] for k in UF.KEYS}, err={k: e[k] for k in UF.KEYS},
                    norms=res.get('norms'), per_dataset=res.get('per_dataset'))
    ox = (res.get('per_dataset') or {}).get('Oxley', {})
    print(f"{tag:24s} N={res['npts']:3d} chi2/N={res['chi2']/res['npts']:.3f} a-b={v['alpha_m_beta']:.2f}({e['alpha_m_beta']:.2f}) "
          f"b={v['betaM1']:.2f}({e['betaM1']:.2f}) gpi={v['gamma_pi']:.2f}  Oxley chi2={ox.get('chi2','-')}  "
          f"normOx={(res.get('norms') or {}).get('Oxley',{}).get('N','-')}")
json.dump(out, open('scenarios/legacy_norms.json', 'w'), indent=1, default=float)
