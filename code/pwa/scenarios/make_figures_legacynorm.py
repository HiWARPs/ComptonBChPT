"""Regenerate the PWA figures with the legacy-set normalisations floated (new PWA baseline).
Old figures (fixed legacy normalisations) archived in figs_v1fixednorm/."""
import sys, runpy
sys.path.insert(0, '.')
import update_fits as UF
_orig = UF.load
LEG = {'Oxley': 0.10, 'Pugh': 0.10, 'Hyman': 0.10, 'Bernardini': 0.10, 'Goldansky': 0.08}
def load_leg(per_beamtime=False):
    d = _orig(per_beamtime=per_beamtime)
    for s, e in LEG.items():
        m = d.dataset == s
        d.loc[m, 'norm_group'] = s; d.loc[m, 'norm_err'] = e
    return d
UF.load = load_leg
runpy.run_path('make_figures.py', run_name='__main__')
