"""Regenerate the PWA figures with the legacy-set normalisations floated (new PWA baseline).
Old figures (fixed legacy normalisations) archived in figs_legacynorm/."""
import sys, runpy
sys.path.insert(0, '.')
import update_fits as UF
_orig = UF.load
LEG = UF.LEGACY_NORM   # Baranov et al. 2001 delta (was 10%, Goldansky 8%)
def load_leg(per_beamtime=False):
    d = _orig(per_beamtime=per_beamtime)
    for s, e in LEG.items():
        m = d.dataset == s
        d.loc[m, 'norm_group'] = s; d.loc[m, 'norm_err'] = e
    return d
UF.load = load_leg
runpy.run_path('make_figures.py', run_name='__main__')
