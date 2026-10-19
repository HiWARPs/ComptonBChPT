"""Scan of the main BChPT fit (also above w = 0.1: weight_scan.py high 0.15,0.2,... 1,2) (scenarios/combined_fit.py) in the weight w of the region above 150 MeV (Delta-region angular
data and forward SRs) up to w = 0.1, with the spin-polarisability priors at their re-examined widths (1 sigma) and doubled
(2 sigma) -- author question, 2026-10-19: "increase the weight of the Delta region to up to w=0.1 and the variation of the
spin polarizabilities in BChPT fit to 2 sigma; what changes?".  For each fit: alpha, beta, R_EM, the spin polarisabilities,
and the unweighted chi^2 of the three parts at the fitted parameters (267 points below threshold with the truncation
covariance; 289 Delta-region angular points with the band; 84 forward-SR amplitude values above threshold with the band).
Writes results/weight_scan.json."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.interpolate import CubicSpline
import combined_fit as CF
from common import build_fit, save_json

WS = [0.01, 0.02, 0.03, 0.05, 0.1]
OUT = {}


def parts(P, g, legs=False):
    """unweighted chi^2 of the three parts at fixed P, g_E (normalisations and forward nuisances re-optimised); legs: include
    the LEGS sets in the Delta part"""
    out = {}
    for part in ('L', 'D'):
        s, C = (CF.sL, CF.C_L) if part == 'L' else (lambda r: (r[0], np.diag(r[1] ** 2)))(CF.rows_delta(360.0, legs))
        idx = CF.T0.subset(s.pid.values)
        tf = lambda P_, ex, idx=idx: CubicSpline(CF.grid, np.array([CF.G[gg].predict_array(P_, idx) for gg in CF.GE_GRID]), axis=0)(g)
        out[part] = build_fit(s, CF.T0, P, [], thcov=C, theory_fn=tf).run(decompose=False).chi2
    yb, Cb, _, tb = CF.fsr_block(360.0, 10.0, 1.0)
    r = yb - tb(P, g)
    out['FSR'] = float(r @ np.linalg.solve(Cb, r))
    return out


if __name__ == '__main__':
    # optional: weight_scan.py <tag> <w1,w2,...> <k1,k2,...>  (prior widths k x re-examined); default: the scan up to 0.1
    tag = sys.argv[1] if len(sys.argv) > 1 else ''
    WS = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else WS
    KS = [float(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1, 2]
    for k, sig in ((k, k * CF.SIG_SPIN) for k in KS):
        for w in WS:
            name = 'w%g_spin%gsig' % (w, k)
            r = CF.fit(name, w=w, sig=sig)
            r['parts'] = parts(np.array(r['P']), r['gE'][0])
            r['prior_sig'] = sig
            OUT[name] = r
            print('   parts (unweighted): below thr %.1f/267, Delta %.1f/289, FSR %.1f/84' % (r['parts']['L'], r['parts']['D'], r['parts']['FSR']))
    save_json(OUT, 'weight_scan%s.json' % (('_' + tag) if tag else ''))
