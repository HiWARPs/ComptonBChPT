"""Unweighted chi^2 parts of the prediction, the main BChPT fit and its with-LEGS variant (paper Sect. fits: overall chi^2
of the 484 values, the LEGS paragraph): below threshold (267 values, truncation covariance), Delta-region angular data
without and with the LEGS sets (band), forward-SR amplitudes above threshold (84, band), and the LEGS sets per set at the
fit parameters (normalisations and forward nuisances re-optimised; n = 1 for LEGS evaluated at a fit that excludes them is
not used here).  Reads results/combined_fit.json; writes results/combined_fit_parts.json (producer saved 2026-10-28; the
2026-10-19/20 file came from an unsaved script with the same definitions as weight_scan.parts)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json
import numpy as np
from scipy.interpolate import CubicSpline
import combined_fit as CF
from weight_scan import parts
from common import build_fit, save_json, PREF_NNLO, RES

J = json.load(open(os.path.join(RES, 'combined_fit.json')))


def legs_sets(P, g):
    sD, sth = CF.rows_delta(360.0, True)
    out = {}
    for ds in CF.LEGS:
        m = (sD.dataset == ds).values
        s = sD[m]; idx = CF.T0.subset(s.pid.values)
        tf = lambda P_, ex, idx=idx: CubicSpline(CF.grid, np.array([CF.G[gg].predict_array(P_, idx) for gg in CF.GE_GRID]), axis=0)(g)
        out[ds] = [round(build_fit(s, CF.T0, P, [], thcov=np.diag(sth[m] ** 2), theory_fn=tf).run(decompose=False).chi2, 1), int(m.sum())]
    return out


OUT = {}
for k, P, g in (('prediction', PREF_NNLO, CF.GE0), ('combined', J['combined']['P'], J['combined']['gE'][0]),
                ('combined_withLEGS', J['combined_withLEGS']['P'], J['combined_withLEGS']['gE'][0])):
    P = np.array(P, float)
    r = parts(P, g); r['D_withLEGS'] = parts(P, g, legs=True)['D']; r['LEGS'] = legs_sets(P, g)
    r['overall_noLEGS'] = r['L'] + r['D'] + r['FSR']
    OUT[k] = r
    print(k, {a: (round(b, 1) if isinstance(b, float) else b) for a, b in r.items()})
save_json(OUT, 'combined_fit_parts.json')
