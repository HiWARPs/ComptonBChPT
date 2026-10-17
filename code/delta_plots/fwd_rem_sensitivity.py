"""fwd_rem_sensitivity.py -- how strongly the forward observables depend on the E2/M1 ratio in the resonance region:
complete NNLO prediction (own static polarisabilities) of dsigma/dOmega_lab(0) and Sigma_2z(0) for 150-400 MeV on the
coupling grid of config.BAND_GRID (g_M = 2.97 +- 0.10, g_E such that R_EM = -2.28 -+ 0.5 %) and for g_E = -1.0 +- 0.5 at
g_M = 2.97.  Output: results/gE_scan/fwd_rem_sensitivity.csv/.md"""
import os, sys
os.environ.setdefault('DP_EXTRA', 'none'); os.environ.setdefault('DP_MAIN', 'pred')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import config as C, runner
E = np.arange(150.0, 400.01, 10.0)
pts = [('fw_%s_%.0f' % (o, e), e, 0.0, o, 'LAB') for o in ('CS', 'S2Z') for e in E]
def pred(gM, gE):
    T = runner.run(pts, gM=gM, gE=gE, name='fwd_gM%.2f_gE%.3f' % (gM, gE), verbose=False)
    v = T.predict_array(T.Pref)
    return v[:len(E)], v[len(E):]
cs0, s0 = pred(C.GM0, C.GE0)
rows = []
for gM, gE in list(C.BAND_GRID) + [(C.GM0, -1.5), (C.GM0, -0.5)]:
    cs, s2 = pred(gM, gE)
    R = 100 * C.C_REM * gE / gM
    for e, a, b, a0, b0 in zip(E, cs, s2, cs0, s0):
        rows.append(dict(gM=gM, gE=round(gE, 3), REM=round(R, 2), E=e, cs=a, s2z=b, dcs_rel=a / a0 - 1, ds2z=b - b0))
t = pd.DataFrame(rows)
OUT = os.path.join(C.RESULTS, 'gE_scan'); os.makedirs(OUT, exist_ok=True)
t.to_csv(os.path.join(OUT, 'fwd_rem_sensitivity.csv'), index=False)
# per 1 % of R_EM at g_M = 2.97 (from g_E = -1.5 / -0.5: Delta R_EM = -+1.14 %)
a = t[(t.gM == C.GM0) & (t.gE == -1.5)].set_index('E'); b = t[(t.gM == C.GM0) & (t.gE == -0.5)].set_index('E')
sens = pd.DataFrame({'E': E, 'S2z0': s0, 'dS2z_per_1pct': 0.5 * (a.ds2z.abs() + b.ds2z.abs()).values / 1.14,
                     'sign(dS2z/dREM)': np.sign((a.ds2z - b.ds2z).values / (a.REM - b.REM).values),
                     'cs0': cs0, 'dcs_rel_per_1pct': 0.5 * (a.dcs_rel.abs() + b.dcs_rel.abs()).values / 1.14})
# band of the figures: envelope over the grid
g = t[t.apply(lambda r: (r.gM, round(r.gE, 3)) in {(x, round(y, 3)) for x, y in C.BAND_GRID}, axis=1)]
env = g.groupby('E').agg(s2z_lo=('s2z', 'min'), s2z_hi=('s2z', 'max'), cs_lo=('cs', 'min'), cs_hi=('cs', 'max'))
sens = sens.set_index('E').join(env)
sens['band_half_s2z'] = 0.5 * (sens.s2z_hi - sens.s2z_lo); sens['band_half_cs_rel'] = 0.5 * (sens.cs_hi - sens.cs_lo) / sens.cs0
pd.set_option('display.width', 200)
with open(os.path.join(OUT, 'fwd_rem_sensitivity.md'), 'w') as f:
    f.write('# Forward observables vs R_EM (NNLO prediction, own static polarisabilities)\n\n')
    f.write(sens.round(4).to_string() + '\n')
print(sens.round(3).to_string())
