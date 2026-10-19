"""Point-by-point comparison of Baranov et al. 2001, Table 2 (db_baranov2001.csv) with the pre-2001 sets of the fit
database as it was before the switch (fit DB rows from db_old.csv / db_klp2018.csv, E <= 150 MeV).
Hallin 93 (c.m. in our files) is compared after the exact c.m. -> lab conversion.
Usage: python compare_baranov2001.py <db_baranov2001.csv> <old fitdb.csv>  (prints a markdown report)"""
import sys
import numpy as np
import pandas as pd

MP = 938.272


def lab_factor(E, thl):
    """dsigma_lab / dsigma_cm at lab angle thl"""
    wl = E / (1 + E * (1 - np.cos(np.radians(thl))) / MP)
    wc = E * MP / np.sqrt(MP ** 2 + 2 * MP * E)
    return (wl / wc) ** 2


b = pd.read_csv(sys.argv[1])
b = b[b.used == 1]
f = pd.read_csv(sys.argv[2], keep_default_na=False)
f = f[f.dataset.isin(b.dataset.unique()) & (f.E <= 150.0001)].copy()
f['val_lab'] = np.where(f.dobs == 'CSCM', f.value * lab_factor(f.E, f.theta_lab), f.value)
f['stat_lab'] = np.where(f.dobs == 'CSCM', f.err_stat * lab_factor(f.E, f.theta_lab), f.err_stat)
used = set()
print('| set | E (B01 / ours) | θ_lab (B01 / ours) | dσ B01 | dσ ours | err B01 | err ours (stat; err_pt) | remark |')
print('|---|---|---|---|---|---|---|---|')
for ds, g in b.groupby('dataset', sort=False):
    h = f[f.dataset == ds]
    for r in g.itertuples():
        c = h[(abs(h.theta_lab - r.theta_lab_deg) < 0.6) & (abs(h.E - r.E_lab_MeV) < 6) & ~h.pid.isin(used)]
        if len(c) > 1:
            c = c.iloc[[np.argmin(abs(c.E - r.E_lab_MeV) + 0.01 * abs(c.val_lab - r.value))]]
        if not len(c):
            print('| %s | %g / – | %g / – | %.1f | – | %.1f | – | **missing in ours** |' % (ds, r.E_lab_MeV, r.theta_lab_deg, r.value, r.stat_err))
            continue
        q = c.iloc[0]
        used.add(q.pid)
        rem = []
        if abs(q.E - r.E_lab_MeV) > 0.05:
            rem.append('energy')
        if abs(q.theta_lab - r.theta_lab_deg) > 0.15:
            rem.append('angle')
        if abs(q.val_lab - r.value) > 0.051:
            rem.append('value %+.1f%%' % (100 * (r.value / q.val_lab - 1)))
        if abs(q.stat_lab - r.stat_err) > 0.051:
            rem.append('error')
        if rem:
            print('| %s | %g / %g | %g / %.2f | %.1f | %.2f | %.1f | %.2f; %.2f | %s |' % (
                ds, r.E_lab_MeV, q.E, r.theta_lab_deg, q.theta_lab, r.value, q.val_lab, r.stat_err, q.stat_lab,
                q.err_pt, ', '.join(rem)))
    print('| %s | N = %d / %d | | | | δ = %.1f%% / norm %.1f%% | | identical points not listed |' % (
        ds, len(g), len(h), 100 * g.norm_err_frac.iloc[0], 100 * float(h.norm.iloc[0]) if len(h) else 0))
rest = f[~f.pid.isin(used)]
for q in rest.itertuples():
    print('| %s | – / %g | – / %.2f | – | %.2f | – | %.2f; %.2f | **not in Baranov 2001** |' % (q.dataset, q.E, q.theta_lab, q.val_lab, q.stat_lab, q.err_pt))
