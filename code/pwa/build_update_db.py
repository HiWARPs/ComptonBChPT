"""Updated proton RCS database for the PWA below 150 MeV.
Sources (read-only): work/data/db_old.csv (primary pre-2015 values, true energies), db_klp2018.csv (only source for
Oxley, Goldansky, Bernardini and the full Pugh/Hyman sets), db_new.csv (A2 2022, HIGS 2022),
db_mornacchi_thesis.csv (per-beamtime A2 variant).
Columns: dataset, texkey, obs (CS = lab dsigma/dOmega, CS_PAR, CS_PERP (lab), S3), E, theta_lab, value, stat, sys,
err (= stat (+) sys, used in chi^2), norm_group, norm_err (Gaussian-equivalent fractional scale), group (old/new)
Output: results/db_update.csv (combined A2) and results/db_update_A2beamtimes.csv (per-beamtime A2)."""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import numpy as np
import pandas as pd
import compton as C

D = (_WORK + '/data/')
EMAX = 150.5


def rows_from(df, dataset, texkey, group, obs=None, norm_group=None, norm=None, E=None):
    out = pd.DataFrame(dict(dataset=dataset, texkey=texkey, obs=df.obs.values if obs is None else obs,
                            E=df.E_lab_MeV.values if E is None else E, theta_lab=df.theta_lab_deg.values,
                            value=df.value.values.astype(float), stat=df.stat_err.values.astype(float),
                            sys=df.sys_err_pt.fillna(0).values.astype(float)))
    out['norm_group'] = norm_group if norm_group else dataset
    out['norm_err'] = norm if norm is not None else df.norm_err_frac.fillna(0).values.astype(float)
    out['group'] = group
    return out


def old_primary():
    o = pd.read_csv(D + 'db_old.csv')
    o = o[(o.obs == 'CS') & (o.E_lab_MeV <= EMAX)]
    parts = []
    x = o[o.dataset == 'OlmosdeLeon2001']
    parts.append(rows_from(x, 'TAPS2001', 'OlmosdeLeon:2001zn', 'old'))
    x = o[(o.dataset == 'Federspiel1991') & (o.E_lab_MeV >= 54.9)]
    parts.append(rows_from(x, 'Federspiel1991', 'Federspiel:1991yd', 'old'))
    x = o[o.dataset == 'Zieger1992']
    parts.append(rows_from(x, 'Zieger1992', 'Zieger:1992jq', 'old'))
    for ds in ('Baranov1974_75_th90', 'Baranov1974_75_th150'):
        x = o[o.dataset == ds]
        parts.append(rows_from(x, ds.replace('1974_75', '1974'), 'Baranov:1974ec', 'old', norm_group='Baranov1974'))
    for ds in ('MacGibbon1995_tagged', 'MacGibbon1995_untagged'):
        x = o[o.dataset == ds]
        parts.append(rows_from(x, ds, 'MacGibbon:1995in', 'old'))
    # Hallin Table II (c.m. frame) -> lab cross sections, errors stat, sys converted with the exact Jacobian
    x = o[(o.dataset == 'Hallin1993_TabII')].copy()
    k = C.kin(x.E_lab_MeV.values, x.theta_cm_deg.values, 'cm')
    lf = C.lab_factor(k)
    h = rows_from(x, 'Hallin1993', 'Hallin:1993ft', 'old')
    h['theta_lab'] = np.degrees(k['thlab'])
    h['value'] *= lf; h['stat'] *= lf; h['sys'] *= lf
    parts.append(h)
    return pd.concat(parts, ignore_index=True)


def old_klp_only():
    """Oxley, Goldansky, Bernardini, Pugh, Hyman from the KLP18 notebook (errors = total errors used there).
    Pugh/Hyman energies: replaced by the true energies of db_old for groups whose 90-deg point matches."""
    k = pd.read_csv(D + 'db_klp2018.csv')
    o = pd.read_csv(D + 'db_old.csv')
    parts = []
    for name, tk in (('Oxley', 'Oxley:1958zz'), ('Goldansky', 'Goldansky:1960zz'), ('Bernardini', 'Bernardini:1960wya')):
        x = k[k.dataset == 'KLP2018_' + name]
        parts.append(rows_from(x, name, tk, 'old', norm=0.0))
    for name, tk, ods in (('Pugh', 'Pugh:1957zz', 'Pugh1957_Frisch1967_legacy_partial'),
                          ('Hyman', 'Hyman:1959zz', 'Hyman1959_legacy_partial')):
        x = k[k.dataset == 'KLP2018_' + name].copy()
        ref = o[o.dataset == ods]
        Enew = x.E_lab_MeV.values.astype(float).copy()
        for Eg in np.unique(x.E_lab_MeV.values):
            g = x[x.E_lab_MeV == Eg]
            g90 = g[abs(g.theta_lab_deg - 90) < 0.5]
            if len(g90):
                m = ref[(abs(ref.value - g90.value.iloc[0]) < 1e-6) & (abs(ref.stat_err - g90.stat_err.iloc[0]) < 1e-6)
                        & (abs(ref.E_lab_MeV - Eg) <= 5)]
                if len(m):
                    Enew[(x.E_lab_MeV == Eg).values] = m.E_lab_MeV.iloc[0]
        parts.append(rows_from(x, name, tk, 'old', E=Enew, norm=0.0))
    return pd.concat(parts, ignore_index=True)


def new_data(per_beamtime=False):
    n = pd.read_csv(D + 'db_new.csv')
    parts = []
    if not per_beamtime:
        x = n[n.dataset == 'A2_Mornacchi2022_CS']
        parts.append(rows_from(x, 'A2_2022_CS', 'A2:2021wmm', 'new'))
        x = n[n.dataset == 'A2_Mornacchi2022_S3']
        parts.append(rows_from(x, 'A2_2022_S3', 'A2:2021wmm', 'new', norm=0.05 / np.sqrt(3)))
    else:
        t = pd.read_csv(D + 'db_mornacchi_thesis.csv')
        for ds in sorted(t.dataset.unique()):
            x = t[t.dataset == ds]
            nm = ds.replace('A2_Mornacchi2021thesis_', 'A2bt_')
            norm = 0.05 / np.sqrt(3) if '_S3_' in ds else None
            parts.append(rows_from(x, nm, 'Mornacchi:2021thesis', 'new', norm=norm))
    x = n[n.dataset == 'HIGS_Li2022_CIRC']
    parts.append(rows_from(x, 'HIGS_2022_circ', 'Li:2022vnz', 'new'))
    x = n[n.dataset == 'HIGS_Li2022_LIN']
    parts.append(rows_from(x, 'HIGS_2022_lin', 'Li:2022vnz', 'new'))
    return pd.concat(parts, ignore_index=True)


def build(per_beamtime=False):
    d = pd.concat([old_primary(), old_klp_only(), new_data(per_beamtime)], ignore_index=True)
    d['err'] = np.hypot(d.stat, d.sys)
    d['theta_cm'] = np.degrees(C.kin(d.E.values, d.theta_lab.values, 'lab')['thcm'])
    return d


if __name__ == '__main__':
    d = build()
    d.to_csv('results/db_update.csv', index=False)
    build(True).to_csv('results/db_update_A2beamtimes.csv', index=False)
    pd.set_option('display.width', 200)
    s = d.groupby(['group', 'dataset', 'obs']).agg(N=('value', 'size'), Emin=('E', 'min'), Emax=('E', 'max'),
                                                   thmin=('theta_lab', 'min'), thmax=('theta_lab', 'max'),
                                                   norm=('norm_err', 'first'))
    print(s.to_string())
    print('total', len(d), ' old', (d.group == 'old').sum(), ' new', (d.group == 'new').sum())
