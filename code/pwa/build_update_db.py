"""Updated proton RCS database for the PWA below 150 MeV.
Sources (read-only): work/data/db_old.csv (TAPS 2001), db_baranov2001.csv (all pre-2001 sets below 150 MeV: Table 2
of Baranov, L'vov, Petrun'kin, Shtarkov 2001, lab cross sections, statistical errors, normalisation uncertainty delta;
since 2026-10-19, replacing db_old/db_klp2018 for these sets), db_new.csv (A2 2022, HIGS 2022),
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
    parts = [rows_from(o[o.dataset == 'OlmosdeLeon2001'], 'TAPS2001', 'OlmosdeLeon:2001zn', 'old')]
    # Hallin Table II at 150 MeV (beyond Baranov et al.'s omega < 150): c.m. -> lab with the exact Jacobian, delta = 4%
    x = o[(o.dataset == 'Hallin1993_TabII') & (o.E_lab_MeV >= 150.0)].copy()
    k = C.kin(x.E_lab_MeV.values, x.theta_cm_deg.values, 'cm')
    lf = C.lab_factor(k)
    h = rows_from(x, 'Hallin1993', 'Hallin:1993ft', 'old', norm=0.04)
    h['theta_lab'] = np.degrees(k['thlab'])
    h['value'] *= lf; h['stat'] *= lf; h['sys'] = 0.0      # Baranov et al. error model for Hallin: stat (+) 4% scale
    parts.append(h)
    return pd.concat(parts, ignore_index=True)


def old_baranov2001():
    """pre-2001 sets from Baranov et al. 2001, Table 2.  The five oldest sets (Oxley, Goldansky, Bernardini, Pugh,
    Hyman) get norm 0 here (fixed-normalisation configuration, main_fits.py); the paper's PWA fit floats them with
    their delta (update_fits.LEGACY_NORM).  Federspiel 1991 only from 55 MeV on (lowest KLP18 bin)."""
    b = pd.read_csv(D + 'db_baranov2001.csv')
    b = b[(b.used == 1) & ~((b.dataset == 'Federspiel1991') & (b.E_lab_MeV < 54.9))].copy()
    b['sys_err_pt'] = 0.0
    name = {'Baranov1974_75_th90': 'Baranov1974_th90', 'Baranov1974_75_th150': 'Baranov1974_th150',
            'Hallin1993_TabII': 'Hallin1993', 'Oxley1958': 'Oxley', 'Goldansky1960': 'Goldansky',
            'Bernardini1960': 'Bernardini', 'Pugh1957': 'Pugh', 'Hyman1959': 'Hyman'}
    parts = []
    for ds in dict.fromkeys(b.dataset):
        x = b[b.dataset == ds]
        nm = name.get(ds, ds)
        parts.append(rows_from(x, nm, x.texkey.iloc[0], 'old',
                               norm=0.0 if nm in ('Oxley', 'Goldansky', 'Bernardini', 'Pugh', 'Hyman') else None))
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
    d = pd.concat([old_primary(), old_baranov2001(), new_data(per_beamtime)], ignore_index=True)
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
