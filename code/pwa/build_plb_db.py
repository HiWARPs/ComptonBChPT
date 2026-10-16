"""Rebuild the PLB 2018 (Krupina-Lensky-Pascalutsa) Table-I database.
Values/errors: from the vector graphics of PLB Fig. 5 (digitize_plb.py), which shows exactly the 138
points of Table I in the 11 energy bins.  Photon energies: matched to db_old.csv (same experiment, angle
and value) where possible; for the sets missing from db_old (Oxley, Goldansky, Bernardini, Hyman/Pugh at
50 and 135 deg) the energy is assigned from the literature / bin (flag E_src).
Output: results/db_plb2018.csv"""
import numpy as np
import pandas as pd
import compton as C
import digitize_plb

F = digitize_plb.main().reset_index(drop=True)
db = pd.read_csv('../data/db_old.csv')
db = db[(db.obs == 'CS') & (db.E_lab_MeV <= 152)].copy()

# convert c.m. (Hallin) to lab, errors stat (+) sys as in the PLB
k = C.kin(db.E_lab_MeV.values, db.theta_cm_deg.values, 'cm')
lf = C.lab_factor(k)
iscm = (db.frame_of_cs == 'cm').values
db['lab_val'] = np.where(iscm, db.value * lf, db.value)
db['lab_err_tot'] = np.where(iscm, np.hypot(db.stat_err, db.sys_err_pt.fillna(0)) * lf, db.stat_err)
db['thl'] = np.where(iscm, np.degrees(k['thlab']), db.theta_lab_deg)

EXPMAP = {'OlmosdeLeon': ['OlmosdeLeon2001'], 'Federspiel': ['Federspiel1991'],
          'MacGibbon': ['MacGibbon1995_tagged', 'MacGibbon1995_untagged'],
          'Baranov': ['Baranov1974_75_th90', 'Baranov1974_75_th150'], 'Zieger': ['Zieger1992'],
          'Hallin': ['Hallin1993_TabII', 'Hallin1993_TabI'], 'Hyman': ['Hyman1959_legacy_partial'],
          'Pugh': ['Pugh1957_Frisch1967_legacy_partial']}
TEXKEY = {'OlmosdeLeon': 'OlmosdeLeon:2001zn', 'Federspiel': 'Federspiel:1991yd', 'MacGibbon': 'MacGibbon:1995in',
          'Baranov': 'Baranov:1974ec', 'Zieger': 'Zieger:1992jq', 'Hallin': 'Hallin:1993ft', 'Hyman': 'Hyman:1959zz',
          'Pugh': 'Pugh:1957zz', 'Oxley': 'Oxley:1958zz', 'Goldansky': 'Goldansky:1960zz', 'Bernardini': 'Bernardini:1960wya'}

rows = []
# per-panel calibration offset of the digitisation, from the OdL points (known exactly)
offs = {}
for b, g in F[F.exp == 'OlmosdeLeon'].groupby('bin'):
    dd = []
    for _, r in g.iterrows():
        m = db[(db.dataset == 'OlmosdeLeon2001') & (abs(db.theta_lab_deg - r.theta_lab) < 0.5) & (abs(db.value - r.value) < 0.02)]
        if len(m) == 1:
            dd.append(r.value - m.value.iloc[0])
    offs[b] = np.median(dd) if dd else 0.0
for i, r in F.iterrows():
    val = r.value - offs.get(r.bin, 0.0)
    E = np.nan; src = 'figure'; dset = ''; norm = np.nan
    for ds in EXPMAP.get(r.exp, []):
        cand = db[(db.dataset == ds) & (abs(db.thl - r.theta_lab) < 1.0) & (abs(db.lab_val - val) / val < 0.03)
                  & (abs(db.E_lab_MeV - r.bin) <= 10)]
        if len(cand):
            j = (abs(cand.lab_val - val) + abs(cand.E_lab_MeV - r.bin) * 1e-3).idxmin()
            E = cand.loc[j, 'E_lab_MeV']; src = 'db_old'; dset = ds; norm = cand.loc[j, 'norm_err_frac']
            break
    if np.isnan(E):
        if r.exp == 'Oxley':
            E = 60.0; src = 'PLB TableI (60 MeV)'
        elif r.exp == 'Goldansky':
            E = 55.0; src = 'PLB TableI (55 MeV)'
        elif r.exp == 'Bernardini':
            E = 120.0 if r.bin == 117 else 139.0; src = 'PLB TableI (120,139 MeV)'
        elif r.exp in ('Hyman', 'Pugh'):
            # energies of the 90-deg points of the same set in the same bin (same photon beam)
            same = [x for x in rows if x['exp'] == r.exp and x['bin'] == r.bin and x['E_src'] == 'db_old']
            if same:
                E = same[0]['E_lab']; src = 'assumed = 90deg point of same bin'
            else:
                E = {69: 65.0, 109: 105.0}.get(r.bin, float(r.bin)); src = 'assumed (bin)'
    rows.append(dict(exp=r.exp, texkey=TEXKEY[r.exp], bin=r.bin, E_lab=E, theta_lab=round(r.theta_lab, 3),
                     value=round(val, 4), err=round(r.err, 4), E_src=src, db_dataset=dset, norm_err_frac=norm,
                     clipped_bar=r.clipped))
D = pd.DataFrame(rows)
# Hyman/Pugh off-90 points processed before their 90-deg partner: second pass
for i, r in D.iterrows():
    if r.E_src == 'assumed (bin)':
        same = D[(D.exp == r.exp) & (D.bin == r.bin) & (D.E_src == 'db_old')]
        if len(same):
            D.loc[i, 'E_lab'] = same.E_lab.iloc[0]; D.loc[i, 'E_src'] = 'assumed = 90deg point of same bin'
D['theta_cm'] = np.round(np.degrees(C.kin(D.E_lab.values, D.theta_lab.values, 'lab')['thcm']), 3)
D = D.sort_values(['bin', 'exp', 'theta_lab']).reset_index(drop=True)
D.to_csv('results/db_plb2018.csv', index=False)
if __name__ == '__main__':
    pd.set_option('display.width', 220)
    print('digitisation offsets per bin [nb/sr]:', {k: round(v, 5) for k, v in offs.items()})
    print(D.groupby('exp').size().to_string())
    print(D.groupby('E_src').size().to_string())
    print(D[D.E_src != 'db_old'].to_string())
    print('N =', len(D))
