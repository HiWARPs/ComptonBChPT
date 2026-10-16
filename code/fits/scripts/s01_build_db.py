"""Step 1: assemble the fit database (db/fitdb.csv, db/fitdb_pb.csv) and the dataset tables."""
import os, sys, re
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
import numpy as np
import pandas as pd
import database
from common import DBDIR, RES, OLD_L, A2, A2PB, HIGS, FWD

os.makedirs(DBDIR, exist_ok=True); os.makedirs(RES, exist_ok=True)
df = database.build(per_beamtime=False)
df.to_csv(os.path.join(DBDIR, 'fitdb.csv'), index=False)
dfp = database.build(per_beamtime=True)
dfp.to_csv(os.path.join(DBDIR, 'fitdb_pb.csv'), index=False)
print('fitdb.csv: %d rows; fitdb_pb.csv: %d rows' % (len(df), len(dfp)))

ORDER = OLD_L + A2 + A2PB + HIGS + FWD + ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_S3', 'LEGS2001_CS',
                                          'Hunger1997_75', 'Hunger1997_90', 'Wissmann1999', 'Camen2002',
                                          'Hallin1993_TabII_D', 'Gray1967', 'Genzel1972', 'Baranov1966', 'Kolb2000',
                                          'Nagashima1964']
REF = {'OlmosdeLeon2001': 'deLeon:2001dnx', 'Federspiel1991': 'Federspiel:1991yd', 'Zieger1992': 'Zieger:1992jq',
       'Baranov1974_75_th90': 'Baranov:1974ec,Baranov:1975ju', 'Baranov1974_75_th150': 'Baranov:1974ec,Baranov:1975ju',
       'Hallin1993_TabII': 'Hallin:1993ft', 'MacGibbon1995_tagged': 'MacGibbon:1995in',
       'MacGibbon1995_untagged': 'MacGibbon:1995in', 'Pugh1957': 'Pugh:1957zz', 'Hyman1959': 'Hyman:1959zz',
       'Goldansky1960': 'Goldansky:1960', 'Oxley1958': 'Oxley:1958zz', 'Bernardini1960': 'Bernardini:1960wya',
       'A2_2022_CS': 'A2:2021wmm', 'A2_2022_S3': 'A2:2021wmm', 'HIGS_2022_CIRC': 'Li:2022vnz', 'HIGS_2022_LIN': 'Li:2022vnz',
       'A2_2018Mar_CS': 'Mornacchi thesis', 'A2_2018Jul_CS': 'Mornacchi thesis', 'A2_2018Mar_S3': 'Mornacchi thesis',
       'A2_2018Jul_S3': 'Mornacchi thesis', 'GHP_fwd_CS': 'Gryniuk:2015eza,Gryniuk:2016gnm',
       'GHP_fwd_S2Z': 'Gryniuk:2016gnm', 'A2_Martel2015_S2x': 'A2:2014iky', 'A2_Paudyal2020_S2z': 'A2:2019bqm',
       'LEGS2001_S3': 'Blanpied:2001ae', 'LEGS2001_CS': 'Blanpied:2001ae', 'Hunger1997_75': 'Hunger:1997pt',
       'Hunger1997_90': 'Hunger:1997pt', 'Wissmann1999': 'Wissmann:1999vi', 'Camen2002': 'Camen:2001st,Kossert:2002ws',
       'Hallin1993_TabII_D': 'Hallin:1993ft', 'Gray1967': 'Gray:1967zzd', 'Genzel1972': 'Genzel:1972ba',
       'Baranov1966': 'Baranov:1966aaa', 'Kolb2000': 'Kolb:2000ix', 'Nagashima1964': 'Nagashima:1964zz'}

def texlabel(x):
    x = re.sub(r'Σ(2x|2z|3)', lambda m: '$\\Sigma_{%s}$' % m.group(1), x)
    for a, b in [('dσ∥,⊥', 'd$\\sigma_{\\parallel,\\perp}$'), ('dσ', 'd$\\sigma$'), ('γ', '$\\gamma$'),
                 ('°cm', '$^\\circ_{\\rm cm}$'), ('°', '$^\\circ$'), ('é', "\\'e"), ('ü', '\\"u'), ('>', '$>$')]:
        x = x.replace(a, b)
    return x


allr = pd.concat([df, dfp[dfp.dataset.isin(A2PB)]])
rows = []
for ds in ORDER:
    d = allr[allr.dataset == ds]
    if not len(d):
        continue
    th = d['theta_lab']
    rows.append(dict(dataset=ds, expt=d['expt'].iloc[0], ref=REF.get(ds, ''), obs='/'.join(sorted(set(d['dobs']))),
                     N=len(d), N150=int((d.E <= 150).sum()), N170=int((d.E <= 170).sum()), N350=int((d.E <= 350).sum()),
                     E='%.1f-%.1f' % (d.E.min(), d.E.max()), th_lab='%.0f-%.0f' % (th.min(), th.max()),
                     frame=d['frame'].iloc[0], norm='%.1f%%' % (100 * d['norm'].iloc[0]) if d['norm'].iloc[0] > 0 else '-',
                     nflag=int(d['flags'].str.contains('GMPF_out').sum()), errmodel=d['prov'].iloc[0]))
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(RES, 'dataset_table.csv'), index=False)
with open(os.path.join(RES, 'dataset_table.md'), 'w') as f:
    f.write('| set | label | obs | N | N(E<=150) | N(E<=170) | N(E<=350) | E_lab [MeV] | theta_lab [deg] | angle frame | norm (1 sigma) | GMPF-pruned | error model |\n')
    f.write('|' + '---|' * 13 + '\n')
    for r in rows:
        f.write('| %s | %s | %s | %d | %d | %d | %d | %s | %s | %s | %s | %d | %s |\n' % (
            r['dataset'], r['expt'], r['obs'], r['N'], r['N150'], r['N170'], r['N350'], r['E'], r['th_lab'],
            r['frame'], r['norm'], r['nflag'], r['errmodel']))
with open(os.path.join(RES, 'dataset_table.tex'), 'w') as f:
    f.write('\\begin{tabular}{llrrrlll}\n\\hline\nExperiment & Obs. & $N_{150}$ & $N_{170}$ & $N_{350}$ & '
            '$E_\\gamma^{\\rm lab}$ [MeV] & $\\theta_{\\rm lab}$ [deg] & norm. \\\\\n\\hline\n')
    for r in rows:
        f.write('%s~\\cite{%s} & %s & %d & %d & %d & %s & %s & %s \\\\\n' % (
            texlabel(r['expt']),
            r['ref'], r['obs'], r['N150'], r['N170'], r['N350'], r['E'].replace('-', '--'), r['th_lab'].replace('-', '--'),
            r['norm'].replace('%', '\\%')))
    f.write('\\hline\n\\end{tabular}\n')
print(tab[['dataset', 'obs', 'N', 'N150', 'N170', 'N350', 'E', 'frame', 'norm', 'nflag']].to_string())
