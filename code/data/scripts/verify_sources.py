#!/usr/bin/env python3
"""Independent re-verification of db_old.csv rows against the raw sources in a second format
(HEPData YAML instead of CSV; arXiv TeX vs HEPData), to guard against parsing errors."""
import csv, os, re, sys
from decimal import Decimal
RAW = sys.argv[2]
db = list(csv.DictReader(open(sys.argv[1])))

def yaml_dep_values(path):
    """crude extraction of dependent variable 'value:' entries and first symerror from HEPData YAML"""
    t = open(path).read()
    dep = t.split('independent_variables:')[0]
    vals = re.findall(r'^\s*value: (\S+)\s*$', dep, flags=re.M)
    return [v for v in vals]

def D(x):
    return Decimal(str(x)).normalize()

checks = [('Hallin1993_HEPData', 'ins369337', range(1, 17), ('Hallin1993_TabI', 'Hallin1993_TabII'), 1),
          ('Hunger1997_HEPData', 'ins458618', range(1, 3), ('Hunger1997_th75cm', 'Hunger1997_th90cm'), 1),
          ('Wissmann1999_HEPData', 'ins514218', range(1, 2), ('Wissmann1999',), 1),
          ('Kossert2003_HEPData', 'ins599960', range(1, 2), ('Camen2002_Kossert2003',), 1),
          ('Genzel1972_HEPData', 'ins77720', range(1, 7), ('Genzel1972',), 1000),
          ('Gray1967_HEPData', 'ins944949', range(3, 4), ('Gray1967',), 1000),
          ('Baranov1975_HEPData', 'ins102948', range(1, 3), ('Baranov1974_75_th90', 'Baranov1974_75_th150'), 1000),
          ('Nagashima1964_HEPData', 'ins1101996', range(1, 2), ('Nagashima1964',), 1000),
          ('Baranov1966b_HEPData', 'ins1392560', range(3, 5), ('Baranov1966',), 1000),
          ('Baranov1966a_HEPData', 'ins1387582', range(1, 2), ('Baranov1966_refined',), 1000)]
for d, ins, tabs, dsets, fac in checks:
    ydir = os.path.join(RAW, d, 'yaml', f'HEPData-{ins}-v1-yaml')
    yv = []
    for i in tabs:
        yv += [D(Decimal(v) * fac) for v in yaml_dep_values(os.path.join(ydir, f'Table{i}.yaml'))]
    dbv = [D(r['value']) for r in db if r['dataset'] in dsets]
    missing = [v for v in dbv if v not in yv]
    print(f'{d}: db rows {len(dbv)}, yaml values {len(yv)}, db values not found in YAML: {missing}')

# MacGibbon: db (from arXiv TeX) vs HEPData CSV
hd = os.path.join(RAW, 'MacGibbon1995_HEPData', 'csv')
hv = []
for i in (1, 2):
    for l in open(os.path.join(hd, f'Table{i}.csv')):
        if l[:1].isdigit():
            p = l.strip().split(',')
            hv.append((D(p[0]), D(p[-3]), D(p[-2])))
dbm = [(D(r['E_lab_MeV']), D(r['value']), D(r['stat_err'])) for r in db if r['dataset'].startswith('MacGibbon')]
print('MacGibbon: db(arXiv TeX) rows', len(dbm), '; not in HEPData:', [x for x in dbm if x not in hv])

# Kossert: db (HEPData) vs arXiv TeX table
import tarfile
tex = None
src = os.path.join(RAW, 'Kossert2003_arXiv', 'source_nucl-ex_0210020.tgz')
with tarfile.open(src) as tf:
    for m in tf.getmembers():
        if m.name.endswith('paper.tex'):
            tex = tf.extractfile(m).read().decode('latin-1')
seg = tex[tex.find('label{csfree_p}'):tex.find('label{table_qffr_proton}')]
tv = [(D(a), D(b), D(c)) for a, b, c in re.findall(r'^\s*([\d.]+) &&\s*([\d.]+) & \$\\pm\$ &\s*([\d.]+)', seg, flags=re.M)]
dbk = [(D(r['E_lab_MeV']), D(r['value']), D(r['stat_err'])) for r in db if r['dataset'] == 'Camen2002_Kossert2003']
print('Kossert/Camen: db(HEPData) rows', len(dbk), '; TeX rows', len(tv), '; db rows not in TeX:', [x for x in dbk if x not in tv])
