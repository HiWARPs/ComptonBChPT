#!/usr/bin/env python3
"""Build db_old_ext450.csv: PLOTTING-ONLY extension of db_old.csv to 400 < E_lab <= 450 MeV.

db_old.csv is cut at E_lab <= 400 MeV (build_db_old.py, EMAX = 400) and feeds the fits; it is NOT modified.
This script runs the very same transcription code (build_db_old.py, all builder functions unchanged) with the
energy cut lifted, checks that the rows with E_lab <= 400 MeV reproduce db_old.csv exactly (string identity of
every field), and writes only the rows with 400 < E_lab <= 450 MeV of the data sets shown in the Delta-region
figures (work/delta_plots), with the same columns as db_old.csv.  The old bremsstrahlung sets whose frames are
unknown (Gray 1967, Genzel 1972, Baranov 1966 [+refined], Nagashima 1964) are excluded (user decision
2026-10-14), as in work/delta_plots/config.EXCLUDE_SETS.

Used only by work/delta_plots (figure of dsigma/dOmega vs E_gamma at fixed angles, axis up to 450 MeV);
not used in any fit or chi^2 table.

Run (from anywhere):  python3 -I build_db_old_ext450.py <raw_old_dir> <db_old.csv> <out_csv>
  e.g. cd work/data && python3 -I scripts/build_db_old_ext450.py raw/old db_old.csv db_old_ext450.csv
"""
import csv
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
RAW, DBOLD, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
ELO, EHI = 400.0, 450.0
KEEP = ('Hunger1997_th75cm', 'Hunger1997_th90cm', 'Wissmann1999', 'Camen2002_Kossert2003', 'Blanpied2001_LEGS',
        'Hallin1993_TabII', 'MacGibbon1995_tagged', 'MacGibbon1995_untagged', 'OlmosdeLeon2001', 'Kolb2000',
        'Barbiellini1968')
EXCLUDED = ('Gray1967', 'Genzel1972', 'Baranov1966', 'Baranov1966_refined', 'Nagashima1964')   # frames unknown

# ---- run build_db_old.py unchanged except for the energy cut
src = open(os.path.join(HERE, 'build_db_old.py')).read()
cut = 'EMAX = 400.0   # keep points with E_lab <= EMAX'
assert src.count(cut) == 1, 'EMAX line of build_db_old.py changed - update this script'
src = src.replace(cut, 'EMAX = %r   # lifted by build_db_old_ext450.py' % EHI)
with tempfile.TemporaryDirectory() as tmp:
    ns = {'__name__': '__build_db_old__', '__file__': os.path.join(HERE, 'build_db_old.py')}
    argv = sys.argv
    sys.argv = ['build_db_old.py', RAW, os.path.join(tmp, 'all.csv')]
    try:
        exec(compile(src, 'build_db_old.py', 'exec'), ns)
    finally:
        sys.argv = argv
    HEADER, rows = ns['HEADER'], ns['rows']

# ---- consistency check: E <= 400 MeV part == db_old.csv (field by field, as strings)
with open(DBOLD, newline='') as f:
    old = list(csv.reader(f))
low = [r for r in rows if float(r[3]) <= ELO]
assert old[0] == HEADER, 'header differs from db_old.csv'
assert [list(map(str, r)) for r in low] == old[1:], 'E<=400 MeV rows do not reproduce db_old.csv'

# ---- the extension
ext = []
for r in rows:
    E = float(r[3])
    if ELO < E <= EHI and r[0] in KEEP:
        r = list(r)
        r[14] = r[14] + '; PLOTTING-ONLY extension 400<E<=450 MeV (build_db_old_ext450.py); not in db_old.csv, not used in fits'
        ext.append(r)
skipped = sorted({r[0] for r in rows if ELO < float(r[3]) <= EHI and r[0] not in KEEP})
with open(OUT, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(HEADER)
    w.writerows(ext)
print('db_old.csv reproduced for E <= %g MeV (%d rows); wrote %d rows with %g < E <= %g MeV to %s'
      % (ELO, len(low), len(ext), ELO, EHI, OUT))
print('sets in the extension:', sorted({r[0] for r in ext}))
print('sets with points in the range but not kept:', skipped, '(excluded: frames unknown)' if set(skipped) <= set(EXCLUDED) else '')
