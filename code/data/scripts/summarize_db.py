#!/usr/bin/env python3
"""Summarise db_old.csv per dataset/observable (markdown table)."""
import csv, sys
from collections import OrderedDict
MP, MPI0 = 938.272, 134.977
ETH = MPI0 + MPI0**2 / (2 * MP)
rows = list(csv.DictReader(open(sys.argv[1])))
g = OrderedDict()
for r in rows:
    g.setdefault((r['dataset'], r['obs']), []).append(r)
print(f'pi0 threshold E_lab = {ETH:.2f} MeV\n')
print('| dataset | texkey | obs | N | E_lab range (MeV) | theta_lab range | theta_cm range | frame | norm | N(E<170) | N(E<=Eth) | N(E<=150) | N(E<=350) |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
tot = [0, 0, 0, 0, 0]
for (ds, obs), rs in g.items():
    E = [float(r['E_lab_MeV']) for r in rs]
    tl = [float(r['theta_lab_deg']) for r in rs]
    tc = [float(r['theta_cm_deg']) for r in rs]
    norms = sorted({r['norm_err_frac'] for r in rs})
    n = [len(rs), sum(e < 170 for e in E), sum(e <= ETH for e in E), sum(e <= 150 for e in E), sum(e <= 350 for e in E)]
    tot = [a + b for a, b in zip(tot, n)]
    print(f"| {ds} | {rs[0]['texkey']} | {obs} | {n[0]} | {min(E):g}-{max(E):g} | {min(tl):.1f}-{max(tl):.1f} | {min(tc):.1f}-{max(tc):.1f} | {rs[0]['frame_of_cs']} | {','.join(x or '-' for x in norms)} | {n[1]} | {n[2]} | {n[3]} | {n[4]} |")
print(f'| TOTAL | | | {tot[0]} | | | | | | {tot[1]} | {tot[2]} | {tot[3]} | {tot[4]} |')
