#!/usr/bin/env python3
"""Cross-check db_old.csv against the author's legacy data embedded in xmgrace files
(dump produced by dump_agr_data.py) and the RCSdata/*.exp files.
Usage: python3 -I crosscheck_agr.py db_old.csv author_agr_embedded_data.txt RCSdata_dir
"""
import csv, math, os, sys

MP = 938.272
db = list(csv.DictReader(open(sys.argv[1])))
blocks = {}
for b in open(sys.argv[2]).read().split('\n\n'):
    lines = b.strip().split('\n')
    if not lines or not lines[0].startswith('# comment='):
        continue
    name = lines[0].split("'")[1]
    data = [list(map(float, l.split())) for l in lines[1:]]
    blocks.setdefault(os.path.basename(name), []).append(data)


def jac_cm_over_lab(E, th_lab):
    """dsigma/dOmega_cm / dsigma/dOmega_lab = (w'_cm / w'_lab)^2"""
    wl = E / (1 + E / MP * (1 - math.cos(math.radians(th_lab))))
    wc = E * MP / math.sqrt(MP * MP + 2 * MP * E)
    return (wc / wl) ** 2


def sel(ds, **kw):
    out = []
    for r in db:
        if r['dataset'] != ds:
            continue
        ok = True
        for k, v in kw.items():
            if abs(float(r[k]) - v) > 1e-6:
                ok = False
        if ok:
            out.append(r)
    return out


def compare(label, legacy, ds, fixkey, fixval, xkey, err_mode='stat'):
    """compare legacy [(x, v, e)] to db rows of dataset ds with r[fixkey]==fixval and x=r[xkey]"""
    rows = sel(ds, **{fixkey: fixval})
    nmatch = nval = 0
    msgs = []
    for x, v, e, *rest in legacy:
        cand = [r for r in rows if abs(float(r[xkey]) - x) < 0.051]
        if not cand:
            msgs.append(f'   {x}: no db point')
            continue
        r = min(cand, key=lambda r: abs(float(r['value']) - v))
        dv = float(r['value'])
        st = float(r['stat_err'] or 0)
        sy = float(r['sys_err_pt'] or 0)
        etot = math.hypot(st, sy)
        if abs(dv - v) < 1e-6:
            nval += 1
            if abs(st - e) < 1e-6:
                nmatch += 1
            else:
                msgs.append(f'   {x}: value ok, error legacy {e} vs db stat {st} (stat+sys {etot:.2f})')
        else:
            msgs.append(f'   {x}: VALUE legacy {v}+-{e} vs db {dv}+-{st}')
    print(f'{label}: {len(legacy)} legacy pts; values equal {nval}; values+errors equal {nmatch}')
    for m in msgs:
        print(m)


def first(name):
    return blocks[name][0] if name in blocks else []


print('=== Olmos de Leon (lab, fixed angle) ===')
for nm, th in (('a059.leon01', 59), ('a085.leon01', 85), ('a133.leon01', 133), ('Olmos_107.dat', 107), ('Olmos_155.dat', 155)):
    compare(nm, first(nm), 'OlmosdeLeon2001', 'theta_lab_deg', th, 'E_lab_MeV')
print('--- e_XXX.leon (fixed energy, x = theta_cm): compare with db lab value x Jacobian ---')
for nm in sorted(k for k in blocks if k.startswith('e_') and k.endswith('.leon')):
    for data in blocks[nm]:
        msgs = []
        for x, v, e in data:
            cand = [r for r in db if r['dataset'] == 'OlmosdeLeon2001' and abs(float(r['theta_cm_deg']) - x) < 0.02]
            if not cand:
                msgs.append(f'{x}:none'); continue
            r = cand[0]
            E = float(r['E_lab_MeV']); th = float(r['theta_lab_deg'])
            j = jac_cm_over_lab(E, th)
            msgs.append(f"E={E} thlab={th}: legacy {v}+-{e} vs db_lab*J {float(r['value']) * j:.2f}+-{float(r['stat_err']) * j:.2f}")
        print(nm, '(', len(data), 'pts )'); [print('   ', m) for m in msgs]

print('\n=== Hallin 1993 ===')
for nm, th in (('a051.sal93', 51), ('a090.sal93', 90), ('a141.sal93', 141)):
    compare(nm, first(nm), 'Hallin1993_TabII', 'theta_cm_deg', th, 'E_lab_MeV')
for nm, E in (('e149.sal93', 149), ('e182.sal93', 182), ('e230.sal93', 230), ('e286.sal93', 286)):
    for data in blocks.get(nm, []):
        compare(nm, data, 'Hallin1993_TabI', 'E_lab_MeV', E, 'theta_cm_deg')
for data in blocks.get('e135.sal93', []):
    print('e135.sal93 (no 135 MeV angular distribution in Hallin Table I; legacy built from Table II neighbours?):', data)
for data in blocks.get('a060.hal93', []):
    print('a060.hal93 (lab 60 deg) vs Hallin TabII theta_cm=67 x (lab/cm) Jacobian:')
    for x, v, e in data:
        r = [r for r in sel('Hallin1993_TabII', theta_cm_deg=67) if abs(float(r['E_lab_MeV']) - x) < 0.01][0]
        j = jac_cm_over_lab(x, float(r['theta_lab_deg']))
        print(f"   E={x}: legacy {v}+-{e};  db cm {r['value']}+-{r['stat_err']} -> lab {float(r['value']) / j:.2f}+-{float(r['stat_err']) / j:.2f} (theta_lab={r['theta_lab_deg']})")

print('\n=== MacGibbon 1995 ===')
for fn, th in (('sal_95_th90.exp', 90), ('sal_95_th135.exp', 135)):
    data = [list(map(float, l.split())) for l in open(os.path.join(sys.argv[3], fn)) if l.strip()]
    leg = [(d[0], d[1], d[3]) for d in data]
    tg = [x for x in leg if x[0] < 100]; ut = [x for x in leg if x[0] > 100]
    compare(fn + ' tagged', tg, 'MacGibbon1995_tagged', 'theta_lab_deg', th, 'E_lab_MeV')
    compare(fn + ' untagged', ut, 'MacGibbon1995_untagged', 'theta_lab_deg', th, 'E_lab_MeV')
d = first('a090.mac95')
compare('a090.mac95 tagged', [x for x in d if x[0] < 100], 'MacGibbon1995_tagged', 'theta_lab_deg', 90, 'E_lab_MeV')
compare('a090.mac95 untagged', [x for x in d if x[0] > 100], 'MacGibbon1995_untagged', 'theta_lab_deg', 90, 'E_lab_MeV')

print('\n=== Baranov 1974/75 ===')
compare('a090.bar74', first('a090.bar74'), 'Baranov1974_75_th90', 'theta_lab_deg', 90, 'E_lab_MeV')

print('\n=== Federspiel 1991 ===')
for fn, th in (('feder_91_th60.exp', 60), ('feder_91_th135.exp', 135)):
    data = [list(map(float, l.split())) for l in open(os.path.join(sys.argv[3], fn)) if l.strip()]
    rows = sel('Federspiel1991', theta_lab_deg=th)
    print(f'{fn} (RCSdata; constant error 1.44 -> apparently digitised from PRL Fig.2) vs db (legacy a0{th}.fed91):')
    for (x, v, _, e), r in zip(data, rows):
        print(f"   E {x} vs {r['E_lab_MeV']}: {v} vs {r['value']}  ratio db/exp = {float(r['value']) / v:.3f}")

print('\n=== Zieger 1992 ===')
data = [list(map(float, l.split())) for l in open(os.path.join(sys.argv[3], 'zieg_92.exp')) if l.strip()]
for x, v, w, e in data:
    r = sel('Zieger1992', E_lab_MeV=x)[0]
    print(f"   E={x}: zieg_92.exp {v}+-{e} (half-width {w}) vs db {r['value']}+-{r['stat_err']}")

print('\n=== LEGS 2001 and MAMI 90 deg ===')
compare('a90.legs01', first('a90.legs01'), 'Blanpied2001_LEGS', 'theta_cm_deg', 90, 'E_lab_MeV')
rows = [r for r in db if r['dataset'] == 'Blanpied2001_LEGS' and r['obs'] == 'S3' and abs(float(r['theta_cm_deg']) - 90) < 1e-6]
print('asym90.legs01 vs LEGS Table III 90 deg:')
for x, v, e in first('asym90.legs01'):
    r = [r for r in rows if abs(float(r['E_lab_MeV']) - x) < 0.01][0]
    flag = '' if (abs(float(r['value']) - v) < 1e-9 and abs(float(r['stat_err']) - e) < 1e-9) else '   <-- DIFF'
    print(f"   {x}: {v}+-{e} vs {r['value']}+-{r['stat_err']}{flag}")
rows286 = [r for r in db if r['dataset'] == 'Blanpied2001_LEGS' and abs(float(r['E_lab_MeV']) - 286.5) < 1e-6]
for nm, obs in (('e_286.legs', 'CS'), ('sigma286.legs01', 'S3')):
    print(nm, 'vs LEGS 286.5 MeV:')
    for x, v, e in first(nm):
        r = [r for r in rows286 if r['obs'] == obs and abs(float(r['theta_cm_deg']) - x) < 1e-6][0]
        flag = '' if (abs(float(r['value']) - v) < 1e-9 and abs(float(r['stat_err']) - e) < 1e-9) else '   <-- DIFF'
        print(f"   {x}: {v}+-{e} vs {r['value']}+-{r['stat_err']}{flag}")
print('asym90.legs96 (legacy) = LEGS 1996 PRL (superseded, not in db):', first('asym90.legs96'))
compare('a90.lags97 (labelled "LEGS 97" in FIGS/*.mgr legends)', first('a90.lags97'), 'Hunger1997_th90cm', 'theta_cm_deg', 90, 'E_lab_MeV')
