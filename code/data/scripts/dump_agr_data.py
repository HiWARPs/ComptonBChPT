#!/usr/bin/env python3
"""Dump all experimental-looking data sets (with symbols, xydy types, data-file comments)
from the author's xmgrace files into one text file, de-duplicated by (comment, data).
Usage: python3 -I dump_agr_data.py out.txt file1 file2 ...
"""
import sys, os, importlib.util
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('parse_agr2', os.path.join(here, 'parse_agr2.py'))
pa = importlib.util.module_from_spec(spec); spec.loader.exec_module(pa)
out = open(sys.argv[1], 'w')
seen = {}
for fn in sys.argv[2:]:
    try:
        meta, sets = pa.parse(fn)
    except Exception as e:
        print('ERR', fn, e); continue
    for (g, s), typ, data in sets:
        md = meta.get((g, s), {})
        c = md.get('comment', '')
        sym = md.get('symbol', '?')
        isdata = typ.startswith('xyd') and (sym not in ('0',) or 'DATA' in c or '.exp' in c or 'leon' in c or 'sal' in c)
        if not isdata: continue
        if 'xcs_' in c or '.out' in c: continue
        key = (c, tuple(tuple(r) for r in data))
        seen.setdefault(key, []).append(f"{os.path.basename(fn)}:G{g}.S{s}")
for (c, data), where in sorted(seen.items(), key=lambda kv: kv[0][0]):
    out.write(f"# comment={c!r} n={len(data)} found_in={where}\n")
    for r in data:
        out.write('  ' + ' '.join(r) + '\n')
    out.write('\n')
print(len(seen), 'unique data sets')
