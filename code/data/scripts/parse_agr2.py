#!/usr/bin/env python3
"""Parse xmgrace .agr and old xmgr .mgr files. Lists data sets with comment and data.
Usage: python3 -I parse_agr2.py file [--data] [--filter substr]"""
import sys, re
def parse(fn):
    meta = {}; sets = []
    lines = open(fn, errors='replace').read().splitlines()
    cur_g = None; curtarget = None; datag = None
    i = 0
    while i < len(lines):
        ln = lines[i]
        m = re.match(r'^@(?:with|WITH) [gG](\d+)', ln)
        if m: cur_g = int(m.group(1)); datag = cur_g
        m = re.match(r'^@\s+s(\d+)\s+(comment|legend)\s+"(.*)"', ln)
        if m and cur_g is not None:
            meta.setdefault((cur_g, int(m.group(1))), {})[m.group(2)] = m.group(3)
        m = re.match(r'^@\s+legend string (\d+) "(.*)"', ln)
        if m and cur_g is not None:
            meta.setdefault((cur_g, int(m.group(1))), {})['legend'] = m.group(2)
        m = re.match(r'^@\s+s(\d+) symbol (\d+)', ln)
        if m and cur_g is not None:
            meta.setdefault((cur_g, int(m.group(1))), {})['symbol'] = m.group(2)
        m = re.match(r'^@(?:target|TARGET) (?:G(\d+)\.)?S(\d+)', ln)
        if m:
            g = int(m.group(1)) if m.group(1) is not None else datag
            curtarget = (g, int(m.group(2)))
        m = re.match(r'^@(?:type|TYPE) (\w+)', ln)
        if m and curtarget is not None:
            typ = m.group(1); data = []; i += 1
            while i < len(lines) and not lines[i].startswith('&'):
                data.append(lines[i].split()); i += 1
            sets.append((curtarget, typ, data)); curtarget = None
        i += 1
    return meta, sets
if __name__ == '__main__':
    fn = sys.argv[1]; showdata = '--data' in sys.argv
    flt = sys.argv[sys.argv.index('--filter') + 1] if '--filter' in sys.argv else None
    meta, sets = parse(fn)
    for (g, s), typ, data in sets:
        md = meta.get((g, s), {})
        c = md.get('comment', '')
        if flt and flt not in c: continue
        print(f"G{g}.S{s} type={typ} n={len(data)} sym={md.get('symbol','?')} comment={c!r} legend={md.get('legend','')!r}")
        if showdata:
            for d in data: print('   ', ' '.join(d))
