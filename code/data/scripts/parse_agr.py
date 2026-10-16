#!/usr/bin/env python3
"""Parse xmgrace .agr/.mgr files: list sets with comment/legend/type and embedded data.
Usage: python3 -I parse_agr.py file.agr [--data] [--filter substr]
Read-only on input; prints to stdout."""
import sys, re

def parse(fn):
    meta = {}  # (g,s) -> dict
    graphs_title = {}
    cur_g = None
    sets = []
    with open(fn, 'r', errors='replace') as f:
        lines = f.readlines()
    i = 0
    curtarget = None
    curtype = None
    while i < len(lines):
        ln = lines[i].rstrip('\n')
        m = re.match(r'^@g(\d+) on', ln)
        if m:
            cur_g = int(m.group(1))
        m = re.match(r'^@with g(\d+)', ln)
        if m:
            cur_g = int(m.group(1))
        m = re.match(r'^@\s+(title|subtitle)\s+"(.*)"', ln)
        if m and cur_g is not None:
            graphs_title.setdefault(cur_g, {})[m.group(1)] = m.group(2)
        m = re.match(r'^@\s+s(\d+)\s+(comment|legend)\s+"(.*)"', ln)
        if m and cur_g is not None:
            meta.setdefault((cur_g, int(m.group(1))), {})[m.group(2)] = m.group(3)
        m = re.match(r'^@\s+s(\d+)\s+(symbol|line type)\s+(\d+)', ln)
        if m and cur_g is not None:
            meta.setdefault((cur_g, int(m.group(1))), {})[m.group(2)] = m.group(3)
        m = re.match(r'^@target G(\d+)\.S(\d+)', ln)
        if m:
            curtarget = (int(m.group(1)), int(m.group(2)))
        m = re.match(r'^@type (\w+)', ln)
        if m and curtarget is not None:
            curtype = m.group(1)
            data = []
            i += 1
            while i < len(lines) and not lines[i].startswith('&'):
                data.append(lines[i].split())
                i += 1
            sets.append((curtarget, curtype, data))
            curtarget = None
        i += 1
    return meta, graphs_title, sets

if __name__ == '__main__':
    fn = sys.argv[1]
    showdata = '--data' in sys.argv
    flt = None
    if '--filter' in sys.argv:
        flt = sys.argv[sys.argv.index('--filter') + 1]
    meta, gt, sets = parse(fn)
    for (g, s), typ, data in sets:
        md = meta.get((g, s), {})
        if flt and flt not in md.get('comment', '') and flt not in typ:
            continue
        print(f"G{g}.S{s} type={typ} n={len(data)} sym={md.get('symbol','?')} comment={md.get('comment','')!r} legend={md.get('legend','')!r} gtitle={gt.get(g,{})}")
        if showdata:
            for d in data:
                print('   ', ' '.join(d))
