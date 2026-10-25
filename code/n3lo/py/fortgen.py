"""fortgen.py -- FORM-optimised Fortran for the projections of one graph (V_1..V_n)."""
import os, re, subprocess, hashlib, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import formtools as ft
FDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'fortran', 'gen')
os.makedirs(FDIR, exist_ok=True)

def tsym(kind, key, uv=False):
    n2, a, b, c = key
    idx = '00' * (n2 // 2) + '1' * a + '2' * b + '3' * c
    pre = {'A': 'aa', 'B': 'bb', 'C': 'cc', 'D': 'dd', 'T': 'tt'}[kind]
    return ('U' if uv else 'T') + pre + (idx if idx else '0')

def optimise(g, kind, nproj, gi=False, tag='', level='O2'):
    """run FORM on the graph with TCF/UCF -> named symbols and O3 optimisation; returns Fortran text"""
    body, ks = g.form_body(gi=gi)
    keys = [(n2, a, b, c) for n2 in range(0, 9, 2) for a in range(6) for b in range(6) for c in range(6) if n2//2*2 + a + b + c <= 6]
    names = sorted({tsym(kind, k) for k in keys} | {tsym(kind, k, True) for k in keys})
    decl = 'Symbols ' + ','.join(names) + ';\n' if kind != 'T' else ''
    ids = ''.join(f'id TCF({k[0]},{k[1]},{k[2]},{k[3]}) = {tsym(kind, k)};\nid UCF({k[0]},{k[1]},{k[2]},{k[3]}) = {tsym(kind, k, True)};\n'
                  for k in keys) if kind != 'T' else 'id TCF(0,0,0,0) = 1;\n'
    src = '#include common.h\n' + decl + body + 'id M^n1? = 1;\n' + ids + '.sort\n'
    src += f'Format {level},stats=off;\nFormat doublefortran;\nFormat 120;\n'
    out = os.path.join(FDIR, g.name + tag + '.f90.inc')
    for j in range(1, nproj + 1):
        src += f'#optimize V{j}\n#write <{out}> "%O"\n#write <{out}> "      W({j}) = %e", V{j}\n#clearoptimize\n'
    src += '.end\n'
    h = hashlib.md5(src.encode()).hexdigest()[:10]
    stamp = out + '.hash'
    if os.path.exists(out) and os.path.exists(stamp) and open(stamp).read() == h:
        return open(out).read()
    if os.path.exists(out): open(out, 'w').close()
    fn = os.path.join(FDIR, g.name + tag + '.frm'); open(fn, 'w').write(src)
    r = subprocess.run(['form', '-q', '-p', ft.FORMDIR, fn], capture_output=True, text=True, cwd=FDIR)
    if r.returncode != 0: raise RuntimeError(r.stdout[-3000:])
    open(stamp, 'w').write(h)
    return open(out).read()
