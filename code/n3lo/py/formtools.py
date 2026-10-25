"""formtools.py -- run FORM on a generated program and turn the printed projections into Python callables."""
import os, re, subprocess, hashlib
FORMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'form')
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'out')

def run_form(name, body, cache=True, header='common.h', bracket=None):
    """body: FORM statements after '#include common.h' (must define KSUB, GAM, call Proj/Finish).
    Returns dict {'V1': text, ...} of the printed expressions (Fortran format, unbracketed)."""
    br = f'Bracket {bracket};\n' if bracket else ''
    src = f'#include {header}\n' + body + f'\nid M^n1? = 1;\n{br}.sort\nFormat doublefortran;\nFormat 250;\nPrint;\n.end\n'
    key = src if header == 'common.h' else src + open(os.path.join(FORMDIR, header)).read()   # header changes invalidate
    h = hashlib.md5(key.encode()).hexdigest()[:10]
    fn = os.path.join(OUTDIR, f'{name}.frm'); fo = os.path.join(OUTDIR, f'{name}.log')
    if cache and os.path.exists(fo) and open(fo).read().startswith(f'#hash {h}\n'):
        txt = open(fo).read()
    else:
        open(fn, 'w').write(src)
        exe = os.environ.get('FORM_BIN', 'form').split()        # e.g. FORM_BIN='tform -w4' for the heavy graphs
        r = subprocess.run(exe + ['-q', '-p', FORMDIR, fn], capture_output=True, text=True, cwd=OUTDIR)
        if r.returncode != 0 or 'rror' in r.stdout:
            raise RuntimeError(f'FORM failed for {name}:\n' + r.stdout[-3000:] + r.stderr[-2000:])
        txt = f'#hash {h}\n' + r.stdout
        open(fo, 'w').write(txt)
    return parse_fortran(txt)

def parse_fortran(txt):
    lines = [l for l in txt.split('\n')[1:] if l.strip()]
    stmts, cur = [], ''
    for l in lines:
        if l.startswith('     &'):
            cur += l[6:].strip()
        else:
            if cur: stmts.append(cur)
            cur = l.strip()
    if cur: stmts.append(cur)
    out = {}
    for s in stmts:
        m = re.match(r'(\w+)\s*=\s*(.*)$', s)
        if not m: continue
        nm, rhs = m.group(1), m.group(2)
        if re.match(nm + r'\s*[+-]', rhs): rhs = re.sub(r'^' + nm + r'\s*', '', rhs)
        out.setdefault(nm, []).append(rhs)
    return {k: '+'.join('(' + r + ')' for r in v) for k, v in out.items()}

_num = re.compile(r'(\d+)\.D0')
def to_python(expr):
    e = _num.sub(r'\1.0', expr)
    e = re.sub(r'TCF\((\d+),(\d+),(\d+),(\d+)\)', r'T[(\1,\2,\3,\4)]', e)
    e = re.sub(r'UCF\((\d+),(\d+),(\d+),(\d+)\)', r'U[(\1,\2,\3,\4)]', e)
    e = re.sub(r'TCFT\((\d+),(\d+),(\d+),(\d+),(\d+)\)', r'T[(\1,\2,\3,\4,\5)]', e)
    e = re.sub(r'UCFT\((\d+),(\d+),(\d+),(\d+),(\d+)\)', r'U[(\1,\2,\3,\4,\5)]', e)
    return e

def compile_proj(d):
    """d: dict V1..V6 -> code objects evaluated with names x,t,mpi2,ixs,iys,iPq,PP,T,U, and couplings."""
    return [compile(to_python(d.get(f'V{j}', '0')), f'V{j}', 'eval') for j in range(1, 7)]
