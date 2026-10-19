"""Build db_baranov2001.csv from Table 2 of P.S. Baranov, A.I. L'vov, V.A. Petrun'kin, L.N. Shtarkov,
"Experimental status of the electric and magnetic polarizabilities of the proton", Fiz. Elem. Chast. Atom. Yadra 32
(2001) 699 [Phys. Part. Nucl. 32 (2001) 376] (raw/old/Baranov2001_PEPAN/), and the full error matrix of the
MacGibbon 1995 untagged cross sections (MacGibbon et al., PRC 52 (1995) 2097, appendix; arXiv TeX source).

Table 2 = all gamma p -> gamma p cross sections below 150 MeV used by Baranov et al. (lab frame, lab angles; the
c.m. values of Zieger 92, Hallin 93 and Bernardini 61 converted to the lab by them), statistical errors, and the
systematic (normalisation) uncertainty delta of each set; corrected values: Federspiel 91 (A.M. Nathan, priv. comm.),
Pugh 57 (final values, D. Frisch, priv. comm., "MIT Summer Study 1967"), Baranov 74 split into 90 deg (delta = 1%) and
150 deg (delta = 3.7%) sets.  Gov56 is listed but was not fitted by Baranov et al. (no systematic error) -> used=0.

Usage: python3 -I build_db_baranov2001.py <table2_pdftotext.txt> <MacGibbon 9507001.tex> <out.csv> <out_cov.csv>"""
import csv
import math
import re
import sys

MP = 938.272
LABEL = {'Gov56': 'Govorkov1956', 'Oxl58': 'Oxley1958', 'Hym59': 'Hyman1959', 'Gol60': 'Goldansky1960',
         'Ber61': 'Bernardini1960', 'Fri67': 'Pugh1957', 'Bar74a': 'Baranov1974_75_th90',
         'Bar74b': 'Baranov1974_75_th150', 'Fed91': 'Federspiel1991', 'Zie92': 'Zieger1992',
         'Hal93': 'Hallin1993_TabII', 'Mac95': 'MacGibbon1995'}
TEXKEY = {'Govorkov1956': '', 'Oxley1958': 'Oxley:1958zz', 'Hyman1959': 'Hyman:1959zz',
          'Goldansky1960': 'Goldansky:1960zz', 'Bernardini1960': 'Bernardini:1960wya', 'Pugh1957': 'Pugh:1957zz',
          'Baranov1974_75_th90': 'Baranov:1974ec', 'Baranov1974_75_th150': 'Baranov:1974ec',
          'Federspiel1991': 'Federspiel:1991yd', 'Zieger1992': 'Zieger:1992jq', 'Hallin1993_TabII': 'Hallin:1993ft',
          'MacGibbon1995_tagged': 'MacGibbon:1995in', 'MacGibbon1995_untagged': 'MacGibbon:1995in'}
NUM = r'\d+(?:,\d+)?'
CS = re.compile(r'(\d+),\s*(\d)\s*±\s*(\d+),\s*(\d)\s*$')


def num(s):
    return float(s.replace(',', '.'))


def th_lab2cm(E, thl):
    b = E / (E + MP)
    c = math.cos(math.radians(thl))
    return math.degrees(math.acos((c - b) / (1 - b * c)))


def parse_table2(fn):
    rows, cur, n, delta, E, nmac = [], None, None, None, None, 0
    for line in open(fn, encoding='utf-8'):
        m = CS.search(line)
        if not m:
            continue
        val, err = num(m.group(1) + ',' + m.group(2)), num(m.group(3) + ',' + m.group(4))
        pre = line[:m.start()]
        lab = re.match(r'\s*(Gov56|Oxl58|Hym59|Gol60|Ber61|Fri67|Bar74a|Bar74b|Fed91|Zie92|Hal93|Mac95)\s*\[\d+\]', pre)
        if lab:
            cur = LABEL[lab.group(1)]
            if cur == 'MacGibbon1995':
                nmac += 1
                cur += '_tagged' if nmac == 1 else '_untagged'
            pre = pre[lab.end():]
            pre = re.sub(r'^[\s∗*]*', ' ', pre)
            tok = re.findall(NUM, pre)
            if cur == 'Govorkov1956':           # delta printed as a dash
                n, delta, tok = int(tok[0]), None, tok[1:]
            else:
                n, delta, tok = int(tok[0]), num(tok[1]), tok[2:]
        else:
            pre = re.sub(r'\([^)]*\)', ' ', pre)
            tok = re.findall(NUM, pre)
        if len(tok) == 2:
            E, th = num(tok[0]), num(tok[1])
        elif len(tok) == 1:
            th = num(tok[0])
        else:
            raise ValueError('cannot parse: %r' % line)
        rows.append(dict(dataset=cur, n_set=n, delta_pct=delta, E=E, th=th, value=val, err=err))
    return rows


def macgibbon_cov(tex):
    """full 5x5 covariance of the untagged cross sections per detector: C = V^-1 diag(dsig'^2) V^-T,
    sigma' = V sigma (MacGibbon et al. 1995, appendix).  Returns {angle: (V, sigp, dsigp)}."""
    s = open(tex, encoding='latin-1').read()
    i0 = s.index('This transformation is written explicitly')
    blk = s[i0:s.index('%====', i0)]
    out = {}
    for ang, a, b in ((135, 'for the 135$^\\circ$ detector as', 'and for the 90$^\\circ$ detector as'),
                      (90, 'and for the 90$^\\circ$ detector as', 'The fitting was done')):
        part = blk[blk.index(a):blk.index(b)]
        V = []
        for ln in part.split('\\nonumber')[:5]:
            rhs = ln.split('=')[-1]
            c = re.findall(r'([+-]?)\s*([0-9.]+)\s*\\:\\sigma_(\d)', rhs.replace('& ', ''))
            row = [0.0] * 5
            for sg, v, k in c:
                row[int(k) - 1] = (-1 if sg == '-' else 1) * float(v)
            V.append(row)
        out[ang] = V
    t = blk[blk.index('The transformed cross sections'):]
    pairs = re.findall(r'=\s*([-0-9.]+)\s*\\pm\s*([0-9]+\.[0-9]+)', t.replace('&', ''))
    out = {135: (out[135], [float(p[0]) for p in pairs[:5]], [float(p[1]) for p in pairs[:5]]),
           90: (out[90], [float(p[0]) for p in pairs[5:10]], [float(p[1]) for p in pairs[5:10]])}
    return out


def inv(A):
    n = len(A)
    M = [list(r) + [1.0 if i == j else 0.0 for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [r[n:] for r in M]


def main(t2, tex, out, outcov):
    rows = parse_table2(t2)
    # self-checks: per-set counts equal the printed n; total fitted = 94 (46 early + 48 late, Baranov et al. Sect. 1)
    cnt = {}
    for r in rows:
        cnt[r['dataset']] = cnt.get(r['dataset'], 0) + 1
    for r in rows:
        assert cnt[r['dataset']] == r['n_set'], (r['dataset'], cnt[r['dataset']], r['n_set'])
    assert sum(v for k, v in cnt.items() if k != 'Govorkov1956') == 94, cnt
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['dataset', 'texkey', 'obs', 'E_lab_MeV', 'theta_lab_deg', 'theta_cm_deg', 'frame_of_cs', 'value',
                    'stat_err', 'norm_err_frac', 'used', 'units', 'source', 'notes'])
        for r in rows:
            used = 0 if r['dataset'] == 'Govorkov1956' else 1
            note = 'Gov56: not fitted by Baranov et al. (no systematic error quoted)' if not used else ''
            if r['dataset'] in ('Zieger1992', 'Hallin1993_TabII', 'Bernardini1960'):
                note = 'c.m. cross section of the original converted to the lab by Baranov et al. (Table 2 footnote)'
            w.writerow([r['dataset'], TEXKEY[r['dataset']], 'CS', '%g' % r['E'], '%g' % r['th'],
                        '%.3f' % th_lab2cm(r['E'], r['th']), 'lab', '%.1f' % r['value'], '%.1f' % r['err'],
                        '' if r['delta_pct'] is None else '%g' % (r['delta_pct'] / 100), used, 'nb/sr',
                        'Baranov, L\'vov, Petrun\'kin, Shtarkov, Fiz.Elem.Chast.Atom.Yadra 32 (2001) 699, Table 2', note])
    # MacGibbon untagged covariance
    mg = macgibbon_cov(tex)
    unt = {r['th']: [x for x in rows if x['dataset'] == 'MacGibbon1995_untagged' and x['th'] == r['th']]
           for r in rows if r['dataset'] == 'MacGibbon1995_untagged'}
    with open(outcov, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['dataset', 'theta_lab_deg', 'E_i_MeV', 'E_j_MeV', 'cov_nb2'])
        for ang, (V, sp, dsp) in mg.items():
            sig = [x['value'] for x in sorted(unt[ang], key=lambda x: x['E'])]
            Es = [x['E'] for x in sorted(unt[ang], key=lambda x: x['E'])]
            # check: V sigma reproduces the published transformed cross sections (rounding of sigma to 0.1 nb/sr)
            for j in range(5):
                assert abs(sum(V[j][i] * sig[i] for i in range(5)) - sp[j]) < 0.1, (ang, j)
            Vi = inv(V)
            C = [[sum(Vi[i][k] * dsp[k] ** 2 * Vi[j][k] for k in range(5)) for j in range(5)] for i in range(5)]
            diag = [math.sqrt(C[i][i]) for i in range(5)]
            print('MacGibbon untagged %d deg: sqrt(diag C) =' % ang, ' '.join('%.2f' % x for x in diag),
                  ' (Table 2 errors:', ' '.join('%.1f' % x['err'] for x in sorted(unt[ang], key=lambda x: x['E'])), ')')
            print('   correlation of adjacent bins:', ' '.join('%+.2f' % (C[i][i + 1] / diag[i] / diag[i + 1]) for i in range(4)))
            for i in range(5):
                for j in range(5):
                    w.writerow(['MacGibbon1995_untagged', ang, Es[i], Es[j], '%.6f' % C[i][j]])
    for k, v in cnt.items():
        print('%-24s %3d' % (k, v))
    print('wrote', out, len(rows), 'rows;', outcov)


if __name__ == '__main__':
    main(*sys.argv[1:5])
