#!/usr/bin/env python3
"""Build db_old.csv: pre-2015 world data on real Compton scattering off the proton,
E_lab ~30-400 MeV, from the raw sources downloaded into work/data/raw/old/.

Run (from anywhere):  python3 -I build_db_old.py <raw_old_dir> <out_csv>

Values are copied exactly as tabulated in the sources (only unit scaling mub->nb or
cm^2->nb where the source uses other units; this is noted per row). Cross sections are
NOT converted between lab and c.m. frames; only the angle in the other frame is computed
from exact two-body Compton kinematics with M_p = 938.272 MeV.
"""
import csv, math, os, re, sys
from decimal import Decimal

MP = 938.272
RAW = sys.argv[1]
OUT = sys.argv[2]
EMAX = 400.0   # keep points with E_lab <= EMAX

HEADER = ['dataset', 'texkey', 'obs', 'E_lab_MeV', 'dE_MeV', 'theta_lab_deg', 'theta_cm_deg',
          'frame_of_cs', 'value', 'stat_err', 'sys_err_pt', 'norm_err_frac', 'units', 'source', 'notes']
rows = []


def beta_cm(E):
    return E / (E + MP)


def lab2cm(E, th_lab):
    b = beta_cm(E)
    c = math.cos(math.radians(th_lab))
    cc = (c - b) / (1.0 - b * c)
    return math.degrees(math.acos(max(-1.0, min(1.0, cc))))


def cm2lab(E, th_cm):
    b = beta_cm(E)
    c = math.cos(math.radians(th_cm))
    cl = (c + b) / (1.0 + b * c)
    return math.degrees(math.acos(max(-1.0, min(1.0, cl))))


def fnum(x):
    """format a computed angle"""
    return f"{x:.3f}"


def scale(s, factor):
    """exact decimal scaling of a numeric string (unit conversion only)."""
    d = (Decimal(str(s)) * Decimal(str(factor))).normalize()
    t = format(d, 'f')
    return t


def add(dataset, texkey, obs, E, th, frame_angle, frame_cs, value, stat, sys_pt='', norm='',
        units='nb/sr', source='', notes='', dE=''):
    E_f = float(E)
    if E_f > EMAX:
        return
    th_f = float(th)
    if frame_angle == 'lab':
        th_lab, th_cm = str(th), fnum(lab2cm(E_f, th_f))
        if th_f in (0.0, 180.0):
            th_cm = str(th)
        anote = 'orig_angle=lab'
    else:
        th_cm, th_lab = str(th), fnum(cm2lab(E_f, th_f))
        if th_f in (0.0, 180.0):
            th_lab = str(th)
        anote = 'orig_angle=cm'
    n = anote + ('; ' + notes if notes else '')
    rows.append([dataset, texkey, obs, str(E), str(dE), th_lab, th_cm, frame_cs, str(value), str(stat),
                 str(sys_pt), str(norm), units, source, n])


def read_hepdata_csv(path):
    """Return list of blocks; each block: dict(qual=..., header=[...], rows=[[...]])"""
    blocks = []
    cur = None
    for line in open(path, encoding='utf-8', errors='replace'):
        line = line.rstrip('\n')
        if not line.strip():
            cur = None
            continue
        if line.startswith('#:'):
            if cur is None or cur.get('header'):
                cur = {'qual': {}, 'header': None, 'rows': []}
                blocks.append(cur)
            m = re.match(r'#: (.+),([^,]*)$', line)
            if m:
                cur['qual'][m.group(1).rstrip(',').strip()] = m.group(2).strip()
            continue
        if cur is None:
            cur = {'qual': {}, 'header': None, 'rows': []}
            blocks.append(cur)
        r = next(csv.reader([line]))
        if cur['header'] is None:
            cur['header'] = r
        else:
            cur['rows'].append(r)
    return [b for b in blocks if b['header'] and b['rows']]


def pct(s):
    s = s.strip()
    return float(s.rstrip('%')) / 100.0 if s.endswith('%') else None


def absfrac(val, frac, nd=3):
    return f"{float(val) * frac:.{nd}f}"


# ---------------------------------------------------------------- Olmos de Leon (TAPS)
def build_odl():
    path = os.path.join(RAW, 'OlmosdeLeon2000_thesis', 'transcription_weighted_mean.csv')
    pts = [l.strip().split(',') for l in open(path) if l[0].isdigit()]
    Es = sorted({float(p[1]) for p in pts})
    width = {}
    for i, E in enumerate(Es):
        if i == 0:
            w = Es[1] - Es[0]
        elif i == len(Es) - 1:
            w = Es[-1] - Es[-2]
        else:
            w = (Es[i + 1] - Es[i - 1]) / 2
        width[E] = f"{w:.1f}"
    for th, E, v, s in pts:
        sys_pt = f"{0.05 * float(v):.3f}"
        flag = '; flagged as outlier (excluded) by Griesshammer2012/McGovern2013/Pasquini2019-SELECTED' if (th, E) == ('133', '108.42') else ''
        add('OlmosdeLeon2001', 'deLeon:2001dnx', 'CS', E, th, 'lab', 'lab', v, s, sys_pt, '0.03',
            source='V. Olmos de Leon, PhD thesis Mainz 2000 (download.uni-mainz.de/fb08-kpha2/publications/phd/2000_olmos.ps.pdf), Tab. A.1-A.5 weighted mean',
            notes='stat as tabulated (incl. ~3% cut-variation component, thesis Sec.4.8.1); sys_err_pt = 5% point-to-point (geometry/simulation) per published EPJA 10,207 as quoted by Griesshammer2012/McGovern2013/Pasquini2019 (not in thesis table); norm 3% = flux 2% (+) target 2% (EPJA); thesis quotes 6% total systematic; values identical to author legacy a0xx.leon01/Olmos_*.dat; dE approx = bin-centre spacing (thesis: 20 tagger channels = 6-10 MeV); theta bins +-10 deg' + flag,
            dE=width[float(E)])


# ---------------------------------------------------------------- Hallin (SAL 1993)
def build_hallin():
    d = os.path.join(RAW, 'Hallin1993_HEPData', 'csv')
    for i in range(1, 17):
        for b in read_hepdata_csv(os.path.join(d, f'Table{i}.csv')):
            h = b['header']
            q = b['qual']
            tabI = i <= 4
            ds = 'Hallin1993_TabI' if tabI else 'Hallin1993_TabII'
            for r in b['rows']:
                if tabI:
                    E = q['E [MeV]']; th = r[0]
                else:
                    E = r[0]; th = q['THETA [deg]']
                v, st, sy = r[1], r[2], r[4]
                note = ('HEPData ins369337 Table %d; ' % i +
                        ('paper Table I (angular distributions, 15 MeV wide brems. bins); NOT used by Griesshammer2012/McGovern2013/Pasquini2019 (not clearly independent of Table II)' if tabI else
                         'paper Table II (excitation functions, 5 MeV bins)') +
                        '; angles are c.m. (HEPData THETA; abstract gives 25<theta_lab<135 which maps to 31-144 c.m.); cross-section frame assumed c.m. (consistent with author legacy a060.hal93 = these values x lab/cm Jacobian); sys_err_pt = per-point systematic as tabulated (paper: correlated flux 3% (+) target 2% = 3.6% plus 1-2% p2p, i.e. contains the norm part - do not add both); LEGS2001 Tab.IX used 3.8% scale')
                if (not tabI) and float(E) > 150:
                    note += '; McGovern2013/Pasquini2019 use only E<=150 MeV (cusp); Griesshammer2012 excluded all Hallin'
                if (not tabI) and q.get('THETA [deg]') == '141.0' and r[0] == '170.0':
                    note += '; flagged as outlier by Griesshammer2012'
                add(ds, 'Hallin:1993ft', 'CS', E, th, 'cm', 'cm', v, st, sy, '0.036',
                    source='HEPData ins369337 (https://www.hepdata.net/record/ins369337) Table %d' % i, notes=note,
                    dE=('15' if tabI else '5'))


# ---------------------------------------------------------------- MacGibbon (SAL 1995)
def build_macgibbon():
    tex = open(os.path.join(RAW, 'MacGibbon1995_arXiv', 'src', '9507001.tex')).read()
    tg = tex[tex.find('label{tab:modeldep}'):]
    t1 = tg[tg.find('Tagged differential'):tg.find('label{tab:cross_tagged}')]
    for line in t1.split('\n'):
        m = re.match(r'^\s*([\d.]+)\s*&\s*\$\s*([\d.]+)\s*\\pm\s*([\d.]+)\$\s*&\s*\$\s*([\d.]+)\s*\\pm\s*([\d.]+)\$', line)
        if m:
            E, v135, e135, v90, e90 = m.groups()
            for th, v, e in (('135', v135, e135), ('90', v90, e90)):
                add('MacGibbon1995_tagged', 'MacGibbon:1995in', 'CS', E, th, 'lab', 'lab', v, e, '', '0.04',
                    source='arXiv:nucl-ex/9507001 source, Table "Tagged differential cross sections" (= HEPData ins397101 T1)',
                    notes='stat errors; systematic (overall) 3.0-4.0% for tagged data (paper Sec. Systematic Errors), 0.04 = upper value, also used by LEGS2001 Tab.IX; tagged bins 8 MeV wide',
                    dE='8')
    t2 = tg[tg.find('Untagged differential'):tg.find('label{tab:cross_untagged}')]
    for line in t2.split('\n'):
        m = re.match(r'^\s*(\d+)--(\d+)\s*&\s*([\d.]+)\s*\$\\pm\$\s*([\d.]+)\s*&\s*([\d.]+)\s*\$\\pm\$\s*([\d.]+)', line)
        if m:
            lo, hi, v135, e135, v90, e90 = m.groups()
            E = f"{(int(lo) + int(hi)) / 2:.0f}"
            for th, v, e in (('135', v135, e135), ('90', v90, e90)):
                add('MacGibbon1995_untagged', 'MacGibbon:1995in', 'CS', E, th, 'lab', 'lab', v, e, '', '0.029',
                    source='arXiv:nucl-ex/9507001 source, Table "Untagged differential cross sections" (= HEPData ins397101 T2)',
                    notes=f'bin {lo}-{hi} MeV; error = diagonal element of 5x5 error matrix (adjacent bins strongly anti-correlated; full matrix not published); overall systematic 2.9% (untagged)',
                    dE=str(int(hi) - int(lo)))


# ---------------------------------------------------------------- Baranov 1974/75 (Moscow)
def build_baranov74():
    d = os.path.join(RAW, 'Baranov1975_HEPData', 'csv')
    for i, th in ((1, '90'), (2, '150')):
        for b in read_hepdata_csv(os.path.join(d, f'Table{i}.csv')):
            for r in b['rows']:
                E, v, e = r[0], r[1], r[2]
                add(f'Baranov1974_75_th{th}', 'Baranov:1974ta;Baranov:1974ec;Baranov:1975ju', 'CS', E, th, 'lab', 'lab',
                    scale(v, 1000), scale(e, 1000), '', '0.08',
                    source=f'HEPData ins102948 (Yad.Fiz.21,689) Table {i}; identical to JETP Lett.19,398 (1974) table (INSPIRE file for ins96876)',
                    notes='converted mub/sr x1000; errors are total (random + reduction + ~1% absolutization systematic; JETP Lett.), placed in stat_err; norm 0.08 = scale uncertainty used for "Moscow 81-110 MeV" in LEGS PRC64 Tab.IX (paper itself quotes ~1%); Baranov2001 compilation treats 90 and 150 deg as separate sets with different systematic errors; Griesshammer2012/McGovern2013 exclude the 150 deg set' + (' (this set)' if th == '150' else ''))


# ---------------------------------------------------------------- Federspiel (Illinois 1991)  [secondary: author legacy]
def build_federspiel():
    leg = open(os.path.join(RAW, 'author_legacy', 'author_agr_embedded_data.txt')).read().split('\n\n')
    for b in leg:
        h = b.split('\n')[0]
        for key, th in (('a060.fed91', '60'), ('a135.fed91', '135')):
            if key in h:
                for l in b.strip().split('\n')[1:]:
                    E, v, e = l.split()[:3]
                    note = ('SECONDARY SOURCE: author legacy file /home/vvp/F90/Compton/DATA2/%s (embedded in LP_ms2/figs/dasfixA.agr etc.); PRL 67,1511 shows cross sections only in Fig.2 (no table); these values are ~5%% below the published figure (consistent with later corrections by A.M. Nathan, priv. comm. 1993, quoted in Baranov et al. 2001 compilation) - verify before use; stat errors; norm 2.0%% (PRL preprint ILL-NPL-91-012), Griesshammer2012 quotes 2.2%%; 4 MeV bins' % key)
                    if th == '135' and E == '44.2':
                        note += '; flagged as outlier by Griesshammer2012/McGovern2013'
                    add('Federspiel1991', 'Federspiel:1991yd', 'CS', E, th, 'lab', 'lab', v, e, '', '0.020',
                        source='author legacy %s (cf. KEK scan of preprint ILL-NPL-91-012, Fig.2)' % key, notes=note, dE='4')


# ---------------------------------------------------------------- Zieger (MAMI 1992) [secondary: author legacy]
def build_zieger():
    for E, v, e, sp in (('98', '18.42', '3.39', 0.064), ('132', '21.82', '1.18', 0.043)):
        add('Zieger1992', 'Zieger:1992jq', 'CS', E, '180', 'lab', 'lab', v, e, f"{float(v) * sp:.3f}", '0',
            source='author legacy Codes/execs/CSeft/RCSdata/zieg_92.exp (= dasfixA.agr G5.S1); PLB 278,34 not accessible',
            notes=f'SECONDARY SOURCE; stat = {float(e) / float(v) * 100:.1f}% matches the 18.4%/5.4% quoted by Griesshammer2012; sys_err_pt = {sp * 100:.1f}% (Griesshammer2012/LEGS2001: no common scale, systematics combined per point); 16 MeV wide brems. bins; recoil proton at 0 deg; cross-section frame taken as lab following author legacy plot (rcs.ps, dsigma/dOmega^lab) - verify',
            dE='16')


# ---------------------------------------------------------------- old partial legacy sets (Hyman 1959, Pugh/Frisch 1967)
def build_old_legacy():
    leg = open(os.path.join(RAW, 'author_legacy', 'author_agr_embedded_data.txt')).read().split('\n\n')
    spec = {'a090.hym59': ('Hyman1959_legacy_partial', 'Hyman:1959zz', 'Hyman et al. PRL 3,93 (1959)'),
            'a090.fri67': ('Pugh1957_Frisch1967_legacy_partial', 'Pugh:1957zz', 'Pugh et al. PR 105,982 (1957) as corrected by D. Frisch 1967 (MIT summer study / priv. comm., cf. Baranov et al. 2001)')}
    for b in leg:
        h = b.split('\n')[0]
        for key, (ds, tk, ref) in spec.items():
            if key in h:
                for l in b.strip().split('\n')[1:]:
                    E, v, e = l.split()[:3]
                    add(ds, tk, 'CS', E, '90', 'lab', 'lab', v, e, '', '',
                        source='author legacy file /home/vvp/F90/Compton/DATA2/%s (embedded in Codes/execs/CSeft/FIGS/dsfixA.mgr, fApolza.mgr)' % key,
                        notes='SECONDARY, PARTIAL (only the 90 deg subset survives in the author files): %s; original numeric tables not publicly accessible; angle assumed lab 90 deg; error type and normalization unknown (Baranov2001 / Griesshammer2012 use these sets with their own systematic assignments)' % ref)


# ---------------------------------------------------------------- LEGS 2001
def build_legs():
    path = os.path.join(RAW, 'Blanpied2001_LEGSweb', 'transcription_TableII_III.csv')
    for l in open(path):
        if not l.startswith('II'):
            continue
        t, E, th, v, e = l.strip().split(',')
        if t == 'II':
            add('Blanpied2001_LEGS', 'Blanpied:2001ae', 'CS', E, th, 'cm', 'cm', v, e, '', '0.02',
                source='PRC 64,025203 Table II, from LEGS-posted manuscript legs_n_delta.pdf (web.archive.org/web/20050118065130/http://www.legs.bnl.gov:80/pub/pi_comp/legs_n_delta.pdf)',
                notes='error = stat (+) point-specific systematic (net measurement error); +2% common scale; Griesshammer2012 notes LEGS > MAMI(Wolf) by up to ~25% at theta_cm>100 deg near peak; tagger resolution 5.4 MeV FWHM')
        else:
            add('Blanpied2001_LEGS', 'Blanpied:2001ae', 'S3', E, th, 'cm', 'n/a', v, e, '', '0',
                units='dimensionless',
                source='PRC 64,025203 Table III, from LEGS-posted manuscript legs_n_delta.pdf (web.archive.org)',
                notes='linear beam asymmetry Sigma (LEGS convention, positive near 90 deg at the Delta peak); error = stat (+) point-specific systematic; no additional scale uncertainty (paper Sec. VI A)')


# ---------------------------------------------------------------- Huenger 1997 (MAMI; includes Peise 1996 75deg and Molinari 1996 90deg)
def build_hunger():
    d = os.path.join(RAW, 'Hunger1997_HEPData', 'csv')
    for i, th, prev in ((1, '75', 'Peise:1996zz'), (2, '90', 'Molinari:1996zw')):
        for b in read_hepdata_csv(os.path.join(d, f'Table{i}.csv')):
            for r in b['rows']:
                E, v, st, sy = r[0], r[1], r[2], r[4]
                f = pct(sy)
                add(f'Hunger1997_th{th}cm', f'Hunger:1997pt;{prev}', 'CS', E, th, 'cm', 'cm', v, st, '', f"{f:.3f}",
                    source=f'HEPData ins458618 (https://www.hepdata.net/record/ins458618) Table {i}',
                    notes=f'c.m. cross section at theta_cm={th} (HEPData RF=CM); stat errors; quoted systematic {sy} placed in norm_err_frac (contains common flux 2% and target 2%; ' +
                          ('energy-dependent: 4.4% <350, 6% 350-380, 8% >380 MeV' if th == '75' else '4% <360, 5% above (pi0 background)') +
                          f'; LEGS2001 used {"4.4" if th == "75" else "5.0"}% scale); NPA 620,385 summarises/supersedes {prev}')


# ---------------------------------------------------------------- Wissmann 1999
def build_wissmann():
    d = os.path.join(RAW, 'Wissmann1999_HEPData', 'csv')
    for b in read_hepdata_csv(os.path.join(d, 'Table1.csv')):
        for r in b['rows']:
            E, v, st = r[0], r[1], r[2]
            add('Wissmann1999', 'Wissmann:1999vi', 'CS', E, '130.7', 'lab', 'lab', v, st, '', '0.03',
                source='HEPData ins514218 (https://www.hepdata.net/record/ins514218) Table 1',
                notes='theta_lab=130.7 (HEPData RF=LAB; LEGS2001 cites "130.7 lab"); cross-section frame NOT stated in HEPData - taken as lab (magnitude ~170 nb/sr at peak matches lab-frame Kossert/Camen data at 136 deg lab; c.m. would be ~1.47x larger) - verify; stat errors; systematic 3% common (flux, target, acceptance)')


# ---------------------------------------------------------------- Camen 2002 / Kossert 2003 free proton
def build_camen():
    d = os.path.join(RAW, 'Kossert2003_HEPData', 'csv')
    for b in read_hepdata_csv(os.path.join(d, 'Table1.csv')):
        for r in b['rows']:
            E, v, st = r[0], r[1], r[2]
            add('Camen2002_Kossert2003', 'Camen:2001st;Kossert:2002ws', 'CS', E, '136.2', 'lab', 'lab', v, st, '', '0.044',
                source='Kossert et al. EPJA 16,259 (2003) appendix table (arXiv:nucl-ex/0210020 source; = HEPData ins599960 Table 1)',
                notes='free-proton data of the Camen et al. PRC 65,032202 experiment (Camen shows them only graphically, shifted to theta_cm=135); theta_lab=136.2; frame lab (Fig.7 of Kossert; Camen Fig.2 c.m. values are 1.2-1.7x larger, consistent with lab->cm Jacobian); stat errors; systematic 4.4%')


# ---------------------------------------------------------------- Kolb 2000 (SAL) free-proton calibration point
def build_kolb():
    import tarfile
    src = os.path.join(RAW, 'Kolb2000_arXiv', 'source_nucl-ex_0003002.tgz')
    with tarfile.open(src) as tf:
        tex = [tf.extractfile(m).read().decode('latin-1') for m in tf.getmembers() if m.name.endswith('paper.tex')][0]
    m = re.search(r'For E\$_\\gamma = (\d+)\$~MeV and \$\\theta_\\gamma=-(\d+)\^\\circ\$,\s*the cross section for elastic Compton scattering\s*from the free proton.*?=\s*([\d.]+)\s*\\pm\s*([\d.]+)\s*\\pm\s*([\d.]+)~\{\\rm nb/sr\}', tex, flags=re.S)
    E, th, v, st, sy = m.groups()
    add('Kolb2000', 'Kolb:2000ix', 'CS', E, th, 'lab', 'lab', v, st, sy, '',
        source='Kolb et al. PRL 85,1388 (2000), Eq.(2), arXiv:nucl-ex/0003002 source',
        notes='free-proton calibration measurement of the SAL quasi-free deuteron experiment (BUNI at -135 deg lab, tagged photons summed up to 260 MeV, quoted at 247 MeV); stat and systematic as given (systematic placed in sys_err_pt, overall for this single point); not in Griesshammer2012 table; frame lab (appears as "SAL-00" lab point in Kossert2003 Fig.7)')


# ---------------------------------------------------------------- old Delta-region data from HEPData
def build_old_delta():
    # Gray & Hanson 1967
    d = os.path.join(RAW, 'Gray1967_HEPData', 'csv')
    for b in read_hepdata_csv(os.path.join(d, 'Table3.csv')):
        for r in b['rows']:
            th, E, lo, hi, v, e = r[0], r[1], r[2], r[3], r[4], r[5]
            add('Gray1967', 'Gray:1967zzd', 'CS', E, th, 'cm', 'cm', scale(v, 1000), scale(e, 1000), '', '',
                source='HEPData ins944949 Table 3', dE=f"{float(hi) - float(lo):g}",
                notes=f'converted mub/sr x1000; bin {lo}-{hi} MeV; angle and cross-section frame not specified in HEPData - assumed c.m. (recoil-proton detection; angle varies with energy) - verify; errors as given (type unknown); Griesshammer2012: some points implausibly low')
    # Genzel 1972
    d = os.path.join(RAW, 'Genzel1972_HEPData', 'csv')
    for i in range(1, 7):
        for b in read_hepdata_csv(os.path.join(d, f'Table{i}.csv')):
            th = b['qual'].get('THETA(P=3,RF=CM) [deg]')
            for r in b['rows']:
                E, v, st, sy = r[0], r[1], r[2], r[4]
                f = pct(sy)
                add('Genzel1972', 'Genzel:1972ba;Genzel:1976wy', 'CS', E, th, 'cm', 'cm', scale(v, 1000), scale(st, 1000),
                    f"{float(v) * 1000 * f:.2f}", '',
                    source=f'HEPData ins77720 (Lett.Nuovo Cim.4,695) Table {i}',
                    notes=f'converted mub/sr x1000; theta c.m. (HEPData RF=CM), cross-section frame assumed c.m.; sys {sy} per point (as given) in sys_err_pt; Griesshammer2012 cite Genzel:1976wy (Z.Phys.A279,399, not on HEPData) with 19 points <=400 MeV = same count as this set; Bonn compilation data')
    # Baranov 1966 (Yad.Fiz. 3,791) tables 3 (215 MeV) and 4 (248 MeV) in mub/sr
    d = os.path.join(RAW, 'Baranov1966b_HEPData', 'csv')
    for i in (3, 4):
        for b in read_hepdata_csv(os.path.join(d, f'Table{i}.csv')):
            for r in b['rows']:
                th, tlo, thi, E, lo, hi, v, e, sy = r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[9]
                th = f"{float(th):g}"
                add('Baranov1966', 'Baranov:1966aaa', 'CS', E, th, 'cm', 'cm', scale(v, 1000), scale(e, 1000), '', f"{pct(sy):.2f}",
                    source=f'HEPData ins1392560 (Sov.J.Nucl.Phys.3,791) Table {i}', dE=f"{float(hi) - float(lo):.1f}",
                    notes=f'converted mub/sr x1000; bin {lo}-{hi} MeV; theta bin {tlo}-{thi}; angle assumed c.m. (Griesshammer2012 quotes 91.9 and 107.2 deg as c.m.), cross-section frame assumed c.m.; 20% photon-spectrum uncertainty common to all points (norm); ' +
                          ('Griesshammer2012: 215 MeV points at 91.9 and 107.2 deg (cm) implausibly low' if i == 3 else 'Table 4 incorporates revised data of Baranov 1962 (JETP 14,1219); overlaps with Baranov:1966oqa refinement - do not use both for the same angles'))
    d = os.path.join(RAW, 'Baranov1966a_HEPData', 'csv')
    for b in read_hepdata_csv(os.path.join(d, 'Table1.csv')):
        for r in b['rows']:
            th, tlo, thi, E, lo, hi, v, e = r[:8]
            th = f"{float(th):g}"
            E = f"{float(E):.1f}"
            ov = float(th) > 80
            add('Baranov1966_refined', 'Baranov:1966oqa', 'CS', E, th, 'cm', 'cm', scale(v, 1000), scale(e, 1000), '', '',
                source='HEPData ins1387582 (Sov.Phys.JETP 23,242) Table 1', dE=f"{float(hi) - float(lo):.1f}",
                notes=f'converted mub/sr x1000; bin {lo}-{hi} MeV; theta bin {tlo}-{thi}; re-analysis of Baranov 1962 (247 MeV); angle/frame assumed c.m.; systematic errors not quoted' +
                      ('; OVERLAPS Baranov1966 Table 4 at the same angle (do not use both)' if ov else '; forward points (with Baranov1966 T3+T4 they make the 12 points of Griesshammer2012 "Moscow 66")'))
    # Nagashima 1964
    d = os.path.join(RAW, 'Nagashima1964_HEPData', 'csv')
    for b in read_hepdata_csv(os.path.join(d, 'Table1.csv')):
        for r in b['rows']:
            E, th, v, e, sy = r[0], r[1], r[2], r[3], r[5]
            add('Nagashima1964', 'Nagashima:1964zz', 'CS', E, th, 'cm', 'cm', scale(v, 1000), scale(e, 1000), '', f"{pct(sy):.2f}",
                source='HEPData ins1101996 (INS report 81, 1964; Bonn compilation) Table 1',
                notes='converted mub/sr x1000; angle ~90 deg c.m. (assumed c.m.), cross-section frame assumed c.m.; 5% systematic as given')
    # Barbiellini 1968 asymmetry
    d = os.path.join(RAW, 'Barbiellini1968_HEPData', 'csv')
    for b in read_hepdata_csv(os.path.join(d, 'Table1.csv')):
        for r in b['rows']:
            th, v, e = r[0], r[1], r[2]
            add('Barbiellini1968', 'Barbiellini:1968rlx', 'S3', '317.5', th, 'cm', 'n/a', v, e, '', '', units='dimensionless',
                source='HEPData ins55887 (Phys.Rev.174,1665) Table 1', dE='35',
                notes='photon asymmetry for E=300-335 MeV at 90 deg c.m. (bin centre 317.5 entered here); SIGN CONVENTION differs from LEGS (LEGS Sigma(90cm,~318MeV) ~ +0.38): value copied as given, do not combine without checking the definition')


build_odl()
build_federspiel()
build_zieger()
build_baranov74()
build_old_legacy()
build_hallin()
build_macgibbon()
build_legs()
build_hunger()
build_wissmann()
build_camen()
build_kolb()
build_old_delta()

with open(OUT, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(HEADER)
    for r in rows:
        w.writerow(r)
print('wrote', len(rows), 'rows to', OUT)
