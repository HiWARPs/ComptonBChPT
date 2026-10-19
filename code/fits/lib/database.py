"""
database.py -- assemble the fit database from the curated CSVs in work/data (never edited).

Unified columns of the fit database (fitdb.csv):
  pid        unique point id
  dataset    independent data set (own normalisation n_k)
  expt       short label for legends
  source     db file of origin
  obs        observable as in the source (CS, S3, S2X, S2Z, CS_PAR, CS_PERP)
  dobs       driver observable code (CS, CSCM, S3, S2X, S2Z, SPAR, SPERP)
  frame      angle frame passed to the driver (LAB or CM)
  E          lab photon energy (MeV, true energy where known)
  theta      angle in `frame` (deg);  theta_lab, theta_cm both given
  value      measured value (nb/sr for cross sections; fraction for asymmetries)
  err_stat, err_syspt, err_pt (= stat (+) uncorrelated syst used in chi^2)
  norm       1-sigma (Gaussian-equivalent) scale uncertainty of the set (fraction); 0 = not floated
  prov       provenance / error-model remarks
  flags      ';'-separated: GMPF_out (outlier/excluded by Griesshammer/McGovern pruning),
             B01 (pre-2001 point from Baranov et al. 2001, Table 2), ...
"""
import os, math
import numpy as np
import pandas as pd
from forward import SOURCE as FWD_SOURCE, SUFFIX as FWD_SUFFIX

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(WORK, 'data')
MP = 938.272
ETHR = 144.69   # pi0 threshold, lab
E2PI = 320.7    # 2 pi(+-) threshold, lab (approx: (2 m_pi+M)^2 - M^2)/2M with m_pi+=139.57)

# --- pre-2001 data below 150 MeV (since 2026-10-19): Table 2 of Baranov, L'vov, Petrun'kin, Shtarkov, Fiz. Elem. Chast.
#     Atom. Yadra 32 (2001) 699 (db_baranov2001.csv): lab cross sections, statistical errors, and the systematic
#     (normalisation) uncertainty delta of each set, incl. the corrections of the original data quoted there.  These
#     rows REPLACE the db_old (and the former db_klp2018) rows of the same sets below 150 MeV (the range of Table 2);
#     the db_old rows from 150 MeV on (Hallin Table II: 150 MeV point, the Delta region) are kept.  The MacGibbon 1995 untagged points carry their full error matrix
#     (db_macgibbon1995_untagged_cov.csv; common.build_fit).
B01_SETS = ('Oxley1958', 'Hyman1959', 'Goldansky1960', 'Bernardini1960', 'Pugh1957', 'Baranov1974_75_th90',
            'Baranov1974_75_th150', 'Federspiel1991', 'Zieger1992', 'Hallin1993_TabII', 'MacGibbon1995_tagged',
            'MacGibbon1995_untagged')

EXPT = {
    'OlmosdeLeon2001': 'A2 01 (Olmos de León)', 'Federspiel1991': 'Illinois 91 (Federspiel)',
    'Zieger1992': 'MAMI 92 (Zieger)', 'Baranov1974_75_th90': 'Moscow 74 90° (Baranov)',
    'Baranov1974_75_th150': 'Moscow 74 150° (Baranov)', 'Hallin1993_TabII': 'SAL 93 (Hallin)',
    'MacGibbon1995_tagged': 'SAL 95 tagged (MacGibbon)', 'MacGibbon1995_untagged': 'SAL 95 untagged (MacGibbon)',
    'Pugh1957': 'MIT 57/67 (Pugh)', 'Hyman1959': 'MIT 59 (Hyman)', 'Goldansky1960': 'Moscow 60 (Goldansky)',
    'Oxley1958': 'Chicago 58 (Oxley)', 'Bernardini1960': 'Illinois 60 (Bernardini)',
    'A2_2022_CS': 'A2 22 dσ (Mornacchi)', 'A2_2022_S3': 'A2 22 Σ3 (Mornacchi)',
    'A2_2018Mar_CS': 'A2 Mar18 dσ', 'A2_2018Jul_CS': 'A2 Jul18 dσ', 'A2_2018Mar_S3': 'A2 Mar18 Σ3', 'A2_2018Jul_S3': 'A2 Jul18 Σ3',
    'HIGS_2022_CIRC': 'HIγS 22 dσ (Li)', 'HIGS_2022_LIN': 'HIγS 22 dσ∥,⊥ (Li)',
    'GHP_fwd_CS': 'forward dσ(0°)', 'GHP_fwd_S2Z': 'forward Σ2z(0°)',   # set keys kept (GHP_*); source: FWD_SOURCE
    'A2_Martel2015_S2x': 'A2 15 Σ2x (Martel)', 'A2_Paudyal2020_S2z': 'A2 20 Σ2z (Paudyal)',
    'LEGS2001_CS': 'LEGS 01 dσ', 'LEGS2001_S3': 'LEGS 01 Σ3', 'Hunger1997_75': 'MAMI 97 75°cm (Hünger)',
    'Hunger1997_90': 'MAMI 97 90°cm (Hünger)', 'Wissmann1999': 'MAMI 99 (Wissmann)', 'Camen2002': 'MAMI 02 (Camen/Kossert)',
    'Hallin1993_TabII_D': 'SAL 93 (Hallin) >150', 'Gray1967': 'Illinois 67 (Gray)', 'Genzel1972': 'Bonn 72 (Genzel)',
    'Baranov1966': 'Moscow 66 (Baranov)', 'Kolb2000': 'SAL 00 (Kolb)', 'Nagashima1964': 'Tokyo 64 (Nagashima)',
}


def th_lab2cm(E, thl):
    b = E / (E + MP)
    c = np.cos(np.radians(thl))
    return np.degrees(np.arccos(np.clip((c - b) / (1 - b * c), -1, 1)))


def th_cm2lab(E, thc):
    b = E / (E + MP)
    c = np.cos(np.radians(thc))
    return np.degrees(np.arccos(np.clip((c + b) / (1 + b * c), -1, 1)))


def _f(x):
    try:
        v = float(x)
        return 0.0 if math.isnan(v) else v
    except Exception:
        return 0.0


def _rows_old():
    d = pd.read_csv(os.path.join(DATA, 'db_old.csv'))
    out = []
    for i, r in d.iterrows():
        ds = r['dataset']
        if ds in ('Pugh1957_Frisch1967_legacy_partial', 'Hyman1959_legacy_partial', 'Hallin1993_TabI',
                  'Baranov1966_refined', 'Barbiellini1968'):
            continue      # duplicates of fuller sets / not independent / sign-convention issue
        if ds in B01_SETS and float(r['E_lab_MeV']) < 150.0:
            continue      # replaced by Baranov et al. 2001, Table 2 (omega < 150 MeV; _rows_b01)
        rr = dict(source='db_old', obs=r['obs'], E=float(r['E_lab_MeV']), theta_lab=float(r['theta_lab_deg']),
                  theta_cm=float(r['theta_cm_deg']), value=float(r['value']), err_stat=_f(r['stat_err']),
                  err_syspt=_f(r['sys_err_pt']), norm=_f(r['norm_err_frac']), prov='', flags=[])
        cmframe = (r['frame_of_cs'] == 'cm')
        orig_cm = str(r['notes']).startswith('orig_angle=cm')
        if r['obs'] == 'CS':
            rr['dobs'] = 'CSCM' if cmframe else 'CS'
        else:
            rr['dobs'] = r['obs']
        rr['frame'] = 'CM' if (cmframe or orig_cm) else 'LAB'
        name = ds
        if ds == 'OlmosdeLeon2001':
            rr['prov'] = 'stat (+) 5% p2p (GMPF/MPG/PPS); norm 3%'
            if abs(rr['E'] - 108.42) < 0.01 and abs(rr['theta_lab'] - 133) < 0.1:
                rr['flags'].append('GMPF_out')
        elif ds == 'Hallin1993_TabII':
            # sys column contains the 3.6% common scale -> keep only the p2p remainder
            s2 = rr['err_syspt'] ** 2 - (0.036 * rr['value']) ** 2
            rr['err_syspt'] = math.sqrt(s2) if s2 > 0 else 0.0
            rr['norm'] = 0.036
            rr['prov'] = 'c.m. angle & c.m. dsigma; stat (+) [sys (-) 3.6% scale]; norm 3.6%'
            if rr['E'] > 150.0:
                rr['flags'].append('GMPF_out')            # MPG 2013: Hallin used only <=150 MeV
                name = 'Hallin1993_TabII_D' if rr['E'] > 170 else 'Hallin1993_TabII'
            if rr['E'] <= 150.0:
                rr['norm'] = 0.04       # the 150-MeV point joins the B01 rows below 150 MeV (delta = 4%); the set's norm in
                                        # the fits is that of its first row (a B01 row); rows above 150 MeV (Delta region) keep 3.6%
            if abs(rr['E'] - 170) < 0.01 and abs(rr['theta_cm'] - 141) < 0.1:
                rr['flags'].append('GMPF_outlier170')
        elif ds == 'Blanpied2001_LEGS':
            name = 'LEGS2001_CS' if r['obs'] == 'CS' else 'LEGS2001_S3'
            rr['prov'] = 'stat (+) pt syst; ' + ('norm 2%' if r['obs'] == 'CS' else 'no scale')
        elif ds == 'Hunger1997_th75cm':
            name = 'Hunger1997_75'; rr['prov'] = 'stat; norm = quoted syst (4.4%/6%/8%)'
        elif ds == 'Hunger1997_th90cm':
            name = 'Hunger1997_90'; rr['prov'] = 'stat; norm = quoted syst (4%/5%)'
        elif ds == 'Wissmann1999':
            rr['prov'] = 'stat; norm 3%; frame lab (inferred)'
        elif ds == 'Camen2002_Kossert2003':
            name = 'Camen2002'; rr['prov'] = 'stat; norm 4.4%; lab (inferred)'
        elif ds == 'Kolb2000':
            rr['prov'] = 'stat (+) syst (3.8) per point; lab frame; not floated'
        elif ds in ('Gray1967', 'Genzel1972', 'Baranov1966', 'Nagashima1964'):
            rr['prov'] = 'as quoted; c.m. frame assumed (not stated in source); norm 5% (20% Baranov66) adopted'
            if rr['norm'] == 0:
                rr['norm'] = 0.05 if ds != 'Baranov1966' else 0.20
        rr['dataset'] = name
        out.append(rr)
    return out


def _rows_b01():
    d = pd.read_csv(os.path.join(DATA, 'db_baranov2001.csv'))
    out = []
    for r in d[d.used == 1].itertuples():
        E, thl = float(r.E_lab_MeV), float(r.theta_lab_deg)
        flags = ['B01']
        if r.dataset in ('Oxley1958', 'Bernardini1960', 'Baranov1974_75_th150') or \
                (r.dataset == 'Federspiel1991' and abs(E - 44.2) < 0.01 and abs(thl - 135) < 0.1):
            flags.append('GMPF_out')          # excluded by the Griesshammer/McGovern pruning (variant 'f')
        out.append(dict(source='db_baranov2001', obs='CS', dobs='CS', frame='LAB', E=E, theta_lab=thl,
                        theta_cm=float(th_lab2cm(E, thl)), value=float(r.value), err_stat=float(r.stat_err),
                        err_syspt=0.0, norm=float(r.norm_err_frac), flags=flags, dataset=r.dataset,
                        prov='Baranov et al. 2001 Tab. 2: lab, stat; norm %.1f%%%s' % (
                            100 * r.norm_err_frac, '; full error matrix' if r.dataset == 'MacGibbon1995_untagged' else '')))
    return out


def _rows_new(per_beamtime=False):
    d = pd.read_csv(os.path.join(DATA, 'db_new.csv'))
    out = []
    for i, r in d.iterrows():
        ds = r['dataset']
        rr = dict(source='db_new', obs=r['obs'], E=float(r['E_lab_MeV']), theta_lab=float(r['theta_lab_deg']),
                  theta_cm=float(r['theta_cm_deg']), value=float(r['value']), err_stat=_f(r['stat_err']),
                  err_syspt=_f(r['sys_err_pt']), norm=_f(r['norm_err_frac']), frame='LAB', flags=[], prov='')
        if ds.startswith('A2_Mornacchi2022'):
            if per_beamtime:
                continue
            if r['obs'] == 'CS':
                rr.update(dataset='A2_2022_CS', dobs='CS', prov='stat (+) p2p syst; norm 3%')
            else:
                rr.update(dataset='A2_2022_S3', dobs='S3', norm=0.05 / math.sqrt(3),
                          prov='stat (+) p2p syst; scale 5% uniform -> rms 2.89%')
        elif ds == 'HIGS_Li2022_CIRC':
            rr.update(dataset='HIGS_2022_CIRC', dobs='CS', prov='stat (+) p2p syst; norm 2.24%')
        elif ds == 'HIGS_Li2022_LIN':
            rr.update(dataset='HIGS_2022_LIN', dobs={'CS_PAR': 'SPAR', 'CS_PERP': 'SPERP'}[r['obs']],
                      prov='stat (+) p2p syst; norm 2.24% (independent of circular)')
        elif ds == 'A2_Martel2015_S2x':
            rr.update(dataset=ds, dobs='S2X', norm=0.0,
                      prov='error = stat (+) beam/target-pol syst (2.1% scale already inside); not floated')
        elif ds == 'A2_Paudyal2020_S2z':
            # Syst. column = correlated: 10.4% scale (+) carbon 3-6% -> float the 10.4%, keep remainder p2p
            s2 = rr['err_syspt'] ** 2 - (0.104 * rr['value']) ** 2
            rr.update(dataset=ds, dobs='S2Z', err_syspt=math.sqrt(s2) if s2 > 0 else 0.0,
                      prov='Rand. (+) [Syst. (-) 10.4% scale]; scale 10.4% floated')
        else:
            raise ValueError(ds)
        out.append(rr)
    if per_beamtime:
        t = pd.read_csv(os.path.join(DATA, 'db_mornacchi_thesis.csv'))
        for i, r in t.iterrows():
            bt = 'Mar' if 'Mar2018' in r['dataset'] else 'Jul'
            isS3 = (r['obs'] == 'S3')
            rr = dict(source='db_mornacchi_thesis', obs=r['obs'], dobs=r['obs'], frame='LAB', E=float(r['E_lab_MeV']),
                      theta_lab=float(r['theta_lab_deg']), theta_cm=float(r['theta_cm_deg']), value=float(r['value']),
                      err_stat=_f(r['stat_err']), err_syspt=_f(r['sys_err_pt']),
                      norm=(0.05 / math.sqrt(3) if isS3 else 0.03), flags=[],
                      dataset='A2_2018%s_%s' % (bt, 'S3' if isS3 else 'CS'),
                      prov='thesis per-beamtime; stat (+) p2p; ' + ('scale 5%/sqrt3' if isS3 else 'norm 3%'))
            out.append(rr)
    return out


def _rows_forward():
    d = pd.read_csv(os.path.join(DATA, 'db_forward_obs%s.csv' % FWD_SUFFIX))
    prov = ('forward amplitudes, ForwardRCS 2026 update' if FWD_SOURCE == '2026' else 'GHP 2016 forward amplitudes') + \
        '; strongly correlated errors (nuisance model)'
    out = []
    for i, r in d.iterrows():
        isS = (r['obs'] == 'S2Z')
        out.append(dict(source='db_forward_obs', obs=r['obs'], dobs=r['obs'], frame='LAB', E=float(r['E_lab_MeV']),
                        theta_lab=0.0, theta_cm=0.0, value=float(r['value']), err_stat=float(r['stat_err']),
                        err_syspt=0.0, norm=0.0, flags=['FWD'],
                        dataset='GHP_fwd_S2Z' if isS else 'GHP_fwd_CS',
                        prov=prov))
    return out


def build(per_beamtime=False):
    rows = _rows_old() + _rows_b01() + _rows_new(per_beamtime) + _rows_forward()
    df = pd.DataFrame(rows)
    df['theta'] = np.where(df['frame'] == 'CM', df['theta_cm'], df['theta_lab'])
    df['err_pt'] = np.sqrt(df['err_stat'] ** 2 + df['err_syspt'] ** 2)
    df['expt'] = [EXPT.get(x, x) for x in df['dataset']]
    df['flags'] = [';'.join(f) for f in df['flags']]
    # deterministic ids
    df = df.sort_values(['dataset', 'dobs', 'E', 'theta'], kind='mergesort').reset_index(drop=True)
    cnt = {}
    pids = []
    for ds in df['dataset']:
        cnt[ds] = cnt.get(ds, 0) + 1
        pids.append('%s_%03d' % (ds, cnt[ds]))
    df.insert(0, 'pid', pids)
    cols = ['pid', 'dataset', 'expt', 'source', 'obs', 'dobs', 'frame', 'E', 'theta', 'theta_lab', 'theta_cm',
            'value', 'err_stat', 'err_syspt', 'err_pt', 'norm', 'flags', 'prov']
    return df[cols]


def points_for_driver(df):
    return [(r.pid, r.E, r.theta, r.dobs, r.frame) for r in df.itertuples()]
