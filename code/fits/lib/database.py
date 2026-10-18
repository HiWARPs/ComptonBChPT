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
             KLP_E (energy reconstructed from the KLP bin label), ...
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

# --- normalisation uncertainties of the 1950s-60s bremsstrahlung sets (not stated in the sources) -----------
# Goldansky: ~8% (MacGibbon 1995, quoted in db_old_README); others: 10% adopted here (assumption,
# documented in the report; the GMPF/MPG/PPS papers do not tabulate them).
OLD_NORM = {'Oxley': 0.10, 'Hyman': 0.10, 'Goldansky': 0.08, 'Bernardini': 0.10, 'Pugh': 0.10}

# --- KLP2018 bin labels -> true energies (Pugh/Hyman), from the 90-deg legacy subsets in db_old
#     (Pugh: 55,75,95,115,125,135;  Hyman: 60,70,78,87,95,106,118,128).  Labels without a 90-deg partner
#     (Pugh 69 -> 65, 109 -> 105) are inferred from the 10-MeV bin pattern of the experiment.
KLP_E = {'Pugh': {59: 55, 69: 65, 79: 75, 99: 95, 109: 105, 117: 115, 127: 125, 135: 135},
         'Hyman': {60: 60, 70: 70, 78: 78, 89: 87, 99: 95, 109: 106, 119: 118, 128: 128}}

EXPT = {
    'OlmosdeLeon2001': 'TAPS 01 (Olmos de León)', 'Federspiel1991': 'Illinois 91 (Federspiel)',
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
        elif ds == 'Federspiel1991':
            rr['prov'] = 'stat; norm 2.0% (legacy values, see db_old_README)'
            if abs(rr['E'] - 44.2) < 0.01 and abs(rr['theta_lab'] - 135) < 0.1:
                rr['flags'].append('GMPF_out')
        elif ds == 'Zieger1992':
            rr['prov'] = 'stat (+) syst per point (6.4%, 4.3%); no common scale (not floated)'
            rr['norm'] = 0.0
        elif ds.startswith('Baranov1974'):
            rr['prov'] = 'error as quoted; norm 8% (LEGS 2001 Tab. IX)'
            if 'th150' in ds:
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
            if abs(rr['E'] - 170) < 0.01 and abs(rr['theta_cm'] - 141) < 0.1:
                rr['flags'].append('GMPF_outlier170')
        elif ds.startswith('MacGibbon'):
            rr['prov'] = 'stat; norm %.1f%%' % (100 * rr['norm'])
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


def _rows_klp():
    d = pd.read_csv(os.path.join(DATA, 'db_klp2018.csv'))
    out = []
    for i, r in d.iterrows():
        nm = r['dataset'].replace('KLP2018_', '')
        if nm not in OLD_NORM:
            continue     # all other KLP sets are taken from db_old with true energies
        E = float(r['E_lab_MeV']); flags = []
        if nm in KLP_E:
            E = float(KLP_E[nm][int(round(E))]); flags.append('KLP_E')
        thl = float(r['theta_lab_deg'])
        rr = dict(source='db_klp2018', obs='CS', dobs='CS', frame='LAB', E=E, theta_lab=thl,
                  theta_cm=float(th_lab2cm(E, thl)), value=float(r['value']), err_stat=float(r['stat_err']),
                  err_syspt=0.0, norm=OLD_NORM[nm], flags=flags,
                  prov='KLP2018 total error; norm %d%% adopted%s' % (round(100 * OLD_NORM[nm]),
                        '; E from legacy 90-deg subset/bin pattern' if nm in KLP_E else '; E = KLP value'))
        rr['dataset'] = {'Oxley': 'Oxley1958', 'Hyman': 'Hyman1959', 'Goldansky': 'Goldansky1960',
                         'Bernardini': 'Bernardini1960', 'Pugh': 'Pugh1957'}[nm]
        if nm in ('Oxley', 'Bernardini'):
            rr['flags'].append('GMPF_out')
        out.append(rr)
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
    rows = _rows_old() + _rows_klp() + _rows_new(per_beamtime) + _rows_forward()
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
