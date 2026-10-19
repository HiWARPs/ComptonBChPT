"""common.py -- shared helpers for the analysis scripts (selections, driver outputs, fit construction)."""
import os, sys, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
FITS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import bchpt, database, forward, theory, fitter  # noqa: E402
from fitter import Chi2Fit, DIR, C_BALDIN, C_GAMMA0, gamma0, gammapi  # noqa: E402,F401

DBDIR = os.path.join(FITS, 'db')
RES = os.path.join(FITS, os.environ.get('FIT_RESULTS', 'results'))
FIGS = os.path.join(FITS, os.environ.get('FIT_FIGS', 'figs'))
PROFILE = bchpt.PROFILE
STYLE = os.path.join(os.path.dirname(FITS), 'figstyle.mplstyle')

OLD_L = ['OlmosdeLeon2001', 'Federspiel1991', 'Zieger1992', 'Baranov1974_75_th90', 'Baranov1974_75_th150',
         'Hallin1993_TabII', 'MacGibbon1995_tagged', 'MacGibbon1995_untagged', 'Pugh1957', 'Hyman1959',
         'Goldansky1960', 'Oxley1958', 'Bernardini1960']
A2 = ['A2_2022_CS', 'A2_2022_S3']
A2PB = ['A2_2018Mar_CS', 'A2_2018Jul_CS', 'A2_2018Mar_S3', 'A2_2018Jul_S3']
HIGS = ['HIGS_2022_CIRC', 'HIGS_2022_LIN']
FWD = ['GHP_fwd_CS', 'GHP_fwd_S2Z']
# MAIN DATABASE (angular sets) = pre-2015 world data + A2 2022.  The HIgS 2022 sets (Li et al.) are EXCLUDED from the
# main database since 2026-10-17 (author decision): they are kept as a test (HIgS alone, HIgS predicted from the main
# fit, parameter-goodness-of-fit) and as the database variant 'j_withHIGS' (= the former main database) in
# scripts/s03_fits_L.py.  MAIN_L_SETS does not include the forward points (FWD) nor the per-beamtime A2 sets (A2PB).
MAIN_L_SETS = OLD_L + A2
MAIN_L_SETS_WITH_HIGS = OLD_L + A2 + HIGS
D_SETS = ['A2_Martel2015_S2x', 'A2_Paudyal2020_S2z', 'LEGS2001_CS', 'LEGS2001_S3', 'Hunger1997_75', 'Hunger1997_90',
          'Wissmann1999', 'Camen2002', 'Hallin1993_TabII', 'Hallin1993_TabII_D', 'OlmosdeLeon2001',
          'Gray1967', 'Genzel1972', 'Baranov1966', 'Kolb2000', 'Nagashima1964']

PREF_NNLO_V1 = np.array([11.246404, 3.870657, -3.335106, 2.942725, 0.178734, 1.146936])
# NNLO static polarisabilities of the base amplitude of the active profile (v2: incl. the physical-mass shift)
PREF_NNLO = bchpt.run([('pref', 100.0, 90.0, 'CS', 'LAB')], 'NNLO', name='pref', verbose=False).Pref.copy()
PAUDYAL = np.array([-2.87, 2.70, -0.85, 2.04])
# Sum-rule constraints (Baldin alpha+beta [1e-4 fm^3]; forward spin polarisability gamma_0 [1e-4 fm^4]); source switch
# FWD_SOURCE (lib/forward.py, environment).  '2026' (default since 2026-10-18): ForwardRCS 2026 update
# (work/data/raw/ForwardRCS2026/unpol_sumrules.json: alpha+beta = 14.0092 +- 0.1814 (exp+fit) +- 0.1022 (data selection);
# pol_sumrules.json 'unconstrained': gamma_0 = -102.479 +- 8.013 (exp+fit) +- 6.248 (PWA model) 1e-6 fm^4), errors added in
# quadrature.  'ghp2016': GHP15 average 14.0(2) and GHP16 gamma_0 = -0.929(105).  (The name GHP_G0 is kept for the consumers.)
if forward.SOURCE == '2026':
    BALDIN = (14.0092, 0.2082)
    GHP_G0 = (-1.0248, 0.1016)
else:
    BALDIN = (14.0, 0.2)
    GHP_G0 = (-0.929, 0.105)
BAL_TXT = '%.1f(%d)' % (BALDIN[0], round(10 * BALDIN[1])) if forward.SOURCE == 'ghp2016' else '%.2f(%d)' % (BALDIN[0], round(100 * BALDIN[1]))
G0_TXT = '%.3f(%d)' % (GHP_G0[0], round(1000 * GHP_G0[1]))
# MAIN DATABASE = MAIN_L_SETS (pre-2015 + A2 2022; HIgS 2022 excluded since 2026-10-17, see above) minus MAIN_DROP.
# MAIN_DROP holds the pid(s) flagged by the 3-sigma rule
# (|pull| > 3 w.r.t. the main-configuration fit F2 + truncation covariance, two passes) applied ONCE to the
# all-points database; it was determined on 2026-10-14 and is held fixed here (s03_fits_L.py verifies it at every
# run and stops if the all-points search gives a different list).  The single point is TAPS (Olmos de Leon 2001)
# at E = 108.42 MeV, theta_lab = 133 deg (pull -4.27).  Every selection made through s03_fits_L.rows_for() drops
# these pids unless keep_outlier=True (variant 'd_allpts'); scripts that select rows directly must drop them
# themselves (see scenarios/nnlo_prediction_chi2.py).
MAIN_DROP = ('OlmosdeLeon2001_029',)
MG_COV = {(int(r.theta_lab_deg), float(r.E_i_MeV), float(r.E_j_MeV)): float(r.cov_nb2) for r in
          pd.read_csv(os.path.join(database.DATA, 'db_macgibbon1995_untagged_cov.csv')).itertuples()}


def main_fit():
    """the main BChPT fit (since 2026-10-19): combined fit of alpha_E1, beta_M1 and R_EM to the data up to 360 MeV with
    the Delta region weighted by 0.01 (scenarios/combined_fit.py, results/combined_fit.json['combined']); the spin
    polarisabilities at NNLO.  cfg = the driver configuration at the fitted g_E (for curves)."""
    J = json.load(open(os.path.join(RES, 'combined_fit.json')))['combined']
    return dict(P=np.array(J['P']), covP=np.array(J['covP']), gE=float(J['gE'][0]), cfg='gE%.4f' % J['gE'][0],
                REM=tuple(J['REM']))


def load_db(per_beamtime=False):
    fn = os.path.join(DBDIR, 'fitdb_pb.csv' if per_beamtime else 'fitdb.csv')
    df = pd.read_csv(fn, keep_default_na=False)
    for c in ['E', 'theta', 'theta_lab', 'theta_cm', 'value', 'err_stat', 'err_syspt', 'err_pt', 'norm']:
        df[c] = df[c].astype(float)
    return df


_OUT = {}


def outputs(df, cfg='NNLO', tag='db'):
    """DriverOutput for all rows of df (cached on disk and in memory)."""
    key = (tag, cfg if isinstance(cfg, str) else json.dumps(cfg, sort_keys=True), len(df))
    if key not in _OUT:
        pts = database.points_for_driver(df)
        _OUT[key] = bchpt.run(pts, cfg, name='%s_%s' % (tag, cfg if isinstance(cfg, str) else 'custom'))
    return _OUT[key]


def base_static(T, cfg):
    """static polarisabilities contained in the base amplitude of configuration cfg.
    v2: the driver's Pref follows the mechanisms switched on (incl. the physical-mass shift) -> Pref."""
    if PROFILE == 'v2':
        return T.Pref.copy()
    c = dict(bchpt.BASE); c.update(bchpt.CONFIGS[cfg] if isinstance(cfg, str) else cfg)
    P = np.zeros(6)
    if c['oneNloop']:
        P += T.comp['piN loops']
    if c['delta_pole']:
        P += T.comp['Delta pole']
    if c['oneDloop']:
        P += T.comp['piDelta loops']
    return P


def select(df, sets, Emax=150.0, Emin=0.0, exclude_flags=(), extra_mask=None):
    m = df['dataset'].isin(sets) & (df['E'] <= Emax + 1e-6) & (df['E'] >= Emin)
    for f in exclude_flags:
        m &= ~df['flags'].str.contains(f, regex=False)
    if extra_mask is not None:
        m &= extra_mask
    return df[m].copy()


def build_fit(df, T, P0, free, constraints=(), fwd_mode='nuisance', thcov=None, theory_fn=None,
              extra=None, extra_prior=None, Toffset=None, append=None):
    """df: selected rows; T: DriverOutput covering (at least) these pids.
    Toffset: shift so that P (total static pols) maps to the driver's P = Pref + d  (P_driver = P - Toffset).
    append: optional extra data block (y, C, label, theory(P, extra)) not tied to rows of df (e.g. the forward
    amplitudes of lib/forward_amp.py), appended to the data vector with its own covariance, no normalisation."""
    idx = T.subset(df['pid'].values)
    y = df['value'].values
    N = len(df)
    C = np.diag(df['err_pt'].values ** 2)
    # MacGibbon 1995 untagged: full 5x5 error matrix per angle (MacGibbon et al., appendix; Baranov et al. 2001, Sect. 4)
    im = np.where((df['dataset'] == 'MacGibbon1995_untagged').values)[0]
    for a in im:
        for b in im:
            if round(df['theta_lab'].iat[a]) == round(df['theta_lab'].iat[b]):
                C[a, b] = MG_COV[(round(df['theta_lab'].iat[a]), df['E'].iat[a], df['E'].iat[b])]
    shifts = None
    isf = df['dataset'].isin(FWD).values
    if isf.any():
        cf, sh, fidx = forward.model(df, mode=fwd_mode)
        C[np.ix_(isf, isf)] = cf[np.ix_(isf, isf)]
        shifts = sh
    if thcov is not None:
        C = C + thcov
    norms = {k: float(df[df.dataset == k]['norm'].iloc[0]) for k in df['dataset'].unique()}
    off = np.zeros(6) if Toffset is None else Toffset
    if theory_fn is None:
        def theory_fn(P, ex):
            return T.predict_array(P - off, idx)
    ds = df['dataset'].values
    if append is not None:
        yb, Cb, lab, tb = append
        nb = len(yb)
        y = np.r_[y, yb]
        C = np.block([[C, np.zeros((N, nb))], [np.zeros((nb, N)), Cb]])
        ds = np.r_[ds, [lab] * nb]
        norms[lab] = 0.0
        if shifts is not None:
            shifts = np.hstack([shifts, np.zeros((len(shifts), nb))])
        tf0 = theory_fn
        theory_fn = lambda P, ex: np.r_[tf0(P, ex), tb(P, ex)]
    return Chi2Fit(y, theory_fn, C, ds, norms, P0, free, constraints=constraints,
                   shifts=shifts, extra=extra, extra_prior=extra_prior)


def fmt(v, e, nd=2):
    return ('%.' + str(nd) + 'f(%.' + str(nd) + 'f)') % (v, e)


def save_json(obj, name):
    os.makedirs(RES, exist_ok=True)

    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        raise TypeError(type(o))
    json.dump(obj, open(os.path.join(RES, name), 'w'), indent=1, default=conv)


def order_values(dfx, s, tag='db', scheme='B'):
    """order-by-order values (y0, y2, y3) for the rows s of database dfx (see theory.py for the schemes)."""
    T3 = outputs(dfx, 'NNLO', tag=tag)
    idx = T3.subset(s['pid'].values)
    y3 = T3.value_file[idx]
    if scheme == 'A':
        return outputs(dfx, 'LO', tag=tag).value_file[idx], outputs(dfx, 'NLO', tag=tag).value_file[idx], y3

    def matched(cfg):
        TT = outputs(dfx, cfg, tag=tag)
        off = base_static(TT, cfg) - TT.Pref
        return TT.predict_array(T3.Pref - off, idx)
    return matched('BornPi0'), matched('NLO'), y3


def truncation_for(dfx, s, tag='db', scheme='B', mode='sqrtQ', Lam=750.0, cbar='rms', corr=True, lE=25.0, lth=40.0):
    y0, y2, y3 = order_values(dfx, s, tag, scheme)
    tr = theory.truncation(y0, y2, y3, s['dobs'].values, s['E'].values, mode=mode, Lam=Lam, cbar=cbar)
    if corr:
        C = theory.covariance(tr['sigma'], s['dobs'].values, s['E'].values, s['theta_cm'].values, lE=lE, lth=lth)
    else:
        C = np.diag(tr['sigma'] ** 2)
    tr['y0'], tr['y2'], tr['y3'] = y0, y2, y3
    return C, tr
