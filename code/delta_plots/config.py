"""
config.py -- parameters of the Delta-region figures (work/delta_plots).

Everything that defines a curve or a band is set here, so that the captions in
work/notes/delta_plots.md can be checked against one file.

Environment (all optional; one command regenerates everything: work/delta_plots/make_all.sh):
  AMPL_EXE        executable used for ALL curves and bands (default: work/code_v2/execs/ampl.exe)
  AMPL_TEMPLATE   optional in.dat template (placeholders as in runner.DEFAULT_TEMPLATE)
  DP_V2OPTS       entries of the &v2 namelist of code_v2 (default 'pimassphys=t, nmassphys=t');
                  dvertex is set by the pipeline itself: dvertex=t for the main curves, dvertex=f for the
                  "without vertex corrections" curve and the delta x (Delta-pole corrections) shift.
                  A v1 executable ignores the &v2 group (then both coincide and are not drawn/used).
  DP_ALPHA, DP_BETA   alphaE1, betaM1 of the main fit [1e-4 fm^3]       (default 11.69, 2.88)
  DP_COV_AB       'c11,c12,c22' fit covariance of (alpha, beta)         (default main fit, rho = -0.83)
  DP_SPINFIT      'a,b,gEE,gMM,gEM,gME' parameter set of the spin-polarisability fit (dash-dotted curves)
  DP_FITS_JSON [+ DP_FIT_KEY]     read alpha, beta, cov from a work/fits fits_L.json (overrides DP_ALPHA...)
  DP_DELTA_JSON [+ DP_SPIN_KEY]   read the spin-fit parameter set from a work/fits delta_region.json
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
ROOT = os.path.dirname(WORK)
CACHE = os.path.join(HERE, 'cache')
RESULTS = os.path.join(HERE, 'results')
PNG = os.path.join(HERE, 'png')
FIGDIR = os.path.join(ROOT, 'paper', 'figs', 'delta')
STYLE = os.path.join(WORK, 'figstyle.mplstyle')
DATA = os.path.join(WORK, 'data')
FITDB = os.path.join(WORK, 'fits', 'db', 'fitdb.csv')

DEFAULT_EXE = os.path.join(WORK, 'code_v2', 'execs', 'ampl.exe')
V1_EXE = os.path.join(WORK, 'code', 'execs', 'ampl.exe')
EXE = os.environ.get('AMPL_EXE', DEFAULT_EXE)
TEMPLATE = os.environ.get('AMPL_TEMPLATE')            # None -> runner.DEFAULT_TEMPLATE
V2OPTS = os.environ.get('DP_V2OPTS', 'pimassphys=t, nmassphys=t')


def v2opts(dvertex):
    """&v2 entries with the vertex-correction switch set explicitly"""
    items = [x.strip() for x in V2OPTS.split(',') if x.strip() and not x.strip().lower().startswith('dvertex')]
    return ', '.join(items + ['dvertex=%s' % ('t' if dvertex else 'f')])


def _vec(name, default):
    v = os.environ.get(name)
    return np.array([float(x) for x in v.split(',')]) if v else np.array(default, float)

MP = 938.272          # proton mass [MeV] for kinematics of the data (as in work/data, work/fits)
MPI = 139.57

# ---------------------------------------------------------------------------------------------
# polarisabilities  P = (alphaE1, betaM1 [1e-4 fm^3]; gE1E1, gM1M1, gE1M2, gM1E2 [1e-4 fm^4])
# ---------------------------------------------------------------------------------------------
# main BChPT fit (fit_report Sect. 0/4.1, F2 + truncation covariance; work/fits/results/fits_L.json
# 'F2_d_all_thB'): alpha = 11.69(47), beta = 2.88(43), rho = -0.83; spin polarisabilities = NNLO values
# of the amplitude (the driver's reference values, i.e. those of the run with the given gM).
ALPHA_MAIN = float(os.environ.get('DP_ALPHA', 11.69))
BETA_MAIN = float(os.environ.get('DP_BETA', 2.88))
_c = _vec('DP_COV_AB', [0.2206, -0.1688, 0.1880])
COV_AB = np.array([[_c[0], _c[1]], [_c[1], _c[2]]])     # fit covariance of (alpha, beta)
# alternatively read the main fit from a fits_L.json-type file of work/fits (rerun with v2):
#   DP_FITS_JSON=work/fits/results/fits_L.json  DP_FIT_KEY=F2_d_all_thB
if os.environ.get('DP_FITS_JSON'):
    import json as _json
    _F = _json.load(open(os.environ['DP_FITS_JSON']))[os.environ.get('DP_FIT_KEY', 'F2_d_all_thB')]
    ALPHA_MAIN, BETA_MAIN = float(_F['P'][0]), float(_F['P'][1])
    COV_AB = np.array(_F['covP'])[:2, :2]
# NNLO spin polarisabilities for gM = 2.97 (analytic, = driver Pref); used only for reporting
SPIN_NNLO = np.array([-3.3351, 2.9427, 0.1787, 1.1469])
# exploratory Delta-region spin-polarisability fit (fit_report Sect. 5.2, bold row:
# A2 Sigma2x+Sigma2z + LEGS Sigma3 below 2pi threshold, gamma_pi = 8.0(1.8), gamma_0 = -0.929(105),
# alpha+beta = 14.0(2), alpha-beta from the main fit; delta_region.json 'S_A2_LEGS_gpi').
# The complete fitted parameter set is used for the dash-dotted curves.
P_SPINFIT = _vec('DP_SPINFIT', [11.2204, 2.8716, -2.9556, 2.5840, -1.3909, 2.7059])
#   or DP_DELTA_JSON=work/fits/results/delta_region.json  DP_SPIN_KEY=S_A2_LEGS_gpi
if os.environ.get('DP_DELTA_JSON'):
    import json as _json
    P_SPINFIT = np.array(_json.load(open(os.environ['DP_DELTA_JSON']))[os.environ.get('DP_SPIN_KEY', 'S_A2_LEGS_gpi')]['P'], float)

# gamma N Delta magnetic coupling and its variation (lighter band)
GM0, GM_LO, GM_HI = 2.97, 2.87, 3.07

# ---------------------------------------------------------------------------------------------
# EFT-truncation band (LMP2015 Eq. (32) type; see delta_plots.md Sect. 2)
#   A_shift = A(P_main) +- delta(omega) * dA_X ,
#   X in {piN loops, piDelta loops, Delta-pole (vertex) corrections, sum of the three}  [LMP2015 Eq. (32)]
#   delta(omega) = omega_cm / LAMBDA_CHI ;  static alphaE1, betaM1 of every dA_X removed (they are fitted)
#   band = pointwise envelope of all shifted observables, combined in quadrature with the
#          propagated (alpha, beta) fit error.
# ---------------------------------------------------------------------------------------------
LAMBDA_CHI = 650.0      # MeV;  Delta/Lambda = 0.45 ~ delta = m_pi/Delta = 0.48 ; alpha_em/Lambda^3 = 2.0e-4 fm^3
SHIFT_SETS = {'piN': ('piN',), 'piD': ('piD',), 'VC': ('VC',), 'sum': ('piN', 'piD', 'VC')}
REMOVE_STATIC_AB = True

# Born + Delta-pole curve: nucleon Born (Dirac + Pauli, s+u) + pi0 t-channel pole + Delta pole (with the
# vertex corrections, i.e. the dressed Delta pole of the NLO Delta-region amplitude), no loops, no contact
# terms, as the blue dashed curves of LMP2015 Figs. 14-18.
TREE_MECHS = ('Born', 'pi0', 'Dpole')

# data sets removed from all Delta-region figures and chi^2 tables
# excluded: frames unknown (user decision 2026-10-14) -- old bremsstrahlung measurements whose lab/c.m. frame of
# the angle and cross section is not established: Gray 1967 (Illinois 67), Genzel 1972 (Bonn 72), Baranov 1966 and
# its 'refined' points (Moscow 66), Nagashima 1964 (Tokyo 64).  Barbiellini 1968 (c.m. stated) is kept.
EXCLUDE_SETS = ('Gray1967', 'Genzel1972', 'Baranov1966', 'Baranov1966_refined', 'Nagashima1964')

# data windows
WIN_THETA = 3.0       # deg, Fig. 1 (vs E at fixed angle), Fig. 3b, Fig. 5
WIN_E = 5.0           # MeV, Fig. 2 and Fig. 3a (vs angle at fixed E): point assigned to nearest panel energy
E_DELTA_MIN = 150.0   # chi^2 table: Delta region E_lab > 150 MeV (as fit_report Sect. 5.1)
E_DELTA_MAX = 350.0   #   ... and E_lab <= 350 MeV (fit_report);  a second column extends to 400 MeV
