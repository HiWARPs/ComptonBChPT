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
  DP_EXTRA        the extra dotted curve of the figures: 'none' (default since 2026-10-17 pm; 'run3' was the default on
                  2026-10-17 am) = no extra curve; 'run3' = the NNLO prediction with the
                  running gammaNDelta couplings and Delta mass of code_v2 (&v2 irun=3: dispersion relations doubly
                  subtracted, values AND slopes fixed at s0 = M_N^2, all six static polarisabilities unchanged;
                  work/notes/delta_running2.md), evaluated at the static polarisabilities of the main curve;
                  'run2' / 'run1' = irun=2 / 1 at the same static polarisabilities (contact terms; exploratory);
                  'novc' = the prediction without the gammaNDelta vertex corrections (the dotted curve until
                  2026-10-17); 'none' = no extra curve.  The chi^2 tables always keep the no-VC column.
  DP_EXTRA_THETA  '1' (default for the running variants) / '0' (default for 'novc', as before): draw the extra
                  curve also in the angular distributions fig_delta_cs_vstheta_a/b
  DP_MAIN         main curve (red solid + dark band): 'pred' (default) = the parameter-free NNLO prediction, all
                  six static polarisabilities at the NNLO values of the amplitude (driver Pref, physical masses);
                  'fit' = alphaE1, betaM1 from the NNLO fit (the main curve before 2026-10-16)
  DP_FITS_JSON [+ DP_FIT_KEY]     fits_L.json of work/fits from which alpha, beta, cov of the 'NNLO fit' curve are
                  read (default: work/fits/results/fits_L.json if present, key F2_d_all_thB)
  DP_ALPHA, DP_BETA, DP_COV_AB    alpha, beta, 'c11,c12,c22' of the NNLO fit, used only when no fits_L.json is found
  DP_SPINFIT / DP_DELTA_JSON [+ DP_SPIN_KEY]   parameter set of the Delta-region spin-polarisability fit; since
                  2026-10-16 no longer drawn, used only as a third parameter set in the validation of the Gram
                  decomposition
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


def v2opts(dvertex, irun=None):
    """&v2 entries with the vertex-correction switch set explicitly.  irun (0-3) overrides any irun/drunning entry of
    DP_V2OPTS (the running of the Delta parameters, code_v2 CHANGES.md Sect. G); irun=None leaves DP_V2OPTS as given,
    i.e. the default amplitude (no running) unless the environment says otherwise (chi2_running.py)."""
    items = [x.strip() for x in V2OPTS.split(',') if x.strip() and not x.strip().lower().startswith('dvertex')]
    if irun is not None:
        items = [x for x in items if not x.lower().startswith(('irun', 'drunning'))] + ['irun=%d' % int(irun)]
    return ', '.join(items + ['dvertex=%s' % ('t' if dvertex else 'f')])


# extra dotted curve of the figures (DP_EXTRA, see the module doc).  'run3' is the variant shown in the paper
# (2026-10-17): same static polarisabilities as the main curve by construction, only the energy dependence of the
# Delta parameters between the threshold and the resonance is added.
EXTRA = os.environ.get('DP_EXTRA', 'none').strip().lower()     # 'none' since 2026-10-17 (user: no running curve in the figures)
assert EXTRA in ('run3', 'run2', 'run1', 'novc', 'none'), 'DP_EXTRA must be run3, run2, run1, novc or none'
EXTRA_IRUN = int(EXTRA[3:]) if EXTRA.startswith('run') else None      # &v2 irun of the extra curve (None: no run)
EXTRA_THETA = os.environ.get('DP_EXTRA_THETA', '0' if EXTRA == 'novc' else '1').strip() == '1'


def _vec(name, default):
    v = os.environ.get(name)
    return np.array([float(x) for x in v.split(',')]) if v else np.array(default, float)

MP = 938.272          # proton mass [MeV] for kinematics of the data (as in work/data, work/fits)
MPI = 139.57

# ---------------------------------------------------------------------------------------------
# polarisabilities  P = (alphaE1, betaM1 [1e-4 fm^3]; gE1E1, gM1M1, gE1M2, gM1E2 [1e-4 fm^4])
# ---------------------------------------------------------------------------------------------
# Main curve (DP_MAIN).  'pred': the parameter-free NNLO prediction -- all six static polarisabilities at the NNLO
# values of the amplitude with the physical masses, i.e. the driver's reference values Pref of the complete run
# (= PREF_NNLO of work/fits/lib/common.py: 11.415, 3.729, -3.246, 2.939, 0.170, 1.133 for gM = 2.97); nothing
# fitted.  'fit': alphaE1, betaM1 of the NNLO fit, spin polarisabilities at the NNLO values (main curve until
# 2026-10-15).  The NNLO fit is drawn in both modes as a thin dashed curve without a band ('NNLO fit').
MAIN = os.environ.get('DP_MAIN', 'pred').strip().lower()
assert MAIN in ('pred', 'fit'), 'DP_MAIN must be pred or fit'

# NNLO fit: alpha, beta (and their covariance) of work/fits/results/fits_L.json 'F2_d_all_thB' (fit_report
# Sect. 0/4.1, F2 + truncation covariance): alpha = 11.59(47), beta = 2.94(43), rho = -0.83.  The file is read by
# default when present (DP_FITS_JSON overrides the path, DP_FIT_KEY the key); the DP_ALPHA/DP_BETA/DP_COV_AB
# values (v1-era defaults 11.69, 2.88) are used only when no file is found -- a warning is printed in that case.
ALPHA_FIT = float(os.environ.get('DP_ALPHA', 11.69))
BETA_FIT = float(os.environ.get('DP_BETA', 2.88))
_c = _vec('DP_COV_AB', [0.2206, -0.1688, 0.1880])
COV_AB = np.array([[_c[0], _c[1]], [_c[1], _c[2]]])     # fit covariance of (alpha, beta)
FITS_JSON = os.environ.get('DP_FITS_JSON') or os.path.join(WORK, 'fits', 'results', 'fits_L.json')
FIT_KEY = os.environ.get('DP_FIT_KEY', 'F2_d_all_thB')
if FITS_JSON.strip().lower() == 'none':          # DP_FITS_JSON=none: use DP_ALPHA, DP_BETA (regression checks)
    FITS_JSON = '/nonexistent'
if os.path.exists(FITS_JSON):
    import json as _json
    _F = _json.load(open(FITS_JSON))[FIT_KEY]
    ALPHA_FIT, BETA_FIT = float(_F['P'][0]), float(_F['P'][1])
    COV_AB = np.array(_F['covP'])[:2, :2]
else:
    print('[config] WARNING: %s not found; NNLO-fit curve uses alpha, beta = %.3f, %.3f from the environment/defaults'
          % (FITS_JSON, ALPHA_FIT, BETA_FIT))
    FITS_JSON = None
# NNLO spin polarisabilities for gM = 2.97 (analytic, = driver Pref); used only for reporting
SPIN_NNLO = np.array([-3.3351, 2.9427, 0.1787, 1.1469])
# exploratory Delta-region spin-polarisability fit (fit_report Sect. 5.2; delta_region.json 'S_A2_LEGS_gpi';
# paper Eq. spinfit).  Since 2026-10-16 this curve is NOT drawn (user decision: the prediction is the main
# result, the gamma fit is dropped); the parameter set serves only as a third point of the validation.
P_SPINFIT = _vec('DP_SPINFIT', [11.2204, 2.8716, -2.9556, 2.5840, -1.3909, 2.7059])
DELTA_JSON = os.environ.get('DP_DELTA_JSON') or os.path.join(WORK, 'fits', 'results', 'delta_region.json')
if os.path.exists(DELTA_JSON):
    import json as _json
    P_SPINFIT = np.array(_json.load(open(DELTA_JSON))[os.environ.get('DP_SPIN_KEY', 'S_A2_LEGS_gpi')]['P'], float)
else:
    DELTA_JSON = None

# gamma N Delta magnetic coupling and its variation (lighter band)
GM0, GM_LO, GM_HI = 2.97, 2.87, 3.07

# ---------------------------------------------------------------------------------------------
# EFT-truncation band (LMP2015 Eq. (32) type; see delta_plots.md Sect. 2)
#   A_shift = A(P_main) +- delta(omega) * dA_X ,
#   X in {piN loops, piDelta loops, Delta-pole (vertex) corrections, sum of the three}  [LMP2015 Eq. (32)]
#   delta(omega) = omega_cm / LAMBDA_CHI ;  band = pointwise envelope of all shifted observables.
#   MAIN = 'pred': nothing is fitted -- the shifts keep their static alphaE1, betaM1 content and no fit error is
#                  added (the band is the pure LMP2015-type envelope);
#   MAIN = 'fit' : static alphaE1, betaM1 of every dA_X removed (they are fitted) and the propagated (alpha, beta)
#                  fit error added in quadrature (band of the figures until 2026-10-15).
# ---------------------------------------------------------------------------------------------
LAMBDA_CHI = 750.0      # MeV; single chiral scale of the paper (of the order of the rho mass; 650 MeV = LMP15 value until 2026-10-16)
SHIFT_SETS = {'piN': ('piN',), 'piD': ('piD',), 'VC': ('VC',), 'sum': ('piN', 'piD', 'VC')}
REMOVE_STATIC_AB = (MAIN == 'fit')
ADD_AB_ERROR = (MAIN == 'fit')

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
