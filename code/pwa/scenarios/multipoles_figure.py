"""Figure of the l=1 non-Born proton Compton multipoles below the pion-production threshold, as dynamical
polarisabilities:  the parameter-free NNLO BChPT prediction (v2 code, static polarisabilities = the code's
physical-mass reference values BV.REF; MAIN curve, solid, dark grey) with the band of the O(p^4) LECs of natural
size, delta alpha_E1 = delta beta_M1 = +-DLEC = +-0.7 (Lensky-Pascalutsa 2010, Eq. 25), propagated through the exact
amplitude (central differences at REF +- 0.7 e_k, added in quadrature; it is drawn in the alpha_E1(w), beta_M1(w)
panels only -- the LECs do not shift the spin multipoles, which is checked at run time); the BChPT main fit
("NNLO fit": alpha_E1, beta_M1 from work/fits/results/fits_L.json['F2_d_all_thB'], spin polarisabilities at NNLO;
dashed, blue); the PWA fit of the paper ("PWA fit": main database = pre-2015 + A2 2022 minus the 3-sigma outlier
TAPS 108.42/133 (update_fits.MAIN_DROP), HIGS 2022 excluded since 2026-10-17 (update_fits.HIGS), legacy
normalisations floated, configuration main_fits.MAIN); and the
published KLP18 Fit 1 and Fit 1'' (PLB 782 (2018) 34, Fig. 4, recovered from the vector graphics).
Until 2026-10-16 the fit was the main (solid) curve, with the band from the (alpha, beta) fit covariance, and the
prediction dashed; that version is kept in work/backups/*_before_predmain.*  The fit-covariance band is still
computed and written to the JSON ('dyn_band_fit').

Conventions (identical for all curves/points; = make_figures.DYN / HMP16 Eq. 3.46 / KLP18 notebook):
  non-Born c.m. multipoles fbar^{l+-}_{rho' rho}(omega) of compton.py (HMP16 Eqs. 3.19-3.22; Phi = T/(8 pi sqrt s)),
  omega = c.m. photon energy = E_lab/sqrt(1 + 2 E_lab/M), plotted against the LAB photon energy E_gamma;
  alpha_E1(w) = [2 f_EE^{1+} + f_EE^{1-}]/w^2,  beta_M1(w) = [2 f_MM^{1+} + f_MM^{1-}]/w^2,
  gamma_E1E1(w) = [f_EE^{1+} - f_EE^{1-}]/w^3,  gamma_M1M1(w) = [f_MM^{1+} - f_MM^{1-}]/w^3,
  gamma_E1M2(w) = 6 f_EM^{1+}/w^3,  gamma_M1E2(w) = 6 f_ME^{1+}/w^3.
Run (one command, from anywhere):
  python3 \
      <release>/code/pwa/scenarios/multipoles_figure.py
Outputs: paper/figs/fig_multipoles.pdf, paper/figs/fig_multipoles_raw.pdf (PNG previews in work/pwa/figs/),
         work/pwa/results/multipoles_figure.json, work/pwa/results/klp18_fig4_points.csv"""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
_ROOT = _os.path.dirname(_WORK)
import os
import sys
import json
import types

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
sys.path.insert(0, PWA)
sys.path.insert(0, HERE)
os.chdir(PWA)

import numpy as np                                    # noqa: E402
import pandas as pd                                   # noqa: E402
import matplotlib                                     # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402
from matplotlib.patches import Patch                  # noqa: E402
from matplotlib.ticker import MultipleLocator         # noqa: E402
import compton as C                                   # noqa: E402
import pwa_model as PM                                # noqa: E402
import update_fits as UF                              # noqa: E402
import plb_inputs as PI                               # noqa: E402
import reproduce_plb as RP                            # noqa: E402
import make_figures as MF                             # noqa: E402  (DYN, dyn_pols, dyn_from_mult)
from main_fits import MAIN                            # noqa: E402
import bchpt_v2_multipoles as BV                      # noqa: E402
import klp18_fig4_digitise as KD                      # noqa: E402

plt.style.use((_WORK + '/figstyle.mplstyle'))
FIGDIR = (_ROOT + '/paper/figs')
PNGDIR = os.path.join(PWA, 'figs')                     # PNG previews (not for the paper)
OUTJSON = os.path.join(PWA, 'results', 'multipoles_figure.json')

# ---------------------------------------------------------------- inputs
# main BChPT fit (F2 + truncation covariance, 3-sigma outlier removed), read from work/fits/results/fits_L.json
# ('F2_d_all_thB'); the literal values are only the fallback if that file is missing
_FL = (_WORK + '/fits/results/fits_L.json')
try:
    _m = json.load(open(_FL))['F2_d_all_thB']
    _c = np.array(_m['covP'])[:2, :2]
    BCHPT_FIT = dict(alphaE1=round(float(_m['P'][0]), 4), betaM1=round(float(_m['P'][1]), 4))
    BCHPT_COV = dict(sa=round(float(np.sqrt(_c[0, 0])), 4), sb=round(float(np.sqrt(_c[1, 1])), 4),
                     rho=round(float(_c[0, 1] / np.sqrt(_c[0, 0] * _c[1, 1])), 4))
    print('BChPT main fit from fits_L.json:', BCHPT_FIT, BCHPT_COV)
except (OSError, KeyError):
    BCHPT_FIT = dict(alphaE1=11.672, betaM1=2.876)
    BCHPT_COV = dict(sa=0.47, sb=0.43, rho=-0.84)
    print('WARNING: fits_L.json not found, using fallback BChPT main fit', BCHPT_FIT)
LEG = UF.LEGACY_NORM   # Baranov et al. 2001 delta (was 10%, Goldansky 8%)
KEYS6 = [k[0] for k in MF.DYN]                         # alphaE1, betaM1, gE1E1, gM1M1, gE1M2, gM1E2
L1 = [('EE', 1, '+'), ('EE', 1, '-'), ('MM', 1, '+'), ('MM', 1, '-'), ('EM', 1, '+'), ('ME', 1, '+')]
EGRID = np.arange(0.5, 160.001, 0.5)
COL = dict(pred='#444444', fit='#1f4e9c', band='#9fb6dd', pwa='#1f4e9c', plb1='#5f5f5f', plb1pp='#e08e0b')
BAND_ALPHA = 0.6
DLEC = 0.7      # natural size of the O(p^4) LECs (LP 2010, Eq. 25): band of the prediction (as in work/fits/scripts/s06_figures.py)
SCALAR = ('alphaE1', 'betaM1')
LABEL = dict(alphaE1=r'$\alpha_{E1}(\omega)$', betaM1=r'$\beta_{M1}(\omega)$', gE1E1=r'$\gamma_{E1E1}(\omega)$',
             gM1M1=r'$\gamma_{M1M1}(\omega)$', gE1M2=r'$\gamma_{E1M2}(\omega)$', gM1E2=r'$\gamma_{M1E2}(\omega)$')
UNIT = dict(alphaE1=3, betaM1=3, gE1E1=4, gM1M1=4, gE1M2=4, gM1E2=4)


def dyn_curves(f, E):
    """dynamical polarisabilities (1e-4 fm^n) from a dict of real non-Born multipoles (fm) on the grid E."""
    return {key: MF.dyn_from_mult(f, E, comb, n) for key, _, comb, n, _ in MF.DYN}


def smooth_lowE(E, y, lo=12.0, hi=30.0, top=70.0, deg=4):
    """Between lo and hi MeV the LoopTools evaluation of the loops loses precision (cancellations ~E^-4;
    after division by omega^3 the noise reaches 0.04 in gamma_M1E2 at 15-30 MeV, < 1e-3 above 30 MeV).
    Replace that interval by a polynomial fitted to E <= lo (contact-term form used by the code below
    elex = 15 MeV) and hi <= E <= top.  Returns the smoothed array, the max |change| in (lo, hi) and the
    max residual of the polynomial on the fitted points."""
    sel = (E <= lo) | ((E >= hi) & (E <= top))
    c = np.polyfit(E[sel], y[sel], deg)
    ys = y.copy()
    rep = (E > lo) & (E < hi)
    ys[rep] = np.polyval(c, E[rep])
    return ys, float(np.max(abs(ys - y)[rep])), float(np.max(abs(np.polyval(c, E[sel]) - y[sel])))


# ---------------------------------------------------------------- BChPT (v2)
def bchpt():
    ref = BV.REF
    from concurrent.futures import ThreadPoolExecutor            # fill the run cache in parallel
    jobs = [(None, False), (dict(alphaE1=ref['alphaE1'] + 1.0), False), (dict(betaM1=ref['betaM1'] + 1.0), False),
            (BCHPT_FIT, False), (None, True)]
    jobs += [(dict(alphaE1=ref['alphaE1'] + sgn * DLEC), False) for sgn in (+1, -1)]      # LEC band of the prediction
    jobs += [(dict(betaM1=ref['betaM1'] + sgn * DLEC), False) for sgn in (+1, -1)]
    with ThreadPoolExecutor(5) as ex:
        list(ex.map(lambda j: BV.multipoles(EGRID, j[0], born_only=j[1]), jobs))
    nb0, born = BV.nonborn(EGRID, None)                                  # NNLO prediction (LECs = 0)
    nba, _ = BV.nonborn(EGRID, dict(alphaE1=ref['alphaE1'] + 1.0))
    nbb, _ = BV.nonborn(EGRID, dict(betaM1=ref['betaM1'] + 1.0))
    nbm, _ = BV.nonborn(EGRID, BCHPT_FIT)                                # main fit, direct run (check)
    assert np.allclose(nb0['_ref'], [ref[k] for k in KEYS6], atol=2e-5), nb0['_ref']    # code's own reference values
    keys = [k for k in nb0 if not str(k).startswith('_')]
    Ja = {k: (nba[k] - nb0[k]).real for k in keys}                       # exact: amplitude linear in the LECs
    Jb = {k: (nbb[k] - nb0[k]).real for k in keys}
    da, db = BCHPT_FIT['alphaE1'] - ref['alphaE1'], BCHPT_FIT['betaM1'] - ref['betaM1']
    f0 = {k: nb0[k].real for k in keys}
    fm = {k: f0[k] + da * Ja[k] + db * Jb[k] for k in keys}
    lin_check = max(float(np.max(abs(fm[k] - nbm[k].real))) / float(np.max(abs(nbm[k].real))) for k in L1)
    d0, dJa, dJb = dyn_curves(f0, EGRID), dyn_curves(Ja, EGRID), dyn_curves(Jb, EGRID)
    # LEC band of the prediction: central differences of the multipoles at REF +- DLEC e_k through the exact
    # amplitude (direct runs), added in quadrature.  The amplitude is linear in the LECs, so ca = DLEC * Ja exactly
    # (checked below, 'lec_lin_check'); the dynamical polarisabilities are linear in the multipoles.
    nbap, _ = BV.nonborn(EGRID, dict(alphaE1=ref['alphaE1'] + DLEC))
    nbam, _ = BV.nonborn(EGRID, dict(alphaE1=ref['alphaE1'] - DLEC))
    nbbp, _ = BV.nonborn(EGRID, dict(betaM1=ref['betaM1'] + DLEC))
    nbbm, _ = BV.nonborn(EGRID, dict(betaM1=ref['betaM1'] - DLEC))
    ca = {k: (nbap[k] - nbam[k]).real / 2 for k in keys}
    cb = {k: (nbbp[k] - nbbm[k]).real / 2 for k in keys}
    lec_lin_check = max(float(np.max(abs(c[k] - DLEC * J[k]))) / max(float(np.max(abs(DLEC * J[k]))), 1e-30)
                        for c, J in ((ca, Ja), (cb, Jb)) for k in L1)
    dca, dcb = dyn_curves(ca, EGRID), dyn_curves(cb, EGRID)
    # how the LEC shifts propagate into the dynamical polarisabilities: d alpha_E1(w)/d alpha_E1 should be 1 and
    # d beta_M1(w)/d alpha_E1 = 0 at all w (and vice versa), the spin multipoles should not move at all
    jac_check = {}
    for k in KEYS6:
        ja, jb = dJa[k], dJb[k]
        jac_check[k] = dict(dJa_max=float(np.max(abs(ja))), dJb_max=float(np.max(abs(jb))),
                            dJa_minus_1_max=float(np.max(abs(ja - 1))), dJb_minus_1_max=float(np.max(abs(jb - 1))),
                            dJa_at_150=float(np.interp(150.0, EGRID, ja)), dJb_at_150=float(np.interp(150.0, EGRID, jb)))
    smooth_info = {}
    for k in KEYS6:
        d0[k], ch, res = smooth_lowE(EGRID, d0[k])
        smooth_info[k] = dict(max_change=ch, max_fit_residual=res)
    sa, sb, rho = BCHPT_COV['sa'], BCHPT_COV['sb'], BCHPT_COV['rho']
    out = dict(E=EGRID, nnlo={}, main={}, band={}, lec={}, raw_main=fm, raw_nnlo=f0, raw_Ja=Ja, raw_Jb=Jb, born=born,
               nb_main_direct=nbm, ref=ref, lin_check=lin_check, lec_lin_check=lec_lin_check, jac_check=jac_check,
               smooth=smooth_info, dlec=DLEC)
    stat_nnlo = dict(ref)
    stat_main = dict(ref, **BCHPT_FIT)
    for k in KEYS6:
        out['nnlo'][k] = np.concatenate([[stat_nnlo[k]], d0[k]])
        main = d0[k] + da * dJa[k] + db * dJb[k]
        out['main'][k] = np.concatenate([[stat_main[k]], main])
        var = (sa * dJa[k])**2 + (sb * dJb[k])**2 + 2 * rho * sa * sb * dJa[k] * dJb[k]
        e0 = {'alphaE1': sa, 'betaM1': sb}.get(k, 0.0)
        out['band'][k] = np.concatenate([[e0], np.sqrt(var)])        # (alpha, beta) covariance of the fit (not drawn)
        l0 = DLEC if k in SCALAR else 0.0                            # static values: shifted by exactly +-DLEC
        out['lec'][k] = np.concatenate([[l0], np.sqrt(dca[k]**2 + dcb[k]**2)])   # LEC band of the prediction (drawn)
    out['Ec'] = np.concatenate([[0.0], EGRID])
    # raw multipole bands: fit covariance (not drawn) and LEC band of the prediction (drawn)
    out['raw_band'] = {k: np.sqrt((sa * Ja[k])**2 + (sb * Jb[k])**2 + 2 * rho * sa * sb * Ja[k] * Jb[k]) for k in keys}
    out['raw_lec'] = {k: np.sqrt(ca[k]**2 + cb[k]**2) for k in keys}
    return out


# ---------------------------------------------------------------- PWA fits
def load_leg(per_beamtime=False):
    d = UF.load(per_beamtime=per_beamtime)
    for s, e in LEG.items():
        m = d.dataset == s
        d.loc[m, 'norm_group'] = s
        d.loc[m, 'norm_err'] = e
    return d


def pwa_fit(**kw):
    cfg = dict(MAIN)
    cfg.update(kw)
    m = UF.build_model(UF.select(load_leg()), **cfg)
    return m, m.fit()


def pwa_points(m, r):
    dp = MF.dyn_pols(m, r)
    st = {k: (float(r.der[k]), float(r.der_err[k])) for k in KEYS6}
    raw = []
    for ib, E in enumerate(m.bins):
        f = m.multipoles_bin(ib, r.p)
        e = m.multipole_cov(ib, r.cov, r.p)
        row = {'E': float(E)}
        row.update({PM.ALLKEYS[j]: (float(f[PM.ALLKEYS[j]]), float(e[j])) for j in range(len(PM.ALLKEYS))})
        raw.append(row)
    return dict(dyn={k: v.tolist() for k, v in dp.items()}, static=st, raw=raw,
                chi2=float(r.chi2), npts=int(r.npts), der={k: float(v) for k, v in r.der.items()},
                der_err={k: float(v) for k, v in r.der_err.items()})


def klp_repro():
    d0 = PI.load_klp2018()
    out = {}
    for tag, which in (('Fit 1', '1'), ("Fit 1''", "1''")):
        m = RP.make_model(RP.select(d0, which))
        r = m.fit()
        out[tag] = pwa_points(m, r)
    return out


# ---------------------------------------------------------------- validation
def validate(B):
    V = {}
    # (1) projection conventions: Born+pi0 multipoles of pwa_gaga (Fortran) vs compton.py (HMP16), 14 multipoles
    Ev = np.array([5., 20., 59., 109., 150.])
    bF = BV.multipoles(Ev, born_only=True)
    xg, wg = np.polynomial.legendre.leggauss(16)
    th = np.degrees(np.arccos(xg))
    dev = []
    for i, e in enumerate(Ev):
        f = C.multipoles_from_pw(C.pw_project(C.phis_from_M(C.born_amplitudes(np.full_like(th, e), th, 'cm')), xg, wg, 5))
        ks = [k for k in bF if not str(k).startswith('_')]
        dev.append(max(float(abs(f[k] - bF[k][i]) / abs(bF[k][i])) for k in ks if abs(bF[k][i]) > 1e-12))
    V['born_projection_rel_dev'] = dict(E=Ev.tolist(), max_rel_dev=dev, n_multipoles=14)
    # (2) observables: direct v2 amplitude (fit driver) vs exact Born + BChPT non-Born multipoles through the
    #     PWA machinery: (a) all multipoles with J <= 5/2, (b) the 10 multipoles kept in the PWA (pwa_model.unit_phis)
    Es, ths = [59., 89., 109., 135., 150.], [30., 60., 90., 120., 150.]
    pts = [(f'{o}_{int(E)}_{int(t)}', E, t, o, 'LAB') for E in Es for t in ths for o in ('CS', 'S3')]
    dr = BV.driver(pts, BCHPT_FIT)
    rows = []
    for E in Es:
        i = int(np.argmin(abs(B['E'] - E)))
        f = {k: B['nb_main_direct'][k][i].real for k in B['nb_main_direct'] if not str(k).startswith('_')}
        k = C.kin(E, np.array(ths), 'lab')
        PB = np.real(C.phis_from_M(C.born_amplitudes(E, np.degrees(k['thcm']), 'cm')))
        Pa = PB + np.real(C.phis_from_pw(C.pw_from_multipoles(f, 5), k['thcm']))
        fv = np.array([f[kk] for kk in PM.ALLKEYS])
        Pb = PB + np.einsum('imk,m->ik', PM.unit_phis(k['thcm']), fv)
        oa, ob = C.obs_from_phis(Pa, k), C.obs_from_phis(Pb, k)
        for j, t in enumerate(ths):
            cs, s3 = dr[f'CS_{int(E)}_{int(t)}'], dr[f'S3_{int(E)}_{int(t)}']
            rows.append(dict(E=E, th=t, cs=cs, s3=s3, dcs_J52=oa['cs_lab'][j] / cs - 1, ds3_J52=oa['S3'][j] - s3,
                             dcs_pwa10=ob['cs_lab'][j] / cs - 1, ds3_pwa10=ob['S3'][j] - s3))
    R = pd.DataFrame(rows)
    V['observables'] = dict(table=R.round(7).to_dict('records'),
                            max_rel_dcs_J52=float(abs(R.dcs_J52).max()), max_abs_dS3_J52=float(abs(R.ds3_J52).max()),
                            max_rel_dcs_pwa10=float(abs(R.dcs_pwa10).max()), max_abs_dS3_pwa10=float(abs(R.ds3_pwa10).max()),
                            max_rel_dcs_J52_below140=float(abs(R[R.E < 140].dcs_J52).max()))
    # (3) cross-check with the PWA projection (compton.pw_project) of the helicity amplitudes of the 2015
    #     code (results/bchpt_multipoles.csv): pwa_gaga in 2015-code mode (no vertex corr., one-sided piDelta
    #     subtraction, isospin-symmetric masses), pure NNLO amplitude (input = the code's reference values)
    v1opt = "dvertex=f, consts='LMP15', iosmode=1"
    Ec = np.array([30., 59., 89., 109., 135.])
    tmp = BV.multipoles(Ec[:1], None, v2=v1opt)
    pref = dict(zip(['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2'], tmp['_ref']))
    nbv1, _ = BV.nonborn(Ec, pref, v2=v1opt)
    csv = pd.read_csv('results/bchpt_multipoles.csv')
    dv = {}
    for k in L1:
        nm = f'f{k[0]}{k[1]}{"p" if k[2] == "+" else "m"}_re'
        ref_2015 = np.interp(Ec, csv.E_lab.values, csv[nm].values)
        dv['%s%d%s' % k] = (nbv1[k].real / ref_2015 - 1).tolist()
    V['lmp2015_projection_crosscheck'] = dict(E=Ec.tolist(), rel_dev=dv, note='pwa_gaga (2015-code mode) vs compton.pw_project of the 2015-code helicity amplitudes (20-pt endpoint-avoiding Gauss)')
    # (4) KLP18 parametrisation of the PWA (pwa_model.Model.multipoles_at) with the BChPT static values and
    #     f^R_i = 0  vs the BChPT curve at low energy; BChPT-implied residual functions must stay finite as E -> 0
    st = dict(BV.REF, **BCHPT_FIT)
    g0 = -(st['gE1E1'] + st['gM1M1'] + st['gE1M2'] + st['gM1E2'])
    ns = types.SimpleNamespace(gnames=['D', 'gEE', 'gME', 'nf', 'ng'], ng=5, rnames=['fR1', 'fR2', 'fR3', 'fR4', 'fR5', 'fR6'],
                               nr=6, fix_gEM=st['gE1M2'], sr_mode='chi2', ab=14.0, sig_ab=0.2, g0=-0.929, sig_g0=0.105,
                               l2='zero', l2_scale=1.0, bchpt=None)
    p = np.zeros(11)
    p[:5] = [st['alphaE1'] - st['betaM1'], st['gE1E1'], st['gM1E2'], (st['alphaE1'] + st['betaM1'] - 14.0) / 0.2, (g0 + 0.929) / 0.105]
    El = np.array([2., 5., 10., 20., 30., 40., 59.])
    fpar = PM.Model.multipoles_at(ns, El, 0, p)
    dpar = dyn_curves({k: np.asarray(v, float) for k, v in fpar.items()}, El)
    comp = {}
    for k in KEYS6:
        comp[k] = dict(param=[float(x) for x in dpar[k]], bchpt=[float(np.interp(e, B['Ec'], B['main'][k])) for e in El])
    V['param_vs_bchpt_lowE'] = dict(E=El.tolist(), values=comp)
    fm = B['raw_main']
    Eb = np.array([2., 5., 10., 20., 30., 40., 59., 69., 79., 89., 99., 109., 117., 127., 135., 143., 150.])
    fr = PM.bchpt_residuals(lambda key, E: np.interp(E, B['E'], fm[key]), Eb,
                            dict(a=st['alphaE1'], b=st['betaM1'], gEE=st['gE1E1'], gMM=st['gM1M1'], gEM=st['gE1M2'], gME=st['gM1E2']))
    V['bchpt_residual_functions'] = dict(E=Eb.tolist(), **{k: v.tolist() for k, v in fr.items()},
                                         note='KLP18 Eq. (3) residuals implied by the BChPT main-fit multipoles; 1e-4 fm^3 (fR1-4), 1e-4 fm^4 (fR5,6)')
    return V


# ---------------------------------------------------------------- figures
def errpoints(ax, x, y, e, col, mk, filled, label=None, ms=3.2, z=4):
    ax.errorbar(x, y, e, ls='none', marker=mk, ms=ms, color=col, mfc=col if filled else 'white', mew=0.7,
                elinewidth=0.7, label=label, zorder=z)


# theory legend entries, as in Figs. cs and sigma3 of the paper (work/fits/scripts/s06_figures.py: TH_LAB, th_handles)
TH_LAB = ['NNLO B$\\chi$PT', 'NNLO fit']      # band = LEC uncertainty, stated in the caption (user, 2026-10-17)


def th_handles():
    return [(Patch(facecolor=COL['band'], alpha=BAND_ALPHA, lw=0), Line2D([], [], color=COL['pred'], lw=1.2)),
            Line2D([], [], color=COL['fit'], ls='--', lw=1.0)]


PANELS = [['alphaE1', 'gE1E1', 'gE1M2'], ['betaM1', 'gM1M1', 'gM1E2']]
YLIM = dict(alphaE1=(8.6, 20.4), betaM1=(0.0, 11.8), gE1E1=(-15.6, 1.6), gM1M1=(-0.6, 6.9), gE1M2=(-0.75, 1.25),
            gM1E2=(-0.4, 4.4))
YTICK = dict(alphaE1=(2, 1), betaM1=(2, 1), gE1E1=(5, 1), gM1M1=(2, 1), gE1M2=(0.5, 0.25), gM1E2=(1, 0.5))
THR = 144.68                                           # pi0 p threshold (lab); pi+ n at 151.4 MeV
OFF = dict(pwa=-2.6, plb1=0.0, plb1pp=2.6)
W2 = 510 / 72.27                                       # double-column width [in]: included without scaling


def fig_main(B, P, K):
    fig, axs = plt.subplots(2, 3, figsize=(W2, 4.45), sharex=True, layout='constrained')
    fig.get_layout_engine().set(hspace=0.0, wspace=0.03, h_pad=0.02, w_pad=0.02)
    for r, row in enumerate(PANELS):
        for c, k in enumerate(row):
            ax = axs[r, c]
            E = B['Ec']
            p, s = B['nnlo'][k], B['lec'][k]
            if k in SCALAR:                                   # LEC band: scalar polarisabilities only
                ax.fill_between(E, p - s, p + s, color=COL['band'], alpha=BAND_ALPHA, lw=0, zorder=1)
            # NNLO fit (dashed, drawn below the prediction: in the spin panels the two coincide)
            ax.plot(E, B['main'][k], '--', color=COL['fit'], lw=1.0, zorder=2)
            ax.plot(E, p, '-', color=COL['pred'], lw=1.2, zorder=2.5)   # NNLO prediction (main curve)
            # KLP18 published (open symbols)
            for tag, key, mk in (('Fit 1', 'plb1', 's'), ("Fit 1''", 'plb1pp', 'D')):
                d = K[(K.fit == tag) & (K.key == k)]
                errpoints(ax, d.E_lab + OFF[key], d.value, d.err, COL[key], mk, False, ms=2.9 if mk == 's' else 2.6, z=3)
            # updated PWA (filled circles), static value at E = 0
            v = np.array(P['dyn'][k])
            errpoints(ax, v[:, 0] + OFF['pwa'], v[:, 1], v[:, 2], COL['pwa'], 'o', True, z=5)
            sv, se = P['static'][k]
            if k == 'gE1M2':
                ax.plot(OFF['pwa'], 0.0, marker='o', ms=3.2, mfc='white', mec=COL['pwa'], mew=0.9, ls='none', zorder=6)
                ax.annotate(r'$\gamma_{E1M2}\equiv0$ fixed' + '\n' + '(PWA, KLP18)', xy=(1.0, -0.04), xytext=(12, -0.66),
                            fontsize=8, color='0.25', va='bottom',
                            arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4', shrinkA=1, shrinkB=1))
            else:
                errpoints(ax, [OFF['pwa']], [sv], [se], COL['pwa'], 'o', True, z=5)
            ax.set_ylim(*YLIM[k])
            ax.yaxis.set_major_locator(MultipleLocator(YTICK[k][0]))
            ax.yaxis.set_minor_locator(MultipleLocator(YTICK[k][1]))
            ax.set_xlim(-7, 157)
            ax.text(0.04, 0.95, LABEL[k], transform=ax.transAxes, va='top', ha='left', fontsize=9)
            ax.set_ylabel(rf'[$10^{{-4}}$ fm$^{UNIT[k]}$]', labelpad=1.5)
            ax.axvline(THR, color='0.55', lw=0.5, ls=':', zorder=0)
            if r == 1:
                ax.set_xlabel(r'$E_\gamma$ [MeV]')
            ax.set_xticks([0, 50, 100, 150])
    h = th_handles() + [
         Line2D([], [], color=COL['pwa'], marker='o', ms=3.2, ls='none'),
         Line2D([], [], color=COL['plb1'], marker='s', ms=2.9, mfc='white', ls='none'),
         Line2D([], [], color=COL['plb1pp'], marker='D', ms=2.6, mfc='white', ls='none')]
    lab = TH_LAB + ['PWA fit', 'KLP18, Fit 1', r"KLP18, Fit 1$''$"]
    fig.legend(h, lab, loc='outside upper center', ncol=5, handletextpad=0.4, columnspacing=1.6)
    for ext in ('pdf', 'png'):
        with plt.rc_context({'savefig.bbox': 'standard'}):          # keep the full 510 pt width (no cropping)
            fig.savefig(os.path.join(FIGDIR if ext == 'pdf' else PNGDIR, f'fig_multipoles.{ext}'), dpi=300)
    plt.close(fig)


RAW = [[('EE', 1, '+'), ('EE', 1, '-'), ('EM', 1, '+')], [('MM', 1, '+'), ('MM', 1, '-'), ('ME', 1, '+')]]
RAWLAB = {k: r'$\bar{f}^{\,%d%s}_{%s}$' % (k[1], k[2], k[0]) for k in L1}


def fig_raw(B, P, Krep):
    fig, axs = plt.subplots(2, 3, figsize=(W2, 4.45), sharex=True, layout='constrained')
    fig.get_layout_engine().set(hspace=0.0, wspace=0.03, h_pad=0.02, w_pad=0.02)
    E = B['E']
    for r, row in enumerate(RAW):
        for c, k in enumerate(row):
            ax = axs[r, c]
            p, s = B['raw_nnlo'][k] * 1e4, B['raw_lec'][k] * 1e4
            ax.fill_between(E, p - s, p + s, color=COL['band'], alpha=BAND_ALPHA, lw=0, zorder=1)
            ax.plot(np.r_[0, E], np.r_[0, B['raw_main'][k] * 1e4], '--', color=COL['fit'], lw=1.0, zorder=2)
            ax.plot(np.r_[0, E], np.r_[0, p], '-', color=COL['pred'], lw=1.2, zorder=2.5)
            for tag, key, mk in (('Fit 1', 'plb1', 's'), ("Fit 1''", 'plb1pp', 'D')):
                rr = Krep[tag]['raw']
                errpoints(ax, [x['E'] + OFF[key] for x in rr], [x[k][0] * 1e4 for x in rr], [x[k][1] * 1e4 for x in rr],
                          COL[key], mk, False, ms=2.9 if mk == 's' else 2.6, z=3)
            rr = P['raw']
            errpoints(ax, [x['E'] + OFF['pwa'] for x in rr], [x[k][0] * 1e4 for x in rr], [x[k][1] * 1e4 for x in rr],
                      COL['pwa'], 'o', True, z=5)
            ax.axhline(0, color='0.5', lw=0.4, zorder=0)
            ax.axvline(THR, color='0.55', lw=0.5, ls=':', zorder=0)
            ax.set_xlim(-5, 157)
            ax.text(0.04, 0.95, RAWLAB[k], transform=ax.transAxes, va='top', ha='left', fontsize=9)
            ax.set_ylabel(r'[$10^{-4}$ fm]', labelpad=1.5)
            if r == 1:
                ax.set_xlabel(r'$E_\gamma$ [MeV]')
            ax.set_xticks([0, 50, 100, 150])
    h = th_handles() + [
         Line2D([], [], color=COL['pwa'], marker='o', ms=3.2, ls='none'),
         Line2D([], [], color=COL['plb1'], marker='s', ms=2.9, mfc='white', ls='none'),
         Line2D([], [], color=COL['plb1pp'], marker='D', ms=2.6, mfc='white', ls='none')]
    lab = TH_LAB + ['PWA fit', 'KLP18 Fit 1 (reanalysis)', r"KLP18 Fit 1$''$ (reanalysis)"]
    fig.legend(h, lab, loc='outside upper center', ncol=5, handletextpad=0.4, columnspacing=1.0)
    for ext in ('pdf', 'png'):
        with plt.rc_context({'savefig.bbox': 'standard'}):
            fig.savefig(os.path.join(FIGDIR if ext == 'pdf' else PNGDIR, f'fig_multipoles_raw.{ext}'), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- main
def tolist(x):
    if isinstance(x, dict):
        return {('%s%d%s' % k if isinstance(k, tuple) else k): tolist(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [tolist(v) for v in x]
    if isinstance(x, np.ndarray):
        return np.real_if_close(x).tolist() if not np.iscomplexobj(x) else [[float(z.real), float(z.imag)] for z in x]
    if isinstance(x, (np.floating, np.integer)):
        return float(x)
    return x


def main():
    K = KD.digitise()
    K.to_csv(os.path.join(PWA, 'results', 'klp18_fig4_points.csv'), index=False)
    # gamma_E1M2 = 0 was fixed in KLP18 (no error bar at E = 0 in their Fig. 4)
    K = pd.concat([K, pd.DataFrame([dict(fit=f, key='gE1M2', E_lab=0.0, value=0.0, err=0.0, clipped=False)
                                    for f in ("Fit 1", "Fit 1''")])], ignore_index=True)
    B = bchpt()
    m, r = pwa_fit()
    P = pwa_points(m, r)
    # must agree with results/main_fits_legacynorm.json['MAIN'] (same selection, same configuration)
    _ref = json.load(open(os.path.join(PWA, 'results', 'main_fits_legacynorm.json')))['MAIN']
    assert abs(P['der']['alpha_m_beta'] - _ref['val']['alpha_m_beta']) < 1e-3 and P['npts'] == _ref['npts'], (P['der'], _ref['val'])
    variants = {}
    for tag, kw in (('bchptctr', dict(prior_center='bchpt')), ('w30', dict(prior_w=(30., 1000., 1000., 30., 3., 3.)))):
        mv, rv = pwa_fit(**kw)
        variants[tag] = pwa_points(mv, rv)
    Krep = klp_repro()
    V = validate(B)
    # PWA vs BChPT at the bins (pulls, PWA errors only)
    pulls = {}
    for k in KEYS6:
        v = np.array(P['dyn'][k])
        bc = np.interp(v[:, 0], B['Ec'], B['main'][k])
        pulls[k] = dict(E=v[:, 0].tolist(), pwa=v[:, 1].tolist(), err=v[:, 2].tolist(), bchpt=bc.tolist(),
                        pull=((v[:, 1] - bc) / v[:, 2]).tolist())
    # published vs reproduced KLP18
    kcmp = {}
    for tag in ("Fit 1", "Fit 1''"):
        kcmp[tag] = {}
        for k in KEYS6:
            d = K[(K.fit == tag) & (K.key == k) & (K.E_lab > 1)].sort_values('E_lab')
            rep = np.array(Krep[tag]['dyn'][k])
            kcmp[tag][k] = dict(E=d.E_lab.tolist(), published=d.value.round(4).tolist(), published_err=d.err.round(4).tolist(),
                                reproduced=rep[:, 1].round(4).tolist(), reproduced_err=rep[:, 2].round(4).tolist())
    fig_main(B, P, K)
    fig_raw(B, P, Krep)
    out = dict(
        conventions=__doc__.split('Run (one command')[0],
        bchpt=dict(settings="work/code_v2 ampl.exe, dvertex=t, pimassphys=t, nmassphys=t, consts='LMP15'; "
                            "baryon.dat M_N = 0.93827 GeV; multipoles from pwa_gaga (16-pt Gauss, J<=5/2) minus Born+pi0",
                   main_fit=BCHPT_FIT, cov=BCHPT_COV, reference_nnlo=B['ref'], linearity_check_rel=B['lin_check'],
                   lowE_smoothing=B['smooth'],
                   dlec=B['dlec'], lec_linearity_check_rel=B['lec_lin_check'], lec_jacobian_check=B['jac_check'],
                   figure='main curve = NNLO prediction (dyn_nnlo) with the LEC band (dyn_band_lec, scalar panels only); '
                          'dashed = NNLO fit (dyn_main); the fit-covariance band (dyn_band_fit) is not drawn',
                   E=B['Ec'], dyn_main=B['main'], dyn_band_fit=B['band'], dyn_nnlo=B['nnlo'], dyn_band_lec=B['lec'],
                   E_raw=B['E'], raw_main_fm={'%s%d%s' % k: B['raw_main'][k] for k in L1},
                   raw_nnlo_fm={'%s%d%s' % k: B['raw_nnlo'][k] for k in L1},
                   raw_band_fit_fm={'%s%d%s' % k: B['raw_band'][k] for k in L1},
                   raw_band_lec_fm={'%s%d%s' % k: B['raw_lec'][k] for k in L1}),
        pwa_main=dict(config='main_fits.MAIN + floating normalisations of the legacy sets ' + str(LEG) +
                             '; main database = pre-2015 + A2 2022 minus ' + str(UF.MAIN_DROP) + ', HIGS 2022 excluded (2026-10-17)', **P),
        pwa_prior_variants=variants,
        klp18_published=dict(source='KLP18 Fig. 4 (arXiv source polarizabilities_NEW_v2.pdf), vector-digitised; errors as published '
                                    '(global parameters only, residual functions fixed)',
                             points=K.round(5).to_dict('records')),
        klp18_published_vs_reproduced=kcmp,
        klp18_reproduced=Krep,
        pwa_vs_bchpt_at_bins=pulls,
        validation=V)
    json.dump(tolist(out), open(OUTJSON, 'w'), indent=1)
    # ---------------- printout
    print('BChPT reference (NNLO) static values:', B['ref'])
    print('linearity check (direct main-fit run vs NNLO + LEC Jacobians), max rel:', B['lin_check'])
    print('LEC band: central differences at REF +- %.1f e_k vs %.1f * Jacobian, max rel:' % (DLEC, DLEC), B['lec_lin_check'])
    print('LEC Jacobians of the dynamical polarisabilities (d/d alpha_E1, d/d beta_M1 of the static shifts):')
    for k in KEYS6:
        j = B['jac_check'][k]
        print('    %-8s max|dJa| = %.2e  max|dJb| = %.2e  max|dJa-1| = %.2e  max|dJb-1| = %.2e  at 150 MeV: %.6f %.6f'
              % (k, j['dJa_max'], j['dJb_max'], j['dJa_minus_1_max'], j['dJb_minus_1_max'], j['dJa_at_150'], j['dJb_at_150']))
    print('LEC band (drawn) at E = 0, 59, 109, 150 MeV:',
          {k: [round(float(np.interp(e, B['Ec'], B['lec'][k])), 4) for e in (0, 59, 109, 150)] for k in KEYS6})
    print('low-E smoothing (12-30 MeV):', {k: round(v['max_change'], 4) for k, v in B['smooth'].items()})
    print('PWA MAIN:', {k: (round(P['der'][k], 3), round(P['der_err'][k], 3)) for k in KEYS6})
    print('validation (1) Born projection max rel dev:', V['born_projection_rel_dev']['max_rel_dev'])
    o = V['observables']
    print('validation (2) dcs/S3 J<=5/2: %.2e / %.2e ; PWA-10: %.2e / %.2e ; J<=5/2 below 140 MeV: %.2e'
          % (o['max_rel_dcs_J52'], o['max_abs_dS3_J52'], o['max_rel_dcs_pwa10'], o['max_abs_dS3_pwa10'], o['max_rel_dcs_J52_below140']))
    print('validation (3) 2015-code cross-check rel dev:', {k: [f'{x:.1e}' for x in v] for k, v in V['lmp2015_projection_crosscheck']['rel_dev'].items()})
    print('validation (4) param(fR=0) vs BChPT:')
    for k, v in V['param_vs_bchpt_lowE']['values'].items():
        print('   ', k, ' '.join(f'{a:.3f}/{b:.3f}' for a, b in zip(v['param'], v['bchpt'])))
    print('    BChPT-implied residuals:')
    for k in ('fR1', 'fR2', 'fR3', 'fR4', 'fR5', 'fR6'):
        print('   ', k, ' '.join(f'{x:6.2f}' for x in V['bchpt_residual_functions'][k]))
    print('pulls PWA vs BChPT main (PWA errors):')
    for k in KEYS6:
        print('   ', k, ' '.join(f'{x:5.1f}' for x in pulls[k]['pull']))
    print('KLP18 published vs reproduced, max |diff| / typical published err:')
    for tag in kcmp:
        for k in KEYS6:
            c = kcmp[tag][k]
            print('   ', tag, k, f'{np.max(np.abs(np.array(c["published"]) - np.array(c["reproduced"]))):.2f}',
                  f'{np.median(c["published_err"]):.2f}')


if __name__ == '__main__':
    main()
