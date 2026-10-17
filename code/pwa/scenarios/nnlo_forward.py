"""Forward Compton amplitudes f(nu), g(nu) of the proton in NNLO BChPT (v2 code, main-fit configuration, static
polarisabilities at the code's own NNLO reference values BV.REF = driver Pref = PREF_NNLO: nothing fitted), for lab
photon energies nu = 50-450 MeV in 5-MeV steps.

Conventions (Gryniuk, Hagelstein, Pascalutsa, PRD 92 (2015) 074031 and PRD 94 (2016) 034043 = GHP):
    T(theta=0) = f(nu) eps'*.eps + i g(nu) sigma.(eps'* x eps),      dsigma/dOmega_lab(0) = |f|^2 + |g|^2,
    Sigma_2z(0) = -2 Re(f g*)/(|f|^2 + |g|^2)   (sigma_3/2 - sigma_1/2 convention),
    Im f = nu sigma_T/(4 pi),  Im g = -nu (sigma_3/2 - sigma_1/2)/(8 pi),
    f -> -alpha/M + (alpha_E1 + beta_M1) nu^2 + ...,   g -> -alpha kappa^2 nu/(2 M^2) + gamma_0 nu^3 + ...
Units: mub GeV (1 fm = 1973.27 mub GeV); the CSV also lists the amplitudes in fm.

Relation to the c.m. helicity amplitudes of compton.py (HMP16: Phi_k = T_{H'H}/(8 pi sqrt s), fm).  At theta = 0 only
the helicity-conserving Phi_1 = T_{-1/2,-1/2} and Phi_5 = T_{3/2,3/2} survive (d^J_{HH'}(0) = delta_{HH'}), and
    f = sqrt(s)/M (Phi_1 + Phi_5)/2,     g = sqrt(s)/M (Phi_1 - Phi_5)/2 ,
i.e. Phi_1 = (M/sqrt s)(f + g) (total helicity 1/2) and Phi_5 = (M/sqrt s)(f - g) (helicity 3/2); the factor
sqrt(s)/M is the lab/c.m. flux factor at theta = 0 (compton.lab_factor = s/M^2).  Established NUMERICALLY with the
exact Born + pi0-pole amplitude of compton.born_amplitudes: f_B = -alpha/M and g_B = -alpha kappa^2 nu/(2 M^2) are
reproduced to machine precision at all energies (sign +1 of both; the opposite sign would flip both).  With the
multipoles, this is the KLP18 relation f = sqrt(s)/M [EE1- + 2 EE1+ + MM1- + 2 MM1+ + 6 EE2- + 9 EE2+ + 6 MM2- + 9 MM2+
+ 18 EE3- + 18 MM3- + ...], g = sqrt(s)/M [EE1- - EE1+ - 6 EM1+ - 6 ME1+ + MM1- - MM1+ + 3 EE2- - 3 EE2+ + 3 MM2-
- 3 MM2+ - 24 EM2+ - 24 ME2+ + 6 EE3- + 6 MM3- + ...]  (pwa_model.forward_amps keeps the l <= 2 EE/MM terms).

Two evaluations of the same NNLO amplitude:
  (1) 'J52' (multipole route, as specified for this scenario):  f = f_Born(exact, closed formulas = compton.born_amplitudes
      at theta = 0) + f_nonBorn from the non-Born multipoles of bchpt_v2_multipoles.nonborn (all 14 multipoles the
      code prints, J <= 5/2, i.e. l <= 3 for the l- and l <= 2 for the l+ multipoles) through pw_from_multipoles ->
      phis_from_pw(theta=0) -> (f, g) as above.  The Fortran projection (pwa_gaga) is dimensioned for J <= 5/2
      (fmEE(3,2) etc.; jmax > 5 would overrun the arrays), so higher partial waves of the non-Born amplitude are
      NOT available: this route carries the truncation error of the non-Born partial-wave series, which grows with
      energy (dsigma(0): < 0.1 % below threshold, 0.5 % around the Delta peak, 1.7 % at 450 MeV; see validation (c)).
  (2) 'exact' (fit driver, no partial-wave truncation; the PRIMARY columns of the CSV): the driver returns every
      observable as an exact quadratic form in the shifts d of the six static polarisabilities.  In the forward
      direction the contact terms are exactly delta f = (d_alpha + d_beta) nu^2 and delta g = -(d_gEE + d_gMM + d_gEM
      + d_gME) nu^3 (lab nu; checked: the quadratic coefficients equal K nu^4 and K nu^6 to all digits, no mixing
      between f and g), so the linear coefficients of dsigma(0) give Re f and Re g, and |f|^2 + |g|^2 (= dsigma(0)/K)
      together with Re(f g*) (= -Sigma_2z dsigma(0)/2K) give Im f and Im g (Im f = nu sigma_T/4pi > 0 and
      |Im f| > |Im g|, which fixes the roots; Im f = Im g = 0 below the pi0 p threshold).
  Both are tabulated; (2) reproduces the driver's dsigma(0) and Sigma_2z(0) to 3e-6 by construction, (1) is the
  cross-check of the PWA machinery (and of the l-truncation).

Validation (numbers in work/notes/nnlo_forward_amplitudes.md, written by this script):
  (a) l-truncation of the Born part: forward amplitude from the Born-only multipoles (J <= 5/2) vs the exact Born;
  (b) low-energy limit of the NNLO non-Born f, g (5-50 MeV) vs the NNLO static values alpha+beta and gamma_0 of
      PREF_NNLO (work/fits/lib/common.py) = BV.REF;
  (c) dsigma(0) and Sigma_2z(0) of route (1) vs the fit driver (work/delta_plots/runner.py AND bchpt_v2_multipoles.driver,
      identical inputs to the multipole runs); and route (1) vs route (2) amplitude by amplitude.
Outputs: work/pwa/results/nnlo_forward_amplitudes.csv, work/pwa/results/nnlo_forward_validation.json,
         work/notes/nnlo_forward_amplitudes.md.
Run:  python3 work/pwa/scenarios/nnlo_forward.py
"""
import os
import sys
import json

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
WORK = os.path.dirname(PWA)
ROOT = os.path.dirname(WORK)
for p in (PWA, HERE, os.path.join(WORK, 'delta_plots')):
    if p not in sys.path:
        sys.path.insert(0, p)
os.chdir(PWA)

import numpy as np                                    # noqa: E402
import pandas as pd                                   # noqa: E402
import compton as C                                   # noqa: E402
import bchpt_v2_multipoles as BV                      # noqa: E402
import runner                                         # noqa: E402  (work/delta_plots: fit driver, cached)

CSV = os.path.join(PWA, 'results', 'nnlo_forward_amplitudes.csv')
VJSON = os.path.join(PWA, 'results', 'nnlo_forward_validation.json')
NOTE = os.path.join(WORK, 'notes', 'nnlo_forward_amplitudes.md')

EGRID = np.arange(50.0, 450.001, 5.0)
ELOW = np.arange(5.0, 50.001, 5.0)                    # low-energy limit (validation b)
UBGEV_PER_FM = 1973.27                                # 1 fm = 1973.27 mub GeV (sumrules.py)
K_NB = 1.0e7                                          # nb/sr per fm^2  (|f|^2 in fm^2 -> nb/sr)
K_UBGEV = 389379.0 / 389.379**2                       # 2.5682 nb/sr per (mub GeV)^2 (forward.py)
MFM = C.MN_GEV / C.HBARC                              # proton mass in fm^-1 (code constants)
E_PI0 = C.MPI0 * 1e3 * (1.0 + C.MPI0 / (2 * C.MN_GEV))  # pi0 p threshold, 144.7 MeV (code masses)
PKEYS = ['alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2']


def nu_fm(E):
    return np.asarray(E, float) / 1e3 / C.HBARC


def born_closed(E):
    """GHP Born forward amplitudes (fm): f_B = -alpha/M, g_B = -alpha kappa^2 nu/(2 M^2)."""
    E = np.asarray(E, float)
    return -C.ALPHA / MFM * np.ones_like(E), -C.ALPHA * C.KAPPA_P**2 * nu_fm(E) / (2 * MFM**2)


def fg_from_phis(P, E):
    """(f, g) [fm] from Phi_1..6 at theta = 0 (last axis), lab energy E [MeV]."""
    sq = np.sqrt(1.0 + 2.0 * np.asarray(E, float) * 1e-3 / C.MN_GEV)
    return sq * (P[..., 0] + P[..., 4]) / 2, sq * (P[..., 0] - P[..., 4]) / 2


def born_exact(E):
    """exact Born + pi0-pole helicity amplitudes of compton.py at theta = 0 -> (f_B, g_B) [fm]"""
    return fg_from_phis(C.phis_from_M(C.born_amplitudes(E, 0.0, 'lab')), E)


def fg_from_multipoles(f, E, Jmax2=5):
    """multipoles (dict key -> array over E, fm) -> (f, g) [fm] through the HMP16 partial-wave machinery."""
    P = C.phis_from_pw(C.pw_from_multipoles(f, Jmax2), 0.0)[..., 0, :]
    return fg_from_phis(P, E)


def zero_below_threshold(x, E):
    x = np.array(x, complex)
    x[np.asarray(E) < E_PI0] = x[np.asarray(E) < E_PI0].real
    return x


def multipole_route(E, pols=None):
    """route (1): exact Born + non-Born multipoles (J <= 5/2).  Returns dict of arrays (fm)."""
    nb, born = BV.nonborn(E, pols)
    keys = [k for k in nb if not str(k).startswith('_')]
    fnb, gnb = fg_from_multipoles({k: nb[k] for k in keys}, E)
    fBm, gBm = fg_from_multipoles({k: born[k] for k in keys}, E)           # Born-only multipoles (validation a)
    fB, gB = born_exact(E)
    fBc, gBc = born_closed(E)
    f_l3 = fg_from_multipoles({k: nb[k] for k in keys if k[1] == 3}, E)     # l = 3 (J = 5/2, '-') content
    f_em2 = fg_from_multipoles({k: nb[k] for k in keys if k[0] in ('EM', 'ME') and k[1] == 2}, E)
    return dict(E=E, f=zero_below_threshold(fB + fnb, E), g=zero_below_threshold(gB + gnb, E), fnb=fnb, gnb=gnb,
                fB=fB, gB=gB, fBc=fBc, gBc=gBc, fBm=fBm, gBm=gBm, f_l3=f_l3[0], g_l3=f_l3[1], g_em2=f_em2[1],
                ref=nb['_ref'], keys=keys)


def driver_route(E, name):
    """route (2): forward amplitudes from the quadratic forms of the fit driver (no partial-wave truncation).
    Returns dict with f, g (fm), the driver's own dsigma(0) [nb/sr], Sigma_2z(0), Pref and exactness checks."""
    E = np.asarray(E, float)
    n = len(E)
    pts = [(f'CS_{i}', e, 0.0, 'CS', 'LAB') for i, e in enumerate(E)] + [(f'S2Z_{i}', e, 0.0, 'S2Z', 'LAB') for i, e in enumerate(E)]
    out = runner.run(pts, name=name, verbose=False)
    assert np.allclose(out.E[:n], E) and np.allclose(out.E[n:], E) and all(out.obs[:n] == 'CS') and all(out.obs[n:] == 'S2Z')
    v = out.predict_array(out.Pref)
    cs, s2z = v[:n], v[n:]
    P = runner.PIDX
    CN, SN, SD = out.CN[:n], out.CN[n:], out.CD[n:]
    c00, c01, c02, c11, c22, c12 = (CN[:, P[k]] for k in ((0, 0), (0, 1), (0, 2), (1, 1), (2, 2), (1, 2)))
    c03, c04, c05, c06, c33, c44, c34, c13 = (CN[:, P[k]] for k in ((0, 3), (0, 4), (0, 5), (0, 6), (3, 3), (4, 4), (3, 4), (1, 3)))
    nu = nu_fm(E)
    u, w = np.sqrt(c11 / K_NB), np.sqrt(c33 / K_NB)            # |d f/d alpha|, |d g/d gamma| in fm per 1e-4 fm^3(4)
    checks = dict(u_over_nu2=float(np.max(abs(u / (nu**2 * 1e-4) - 1))), w_over_nu3=float(np.max(abs(w / (nu**3 * 1e-4) - 1))),
                  alpha_beta_sym=float(max(np.max(abs(c01 / c02 - 1)), np.max(abs(c11 / c22 - 1)), np.max(abs(c12 / (2 * c11) - 1)))),
                  spin_sym=float(max(np.max(abs(c03 / c04 - 1)), np.max(abs(c05 / c03 - 1)), np.max(abs(c06 / c03 - 1)),
                                     np.max(abs(c33 / c44 - 1)), np.max(abs(c34 / (2 * c33) - 1)))),
                  fg_mixing=float(np.max(abs(c13) / np.sqrt(c11 * c33))),
                  s2z_den_vs_cs=float(np.max(abs(SD - CN) / abs(c00[:, None]))))
    Ref = c01 / (2 * K_NB * u)                                   # dsigma = K |f + u d|^2 + ... -> coefficient of d
    Reg = -c03 / (2 * K_NB * w)                                  # delta g = -w d_gamma
    Refg = -0.5 * (SN[:, 0] / SD[:, 0]) * c00 / K_NB             # -2 K Re(f g*) = Sigma_2z * dsigma
    S = c00 / K_NB - Ref**2 - Reg**2                             # Im f^2 + Im g^2
    Pp = Refg - Ref * Reg                                        # Im f Im g
    disc = S**2 - 4 * Pp**2                                      # = (Im f^2 - Im g^2)^2 >= 0 (rounding: ~0 just above threshold)
    above = E >= E_PI0
    assert np.all(disc[above] > -1e-3 * S[above]**2), 'inconsistent imaginary parts'
    Imf = np.sqrt(np.maximum((S + np.sqrt(np.maximum(disc, 0))) / 2, 0))      # larger root: |Im f| > |Im g|
    Img = np.where(Imf > 0, Pp / np.where(Imf > 0, Imf, 1.0), 0.0)
    f = Ref + 1j * np.where(above, Imf, 0.0)
    g = Reg + 1j * np.where(above, Img, 0.0)
    rec = K_NB * (abs(f)**2 + abs(g)**2)
    checks['cs_reconstruction'] = float(np.max(abs(rec / cs - 1)))
    checks['s2z_reconstruction'] = float(np.max(abs(-2 * np.real(f * np.conj(g)) * K_NB / rec - s2z)))
    return dict(E=E, f=f, g=g, cs=cs, s2z=s2z, Pref=out.Pref, checks=checks)


def observables(f, g):
    """dsigma/dOmega_lab(0) [nb/sr] and Sigma_2z(0) from f, g in fm."""
    den = abs(f)**2 + abs(g)**2
    return K_NB * den, -2 * np.real(f * np.conj(g)) / den


def lowE_limit(E, fnb, gnb):
    """(alpha+beta)(nu) = fnb/nu^2 and gamma_0(nu) = gnb/nu^3 in 1e-4 fm^3, 1e-4 fm^4, and their nu -> 0 extrapolation
    (least squares a0 + a1 nu^2 + a2 nu^4 on the points outside the 12-30 MeV window, where the LoopTools evaluation of
    the loops loses precision -- the same window as multipoles_figure.smooth_lowE; the excluded points are reported)."""
    nu = nu_fm(E)
    ab = np.real(fnb) / nu**2 * 1e4
    g0 = np.real(gnb) / nu**3 * 1e4
    A = np.stack([np.ones_like(nu), nu**2, nu**4], -1)
    use = (E <= 12.0) | (E >= 30.0)
    cab = np.linalg.lstsq(A[use], ab[use], rcond=None)[0]
    cg0 = np.linalg.lstsq(A[use], g0[use], rcond=None)[0]
    return dict(E=E.tolist(), ab=ab.tolist(), g0=g0.tolist(), ab_extrap=float(cab[0]), g0_extrap=float(cg0[0]),
                ab_resid=float(np.max(abs(A @ cab - ab)[use])), g0_resid=float(np.max(abs(A @ cg0 - g0)[use])),
                ab_noise_12_30=float(np.max(abs(A @ cab - ab)[~use])), g0_noise_12_30=float(np.max(abs(A @ cg0 - g0)[~use])),
                fit_points=E[use].tolist())


def table():
    """DataFrame of the forward amplitudes on EGRID (reads the CSV if it exists, else computes and writes it)."""
    if os.path.exists(CSV):
        return pd.read_csv(CSV, comment='#')
    main()
    return pd.read_csv(CSV, comment='#')


def fmt_rows(rows, header):
    w = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    line = lambda r: '| ' + ' | '.join(str(x).rjust(w[i]) for i, x in enumerate(r)) + ' |'
    return '\n'.join([line(header), '|' + '|'.join('-' * (x + 2) for x in w) + '|'] + [line(r) for r in rows])


def main():
    E = EGRID
    M1 = multipole_route(E)
    D = driver_route(E, 'nnlo_forward')
    dr = BV.driver([(f'CS_{i}', e, 0.0, 'CS', 'LAB') for i, e in enumerate(E)] + [(f'S2Z_{i}', e, 0.0, 'S2Z', 'LAB') for i, e in enumerate(E)])
    csB = np.array([dr[f'CS_{i}'] for i in range(len(E))]); s2B = np.array([dr[f'S2Z_{i}'] for i in range(len(E))])
    try:
        sys.path.insert(0, os.path.join(WORK, 'fits', 'lib'))
        from common import PREF_NNLO
    except Exception as ex:                                       # noqa: BLE001
        print('WARNING: PREF_NNLO not imported (%s); using the driver Pref' % ex)
        PREF_NNLO = D['Pref']
    ref = np.array([BV.REF[k] for k in PKEYS])
    assert np.allclose(D['Pref'], ref, atol=2e-5) and np.allclose(PREF_NNLO, ref, atol=2e-5) and np.allclose(M1['ref'], ref, atol=2e-5), (D['Pref'], PREF_NNLO, M1['ref'], ref)
    ab_nnlo, g0_nnlo = ref[0] + ref[1], -np.sum(ref[2:])
    U = UBGEV_PER_FM

    # ---------- mapping / normalisation: exact Born at theta = 0 vs the GHP closed formulas
    map_dev = dict(f=float(np.max(abs(M1['fB'] / M1['fBc'] - 1))), g=float(np.max(abs(M1['gB'] / M1['gBc'] - 1))),
                   imag=float(max(np.max(abs(M1['fB'].imag)), np.max(abs(M1['gB'].imag)))))
    # ---------- (a) l-truncation of the Born part (Born-only multipoles, J <= 5/2, vs exact Born)
    aF = (M1['fBm'].real / M1['fBc'] - 1); aG = (M1['gBm'].real / M1['gBc'] - 1)
    # ---------- (b) low-energy limit
    L1 = multipole_route(ELOW)
    LD = driver_route(ELOW, 'nnlo_forward_lowE')
    fBl, gBl = born_closed(ELOW)
    lowE = dict(J52=lowE_limit(ELOW, L1['fnb'], L1['gnb']), exact=lowE_limit(ELOW, LD['f'] - fBl, LD['g'] - gBl),
                nnlo=dict(alpha_p_beta=float(ab_nnlo), gamma0=float(g0_nnlo)))
    # ---------- (c) observables of route (1) vs the driver, and route (1) vs route (2)
    cs1, s21 = observables(M1['f'], M1['g'])
    cs2, s22 = observables(D['f'], D['g'])
    rel_cs = cs1 / D['cs'] - 1; d_s2 = s21 - D['s2z']
    sel = {'<=150': E <= 150, '150-300': (E > 150) & (E <= 300), '300-400': (E > 300) & (E <= 400), '400-450': E > 400}
    c_summary = {k: dict(max_rel_dcs=float(np.max(abs(rel_cs[m]))), max_abs_dS2z=float(np.max(abs(d_s2[m]))),
                         max_df_ubgev=float(np.max(abs(M1['f'] - D['f'])[m]) * U), max_dg_ubgev=float(np.max(abs(M1['g'] - D['g'])[m]) * U))
                 for k, m in sel.items()}
    E_1pct = float(E[np.where(abs(rel_cs) > 0.01)[0][0]]) if np.any(abs(rel_cs) > 0.01) else None
    V = dict(mapping=dict(statement='f = sqrt(s)/M (Phi1 + Phi5)/2, g = sqrt(s)/M (Phi1 - Phi5)/2 [fm]; Phi_k of compton.py at theta = 0',
                          max_rel_dev_vs_closed_born=map_dev),
             a_born_truncation=dict(E=E.tolist(), rel_dev_f=aF.tolist(), rel_dev_g=aG.tolist(),
                                    note='forward f, g from the 14 Born-only multipoles (J <= 5/2) over the exact Born, minus 1'),
             b_lowE=lowE,
             c_vs_driver=dict(E=E.tolist(), rel_dcs_J52=rel_cs.tolist(), dS2z_J52=d_s2.tolist(), summary=c_summary,
                              first_E_rel_dcs_above_1pct=E_1pct,
                              runner_vs_BVdriver=dict(max_rel_dcs=float(np.max(abs(D['cs'] / csB - 1))), max_abs_dS2z=float(np.max(abs(D['s2z'] - s2B)))),
                              exact_route_checks=D['checks'],
                              l3_content_ubgev=dict(f=float(np.max(abs(M1['f_l3'])) * U), g=float(np.max(abs(M1['g_l3'])) * U)),
                              em2_me2_in_g_ubgev=dict(max_all=float(np.max(abs(M1['g_em2'])) * U), at_150=float(abs(M1['g_em2'][E == 150][0]) * U))),
             pref=dict(driver=D['Pref'].tolist(), bv_ref=ref.tolist(), pref_nnlo=np.asarray(PREF_NNLO).tolist()))
    json.dump(V, open(VJSON, 'w'), indent=1)

    # ---------- table
    T = pd.DataFrame(dict(E_MeV=E,
                          Re_f=D['f'].real * U, Im_f=D['f'].imag * U, Re_g=D['g'].real * U, Im_g=D['g'].imag * U,
                          dsigma0_nb=cs2, Sigma2z0=s22,
                          Re_f_J52=M1['f'].real * U, Im_f_J52=M1['f'].imag * U, Re_g_J52=M1['g'].real * U, Im_g_J52=M1['g'].imag * U,
                          dsigma0_nb_J52=cs1, Sigma2z0_J52=s21,
                          dsigma0_nb_driver=D['cs'], Sigma2z0_driver=D['s2z'],
                          Re_f_fm=D['f'].real, Im_f_fm=D['f'].imag, Re_g_fm=D['g'].real, Im_g_fm=D['g'].imag,
                          f_Born=M1['fBc'] * U, g_Born=M1['gBc'] * U))
    with open(CSV, 'w') as fh:
        fh.write('# NNLO BChPT forward Compton amplitudes of the proton (v2 code, main-fit configuration, static polarisabilities = NNLO\n'
                 '# reference values %s; nothing fitted).  GHP conventions, mub GeV (1 fm = 1973.27 mub GeV).\n'
                 '# Re_f..Im_g, dsigma0_nb, Sigma2z0: exact forward amplitudes from the quadratic forms of the fit driver (no partial-wave\n'
                 '# truncation).  *_J52: exact Born + non-Born multipoles J <= 5/2 through the PWA machinery (l-truncated above ~300 MeV).\n'
                 '# *_driver: dsigma(0) [nb/sr] and Sigma_2z(0) as returned by the driver.  *_fm: primary amplitudes in fm.  f_Born, g_Born:\n'
                 '# -alpha/M, -alpha kappa^2 nu/(2 M^2) [mub GeV].  Im parts vanish below the pi0 p threshold (144.7 MeV).\n'
                 '# Written by work/pwa/scenarios/nnlo_forward.py\n' % ', '.join(f'{k}={v:.5f}' for k, v in zip(PKEYS, ref)))
        T.to_csv(fh, index=False, float_format='%.7g')

    # ---------- note
    i_sel = [int(np.argmin(abs(E - e))) for e in (50, 100, 150, 200, 250, 300, 330, 350, 400, 450)]
    rows_c = [[f'{E[i]:.0f}', f'{cs1[i]:.4f}', f'{D["cs"][i]:.4f}', f'{100 * rel_cs[i]:+.3f}', f'{s21[i]:+.4f}', f'{D["s2z"][i]:+.4f}', f'{d_s2[i]:+.5f}'] for i in i_sel]
    rows_fg = [[f'{E[i]:.0f}', f'{M1["f"][i].real * U:+.4f}', f'{D["f"][i].real * U:+.4f}', f'{M1["f"][i].imag * U:+.4f}', f'{D["f"][i].imag * U:+.4f}',
                f'{M1["g"][i].real * U:+.4f}', f'{D["g"][i].real * U:+.4f}', f'{M1["g"][i].imag * U:+.4f}', f'{D["g"][i].imag * U:+.4f}'] for i in i_sel]
    rows_a = [[f'{E[i]:.0f}', f'{100 * aF[i]:+.3f}', f'{100 * aG[i]:+.3f}'] for i in i_sel]
    rows_b = [[f'{ELOW[i]:.0f}', f'{lowE["J52"]["ab"][i]:.4f}', f'{lowE["exact"]["ab"][i]:.4f}', f'{lowE["J52"]["g0"][i]:+.4f}', f'{lowE["exact"]["g0"][i]:+.4f}'] for i in range(len(ELOW))]
    txt = f"""# NNLO BChPT forward amplitudes f(nu), g(nu) of the proton, 50-450 MeV: validation

Script: `work/pwa/scenarios/nnlo_forward.py` (run with the project venv).  Outputs: `work/pwa/results/nnlo_forward_amplitudes.csv`
(table, mub GeV and fm), `work/pwa/results/nnlo_forward_validation.json` (all validation arrays), this note.
Amplitude: v2 code `work/code_v2/execs/ampl.exe`, main-fit configuration (dvertex=t, pimassphys=t, nmassphys=t, consts='LMP15',
g_M = 2.97, g_E = -1), static polarisabilities at the code's NNLO reference values (driver Pref = `bchpt_v2_multipoles.REF` =
`PREF_NNLO`): alpha_E1, beta_M1, gamma_E1E1, gamma_M1M1, gamma_E1M2, gamma_M1E2 = {', '.join(f'{x:.5f}' for x in ref)}
(alpha+beta = {ab_nnlo:.4f}, gamma_0 = {g0_nnlo:.4f}); nothing fitted.  GHP conventions (PRD 92, 074031; PRD 94, 034043):
T(0) = f eps'*.eps + i g sigma.(eps'* x eps), dsigma/dOmega_lab(0) = |f|^2 + |g|^2, Sigma_2z(0) = -2 Re(f g*)/(|f|^2+|g|^2)
(sigma_3/2 - sigma_1/2 convention), Im f = nu sigma_T/(4 pi), Im g = -nu Delta sigma/(8 pi).  Units mub GeV, 1 fm = 1973.27 mub GeV.

## Mapping of the helicity amplitudes onto f, g

At theta = 0 only Phi_1 = T_{{-1/2,-1/2}}/(8 pi sqrt s) and Phi_5 = T_{{3/2,3/2}}/(8 pi sqrt s) of `compton.py` survive
(d^J_{{HH'}}(0) = delta_{{HH'}}).  Numerically, with the exact Born + pi0-pole amplitude `compton.born_amplitudes(E, 0, 'lab')`:

    f = sqrt(s)/M (Phi_1 + Phi_5)/2,    g = sqrt(s)/M (Phi_1 - Phi_5)/2        [fm]

reproduces f_B = -alpha/M and g_B = -alpha kappa^2 nu/(2 M^2) at all energies 50-450 MeV: max |f/f_B - 1| = {map_dev['f']:.1e},
max |g/g_B - 1| = {map_dev['g']:.1e}, imaginary parts {map_dev['imag']:.1e} fm.  (Phi_1 = (M/sqrt s)(f+g) is the total-helicity-1/2 and
Phi_5 = (M/sqrt s)(f-g) the helicity-3/2 amplitude; the lab/c.m. factor sqrt(s)/M is the square root of `compton.lab_factor` at
theta = 0; the pi0 pole vanishes at t = 0.)  With the multipoles inserted, the same chain gives the KLP18 relations of
`pwa_model.forward_amps` (l <= 2 EE/MM terms) exactly, plus the terms -24 (f_EM^2+ + f_ME^2+) in g and the l = 3 terms
18 (f_EE^3- + f_MM^3-) in f, 6 (f_EE^3- + f_MM^3-) in g.

## Two evaluations

1. **J52** (multipole route, as specified): exact Born (closed formulas) + non-Born multipoles of `bchpt_v2_multipoles.nonborn`
   (all 14 multipoles printed by the code: J <= 5/2) -> `pw_from_multipoles` -> `phis_from_pw(theta=0)` -> (f, g).  The Fortran
   projection `pwa_gaga` is dimensioned for J <= 5/2 (`fmEE(3,2)`, `fpEE(3,2)`, `fpEM(3,2)`; `jmax > 5` would overrun the arrays), so
   this route truncates the non-Born partial-wave series at J = 5/2.
2. **exact** (primary columns of the CSV): the fit driver returns dsigma(0) and Sigma_2z(0) as exact quadratic forms in the shifts d of
   the six static polarisabilities.  In the forward direction the O(p^4) contact terms are exactly delta f = (d_alpha + d_beta) nu^2 and
   delta g = -(d_gEE + d_gMM + d_gEM + d_gME) nu^3 with the lab nu: the quadratic coefficients equal K nu^4 and K nu^6 to
   {max(D['checks']['u_over_nu2'], D['checks']['w_over_nu3']):.0e} (relative), alpha/beta and the four gamma's enter symmetrically to
   {max(D['checks']['alpha_beta_sym'], D['checks']['spin_sym']):.0e}, and there is no f-g mixing ({D['checks']['fg_mixing']:.0e}).  Hence the linear
   coefficients give Re f and Re g exactly, and |f|^2 + |g|^2 = dsigma(0)/K with Re(f g*) = -Sigma_2z dsigma(0)/(2K) give Im f, Im g
   (Im f > 0 and |Im f| > |Im g| since sigma_T > |Delta sigma|/2 fix the roots; Im f = Im g = 0 below the pi0 p threshold, {E_PI0:.1f} MeV).
   Reconstruction of the driver's dsigma(0) and Sigma_2z(0) from these f, g: {D['checks']['cs_reconstruction']:.1e} (relative) and
   {D['checks']['s2z_reconstruction']:.1e} (absolute).

## (a) l-truncation of the Born part

Forward f, g from the 14 Born-only multipoles (J <= 5/2, same projection) relative to the exact Born, in per cent
(this is why the Born part is taken exactly; the non-Born part is better convergent, see (c)):

{fmt_rows(rows_a, ['E [MeV]', 'f: 100(f_mult/f_B - 1)', 'g: 100(g_mult/g_B - 1)'])}

With J <= 3/2 only (l <= 2 in f, no l = 3), the Born g is off by 8-9 % and f by up to 4 % at 150-450 MeV.

## (b) Low-energy limit

Non-Born f/nu^2 and g/nu^3 (10^-4 fm^3, 10^-4 fm^4) of both routes at 5-50 MeV; NNLO static values alpha+beta = {ab_nnlo:.4f},
gamma_0 = {g0_nnlo:.4f}.  Extrapolation nu -> 0 (least squares a0 + a1 nu^2 + a2 nu^4 on the points outside 12-30 MeV):
alpha+beta -> {lowE['J52']['ab_extrap']:.4f} (J52) / {lowE['exact']['ab_extrap']:.4f} (exact), gamma_0 -> {lowE['J52']['g0_extrap']:.4f} (J52) /
{lowE['exact']['g0_extrap']:.4f} (exact); max fit residuals {lowE['J52']['ab_resid']:.1e}, {lowE['J52']['g0_resid']:.1e} (J52).  Direct values at
5 MeV: {lowE['J52']['ab'][0]:.4f} and {lowE['J52']['g0'][0]:.4f}, i.e. {lowE['J52']['ab'][0] - ab_nnlo:+.4f} and {lowE['J52']['g0'][0] - g0_nnlo:+.4f} from the static values (the size of
the nu^2 terms of the sum rules, c_4 nu^2 and gamma_0-bar nu^2, at 5 MeV).

{fmt_rows(rows_b, ['E [MeV]', '(a+b)(nu) J52', '(a+b)(nu) exact', 'g0(nu) J52', 'g0(nu) exact'])}

(The rise of (alpha+beta)(nu) is the nu^2 term of the sum rule, c_4 ~ 6 x 10^-4 fm^5 in GHP; the 15-25 MeV values of gamma_0(nu) carry
the LoopTools cancellation noise noted in `work/notes/multipoles_figure.md`: up to {lowE['J52']['g0_noise_12_30']:.2f} off the smooth fit there
(J52; {lowE['exact']['g0_noise_12_30']:.2f} exact), {lowE['J52']['ab_noise_12_30']:.3f} in (alpha+beta)(nu); the fits exclude this window.)

## (c) Route (1) vs the fit driver at theta = 0 (dsigma/dOmega_lab in nb/sr)

`runner.run` (work/delta_plots; template with M_N = 0.9383 GeV in baryon.dat) and `bchpt_v2_multipoles.driver` (identical inputs to the
multipole runs) agree with each other to {V['c_vs_driver']['runner_vs_BVdriver']['max_rel_dcs']:.1e} (dsigma) and
{V['c_vs_driver']['runner_vs_BVdriver']['max_abs_dS2z']:.1e} (Sigma_2z); both have Pref = REF.

{fmt_rows(rows_c, ['E [MeV]', 'dsigma J52', 'dsigma driver', 'dev. [%]', 'S2z J52', 'S2z driver', 'diff'])}

Maximal deviations of route (1) from the driver: {'; '.join(f"{k} MeV: dsigma {100 * v['max_rel_dcs']:.2f} %, Sigma_2z {v['max_abs_dS2z']:.4f}" for k, v in c_summary.items())}.
The 1 % level is exceeded from E = {E_1pct} MeV on.  **This is the J <= 5/2 truncation of the non-Born partial-wave series, not an error of
the mapping**: the exact route (2) reproduces the driver to 1e-6, and the difference between the two routes, amplitude by amplitude, is

{fmt_rows(rows_fg, ['E [MeV]', 'Re f J52', 'Re f exact', 'Im f J52', 'Im f exact', 'Re g J52', 'Re g exact', 'Im g J52', 'Im g exact'])}

(max |f_J52 - f_exact| = {', '.join(f"{v['max_df_ubgev']:.3f} ({k})" for k, v in c_summary.items())} mub GeV; for g
{', '.join(f"{v['max_dg_ubgev']:.3f} ({k})" for k, v in c_summary.items())} mub GeV.)  The l = 3 (J = 5/2, '-') multipoles alone contribute up to
{V['c_vs_driver']['l3_content_ubgev']['f']:.3f} mub GeV to f and {V['c_vs_driver']['l3_content_ubgev']['g']:.3f} mub GeV to g (at 450 MeV); the J >= 7/2 remainder is of
the same size as the J = 5/2 term there (geometric convergence with ratio ~0.5), i.e. the series is not exhausted at J = 5/2 above ~300 MeV.  Below the pion threshold the two routes agree to better than
{1e3 * c_summary['<=150']['max_df_ubgev']:.1f} x 10^-3 mub GeV (f) and {1e3 * c_summary['<=150']['max_dg_ubgev']:.1f} x 10^-3 mub GeV (g).
Side remark for the PWA: the mixed l = 2 multipoles f_EM^2+, f_ME^2+, which the KLP18 parametrisation does not contain, enter g with the
coefficient -24 sqrt(s)/M; in NNLO BChPT this term is {1e3 * V['c_vs_driver']['em2_me2_in_g_ubgev']['at_150']:.2f} x 10^-3 mub GeV at 150 MeV
(GHP error of Re g there: 0.023 mub GeV), negligible for the fits.

## Decision

The CSV's primary columns (Re_f, Im_f, Re_g, Im_g, dsigma0_nb, Sigma2z0) and the figure `paper/figs/pwa_forward.pdf` use route (2),
the exact forward amplitudes of the same NNLO prediction; route (1) is kept in the columns `*_J52` as the check of the PWA machinery
(and would be the only choice if only the multipoles were available).  Below 300 MeV the two are indistinguishable on the figure.
"""
    open(NOTE, 'w').write(txt)
    print(txt)
    print('written:', CSV, VJSON, NOTE)
    return T, V


if __name__ == '__main__':
    main()
