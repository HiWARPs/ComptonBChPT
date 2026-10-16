"""Validation of the PWA machinery against the Fortran BChPT code (step 1 of the task).
Writes validation/validation_log.txt"""
import numpy as np
import compton as C
import fortran_runner as FR

out = open('validation/validation_log.txt', 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    out.write(s + '\n')


HBARC_MEV = C.HBARC * 1e3
MN = C.MN_GEV * 1e3
NG = 64
xg, wg = np.polynomial.legendre.leggauss(NG)
thg = np.degrees(np.arccos(xg))


def mult_of(Mgrid, Jmax2=31):
    P = C.phis_from_M(Mgrid)
    pw = C.pw_project(P, xg, wg, Jmax2)
    return C.multipoles_from_pw(pw), pw


# ---------- 1. Python Born / pi0 vs Fortran ----------
log('=== 1. exact Born and pi0-pole amplitudes: Python (Feynman rules) vs Fortran (A_i route)')
Es = np.array([10., 59., 89., 109., 135., 150.])
for E in Es:
    Mf, _ = FR.run(E, thg, 'CM', 'bchpt')
    Mb = C.born_amplitudes(E, thg, 'cm', pi0=False)
    Mp = C.born_amplitudes(E, thg, 'cm', born=False)
    log(f'  E={E:6.1f}  max rel dev Born {np.max(abs(Mb-Mf[:,0]))/np.max(abs(Mf[:,0])):.1e}'
        f'   pi0 {np.max(abs(Mp-Mf[:,1]))/np.max(abs(Mf[:,1])):.1e}')

# ---------- 2. low-energy theorem ----------
log('=== 2. LET: Born multipoles at E_lab=0.05 MeV (units alpha/M); expect fEE1-=fEE1+=-1/3, others 0')
E = 0.05
Mb = C.born_amplitudes(E, thg, 'cm', pi0=True)
f, pw = mult_of(Mb, 9)
aM = C.ALPHA * C.FM_PER_M
for key in [('EE', 1, '-'), ('EE', 1, '+'), ('MM', 1, '-'), ('MM', 1, '+'), ('EM', 1, '+'), ('ME', 1, '+'),
            ('EE', 2, '-'), ('EE', 2, '+'), ('MM', 2, '-'), ('MM', 2, '+')]:
    log(f'  f{key[0]}{key[1]}{key[2]} = {np.real(f[key])/aM:+.6f}')

# ---------- 3. resummation of partial waves ----------
log('=== 3. partial-wave resummation (theta_cm = 20..160 deg)')
log('   NOTE: the NNLO loop amplitudes of ampl.f90 are numerically unstable within ~2-3 deg of theta_cm=0 (piN loops)')
log('   and theta_cm=180 deg (piN and piDelta loops); the non-Born part is therefore projected with an')
log('   endpoint-avoiding 20-point Gauss-Legendre rule (extreme nodes at 7.2/172.8 deg).')
xg20, wg20 = np.polynomial.legendre.leggauss(20)
thg20 = np.degrees(np.arccos(xg20))
thc = np.array([20., 45., 70., 90., 110., 135., 150., 160.])
for E in [59., 89., 109., 135., 150.]:
    k = C.kin(E, thc, 'cm')
    Md, _ = FR.run(E, thc, 'CM', 'bchpt')
    od = C.obs_from_M(Md.sum(axis=1), k)
    PB = C.phis_from_M(C.born_amplitudes(E, thc, 'cm'))
    # (a) Born+pi0 alone, resummed from partial waves (64-pt grid, exact amplitude)
    Pg = C.phis_from_M(C.born_amplitudes(E, thg, 'cm'))
    oB = C.obs_from_phis(PB, k)
    msg = []
    for Jm2 in (9, 15, 21):
        Pr = C.phis_from_pw(C.pw_project(Pg, xg, wg, Jm2), np.radians(thc))
        o = C.obs_from_phis(Pr, k)
        msg.append(f'J<={Jm2}/2: {np.max(abs(o["cs_cm"]/oB["cs_cm"]-1)):.1e}')
    log(f'  E={E:5.1f} Born+pi0 resummed vs exact, max rel dcs: ' + ', '.join(msg))
    # (b) full = exact Born+pi0 + projected non-Born, vs direct Fortran
    Mf, _ = FR.run(E, thg20, 'CM', 'bchpt')
    Pnb = C.phis_from_M(Mf[:, 2:].sum(axis=1))
    msg = []
    for Jm2 in (3, 5, 9, 13):
        Pr = C.phis_from_pw(C.pw_project(Pnb, xg20, wg20, Jm2), np.radians(thc)) + PB
        o = C.obs_from_phis(Pr, k)
        msg.append(f'J<={Jm2}/2: {np.max(abs(o["cs_cm"]/od["cs_cm"]-1)):.1e}/{np.max(abs(o["S3"]-od["S3"])):.1e}')
    log(f'  E={E:5.1f} Born exact + non-Born(J<=Jmax) vs direct BChPT, max rel dcs / max |dS3|: ' + ', '.join(msg))
    # (c) PLB truncation: ten multipoles (l=1 and l=2 EE,MM)
    f, _ = (C.multipoles_from_pw(C.pw_project(Pnb, xg20, wg20, 9)), None)
    keep = [('EE', 1, '-'), ('EE', 1, '+'), ('MM', 1, '-'), ('MM', 1, '+'), ('EM', 1, '+'), ('ME', 1, '+'),
            ('EE', 2, '-'), ('EE', 2, '+'), ('MM', 2, '-'), ('MM', 2, '+')]
    f10 = {kk: f[kk] for kk in keep}
    Pr = C.phis_from_pw(C.pw_from_multipoles(f10, 5), np.radians(thc)) + PB
    o = C.obs_from_phis(Pr, k)
    f6 = {kk: f[kk] for kk in keep[:6]}
    Pr6 = C.phis_from_pw(C.pw_from_multipoles(f6, 3), np.radians(thc)) + PB
    o6 = C.obs_from_phis(Pr6, k)
    log(f'  E={E:5.1f} PLB truncation (10 multipoles) vs direct: max rel dcs {np.max(abs(o["cs_cm"]/od["cs_cm"]-1)):.1e}, '
        f'max |dS3| {np.max(abs(o["S3"]-od["S3"])):.1e};  l=1 only: {np.max(abs(o6["cs_cm"]/od["cs_cm"]-1)):.1e} / '
        f'{np.max(abs(o6["S3"]-od["S3"])):.1e}')
for E in [109.]:
    f, pw = mult_of(C.born_amplitudes(E, thg, 'cm'), 31)
    pw2 = C.pw_from_multipoles(f, 31)
    dev = max(np.max(abs(pw2[J]-pw[J])) for J in pw)
    log(f'  multipole <-> partial-wave round trip (Born, J<=31/2) max dev {dev:.1e} fm')

# ---------- 4. observables vs Fortran's own observables (run_points / SUBROUTINE observables) ----------
log('=== 4. observables: Python formulas vs Fortran fit-driver (pointsmode) at lab kinematics')
import os, subprocess, tempfile
pts = [(59., 60.), (89., 90.), (109., 133.), (135., 155.), (92.2, 35.), (129.5, 145.), (150., 120.)]
txt = open(os.path.join(FR.BASE, 'fortran', 'runs', 'bchpt', 'in.dat')).read().replace(
    "helmode=t, pointsfile='pts.txt', outfile='hel.out', pointsmode=f", "helmode=f, pointsfile='pp.txt', outfile='pp.out', pointsmode=t")
with tempfile.TemporaryDirectory() as td:
    open(os.path.join(td, 'in.dat'), 'w').write(txt)
    with open(os.path.join(td, 'pp.txt'), 'w') as fh:
        for i, (e, t) in enumerate(pts):
            for ob in ('CS', 'S3', 'SPAR', 'SPERP', 'CSCM'):
                fh.write(f'p{i} {e} {t} {ob} LAB\n')
    os.symlink(FR.BARYON, os.path.join(td, 'baryon.dat'))
    subprocess.run([FR.EXE], cwd=td, capture_output=True)
    rows = [l.split() for l in open(os.path.join(td, 'pp.out')) if not l.startswith('#')]
fort = {(r[0], r[1]): float(r[5]) for r in rows}
for i, (e, t) in enumerate(pts):
    Md, _ = FR.run(e, t, 'LAB', 'bchpt')
    k = C.kin(e, t, 'lab')
    o = C.obs_from_M(Md.sum(axis=1), k)
    log(f'  E={e:6.1f} th_lab={t:5.1f}: CS_lab {o["cs_lab"][0]:8.4f} vs {fort[("p%d"%i,"CS")]:8.4f} | '
        f'CS_cm {o["cs_cm"][0]:8.4f} vs {fort[("p%d"%i,"CSCM")]:8.4f} | S3 {o["S3"][0]:+.5f} vs {fort[("p%d"%i,"S3")]:+.5f} | '
        f'par {o["cs_par_lab"][0]:8.4f} vs {fort[("p%d"%i,"SPAR")]:8.4f} | perp {o["cs_perp_lab"][0]:8.4f} vs {fort[("p%d"%i,"SPERP")]:8.4f}')

# ---------- 5. LEX: non-Born multipoles of the contact (polarisability) amplitude ----------
log('=== 5. LEX check of the PLB parametrisation (Eq. nonBorn_multipoles_expansion), R=(E/M)^2 f^R residual')
log('   contact (LEC) amplitude of ampl.f90 for a shift of one static polarisability at a time')
log('   [NB: onlyLEX=t does NOT work in this code version: the NNLO-completeness check that follows the onlyLEX')
log('    block zeroes all LEC shifts, so onlyLEX gives pure Born+pi0. We use the complete NNLO run with shifted')
log('    polarisabilities instead; component 6 (contact) then carries exactly the shift.]')
log('   printed: f^R_i(E) = [fbar/(E^k M/sqrt s) - (static + recoil terms)] (M/E)^2 -> should stay O(1)')


def plb_terms(E_fm, pol):
    """bracket terms of the PLB parametrisation without residuals; pol dict in 1e-4 units -> fm units."""
    a = pol.get('a', 0) * 1e-4; b = pol.get('b', 0) * 1e-4
    gEE = pol.get('gEE', 0) * 1e-4; gMM = pol.get('gMM', 0) * 1e-4
    gEM = pol.get('gEM', 0) * 1e-4; gME = pol.get('gME', 0) * 1e-4
    Mfm = MN / HBARC_MEV
    E = E_fm
    return {('EE', 1, '+'): a / 3 + E / 3 * ((-a + b) / Mfm + gEE),
            ('EE', 1, '-'): a / 3 + E / 3 * ((-a + b) / Mfm - 2 * gEE),
            ('MM', 1, '+'): b / 3 + E / 3 * ((-b + a) / Mfm + gMM),
            ('MM', 1, '-'): b / 3 + E / 3 * ((-b + a) / Mfm - 2 * gMM),
            ('EM', 1, '+'): gEM / 6 + E / 6 * ((-6 * gEM + 3 * gME + 3 * gMM) / (4 * Mfm) - b / (8 * Mfm**2)),
            ('ME', 1, '+'): gME / 6 + E / 6 * ((-6 * gME + 3 * gEM + 3 * gEE) / (4 * Mfm) - a / (8 * Mfm**2))}


Mfm = MN / HBARC_MEV
cases = [('a', dict(alphaE1fm=10.0, betaM1fm=0, gamEEfm=0, gamMMfm=0, gamEMfm=0, gamMEfm=0), {'a': 10}),
         ('b', dict(alphaE1fm=0, betaM1fm=10.0, gamEEfm=0, gamMMfm=0, gamEMfm=0, gamMEfm=0), {'b': 10}),
         ('gEE', dict(alphaE1fm=0, betaM1fm=0, gamEEfm=5.0, gamMMfm=0, gamEMfm=0, gamMEfm=0), {'gEE': 5}),
         ('gMM', dict(alphaE1fm=0, betaM1fm=0, gamEEfm=0, gamMMfm=5.0, gamEMfm=0, gamMEfm=0), {'gMM': 5}),
         ('gEM', dict(alphaE1fm=0, betaM1fm=0, gamEEfm=0, gamMMfm=0, gamEMfm=5.0, gamMEfm=0), {'gEM': 5}),
         ('gME', dict(alphaE1fm=0, betaM1fm=0, gamEEfm=0, gamMMfm=0, gamEMfm=0, gamMEfm=5.0), {'gME': 5})]
keys = [('EE', 1, '+'), ('EE', 1, '-'), ('MM', 1, '+'), ('MM', 1, '-'), ('EM', 1, '+'), ('ME', 1, '+')]
for name, pols, pd in cases:
    log(f'  -- only {name} = {list(pd.values())[0]}')
    for E in (2., 10., 30., 60., 100., 150.):
        pp = {kk: FR.NNLO[kk] + v for kk, v in pols.items()}
        Mf, _ = FR.run(E, thg20, 'CM', 'bchpt', pp)
        f = C.multipoles_from_pw(C.pw_project(C.phis_from_M(Mf[:, 5]), xg20, wg20, 9))
        k = C.kin(E, 90., 'cm')
        Efm = E / HBARC_MEV
        pref = (1.0 / np.sqrt(k['s']))   # M/sqrt(s)
        terms = plb_terms(Efm, pd)
        row = []
        for key in keys:
            n = 3 if key[0] in ('EM', 'ME') else 2
            br = np.real(f[key]) / (Efm**n * pref)
            fR = (br - terms[key]) / (Efm / Mfm)**2 / 1e-4
            row.append(f'{fR:+9.3f}')
        log(f'     E={E:6.1f}  fR(EE+,EE-,MM+,MM-,EM,ME) [1e-4 fm^3|fm^4] = ' + ' '.join(row))

# ---------- 6. forward amplitudes: multipole sum rules vs direct; sign of g ----------
log('=== 6. forward amplitudes f, g: direct (theta=0) vs multipole sums (PLB Eq. A_1_A_3_multipoles)')
for E in (20., 59., 109., 150.):
    k0 = C.kin(E, 1e-4, 'cm')
    # direct: T = f eps'*.eps + g i (eps'* x eps).sigma  ->  Phi1 = (M/sqrt s)(f+g)?, Phi5 = (M/sqrt s)(f-g)?
    Mb = C.born_amplitudes(E, [1e-4], 'cm')
    xg20, wg20 = np.polynomial.legendre.leggauss(20)
    thg20 = np.degrees(np.arccos(xg20))
    Pnb20 = C.phis_from_M(FR.run(E, thg20, 'CM', 'bchpt')[0][:, 2:].sum(axis=1))
    PBf = C.phis_from_M(Mb[0])
    Pfull0 = PBf + C.phis_from_pw(C.pw_project(Pnb20, xg20, wg20, 9), np.array([1e-6]))[0]
    for lab, P in (('Born+pi0', PBf), ('BChPT', Pfull0)):
        sq = np.sqrt(k0['s'])
        fdir = sq * (P[0] + P[4]) / 2
        gdir = sq * (P[0] - P[4]) / 2
        if lab == 'Born+pi0':
            Mg = C.born_amplitudes(E, thg, 'cm')
        else:
            Mg = None
        if Mg is not None:
            f, _ = mult_of(Mg, 41)
        else:
            fB, _ = mult_of(C.born_amplitudes(E, thg, 'cm'), 41)
            fN = C.multipoles_from_pw(C.pw_project(Pnb20, xg20, wg20, 9))
            f = {kk: fB.get(kk, 0) + fN.get(kk, 0) for kk in set(fB) | set(fN)}
        fs = 0; gs = 0
        for L in range(0, 20):
            fs += (L + 1)**2 * ((L + 2) * (f.get(('EE', L + 1, '-'), 0) + f.get(('MM', L + 1, '-'), 0))
                                + L * (f.get(('EE', L, '+'), 0) + f.get(('MM', L, '+'), 0)))
            gs += (L + 1) * ((L + 2) * (f.get(('EE', L + 1, '-'), 0) + f.get(('MM', L + 1, '-'), 0))
                             - L * (f.get(('EE', L, '+'), 0) + f.get(('MM', L, '+'), 0))
                             - 2 * L * (L + 2) * (f.get(('EM', L, '+'), 0) + f.get(('ME', L, '+'), 0)))
        fs *= sq / 2; gs *= sq / 2
        nu = E / HBARC_MEV
        gB = -C.ALPHA * C.KAPPA_P**2 * nu / (2 * Mfm**2)
        log(f'  E={E:6.1f} {lab:9s}: f_dir={np.real(fdir):+.6e} f_mult={np.real(fs):+.6e} | g_dir={np.real(gdir):+.6e} '
            f'g_mult={np.real(gs):+.6e} | -alpha kappa^2 nu/2M^2={gB:+.6e}  -alpha/M={-C.ALPHA/Mfm:+.6e} [fm]')
out.close()
