"""Validation of the double-polarisation asymmetries Sigma_2x, Sigma_2z of compton.py (S2_from_M, S2_from_phis; keys
S2X, S2Z of obs_from_M / obs_from_phis, added 2026-10-17) against the Fortran fit driver (work/code_v2/execs/ampl.exe,
observables S2X, S2Z, LAB frame = the conventions of the A2 data, driver README Sect. 3).

  (1) Born + pi0-pole amplitude: compton.born_amplitudes(E, theta, 'lab') -> S2_from_M / S2_from_phis
      vs the driver with loops and Delta switched off (runner.run(points, mechs=('Born', 'pi0'))).  Must agree to ~1e-13.
  (2) NNLO prediction: the v2 multipoles (scenarios/bchpt_v2_multipoles.py, J <= 5/2) -> partial waves -> Phi_1..6 ->
      the same helper, vs the driver's S2X/S2Z of the complete NNLO amplitude at Pref (runner.run default mechs):
      (2a) all multipoles of the full amplitude (Born included, l-truncated);
      (2b) non-Born multipoles + the exact Born + pi0 amplitude (as the PWA model, pwa_model.Model.curve, does).
      Agreement at the 1e-3 level (l-truncation) is expected for (2a); (2b) should be better.
Points: E_lab = 80, 100, 140 MeV, theta_lab = 30, 60, 90, 120, 150 deg.
Run from anywhere:  .venv/bin/python work/pwa/validation/validate_S2.py   (writes validate_S2.log next to this file;
driver runs are cached in work/delta_plots/cache and work/pwa/fortran/runs/cache_v2)."""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PWA = os.path.dirname(HERE)
WORK = os.path.dirname(PWA)
sys.path.insert(0, PWA); sys.path.insert(0, os.path.join(PWA, 'scenarios')); sys.path.insert(0, os.path.join(WORK, 'delta_plots'))
os.chdir(os.path.join(WORK, 'delta_plots'))           # runner/config resolve their paths relative to the module, cwd only for logs
import compton as C                                   # noqa: E402
import runner                                         # noqa: E402
import bchpt_v2_multipoles as V2                      # noqa: E402

ES = (80.0, 100.0, 140.0)
THS = (30.0, 60.0, 90.0, 120.0, 150.0)
OBS = ('S2X', 'S2Z')
LOG = []


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s); LOG.append(s)


def driver(mechs, name):
    pts = [('%s_%g_%g' % (o, E, th), E, th, o, 'LAB') for E in ES for th in THS for o in OBS]
    out = runner.run(pts, mechs=mechs, name=name, verbose=False)
    val = dict(zip(out.ids, out.predict_array(out.Pref)))
    assert np.allclose(out.value_file, out.predict_array(out.Pref), atol=1e-12)
    return out.Pref, val


def table(title, drv, cols):
    """cols: dict label -> dict (o, E, th) -> value"""
    say('\n' + title)
    say('%6s %6s | ' % ('E', 'theta') + ' | '.join('%s: %10s %s' % (o, 'driver', ' '.join('%10s' % c for c in cols)) for o in OBS))
    dmax = {(o, c): 0.0 for o in OBS for c in cols}
    for E in ES:
        for th in THS:
            row = '%6.0f %6.0f | ' % (E, th)
            parts = []
            for o in OBS:
                d = drv['%s_%g_%g' % (o, E, th)]
                parts.append('%s: %10.6f %s' % (o, d, ' '.join('%10.6f' % cols[c][(o, E, th)] for c in cols)))
                for c in cols:
                    dmax[(o, c)] = max(dmax[(o, c)], abs(cols[c][(o, E, th)] - d))
            say(row + ' | '.join(parts))
    for o in OBS:
        say('  max |%s(helper) - %s(driver)|: ' % (o, o) + ', '.join('%s %.2e' % (c, dmax[(o, c)]) for c in cols))
    return dmax


# ---------------------------------------------------------------- (1) Born + pi0 pole
pref_b, drv_b = driver(('Born', 'pi0'), 's2check_born')
say('(1) Born + pi0 pole: driver Pref =', pref_b, '(no loops, no Delta: no contact terms)')
assert np.allclose(pref_b, 0.0)
colM, colP = {}, {}
for E in ES:
    for th in THS:
        M = C.born_amplitudes(E, th, 'lab')
        sM = C.S2_from_M(M); sP = C.S2_from_phis(C.phis_from_M(M))
        for o in OBS:
            colM[(o, E, th)] = float(sM[o]); colP[(o, E, th)] = float(sP[o])
d1 = table('(1) Born + pi0: driver vs compton.py (from M: bilinear forms of ampl.f90; from Phi: HMP16 Eq. 3.50)',
           drv_b, {'from_M': colM, 'from_Phi': colP})
assert max(d1.values()) < 1e-10, d1

# ---------------------------------------------------------------- (2) NNLO prediction
pref_n, drv_n = driver(runner.MECHS, 's2check_nnlo')
say('\n(2) NNLO: driver Pref =', pref_n, ' (v2 multipoles REF =', [V2.REF[k] for k in ('alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2')], ')')
assert np.allclose(pref_n, [V2.REF[k] for k in ('alphaE1', 'betaM1', 'gE1E1', 'gM1M1', 'gE1M2', 'gM1E2')], atol=2e-5)
colA, colB, colA_M = {}, {}, {}
for E in ES:
    full = V2.multipoles([E])                                  # full amplitude, J <= 5/2 (jmax=5)
    born = V2.multipoles([E], born_only=True)
    fd = {k: complex(np.asarray(v).ravel()[0]) for k, v in full.items() if not str(k).startswith('_')}
    nb = {k: fd[k] - complex(np.asarray(born[k]).ravel()[0]) for k in fd}
    th = np.array(THS)
    k = C.kin(E, th, 'lab')
    Pfull = C.phis_from_pw(C.pw_from_multipoles(fd, 5), k['thcm'])                      # (2a) l-truncated full amplitude
    Pnb = C.phis_from_pw(C.pw_from_multipoles(nb, 5), k['thcm'])
    PB = C.phis_from_M(C.born_amplitudes(E, np.degrees(k['thcm']), 'cm'))               # exact Born + pi0
    sA = C.S2_from_phis(Pfull); sB = C.S2_from_phis(PB + Pnb); sAM = C.S2_from_M(C.M_from_phis(Pfull))
    for i, t in enumerate(THS):
        for o in OBS:
            colA[(o, E, t)] = float(sA[o][i]); colB[(o, E, t)] = float(sB[o][i]); colA_M[(o, E, t)] = float(sAM[o][i])
d2 = table('(2) NNLO prediction: driver vs multipoles (2a: full amplitude l<=3 truncated, via Phi and via M; 2b: non-Born multipoles + exact Born)',
           drv_n, {'2a_full': colA, '2a_full_M': colA_M, '2b_nB+Born': colB})
assert max(d2.values()) < 5e-3, d2
say('\nPASSED: Born check < 1e-10; NNLO check < 5e-3 (l-truncation)')
open(os.path.join(HERE, 'validate_S2.log'), 'w').write('\n'.join(LOG) + '\n')
print('written', os.path.join(HERE, 'validate_S2.log'))
