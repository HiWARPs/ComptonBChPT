"""pi0ff_multipoles.py -- shift of the dynamical dipole polarisabilities (l = 1 non-Born multipoles of the paper's
Fig. 'multipoles') caused by the pi0 gamma gamma transition form factor of the proton triangle loop on the pi0-pole
vertex (work/code_v2 &v2 pi0ff=t, driver/README.md; work/notes/pi0_formfactor.md).

Projection: bv2mult.py (copy of work/pwa/scenarios/bchpt_v2_multipoles.py: main program &multipoles pwa=t, J <= 5/2,
16-point Gauss rule, printed multipoles -> fm).  Conventions of the paper figure (work/pwa/make_figures.DYN):
  alpha(w) = [2 f_EE^{1+} + f_EE^{1-}]/w^2,  beta(w) = [2 f_MM^{1+} + f_MM^{1-}]/w^2,
  gE1E1(w) = [f_EE^{1+} - f_EE^{1-}]/w^3, gM1M1(w) = [f_MM^{1+} - f_MM^{1-}]/w^3, gE1M2(w) = 6 f_EM^{1+}/w^3,
  gM1E2(w) = 6 f_ME^{1+}/w^3,  w = c.m. photon energy = E_lab/sqrt(1+2E_lab/M) (fm^-1), f in fm, result in 1e-4 fm^n.
Non-Born = full - (Born + POINTLIKE pi0 pole), the convention of the figure, for both pi0ff=f and pi0ff=t.  The shift
is then the projection of [F(t) - 1] x (pi0 pole) onto the l = 1 multipoles; it does not depend on the LECs (the pole
is LEC independent) and it is checked here against the difference of two full-amplitude runs.  If the form factor is
instead regarded as part of the Born + pi0 reference, the non-Born multipoles are unchanged identically.
Energies: 59-150 MeV (figure range above the lowest PWA bin) and 200-300 MeV (real parts; above threshold the
multipoles are complex, the pole is real).
Run:  python3 pi0ff_multipoles.py
Outputs (this directory): pi0ff_multipoles.json, pi0ff_multipoles.md; cache_mult/ (run cache)."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import bv2mult as BV

HB = 197.327                         # make_figures.HB (MeV fm)
MN_MEV = 0.93827 * 1e3               # compton.MN_GEV
DYN = [('alphaE1', [('EE', 1, '+', 2), ('EE', 1, '-', 1)], 2),
       ('betaM1', [('MM', 1, '+', 2), ('MM', 1, '-', 1)], 2),
       ('gE1E1', [('EE', 1, '+', 1), ('EE', 1, '-', -1)], 3),
       ('gM1M1', [('MM', 1, '+', 1), ('MM', 1, '-', -1)], 3),
       ('gE1M2', [('EM', 1, '+', 6)], 3),
       ('gM1E2', [('ME', 1, '+', 6)], 3)]
V2F = BV.V2OPTS
V2T = BV.V2OPTS + ', pi0ff=t'
FITS_JSON = os.path.join(HERE, 'pi0ff_fits.json')


def dyn(f, E, comb, n):
    w = (E / HB) / np.sqrt(1 + 2 * E / MN_MEV)
    return sum(c * f[(a, b, s)] for a, b, s, c in comb) / w ** n * 1e4


def dyn_all(f, E):
    return {k: dyn(f, E, comb, n) for k, comb, n in DYN}


def run(E, label):
    full_f = BV.multipoles(E, None, v2=V2F); full_t = BV.multipoles(E, None, v2=V2T)
    born_f = BV.multipoles(E, None, born_only=True, v2=V2F); born_t = BV.multipoles(E, None, born_only=True, v2=V2T)
    keys = [k for k in full_f if not str(k).startswith('_')]
    nb_f = {k: full_f[k] - born_f[k] for k in keys}                   # non-Born, pi0ff=f  (= paper figure, NNLO pols)
    nb_t = {k: full_t[k] - born_f[k] for k in keys}                   # non-Born with pointlike reference, pi0ff=t
    d_pole = {k: born_t[k] - born_f[k] for k in keys}                 # [F-1] x pole, projected (exact, no loops)
    d_full = {k: full_t[k] - full_f[k] for k in keys}                 # same from the full runs (LoopTools jitter)
    pole = None
    D0 = dyn_all({k: nb_f[k].real for k in keys}, E)
    Dt = dyn_all({k: nb_t[k].real for k in keys}, E)
    Dp = dyn_all({k: d_pole[k].real for k in keys}, E)
    Df = dyn_all({k: d_full[k].real for k in keys}, E)
    jitter = max(float(np.max(np.abs(Dp[k] - Df[k]))) for k in Dp)
    out = dict(E=E.tolist(), dyn_f={k: D0[k].tolist() for k in D0}, dyn_t={k: Dt[k].tolist() for k in Dt},
               shift={k: Dp[k].tolist() for k in Dp}, shift_from_full_runs={k: Df[k].tolist() for k in Df},
               jitter_check=jitter)
    print('=== %s: E = %g-%g MeV (%d energies); shift from (Born+pi0)(pi0ff) - (Born+pi0) vs full(pi0ff) - full: max |diff| = %.1e' % (
        label, E[0], E[-1], len(E), jitter))
    print('%-8s %12s %12s %12s %12s %12s' % ('pol', 'max|shift|', 'at E', 'shift(Emax)', 'value(Emax)', 'rel(Emax)'))
    summ = {}
    for k in D0:
        i = int(np.argmax(np.abs(Dp[k])))
        summ[k] = dict(max_abs_shift=float(np.abs(Dp[k][i])), at_E=float(E[i]), shift_Emax=float(Dp[k][-1]),
                       value_Emax=float(D0[k][-1]), value_E0=float(D0[k][0]), rel_Emax=float(Dp[k][-1] / D0[k][-1]) if D0[k][-1] != 0 else None)
        print('%-8s %12.5f %12.1f %+12.5f %12.4f %12.2e' % (k, summ[k]['max_abs_shift'], summ[k]['at_E'], summ[k]['shift_Emax'],
                                                          summ[k]['value_Emax'], summ[k]['rel_Emax'] or 0))
    out['summary'] = summ
    return out


if __name__ == '__main__':
    E1 = np.arange(59.0, 150.001, 1.0)
    E2 = np.arange(200.0, 300.001, 5.0)
    R = {'below': run(E1, 'below threshold'), 'delta': run(E2, 'Delta region')}
    # the pi0-pole dynamical polarisabilities themselves (for scale): pole = (Born+pi0) - Born(anomaly off)
    # anomaly off: run the Born-only configuration with anomaly=f through the INDAT template
    tmpl = BV.INDAT
    BV.INDAT = tmpl.replace('&incl anomaly=t,', '&incl anomaly=f,')
    born_noan = BV.multipoles(E1, None, born_only=True, v2=V2F)
    BV.INDAT = tmpl
    born_f = BV.multipoles(E1, None, born_only=True, v2=V2F)
    keys = [k for k in born_f if not str(k).startswith('_')]
    pole = dyn_all({k: (born_f[k] - born_noan[k]).real for k in keys}, E1)
    R['pole_dyn_below'] = {k: v.tolist() for k, v in pole.items()}
    print('=== pi0-pole contribution to the dynamical polarisabilities (pointlike), 59 and 150 MeV:')
    for k in pole:
        print('   %-8s %9.4f %9.4f' % (k, pole[k][0], pole[k][-1]))
    json.dump(R, open(os.path.join(HERE, 'pi0ff_multipoles.json'), 'w'), indent=1)
    L = ['# pi0ff: shift of the dynamical dipole polarisabilities (l = 1 non-Born multipoles)', '',
         'Produced by `pi0ff_multipoles.py` (data in `pi0ff_multipoles.json`).  Non-Born = full - (Born + pointlike pi0 pole)',
         'in both cases (paper convention); the shift is the l = 1 projection of [F(t) - 1] x (pi0 pole), independent of the LECs.',
         'Units: 1e-4 fm^3 (alpha, beta), 1e-4 fm^4 (gamma).  Values at the NNLO polarisabilities.', '']
    for lab, key in [('59-150 MeV (below threshold)', 'below'), ('200-300 MeV (Delta region, real parts)', 'delta')]:
        s = R[key]['summary']
        L += ['## %s' % lab, '', '| polarisability | max abs shift | at E [MeV] | shift at E_max | value at E_max (pi0ff=f) | relative |', '|---|---|---|---|---|---|']
        for k in s:
            L.append('| %s | %.5f | %.0f | %+.5f | %.4f | %.1e |' % (k, s[k]['max_abs_shift'], s[k]['at_E'], s[k]['shift_Emax'], s[k]['value_Emax'], s[k]['rel_Emax'] or 0))
        L += ['', 'Check: shift from the two Born+pi0 runs vs from the two full runs agree to %.1e (LoopTools history jitter).' % R[key]['jitter_check'], '']
    L += ['## pi0-pole (pointlike) dynamical polarisabilities, for scale', '', '| polarisability | 59 MeV | 150 MeV |', '|---|---|---|']
    for k in pole:
        L.append('| %s | %.4f | %.4f |' % (k, pole[k][0], pole[k][-1]))
    L += ['', 'Shift at selected energies (below threshold):', '', '| E [MeV] | ' + ' | '.join(k for k, _, _ in DYN) + ' |', '|---|' + '---|' * 6]
    for Ev in [59.0, 80.0, 100.0, 120.0, 140.0, 150.0]:
        i = int(np.argmin(np.abs(E1 - Ev)))
        L.append('| %g | ' % Ev + ' | '.join('%+.5f' % R['below']['shift'][k][i] for k, _, _ in DYN) + ' |')
    L += ['', '| E [MeV] | ' + ' | '.join(k for k, _, _ in DYN) + ' |', '|---|' + '---|' * 6]
    for Ev in [200.0, 250.0, 300.0]:
        i = int(np.argmin(np.abs(E2 - Ev)))
        L.append('| %g | ' % Ev + ' | '.join('%+.5f' % R['delta']['shift'][k][i] for k, _, _ in DYN) + ' |')
    open(os.path.join(HERE, 'pi0ff_multipoles.md'), 'w').write('\n'.join(L) + '\n')
    print('written pi0ff_multipoles.json / .md')
