"""
theory_dp.py -- BChPT curves and uncertainty bands for the Delta-region figures.

Gram-matrix decomposition
-------------------------
Every observable is a real quadratic form Q in the helicity amplitudes (cross sections) or a ratio N/D of
two such forms (asymmetries), and the amplitude is linear in the mechanisms and contact terms,
    A = sum_i c_i M_i + sum_k d_k C_k ,
    M = (Born, pi0, piN, D1, piD, D0)
        D1 = Delta pole with the gamma-N-Delta vertex corrections (code_v2 dvertex=t; the main amplitude)
        D0 = Delta pole without them (dvertex=f, = v1),
    C_k = O(p^4) contact term per unit shift of the k-th dipole polarisability (d = P - Pref).
Driver runs of every single mechanism and every pair (dvertex=t: 5 + 10 runs; dvertex=f: the D0 single and
the 4 pairs with D0) give by polarisation the real Gram matrix H (12 x 12) of each form at each point:
    H_ii = Q(M_i),  H_ij = [Q(M_i + M_j) - Q(M_i) - Q(M_j)]/2,  H_{i,C_k} = c_0k(run i)/2,  H_{C_k C_l} = c_kl.
The single entry not obtainable from runs, <D1, D0>, is set to [Q(D1) + Q(D0)]/2, i.e. Q(V) := 0 for the
vertex-correction amplitude V = D1 - D0.  This defines the delta x (Delta-pole corrections) shift of the band:
for every quadratic form  Q(A + l V) -> Q(A) + l [Q(A) - Q(A - V)] (+ the exact cross terms l^2 2<L,V> with the
loops in the 'sum' shift), i.e. the shift is delta times the actual change of the observable's numerator and
denominator caused by the vertex corrections (linear in the vertex-correction effect).  The amplitude-level
shift of LMP2015 would add +-l(1+l) Q(V); V is mostly 90 deg out of phase with the Delta pole (Watson), so it
changes observables much less than |V| suggests -- the linear prescription measures that change directly.
Every other curve (main, without VC, tree, spin fit, loop shifts) is exact.
Then Q(c, d) = v^T H v with v = (c, d) for ANY linear combination, e.g. the LMP2015 Eq. (32) shifted
amplitudes A +- delta dA_X.  validate() compares with direct runs (complete amplitude with and without the
vertex corrections, tree level, a 3-mechanism configuration with contact terms).

Curves (config.py):
  main   : A = Born + pi0 + piN + D1 + piD, alphaE1, betaM1 = main fit, spin polarisabilities = NNLO (Pref)
  noVC   : same with D0 instead of D1 ("without vertex corrections", = v1 amplitude with the v2 masses)
  trunc  : envelope over X in {piN, piD, VC = D1 - D0, sum} and +- of O(A +- delta(omega_cm) dA_X),
           static alphaE1, betaM1 of the loop shifts removed (they are fitted); delta = omega_cm / Lambda_chi
  ab     : (alpha, beta) fit error, linear propagation with the fit covariance
  band   : main -+ sqrt(trunc_side^2 + ab^2)
  gM     : envelope of the complete amplitude with gM = 2.87, 3.07 (alpha, beta fixed; spin pols = NNLO of gM)
  tree   : Born + pi0 + D1 (no loops, no contact terms)
  spin   : complete amplitude at the parameter set of the Delta-region spin-polarisability fit
"""
import itertools
import numpy as np
import config as C
import runner
from runner import MECHS

MB = ('Born', 'pi0', 'piN', 'D1', 'piD', 'D0')     # Gram basis (mechanisms)
NB = len(MB)
NT = NB + 6
IDX = {m: i for i, m in enumerate(MB)}
RUNMECH = {'Born': 'Born', 'pi0': 'pi0', 'piN': 'piN', 'D1': 'Dpole', 'piD': 'piD', 'D0': 'Dpole'}
ASYM = ('S3', 'S2X', 'S2Z')
COMPKEY = {'piN': 'piN loops', 'D1': 'Delta pole', 'D0': 'Delta pole', 'piD': 'piDelta loops'}
PIDX = runner.PIDX


def omega_cm(E):
    """c.m. photon energy [MeV] for lab photon energy E [MeV]"""
    E = np.asarray(E, float)
    return E * C.MP / np.sqrt(C.MP ** 2 + 2 * C.MP * E)


def delta_of(E):
    return omega_cm(E) / C.LAMBDA_CHI


def cvec(**kw):
    c = np.zeros(NB)
    for k, v in kw.items():
        c[IDX[k]] = v
    return c


C_MAIN = cvec(Born=1, pi0=1, piN=1, D1=1, piD=1)
C_NOVC = cvec(Born=1, pi0=1, piN=1, D0=1, piD=1)
C_TREE = cvec(**{({'Dpole': 'D1'}.get(m, m)): 1 for m in C.TREE_MECHS})
SHIFT_VEC = {'piN': cvec(piN=1), 'piD': cvec(piD=1), 'VC': cvec(D1=1, D0=-1)}


def predict_d(T, d):
    """driver observable at contact shifts d (relative to the run's own reference = no contact terms at d=0)"""
    return T.predict_array(T.Pref + np.asarray(d, float))


class Theory:
    def __init__(self, points, tag='all', verbose=True):
        self.points = list(points)
        self.pid = np.array([p[0] for p in self.points])
        self.E = np.array([float(p[1]) for p in self.points])
        self.obs = np.array([p[3] for p in self.points])
        self.isasym = np.isin(self.obs, ASYM)
        self.tag = tag
        n = len(self.points)
        run = lambda mechs, nm, dv=True: runner.run(self.points, mechs=mechs, name='%s_%s' % (tag, nm), verbose=verbose, dvertex=dv)
        R = {}                                   # Gram-basis runs: key = frozenset of basis labels
        for b in MB:
            R[frozenset([b])] = run((RUNMECH[b],), b, b != 'D0')
        for a, b in itertools.combinations(MB, 2):
            if {a, b} == {'D1', 'D0'}:
                continue
            R[frozenset([a, b])] = run((RUNMECH[a], RUNMECH[b]), a + '+' + b, not ('D0' in (a, b)))
        self.full = run(MECHS, 'NNLO_VC', True)
        self.full_noVC = run(MECHS, 'NNLO_noVC', False)
        self.Pref = self.full.Pref.copy()
        for T in R.values():
            assert np.array_equal(T.ids, self.full.ids)
        # static polarisabilities of the mechanisms: v2 sets Pref = static pols of the mechanisms switched on
        # (incl. physical-mass shifts); v1 writes the complete NNLO values for every configuration -> header.
        self.static = np.zeros((NB, 6))
        v2pref = not np.allclose(R[frozenset(['piN'])].Pref, self.Pref)
        for b in ('piN', 'D1', 'piD', 'D0'):
            self.static[IDX[b]] = R[frozenset([b])].Pref if v2pref else self.full.comp[COMPKEY[b]]
        self.v2pref = v2pref
        self.H = {}
        for form in ('N', 'D'):
            get = (lambda T: T.CN) if form == 'N' else (lambda T: T.CD)
            H = np.zeros((n, NT, NT))
            Q1 = {b: get(R[frozenset([b])])[:, PIDX[0, 0]] for b in MB}
            for b in MB:
                i = IDX[b]
                H[:, i, i] = Q1[b]
                for k in range(1, 7):
                    H[:, i, NB + k - 1] = H[:, NB + k - 1, i] = 0.5 * get(R[frozenset([b])])[:, PIDX[0, k]]
            for key, T in R.items():
                if len(key) == 2:
                    a, b = sorted(key, key=MB.index)
                    H[:, IDX[a], IDX[b]] = H[:, IDX[b], IDX[a]] = 0.5 * (get(T)[:, PIDX[0, 0]] - Q1[a] - Q1[b])
            # <D1, D0>: Q(D1 - D0) neglected (see module doc)
            H[:, IDX['D1'], IDX['D0']] = H[:, IDX['D0'], IDX['D1']] = 0.5 * (Q1['D1'] + Q1['D0'])
            cf = get(self.full)
            for k in range(1, 7):
                for l in range(k, 7):
                    H[:, NB + k - 1, NB + l - 1] = H[:, NB + l - 1, NB + k - 1] = cf[:, PIDX[k, l]] * (1.0 if k == l else 0.5)
            self.H[form] = H
        self.H['D'][~self.isasym] = 0.0           # cross sections: denominator == 1
        self.hasVC = bool(np.max(np.abs(self.H['N'][:, IDX['D1'], IDX['D1']] - self.H['N'][:, IDX['D0'], IDX['D0']])
                                 / np.maximum(np.abs(self.H['N'][:, IDX['D1'], IDX['D1']]), 1e-30)) > 1e-9)
        self._gm = None

    # ------------------------------------------------------------------ evaluation
    def nd(self, c, d):
        n = len(self.points)
        c = np.broadcast_to(np.asarray(c, float), (n, NB))
        d = np.broadcast_to(np.asarray(d, float), (n, 6))
        v = np.concatenate([c, d], axis=1)
        num = np.einsum('pi,pij,pj->p', v, self.H['N'], v)
        den = np.where(self.isasym, np.einsum('pi,pij,pj->p', v, self.H['D'], v), 1.0)
        return num, den

    def evaluate(self, c, d):
        num, den = self.nd(c, d)
        return num / den

    def d_of(self, P):
        return np.asarray(P, float) - self.Pref

    def P_main(self):
        return np.r_[C.ALPHA_MAIN, C.BETA_MAIN, self.Pref[2:]]

    def main(self):
        return self.evaluate(C_MAIN, self.d_of(self.P_main()))

    def noVC(self):
        return self.evaluate(C_NOVC, self.d_of(self.P_main()))

    def tree(self):
        return self.evaluate(C_TREE, np.zeros(6))

    def spinfit(self):
        return self.evaluate(C_MAIN, self.d_of(C.P_SPINFIT))

    def ab_error(self, eps=0.01):
        P0 = self.P_main()
        g = []
        for k in range(2):
            dp = np.zeros(6); dp[k] = eps
            g.append((self.evaluate(C_MAIN, self.d_of(P0 + dp)) - self.evaluate(C_MAIN, self.d_of(P0 - dp))) / (2 * eps))
        g = np.array(g)
        return np.sqrt(np.einsum('ip,ij,jp->p', g, C.COV_AB, g))

    def truncation(self):
        """LMP2015 Eq. (32)-type envelope; returns (lo, hi, dict of the individual shifted curves)"""
        dlt = delta_of(self.E)[:, None]
        d0 = self.d_of(self.P_main())
        O0 = self.evaluate(C_MAIN, d0)
        lo, hi = O0.copy(), O0.copy()
        parts = {}
        for name, comps in C.SHIFT_SETS.items():
            comps = [x for x in comps if x != 'VC' or self.hasVC]
            if not comps:
                continue
            e = np.sum([SHIFT_VEC[x] for x in comps], axis=0)
            s = np.zeros(6)
            for x in comps:
                if x in ('piN', 'piD'):
                    s += self.static[IDX[x]]
            sab = np.r_[s[:2], np.zeros(4)] if C.REMOVE_STATIC_AB else np.zeros(6)
            for sg in (+1, -1):
                O = self.evaluate(C_MAIN[None, :] + sg * dlt * e[None, :], d0[None, :] - sg * dlt * sab[None, :])
                parts['%s%+d' % (name, sg)] = O
                lo = np.minimum(lo, O); hi = np.maximum(hi, O)
        return lo, hi, parts

    def gM_runs(self):
        if self._gm is None:
            out = {}
            for g in (C.GM_LO, C.GM_HI):
                T = runner.run(self.points, gM=g, name='%s_NNLO_VC_gM%.2f' % (self.tag, g))
                out[g] = T.predict_array(np.r_[C.ALPHA_MAIN, C.BETA_MAIN, T.Pref[2:]])
                out['Pref%.2f' % g] = T.Pref
            self._gm = out
        return self._gm

    def curves(self):
        main = self.main()
        tlo, thi, parts = self.truncation()
        ab = self.ab_error()
        lo = main - np.sqrt((main - tlo) ** 2 + ab ** 2)
        hi = main + np.sqrt((thi - main) ** 2 + ab ** 2)
        gm = self.gM_runs()
        glo = np.minimum.reduce([main, gm[C.GM_LO], gm[C.GM_HI]])
        ghi = np.maximum.reduce([main, gm[C.GM_LO], gm[C.GM_HI]])
        out = dict(main=main, lo=lo, hi=hi, trunc_lo=tlo, trunc_hi=thi, ab=ab, gm_lo=glo, gm_hi=ghi,
                   gm_287=gm[C.GM_LO], gm_307=gm[C.GM_HI], tree=self.tree(), spin=self.spinfit(),
                   flag=self.full.flag.copy())
        for k, v in parts.items():
            out['shift_' + k] = v
        if self.hasVC:
            out['noVC'] = self.noVC()
        num, den = self.nd(C_MAIN, self.d_of(self.P_main()))
        out['num'], out['den'] = num, den
        return out

    # ------------------------------------------------------------------ validation
    def validate(self):
        rep = {}
        for lab, P in [('complete (VC), P=Pref', self.Pref), ('complete (VC), P=main fit', self.P_main()),
                       ('complete (VC), P=spin fit', C.P_SPINFIT)]:
            rep[lab] = self._dev(self.evaluate(C_MAIN, self.d_of(P)), self.full.predict_array(P))
        dm = self.d_of(self.P_main())
        rep['complete without VC, P=main fit'] = self._dev(self.evaluate(C_NOVC, dm), predict_d(self.full_noVC, dm))
        Tt = runner.run(self.points, mechs=('Born', 'pi0', 'Dpole'), name='%s_check_tree' % self.tag)
        rep['tree Born+pi0+D1'] = self._dev(self.tree(), predict_d(Tt, np.zeros(6)))
        Tn = runner.run(self.points, mechs=('Born', 'pi0', 'piN'), name='%s_check_NLO' % self.tag)
        rep['Born+pi0+piN, d=main fit'] = self._dev(self.evaluate(cvec(Born=1, pi0=1, piN=1), dm), predict_d(Tn, dm))
        Tq = runner.run(self.points, mechs=('Born', 'piN', 'Dpole', 'piD'), name='%s_check_noPi0_noVC' % self.tag, dvertex=False)
        rep['Born+piN+D0+piD, d=main fit'] = self._dev(self.evaluate(cvec(Born=1, piN=1, D0=1, piD=1), dm), predict_d(Tq, dm))
        rep['static pols per mechanism (piN, D1, piD)'] = {b: [round(float(x), 4) for x in self.static[IDX[b]]] for b in ('piN', 'D1', 'piD')}
        rep['Pref (complete)'] = [round(float(x), 6) for x in self.Pref]
        rep['vertex corrections present'] = self.hasVC
        return rep

    def _dev(self, g, r):
        cs = ~self.isasym
        out = {}
        if cs.any():
            out['max rel. dev. (cross sections)'] = float(np.max(np.abs(g[cs] - r[cs]) / np.maximum(np.abs(r[cs]), 1e-12)))
        if self.isasym.any():
            out['max abs. dev. (asymmetries)'] = float(np.max(np.abs(g[self.isasym] - r[self.isasym])))
        return out
