"""Empirical non-Born forward amplitudes fbar(nu), gbar(nu) from the table of Re f, Re g at 12 energies 50-150 MeV.
Source switch FWD_SOURCE (environment): '2026' (default since 2026-10-18) = the 2026 update of the ForwardRCS repository
(github.com/HiWARPs/ForwardRCS; work/data/db_forward_amp.csv built from work/data/raw/ForwardRCS2026/forward_points.json,
sum-rule coefficients from unpol_sumrules.json / pol_sumrules.json of the same directory); 'ghp2016' = Gryniuk-
Hagelstein-Pascalutsa (PRD 92 (2015) 074031, PRD 94 (2016) 034043), table work/data/db_forward_amp_ghp2016.csv.
Born parts subtracted with GHP constants: f_B = -alpha/M, g_B = -alpha kappa^2 nu/(2M^2) (exact for the
Born+pi0 amplitude used here, verified numerically in validate_machinery.py).
Interpolation (ghp2016; 2026: see ForwardAmps.__init__, the 5-MeV grid of the repository): for nu <= 135 MeV a smooth fit with the low-energy (sum-rule) coefficients
  fbar = (alpha+beta) nu^2 + c4 nu^4 + c6 nu^6 + nu^8 (a8 + a10 nu^2)
  gbar = gamma0 nu^3 + gamma0bar nu^5 + nu^7 (b7 + b9 nu^2)
with alpha+beta, c4, c6, gamma0, gamma0bar of the source: 2026 = 14.009, 5.950, 4.293, -1.0248, 0.4484 (ForwardRCS 2026,
gamma's 'unconstrained', i.e. without the GDH sum rule imposed); ghp2016 = 14.0, 6.04, 4.39 (GHP 2015, fits I/II averaged),
-0.929, 0.484 (GHP 2016),
the remaining coefficients fitted to the table points; above 135 MeV (pion-threshold cusp region)
a monotone cubic (PCHIP) through the table points 130-150 MeV joined continuously.
Errors: err(nu) = total error of the table (PCHIP); split(nu) = (coherent f, coherent g, uncorrelated f, uncorrelated g)
for the nuisance model of pwa_model (one nuisance per amplitude): 2026: uncorrelated = err_unc column (remainder after the
leading exp+fit eigenmode (+) systematic mode), coherent = sqrt(total^2 - unc^2); ghp2016: coherent = total, uncorrelated =
0.0029 mub GeV (the 2-decimal rounding; the model used until 2026-10-18).
Constraint values of pwa_model (sr_mode 'chi2' centres and widths of alpha+beta and gamma0, units 1e-4 fm^3, fm^4):
AB_SR, SIG_AB, G0_SR, SIG_G0 (2026: 14.009(208), -1.025(102), exp+fit (+) data selection resp. PWA model in quadrature).
Units: fm (amplitudes), MeV (nu)."""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import os
import json
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator, CubicSpline

HBARC_MEV = 197.3269804
UBGEV_PER_FM = 1973.27          # 1 fm = 1973.27 mub GeV
MP_GHP = 938.272
ALPHA = 1 / 137.035999
KAPPA = 1.792847
DATA = (_WORK + '/data')
SOURCE = os.environ.get('FWD_SOURCE', '2026')
assert SOURCE in ('2026', 'ghp2016'), SOURCE
R_ROUND_GHP = 0.0029 / UBGEV_PER_FM          # fm: uncorrelated (rounding) part of the GHP table
if SOURCE == 'ghp2016':
    AB, C4, C6 = 14.0e-4, 6.04e-4, 4.39e-4       # fm^3, fm^5, fm^7
    G0, G0B = -0.929e-4, 0.484e-4                 # fm^4, fm^6
    AB_SR, SIG_AB, G0_SR, SIG_G0 = 14.0, 0.2, -0.929, 0.105
    TABLE = DATA + '/db_forward_amp_ghp2016.csv'
else:
    _u = json.load(open(DATA + '/raw/ForwardRCS2026/unpol_sumrules.json'))
    _p = json.load(open(DATA + '/raw/ForwardRCS2026/pol_sumrules.json'))
    AB, C4, C6 = [_u['sum_rules'][k]['regions']['total'][0] * 1e-4 for k in ('2', '4', '6')]
    G0, G0B = [_p['unconstrained'][k]['value'] * 1e-6 for k in ('gamma0', 'gamma0bar')]
    AB_SR, G0_SR = AB * 1e4, G0 * 1e4
    SIG_AB = float(np.hypot(_u['sum_rules']['2']['regions']['total'][1], _u['selection_systematic']['2']['half_difference']))
    SIG_G0 = float(np.hypot(_p['unconstrained']['gamma0']['err'], _p['pwa_model_systematic']['unconstrained']['gamma0'])) / 100
    TABLE = DATA + '/db_forward_amp.csv'


def load_table(path=TABLE, unc=False):
    d = pd.read_csv(path)
    f = d[d.obs == 'FWD_REF'].sort_values('E_lab_MeV')
    g = d[d.obs == 'FWD_REG'].sort_values('E_lab_MeV')
    nu = f.E_lab_MeV.values.astype(float)
    Mfm = MP_GHP / HBARC_MEV
    nufm = nu / HBARC_MEV
    fB = -ALPHA / Mfm
    gB = -ALPHA * KAPPA**2 * nufm / (2 * Mfm**2)
    fbar = f.value.values / UBGEV_PER_FM - fB
    gbar = g.value.values / UBGEV_PER_FM - gB
    out = (nu, fbar, f.stat_err.values / UBGEV_PER_FM, gbar, g.stat_err.values / UBGEV_PER_FM)
    if unc:
        if 'err_unc' in d:
            out += (f.err_unc.values / UBGEV_PER_FM, g.err_unc.values / UBGEV_PER_FM)
        else:
            out += (np.full(len(nu), np.nan), np.full(len(nu), np.nan))
    return out


class ForwardAmps:
    def __init__(self, path=None):
        nu, fb, fe, gb, ge, fu, gu = load_table(path, unc=True) if path else load_table(unc=True)
        self.nu, self.fb, self.fe, self.gb, self.ge = nu, fb, fe, gb, ge
        self.has_unc = bool(np.all(np.isfinite(fu)))
        x = nu / HBARC_MEV
        lo = nu <= 135.0
        # least-squares for the remainder coefficients (weights 1/err^2)
        Af = np.stack([x**8, x**10], -1)[lo]
        rf = (fb - AB * x**2 - C4 * x**4 - C6 * x**6)[lo]
        self.af = np.linalg.lstsq(Af / fe[lo, None], rf / fe[lo], rcond=None)[0]
        Ag = np.stack([x**7, x**9], -1)[lo]
        rg = (gb - G0 * x**3 - G0B * x**5)[lo]
        self.ag = np.linalg.lstsq(Ag / ge[lo, None], rg / ge[lo], rcond=None)[0]
        hi = nu >= 130.0
        # cusp region: PCHIP through the (smooth value at 130,135) + table at 145,150
        nn = nu[hi]
        fv = np.where(nn <= 135, self._fs(nn / HBARC_MEV), fb[hi])
        gv = np.where(nn <= 135, self._gs(nn / HBARC_MEV), gb[hi])
        self.fh = PchipInterpolator(nn, fv)
        self.gh = PchipInterpolator(nn, gv)
        self.ferr = PchipInterpolator(nu, fe)
        self.gerr = PchipInterpolator(nu, ge)
        if self.has_unc:
            self.func = PchipInterpolator(nu, fu)
            self.gunc = PchipInterpolator(nu, gu)
        # 2026: the repository tabulates Re f, Re g on a 5-MeV grid (unpol_f.txt, pol_g.txt; identical to the 12 table
        # values at the table energies).  The LEX representation above misses the approach to the pi0 cusp by up to
        # 0.0016 (f) / 0.0029 (g) mub GeV at 117-135 MeV, i.e. by more than the 2026 uncorrelated errors; hence fbar, gbar
        # interpolate the grid instead: cubic spline for nu <= 140 MeV, PCHIP through 130...150 MeV across the cusp.
        self.grid = None
        if SOURCE == '2026' and path is None:
            def rd(fn):
                rows = [l.split()[:2] for l in open(DATA + '/raw/ForwardRCS2026/' + fn) if l.strip() and l.split()[0][0].isdigit()]
                a = np.array(rows, float)
                return 1000 * a[:, 0], a[:, 1]
            nf_, f_ = rd('unpol_f.txt'); ng_, g_ = rd('pol_g.txt')
            m = (nf_ <= 150.01); assert np.allclose(ng_[m], nf_[m])
            ng = nf_[m]; Mfm = MP_GHP / HBARC_MEV
            fbg = f_[m] / UBGEV_PER_FM + ALPHA / Mfm
            gbg = g_[:m.sum()] / UBGEV_PER_FM + ALPHA * KAPPA**2 * (ng / HBARC_MEV) / (2 * Mfm**2)
            lo, hi = ng <= 140.01, ng >= 129.99
            self.grid = (CubicSpline(ng[lo], fbg[lo]), CubicSpline(ng[lo], gbg[lo]), PchipInterpolator(ng[hi], fbg[hi]),
                         PchipInterpolator(ng[hi], gbg[hi]))
            assert np.allclose(self.fbar(nu), fb, atol=1e-9 / UBGEV_PER_FM) and np.allclose(self.gbar(nu), gb, atol=1e-9 / UBGEV_PER_FM)

    def _fs(self, x):
        return AB * x**2 + C4 * x**4 + C6 * x**6 + self.af[0] * x**8 + self.af[1] * x**10

    def _gs(self, x):
        return G0 * x**3 + G0B * x**5 + self.ag[0] * x**7 + self.ag[1] * x**9

    def fbar(self, nu):
        nu = np.asarray(nu, float)
        if self.grid is not None:
            return np.where(nu <= 140.0, self.grid[0](np.clip(nu, 0, 140)), self.grid[2](np.clip(nu, 140, 150)))
        return np.where(nu <= 135.0, self._fs(nu / HBARC_MEV), self.fh(np.clip(nu, 130, 150)))

    def gbar(self, nu):
        nu = np.asarray(nu, float)
        if self.grid is not None:
            return np.where(nu <= 140.0, self.grid[1](np.clip(nu, 0, 140)), self.grid[3](np.clip(nu, 140, 150)))
        return np.where(nu <= 135.0, self._gs(nu / HBARC_MEV), self.gh(np.clip(nu, 130, 150)))

    def err(self, nu):
        nu = np.asarray(nu, float)
        return self.ferr(np.clip(nu, 50, 150)), self.gerr(np.clip(nu, 50, 150))

    def split(self, nu):
        """(coherent f, coherent g, uncorrelated f, uncorrelated g) errors at nu [fm] for the nuisance model; without
        an err_unc column (GHP table): coherent = total, uncorrelated = the rounding R_ROUND_GHP (model until 2026-10-18)"""
        ef, eg = self.err(nu)
        if not self.has_unc:
            r = np.full(np.shape(ef), R_ROUND_GHP)
            return ef, eg, r, r
        x = np.clip(np.asarray(nu, float), 50, 150)
        uf, ug = self.func(x), self.gunc(x)
        return np.sqrt(np.maximum(ef**2 - uf**2, 0)), np.sqrt(np.maximum(eg**2 - ug**2, 0)), uf, ug


if __name__ == '__main__':
    F = ForwardAmps()
    print('nu  fbar_tab  fbar_smooth  err  |  gbar_tab gbar_smooth err   [mub GeV]')
    for i, n in enumerate(F.nu):
        print(f'{n:5.0f} {F.fb[i]*UBGEV_PER_FM:8.4f} {F.fbar(n)*UBGEV_PER_FM:8.4f} {F.fe[i]*UBGEV_PER_FM:6.4f} | '
              f'{F.gb[i]*UBGEV_PER_FM:8.4f} {F.gbar(n)*UBGEV_PER_FM:8.4f} {F.ge[i]*UBGEV_PER_FM:6.4f}')
    print('remainder coefs f:', F.af, ' g:', F.ag)
    for n in (59, 69, 79, 89, 99, 109, 117, 127, 135, 143, 149, 150):
        print(n, f'{F.fbar(n)*UBGEV_PER_FM:.4f}', f'{F.gbar(n)*UBGEV_PER_FM:.4f}')
