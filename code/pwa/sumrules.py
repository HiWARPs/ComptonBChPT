"""Empirical non-Born forward amplitudes fbar(nu), gbar(nu) of Gryniuk-Hagelstein-Pascalutsa
(PRD 92 (2015) 074031, PRD 94 (2016) 034043), from their table of Re f, Re g (work/data/db_forward_amp.csv).
Born parts subtracted with GHP constants: f_B = -alpha/M, g_B = -alpha kappa^2 nu/(2M^2) (exact for the
Born+pi0 amplitude used here, verified numerically in validate_machinery.py).
Interpolation: for nu <= 135 MeV a smooth fit with the low-energy (sum-rule) coefficients
  fbar = (alpha+beta) nu^2 + c4 nu^4 + c6 nu^6 + nu^8 (a8 + a10 nu^2)
  gbar = gamma0 nu^3 + gamma0bar nu^5 + nu^7 (b7 + b9 nu^2)
with alpha+beta=14.0, c4=6.04, c6=4.39 (GHP 2015, fits I/II averaged) and gamma0=-0.929, gamma0bar=0.484 (GHP 2016),
the remaining coefficients fitted to the table points; above 135 MeV (pion-threshold cusp region)
a monotone cubic (PCHIP) through the table points 130-150 MeV joined continuously.
Units: fm (amplitudes), MeV (nu)."""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

HBARC_MEV = 197.3269804
UBGEV_PER_FM = 1973.27          # 1 fm = 1973.27 mub GeV
MP_GHP = 938.272
ALPHA = 1 / 137.035999
KAPPA = 1.792847
AB, C4, C6 = 14.0e-4, 6.04e-4, 4.39e-4       # fm^3, fm^5, fm^7
G0, G0B = -0.929e-4, 0.484e-4                 # fm^4, fm^6


def load_table(path=(_WORK + '/data/db_forward_amp.csv')):
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
    return nu, fbar, f.stat_err.values / UBGEV_PER_FM, gbar, g.stat_err.values / UBGEV_PER_FM


class ForwardAmps:
    def __init__(self, path=None):
        nu, fb, fe, gb, ge = load_table(path) if path else load_table()
        self.nu, self.fb, self.fe, self.gb, self.ge = nu, fb, fe, gb, ge
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

    def _fs(self, x):
        return AB * x**2 + C4 * x**4 + C6 * x**6 + self.af[0] * x**8 + self.af[1] * x**10

    def _gs(self, x):
        return G0 * x**3 + G0B * x**5 + self.ag[0] * x**7 + self.ag[1] * x**9

    def fbar(self, nu):
        nu = np.asarray(nu, float)
        return np.where(nu <= 135.0, self._fs(nu / HBARC_MEV), self.fh(np.clip(nu, 130, 150)))

    def gbar(self, nu):
        nu = np.asarray(nu, float)
        return np.where(nu <= 135.0, self._gs(nu / HBARC_MEV), self.gh(np.clip(nu, 130, 150)))

    def err(self, nu):
        nu = np.asarray(nu, float)
        return self.ferr(np.clip(nu, 50, 150)), self.gerr(np.clip(nu, 50, 150))


if __name__ == '__main__':
    F = ForwardAmps()
    print('nu  fbar_tab  fbar_smooth  err  |  gbar_tab gbar_smooth err   [mub GeV]')
    for i, n in enumerate(F.nu):
        print(f'{n:5.0f} {F.fb[i]*UBGEV_PER_FM:8.4f} {F.fbar(n)*UBGEV_PER_FM:8.4f} {F.fe[i]*UBGEV_PER_FM:6.4f} | '
              f'{F.gb[i]*UBGEV_PER_FM:8.4f} {F.gbar(n)*UBGEV_PER_FM:8.4f} {F.ge[i]*UBGEV_PER_FM:6.4f}')
    print('remainder coefs f:', F.af, ' g:', F.ag)
    for n in (59, 69, 79, 89, 99, 109, 117, 127, 135, 143, 149, 150):
        print(n, f'{F.fbar(n)*UBGEV_PER_FM:.4f}', f'{F.gbar(n)*UBGEV_PER_FM:.4f}')
