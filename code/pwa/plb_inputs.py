"""Inputs of the original KLP18 analysis (Krupina's notebook PWA/Final_NB/MEX_analysis.nb):
energy bins, forward amplitudes Re f (average of GHP fits I and II) and Re g at the bins, data loader."""
import os as _os  # [release] paths relative to this file: _WORK = <release>/code (the former work/), _ROOT = <release>
_WORK = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..'))
_ROOT = _os.path.dirname(_WORK)
import numpy as np
import pandas as pd

BINS = [59, 69, 79, 89, 99, 109, 117, 127, 135, 143, 150]
# notebook values, mub GeV  (A1fit1, A1fit2, A3)
NB = {59: (-2.7567, -2.7648, -0.31606775),
      69: (-2.6527, -2.6637, -0.37052864),
      79: (-2.5281, -2.5425, -0.42553081),
      89: (-2.3803, -2.3986, -0.48081548),
      99: (-2.2057, -2.2283, -0.53578031),
      109: (-1.9988, -2.0260, -0.58913038),
      117: ((-1.8808 - 1.7512) / 2, (-1.9104 - 1.7831) / 2, (-0.61442061 - 0.63803654) / 2),
      127: ((-1.6077 - 1.4475) / 2, (-1.6419 - 1.4836) / 2, (-0.65900488 - 0.67570786) / 2),
      135: (-1.2654, -1.3031, -0.68524228),
      143: ((-1.0527 - 0.7885) / 2, (-1.0910 - 0.8254) / 2, (-0.68171866 - 0.64996667) / 2),
      150: (-0.3937, -0.4248, -0.52176691)}
HBARC = 197.327
UB = 1973.27
MNB, ALNB, KANB = 938.27, 1 / 137.036, 1.7928474


class PLBForward:
    """fbar, gbar (fm) at the KLP18 bins, Born parts subtracted at the bin energy as in the notebook."""
    def __init__(self):
        self.E = np.array(BINS, float)
        Mfm = MNB / HBARC
        f = np.array([(NB[b][0] + NB[b][1]) / 2 for b in BINS]) / UB
        g = np.array([NB[b][2] for b in BINS]) / UB
        self.fb = f + ALNB / Mfm
        self.gb = g + ALNB * KANB**2 * (self.E / HBARC) / (2 * Mfm**2)

    def _look(self, arr, E, n):
        """exact values at the bins; between bins (plotting only) interpolate arr/E^n linearly."""
        E = np.asarray(E, float)
        i = np.argmin(abs(E[..., None] - self.E), axis=-1)
        exact = abs(self.E[i] - E) <= 0.5
        interp = np.interp(E, self.E, arr / self.E**n) * E**n
        return np.where(exact, arr[i], interp)

    def fbar(self, E):
        return self._look(self.fb, E, 2)

    def gbar(self, E):
        return self._look(self.gb, E, 3)

    def err(self, E):
        return 0.01 / UB + 0 * np.asarray(E, float), 0.01 / UB + 0 * np.asarray(E, float)


def load_klp2018(path=(_WORK + '/data/db_klp2018.csv')):
    d = pd.read_csv(path)
    out = pd.DataFrame(dict(dataset=d.dataset.str.replace('KLP2018_', ''), obs='CS', E=d.E_lab_MeV.astype(float),
                            theta_lab=d.theta_lab_deg.astype(float), value=d.value.astype(float),
                            err=np.hypot(d.stat_err.astype(float), d.sys_err_pt.fillna(0).astype(float))))
    bins = np.array(BINS, float)
    out['bin'] = bins[np.argmin(abs(out.E.values[:, None] - bins[None, :]), axis=1)]
    return out
