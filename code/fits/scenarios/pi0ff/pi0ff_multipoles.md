# pi0ff: shift of the dynamical dipole polarisabilities (l = 1 non-Born multipoles)

Produced by `pi0ff_multipoles.py` (data in `pi0ff_multipoles.json`).  Non-Born = full - (Born + pointlike pi0 pole)
in both cases (paper convention); the shift is the l = 1 projection of [F(t) - 1] x (pi0 pole), independent of the LECs.
Units: 1e-4 fm^3 (alpha, beta), 1e-4 fm^4 (gamma).  Values at the NNLO polarisabilities.

## 59-150 MeV (below threshold)

| polarisability | max abs shift | at E [MeV] | shift at E_max | value at E_max (pi0ff=f) | relative |
|---|---|---|---|---|---|
| alphaE1 | 0.00000 | 59 | +0.00000 | 16.9924 | 0.0e+00 |
| betaM1 | 0.00000 | 59 | +0.00000 | 9.1771 | 2.7e-16 |
| gE1E1 | 0.01239 | 150 | -0.01239 | -10.8515 | 1.1e-03 |
| gM1M1 | 0.01239 | 150 | +0.01239 | 3.9897 | 3.1e-03 |
| gE1M2 | 0.01328 | 150 | -0.01328 | 0.4219 | -3.1e-02 |
| gM1E2 | 0.01328 | 150 | +0.01328 | 1.5129 | 8.8e-03 |

Check: shift from the two Born+pi0 runs vs from the two full runs agree to 1.1e-13 (LoopTools history jitter).

## 200-300 MeV (Delta region, real parts)

| polarisability | max abs shift | at E [MeV] | shift at E_max | value at E_max (pi0ff=f) | relative |
|---|---|---|---|---|---|
| alphaE1 | 0.00000 | 240 | +0.00000 | 6.2972 | 3.1e-16 |
| betaM1 | 0.00000 | 225 | -0.00000 | 15.4716 | -9.9e-17 |
| gE1E1 | 0.01352 | 260 | -0.01346 | 2.4598 | -5.5e-03 |
| gM1M1 | 0.01352 | 260 | +0.01346 | 5.0697 | 2.7e-03 |
| gE1M2 | 0.01403 | 235 | -0.01383 | 0.3966 | -3.5e-02 |
| gM1E2 | 0.01403 | 235 | +0.01383 | -1.6167 | -8.6e-03 |

Check: shift from the two Born+pi0 runs vs from the two full runs agree to 8.6e-15 (LoopTools history jitter).

## pi0-pole (pointlike) dynamical polarisabilities, for scale

| polarisability | 59 MeV | 150 MeV |
|---|---|---|
| alphaE1 | -0.0000 | 0.0000 |
| betaM1 | 0.0000 | -0.0000 |
| gE1E1 | 7.2345 | 2.7567 |
| gM1M1 | -7.2345 | -2.7567 |
| gE1M2 | 6.7001 | 2.2272 |
| gM1E2 | -6.7001 | -2.2272 |

Shift at selected energies (below threshold):

| E [MeV] | alphaE1 | betaM1 | gE1E1 | gM1M1 | gE1M2 | gM1E2 |
|---|---|---|---|---|---|---|
| 59 | +0.00000 | -0.00000 | -0.00619 | +0.00619 | -0.00711 | +0.00711 |
| 80 | -0.00000 | +0.00000 | -0.00849 | +0.00849 | -0.00955 | +0.00955 |
| 100 | -0.00000 | +0.00000 | -0.01011 | +0.01011 | -0.01118 | +0.01118 |
| 120 | +0.00000 | -0.00000 | -0.01127 | +0.01127 | -0.01229 | +0.01229 |
| 140 | -0.00000 | +0.00000 | -0.01209 | +0.01209 | -0.01302 | +0.01302 |
| 150 | +0.00000 | +0.00000 | -0.01239 | +0.01239 | -0.01328 | +0.01328 |

| E [MeV] | alphaE1 | betaM1 | gE1E1 | gM1M1 | gE1M2 | gM1E2 |
|---|---|---|---|---|---|---|
| 200 | +0.00000 | +0.00000 | -0.01327 | +0.01327 | -0.01394 | +0.01394 |
| 250 | +0.00000 | +0.00000 | -0.01351 | +0.01351 | -0.01402 | +0.01402 |
| 300 | +0.00000 | -0.00000 | -0.01346 | +0.01346 | -0.01383 | +0.01383 |
