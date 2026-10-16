# pi0 gamma gamma transition form factor on the pi0-pole graph: impact on the analysis (2026-10-14)

Scenario directory for the switch `&v2 pi0ff` of `work/code_v2` (commit bfe883f; test g in commit 7531928).
Physics, derivation and the proposed paper paragraph: `work/notes/pi0_formfactor.md`.  Nothing in
`work/fits/results`, `work/fits/cache` (pi0ff=f entries) or `work/pwa` is modified; the pi0ff=t driver runs are new
cache entries `work/fits/cache/pi0ff_*` and `./cache_mult/` (multipole projections).

| file | content |
|---|---|
| `pi0ff_fits.py` -> `pi0ff_fits.json`, `pi0ff_fits.md`, `pi0ff_fits.log` | (1) NNLO fits F2_d_all_thB / F1_d_all_thB with and without the form factor (truncation covariance recomputed), (2) Delta-region spin fit S_A2_LEGS_gpi, (3) chi^2 per point in the Delta region, (5) observables |
| `pi0ff_multipoles.py` (+ `bv2mult.py`, copy of `work/pwa/scenarios/bchpt_v2_multipoles.py`) -> `pi0ff_multipoles.json`, `pi0ff_multipoles.md` | (4) shift of the dynamical dipole polarisabilities (l = 1 non-Born multipoles) at 59-150 and 200-300 MeV |

Units: alpha, beta in 1e-4 fm^3, gamma in 1e-4 fm^4.  F(t) - 1 = -0.0062 at the largest |t| of the data below
threshold (TAPS, 149.9 MeV, theta_lab = 155 deg), -0.020 at 300 MeV / 180 deg.

## Headline numbers (pi0ff=t minus pi0ff=f; the pi0ff=f runs reproduce `results/*.json` exactly)

| quantity | pi0ff=f (paper) | pi0ff=t | shift |
|---|---|---|---|
| main NNLO fit F2_d_all_thB: alpha_E1 | 11.591(47) | 11.599(47) | +0.0073 |
| beta_M1 | 2.941(43) | 2.933(43) | -0.0077 |
| alpha-beta | 8.651(86) | 8.666(86) | +0.015 |
| chi^2 (263 dof) | 219.69 | 219.64 | -0.050 |
| Baldin fit F1_d_all_thB: alpha_E1, beta_M1 | 11.336(45), 2.827(43) | 11.343(45), 2.819(43) | +0.0074, -0.0076 (chi^2 -0.053) |
| same without truncation covariance (F2_d_all / F1_d_all) | | | +0.013, -0.013 (chi^2 -0.14) |
| truncation sigma_th recomputed with pi0ff=t | | | max relative change 1.5e-3 (no effect on the fit: the `_covF` variant gives the same shifts) |
| Delta-region spin fit S_A2_LEGS_gpi: gE1E1, gM1M1, gE1M2, gM1E2 | -3.094(563), 2.686(417), -1.256(767), 2.604(439) | -3.085, 2.683, -1.261, 2.604 | +0.008, -0.003, -0.005, -0.001 |
| gamma_pi, chi^2 (76 dof) | 9.640(1.63), 85.63 | 9.633(1.63), 85.59 | -0.007, -0.035 |
| Delta-region chi^2 (317 points, n = 1 / n fitted), sum over sets | | | -5.7 / -3.1 (LEGS dsigma -3.7/-2.2, SAL -2.5/-0.8, LEGS Sigma_3 +0.45; all others < 0.4) |
| dynamical polarisabilities, 59-150 MeV: max shift | | | alpha, beta: 0; gE1E1, gM1M1: 0.012 (150 MeV); gE1M2, gM1E2: 0.013 |
| 200-300 MeV | | | 0; 0.0135; 0.0140 |
| dsigma/dOmega_lab, relative change at 150 / 180 deg | | | 100 MeV +3.0e-4/+3.1e-4; 140 MeV +1.1e-3/+1.1e-3; 150 MeV +1.4e-3/+1.5e-3; 250 MeV +1.9e-3/+2.1e-3; 300 MeV +9.8e-4/+1.0e-3 (zero at 0 deg) |
| Sigma_3, change (max over angles) | | | 100 MeV 9e-5; 140 MeV 2.4e-4; 250 MeV 2e-5; 300 MeV -1.1e-4 |

Conclusion: the proton-loop form factor of the pi0 gamma gamma vertex is a < 2 % correction to the pi0-pole graph
in the kinematics of the database; it leaves the static polarisabilities unchanged exactly and shifts every extracted
quantity by two orders of magnitude less than its uncertainty.  The pointlike vertex is kept in the paper.
