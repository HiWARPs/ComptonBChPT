# Fit driver for the NNLO BChPT proton-Compton code (`ampl.exe`), version 2

The driver is an optional mode of the program (`ampl.f90`). If the namelist `&fitdriver` is
absent from `in.dat`, `ampl.exe` runs its original main loop.  **Version 2** (`work/code_v2`,
see below) adds bug fixes, robust evaluation at 0/180 deg and at low energy, an
order switch, the gamma-N-Delta vertex corrections of LMP2015 (on by default), optional
physical pion/neutron masses in the piN loops, optional running Delta parameters and a PDG-2024
parameter set, all controlled by the optional namelist `&v2` (keys below).  With
`dvertex=f` (and `v2="iosmode=1,loopsafe=f"`) the amplitude reproduces the code of LMP2015 to <1e-6.

## 1. Build

```
cd work/code
make -f makefile.gfortran-native      # -> execs/ampl.exe   (primary build, arm64, LoopTools 2.16)
make -f makefile.ifort                # -> execs/ampl_ifort.exe (optional x86_64/Rosetta cross-check)
```

## 2. Quick start

```
cd work/code/driver
./run_points.sh test_points.txt                         # NNLO prediction, writes test_points.out
./run_points.sh test_points.txt out.txt alpha=10.65 beta=3.15 gMM=2.2 selftest=t
python3 compton_qf.py test_points.out --demo-fit        # needs numpy (+scipy for the fit); use ../../../.venv
python3 test_quadratic_exactness.py                     # end-to-end test vs. the original main program
```

`run_points.sh POINTS [OUT] [key=value ...]` keys: `alpha beta` (1e-4 fm^3), `gEE gMM gEM gME`
(1e-4 fm^4) = gamma_E1E1, gamma_M1M1, gamma_E1M2, gamma_M1E2; `selftest=t`; `onlylex=t`;
`ref="refalphaE1=..,refbetaM1=..,refgamEE=..,refgamMM=..,refgamEM=..,refgamME=.."`; `exe=`; `keep=1`.
Unspecified polarisabilities are set to the reference values `Pref` (= static polarisabilities of
the base amplitude, i.e. the pure prediction at the chosen order).

v2 keys (namelist `&v2` of `ampl.exe`):

| key | values (default first) | meaning |
|---|---|---|
| `order=` | `3.5`, `3`, `2` | chiral order of the base amplitude (LMP2015 Table 1, low-energy counting): `2` = O(p^2) nucleon Born graphs with charge couplings (kappa=0); `3` = O(p^3): + anomalous magnetic moments in the Born graphs, pi0 (WZW) pole, piN loops; `3.5` = O(p^7/2) = NNLO: + Delta pole (s+u) and piDelta loops (the content of the 2015 code). `Pref` follows the order (0 / piN-loop / NNLO values). For Bayesian truncation-error estimates run the same points with order=2,3,3.5. |
| `dvertex=` | `t`, `f` | imaginary parts of the one-loop gamma-N-Delta vertex (LMP2015 Fig. 4; NLO in the Delta region; zero below the pion threshold). `f` = 2015 code. |
| `drunning=` | `f`, `t`, `0`, `1`, `2`, `3` | running Re gM(s), Re gE(s) (dispersion relation of the vertex Im parts), M_Delta(s), 1+Re tau_R(s) (LMP2015 App. D; `work/notes/delta_running2.md`). `t` = `1`: once subtracted at the pole (LMP2015 Eq. D5; changes the static polarisabilities, beta_M1 3.87 -> 6.12). `2`: doubly subtracted (pole and static point s0 = M_N^2): alpha_E1, beta_M1, gamma_E1M2, gamma_M1E2 unchanged, gamma_E1E1 -0.05, gamma_M1M1 -0.38 (slopes at s0). `3`: also vanishing slopes at s0, all six static polarisabilities unchanged (exploratory). Integers set `&v2 irun`; `Pref` follows (numerical extraction). |
| `consts=` | `LMP15`, `PDG24` | parameter set; PDG24: gA=1.2754, f_pi=92.07 MeV, m_pi+-, m_pi0, M_p, kappa_p (CODATA 2022), hbar c, h_A refitted to Gamma_Delta=117 MeV (2.873); gM, gE unchanged. |
| `v2="..."` | | further `&v2` entries, comma separated, e.g. `v2="pimassphys=t,nmassphys=t"`: `pimassphys` (pi0 p loops with m_pi0), `nmassphys` (pi+ n loops with the neutron mass: physical cusp at 151.4 MeV), `hA=` (pi N Delta coupling), `vxmode=` (2 Watson [default], 0 PV2006 on-shell mapping, 1 PV2006 with M_Delta->W), `pi0ff=t` (pi0 gamma gamma transition form factor of the proton triangle loop on the pi0-pole vertex, F(t) = -2M^2 C0(0,0,t;M^2,M^2,M^2), F(0)=1; static polarisabilities unchanged), `loopsafe=f` (raw loop evaluation), `elex=` (15 MeV), `tausafe=` (4e-4 M^2), `iosmode=1` (one-sided piDelta subtraction of the 2015 code, diagnostics). |

Runtime: **1.4 ms per point** (gfortran -O3, Apple M3 Max; 2000 random points in 2.9 s including
start-up; ~3 ms/point for the ifort build under Rosetta).  Points within the "Gram-unsafe" zones
near 0/180 deg cost 4x (see Sect. 6).

## 3. Input file (points)

One point per line, `#` or `!` starts a comment line:

```
id   E_lab[MeV]   theta[deg]   OBS   [LAB|CM]
```

| OBS     | meaning (proton target)                                                                 | units  |
|---------|-----------------------------------------------------------------------------------------|--------|
| `CS`    | unpolarised dσ/dΩ in the **lab frame** (dΩ_lab)                                          | nb/sr  |
| `CSCM`  | unpolarised dσ/dΩ in the c.m. frame                                                      | nb/sr  |
| `SPAR`  | lab dσ/dΩ for a 100 % linearly polarised beam, polarisation **in** the scattering plane   | nb/sr  |
| `SPERP` | same, polarisation **perpendicular** to the scattering plane                            | nb/sr  |
| `S3`    | Σ3 = (dσ∥ − dσ⊥)/(dσ∥ + dσ⊥)                                                             | 1      |
| `S2X`   | Σ2x = (σ_R − σ_L)/(σ_R + σ_L), target spin along +x (in plane, ⊥ beam, towards the scattered photon) | 1 |
| `S2Z`   | Σ2z = (σ_R − σ_L)/(σ_R + σ_L), target spin along +z (beam direction)                    | 1      |

R/L = photon helicity +1/−1. These are the conventions of Lensky–McGovern–Pascalutsa (EPJC 75, 604)
and Babusci et al. (PRC 58, 1013), as used by A2@MAMI (Martel 2015, Paudyal 2020, Sokhoyan 2017,
Mornacchi 2022).  Verified here: the code's Σ3, Σ2x, Σ2z reproduce the signs and shapes of the
LMP2015 curves (which are compared there with LEGS and A2 data); Σ3 → −1 at 90° for E→0 (Thomson).
Note: unpolarised dσ = (SPAR+SPERP)/2.  (LP2010, Fig. 9 / Eq. 33 use instead dσ = dσx+dσy, i.e.
their dσx,y are half of SPAR/SPERP.)

The angle column is the **lab** angle unless the 5th column is `CM`.  Asymmetries are
frame-independent ratios, but the angle must be given in the frame the data use (check each
data set: e.g. the A2 Σ2x data were compared with c.m. angle in the published LMP2015 Fig. 17 but
with the lab angle in the draft figure).  Energies are always the lab photon energy.

## 4. Output file

Header lines (`#`) contain the reference values `Pref`, the analytic NNLO values with their
πN/Δ-pole/πΔ breakdown, the hard-coded values of `ampl.f90`, `Puser`, flags and constants.
Each data line:

```
id  obs  E_lab  theta_lab  theta_cm  value(Puser)  flag  c_num[28]  c_den[28]
```

The observable is an **exact** ratio of quadratic forms in the shifts d_k = P_k − Pref_k:

```
O(P) = N(d)/D(d),   N(d) = Σ_{a≤b} c_num[ab] d_a d_b,   d_0 = 1,
(a,b) ordered (0,0),(0,1),…,(0,6),(1,1),(1,2),…,(6,6)
P = (αE1, βM1 [1e-4 fm^3];  γE1E1, γM1M1, γE1M2, γM1E2 [1e-4 fm^4])
```

For `CS, CSCM, SPAR, SPERP`: D ≡ 1 (`c_den = [1,0,…,0]`).  For `S3, S2X, S2Z`: D = unpolarised
cross section (lab), N = (dσ∥−dσ⊥)/2, (σ_R−σ_L)/2 resp.  γ0 and γπ follow as
γ0 = −(γE1E1+γM1M1+γE1M2+γM1E2), γπ = γM1M1+γM1E2−γE1E1−γE1M2.
`flag=1` marks points within the Gram-singular zones near 0/180 deg where the *loop amplitudes*
were extrapolated in t from the physical side (Sect. 6); `flag=2` marks E_lab < `elex` (15 MeV),
where the loops are replaced by their static-polarisability (contact-term) form.

`compton_qf.py` reads this file; `DriverOutput.predict_array(P)` evaluates all points for any
parameter set (two matrix products), `chi2_fit()` is a minimal least-squares example with optional
floating normalisations and a Baldin constraint.

## 5. How the polarisabilities enter (physics)

The amplitude is written in the covariant, gauge-invariant basis of LP2010/LMP2015 Eq. (5–7),
T = Σ_{i=1..8} A_i(s,t) ū' O_i^{μν} u ℰ'*_μ ℰ_ν with ℰ = ε − (P·ε)/(P·q) q (every structure is
gauge invariant by construction).  The NNLO amplitude A_i^NNLO (Born + π0 pole + πN loops +
Δ pole + πΔ loops) has fixed polarisabilities P^NNLO.  Shifts are introduced as O(p^4) contact
(LEC) terms, `SUBROUTINE add_LECs` (identical code used by the main program):

```
A1 += ½(δα+δβ)(x²+y²) + ½tδβ + ½ξ²t(δγMM−δγME)      A2 += −δβ − ⅛(δα+δβ)t − ξ²(δγMM+δγEM)
A3 += ξ³δγ0 − ½ξt(δγMM+δγME)                        A4 += (1−t/4)ξ(δγME−δγMM)
A5 += ξδγMM + ⅛ξt(δγEM+δγME)                        A6 += ξδγEM − ⅛ξt(δγEM+δγME)
A7 += ½ξ(δγMM+δγEM)                                  A8 += −½ξ²δγ0
x=(s−M²)/2, y=(u−M²)/2, ξ=(x−y)/2=ν M,  δγ0 = −(δγEE+δγMM+δγEM+δγME)
```

(nucleon-mass units, δP in units α_em/M^3 or α_em/M^4).  These terms are crossing symmetric
(A1,2,8 even, A3..7 odd under s↔u) and shift **exactly** the six static (Breit-frame / L'vov)
dipole polarisabilities and nothing else at that order.  Verified (i) analytically: the ω³ terms
reproduce the Breit-frame amplitudes of LMP2015 Eq. (16); (ii) symbolically via the L'vov
amplitudes (LMP2015 Eqs. A2, B3): δα→α, δβ→β, δγ→γ with one common normalisation, no mixing;
(iii) numerically: ∂Σ2x/∂γ, ∂Σ2z/∂γ, ∂Σ3/∂β tend to the LEX (Babusci) coefficients as E→0.
Higher-order (quadrupole/dispersive) polarisabilities receive the usual 1/M recoil admixtures of a
covariant contact term.  Since the helicity amplitudes are linear in δP, every cross section is a
quadratic form and every asymmetry a ratio of two quadratic forms – exactly.  Exactness tests:
`selftest=t` (direct full-amplitude evaluation vs. quadratic form, ≤2e-10, limited by LoopTools
jitter) and `test_quadratic_exactness.py` (original main program vs. driver for random and very
large shifts, ≤6e-8 = 8-digit output precision of the main program).

**Reference values.**  `Pref` defaults to the *analytic* NNLO values (LMP2015 App. C) evaluated
with the constants actually used (computed in Fortran, `nnlo_static_pols`):
αE1=11.2464, βM1=3.8707, γE1E1=−3.3351, γM1M1=2.9427, γE1M2=0.1787, γM1E2=1.1469
(γ0=−0.933, γπ=7.246).  v2: for `order=3` / `2` the piN-loop values / zero; with physical masses
or running Delta parameters the analytic values plus a numerically extracted shift
(`piN_static_numeric`, validated to 1e-3 against the analytic formulas).  The main program (no
`&fitdriver`) uses the same references (the 2015 code had them hard-coded).  The amplitude's own static polarisabilities agree with these; the hard-coded constants of the main program (11.2461, 3.87054, −3.34375, 2.95934,
0.185504, 1.14328) differ by ≤0.017 (γM1M1), so with `Pref` = analytic the labels P are the true
total polarisabilities of the amplitude.  With `onlylex=t` the base amplitude is Born + π0 pole
and `Pref` = 0 (LEX-type calculation).

## 6. Numerical limits (read before fitting)

* **Low energies (v2).**  Below `elex` = 15 MeV the loop amplitudes (which LoopTools cannot
  evaluate there: 1.8e6 nb/sr at 0.5 MeV in the 2015 code) are replaced by contact terms carrying their
  analytic static dipole polarisabilities; the neglected O(ω⁴) terms are ≤2e-5 of dσ at 15 MeV.
  Above 15 MeV the raw loops are accurate to ≲1e-5 (gfortran vs ifort ≤4e-7 at 20 MeV).
* **0° and 180° (v2).**  The box reduction is singular at exactly forward/backward scattering and
  noisy for |t| or |η| = (M⁴−su)/M² small.  For |t| or |η| < τ = min(4e-4 M², 0.2|t_min|) the loop
  amplitudes A_i are extrapolated (cubic, at fixed s) from four nodes on the physical side;
  θ = 0 and 180° are evaluated exactly in this way (forward Im A1, Im A3 reproduce the optical
  theorem with Born σ_tot to 1e-5).  Accuracy ~1e-5 (compiler/τ variations); the observable-level
  extrapolation of the 2015 code was off by up to 5.6e-4.
* **Δ region.**  v2 includes the γNΔ vertex corrections (default): Σ2x, Σ3 agree with the LMP2015
  curves to ≤0.002/0.005, dσ/dΩ (Fig. 14) to 0.4% (median).  The running Δ parameters of App. D
  are optional (`drunning`).
* **π thresholds.**  By default the piN loops use the charged pion mass and M_n = M_p (cusp at
  149.95 MeV); `v2="pimassphys=t,nmassphys=t"` puts the π0 p and π+ n thresholds at their physical
  positions (144.7 and 151.4 MeV) – recommended for data at 140–155 MeV (forward dσ at 150 MeV:
  0.18 → 0.91 nb/sr).  This is an isospin-breaking (higher-order) refinement which shifts the
  static polarisabilities (αE1 +0.17, βM1 −0.14 with the neutron mass).
* `proton=f` is **not** a neutron calculation (πN loops are proton-only); the heavy-baryon branch
  (`relchpt=f`) is legacy (loop integrals vanish above threshold).

## 7. Files

| file | purpose |
|---|---|
| `run_points.sh` | wrapper: builds a temporary run directory, writes `in.dat`, runs `ampl.exe` |
| `in_driver.template`, `baryon.dat` | input template (current constants) and baryon table |
| `test_points.txt`, `test_points.out` | small example (CS, S3, S2X, S2Z, SPAR/SPERP, lab/c.m., one 0° point) |
| `compton_qf.py` | reader / evaluator / minimal χ² fitter |
| `test_quadratic_exactness.py` | end-to-end exactness test against the original main program |
| `static_pols.py` | extraction of static polarisabilities from the amplitude (`&fitdriver staticmode=t`) |
| `analytic_pols.py` | LMP2015 App. C analytic static polarisabilities (any constants) |
| `../tests/run_regression.sh` | v2 regression suite (static pols, Thomson/LEX, unitarity, gfortran vs ifort, driver exactness, pi0ff, running Delta parameters) |
| `../tools/` | Born photoproduction multipoles, PWA runner, unitarity checks, PV2006 vertex functions, running Delta prototype (`running_delta.py`, schemes 1-3; `running_table.f90`; `check_unitarity_running.py`), LMP2015 figure validation (`validate_lmp2015.py`, `validate_fig14.py`) |

Other `&fitdriver` options: `staticmode=t` (dump A_i by mechanism on a small (ν,t) grid to
`static_As.out` + diagnostics: UV/Δ-state dependence, reproducibility, Thomson limit),
`ampsfile='xt.txt'` (dump A_i by mechanism at user (x,t) pairs to `amps.out`).
