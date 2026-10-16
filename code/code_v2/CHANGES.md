# ComptonBChPT code — version 2 (work/code_v2): changes with respect to v1

v1 = audited working copy `work/code` (= `<original Fortran source (SRC_MOD), not included>` + fit driver), see
`work/notes/code_audit.md` (bug labels B1–B10, S1–S9 refer to its Sect. 5).
History: `git log` in this directory (one commit per change; the first commit is the v1 state).
Regression suite: `tests/run_regression.sh [--build]` (16 s; all PASS at the final commit).
Physics validation against the published curves: `tools/validate_lmp2015.py`, `tools/validate_fig14.py`.

Default behaviour of v2 = v1 + bug fixes + robust evaluation + **γNΔ vertex corrections on**.
v1 is reproduced to < 1e-6 (driver grid 20–350 MeV, all angles away from 0/180°) with
`dvertex=f` and `v2="iosmode=1,loopsafe=f"` (test d).  All new options live in the optional
namelist `&v2` (absent ⇒ defaults); `driver/run_points.sh` passes them (`order=`, `dvertex=`,
`drunning=`, `consts=`, `v2="..."`).  Input/output formats of the driver are unchanged
(additional `#` header lines only; flag column now 0/1/2).

Units below: α, β in 1e-4 fm³, γ in 1e-4 fm⁴; dσ relative, asymmetries absolute.

---

## A. Bug fixes

### A1. πΔ loops: on-shell subtraction point (audit B1) — `Cov_1loop_As.F90:1551–1553, 2163ff`
*Problem.* Graph 4 of `CorrIsoCovOneDeltaLoop` evaluated its on-shell γNN-vertex constants at
p² = M² + 1e-4 (one-sided offset; the LoopTools C-reduction is singular at q = 0): A1(0,0) =
1.4e-5 (Born −1), Ξ- and μ-dependent — a Thomson-limit violation and residual scale dependence.
*Fix.* New `SUBROUTINE piD_onshell_consts`: symmetric average of p1p2s = M² ± ε (error O(ε²)),
results cached (kinematics independent; also saves ~40 LoopTools calls and two cache flushes per
call). `iosmode=1` restores v1 for diagnostics.
*Effect.* Ξ/μ dependence of A1(0,0) gone (spread 0 over Ξ = 0,1,10, μ = 1,4; v1: 4.9e-5).
Static polarisabilities extracted from the amplitude now equal the analytic LMP2015 App. C values:
α(πΔ) 4.4525 → **4.4253** (analytic 4.4253), γE1E1(πΔ) 0.4307 → 0.4343, γM1M1 −0.2186 → −0.2229,
γπ −0.873 → −0.882. Observables: dσ ≤ 6.5e-5 (max near 145 MeV), asymmetries ≤ 1e-5.
(The unused routines with the same offset, l. ~887/2446 of v1, are left as they are.)

### A2. Reference polarisabilities computed, not hard-coded (audit B4) — `ampl.f90:1536 (ref_pols_v2), 1563 (set_LECs_v2)`
v1 shifted the input polarisabilities of the main program by hard-coded NNLO values
(11.2461, 3.87054, −3.34375, 2.95934, 0.185504, 1.14328), off from the true values of the
amplitude by ≤ 3e-4 (α, β) and −0.0086/+0.0166/+0.0068/−0.0037 (γE1E1/γM1M1/γE1M2/γM1E2).
v2 computes them analytically (LMP2015 App. C, `nnlo_static_pols`) from the constants in use and
for the mechanisms switched on (order switch), plus numerically extracted shifts for the options
that the closed formulas do not cover (physical masses, running Δ). Default values:
**11.2464, 3.8707, −3.3351, 2.9427, 0.1787, 1.1469** (γ0 −0.933, γπ 7.246). The driver already used
them; in the main program input γ’s now mean total γ’s (the stored `execs/in.dat`, which passes
the old hard-coded γ’s, therefore runs with γ shifted by the differences above; α = 10, β = 4 there
is a user choice, not NNLO).

### A3. LoopTools UV-parameter state leak (audit B5) — `Cov_1loop_As.F90:1300ff, 2141`
`CorrIsoCovOneDeltaLoop` left Δ_UV = 1 for all later πN-loop evaluations; it now restores the
entry value (`getdelta`). Effect ≤ 4e-9.  Δ mass inside the πΔ loops taken from `baryon.dat`
(audit S1, COMMON /v2loopmass/; identical value 1232/938.27).

### A4. Theory-error band (audit B7) — `ampl.f90:2390`
`RTOT_ERR` used √((1+z)² + (2z)²); the LEX gives ∂σ/∂α ∝ 1+z², ∂σ/∂β ∝ 2z ⇒ √((1+z²)² + (2z)²).
Band unchanged forward and at 90°, ×1.43 at 133°, ×1.41 at 180° (only the printed error column).

### A5. `varyscalarPs` (audit B8) — `ampl.f90:441, 2303`
β scan written to unit 335 (was 325, overwriting the α scan); integer DO loop instead of a REAL
loop variable; no LECs for incomplete calculations.

### A6. Single-precision `acos` (audit B10) — `ampl.f90:516`
θ_cm = acos(real(cos θ)) → `dacos` (clamped). Changes main-program outputs at ≤ 1e-7.

### A7. Robust evaluation near 0°/180° and at low energy (audit B2, B3) — `ampl.f90:876 (loop_raw), 898 (loop_amps)`
* **Gram-singular zones.** LoopTools’ box reduction fails at t = 0 / η = (M⁴−su)/M² = 0 (NaN /
  1e55 in v1) and is noisy close to them. For |t| or |η| < τ = min(4e-4 M², 0.2|t_min|) the loop
  amplitudes A_i(s,t) are extrapolated (cubic, fixed s) from four nodes **on the physical side**
  (LoopTools returns garbage for t > 0 or η < 0 above threshold — tested). θ = 0 and 180° are thus
  evaluated exactly (only the forward amplitudes contribute; Σ3, Σ2x vanish there identically).
  Checks: forward Im A1, Im A3 at t = 0 reproduce ν σ_tot/4πα and ν(σ_1/2−σ_3/2)/8πα of
  tree-level photoproduction to 1e-5 (test c); stability against τ, node choice and compiler
  ~1e-5; v1’s observable-level extrapolation (driver only) was off by up to 5.6e-4 (CS 120 MeV, 180°)
  and gave Σ2z(0°) only after extrapolation.
* **Low energies.** Below `elex` = 15 MeV the loops are replaced by contact terms carrying their
  analytic static dipole polarisabilities (`add_LECs`); neglected O(ω⁴) terms ≤ 2e-5 of dσ at 15 MeV
  (v1: 0.27 % error at 5 MeV, 1.8e6 nb/sr at 0.5 MeV at 150°). NNLO vs Born+π0+LEX: 3e-11 (0.5 MeV),
  6e-6 (10 MeV), 1.4e-4 (20 MeV) — the O(ω⁴) Δ-pole/loop terms (test b).
* Driver: flag 1 = extrapolated loop amplitudes, 2 = low-energy replacement; the v1 observable
  extrapolation was removed. `loopsafe=f` restores the raw v1 evaluation.

### A8. Documented, not fixed
* `proton=f` (audit B6): Born and π0 pole are switched to the neutron, but the πN loops have no
  isospin argument (proton loops; α^πN 6.93 instead of 9.4 for the neutron) and the references are
  proton values ⇒ **not** a neutron calculation. The program now prints a warning.
* Heavy-baryon branch `relchpt=f` (audit B9): `Loops.f90` uses the unset `np` (loop integrals vanish
  above threshold), local hA = 2.7 in `piDloop`, debug output to fort.8/11; legacy only, warning.
* New finding **B11 (fixed)**: with `onlyLEX=t` the main program zeroed the polarisability LECs
  again in the completeness check (input polarisabilities ignored, except in varyscalarPs mode).
* S7: tolerant loop ends (`Emin ≤ Emax + 1e-6 Edel`).  S5: per-call `write` in `s_Delta2` removed.

## B. Order switch (`&v2 order`, `run_points.sh order=`) — `ampl.f90:1511 (set_order_v2), 1536`
LMP2015 Table 1, low-energy counting:
* `order=2`  O(p²): nucleon Born graphs, charge couplings only (κ = 0);
* `order=3`  O(p³): + anomalous magnetic moments in the Born graphs, π0 (WZW) pole, covariant πN loops;
* `order=3.5` O(p^{7/2}) = NNLO (default): + Δ pole (s + u) and πΔ loops (= v1 content).
`Pref` follows the order (0 / πN-loop values 6.931, −1.833, −3.370, −0.120, 0.548, 0.928 / NNLO);
the main program runs pure order-by-order predictions (no LECs) for order < 3.5. order=3.5 is
identical to the default (test d); with `dvertex=f` it reproduces v1 (< 1e-6).
Example (lab): dσ(79 MeV, 60°) = 13.57 / 13.68 / 11.64 nb/sr; Σ3(90 MeV, 90°) = −0.992 / −0.807 /
−0.770; Σ2z(250 MeV, 0°) = 0 / −0.944 / 0.713 for order 2 / 3 / 3.5 (at order 2 the Dirac-Born
Σ2z vanishes forward and at θ_lab = 90°, as in the Lipps–Tolhoek formula); 0.713 is the order-3.5
value without vertex corrections (0.745 with them, the default).

## C. Δ-region upgrade: γNΔ vertex corrections (LMP2015 Sect. 2.2, Fig. 4)
`As.f90:129ff (s_Delta2), 384 (gnd_vertex_im), 644 (born_p33)`; option `dvertex` (default `t`).
* **Implementation.** In the s-channel Δ pole the couplings become complex, gM(s) = gM + i Im gM(W),
  gE(s) = gE + i Im gE(W), above the πN threshold; the amplitude stays quadratic in g(s) (both
  vertices at the same s), so Im[g(s)²/D] supplies the Born–Δ interference 2k M_B Re M_Δ required by
  unitarity. Below threshold the Im parts vanish ⇒ **no change of any static polarisability or of
  any observable below 150 MeV** (verified, test d).
* **Im parts.** Three equivalent-at-the-peak constructions (`vxmode`):
  0 — literal transcription of Pascalutsa–Vanderhaeghen PRD 73 (2006) 034003 App. B (Im g1, g2 of
  the photon-on-pion and photon-on-nucleon loops), real-photon limit taken numerically in quad
  precision (1/Q² terms cancel; = 60-digit mpmath to 4e-13), on-shell mapping to gM, gE;
  1 — same with M_Δ → W in the mapping;
  **2 (default) — unitarity/Watson construction**: the code’s s-channel Δ pole has, exactly, the J=3/2
  multipole f_MM^{1+} = −a(W) G_M*(W)²/D(W) with a(W) = α ω²(W+M)²(W+M_Δ)/(8(M+M_Δ)²) and
  G_M* = gM − κ(W) gE, κ = (W−M)/(2(W+M)) (identified from the code’s partial waves to 2e-8);
  Watson’s theorem then fixes Im G_M* = M1+^B √(k M_Δ Γ/a), Im gE = (E1+^B/κ) √(k M_Δ Γ/a) with the
  LO Born multipoles of I = 3/2 (`born_p33`, Fortran port of `tools/born_multipoles.py`, which is
  verified by gauge invariance, the LETs E0+(π+n) = +26.6, E0+(π0p) = −2.37, and by reproducing the
  Im parts of the code’s πN loops in all multipoles to 1e-6).
  At W = M_Δ all agree: Im gM = 0.2166, Im gE = −2.201 (gM = 2.97, gE = −1). Off the peak Im gE
  differs by up to ~15 % between the variants; mode 2 reproduces LMP2015 best (below) and makes the
  complete amplitude satisfy P33 unitarity, Im f_MM^{1+} = k(|M_B,1/2|² + |M_B,3/2 + M_Δ|²), to 1e-9.
* **Validation vs LMP2015** (draft vector figures = published curves; `tools/validate_lmp2015.py`):

  | quantity | v1 (no vertex corr.) | v2 (default) |
  |---|---|---|
  | Σ2x vs θ_cm, 224/245/265/286/310/334 MeV, max \|Δ\| to band mid-line | 0.006/0.012/0.018/0.026/0.037/0.057 | **0.002/0.002/0.001/0.002/0.001/0.002** (all inside the bands) |
  | Σ3, same energies | ≤0.004, 0.013 at 334 MeV | ≤ 0.005 |
  | Σ2z, same energies | ≤ 0.023 | ≤ 0.016 (bands ±0.03–0.07) |
  | Im f_MM^{1+}, f_ME^{1+}, f_EE^{1+}, f_EM^{1+} (Fig. 10, red solid), max dev./max | (equal to LMP’s red dashed “without” curves to 0.1 %) | 0.1 %, 0.5 %, 0.5 %, 0.1 % |
  | dσ/dΩ_lab 150–400 MeV, 6 angles (published Fig. 14, digitised, 143 points) | median 1.8 %, max 14 % | **median 0.44 %, 90 % ≤ 1.6 %**, max 3.7 % (steep flanks) |

  With vxmode = 0 the multipoles agree to 0.3/2.7/2.4/0.5 % and Σ2x to 0.004.
  Effect on observables (v2 − v1, θ_lab 10–170°): up to 250 MeV dσ ≤ 2.3 %, Σ3 ≤ 0.004, Σ2x ≤ 0.014,
  Σ2z ≤ 0.038; up to 325 MeV dσ ≤ 6.3 %, Σ3 ≤ 0.011, Σ2x ≤ 0.047; at 400 MeV dσ up to −16 %
  (forward), Σ3 0.066, Σ2x 0.15, Σ2z 0.08. Forward Σ2z(250 MeV, 0°) 0.713 → 0.745. Δ-region data (no fit, statistical+point errors, no
  normalisation): A2 Σ2x (4 pts) χ² 2.5 → 2.5; A2 Σ2z (10) 13.9 → 12.1; LEGS Σ3 (79) 146.6 → 141.8;
  LEGS dσ (77) 2761 → 2889 (dominated by normalisation; LMP2015 note that gM = 2.97 is probably
  slightly too large).
* **Optional higher order: `drunning`** (`As.f90:490 delta_running, 556 delta_selfenergy`; default off).
  Running Re gM(s), Re gE(s) from the once-subtracted dispersion relation (LMP2015 Eq. D5) of the
  vertex Im parts (singularity-subtracted principal values, 400-node mapped quadrature to s → ∞),
  running M_Δ(s) = M_Δ + Re σ_R(s) and 1 + Re τ_R(s) from the πN and πΔ self-energies (D1–D4, Ξ = 0,
  scale M_N), applied to the s- and u-channel Δ pole (numerator mass, couplings, denominator, Z).
  Reproduces LMP2015 Fig. 5: M_Δ(s)/M_N (1.2156 at √s = 0.94 → 1.3131 at the peak) and 1 + Re τ_R
  (1.189 → 1) exactly as read off; gM(s) within 0.03, gE(s) within 0.14 (vxmode 2; 0.3 with mode 0).
  Being beyond the working order it changes the static polarisabilities strongly: Δ-pole
  β_M1 7.15 → 9.35, total (α, β, γE1E1, γM1M1, γE1M2, γM1E2) = (11.28, 6.12, −3.50, 2.54, 0.07, 1.19)
  (the reference values follow automatically, numerical extraction). For studies only.

## Physical pion masses per isospin channel (coordinator request) — `Cov_1loop_As.F90:18ff (CovOneLoop)`
Options `pimassphys` and `nmassphys` (both default `f` = v1).
* Channel assignment of the covariant O(p³) graphs (proton): photons on the charged pion
  (boxes, reducible pion-photon graphs, γγππ seagull) → π+n with m_π±; photons on the nucleon
  (reducible and box graphs) → π0p with m_π0; nucleon self-energy split with isospin weights 2:1;
  isoscalar ππ graphs → charged pions. With `nmassphys` the π+n loops use M_n for the internal
  nucleon (cusp at the physical 151.4 MeV instead of 149.95 MeV). Each channel is renormalised on
  shell with its own masses (separately gauge invariant). `pichan=2|3` keeps one channel (diagnostics).
* Checks: per-channel unitarity vs tree-level γp → π0p (m_π0) and π+n (m_π+) in all multipoles
  (1e-6–1e-13 above threshold; Im parts vanish below each threshold, e.g. π+n at 147 MeV);
  Thomson limit (A1(0,0) change < 1e-5); equal masses reproduce the single-mass code to 1.5e-8
  (2e-6 at extrapolated 0/180° points; LoopTools cache-history jitter).
* Static polarisabilities (shift of the reference values; π0 mass alone / plus neutron mass):
  Δα = +0.003 / **+0.169**, Δβ = +0.003 / **−0.142**, ΔγE1E1 = +0.013 / +0.089, ΔγM1M1 = −0.001 / −0.004,
  ΔγE1M2 = 0.000 / −0.009, ΔγM1E2 = 0.000 / −0.014 (numerical extraction, accuracy ~1e-3).
* dσ/dΩ_lab and Σ3 at 100–150 MeV (π0 mass only | π0 + neutron mass): dσ ≤ 0.6 % | 0.3–0.9 % at
  100–140 MeV, up to +3.6 %/−1.1 % at 145 MeV and +11 %/−5.6 % at 150 MeV (cusp region);
  Σ3 ≤ 0.0014 | ≤ 0.005 below 140 MeV, up to −0.10 at 150 MeV.
* Forward dσ(0°) vs GHP sum-rule points (`work/data/db_forward_obs.csv`), pulls (theory − GHP)/σ:
  145 MeV 2.168 → 2.382 nb/sr (GHP 2.770 ± 0.164; pull −3.7 → −2.4); 150 MeV 0.179 → **0.906**
  (GHP 1.116 ± 0.088; pull −10.7 → −2.4). The π0 mass alone does not help (0.183 at 150 MeV): the
  150 MeV point lies 1.4 MeV below the physical π+n threshold but above the isospin-symmetric one.
  50–135 MeV: pulls −1.2 … −3.8 unchanged (the base NNLO prediction at Pref; fits shift them).
* `hA` is now an `&v2` variable (default 2.85 as in v1; PDG24 refits it unless given).

## F. π⁰γγ transition form factor on the π⁰-pole graph (`&v2 pi0ff`, default off) — `ampl.f90` (`pi0_formfactor`, computeAs, comp_bkmAs, components_As)
Author's suggestion: the pointlike Wess–Zumino–Witten vertex of the t-channel π⁰ pole (graph 5 of the paper's
Fig. 1) is multiplied by the π⁰ → γγ transition form factor generated by the proton triangle loop
(Steinberger-type π⁰ → p p̄ → γγ; the neutron is neutral and does not contribute),
F(t) = −2 M_p² C0(0, 0, t; M_p², M_p², M_p²) (LoopTools convention, both photons real, pion momentum² = t),
normalised to F(0) = 1: the anomaly strength at t = 0 is kept, only the t dependence of the vertex is new.
Closed form (`REAL*8 FUNCTION pi0_formfactor(t)`, nucleon-mass units): t < 0: F = (4M²/(−t)) arsinh²(√(−t)/2M);
0 < t < 4M²: F = (4M²/t) arcsin²(√t/2M); F ≈ 1 + t/12M² + t²/90M⁴. Checked against LoopTools' C0 to 1e-15 at
−2 ≤ t ≤ 1 GeV² (`tools/check_pi0ff.F90`). `pi0ff=t` multiplies gg by F(t) in the three places where the pole
enters (computeAs A(8); comp_bkmAs A(3,5,6) of the HB branch; components_As Ac(8,2)); the closed-form LEX
expressions of `LEXobservables` (main-program LEX output only) keep the pointlike vertex. The driver header lists
`pi0ff`; `run_points.sh ... v2="pi0ff=t"`.
* **Static polarisabilities unchanged** (exactly: g F(t)/(t − m²) = g/(t − m²) + g [F(t) − 1]/(t − m²) and the second
  term vanishes at t = 0, i.e. the modification starts at O(ω²) relative to the pole, beyond the dipole
  polarisabilities): numerical extraction per mechanism (`staticmode`, `driver/static_pols.py`) gives the anomaly
  contribution γE1E1 = −γM1M1 = γE1M2 = −γM1E2 = 11.488 (γπ(π⁰) = −45.95) with and without the form factor,
  difference 1.5e-7; all other components bit-identical. `Pref` (analytic) is unaffected by construction.
* Size: F(t) − 1 = −0.0062 at the largest |t| of the database below threshold (TAPS 149.9 MeV, θ_lab = 155°,
  t = −0.066 GeV²), −0.0064 at 150 MeV/180°, −0.020 at 300 MeV/180°. Observables (main-fit polarisabilities):
  dσ(140 MeV, 180°) +0.11 %, dσ(300 MeV, 180°) +0.10 %, Σ2z(300, 180°) +7e-4, Σ3 ≤ 3e-7 (vanishes at 0/180°);
  forward observables unchanged (t = 0). Fit impact (`work/notes/pi0_formfactor.md`, `work/fits/scenarios/pi0ff/`):
  main NNLO fit F2_d_all_thB (truncation covariance recomputed with the form factor) Δα_E1 = +0.007, Δβ_M1 = −0.008,
  Δχ² = −0.05 (errors 0.47/0.43 unchanged); Baldin fit F1_d_all_thB +0.007/−0.008; Δ-region spin fit (S_A2_LEGS_gpi)
  γ's shift by ≤ 0.008, γπ by −0.007, χ² by −0.04; dynamical polarisabilities (l = 1 non-Born multipoles, pointlike
  Born+π⁰ reference) shift by ≤ 0.013 at 150 MeV (α(ω), β(ω) not at all); Δ-region χ² (317 points) −5.7 (n = 1).
  Negligible; the pointlike vertex is kept in the paper.
* `pi0ff=f` (default) is bit-identical to the previous state (driver outputs identical to the cached pipeline
  runs; regression suite 0 FAIL; flags and the exact 0°/180° evaluation unchanged, the pole is evaluated at the
  exact t and is not part of the loop extrapolation). Regression test g covers the switch.

## D. Constants (`&v2 consts='PDG24'`) — `ampl.f90:1468 (set_constants_v2), 1501 (set_hA_from_width)`
PDG 2024 / CODATA 2022: gA = 1.2754, f_π = 130.2/√2 = 92.07 MeV, m_π± = 139.570, m_π0 = 134.977,
M_p = 938.272 MeV, κ_p = 1.7928473, ħc = 197.327 MeV fm, α = 1/137.036, M_Δ = 1232 MeV, h_A refitted
to Γ_Δ = 117 MeV ⇒ h_A = 2.873 (LMP15: 2.85, Γ = 114.8 MeV); gM, gE unchanged.
Static polarisabilities LMP15 → PDG24 (πN / πΔ / Δ-pole shifts in brackets):
α 11.246 → **11.414** (+0.081/+0.086/0.000), β 3.871 → **3.821** (−0.021/−0.028/0.000),
γE1E1 −3.335 → −3.366, γM1M1 2.943 → 2.937, γE1M2 0.179 → 0.186, γM1E2 1.147 → 1.154,
γ0 −0.933 → −0.912, γπ 7.246 → 7.271.
Observables: dσ −0.3 % (80 MeV) … −0.7 % (140 MeV, 150°), −0.2 % at 150 MeV/90°, −3.2/−3.6 % at the
Δ peak (320 MeV, 90°/135°: larger width); asymmetries ≤ 0.003 (Σ3(90 MeV, 90°) +0.0003,
Σ2x(288, 90°) +0.0017, Σ2z(250, 0°) −0.0034).

## E. Other
* Version control: git repository in `work/code_v2` (no remote); first commit = v1 state.
* Speed: 1.57 ms/point (2000 random points, gfortran; v1 1.45 ms/point): the cached πΔ on-shell
  constants compensate the Watson Born projection and the 0/180° extrapolations.
* gfortran-16 and ifort 2022 build cleanly (`makefile.gfortran-native`, `makefile.ifort`);
  results agree to ≤ 4e-7 (regular points, ≥ 20 MeV) and ≤ 6e-6 (extrapolated 0/180° points).
* `tests/`: `run_regression.sh`, `regression.py`, `cmp_outputs.py`, `ref_bin/` (v1 executables,
  not in git). `tools/`: Python cross-checks (Born multipoles, PWA runner, unitarity, PV2006 vertex,
  running Δ, LMP2015 figure digitisers).
* Not done: the off-shell Δ P_− components of the vertex (beyond NLO); isospin-breaking masses in
  the πΔ loops and in the Δ width; a neutron version of the πN loops; COLLIER/quad-precision loop
  integrals (the extrapolation/LEX replacement makes them unnecessary for observables).
