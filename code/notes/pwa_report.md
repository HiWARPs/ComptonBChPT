# Partial-wave analysis of proton Compton scattering below the pion threshold: re-implementation of Krupina–Lensky–Pascalutsa (PLB 782 (2018) 34) and update with A2/MAMI 2022 and HIγS 2022

Working directory: `work/pwa/` (all code, results and figures). Date: 2026-10-14.
Units: scalar polarisabilities in 1e-4 fm³, spin polarisabilities in 1e-4 fm⁴. γπ never includes the π⁰ pole (the π⁰ pole is part of the exactly-treated "Born" amplitude, as in KLP18).

---

## 0. Executive summary

* **Machinery.** A new Python PWA code was written and validated against the Fortran NNLO BχPT code.
  * Born+π⁰ helicity amplitudes are computed from Feynman rules. They agree with the Fortran code to 5e-13.
  * The multipole ↔ partial-wave ↔ helicity transformations satisfy the low-energy theorem f_EE^{1±}(0) = −α/3M, and the forward sum-rule identities hold exactly.
  * Observables (dσ/dΩ in the lab and c.m. frames, dσ∥, dσ⊥, Σ3) agree with the Fortran code to all printed digits.
  * Resumming the projected BχPT amplitude reproduces the Fortran dσ/dΩ to ≤4e-5 and Σ3 to ≤1e-5 (J ≤ 9/2).
  * The KLP18 parametrisation (static + recoil terms) is the correct low-energy limit of the Fortran contact amplitude, to O(E²).
* **Reproduction of KLP18.**
  * The original 138-point database was rebuilt in two independent ways:
    * the author's notebook `datafull` (`db_klp2018.csv`);
    * digitisation of the vector graphics of PLB Fig. 5. Both agree to 1e-4 nb/sr.
  * The original analysis notebook (`PWA/Final_NB/MEX_analysis.nb`) was read to fix the exact procedure:
    * data are evaluated at the bin energies;
    * the GHP forward amplitudes are taken at the nearest tabulated energy;
    * a Monte-Carlo sweep is followed by "stage-2" minimisation with the residual functions fixed.
  * The published Born and BχPT curves of PLB Fig. 5 are reproduced to 1e-4.
  * **Fit 1 is reproduced:** α_E1 = 12.3, β_M1 = 1.7, γE1E1 = −3.0, γM1M1 = 1.4, γM1E2 = 2.5, χ²/N = 1.357. KLP18 has 12.2, 1.8, −3.1, 1.6, 2.5 and 1.35. The stage-2 errors agree as well. The reproduced Fit-1 angular distributions lie within 0.1–3 % of the published curves.
  * **The "refined-database" fits are reproduced qualitatively but not to the last digit.**

    | fit | α−β (this work) | α−β (KLP18) |
    |---|---|---|
    | Fit 1′ | 8.3–8.6 | 7.2 |
    | Fit 1″ | 7.4–7.7 | 6.4 |

    * The difference is inside the intrinsic seed-to-seed scatter of the KLP18 Monte-Carlo procedure (0.3–0.8 in α−β).
    * The notebook's own stored Fit 1′ result (α−β = 8.62) coincides with ours.
    * The published uncertainties on α, β and α−β carry a spurious factor cos(π/4) (see §3.4).
* **Methodological findings.**
  * In the KLP18 parametrisation the static polarisabilities are **exactly degenerate** with the per-bin residual functions. They are fixed only by the "natural-size" regularisation of the residuals.
  * The data do not constrain f^R_4, f^R_5 or f^R_6 at all: their posterior equals the prior.
  * The NNLO BχPT residual functions are f^R_4 ≈ 40–50 and f^R_6 ≈ 4.5. These are 4–5 times the KLP18 widths (10, 1).
  * Consequences: γM1M1 and γπ are regularisation-driven. With BχPT-centred priors, γM1M1 moves from 0.9 to 2.6 and γπ from 6.8 to 9.4 at identical χ².
  * The forward amplitudes are essential. Without them the error on α−β grows from 1.0 to 4.5.
* **Updated database.** 254 points at ≤ 150 MeV: 142 pre-2015 points plus 112 new ones (A2 2022 dσ 60 + Σ3 36; HIγS 2022 16).
* **MAIN fit (user decision of 2026-10-14).** All data, floating normalisations, and the GHP forward amplitudes Re f, Re g as χ² data at the bin centres. Their correlated uncertainties enter through two nuisance parameters n_f, n_g, which also move α+β and γ0. Result (χ²/N = 230.6/254 = 0.91):

  | quantity | value | stat | method | total |
  |---|---|---|---|---|
  | **α_E1 − β_M1** | **9.90** | 1.08 | 0.97 | **1.45** |
  | **β_M1** | **2.14** | 0.53 | 0.49 | **0.73** |
  | **α_E1** | **12.04** | 0.56 | 0.49 | **0.75** |
  | **γπ** (non-pole) | **6.43** | 1.67 | 2.70 | **3.17** |
  | α_E1 + β_M1 (posterior) | 14.18 | 0.17 | | |
  | γ0 (posterior) | −0.866 | 0.098 | | |

  "stat" includes the normalisations, the forward-point nuisances and the residual-function priors; "method" is the quadrature sum of the variants in §6.2. The spin polarisabilities are γE1E1 = −2.78(82), γM1M1 = 0.89(37) and γM1E2 = 2.75(69), with γE1M2 = 0 fixed.

  Database variants, all in the MAIN configuration:

  | variant | α−β | β_M1 |
  |---|---|---|
  | 3σ rule | 8.84(1.11) | 2.67 |
  | McGovern pruning | 8.03(1.15) | 3.07 |
  | modern sets only | 9.15(1.20) | 2.52 |
  | without pre-1961 legacy sets | 9.26(1.19) | 2.46 |
  | without HIγS | 9.40(1.13) | 2.38 |
  | without TAPS | 9.79(1.16) | 2.19 |
  | A2 + HIγS only | 8.94(1.41) | 2.62(69) |
  | A2 only | 7.31(1.63) | 3.41 |
  | without A2 | 11.74(1.41) | 1.12 |

  The previous main (sum rules as exact constraints) is kept as a variant: 9.91(99), β 2.04, γπ 6.82. It differs from MAIN by < 0.1 in α−β and β.
* **Outliers now.**
  * In the MAIN fit only two points exceed 3σ: Oxley 60 MeV/150° (−3.8σ) and TAPS 108.4 MeV/133° (−3.7σ).
  * The TAPS point has lost its leverage: removing it changes α−β by −0.27, compared with ≈ −1 in the 2018 database.
  * The largest remaining sensitivity comes from the 1958 Oxley set at 60 MeV. Removing all 4 points lowers α−β by 1.27; the 150° point alone by 0.77. Goldansky 1960 (55 MeV) pulls the other way (+0.66).
  * The A2, HIγS and TAPS sets are mutually compatible:

    | set | χ²/N | normalisation |
    |---|---|---|
    | A2 dσ | 50.4/60 | 0.996(13) |
    | A2 Σ3 | 25.2/36 | 1.013(20) |
    | HIγS | 15.3/16 | 0.980(19), 0.997(19) |
    | TAPS | 67.4/55 | 1.012(14) |

  * Parameter-goodness-of-fit Δχ²: TAPS-vs-rest 2.4, HIγS-vs-rest 1.9, A2-vs-rest 11.2, new-vs-old 8.4. The A2 tension is with the old non-TAPS data, which prefer α−β ≈ 11.7 against A2 ≈ 7.3.
* **Physics conclusion.** The new data shrink the database-induced spread of β_M1 from 1.4–3.5 (2018-type database, same configuration) to about 2.1–3.1. They favour β_M1 ≈ 2.4–2.7, between the DR (1.8) and χPT (3.2) extractions and consistent with PDG 2026 (2.31 ± 0.29). The DR-vs-χPT tension is reduced but not removed by the PWA. What remains is driven by:
  * pre-1961 data at 55–60 MeV;
  * the higher-order energy dependence of the M1 and spin multipoles, which the data cannot constrain.

  Precision data at E_γ ≈ 50–80 MeV would remove the first source.

---

## 1. Method and implementation

### 1.1 Amplitudes and multipoles (`compton.py`)

**Helicity conventions.**
* The c.m. helicity amplitudes are T_{H′H}, with H = λγ − λN, normalised as Φ = T/(8π√s). This follows HMP16 Eqs. (3.8)–(3.10).
* The internal array M[i,j,h′,h] uses exactly the index order and phases of `ampl.f90`. The mapping to Φ1…Φ6 was verified to satisfy all P,T relations of HMP16 Eq. (3.8) to 1e-16.

**Born and π⁰ pole.**
* The Born amplitude (nucleon s- and u-poles, Dirac+Pauli couplings) and the π⁰ t-channel pole are computed from Feynman rules with explicit Dirac algebra.
* Constants are those of the Fortran code: M = 938.27 MeV, κ = 1.793, g_A = 1.27, f_π = 92.21 MeV, m_π⁰ = 134.98 MeV, α = 1/137.036.

**Partial waves and multipoles.**
* Partial waves use T = Σ(2J+1) T^J d^J_{HH′}(θ), with the Edmonds d-functions of the Fortran routine `dfun`.
* Multipoles follow HMP16 Eqs. (3.19) and (3.21), identical to `pwa_gaga` in `ampl.f90`. The inverse is obtained by exact matrix inversion per J.

**Observables.**
* The unpolarised cross section is dσ/dΩ_cm = ½|Φ1|²+½|Φ2|²+|Φ3|²+|Φ4|²+½|Φ5|²+½|Φ6|², and dσ/dΩ_lab = (ν′/ν)² (s/M²) dσ/dΩ_cm.
* The beam asymmetry is Σ3 = (dσ∥ − dσ⊥)/(dσ∥ + dσ⊥), with Σ3·dσ = −Re[(Φ1+Φ5)Φ3* + (Φ2−Φ6)Φ4*].
* dσ∥ and dσ⊥ are the cross sections for 100 % linear polarisation in and out of the scattering plane, as in HIγS 2022.

### 1.2 PWA model (`pwa_model.py`)

The model follows KLP18 Eq. (nonBorn_multipoles_expansion).

**Multipoles.**
* Born+π⁰ is treated exactly (all J).
* The six non-Born l = 1 multipoles are written as E²(M/√s)[static + recoil + (E/M)² f^R_i] (E³ for EM/ME), with residual functions f^R_i constant within each energy bin.
* The l = 2 multipoles f^{2±}_{EE,MM} are fixed to NNLO BχPT (projected from the Fortran code with an endpoint-avoiding Gauss rule, §2) or set to zero. Our BχPT l = 1,2 multipoles agree with the numbers in Krupina's notebook to ≤ 0.1 % (l = 2: ≤ 1 %).

**Forward sum rules.**
* They are either exact linear constraints eliminating f^R_2 and f^R_3 at each energy (KLP18), or χ² points with correlated nuisances (the MAIN treatment of the update, §6; first study in §6.8.6).
* α+β = 14.0 and γ0 = −0.929 eliminate α+β and γM1M1.
* Global parameters are α−β, γE1E1, γM1E2 (and γE1M2, fixed to 0 as in KLP18 Fit 1).

**Fit.** The model is linear in all parameters at the amplitude level: Φ = u + V·p (exactly affine). The χ² is minimised with an analytic Jacobian (Levenberg–Marquardt).

**Regularisation of the residual functions (§4).**
* Gaussian priors are used: mean 0 and widths (10, 10, 1, 1) for (f^R_1, f^R_4, f^R_5, f^R_6). These are exactly the KLP18 Monte-Carlo widths (footnote 2).
* This "MAP" fit replaces the KLP18 random sweep, which is also implemented (`reproduce_plb.mc_fit`).

**Errors.**
* "full": the inverse Hessian of χ² + priors, i.e. marginalised over residuals and normalisations.
* "stage-2": residuals and normalisations held fixed (the KLP18 convention).

**Options.**
* Floating normalisations per data set: factor 1 + σ_k n_k, with penalty n_k².
* Prior centring on the BχPT residual functions.
* Prior width scaling and a smoothness prior.
* Data evaluated at their actual energy, with bin-constant f^R; this is the default for the update, while KLP18 used the bin energy.
* Forward amplitudes as χ² data with correlated nuisances.

### 1.3 Deviations from KLP18 (update only; the reproduction follows KLP18)

1. **Energies.** Data are evaluated at their actual energies, not at the bin centre. The effect on the main fit is −0.46 in α−β relative to bin-centre evaluation.
2. **Sum-rule inputs.**
   * A smooth representation of the GHP table is used: the low-energy sum-rule series (α+β = 14.0, c4 = 6.04, c6 = 4.39; γ0 = −0.929, γ̄0 = 0.484) plus two fitted higher terms, with PCHIP in the threshold region 130–150 MeV. It reproduces the table within its rounding (`sumrules.py`).
   * KLP18 instead used the table value at the nearest tabulated energy, e.g. 60 MeV for the 59 MeV bin.
3. **Statistics.** A MAP fit with KLP18's widths as priors replaces the random sweep. Errors are given both marginalised ("full") and stage-2.
4. **Database.**
   * Primary-source values from `db_old.csv` are used where available.
   * For OdL the errors are stat ⊕ 5 % point-to-point (KLP18 Fit 1 used stat only; this variant is `old_raw`/`all_raw`).
   * Hallin Table II (13 points, 130–150 MeV) replaces KLP18's mixture of Table I at 149 MeV and Table II.
   * Federspiel adds the 55 and 65.8 MeV points (8 instead of 4).
   * Floating normalisations are optional.

---

## 2. Validation (`validate_machinery.py` → `pwa/validation/validation_log.txt`)

| test | result |
|---|---|
| Born and π⁰ helicity amplitudes, Python (Feynman rules) vs Fortran (A_i → helicity), 10–150 MeV, 64 angles | max rel. dev. 5e-13 |
| LET: Born multipoles at E = 0.05 MeV (units α/M) | f_EE^{1−} = −0.33328, f_EE^{1+} = −0.33333; other multipoles ≤ 1.5e-4 (O(ω/M)) |
| Born+π⁰ partial-wave resummation, J ≤ 21/2, 59–150 MeV | ≤ 1.5e-6 |
| Born exact + projected BχPT non-Born, J ≤ 9/2, vs direct Fortran (20–160° c.m.) | dσ ≤ 4e-5, Σ3 ≤ 1e-5 |
| Truncation: KLP18 ten multipoles (l = 1 and EE/MM l = 2) vs full BχPT | 59 / 109 / 150 MeV: dσ 2e-5 / 4e-4 / 1.7e-3; Σ3 3e-5 / 7e-4 / 5e-3. With l = 1 only: up to 2.3 % (dσ), 0.026 (Σ3) |
| Observables (CS_lab, CS_cm, Σ3, dσ∥, dσ⊥): Python formulas vs Fortran fit-driver | identical to all printed digits (7 kinematics) |
| KLP18 parametrisation vs Fortran contact amplitude (one polarisability at a time) | implied f^R stays O(1)·pol as E→0, i.e. static + recoil terms correct to O(E²) |
| Forward amplitudes: (√s/2)(Φ1 ± Φ5)(θ=0) vs KLP18 multipole sums (Born+π⁰ and full BχPT) | identical (1e-7 rel.); g_Born = −ακ²ν/2M² exactly (sign of g as in GHP; Σ2z(0°) = −2Re(fg*)/(|f|²+|g|²) consistent) |
| Born+π⁰ and BχPT dσ/dΩ_lab vs curves digitised from PLB Fig. 5 (all 11 bins) | ≤ 1e-4 (≤ 3e-3 for BχPT at forward angles, 127–135 MeV); the panel labelled "149 MeV" is computed at 150 MeV |

**Issues found in the Fortran code.**
* **Loop amplitudes near 0° and 180°.** The NNLO loop amplitudes (πN and πΔ) are numerically unstable within about 2–3° of θ_cm = 0° (πN) and 180° (both). Values blow up by up to 1e9×. The check of the 2015 code found the same.
* **Consequences.** Partial-wave projections must use endpoint-avoiding quadrature, and the BχPT dσ/dΩ at 180° (Zieger) must be taken from the resummed partial waves (done).
* **`onlyLEX=t` does not work in the main program and `components_As`.** The NNLO-completeness check that follows the onlyLEX block zeroes all LEC shifts, so onlyLEX returns pure Born+π⁰. The fit driver `run_points` is not affected. LEX tests were therefore done with the full NNLO run and shifted polarisabilities.

**Fortran usage.** My copy is `work/pwa/fortran/` (a copy of the 2015 code used for the validation; the original is not modified). One addition: `run_helpoints` (namelist `&fitdriver helmode=t`) dumps the c.m. helicity amplitudes per mechanism. The interface is `fortran_runner.py`, with results cached.

---

## 3. Reproduction of KLP18 (`reproduce_plb.py`, `results/plb_reproduction.json`)

### 3.1 Database (KLP18 Table I)

**Sources.**
* `work/data/db_klp2018.csv`, which is the author's `datafull` list: 138 points, OdL 55, Oxley 4, Hyman 12, Goldansky 5, Pugh 16, Federspiel 4, MacGibbon 18, Baranov 7, Zieger 2, Bernardini 2, Hallin 13.
* Independently, my vector digitisation of PLB Fig. 5 (`digitize_plb.py` → `results/plb_fig5_points.csv`) recovers the same 138 points to ≤ 1e-3 nb/sr in value and error.

**Notable details.**
* Energies are the data-group labels (e.g. Baranov 80.9 → 79, Hyman 87 → 89).
* Hallin = Table II (9 points) + Table I at 149 MeV (4 points), converted c.m. → lab with stat ⊕ sys errors.
* OdL errors are statistical only.

### 3.2 Procedure (from `PWA/Final_NB/MEX_analysis.nb`)

**Kinematics.** All data of a bin are evaluated at the bin energy (59, 69, 79, 89, 99, 109, 117, 127, 135, 143, 150 MeV). Groups go to the nearest bin; Bernardini 139 MeV goes to 135.

**Sum-rule inputs.**
* Re f is the average of GHP fits I and II; Re g is from GHP.
* Values are taken at the tabulated energies 60, 70, 80, 90, 100, 110 MeV for the bins 59…109 ("although given at 80 MeV we use them for 79 MeV").
* Averages of neighbouring tabulated energies are used for 117, 127 and 143 MeV.
* The constraints solve for f^R_2 and f^R_3 (notebook `tempcons1v6`/`tempcons2v6`).

**Monte Carlo sweep and stage 2.**
* The MC sweep uses 500 000 draws in each of 4 parallel kernels, over 3 globals and 44 residuals, keeping the best χ².
* `FindMinimum` then varies the 3 globals with the residuals fixed.
* The "Hessian" is ½∂²χ² of this stage.

**Error propagation (the cos(π/4) factor).**
* The notebook sets `erraE1 = unc(α−β)/2·Cos[Pi/4]`.
* Since α = (14.0 + (α−β))/2, the correct value is unc(α−β)/2.
* The published σ(α), σ(β) are therefore too small by √2 even within the stage-2 definition, and Table III's σ(α−β) is also below the stage-2 value. For Fit 1, stage-2 gives σ(α−β) = 0.82, compared with ±0.4 in the PLB.

### 3.3 Results

The MAP fit uses the KLP18 priors. "MC" is the KLP18 procedure (3–6 seeds × 4·10⁵ draws; mean and spread in brackets). For PLB the published values are given; the χ²/N column in the PLB tables is consistent with χ²/138 for all fits (e.g. Fit 1′: 136/138 = 0.99).

| Fit | α_E1 | β_M1 | γE1E1 | γM1M1 | γM1E2 | α−β | γπ | χ²/N |
|---|---|---|---|---|---|---|---|---|
| **1 (MAP)** | 12.32 | 1.68 | −2.97 | 1.43 | 2.46 | 10.63 | 6.87 | 1.357 |
| 1 stage-2 err | 0.41 | 0.41 | 0.65 | 0.30 | 0.71 | 0.82 | 1.29 | |
| 1 full err | 0.51 | 0.51 | 0.96 | 0.41 | 0.85 | 1.03 | 1.92 | |
| 1 (MC) | 12.15 (0.21) | 1.85 (0.21) | −2.58 (0.56) | 1.29 (0.17) | 2.21 (0.40) | 10.30 (0.41) | 6.09 (1.12) | 1.35–1.37 |
| 1 PLB | 12.2(3) | 1.8(3) | −3.1(7) | 1.6(3) | 2.5(7) | 10.5(4) | 7.2(1.0) | 1.35 |
| 0 (γE1M2 free; MAP) | 12.37 | 1.63 | −2.29 | 1.74 | 2.44 (γE1M2 −0.97(3.85)) | 10.74 | 7.44 | 1.357 |
| 0 PLB | 12.2 | 1.8 | −1.6 | 1.8 | 2.0 (γE1M2 −1.3(3.7)) | | | 1.35 |
| 1_3σ (PLB's 4 points; MAP / MC) | 12.30 / 11.98 | 1.70 / 2.02 | −3.03 / −2.14 | 1.47 / 1.16 | 2.48 / 1.91 | 10.59 / 9.97 | 6.98 / 5.21 | 0.977 |
| 1_3σ PLB | 11.8 | 2.2 | −2.7 | 1.5 | 2.2 | | | 0.97 |
| 1′ (MAP / MC) | 11.13 / 11.31 | 2.87 / 2.69 | −1.65 / −2.15 | 1.05 / 1.20 | 1.52 / 1.88 | 8.25 / 8.62 (0.44) | 4.23 / 5.24 | 1.015 (0.985 per 138) |
| 1′ PLB | 10.6 | 3.4 | −1.0 | 1.0 | 1.0 | 7.2(6) | 3.0(1.1) | 0.99 |
| 1″ (MAP / MC) | 10.87 / 10.70 (0.16) | 3.13 / 3.30 | −2.04 / −1.85 | 0.91 / 0.89 | 2.05 / 1.89 | 7.74 / 7.40 (0.32) | 5.00 / 4.63 | 0.688 (0.634 per 138) |
| 1″ PLB | 10.2(4) | 3.8(4) | −1.2(8) | 0.6(3) | 1.6(8) | 6.4(6) | 3.5(1.2) | 0.62 |
| 2 (l=2 → 0; MAP / MC) | 12.26 / 12.20 | 1.74 / 1.80 | −3.05 | 1.29 | 2.69 | 10.51 | 7.02 | 1.363 |
| 2 PLB | 11.7 | 2.3 | −2.6 | 1.1 | 2.4 | | | 1.35 |
| 2″ (MAP / MC) | 10.82 / 10.68 | 3.18 / 3.32 | −2.13 | 0.75 | 2.30 | 7.64 | 5.18 | 0.697 |
| 2″ PLB | 10.8(4) | 3.2(4) | −1.9(8) | 0.7(3) | 2.2(8) | | | 0.69 |

### 3.4 Diagnosis

1. **Fit 1, Fit 0 and Fit 2″ are reproduced** to ≤ 0.2 in α and β and ≤ 0.5 in the γ's. The angular distributions match the published Fig. 5 curves:

   | bin [MeV] | 59 | 79 | 109 | 135 | 150 |
   |---|---|---|---|---|---|
   | max deviation, Fit 1 (%) | 0.1 | 0.6 | 1.0 | 2.9 | 2.1 |

   For Fit 1″ the deviation is ≤ 3.4 % (see `figs/fig_klp18_reproduction.pdf`).
2. **Fits 1′, 1″, 1_3σ and 2 are reproduced qualitatively.** The pattern is that pruning raises β by about 1.5 and lowers α−β by about 3. The published central values of 1′ and 1″ lie 1–2.5 MC-spreads below ours.
   * The MC+stage-2 procedure is intrinsically seed-dependent: the spread over seeds is 0.2–0.8 in α−β and up to 2 in γπ.
   * The notebook stores the Fit 1′ result of an earlier run (March 2017): **α−β = 8.62, γE1E1 = −2.06, γM1E2 = 1.78, χ² = 135.8**. This is equal to our MC mean (8.62) and close to our MAP value (8.25), but differs from the published 7.2.
   * **The published 1′/1″ numbers therefore come from a particular, more extreme MC run. Our result is the reproducible one.**
3. **The l = 2 multipoles are not the cause.** KLP18 found sizeable shifts of α between Fit 1 and Fit 2 (−0.5) and between Fit 1″ and Fit 2″ (+0.6). These come from MC noise, not from the l = 2 multipoles. In our deterministic fits the effect of setting l = 2 to zero is ≤ 0.06 in α, as expected from the BχPT truncation test (≤ 0.2 % in dσ).
4. **Errors.**
   * The published σ(α), σ(β) and σ(α−β) are smaller than the stage-2 errors (factor cos(π/4), §3.2).
   * The stage-2 errors themselves ignore the residual-function freedom. The full (marginalised) errors are 1.2–1.5× the stage-2 ones even with the narrow KLP18 priors.
5. **Fit 1 outliers (> 3σ)** in our reproduction:

   | data set | E [MeV] | θ | pull |
   |---|---|---|---|
   | Oxley | 60 | 150° | −4.0 |
   | OdL | 109 | 133° | −4.3 |
   | OdL | 89 | 155° | +3.7 |
   | Goldansky | 55 | 150° | +3.1 |
   | OdL | 79 | 107° | +3.08 |

   The first four are exactly KLP18's four outliers (Fig. 2 caption). The fifth is marginal.

---

## 4. Methodological findings on the KLP18 parametrisation

1. **Exact degeneracy.**
   * At each bin energy the four free residuals (f^R_1, f^R_4, f^R_5, f^R_6) span the full four-dimensional space of l = 1 multipoles allowed by the two forward constraints.
   * Any change of the static polarisabilities can be absorbed by the residuals of every bin. The data χ² is exactly invariant along four directions.
   * The static polarisabilities are therefore fixed **only** by the regularisation of the residual functions: the Gaussian widths of the KLP18 MC, or our priors.
   * The data determine the multipoles at the bin energies. The statics are an extrapolation.
2. **The data do not constrain f^R_4, f^R_5 or f^R_6.** In the main updated fit their posteriors equal the priors: 0 ± 10, 0 ± 1, 0 ± 1 (`results/residual_functions.txt`). f^R_1 is weakly constrained (±5–10).
3. **The natural-size widths are too narrow.** The NNLO BχPT multipoles imply these residual functions in the KLP18 parametrisation:

   | residual | 59–150 MeV |
   |---|---|
   | f^R_1 | 8 → 11 |
   | f^R_3 | 36 → 40 |
   | f^R_4 | 52 → 38 (Δ-driven rise of β_M1(ω)) |
   | f^R_5 | −3 → −2 |
   | f^R_6 | 4.9 → 4.4 |
   | f^R_2 | −3 → 241 (π-threshold cusp) |

   The sum rules force the fitted f^R_2 = 34–175 and f^R_3 = 50–80, so those two are far outside any "natural" width of 10. f^R_4 and f^R_6 are 4–5σ outside the KLP18 priors.

   The narrow f^R_4 prior forces the energy rise of the M1 multipole into the static recoil term −(2/3)E·γM1M1. This biases γM1M1 low: 0.6–1.6 in all KLP18-type fits, compared with about 3 in BχPT and DR.
4. **Regularisation scan** (PLB database, Fit 1 / Fit 1″; α−β, γπ):

   | priors | Fit 1 | Fit 1″ |
   |---|---|---|
   | KLP18 (0 ± 10/1) | 10.6(1.0), 6.9(1.9) | 7.7(1.3), 5.0(2.3) |
   | widths ×3 | 10.3(1.5), 6.8(3.8) | 6.9(1.8), 3.5(4.5) |
   | widths ×5 | 10.6(1.9), 7.4(5.6) | 6.7(2.2), 3.1(6.6) |
   | BχPT-centred (±10/1) | 11.6(1.0), 9.5(1.9) | 8.8(1.3), 7.8(2.3) |
   | smoothness prior only | 8.8–9.4 (±2.3), ≈0–2 (±8.5) | 5.8–6.0 (±2.8), 2–2.5 (±10) |

   **The database effect on α−β is robust:** it shifts by 3 between Fit 1 and Fit 1″ for every regularisation. **γπ and γM1M1 are not robust.**
5. **The sum-rule information is essential.**
   * With the forward amplitudes dropped (α+β and γ0 kept as Gaussian priors), α−β = 6.7 ± 4.5.
   * Imposing KLP18-size priors on f^R_2 and f^R_3 instead of the sum rules biases the fit badly: γπ = 14, and the α+β nuisance runs to +3σ.

---

## 5. Updated database (`build_update_db.py` → `results/db_update.csv`; per-beam-time A2 variant `db_update_A2beamtimes.csv`)

**Error model.**
* For χ²: error = stat ⊕ point-to-point sys.
* Normalisation: Gaussian scale σ_k per set. The A2 Σ3 scale is 5 % uniform → 2.89 % rms.
* Sets with unknown or no common scale are not floated: Oxley, Goldansky, Bernardini, Pugh, Hyman, Zieger.

The last two columns are the χ²/N and normalisation factor in the MAIN fit (§6).

| set | source | obs | N | E [MeV] | θ_lab | scale | χ²/N (MAIN) | norm (MAIN) |
|---|---|---|---|---|---|---|---|---|
| TAPS 2001 (Olmos de León) | db_old (thesis, stat ⊕ 5 % p2p) | dσ | 55 | 58.9–149.9 | 59–155 | 3 % | 67.4/55 | 1.012(14) |
| Federspiel 1991 | db_old (author legacy) | dσ | 8 | 55–70 | 60, 135 | 2 % | 2.3/8 | 0.998(17) |
| Zieger 1992 | db_old | dσ | 2 | 98, 132 | 180 | – | 1.9/2 | – |
| Baranov 1974 90° / 150° | db_old (HEPData) | dσ | 3 / 4 | 81–111 | 90, 150 | 8 % | 6.2/3, 8.6/4 | 0.989(23) |
| MacGibbon 1995 tagged / untagged | db_old | dσ | 8 / 10 | 73–145 | 90, 135 | 4 % / 2.9 % | 5.2/8, 3.0/10 | 1.007, 1.004 |
| Hallin 1993 (Table II, c.m.→lab) | db_old (HEPData) | dσ | 13 | 130–150 | 45–136 | 3.6 % | 7.1/13 | 1.041(27) |
| Oxley 1958 | db_klp2018 | dσ | 4 | 60 | 70–150 | – | **21.2/4** | – |
| Goldansky 1960 | db_klp2018 | dσ | 5 | 55 | 75–150 | – | 9.7/5 | – |
| Bernardini 1960 | db_klp2018 | dσ | 2 | 120, 139 | 134 | – | 0.9/2 | – |
| Pugh 1957 / Frisch 1967 | db_klp2018 (E from db_old for 90° groups) | dσ | 16 | 55–135 | 50–135 | – | 4.5/16 | – |
| Hyman 1959 | db_klp2018 (E from db_old) | dσ | 12 | 60–128 | 50, 90 | – | 1.8/12 | – |
| **A2 2022 (Mornacchi)** | db_new | dσ_lab | 60 | 92.2–135.3 | 35–145 | 3 % | 50.4/60 | 0.996(13) |
| **A2 2022** | db_new | Σ3 | 36 | 92.2, 108.5, 129.5 | 35–145 | 2.9 % | 25.2/36 | 1.013(20) |
| **HIγS 2022 (Li), circular** | db_new | dσ_lab | 8 | 81.3 | 55–125 | 2.24 % | 7.1/8 | 0.980(19) |
| **HIγS 2022, linear** | db_new | dσ∥ (5), dσ⊥ (3) | 8 | 83.4 | 55–125 | 2.24 % | 8.2/8 | 0.997(19) |

Totals: 254 points (142 old, 112 new). The per-beam-time A2 variant has 350 points, with separate March/July normalisations for dσ and Σ3.

**Binning.**
* Default: the KLP18 bins (59, 69, 79, 89, 99, 109, 117, 127, 135, 143, 150 MeV), with the data at their own energies.
* The new data sit at most 4.4 MeV from a bin centre: HIγS → 79; A2 92.2 → 89; 103.2 → 99; 113.6 → 117; 124.1 → 127; 135.3 → 135; Σ3 108.5 → 109; Σ3 129.5 → 127.
* Σ3 and dσ at the same energy thus constrain the same residual functions.
* Variant B re-centres the bins on the new data: 59, 69, 81, 92, 103, 113, 124, 135, 143, 150. The effect is −0.05 in α−β (MAIN).

**BχPT NNLO prediction (no fit; α = 11.2, β = 3.9)** — `results/chi2_bchpt_prediction.csv`:

| data set | χ²/N |
|---|---|
| A2 dσ | 62.5/60 |
| A2 Σ3 | 27.0/36 |
| HIγS | 13.9/8, 10.3/8 |
| TAPS | 80/55 |
| Oxley | 32/4 |
| Baranov 150° | 37/4 |
| total | 303/254 |

Born+π⁰ alone gives χ² ≈ 10⁴.

---

## 6. Updated fits

* Scripts: `main_fits.py` → `results/main_fits.json` (MAIN and its variants); the earlier exact-constraint fits come from `update_fits.py` → `results/update_fits.json`.
* Tables: `results/tables.md`, `results/tables.tex`.

**MAIN configuration.**
* All data, with floating normalisations (penalty (N−1)²/σ²).
* The GHP forward amplitudes Re f, Re g are χ² data at the 11 bin centres, using the smooth GHP representation (§1.3).
* Correlated uncertainties: the data are shifted by n_f σ_f(ν) and n_g σ_g(ν), with penalty n_f² + n_g² and an uncorrelated (rounding) error of 0.0029 μb·GeV.
* α+β = 14.0 + 0.2 n_f and γ0 = −0.929 + 0.105 n_g move with the nuisances, so Baldin and γ0 are no longer fixed inputs.
* f^R_2 and f^R_3 are free per bin (prior width 1000); KLP18 widths (10, 10, 1, 1) are used for f^R_1, f^R_4, f^R_5, f^R_6.
* γE1M2 = 0, l = 2 multipoles from BχPT, data at their actual energies, KLP18 bins.

All errors are "full": marginalised over residuals, normalisations and nuisances.

### 6.1 MAIN result and database variants (all in the MAIN configuration)

| fit | N | χ²/N | α_E1 | β_M1 | α−β | γπ | γE1E1 | γM1M1 | γM1E2 | α+β | γ0 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **MAIN** | 254 | 0.91 | **12.04(56)** | **2.14(53)** | **9.90(1.08)** | **6.43(1.67)** | −2.78(82) | 0.89(37) | 2.75(69) | 14.18(17) | −0.866(98) |
| 3σ outliers removed (Oxley 60/150°, TAPS 108.4/133°) | 252 | 0.80 | 11.51(58) | 2.67(55) | 8.84(1.11) | 5.51(1.76) | −2.30(87) | 0.85(38) | 2.36(72) | 14.19 | −0.900 |
| McGovern pruning | 243 | 0.73 | 11.10(59) | 3.07(56) | 8.03(1.15) | 5.31(1.82) | −2.19(90) | 0.82(38) | 2.30(74) | 14.16 | −0.930 |
| modern sets only (TAPS, MacGibbon, Hallin, Federspiel, A2, HIγS) | 206 | 0.86 | 11.67(62) | 2.52(59) | 9.15(1.20) | 6.14(1.78) | −2.61(88) | 0.93(38) | 2.59(74) | 14.19 | −0.912 |
| without pre-1961 legacy sets (Oxley, Goldansky, Bernardini, Pugh, Hyman) | 215 | 0.89 | 11.72(61) | 2.46(58) | 9.26(1.19) | 5.84(1.78) | −2.46(88) | 0.91(38) | 2.47(74) | 14.18 | −0.916 |
| A2 2022 + HIγS 2022 only | 112 | 0.79 | 11.56(73) | 2.62(69) | 8.94(1.41) | 7.25(2.23) | −3.18(1.11) | 1.15(51) | 2.92(88) | 14.18 | −0.895 |
| A2 2022 only | 96 | 0.74 | 10.72(84) | 3.41(80) | 7.31(1.63) | 5.88(2.51) | −2.48(1.25) | 1.15(54) | 2.25(98) | 14.14 | −0.925 |
| without HIγS 2022 | 238 | 0.91 | 11.78(59) | 2.38(56) | 9.40(1.13) | 5.99(1.72) | −2.56(85) | 0.88(37) | 2.55(72) | 14.16 | −0.873 |
| without TAPS 2001 | 199 | 0.82 | 11.98(60) | 2.19(57) | 9.79(1.16) | 7.07(1.90) | −3.11(94) | 0.95(42) | 3.01(76) | 14.16 | −0.849 |
| without A2 2022 | 158 | 0.94 | 12.85(72) | 1.12(70) | 11.74(1.41) | 7.19(2.15) | −3.14(1.06) | 1.19(49) | 2.87(1.01) | 13.97 | −0.917 |
| pruned and without HIγS | 227 | 0.70 | 10.60(63) | 3.53(60) | 7.07(1.22) | 4.39(1.94) | −1.72(96) | 0.78(40) | 1.89(78) | 14.14 | −0.949 |
| pre-2015 data only | 142 | 0.94 | 12.57(77) | 1.37(76) | 11.20(1.52) | 6.63(2.25) | −2.86(1.11) | 1.15(49) | 2.62(1.05) | 13.94 | −0.916 |
| pre-2015 data only, pruned | 131 | 0.64 | 10.49(86) | 3.49(85) | 7.00(1.70) | 3.61(2.69) | −1.33(1.34) | 0.61(55) | 1.67(1.18) | 13.98 | −0.956 |

The no-HIγS rows reproduce exactly the scenario results in `work/pwa/scenarios/no_higs.json` (`no_higs_fwd.py`): noHIGS 9.40(1.13), pruned+noHIGS 7.07(1.22).

The HIγS data raise α−β by 0.5 and lower β by 0.24. Their own preference, not meaningful alone with 16 points at one energy, is for a very small β, as in the HIγS paper's own fit (β = 0.2 ± 1.2).

### 6.2 Uncertainty budget of MAIN

`make_tables.py`, Table S; numbers also in `results/main_total_errors.json`. Shift = variant − MAIN.

| source | Δ(α−β) | Δγπ | Δβ_M1 | Δα_E1 |
|---|---|---|---|---|
| statistical (norms, forward nuisances n_f/n_g, residual priors) | 1.08 | 1.67 | 0.53 | 0.56 |
| residual priors ×3 wider | +0.40 | +0.20 | −0.21 | +0.19 |
| residual priors centred on BχPT residual functions | +0.82 | +2.60 | −0.41 | +0.42 |
| l = 2 multipoles BχPT → 0 | +0.07 | +0.23 | −0.06 | +0.01 |
| γE1M2 free (fit: −0.43) | +0.06 | +0.29 | −0.03 | +0.03 |
| binning re-centred on A2/HIγS energies | −0.05 | +0.24 | +0.04 | −0.01 |
| energy mode: bin-centre evaluation (KLP18) vs actual energies | +0.30 | +0.48 | −0.16 | +0.14 |
| A2 per beam time (4 scales) vs combined | −0.03 | −0.15 | +0.01 | −0.02 |
| uncorrelated forward-point errors vs correlated nuisances | −0.01 | +0.16 | −0.08 | −0.09 |
| method/theory total | 0.97 | 2.70 | 0.49 | 0.49 |
| **total** | **1.45** | **3.17** | **0.73** | **0.75** |

Variants not added to the budget (shift relative to MAIN):

| variant | Δ(α−β) | Δγπ | Δβ_M1 |
|---|---|---|---|
| exact sum-rule constraints (previous main) | +0.01 | +0.40 | −0.09 |
| no floating normalisations | −0.22 | −0.42 | +0.11 |
| 3σ outliers removed | −1.07 | −0.92 | +0.53 |
| McGovern pruning | −1.88 | −1.11 | +0.93 |
| modern sets only | −0.76 | −0.29 | +0.38 |
| without HIγS | −0.51 | −0.44 | +0.24 |
| without legacy sets | −0.64 | −0.59 | +0.32 |

**Sum-rule quantities (posteriors).** α+β = 14.18 ± 0.17 and γ0 = −0.866 ± 0.098, with n_f = +0.92 and n_g = +0.60. The data pull the forward amplitudes about 1σ up within the GHP band; the effect on α−β and β is < 0.1. With uncorrelated forward errors the posteriors are α+β = 14.01 ± 0.20 and γ0 = −0.928 ± 0.105.

**Result:**
* α_E1 − β_M1 = 9.9 ± 1.1 (stat) ± 1.0 (method), total ±1.45;
* β_M1 = 2.14 ± 0.53 ± 0.49 (±0.73);
* α_E1 = 12.04 ± 0.56 ± 0.49 (±0.75);
* γπ = 6.4 ± 1.7 ± 2.7 (±3.2);
* ρ(α−β, γπ) = 0.72.

**Database range** (3σ / pruned / modern / legacy-free): α−β = 8.0 … 9.9, β_M1 = 2.1 … 3.1. The extreme subsets are A2 only (7.3, β = 3.4) and everything except A2 (11.7, β = 1.1).

### 6.3 Outliers and influence (MAIN configuration)

**3σ rule.**
* Two points are flagged: Oxley {60 MeV, 150°} (pull −3.83) and TAPS {108.4 MeV, 133°} (−3.71).
* After removing them no further point exceeds 3σ; the iterated rule converges at once.
* TAPS {89, 155°}, TAPS {79, 107°} and Goldansky {55, 150°}, outliers in 2018, are now below 3σ.

**Single removals** (Table I in `tables.md`):

| removed | α−β | β_M1 | γπ |
|---|---|---|---|
| none (MAIN) | 9.90 | 2.14 | 6.43 |
| Oxley (4 points) | 8.63 | 2.78 | 5.04 |
| Oxley 60/150° only | 9.13 | 2.53 | 5.59 |
| TAPS 108.4/133° | 9.63 | 2.28 | 6.37 |
| Baranov 150° | 9.70 | 2.24 | 6.71 |
| Bernardini | 9.89 | 2.15 | 6.44 |
| Goldansky | 10.56 | 1.81 | 7.11 |
| Oxley + Goldansky | 9.29 | 2.45 | 5.78 |
| A2 Σ3 | 9.92 | 2.14 | 6.38 |

**The leverage has moved.** In 2018 the decisive point was TAPS 108 MeV/133°, which motivated KLP18's call for 109 MeV backward data. That region is now covered by A2 2022. The remaining database sensitivity sits at 55–60 MeV, in the 1958–60 Oxley and Goldansky data.

A2 Σ3 has almost no influence on α−β (+0.02). Its effect enters through β_M1(ω) at ≥ 92 MeV, where the residual functions absorb it.

### 6.4 Consistency between data sets (MAIN configuration)

**χ² per set** (§5 table):
* All new sets: χ²/N = 0.7–1.0.
* TAPS: 1.23 (5 % point-to-point errors).
* Inconsistent old sets: Oxley (5.3/point), Goldansky (1.9/point) and Baranov 150° (2.2/point).

**Normalisations** all lie within 1.5σ of 1, the largest being Hallin with 1.041(27).

**Parameter-goodness-of-fit** (data χ² only; the subset and the rest are fitted separately in the MAIN configuration):

| subset vs rest | N | Δχ² | α−β subset / rest | β subset / rest |
|---|---|---|---|---|
| TAPS | 55 / 199 | 2.4 | 10.8(2.4) / 9.8(1.2) | 1.6 / 2.2 |
| HIγS | 16 / 238 | 1.9 | (19 ± 6) / 9.4(1.1) | – / 2.4 |
| A2 | 96 / 158 | 11.2 | 7.3(1.6) / 11.7(1.4) | 3.4 / 1.1 |
| new (A2+HIγS) vs old | 112 / 142 | 8.4 | 8.9(1.4) / 11.2(1.5) | 2.6 / 1.4 |

* With ≥ 3 shared global parameters (plus residuals and nuisances), Δχ² = 11.2 corresponds to p ≈ 0.01. This is a ≈ 2.5σ tension between A2 2022 and the rest of the database, driven by α−β (β).
* Within the old data, the tension is carried mainly by the sets that disagree with A2 at backward angles: Oxley, Goldansky, Baranov 150°.
* TAPS and HIγS are fully compatible with the rest.

### 6.5 Forward points and their impact

The sign conventions were verified in §2:
* g = (√s/2M)(Φ1 − Φ5)(0°) → −ακ²ν/2M²;
* the KLP18 multipole sums for f and g hold identically in the Fortran conventions;
* GHP's Σ2z(0°) = −2Re(fg*)/(|f|²+|g|²) corresponds to the helicity-3/2 ↔ f − g assignment.

| treatment | α−β | γπ | β | α+β | γ0 | n_f, n_g | χ²/N |
|---|---|---|---|---|---|---|---|
| **MAIN** (χ² points, correlated nuisances, floating norms) | 9.90(1.08) | 6.43(1.67) | 2.14(53) | 14.18 | −0.866 | 0.92, 0.60 | 0.908 |
| χ² points, uncorrelated table errors | 9.89(1.05) | 6.58(1.68) | 2.06(54) | 14.01 | −0.928 | – | 0.903 |
| exact constraints (KLP18 treatment; previous main) | 9.91(99) | 6.82(1.60) | 2.04(50) | 14.00 | −0.929 | – | 0.923 |
| no forward points (α+β, γ0 priors only), no norms | 6.7(4.5) | 3.6(2.3) | 3.6(2.3) | 14.00 | −0.930 | – | 0.868 |
| no forward points, KLP18-size priors on f^R_2,3, no norms | 9.4(0.6) | 14.2(1.1) | 2.6(0.3) | 14.60 | −0.885 | n_f = 3.0 | 1.127 |

**Impact of the forward points.**
* **On α−β.** They shrink its uncertainty by a factor ≈ 4–5. The correlated treatment versus exact constraints changes α−β by < 0.02 and γπ by −0.4.
* **On the residual functions** (`results/residual_functions_MAIN.txt`). f^R_2 and f^R_3 are now determined by the forward points with finite errors:

  | E [MeV] | f^R_2 | f^R_3 |
  |---|---|---|
  | 59 | 122 ± 93 | 48 ± 55 |
  | 89 | 55 ± 20 | 63 ± 13 |
  | 150 | 179 ± 12 | 80 ± 8 |

  BχPT gives f^R_2 = −3 → 241 and f^R_3 = 36–40 over the same range. At low energy the (E/M)² suppression makes them weakly determined, without consequences for the multipoles.
* **Unconstrained residuals.** f^R_4, f^R_5 and f^R_6 remain prior-dominated (0 ± 10, 0 ± 1, 0 ± 1); f^R_1 is weakly constrained (±5–10).
* **Figure.** `figs/fig_forward.pdf` shows dσ/dΩ_lab(0°) and Σ2z(0°) from GHP, GHP shifted by the fitted nuisances, and the PWA (MAIN, pruned) at the bin centres. The fit sits on the shifted GHP curve. The shift is −0.5 % to −2.5 % in dσ(0°) (largest near threshold) and ≤ 0.005 in Σ2z(0°), well inside the GHP band.

### 6.6 Mornacchi per-beam-time data (variant)

With separate scales for March/July and dσ/Σ3: α−β = 9.88(1.02), β = 2.15(51), γπ = 6.28(1.58). This is identical to the combined set within 0.03; the scales are within 2 % of unity.

### 6.7 Dynamical polarisabilities (MAIN)

`figs/fig_dynpols.pdf`; values in `results/dynamical_pols.json`. Definitions follow HMP16 Eq. (3.46), with ω the c.m. photon energy.

* **α_E1(ω):** follows the BχPT shape.
* **β_M1(ω):** MAIN lies 1–1.5 below BχPT at 60–110 MeV, while pruned and KLP18 Fit 1″ are closer to BχPT. All fits converge at 130–150 MeV (Δ rise).
* **γE1E1(ω):** compatible with BχPT throughout.
* **γM1M1(ω):** about 1.5 below BχPT at low energy (regularisation bias, §4.3).
* **γE1M2(ω) and γM1E2(ω):** their errors reflect the priors on f^R_5 and f^R_6. γM1E2(ω) is about 1.5 above BχPT.

DR curves are not available numerically.

### 6.8 Variant kept from the previous version: sum rules as exact constraints (`update_fits.py`)

In this treatment (the KLP18 one, previously the main fit), errors include the normalisations and residual priors, and α+β = 14.0 and γ0 = −0.929 are fixed. The results:

#### 6.8.1 Fits with exact constraints

| fit | N | χ²/N | α_E1 | β_M1 | α−β | γπ | γE1E1 | γM1M1 | γM1E2 |
|---|---|---|---|---|---|---|---|---|---|
| pre-2015 data, TAPS stat only (≈ KLP18 Fit 1) | 142 | 1.29 | 11.89(52) | 2.11(52) | 9.78(1.04) | 5.4(2.0) | −2.3(1.0) | 1.35(45) | 1.8(0.9) |
| pre-2015 data, pruned (≈ KLP18 Fit 1″) | 131 | 0.66 | 10.49(65) | 3.51(65) | 6.98(1.29) | 3.6(2.5) | −1.3(1.2) | 0.76(52) | 1.5(1.0) |
| TAPS 2001 only | 55 | 1.21 | 11.41(82) | 2.59(82) | 8.81(1.64) | 3.6(2.9) | −1.4(1.5) | 1.2(0.6) | 1.1(1.3) |
| **A2 2022 + HIγS 2022 only** | 112 | 0.87 | 11.56(45) | **2.44(45)** | **9.11(90)** | 6.5(1.6) | −2.8(0.8) | 0.60(35) | 3.1(0.6) |
| A2 2022 only (dσ+Σ3) | 96 | 0.81 | 11.22(48) | 2.78(48) | 8.44(96) | 5.9(1.6) | −2.5(0.8) | 0.52(36) | 2.9(0.6) |
| all data | 254 | 0.96 | 11.92(36) | 2.08(36) | 9.84(72) | 6.4(1.3) | −2.7(0.6) | 0.71(29) | 3.0(0.5) |
| **all data, floating norms [previous main]** | 254 | 0.92 | **11.96(50)** | **2.04(50)** | **9.91(99)** | **6.8(1.6)** | −2.95(80) | 0.87(36) | 3.00(65) |
| previous main, 3σ outliers removed (Oxley 60/150°, TAPS 108/133°) | 252 | 0.81 | 11.35(52) | 2.65(52) | 8.70(1.04) | 5.8(1.7) | −2.4(0.9) | 0.81(37) | 2.5(0.7) |
| previous main, McGovern pruning (no Oxley, Bernardini, Baranov 150°, TAPS 108/133°) | 243 | 0.74 | 10.89(54) | 3.11(54) | 7.79(1.08) | 5.4(1.8) | −2.2(0.9) | 0.79(38) | 2.4(0.7) |
| previous main, "modern" sets only (TAPS, MacGibbon, Hallin, Federspiel, A2, HIγS) | 206 | – | 11.45(57) | 2.55(57) | 8.89(1.14) | 6.3(1.8) | | | |
| all data (no floating norms) without TAPS | 199 | 0.88 | 12.01(39) | 1.99(39) | 10.02(78) | 6.8(1.4) | | | |
| all data (no floating norms) without A2 | 158 | 0.99 | 12.61(54) | 1.39(54) | 11.22(1.09) | 7.3(2.0) | | | |
| A2 per beam time (4 separate scales) | 350 | 0.88 | 11.96(47) | 2.04(47) | 9.93(94) | 6.7(1.5) | | | |
| γE1M2 free | 254 | 0.96 | 11.88(46) | 2.12(46) | 9.76(93) | 6.0(3.0) | −3.1(2.2) | 0.55(1.17) | 2.9(0.6), γE1M2 = 0.5(2.1) |

#### 6.8.2 Uncertainty budget of the previous main (exact constraints)

`make_tables.py`; total in `results/main_total_errors.json`.

| source | Δ(α−β) | Δγπ | Δβ_M1 |
|---|---|---|---|
| statistical (incl. normalisations and residual priors) | 0.99 | 1.60 | 0.50 |
| residual priors ×3 wider | +0.42 | +0.12 | −0.21 |
| residual priors centred on BχPT | +0.80 | +2.61 | −0.40 |
| α+β = 14.0 ± 0.2 | ∓0.12 | ±0.20 | ±0.16 |
| γ0 = −0.929 ± 0.10 | ±0.09 | ±0.01 | ∓0.04 |
| forward amplitudes as χ² points (correlated nuisances) instead of constraints | −0.01 | −0.40 | +0.09 |
| l = 2 BχPT vs 0 | 0.00 | +0.18 | 0.00 |
| γE1M2 free | −0.08 | −0.39 | +0.04 |
| binning (variant B) | −0.06 | +0.11 | +0.03 |
| A2 per beam time vs combined | +0.02 | −0.17 | −0.01 |
| **total** | **1.35** | **3.13** | **0.70** |
| *database: 3σ outliers removed* | *−1.21* | *−1.06* | *+0.61* |
| *database: McGovern pruning* | *−2.12* | *−1.39* | *+1.06* |

**Result (previous main, exact constraints):**
* α_E1 − β_M1 = 9.9 ± 1.0 (stat) ± 0.9 (method), α_E1 = 12.0 ± 0.7, β_M1 = 2.0 ± 0.5 ± 0.5;
* γπ = 6.8 ± 1.6 ± 2.7;
* correlation ρ(α−β, γπ) = 0.74.

**Database range:** α−β = 7.8 … 9.9, β_M1 = 2.0 … 3.1. The 2018 database gave 7.0 … 10.6 and 1.7 … 3.5 in the same approach.

#### 6.8.3 Outliers and influence (exact constraints, floating norms)

* **3σ rule.** Two points are flagged: Oxley {60 MeV, 150°} (pull −3.9) and TAPS {108.4 MeV, 133°} (−3.7). After removing them no further point exceeds 3σ (the iterated rule converges at once). TAPS {89, 155°} and Goldansky {55, 150°}, outliers in 2018, are now below 3σ.
* **Single removals:**

  | removed | α−β | β |
  |---|---|---|
  | Oxley (4 points) | 8.34 | 2.83 |
  | Oxley 60/150° alone | 8.95 | |
  | TAPS 108/133° | 9.68 | |
  | Baranov 150° | 9.73 | |
  | Bernardini | 9.89 | |
  | Goldansky | 10.62 | 1.69 |
  | Oxley + Goldansky | 9.06 | |
  | all pre-1961 legacy sets | 9.01 | |
  | all TAPS points at θ > 130° | 9.58 | |

  (`results/influence.json`)
* **The leverage has moved.** In 2018 the decisive point was TAPS 108 MeV/133°, which motivated KLP18's call for 109 MeV backward-angle data. The A2 2022 data now cover that region, and the point has little influence left. **The remaining database sensitivity sits at 55–60 MeV, in the 1958–60 Oxley and Goldansky data.**

#### 6.8.4 Consistency between data sets (exact constraints, no floating norms)

* **χ² per set in the global fit** (§5 table): all new sets have χ²/N ≈ 0.7–1.0. TAPS has 1.23 with 5 % point-to-point errors (2.0 with statistical errors only). Oxley (5.6/point) and Baranov 150° (2.2/point) are the inconsistent old sets.
* **Floating normalisations** all lie within 1–2σ of 1: A2 dσ 0.991(13), A2 Σ3 1.014(20), HIγS 0.979(19)/0.996(19), TAPS 1.009(13), Hallin 1.037(27).
* **Separate fits** (α−β, γπ):

  | subset | α−β | γπ |
  |---|---|---|
  | TAPS only | 8.8(1.6) | 3.6(2.9) |
  | A2 only | 8.4(1.0) | 5.9(1.6) |
  | A2 + HIγS | 9.1(0.9) | 6.5(1.6) |
  | all except A2 | 11.2(1.1) | 7.3(2.0) |
  | all except TAPS | 10.0(0.8) | 6.8(1.4) |

  HIγS alone (16 points at one energy) cannot fix four global parameters; its consistency is shown by its χ² in the global fit (15.6/16 with floated scales).
* **Parameter-goodness-of-fit Δχ²** = χ²(all) − χ²(subset) − χ²(rest), data χ² only:

  | subset vs rest | Δχ² |
  |---|---|
  | TAPS | 1.3 |
  | HIγS | 4.1 |
  | new (A2 + HIγS) vs old | 7.0 |
  | A2 | 8.9 |

  At least 3 global parameters (plus shared residuals) are common, so the A2-vs-rest value corresponds to p ≳ 0.03. It is driven by α−β: A2 prefers about 8.4, the old non-TAPS sets about 11.2.
* **Interpretation.** A2 2022, HIγS 2022 and TAPS 2001 are mutually consistent within the PWA. The tension is between the new data and the pre-1975 data, especially Oxley.

#### 6.8.5 Mornacchi per-beam-time data (exact constraints)

With separate scales for March/July and dσ/Σ3: α−β = 9.93(94) and γπ = 6.65(1.51), identical to the combined set (Δ = +0.02 in α−β). The fitted scales are within 2 % of unity.

#### 6.8.6 First study of the forward points as χ² data (`all_std_*fwdchi2*`)

**Implementation.**
* The values are the smooth GHP Re f, Re g evaluated at the bin centres, with f^R_2 and f^R_3 now free per bin (wide prior).
* Correlated treatment: shift = n_f σ_f(ν) + n_g σ_g(ν), with penalty n_f² + n_g² and uncorrelated rounding error 0.0029 μb·GeV.
* α+β = 14.0 + 0.2 n_f and γ0 = −0.929 + 0.105 n_g move coherently with the forward amplitudes, since they come from the same photoabsorption fits.
* Variant: uncorrelated table errors.

**Sign conventions verified (§2).**
* g = (√s/2M)(Φ1 − Φ5)(0°) → −ακ²ν/2M² at low energy.
* The KLP18 multipole sums for f and g hold identically in the Fortran conventions.
* The GHP convention T = f ε′*·ε + g iσ·(ε′*×ε) gives total-helicity-3/2 ↔ f−g, i.e. Σ2z(0°) = −2Re(fg*)/(|f|²+|g|²).

**Results:**

| treatment | α−β | γπ | β | α+β | γ0 | n_f, n_g | χ²/N |
|---|---|---|---|---|---|---|---|
| exact constraints (KLP18) | 9.84(72) | 6.41(1.28) | 2.08(36) | 14.00 | −0.929 | – | 0.957 |
| χ² points, correlated nuisances | 9.68(79) | 6.00(1.35) | 2.25(40) | 14.18 | −0.860 | 0.92, 0.66 | 0.941 |
| χ² points, uncorrelated | 9.78(77) | 6.17(1.37) | 2.12(41) | 14.01 | −0.925 | – | 0.937 |
| no forward points, wide f^R_2,3 | 6.7(4.5) | 3.6(2.3) | 3.6(2.3) | 14.00 | −0.930 | – | 0.868 |
| no forward points, KLP18-size priors on f^R_2,3 | 9.4(6) | 14.2(1.1) | 2.6(3) | 14.60 | −0.885 | 3.0 | 1.127 |

With floating normalisations the correlated treatment gives 9.90(1.08) and 6.43(1.67); with pruning, 8.03(1.15) and 5.31(1.82).

**Impact of the forward points.**
* They are the main constraint on f^R_2 and f^R_3. The constrained and χ² treatments agree, so the residual functions in the two treatments agree as well (`results/residual_functions.txt`; f^R_2 at 59 MeV differs, 34 vs 122, because the low-energy forward points are only weakly constraining there).
* They shrink the α−β uncertainty by a factor ≈ 5. Their correlated uncertainty shifts α−β by < 0.2 and γπ by ≈ 0.4.
* The data mildly prefer the upper side of the GHP band (n_f ≈ +0.9, i.e. α+β ≈ 14.2).


---

## 7. Physics conclusions (MAIN configuration)

1. **The KLP18 picture is confirmed and quantified.** With the 2018-type database the PWA gives β_M1 = 1.4 (all pre-2015 data) to 3.5 (pruned). The size of this database effect is robust against the regularisation of the residual functions and against the treatment of the forward amplitudes.
2. **The new A2 2022 and HIγS 2022 data reduce, but do not remove, the database dependence.**
   * The MAIN fit gives α_E1 − β_M1 = 9.9 ± 1.45 and β_M1 = 2.14 ± 0.73 (total).
   * Reasonable database refinements (3σ rule, pruning, modern sets, legacy-free) move these to α−β = 8.0–9.3 and β_M1 = 2.5–3.1.
   * The new data alone give β_M1 = 2.6 ± 0.7 (stat) and α−β = 8.9 ± 1.4.
   * Overall β_M1 ≈ 2.1–3.1, intermediate between the DR (PPS 2019: 1.8) and χPT (MPG 2013: 3.15; BχPT prediction 3.9) extractions. This is consistent with PDG 2026 (2.31 ± 0.29) and with A2's own fits (HDPV 2.79, BχPT 3.28).
   * **The DR-vs-χPT tension in β_M1 is therefore not resolved by the PWA alone.** The data sit in between, and the PWA uncertainty (±0.7, plus up to +0.9 from database choices) covers both ends.
3. **Data-set consistency.** A2 2022, HIγS 2022 and TAPS 2001 are mutually consistent: χ²/N ≈ 0.7–1.2, normalisations within 1–2 %. There is a ≈ 2.5σ tension (Δχ² = 11) between A2 2022 and the remaining pre-2015 non-TAPS data in α−β (7.3 vs 11.7). Including HIγS raises α−β by 0.5.
4. **Outliers.** The 3σ rule flags only Oxley {60, 150°} and TAPS {108, 133°}.
   * The 2018 key point TAPS 108 MeV/133° is no longer influential (−0.27).
   * The remaining sensitivity comes from the 1958–60 data at 55–60 MeV (Oxley −1.27; Goldansky +0.66).
   * **Recommendation:** new precision dσ/dΩ (and Σ3) at 50–80 MeV, especially backward. Here α−β separates cleanly from the spin polarisabilities. The KLP18 recommendation of ≈ 109 MeV backward data has meanwhile been met by A2 2022.
5. **Forward amplitudes.** Treating them as χ² data with correlated nuisances is equivalent to the KLP18 exact constraints for the scalar polarisabilities (Δ < 0.1) and gives posteriors α+β = 14.18 ± 0.17 and γ0 = −0.87 ± 0.10. They remain indispensable: without them σ(α−β) ≈ 4.5.
6. **Spin polarisabilities.**
   * γπ (non-pole) = 6.4 ± 1.7 (stat) ± 2.7 (method). This agrees with DR (≈ 8), BχPT (7.2 ± 1.7) and the bootstrap (9.9 ± 2.0).
   * γM1M1 and γπ are controlled by the assumed size of the unconstrained residual energy dependence. With KLP18's priors γM1M1 ≈ 0.9 comes out low; with BχPT-centred priors it is 2.6 ± 0.4, at identical χ².
   * Static spin polarisabilities cannot be extracted reliably from sub-threshold unpolarised and Σ3 data alone.

---

## 8. Figures (`work/pwa/figs/`, shared style `work/figstyle.mplstyle`; "MAIN" = §6 configuration)

| file | size | content |
|---|---|---|
| `fig_klp18_reproduction.pdf` | 3.4″ | 59 and 109 MeV: KLP18 data, published Fit 1 / Fit 1″ curves (thick, digitised from the PLB PDF) vs this reproduction (dashed); BχPT from the PLB figure and from the Fortran code |
| `fig_pulls.pdf` | 3.4″ | χ²_i per data point vs E_γ, with 3σ line and outliers labelled. Top: KLP18 database, Fit 1 reproduced (analogue of PLB Fig. 2). Bottom: updated database, MAIN fit |
| `fig_cs_A2_HIGS.pdf` | 7.0″ | dσ/dΩ_lab at the HIγS (81.3; 83.4 ∥/⊥) and A2 (92.2–135.3 MeV) energies. Data: A2 2022, HIγS 2022, TAPS 2001 within 5 MeV. Bands: MAIN and pruned fits (1σ, full covariance). Curves: KLP18 Fit 1, BχPT NNLO prediction, Born+π⁰ |
| `fig_S3_A2.pdf` | 7.0″ | Σ3 at 92.2, 108.5, 129.5 MeV with the A2 2022 data, MAIN and pruned bands, KLP18 Fit 1 and 1″ (predictions), BχPT and Born |
| `fig_polplane.pdf` | 3.4″ | (α−β, γπ) plane: 1σ ellipses (full errors) for KLP18 Fit 1/1″ (reproduced) and this work (MAIN, pruned, A2+HIγS only); published KLP18 points, BχPT (LMP15), DR (OdL01/BGLMN98), DR bootstrap (MRPP22), PDG 2026 α−β band |
| `fig_dynpols.pdf` | 7.0″ | dynamical polarisabilities per bin (MAIN, pruned, KLP18 1/1″ reproduced) vs BχPT NNLO, with static values at E = 0 (BχPT stars) |
| `fig_forward.pdf` | 7.0″ | forward observables dσ/dΩ_lab(0°) = \|f\|²+\|g\|² and Σ2z(0°) = −2fg/(f²+g²) vs E_γ. Shown: GHP table points and smooth representation, GHP shifted by the fitted nuisances (n_f = 0.92, n_g = 0.60), and PWA MAIN/pruned at the bin centres with errors. The right panels show the differences from GHP |

---

## 9. LaTeX tables

The full code is in `work/pwa/results/tables.tex`:
* R1: reproduction;
* D: database, with χ² from MAIN;
* U-old: exact-constraint variant;
* M: MAIN and database variants.

Excerpt of Table M (errors "full"):

```latex
\begin{tabular}{l r c cccc ccc cc}
\hline
Fit & $N$ & $\chi^2/N$ & $\alpha_{E1}$ & $\beta_{M1}$ & $\alpha_{E1}-\beta_{M1}$ & $\gamma_\pi$ & $\gamma_{E1E1}$ & $\gamma_{M1M1}$ & $\gamma_{M1E2}$ & $\alpha_{E1}+\beta_{M1}$ & $\gamma_0$\\\hline
MAIN: all data & 254 & 0.91 & $12.0\pm0.6$ & $2.1\pm0.5$ & $9.9\pm1.1$ & $6.4\pm1.7$ & $-2.8\pm0.8$ & $0.9\pm0.4$ & $2.7\pm0.7$ & $14.2\pm0.2$ & $-0.87\pm0.10$\\
$3\sigma$ outliers removed & 252 & 0.80 & $11.5\pm0.6$ & $2.7\pm0.6$ & $8.8\pm1.1$ & $5.5\pm1.8$ & $-2.3\pm0.9$ & $0.8\pm0.4$ & $2.4\pm0.7$ & $14.2\pm0.2$ & $-0.90\pm0.10$\\
McGovern pruning & 243 & 0.73 & $11.1\pm0.6$ & $3.1\pm0.6$ & $8.0\pm1.1$ & $5.3\pm1.8$ & $-2.2\pm0.9$ & $0.8\pm0.4$ & $2.3\pm0.7$ & $14.2\pm0.2$ & $-0.93\pm0.10$\\
modern sets only & 206 & 0.86 & $11.7\pm0.6$ & $2.5\pm0.6$ & $9.2\pm1.2$ & $6.1\pm1.8$ & $-2.6\pm0.9$ & $0.9\pm0.4$ & $2.6\pm0.7$ & $14.2\pm0.2$ & $-0.91\pm0.10$\\
A2 2022 + HI$\gamma$S 2022 only & 112 & 0.79 & $11.6\pm0.7$ & $2.6\pm0.7$ & $8.9\pm1.4$ & $7.2\pm2.2$ & $-3.2\pm1.1$ & $1.1\pm0.5$ & $2.9\pm0.9$ & $14.2\pm0.2$ & $-0.90\pm0.10$\\
without HI$\gamma$S & 238 & 0.91 & $11.8\pm0.6$ & $2.4\pm0.6$ & $9.4\pm1.1$ & $6.0\pm1.7$ & $-2.6\pm0.9$ & $0.9\pm0.4$ & $2.5\pm0.7$ & $14.2\pm0.2$ & $-0.87\pm0.10$\\
without A2 & 158 & 0.94 & $12.9\pm0.7$ & $1.1\pm0.7$ & $11.7\pm1.4$ & $7.2\pm2.2$ & $-3.1\pm1.1$ & $1.2\pm0.5$ & $2.9\pm1.0$ & $14.0\pm0.2$ & $-0.92\pm0.10$\\
\hline
\end{tabular}
```

KLP18 reproduction excerpt (stage-2 errors, comparable to the PLB):

```latex
\begin{tabular}{l|ccccccc|c}
\hline
Fit & $\alpha_{E1}$ & $\beta_{M1}$ & $\gamma_{E1E1}$ & $\gamma_{M1M1}$ & $\gamma_{M1E2}$ & $\alpha_{E1}-\beta_{M1}$ & $\gamma_\pi$ & $\chi^2/N$\\\hline
Fit 1 (this repr.) & $12.3\pm0.4$ & $1.7\pm0.4$ & $-3.0\pm0.6$ & $1.4\pm0.3$ & $2.5\pm0.7$ & $10.6\pm0.8$ & $6.9\pm1.3$ & $1.36$\\
Fit 1 (KLP18) & $12.2\pm0.3$ & $1.8\pm0.3$ & $-3.1\pm0.7$ & $1.6\pm0.3$ & $2.5\pm0.7$ & $10.5\pm0.4$ & $7.2\pm1.0$ & $1.35$\\\hline
Fit 1$''$ (this repr.) & $10.9\pm0.5$ & $3.1\pm0.5$ & $-2.0\pm0.8$ & $0.9\pm0.3$ & $2.1\pm0.8$ & $7.7\pm1.1$ & $5.0\pm1.6$ & $0.69$\\
Fit 1$''$ (KLP18) & $10.2\pm0.4$ & $3.8\pm0.4$ & $-1.2\pm0.8$ & $0.6\pm0.3$ & $1.6\pm0.8$ & $6.4\pm0.6$ & $3.5\pm1.2$ & $0.62$\\\hline
\end{tabular}
```

---

## 10. Limitations and caveats

* **Legacy data.**
  * Oxley, Goldansky, Bernardini and the full Pugh and Hyman sets exist only in Krupina's notebook. They are secondary sources with total errors as used in 2018.
  * Energies are approximate for some Pugh/Hyman points (bin labels, or the energy of the 90° point of the same group).
  * The Federspiel values are the author's legacy numbers, 3–5 % below the published figure (`db_old_README.md`).
  * These are exactly the sets that drive the remaining spread.
* **The static polarisabilities rest on the regularisation of the residual functions** (§4). The "statistical" errors include the priors; they would grow without them. The prior-centring systematic (BχPT-centred vs zero) is the dominant method uncertainty for γπ and γM1M1.
* **Fixed spin polarisability.** γE1M2 = 0 throughout, as in KLP18; freeing it costs little in α−β but inflates the γ errors.
* **Sum-rule quantities.** In MAIN, α+β and γ0 follow the GHP forward-amplitude nuisances (posteriors 14.18 ± 0.17 and −0.87 ± 0.10). The coherent coupling of the Baldin/γ0 uncertainties to the forward-amplitude uncertainties (one nuisance each) is a modelling assumption; the uncorrelated variant gives the same scalar polarisabilities.
* **Truncation of the multipole expansion.** Non-Born J = 5/2 EM/ME and l = 3 multipoles are neglected. In BχPT this matters at ≤ 0.2 % in dσ and ≤ 0.005 in Σ3 at 135–150 MeV.
* **Near threshold.** Real non-Born multipoles are assumed, including in the 143 and 150 MeV bins just above the π⁰ threshold, as in KLP18.
* **Systematics model of the new data.** A2 Σ3: the uniform 5 % scale is treated as a Gaussian (2.9 %). HIγS: the linear polarisation is taken as 100 %. TAPS: a 5 % point-to-point error is added, following the published paper.
* **KLP18 reproduction.** It used the MAP fit plus a literal MC emulation. The exact published 1′/1″ numbers depend on the random seeds of the original runs, which are not recoverable.
* **No DR comparison curves.** The DR dynamical polarisability curves (Hildebrandt et al.) are not available numerically.
* **Sum-rule inputs.** The GHP photoabsorption parametrisation is not fully published (the low-energy region uses MAID/SAID), so the dispersion integrals were not recomputed. The GHP table was used (smooth interpolation for the update; the exact notebook values for the reproduction). It is consistent with α+β = 14.0, c4 = 6.04, γ0 = −0.929 and γ̄0 = 0.484 within its rounding.

---

## 11. Files

**`work/pwa/`**

* **Machinery**
  * `compton.py`: kinematics, Born+π⁰ (Feynman rules), d-functions, partial waves, multipoles, observables.
  * `fortran_runner.py`: interface to my Fortran copy (`fortran/`; helmode added in `ampl.f90`).
  * `bchpt_multipoles.py` → `results/bchpt_multipoles.csv`.
  * `sumrules.py`: smooth GHP forward amplitudes.
  * `validate_machinery.py` → `validation/validation_log.txt`.
* **KLP18 reproduction**
  * `plb_inputs.py`: KLP18 bins, notebook sum-rule inputs, loader of `db_klp2018.csv`.
  * `pdfvec.py`, `digitize_plb.py`: PDF vector digitisation of PLB Fig. 5 (data points and curves).
  * `build_plb_db.py` → `results/db_plb2018.csv` (independent cross-check of `db_klp2018.csv`).
  * `nb2txt.py`: reader for the Mathematica notebooks.
  * `reproduce_plb.py` → `results/plb_reproduction.json`, `plb_fit1_pulls.csv`, `plb_fit1(pp)_{p,cov}.npy`.
* **Update**
  * `pwa_model.py`: model and fitter.
  * `build_update_db.py` → `results/db_update.csv` and `results/db_update_A2beamtimes.csv`.
  * `main_fits.py` → `results/main_fits.json`: **MAIN** (forward points as χ² data) with budget, database variants, influence and PG tests.
  * `update_fits.py` → `results/update_fits.json`: exact-constraint variants and the first studies.
  * Per-fit pulls, parameters and covariances are in `results/update/` (prefix `M_` for MAIN-configuration fits).
  * `scenarios/no_higs*.py` → `scenarios/no_higs.json`: no-HIγS scenarios (coordinator).
  * Also in `results/`: `influence.json` (exact-constraint config), `residual_functions_MAIN.txt`, `residual_functions.txt` (exact constraints), `chi2_bchpt_prediction.csv`, `dynamical_pols.json` (MAIN).
* **Output**
  * `make_tables.py` → `results/tables.md`, `results/tables.tex`, `results/main_total_errors.json`.
  * `make_figures.py` → `figs/*.pdf`.

**Reproduction.** Run with `.venv/bin/python` from `work/pwa/` in this order:
1. `validate_machinery.py`
2. `bchpt_multipoles.py`
3. `build_update_db.py`
4. `reproduce_plb.py` (about 2 min, MC)
5. `update_fits.py` (seconds)
6. `main_fits.py` (seconds)
7. `make_tables.py`
8. `make_figures.py`

---

## 12. Addendum (2026-10-14, later): configuration used in the paper

The paper's PWA **main fit floats the normalisations of all data sets, including the pre-1961 sets** (10 %, 8 % for
Goldansky), exactly as the BχPT fits do. It is produced by `scenarios/main_fits_legacynorm.py` and stored in
`results/main_fits_legacynorm.json` (key `MAIN`: α−β = 9.34(1.17), β = 2.42(57), α = 11.76(61), α+β = 14.17(17),
χ²/N = 0.834 for N = 254). The fixed-legacy-norm fit presented as MAIN in Sections 5–10 above (9.90/2.14) appears in
the paper only as the last row of its PWA table ("normalizations of the oldest data fixed").

The "method" uncertainty quoted in the paper (±1.0 on α−β, ±0.5 on β, ±0.5 on α) is the quadrature sum of the
variants `w30`, `bchptctr`, `l2zero`, `gEMfree`, `binsB`, `binmode`, `A2beamtimes` and `fwd_uncorr` in
`main_fits_legacynorm.json` (0.98, 0.49, 0.51). Figures `paper/figs/pwa_{pulls,forward}.pdf` come from
`scenarios/make_figures_legacynorm.py`; `paper/figs/fig_multipoles.pdf` from `scenarios/multipoles_figure.py`
(see `paper/figs/collect_figs.sh`). The main fit was re-run on 2026-10-14 (check) and is bit-identical.

**Main database (2026-10-14, later; `work/notes/pwa_mainfit_switch.md`).** By decision of the author both analyses
now use the same main database: all data at E ≤ 150 MeV **minus the single 3σ outlier TAPS 108.42 MeV/133°**
(`update_fits.MAIN_DROP`; BχPT pid `OlmosdeLeon2001_029`), which `update_fits.select` drops from every selection
unless `keep_outlier=True`. The paper's PWA fit (`main_fits_legacynorm.json`, key `MAIN`) is therefore N = 253,
χ²/N = 0.787, α−β = 9.01 ± 1.17 (stat) ± 1.02 (method), β = 2.58 ± 0.58 ± 0.51, α = 11.59 ± 0.61 ± 0.52,
α+β = 14.17(17), γ0 = −0.904(100), γπ = 5.76(1.78); the method uncertainty (same eight variants, key `method`) is
recomputed by the script. All variants are subsets of the main database; `db_allpts` keeps the point (N = 254: the
former MAIN, 9.34/2.42), `MAIN_DROP_check` verifies that the iterated 3σ rule on all points flags exactly this point,
and `db_3sigma` (3σ rule rerun on the main database) flags nothing and equals MAIN. The fixed-legacy-norm fit
(`main_fits.json`, Sections 5–10 above) is on the main database as well (9.58(1.09)/2.33; the numbers printed in
Sections 5–10 are the old all-points ones). Tables: `results/tables.md` now has Tables P-D…P-C for the paper
configuration and D…C for the fixed-legacy-norm one; `results/main_total_errors.json` holds the paper's errors.
Figures: `pwa_pulls.pdf` (lower panel "Updated database, PWA fit", excluded point drawn open with its χ² = 17.2 at
the fitted parameters), `pwa_forward.pdf` (legend "PWA fit"), `fig_multipoles.pdf` (legend "NNLO fit"/"PWA fit")
regenerated with unchanged layout. The χ² comparison with the BχPT Baldin fit on the 197 modern points (158.1 vs
165.3) is in `results/modern_chi2_comparison.md` (`scenarios/modern_chi2_comparison.py`). The KLP18 reproduction
(Sec. 3) is unaffected (rerun, bit-identical). Every OLD → NEW number of the paper is listed in
`work/notes/pwa_mainfit_switch.md`.

## 13. Addendum (2026-10-17): HIγS 2022 dropped from the main database

By decision of the author (2026-10-17) the HIγS 2022 sets (`HIGS_2022_circ`, `HIGS_2022_lin`, 16 points) are no
longer part of the main database of either analysis. `update_fits.select` now drops them from every selection
unless `with_higs=True` (or they are named explicitly in `sets`); the paper's PWA fit (`main_fits_legacynorm.json`,
key `MAIN`) is therefore N = 237, χ²/N = 0.771, α−β = 8.15 ± 1.26 (stat) ± 1.05 (method), β = 2.99 ± 0.62 ± 0.53,
α = 11.14 ± 0.65 ± 0.53, α+β = 14.14(17), γ0 = −0.921(101), γπ = 4.97(1.89), n_f = 0.68, n_g = 0.07 (= the former
variant `db_noHIGS`). The former main fit (N = 253, 9.01/2.58) is kept as the variant **`db_withHIGS`**
(`db_pruned_withHIGS`, `all_withHIGS[_norm]` in `update_fits.json`); `db_new_only` (A2 + HIγS) and `new_only_HIGS`
remain as tests; the PG tests of the HIγS-containing subsets use the with-HIγS fit as "all". MAIN_DROP is re-derived
on the new database (verified; nothing else flagged). Every OLD → NEW number of the paper (Table pwa, Eq. pwamain,
the pulls/forward/multipole figures, the modern-χ² comparison: 181 points, PWA 141.2 vs BχPT 145.1) is listed in
`work/notes/main_without_HIGS.md`; the numbers of Sections 1–12 above refer to the former main database.

## 14. Addendum (2026-10-18): forward points from the 2026 ForwardRCS evaluation; team audit

With FWD_SOURCE=2026 (see `work/notes/forward2026_update.md`) the paper's PWA fit (`main_fits_legacynorm.json`, key `MAIN`) is
N = 237, χ²/N = 0.772, α−β = 8.05 ± 1.22 (stat) ± 1.02 (method), β = 3.08 ± 0.61 ± 0.52, α = 11.13 ± 0.62 ± 0.52,
α+β = 14.21(18), γ0 = −1.008(101), γπ = 5.25(1.89), n_f = 0.95, n_g = 0.16; the modern-χ² comparison is 141.5 (PWA) vs
145.0 (BχPT with Baldin, no truncation covariance) for 181 points. The numbers of Sect. 13 are superseded by these.
`results/influence.json` (Oct 14, GHP2016 source) was stale and has been moved to `work/_to_delete/stale_results_2026-10-18/`;
the live influence table is `main_fits_legacynorm.json['influence']`. Team audit: `work/notes/audit_2026-10-18_team.md`.


## 15. Addendum (2026-10-19): pre-2001 data from Baranov et al. 2001, Table 2
See `work/notes/baranov2001_database.md`. PWA MAIN: α−β = 8.22(1.12)(1.03), β = 2.99(56)(52), α = 11.21(57)(52), χ² = 174.4/237;
legacy normalisations = the review's δ (`update_fits.LEGACY_NORM`); MacGibbon untagged with its full error matrix (whitening in `pwa_model`).
