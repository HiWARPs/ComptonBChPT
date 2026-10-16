# Proton polarisabilities from the updated world Compton database in NNLO BχPT — fit report

*Analysis for the single-author paper (V. Pascalutsa). Date: 2026-10-14.*

- **Baseline: code v2** (`work/code_v2`), which includes the γNΔ vertex corrections and physical pion and nucleon
  masses in the πN loops. The new results are in Part V.
- **Appendix A** keeps the original v1 analysis (`work/code`) unchanged. Its numbers are archived in
  `work/fits/results_v1/` and its figures in `work/fits/figs_v1/`; a verbatim copy of that report is
  `results_v1/fit_report_v1.md`.
- **Reproduction:** `cd work/fits && make` (v2 profile, about 40 s with a filled cache).
- **Tables:** all v2 tables are in `results/tables.md` and `results/tables.tex`. The v1 → v2 comparison is in
  `results/v1_v2_comparison.md`.
- **MAIN-FIT SWITCH (2026-10-14, after the v2 baseline below was written): the main fit is now the one with the
  3σ outlier removed.** The main database ('d_all' in all keys) = all region-L data minus `MAIN_DROP =
  ('OlmosdeLeon2001_029',)` (TAPS, 108.42 MeV, 133°, pull −4.27); the all-points database is the variant 'd_allpts'.
  The new main result is α = 11.59(47)(29), β = 2.94(43)(27), χ² = 219.7/263 (Sect. V.9); the numbers quoted in
  Sects. V.1–V.6 are those of the all-points fit (now 'd_allpts', formerly 'd_all') and are superseded where they
  differ — the complete OLD → NEW list is in `work/notes/mainfit_switch.md`.

Units: α, β in 10⁻⁴ fm³; γ in 10⁻⁴ fm⁴.

---

# Part V — v2 results (baseline)

## V.0 What changed with respect to v1

### Code (main configuration)

- **Settings.** `dvertex=t` (imaginary parts of the γNΔ vertex loops, LMP2015 Fig. 4), `pimassphys=t`,
  `nmassphys=t`, `consts='LMP15'`. The physical masses mean π⁰p loops with m_π⁰ and π⁺n loops with m_π± and M_n,
  so the cusps sit at the physical thresholds, 144.7 and 151.4 MeV.
- **Bug fixes.**
  - The πΔ Thomson-limit leak is fixed.
  - θ = 0° and 180° are now exact: the loop amplitudes are extrapolated in t at fixed s. v1 extrapolated the
    observable instead, which was off by up to 5.6e-4.
  - Analytic reference polarisabilities.
  - LoopTools state leak fixed.
- **New options** used here: the order switch (`order=2|3|3.5`), h_A as a namelist variable, and PDG-2024
  constants.

### NNLO static polarisabilities of the v2 base amplitude

The physical masses add a numerically extracted shift to the πN-loop part:
(α, β, γE1E1, γM1M1, γE1M2, γM1E2) = (11.415, 3.729, −3.246, 2.939, 0.170, 1.133), with γ0 = −0.996 and
γπ = 7.15. The v1 values were 11.246, 3.871, −3.335, 2.943, 0.179, 1.147, with γ0 = −0.933. In all F1/F2 fits the
spin polarisabilities are fixed at these v2 values.

### Pipeline

- `FIT_PROFILE=v2` (default) selects the v2 executable and settings (`lib/bchpt.py`). `FIT_PROFILE=v1` reproduces
  the v1 configuration.
- Database, χ² machinery, nuisance model, database variants and truncation scheme are unchanged. The fixed
  references (`PREF_NNLO`) and the base static polarisabilities are now taken from the driver's `Pref`, which in
  v2 follows the mechanisms switched on.

### Forward points

The main fit now uses the GHP points up to **145 MeV** (11 dσ(0°) + 11 Σ2z(0°), N = 284). Variants: ≤ 135 MeV
(without 145 and 150) and ≤ 150 MeV.

At 145 MeV the v2 NNLO prediction is 2.38 nb/sr against GHP 2.77(16) (v1: 2.17). In the main fit the 145 MeV points
have pulls of +0.05 (dσ) and −0.7 (Σ2z), relative to exp ⊕ truncation errors.

### Truncation errors from the v2 order switch

The scheme-B coefficients were recomputed with `order=2|3` (LO/NLO) and the Born + π⁰ configuration (flags at
order 3.5).

- **Consistency with v1's flag-based orders.** `order=2` ≡ my LO flags and `order=3` ≡ my NLO flags, to all
  printed digits (8 test kinematics). `Pref` follows the order.
- **Size compared with v1.** The ratio σ_th(v2)/σ_th(v1) has median 0.995 over the 284 fit points. Per-set medians
  are 0.986–1.005; individual points range 0.91–1.09, at the cusp. corr(c3^v1, c3^v2) = 0.9997.
- **GP length scales.** Unchanged: ℓ_E = 12–27 MeV, ℓ_θ = 25–57°.

### Validation at θ = 0 (`results/forward_validation.md`)

- **Born + π⁰.** Agrees with the analytic forward amplitudes to 2e-10.
- **NNLO, driver value vs independent small-angle extrapolation.** The comparison uses a cubic in t through the
  8 smallest directly evaluated angles, 7–40°. Differences are ≤ 3e-5 (dσ) and ≤ 1e-4 (Σ2z) for E ≤ 145 MeV, and
  1e-4 / 3e-4 at 150 MeV.
- **NNLO prediction vs GHP, pulls (theory − GHP)/σ.** dσ(0°): −1.2 … −3.8. Σ2z(0°): −0.9 … +0.5. The fit absorbs
  this through α+β (Sect. V.1).

## V.1 Main result (F2 + truncation covariance, N = 284)

| | αE1 | βM1 | α+β | α−β | χ²/dof (p) |
|---|---|---|---|---|---|
| **v2 main** | **11.67 ± 0.47 (fit) ± 0.29 (syst) = 11.7 ± 0.55** | **2.88 ± 0.43 (fit) ± 0.27 (syst) = 2.9 ± 0.5** | 14.55(26) | 8.80(86) | 241.6/264 (0.84) |
| v1 main (Appendix A) | 11.69 ± 0.47 ± 0.31 | 2.88 ± 0.43 ± 0.27 | 14.57(27) | 8.81(86) | 241.2/262 (0.82) |

- **Correlation.** ρ(α, β) = −0.84.
- **Nuisance pulls.** λ_f = −0.61(42), λ_g = −0.01(17). Normalisation pulls are within ±1.4.
- **Spin polarisabilities.** γ0 = −1.00 is the fixed NNLO input (v2 physical masses).
- **F3 (main + γM1M1, γE1E1 free).** γM1M1 = 3.03(49), γE1E1 = −3.25(56), γ0 = −1.08(29) (fixed by the forward
  Σ2z(0°) points), γπ = 7.25(1.01). With γπ free at fixed γ0 (F2π): γπ = 7.32(97).
- **F1 (Baldin 14.0(2) penalty, no forward points) + trunc. cov.** α = 11.41 ± 0.45 ± 0.30, β = 2.76 ± 0.43 ± 0.27,
  χ² = 236.9/245. Without the truncation covariance (LM14-comparable): 10.87(29), 3.40(24).

### Error budget (v2 main)

| source | δα | δβ | v1 δα | v1 δβ |
|---|---|---|---|---|
| statistics (stat ⊕ p2p syst) | 0.12 | 0.12 | 0.12 | 0.12 |
| floating normalisations | 0.30 | 0.28 | 0.30 | 0.28 |
| forward points (λ_f, λ_g) | 0.16 | 0.06 | 0.17 | 0.06 |
| EFT truncation (covariance in χ²) | 0.30 | 0.30 | 0.30 | 0.30 |
| *total fit error* | *0.47* | *0.43* | *0.47* | *0.43* |
| truncation-model choice (Λ_b, δ, c̄, ℓ) | 0.11 | 0.07 | 0.15 | 0.07 |
| database selection (half range d, e, f, g, h, F4) | 0.14 | 0.13 | 0.14 | 0.12 |
| spin polarisabilities (γπ free) | 0.20 | 0.17 | 0.20 | 0.16 |
| gM ± 0.1 | 0.09 | 0.15 | 0.08 | 0.15 |
| **constants → PDG24** (gA 1.2754, fπ 92.07, masses, h_A refitted to Γ = 117 MeV) | 0.02 | 0.01 | 0.01 | 0.01 |
| **h_A ± (Γ_Δ 114–120 MeV ⇒ h_A 2.84–2.91)** | 0.00 | 0.00 | (est.) 0.00 | 0.00 |
| forward-error model (½ \|nuisance − uncorrelated\|) | 0.07 | 0.02 | 0.08 | 0.02 |
| **total** | **0.55** | **0.51** | 0.56 | 0.51 |

Notes on the budget:

- **h_A and the PDG24 constants.** These now run as genuine variants: h_A = 2.84 and 2.91, and `consts=PDG24`.
  Their effect on the fitted α, β is ≤ 0.02. The static polarisabilities change (PDG24: Pref α 11.58, β 3.68;
  h_A = 2.91: α 11.60, β 3.67), but that change is absorbed by the fitted O(p⁴) LECs. Only the residual energy
  dependence matters, and it is negligible below 150 MeV.
- **Conservative next-order estimate.** It is unchanged: δ × |NNLO − NLO⁺| = 0.46 × (0.82, 1.68) = ±0.38 (α) and
  ±0.77 (β), where NLO⁺ is the NLO + δα, δβ fit with the spin polarisabilities matched to NNLO (10.86/4.56).

## V.2 v1 → v2 comparison of the headline numbers

(Full table: `results/v1_v2_comparison.md`.)

| fit | v1 α / β | v2 α / β | Δ(v2−v1) |
|---|---|---|---|
| **main F2 + trunc.** (v1 fwd ≤ 135, v2 fwd ≤ 145 MeV) | 11.69(47) / 2.88(43) | 11.67(47) / 2.88(43) | −0.02 / 0.00 |
| F2 + trunc., fwd ≤ 135 in both | 11.69 / 2.88 | 11.70 / 2.88 | +0.00 / +0.01 |
| F2 + trunc., fwd ≤ 150 | 12.24 / 3.09 | 12.73 / 3.26 | +0.49 / +0.17 |
| F2 + trunc., fwd uncorrelated | 11.53 / 2.83 | 11.53 / 2.84 | 0.00 / +0.01 |
| F2, no trunc. (fwd ≤ 135 in both) | 11.16 / 3.38, χ² 616/262 | 11.21 / 3.43, χ² 478/262 | +0.05 / +0.05 |
| F1 + trunc. | 11.41 / 2.75 | 11.41 / 2.76 | 0.00 / +0.01 |
| F1, no trunc. | 10.89 / 3.36 | 10.87 / 3.40 | −0.03 / +0.04 |
| F3: γM1M1, γE1E1, γ0, γπ | 3.00(48), −3.25(55), −1.07(30), 7.22(99) | 3.03(49), −3.25(56), −1.08(29), 7.25(1.01) | ≤ 0.03 |
| main, (b) A2 only | 10.95 / 3.52 | 10.89 / 3.55 | −0.05 / +0.03 |
| main, (c) HIγS only | 13.94 / 0.46 | 13.94 / 0.42 | 0.00 / −0.04 |
| main, (i) − HIγS | 11.31 / 3.20 | 11.29 / 3.20 | −0.02 / 0.00 |
| main, (e) − TAPS | 11.52 / 3.03 | 11.49 / 3.04 | −0.03 / +0.01 |
| main, (f) GMPF pruning | 11.47 / 3.08 | 11.45 / 3.08 | −0.02 / 0.00 |
| Mornacchi benchmark (thesis 10.66(15), 3.40(20)) | 10.81 / 3.32, χ² 153.9 | 10.85 / 3.33, χ² 153.8 | +0.04 / +0.01 |
| LM14-type (10.6(5), 3.2(5)) | 10.48 / 3.40 | 10.47 / 3.42 | −0.01 / +0.02 |
| PG test HIγS vs rest (main), χ²_PG (p) | 5.31 (0.070) | 5.44 (0.066) | |

**All headline numbers are stable to ≤ 0.05 in α and β.** Neither the vertex corrections (no effect below threshold
by construction) nor the physical masses nor the bug fixes change the region-L results.

- **Effect of the physical masses on the fit.** Comparing with a single pion mass ("v1 physics" inside v2):
  ±0.001 with the truncation covariance; −0.03/+0.04 for F1 without it.
- **Effect on the cross sections.** The physical masses change dσ at 100–140 MeV by only 0.3–0.9%, and the
  truncation covariance down-weights the cusp region anyway.

## V.3 Database variants and the HIγS test (v2 main configuration, + trunc.)

| variant | N | α | β | α−β | χ²/dof | Δα | Δβ |
|---|---|---|---|---|---|---|---|
| (a) pre-2015 only | 172 | 11.74(60) | 2.56(60) | 9.18(1.16) | 145.2/156 | +0.07 | −0.32 |
| (b) A2 2022 only | 118 | 10.89(64) | 3.55(58) | 7.35(1.19) | 78.5/112 | −0.78 | +0.67 |
| (c) HIγS 2022 only | 38 | 13.94(1.17) | 0.42(1.18) | 13.52(2.34) | 21.4/32 | +2.27 | −2.46 |
| **(d) all (main)** | 284 | 11.67(47) | 2.88(43) | 8.80(86) | 241.6/264 | 0 | 0 |
| (e) all − TAPS | 229 | 11.49(49) | 3.04(46) | 8.45(92) | 166.4/210 | −0.18 | +0.16 |
| (f) GMPF/McGovern pruning | 272 | 11.45(47) | 3.08(44) | 8.37(87) | 195.5/255 | −0.23 | +0.20 |
| (g) A2 per beamtime | 380 | 11.73(46) | 2.82(43) | 8.90(85) | 314.6/358 | +0.06 | −0.05 |
| (h) 3σ rule (removes TAPS 108.4 MeV/133° only) | 283 | 11.59(47) | 2.94(43) | 8.65(86) | 219.7/263 | −0.08 | +0.06 |
| **(i) all − HIγS 2022** | 268 | 11.29(50) | 3.20(46) | 8.09(92) | 222.5/250 | **−0.38** | **+0.32** |
| F4: E ≤ 170 MeV | 294 | 11.63(47) | 2.92(43) | – | 252.8/274 | −0.04 | +0.04 |

**HIγS test (v2).** The PG test gives Δχ²_PG = χ²(all) − χ²(HIγS alone, α and β free) − χ²(rest):

| configuration | χ²_PG (dof 2) | p |
|---|---|---|
| main + trunc. | 5.44 | 0.066 |
| main, no trunc. | 6.45 | 0.040 |
| F1 + trunc. | 6.04 | 0.049 |
| F1, no trunc. | 8.38 | 0.015 |

- **Fixed-sum version** (α+β = 14.0 exactly, 1 dof): 7.90 (p = 0.005) without and 5.13 (p = 0.024) with the
  truncation covariance.
- **Prediction of the HIγS points.** The 16 HIγS points predicted by the no-HIγS main fit have χ² = 22.9 (n = 1) or
  20.1 (n fitted: 0.971, 0.987). The largest pulls are −2.6 (CIRC 125°, 11.46 nb/sr) and +2.1 (LIN dσ∥ 90°, 3.99).
- **HIγS alone.** α = 15.18(1.74), β = 2.02(1.97) with α, β free (Li et al.: 15.4(1.8), 2.1(2.0)); 13.79/0.23 with
  fixed sum (Li et al.: 13.8/0.2).
- **Verdict.** Unchanged from v1: a mild tension of about 2σ.

## V.4 E_max scan: does the β drift shrink with the physical cusps? — **No.**

| E_max [MeV] | 90 | 100 | 110 | 120 | 130 | 140 | 150 | 170 |
|---|---|---|---|---|---|---|---|---|
| β, F1 no trunc., v1 | 1.34(69) | 1.79(56) | 2.01(44) | 2.43(36) | 3.11(29) | 3.37(24) | 3.36(24) | 3.45(24) |
| β, F1 no trunc., **v2** | 1.32(69) | 1.78(56) | 2.00(44) | 2.43(36) | 3.13(28) | 3.40(24) | 3.40(24) | 3.49(23) |
| β, main + trunc., v1 | 1.52(71) | 1.94(60) | 2.12(52) | 2.48(49) | 2.73(45) | 2.91(44) | 2.88(43) | 2.94(43) |
| β, main + trunc., **v2** | 1.49(71) | 1.91(60) | 2.10(52) | 2.46(49) | 2.73(45) | 2.91(43) | 2.88(43) | 2.92(43) |

- **Where the drift comes from.** Its rise of about 1.4 units in β between E_max = 110 and 140 MeV happens
  *below* both thresholds. It is driven by the 120–140 MeV data (Hallin, MacGibbon untagged, TAPS 126–143 MeV),
  not by the cusp.
- **Why the physical cusps do not help.** They modify the cross sections there by < 1% (more only above
  144.7 MeV), so the drift is unchanged.
- **Interpretation.** The drift signals the NNLO energy dependence beyond the static polarisabilities
  (higher-order effects on the approach to threshold), not the cusp position. The truncation covariance remains
  essential.

## V.5 Forward points with the physical cusps

| forward treatment (main + trunc.) | α | β | α+β | χ²/dof |
|---|---|---|---|---|
| nuisance, E ≤ 145 MeV (**main**) | 11.67(47) | 2.88(43) | 14.55(26) | 241.6/264 |
| nuisance, E ≤ 135 MeV (without 145 and 150) | 11.70(47) | 2.88(43) | 14.58(27) | 241.0/262 |
| nuisance, E ≤ 150 MeV | 12.73(43) | 3.26(43) | 15.99(6) | 275.0/266 |
| uncorrelated, E ≤ 145 | 11.53(46) | 2.84(43) | 14.37(23) | 237.5/266 |
| uncorrelated, E ≤ 150 | 11.59(45) | 2.86(43) | 14.44(22) | 242.5/268 |
| dσ(0°) only | 11.63(47) | 2.86(43) | 14.49(28) | 237.1/253 |

- **145 MeV point.** It is now usable and changes nothing.
- **150 MeV point.** 1.4 MeV below the physical π⁺n threshold, it is still incompatible with the correlated
  nuisance model. It pulls λ_f to −2.3 and α+β to 16.0, with Δχ² = +33 for 2 points. In v2 the NNLO dσ(0°, 150)
  is 0.91 against GHP 1.12(9); in v1 it was 0.18. That is much closer, but within the fully correlated model the
  steep threshold dependence still gives the point excessive leverage. It is kept as a variant only.
- **Without the truncation covariance**, the forward Σ2z(0°) shape problem is reduced but persists: the
  Σ2z(0°) χ² for the same ≤ 135 MeV points falls from 328 (v1) to 198 (v2), still unacceptable.

## V.6 Δ region with γNΔ vertex corrections

### χ²/N of the parameter-free NNLO prediction (normalisation fitted)

χ² with n = 1 in brackets.

| set | N | v1 | v2 without vertex corr. | **v2** |
|---|---|---|---|---|
| A2 Σ2x (Martel) | 4 | 0.63 (2.5) | 0.64 (2.6) | 0.63 (2.5) |
| A2 Σ2z (Paudyal) | 10 | 1.42 (28.2) | 1.44 (28.5) | 1.48 (**23.6**) |
| LEGS Σ3 (all / < 2π) | 79 / 65 | 1.86 / 2.09 | 1.91 / 2.15 | 1.85 / 2.08 |
| LEGS dσ | 77 | 6.97 (624) | 7.13 (646) | 7.56 (684) |
| MAMI 97 75° / 90° (Hünger) | 12 / 9 | 1.54 / 2.27 | 1.55 / 2.30 | 1.51 / 2.24 |
| MAMI 99 (Wissmann) | 8 | 1.13 | 1.14 | 1.12 |
| MAMI 02 (Camen/Kossert) | 7 | 0.87 (6.8) | 0.87 (6.6) | 0.66 (**4.6**) |
| SAL 93 (Hallin) > 150 MeV | 70 | 3.30 | 3.34 | 3.35 |
| Bonn 72 / Illinois 67 / Moscow 66 | 14 / 5 / 10 | 2.9 / 4.9 / 12.6 | 2.9 / 4.9 / 12.6 | 2.9 / 4.9 / 12.9 |

The vertex corrections improve Paudyal Σ2z (with n = 1: 28.5 → 23.6), Camen dσ (6.6 → 4.6) and LEGS Σ3 slightly
(151 → 146). They worsen LEGS dσ (646 → 684). All other changes are ≤ 2%. This is consistent with CHANGES.md
(Sect. C): the vertex corrections mostly fix Σ2x/Σ2z shapes at 300+ MeV, where there are few data.

### gM fitted to single data sets

Spin polarisabilities at NNLO; α, β constrained to the main fit; prior 2.97(50).

| data | v1 | v2 without vertex corr. | **v2** |
|---|---|---|---|
| A2 Σ2x + Σ2z | 3.228(177) | 3.221(177) | 3.221(164) |
| LEGS Σ3 (< 2π) | 3.406(52) | 3.407(52) | **3.352(47)** |
| MAMI dσ (Hünger, Wissmann, Camen) | 2.912(26) | 2.912(26) | 2.925(24) |
| LEGS dσ | 2.819(19) | 2.819(19) | **2.857(17)** |
| Hallin dσ (> 170 MeV) | 2.920(42) | 2.913(43) | 2.919(41) |
| spread (max − min) | 0.59 (2.82–3.41) | 0.59 | **0.50 (2.86–3.35)** |

**Answer: only marginally.** The vertex corrections move the two extremes towards each other by 0.04–0.05 each:
the spread drops from 0.59 to 0.50, still about 8σ between LEGS Σ3 (3.35(5)) and the MAMI cross sections
(2.925(24)). The incompatibility of
the per-set gM values is therefore not an artefact of the missing vertex corrections.

- **Agreement among the cross sections.** The MAMI and Hallin cross sections agree on gM = 2.92–2.93. In the gM
  scan of the parameter-free χ², the MAMI dσ and A2 Σ2z minima now lie at 2.90–3.00.
- **What disagrees.** LEGS Σ3 (and, weakly, A2 Σ2x) want a larger gM, while LEGS dσ wants a smaller one with a poor
  χ². This mirrors the LEGS–MAMI inconsistency and possibly missing higher-order Δ effects, e.g. running Δ
  parameters (`drunning`, beyond the working order).

### Exploratory spin-polarisability fit (A2 Σ2x+Σ2z + LEGS Σ3 < 2π, γπ = 8.0(1.8), γ0 = −0.929(105))

| | γE1E1 | γM1M1 | γE1M2 | γM1E2 | γπ | χ²/dof |
|---|---|---|---|---|---|---|
| v1 | −2.96(52) | 2.58(41) | −1.39(76) | 2.71(43) | 9.6(1.6) | 84.6/76 |
| v2 without vertex corr. | −2.86(52) | 2.58(41) | −1.47(77) | 2.69(43) | 9.6(1.6) | 84.2/76 |
| **v2 (vertex corr.)** | **−3.10(56)** | **2.69(42)** | **−1.27(77)** | **2.62(44)** | 9.7(1.6) | 85.9/76 |
| v2, + gM floated (prior 2.97(30)) | −3.03(62) | 2.86(77) | −1.32(81) | 2.43(77) | 9.6(1.6) | 85.9/76, gM = 3.04(26) |
| v2, no γπ constraint | −4.19(74) | 4.52(89) | −3.93(1.35) | 4.55(93) | 17.2(3.6) | 80.9/75 |
| v2, A2 Σ2x+Σ2z only (γE1E1, γM1M1) | −3.52(1.16) | 3.13(1.16) | 0.17* | 1.13* | 7.6(2.3) | 18.7/12 |

For comparison, Paudyal 2020 BχPT gives −2.65(43), 2.43(42), −1.32(72), 2.47(42), and MRPP 2022 gives −3.0(6),
3.7(5), −1.2(1.0), 2.0(7).

- **Shifts.** The vertex corrections move the γ's by ≤ 0.24, well within errors. The A2-only fits are unchanged.
- **Without the γπ constraint.** LEGS Σ3 drives γπ even higher with the vertex corrections (14.6 → 17.2).
- **Δ-region dσ.** It is still not used in the spin fits (LEGS–MAMI conflict).

## V.7 What changed qualitatively from v1 to v2

1. **Nothing in region L.** α, β, α±β, the error budget, the database/HIγS/TAPS conclusions, the outlier list and
   the benchmarks all agree to ≤ 0.05. The main result is α = 11.7 ± 0.55, β = 2.9 ± 0.5.
2. **The 145 MeV forward points are now usable**, thanks to the physical cusps. The 150 MeV point is still not.
3. **The β drift with the energy cut is a genuine NNLO-truncation effect, not a cusp artefact** (V.4). This
   strengthens the case for including the truncation covariance.
4. **h_A and the PDG24 constants are now explicit variants.** Their effect on the fitted α, β is ≤ 0.02.
5. **Δ region.** The vertex corrections improve Σ2z and the MAMI 136° cross sections, and slightly narrow the gM
   spread (0.59 → 0.50). The per-set gM tension (LEGS Σ3 vs cross sections) remains. The spin polarisabilities from
   A2 + LEGS Σ3 + γπ change by ≤ 0.24.
6. **The NNLO γ0 changes from −0.93 to −1.00.** This comes from the physical-mass shift of γE1E1 (+0.09). The
   forward-Σ2z fit (F3) gives γ0 = −1.08(29) in both versions.

## V.8 Files (v2)

- **Results (`work/fits/results/`).**
  - `fits_L.json` and `delta_region.json`;
  - `tables.md` and `tables.tex`;
  - `error_budget.json`;
  - `v1_v2_comparison.md`;
  - `forward_validation.md`;
  - `dataset_table.*`.
- **Figures.** `work/fits/figs/*.pdf`, same names and captions as Appendix A, Sect. 6, now from v2. The Δ-region
  figures (fig5a–c) additionally show the NNLO curve without vertex corrections, dotted. In fig7 the 145 MeV
  forward points are "used".
- **Archive.** v1 results in `work/fits/results_v1/` (with `fit_report_v1.md`) and v1 figures in
  `work/fits/figs_v1/`.
- **Code.** `lib/bchpt.py` (profiles `v1`/`v2`), `Makefile` (`FIT_PROFILE`, `FIT_RESULTS`, `FIT_FIGS`; target
  `compare`) and `scripts/s08_compare_v1_v2.py`.

### Remaining caveats

- The vertex corrections are included, but the γNΔ coupling is not refitted in region L. Running Δ parameters and
  isospin-breaking masses in the πΔ loops and in the Δ width are not included.
- All database caveats of Appendix A, Sect. 8 still apply.

## V.9 Main-fit switch (2026-10-14): the main database excludes the 3σ outlier

**Decision (author).** The main fit is "all data, forward points as data, truncation covariance, 3σ outlier removed",
formerly the variant `F2_h_3sigma_thBmain`.

**Implementation (minimal disruption; key names unchanged).**

- `lib/common.py`: `MAIN_DROP = ('OlmosdeLeon2001_029',)`, the pid(s) flagged by the two-pass 3σ rule (|pull| > 3,
  pull = r/√(err_pt² + σ_th²)) around the all-points fit F2 + truncation covariance. Determined once and held fixed;
  `s03_fits_L.py` re-derives it from the all-points fit `F2_d_allpts_thB` at every run and stops if the list differs
  (log: `results/s03_log.txt`, "MAIN_DROP: 3-sigma outliers of the all-points fit"). The point is the TAPS cross
  section at 108.42 MeV, 133° (10.84(94) nb/sr, pull −4.27); the next largest pulls are Federspiel 1991 #6 (+2.70)
  and A2 Σ3 #25 (+2.50).
- `s03_fits_L.py`: the drop is applied in `rows_for()` (switch `keep_outlier`), i.e. to every fit built through
  `do_fit()`: the whole 'd_all' family (d_all, d_all_thB, truncation-model, order/parameter/spin, forward-point, F4,
  E_max scans, F2pi/F3*/F3g0, fixsum, HIGS test, single-observable fits) with unchanged keys, and every subset that
  contains the TAPS set (a_pre2015, g_A2perbeam, i_noHIGS, TAPSonly, pre2015+HIGS, pre2015+A2; the GMPF-pruned
  sets never contained the point, which carries the flag GMPF_out, so (f), F4_pruned and LM14_type are unchanged).
  Subsets without TAPS (b, c, e, A2+HIGS, Mornacchi benchmark) are unchanged. **Choice:** every subset is a subset of
  the main database, so that the database variants measure the effect of the data selection alone.
- New variant **(d') 'd_allpts'** = all points incl. the outlier (`F1/F2_d_allpts[_thB]`); it is in the variants
  tables of `s07` and in the database half-range of the error budget (d, d', e, f, g, h, F4). Its F2 + trunc. row
  equals the former main fit.
- The outlier variants `h_3sigma`, `h_3sigma_thBmain` now run the search on the main database: nothing else is
  flagged (|pull| > 2.5 only: Federspiel #6 +2.70, A2 Σ3 #25 +2.51), so they are copies of the d_all fits with
  `removed = []`; `s07` handles that.
- Downstream: `s05`–`s07`, `delta_plots` (`DP_FITS_JSON`/`DP_DELTA_JSON`), `radcorr/run_radcorr_fit.py` (uses
  `s03.rows_for`), `scenarios/asym_role.py` (uses `s03.do_fit`) follow automatically; `sigma_tot/sigma_tot.py` and
  `pwa/scenarios/multipoles_figure.py` now read α, β, cov from `fits_L.json['F2_d_all_thB']` instead of literals;
  `scenarios/nnlo_prediction_chi2.py` uses the main database (TAPS 54.4/54; 78.6/55 with all points).

**Main result (N = 283 = 261 angular + 22 forward; dof 263).**

| | αE1 | βM1 | α+β | α−β | χ²/dof (p) |
|---|---|---|---|---|---|
| **main (3σ outlier removed)** | **11.59 ± 0.47 (fit) ± 0.29 (syst) = 11.6 ± 0.55** | **2.94 ± 0.43 (fit) ± 0.27 (syst) = 2.9 ± 0.5** | 14.53(26) | 8.65(86)(56) | 219.7/263 (0.98) |
| all points incl. the outlier (d', former main) | 11.67(47)(29) | 2.88(43)(27) | 14.55(26) | 8.80(86) | 241.6/264 (0.84) |

ρ(α, β) = −0.84; λ_f = −0.58(42), λ_g = −0.01(17); normalisation pulls within ±1.43 (Moscow 74 150°); TAPS
χ²/N = 51.9/54 = 0.96 (was 1.33 with the outlier). Error budget: fit components unchanged (0.12/0.12, 0.30/0.28,
0.16/0.06, 0.30/0.30); truncation-model 0.13/0.09 (was 0.11/0.07), database half range 0.11/0.10 (was 0.14/0.13),
the rest unchanged; total 0.55/0.51. Baldin variant: 11.34(45), 2.83(43), 214.8/244; F3: γM1M1 = 3.07(49),
γE1E1 = −3.28(56), γ0 = −1.09(29), γπ = 7.31(1.01). Variants (F2 + trunc.): pre-2015 11.49(60)/2.82(60); no HIγS
11.19(50)/3.28(46) (−0.40/+0.33); A2 per beam time 11.66(46)/2.88(43); E ≤ 170 11.55(47)/2.99(43); no TAPS
+0.09 in β (was +0.16), pruning +0.14 (was +0.20); keeping the outlier +0.08/−0.06. PG test (HIγS) 5.84 (p = 0.054),
fixed sum 5.52 (p = 0.019). E_max scan: β(≤110, + trunc.) = 2.31(52); nested drift 110 → 150 MeV 2.5σ (no trunc.,
was 3.4σ) → 2.2σ (with, was 2.7σ). Δ-region spin fit (α−β = 8.65(86) prior): γ's unchanged within 0.02,
γπ = 9.64(1.63), χ² = 85.6/76. Everything else (A2-only, HIγS-only, benchmarks, g_M per set, σ_tot, radiative
corrections, multipoles) is unchanged or changes in the last digit; see `work/notes/mainfit_switch.md`.

**Figures.** `fig1a_cs_A2_HIGS` (paper Fig. cs): common y-range 4.5–24.5 nb/sr on all panels (drawn range
4.67–20.95, nothing clipped), MIT 57 (Pugh) omitted from the plot only, GHP forward point (spline-interpolated to
the panel energy, quoted error) as an open star at θ = 0. `fig7_forward` (paper Fig. fwd): 2 × 2 panels, absolute
and relative to the empirical (spline through the GHP amplitudes) forward dσ(0) and Σ2z(0), with the main fit ±
truncation (solid), the main fit without truncation covariance (dotted), the Baldin fit (dash-dotted), the NNLO
prediction (dashed) and the empirical values shifted by the fitted n_f, n_g. `fig3_higs_linear`: label "main fit".

---

# Appendix A — v1 analysis (archived; code `work/code`; numbers in `results_v1/`, figures in `figs_v1/`)

*The text below is the v1 report as delivered. Its pipeline description (Sects. 1–3: database, method, fits,
truncation scheme) applies unchanged to v2, except the items listed in V.0. File paths `results/…` and `figs/…`
in this appendix refer to the v1 versions, now in `results_v1/…` and `figs_v1/…`.*

## 0. Executive summary

**Main configuration (user decision).** The main fit (F2) uses all region-L data with E_lab ≤ 150 MeV: the pre-2015
sets, A2 2022 dσ and Σ3, and HIγS 2022. It adds the GHP forward "data" dσ/dΩ_lab(0°) and Σ2z(0°) for 50–135 MeV as
χ² data. Their correlated errors enter through two common nuisance parameters λ_f, λ_g for the Re f and Re g
uncertainties. There is no Baldin penalty, so there is no double counting with the forward points. The spin
polarisabilities are fixed at their NNLO values, every data set has a floating normalisation, and the χ² includes
the EFT truncation covariance ("scheme B", Sect. 3.4).

| | αE1 | βM1 | α+β | α−β | χ²/dof (p) |
|---|---|---|---|---|---|
| **Main fit (F2 + trunc. cov.)** | **11.69 ± 0.47 (fit) ± 0.31 (syst) = 11.7 ± 0.6** | **2.88 ± 0.43 (fit) ± 0.27 (syst) = 2.9 ± 0.5** | 14.57(27) | 8.81(86) | 241.2/262 (0.82) |
| F2 without trunc. cov. | 11.16(25) | 3.38(24) | 14.54(7) | 7.79(48) | 616.4/262 (≈0) |
| F1 (Baldin 14.0(2)) + trunc. cov. | 11.41 ± 0.45 ± 0.28 | 2.75 ± 0.43 ± 0.29 | 14.16(19) | 8.66(86) | 236.1/245 (0.65) |
| F1 without trunc. cov. (Lensky–McGovern-type) | 10.89(29) | 3.36(24) | 14.25(18) | 7.53(50) | 244.9/245 (0.49) |
| F3 = main + γM1M1, γE1E1 free | 11.67(52) | 2.88(46) | 14.55(27) | 8.79(94) | 241.0/260 (0.80) |

Full error budget of the main fit:

- α = 11.69 ± 0.12 (stat) ± 0.30 (normalisations) ± 0.17 (forward-point nuisances) ± 0.30 (EFT truncation) ± 0.31 (model).
- β = 2.88 ± 0.12 ± 0.28 ± 0.06 ± 0.30 ± 0.26.
- The model item combines the truncation-model choice, the database-selection spread, the spin polarisabilities,
  gM ± 0.1, the constants and the forward-error model.

What the main fit implies:

- α+β = 14.57 ± 0.27. This is 1.7σ above the Baldin value 14.0(2). The forward dσ(0°) points alone give
  14.10(9), which does agree with Baldin; the angular data pull the sum up.
- γ0 is not fitted in F2; it is fixed at −0.933, consistent with GHP −0.929(105).
- When the spin polarisabilities are freed (F3), the forward Σ2z(0°) points give γ0 = −1.07(30) and
  γπ = 7.2(1.0); also γM1M1 = 3.00(48) and γE1E1 = −3.25(55). These are consistent with the NNLO values.

Nuisance pulls in the main fit:

- λ_f = −0.65 ± 0.43 and λ_g = −0.05 ± 0.17 (prior N(0,1)), so the forward amplitudes are shifted by well under 1σ.
- Normalisation pulls (n_k − 1)/δ_k are all within ±1.4 (Table 4.2).

### Key findings

1. **Theory (truncation) uncertainty matters.** Without the truncation covariance, the χ²/dof of the main
   configuration is 2.35. The excess comes entirely from the shape of Σ2z(0°) at 120–135 MeV, where the data
   uncertainty is far below the expected size of higher-order corrections. With the covariance, χ²/dof = 0.92.
   - Including it moves the central values by Δα = +0.53 and Δβ = −0.50, with the same sign in F1 and F2.
   - The fitted β drifts strongly with the upper energy cut when no truncation error is included: β = 1.8(6) at
     E ≤ 100 MeV, 2.4(4) at E ≤ 120 MeV and 3.4(2) at E ≤ 150 MeV (Fig. 8). The truncation covariance tempers the
     drift to 1.9(6) → 2.9(4), within errors.
   - The NNLO energy dependence of the 120–150 MeV data (cusp region) is therefore the leading theory issue.
2. **Database.** The full sub-threshold database (262 angular points) is described with χ²/N ≈ 0.9–1.0, with no
   pruning needed.
   - McGovern/GMPF pruning shifts the main fit by Δα = −0.23, Δβ = +0.20. Removing TAPS shifts it by −0.18/+0.15.
   - The 3σ rule removes a single point, TAPS (108.4 MeV, 133°), which is the known GMPF outlier.
   - TAPS has the largest χ²/N, 1.32 for 55 points. Fitted alone it gives β = 1.8(8).
3. **A2 vs HIγS.** In the main configuration:
   - A2 2022 alone gives α = 10.95(65), β = 3.52(58). HIγS alone gives α = 13.94(1.18), β = 0.46(1.18),
     reproducing Li et al.'s HBχPT result 13.8(1.2)/0.2(1.2) with our BχPT. The two differ by Δβ = 3.1 ± 1.3 (2.3σ).
   - Dropping HIγS from the main fit gives α = 11.31(50), β = 3.20(46), a shift of −0.38/+0.32.
   - Parameter-goodness-of-fit test: χ²_PG = 5.3 (dof 2, p = 0.07) with the truncation covariance; 7.2 (p = 0.03)
     without it.
   - The 16 HIγS points predicted by the no-HIγS fit have χ² = 22.7 (n = 1) or 19.9 (n fitted). The largest pulls
     are −2.6 (CIRC 125°) and +2.1 (LIN dσ∥ 90°).
   - The tension is mild, at the 2σ level.
4. **Benchmarks.**
   - Mornacchi's thesis BχPT fit of the per-beamtime A2 data (Baldin 13.8(4), Paudyal spin polarisabilities) is
     reproduced: our α = 10.81(16)_stat, β = 3.32(21)_stat, χ² = 153.9/187 against his 10.66(15), 3.40(20),
     154.6/187.
   - A Lensky–McGovern-2014-type fit (pre-2015 data, GMPF pruning, E ≤ 170 MeV, Baldin 13.8(4), no truncation
     covariance) gives α = 10.48(45), β = 3.40(43), against LM14's 10.6(5), 3.2(5).
5. **Δ region.** The parameter-free NNLO amplitude describes the polarised A2 data and the MAMI cross sections
   reasonably well:
   - χ²/N is 0.6 for Σ2x (Martel), 1.4 for Σ2z (Paudyal, with fitted scale n = 0.90), and 0.9–2.3 for the MAMI
     dσ sets.
   - It describes LEGS dσ (χ²/N ≈ 7), LEGS Σ3 (1.9) and the old sets poorly.
   - Data-driven γNΔ couplings disagree among themselves: MAMI dσ prefer gM = 2.91(3), LEGS Σ3 prefers 3.41(5),
     and A2 Σ2x+Σ2z prefer 3.23(18). Some of this reflects the missing γNΔ vertex corrections (code limitation)
     and the known LEGS–MAMI inconsistency.
   - The exploratory spin-polarisability fit (A2 Σ2x+Σ2z + LEGS Σ3 below the 2π threshold, γπ = 8.0(1.8),
     γ0 = −0.929(105)) gives γE1E1 = −2.96(52), γM1M1 = 2.58(41), γE1M2 = −1.39(76), γM1E2 = 2.71(43). This is
     almost identical to Paudyal et al.'s BχPT fit (−2.65(43), 2.43(42), −1.32(72), 2.47(42)).
   - Without the γπ constraint, LEGS Σ3 drives γπ up to 14.8(3.1).
6. **Size of the theory error.**
   - The Bayesian truncation covariance contributes ±0.30 to α and to β; varying the truncation model adds
     ±0.15 (α) and ±0.07 (β).
   - A cruder "next-order shift" estimate, δ × |fit(NNLO) − fit(NLO+LECs)| with δ = 0.46, gives ±0.4 (α) and
     ±0.8 (β). It should be read as a conservative upper bound.
   - Parameter uncertainties are small: gM ± 0.1 gives ±0.08 (α) and ±0.15 (β); the PDG gA and fπ shift results
     by ≤0.01; hA is negligible.

---

## 1. Pipeline and code interface

- **Location.** All analysis files are in `work/fits/`. Nothing in `work/code` or `work/data` was modified.
- **Running.** `make` (or `./run_all.sh [path/to/ampl.exe]`) runs:
  - `s01` builds the database;
  - `s02` runs the driver for all configurations;
  - `s00` validates the forward (θ = 0) evaluation;
  - `s03` runs the region-L fits;
  - `s04` fits the GP length scales of the truncation coefficients;
  - `s05` runs the Δ-region analysis;
  - `s06` makes the figures;
  - `s07` writes the tables and error budgets.
- **Executable-agnostic.** The executable comes from `$AMPL_EXE` (default `work/code/execs/ampl.exe`); the in.dat
  template from `$AMPL_TEMPLATE` (optional, same `@KEY@` placeholders as `lib/bchpt.py`).
  - Driver outputs are cached in `work/fits/cache/`. The cache key is a hash of the point list, in.dat, baryon.dat
    and the executable path, size and modification time. Rerunning with `code_v2` (`make AMPL_EXE=...`) therefore
    never mixes results.
  - If `code_v2` exposes hA or other new namelist variables, add them to the template and to `BASE`/`CONFIGS` in
    `lib/bchpt.py`.
- **Library** (`work/fits/lib`):
  - `bchpt.py`: runner with configurations NNLO, LO, NLO, BornPi0, NNLO_noPiD, NNLO_noDpole, gM grid, gA, fπ.
    gM enters via `&Delta dng1`. The Δ pole is switched on/off through `baryon.dat`. Order-by-order runs use
    `&incl oneNloop/oneDloop`, `&nucleon nuclamm` and `anomaly`.
  - `database.py`: database assembly.
  - `forward.py`: GHP forward-point error model.
  - `theory.py`: truncation errors.
  - `fitter.py`: χ² machinery.
  - `common.py`: shared helpers.
- **Polarisabilities as fit parameters.** In the driver, P = Pref + d, with Pref the analytic NNLO values
  (11.2464, 3.8707, −3.3351, 2.9427, 0.1787, 1.1469).
  - For incomplete configurations (NLO, no πΔ loops), the total static polarisabilities are P = (static
    polarisabilities of the components switched on, read from the driver header) + d. So order-by-order fits always
    quote total α and β.
  - With modified gM or constants, Pref is the analytic NNLO value for those parameters, so P is again the total
    static polarisability.

---

## 2. Database (`work/fits/db/fitdb.csv`, 624 points; `fitdb_pb.csv` with the A2 per-beamtime tables)

The full table with error models is `results/dataset_table.md`; the LaTeX version is `results/dataset_table.tex`.

### 2.1 Region L (E ≤ 150 MeV; variant E ≤ 170 MeV)

**Sources and choices**

- **True energies.** Taken from db_old/db_new wherever available.
- **Pugh 1957 and Hyman 1959 (full sets).** Taken from db_klp2018, where the energies are bin labels. The true
  energies were reconstructed from the 90° legacy subsets in db_old (Pugh 55, 75, 95, 115, 125, 135; Hyman 60, 70,
  78, 87, 95, 106, 118, 128 MeV). The two Pugh labels without a 90° partner (69 → 65 MeV, 109 → 105 MeV) are
  inferred from the experiment's 10-MeV bin pattern. These points carry the flag `KLP_E`.
- **Oxley (60 MeV), Goldansky (55 MeV), Bernardini (120, 139 MeV).** Taken from db_klp2018 with the stored energies.
  For Bernardini these are bin labels, and KLP converted c.m. → lab at the label energy; the true energies are not
  available (caveat).
- **Duplicates avoided.**
  - The partial Pugh/Hyman legacy subsets in db_old are not used.
  - Hallin Table I is not used (not independent of Table II).
  - Baranov 1966 "refined" is not used (overlaps Baranov 1966 T4).
  - Of LEGS, only the 2001 data are used.
  - The A2 data enter either combined (PRL) or per beamtime (thesis), never both.
- **Frames.**
  - Hallin and all old Δ-region sets: c.m. angle and c.m. dσ (driver `CSCM` with `CM` angles).
  - LEGS 2001: c.m.
  - A2, HIγS, TAPS, MacGibbon, Federspiel, Baranov 1974, Zieger, Wissmann, Camen/Kossert: lab.
  - Martel Σ2x and Paudyal Σ2z: lab angles. This was checked: with c.m. angles, χ² rises from 2.5 to 11.9 (Σ2x)
    and from 28 to 56 (Σ2z).
- **HIγS linear-beam data.** `CS_PAR` and `CS_PERP` map to the driver's `SPAR` and `SPERP` (lab, 100%
  polarisation). Circular and linear beams are two independent normalisation sets.

**Error models (all σ_pt = stat ⊕ uncorrelated point-to-point syst)**

| set | point error | normalisation |
|---|---|---|
| TAPS (Olmos de León) | stat ⊕ 5% (GMPF/MPG/PPS) | 3% |
| Federspiel | stat | 2% |
| Zieger | stat ⊕ 6.4%/4.3% | none (not floated) |
| Baranov 1974/75, 90° and 150° (two sets) | as quoted | 8% |
| Hallin TabII | stat ⊕ (sys ⊖ 3.6%); the sys column contains the 3.6% scale, which is not double counted | 3.6% |
| MacGibbon tagged / untagged | stat | 4% / 2.9% |
| Pugh, Hyman, Oxley, Bernardini | KLP total errors | **10% adopted** (not stated in the sources; assumption) |
| Goldansky | KLP total errors | 8% (MacGibbon 1995) |
| A2 dσ | stat ⊕ p2p | 3% |
| A2 Σ3 | stat ⊕ p2p | 5% uniform scale → Gaussian rms 2.89% |
| HIγS | stat ⊕ p2p | 2.24% per set |

**Counts.** N(≤150) = 150 pre-2015 points (identical to the PPS-2019 "FULL" count), 96 A2 points and 16 HIγS
points, i.e. 262 angular points. Adding the 20 GHP forward points (E ≤ 135 MeV) gives N = 282 in the main fit.

**GMPF/McGovern pruning flag.** The `GMPF_out` flag marks:

- Oxley and Bernardini (whole sets);
- Baranov 150°;
- TAPS (108.42 MeV, 133°);
- Federspiel (44.2 MeV, 135°);
- Hallin above 150 MeV.

### 2.2 GHP forward points (db_forward_obs / db_forward_amp)

- **Definitions.** dσ/dΩ_lab(0) = K(f² + g²) with K = 2.5682 nb/(μb GeV)², and Σ2z(0) = −2fg/(f² + g²).
- **Sign convention.** The convention is σ3/2 − σ1/2. Analytically this coincides with the driver's
  S2Z = (σR − σL)/(σR + σL) for target spin +z: R photon + target +z gives J_z = 3/2. Numerically, the Born forward
  limit has the same sign in both.
- **Error propagation.** The propagated errors from (δf, δg) reproduce the tabulated ones to 0.3%.
- **Correlation model (main).** Each δf_i and δg_i is split into:
  - an uncorrelated rounding part, r = 0.005/√3 μb GeV, since the table is given to two decimals. It produces a
    2×2 (dσ, Σ2z) covariance block at each energy;
  - a fully correlated part, √(δ² − r²), multiplied by the common nuisances λ_f and λ_g, which have N(0,1) priors.
- **Variant.** Uncorrelated diagonal errors as listed.
- **Energy range.** The 145 and 150 MeV points are **not used** in the main fit (E ≤ 135 MeV). The code uses one
  (charged) pion mass in the loops, so its threshold cusp sits at 149.95 MeV, while the physical π⁰p threshold is
  at 144.7 MeV. At 150 MeV the NNLO dσ(0) is 0.18 nb/sr against GHP's 1.12; the t dependence is non-analytic
  there. Including these two points (variant) gives α+β = 15.3, driven by this cusp artefact.
- **Accuracy of the θ = 0 evaluation** (`results/forward_validation.md`). The driver obtains θ = 0 by a cubic
  extrapolation in t (flag 1).
  - **Born + π⁰ check.** At θ = 0 the driver agrees with the analytic forward Born amplitudes
    f = −α/M, g = −ακ²ν/2M² to 2e-10 for 50–150 MeV. The π⁰ pole vanishes forward.
  - **NNLO check.** For the NNLO and strongly shifted polarisabilities, the driver value agrees with an independent
    cubic-in-t extrapolation of direct points at θ_lab = 7–24° to ≤1e-4. That is about 100× below the data
    errors for E ≤ 145 MeV.
  - **No fallback needed.** Computing f and g from small-angle amplitudes was therefore unnecessary.

### 2.3 Region D (150 < E ≤ 350 MeV)

The following sets are used, with their quoted normalisations as floating scales:

- Martel Σ2x: 4 points. The error already contains the 2.1% scale, so it is not floated.
- Paudyal Σ2z: 10 points. The "Rand." error is used, plus the part of "Syst." beyond the 10.4% scale; the 10.4%
  scale is floated.
- LEGS 2001: Σ3 (79 points) and dσ (77 points).
- Hünger 75° c.m. and 90° c.m.
- Wissmann (130.7° lab) and Camen/Kossert (136.2° lab).
- Hallin TabII above 150 MeV.
- TAPS at 157–164 MeV.
- Gray, Genzel, Baranov 1966 (c.m. assumed), Kolb and Nagashima.

Wolf/Galler (LARA) are not available as numbers.

---

## 3. Method

### 3.1 χ²

- **Definition.**
  χ² = (y − m)ᵀ C⁻¹ (y − m) + Σ_k λ_k² + λ_f² + λ_g² + Σ constraints, with
  m_i = n_k(i) t_i(P) + λ_f s_f,i + λ_g s_g,i and n_k = 1 + δ_k λ_k.
- **Normalisations.** The floating normalisation multiplies the theory, so there is no d'Agostini bias. For
  asymmetries it is a scale factor (A2 Σ3, Paudyal Σ2z).
- **Covariance.** C = diag(σ_pt²) ⊕ the forward rounding blocks ⊕ (optionally) the truncation covariance.
- **Degrees of freedom.** dof = N_points + N_constraints − N_free − N_norms − N_nuisances (GMPF/MPG convention).
  p-values are χ² tail probabilities.
- **Minimisation.** Levenberg–Marquardt on the whitened residuals; the parameter covariance comes from a
  central-difference Jacobian.
- **Polarisability shifts.** These enter exactly as ratios of quadratic forms, so there is no linearisation.
- **Error decomposition by sequential release.**
  - "stat": all nuisances fixed at their best-fit values;
  - "+norm": normalisations released;
  - "+forward/Baldin": the remaining constraint released;
  - each item is the quadrature difference between successive steps.
  - In fits with the truncation covariance, the pure experimental "stat" is taken from the same fit without the
    covariance. The truncation item is then √(stat⊕th² − stat²).

### 3.2 Fits

- **F1.** α, β free; Baldin α+β = 14.0(2) (GHP 2015) as a penalty.
- **F2 (main).** α, β free; GHP forward points; no Baldin.
- **F3 and its variants.**
  - F3: F2 + γM1M1 and γE1E1 free (γ0 constrained by Σ2z(0°)).
  - F2π: F2 + γπ free at fixed γ0, along the direction δγM1M1 = −δγE1E1 = ½.
  - F3g0: forward dσ(0) + the GHP γ0 = −0.929(105) penalty (instead of Σ2z(0)), γM1M1 and γE1E1 free.
  - F3b: F1 + γπ free.
  - F3c: F1 + γM1M1 free (MPG style).
- **F4.** E ≤ 170 MeV (TAPS up to 164 MeV, Hallin ≤ 150 MeV as in MPG Set III; variant with Hallin ≤ 170 MeV).
- **Database variants (a)–(i).**
  - (a) pre-2015 only;
  - (b) A2 2022 only;
  - (c) HIγS 2022 only;
  - (d) all;
  - (e) all − TAPS;
  - (f) GMPF/McGovern pruning;
  - (g) A2 per beamtime;
  - (h) 3σ outliers removed (two passes; with the truncation covariance the pull uses √(σ_pt² + σ_th²));
  - (i) all − HIγS (requested; both HIγS sets removed).

### 3.3 Order-by-order observables from the code

| order | content | how it is run |
|---|---|---|
| LO, O(p²) | Born without κ (Thomson/Dirac) | `nuclamm=f, anomaly=f, oneNloop=f, oneDloop=f`, Δ pole off in baryon.dat |
| NLO, O(p³) | + Pauli (κ) Born terms + π⁰ pole + πN loops | `oneDloop=f`, Δ pole off |
| NNLO, O(p^{7/2}) | + Δ pole + πΔ loops | default |

Additional configurations:

- BornPi0 (Born with κ + π⁰ pole): the reference amplitude of scheme B.
- NNLO without πΔ loops: diagnostic, and used for the hA estimate.

The driver evaluates the base amplitude with all contact terms zero and adds the O(p⁴) contact terms (shifts d) on
top. Hence any configuration can be evaluated at prescribed total static polarisabilities.

### 3.4 EFT truncation error (`lib/theory.py`)

**Expansion.**

- y = y_ref (c0 + c2 δ² + c3 δ³), with c2 = (y_NLO − y_LO)/(y_ref δ²) and c3 = (y_NNLO − y_NLO)/(y_ref δ³).
- The truncation error is σ_th = y_ref · c̄ · δ⁴/√(1−δ²). This sums all omitted orders δ⁴, δ⁵, … in quadrature
  (Furnstahl et al. 2015; Melendez et al. 2019).
- c̄² = (c2² + c3²)/2 by default; variant c̄ = max(|c2|, |c3|).
- δ(E) = √(max(mπ, E_lab)/Λ_b) with Λ_b = 650 MeV, i.e. δ = 0.463 below mπ. This is the δ counting:
  δ² ≈ mπ/Λ and δ ≈ mπ/Δ = 0.48.
- Variants: Λ_b = 500 and 800 MeV; fixed δ = 0.48.
- y_ref = |y_LO| for cross sections and 1 for asymmetries.

**Two schemes.**

- **Scheme A (naive).** LO, NLO and NNLO are each evaluated with their own static polarisabilities. At NLO,
  β = −1.8 instead of 3.9, so the order steps are dominated by the static polarisabilities. Those are precisely
  the O(p⁴) LEC directions that the fit determines, so scheme A double counts.
  - It gives σ_th ≈ 5% of dσ and up to 15× the forward-point errors.
  - It is shown only as a variant.
- **Scheme B (LEC-matched, default).** All orders are evaluated with the same six static polarisabilities (the
  NNLO values) via the contact terms. Born (with κ) + π⁰ is treated as exact.
  - The steps then measure only the energy dependence beyond the static polarisabilities (dispersive, quadrupole,
    cusp, Δ tail), which the fit cannot absorb.
  - Typical σ_th/dσ: 0.2% at ≤ 80 MeV, 0.7% at 80–100 MeV, 1.5% at 100–120 MeV, 3% at 120–140 MeV and 4% at
    140–150 MeV. Relative to the point errors that is a median of 0.38 for A2 dσ, 0.17 for TAPS, 0.06 for A2 Σ3
    and ~1.0 / 0.5 for the forward dσ(0) / Σ2z(0).

**Correlations.**

- Squared-exponential kernel in (E_lab, θ_cm) within one observable family (dσ, σ∥, σ⊥, Σ3, Σ2x, Σ2z), no
  correlation across families.
- Default lengths ℓ_E = 25 MeV, ℓ_θ = 40°, from a GP marginal-likelihood fit to the scheme-B coefficients on a grid
  (`s04_theory_lengths.py`), which gives ℓ_E = 12–26 MeV and ℓ_θ = 25–58°.
- The coefficients grow roughly as E³–E⁴ below mπ, so a stationary GP is only approximate. Variants:
  ℓ = 12/25, 50/50 and 100/100 (MeV/deg), and diagonal.

**Effect on α, β.** The truncation covariance is either included in χ² (fits labelled "+ trunc.") or not.

- Its contribution to the error is ±0.30 for both α and β.
- Its central-value shift is +0.53 (α) and −0.50 (β).

---

## 4. Results in region L

### 4.1 Main fits

(The LaTeX version of all tables is in `results/tables.tex`.)

| fit | N | αE1 | βM1 | α+β | α−β | χ²/dof | p | γM1M1 | γE1E1 | γ0 | γπ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **MAIN: F2 + trunc. cov.** | 282 | 11.69(47) | 2.88(43) | 14.57(27) | 8.81(86) | 241.2/262 | 0.82 | 2.94* | −3.34* | −0.93* | 7.25* |
| F2, no trunc. cov. | 282 | 11.16(25) | 3.38(24) | 14.54(7) | 7.79(48) | 616.4/262 | 0.00 | * | * | * | * |
| F2, forward uncorrelated + trunc. | 282 | 11.53(46) | 2.83(43) | 14.36(23) | 8.70(86) | 237.4/264 | 0.88 | * | * | * | * |
| F2, forward uncorrelated, no trunc. | 282 | 11.02(25) | 3.37(24) | 14.39(9) | 7.65(48) | 266.3/264 | 0.45 | * | * | * | * |
| F2π: + γπ free (fixed γ0) + trunc. | 282 | 11.71(51) | 2.86(46) | 14.58(27) | 8.85(93) | 241.2/261 | 0.81 | 2.99(48) | −3.38(48) | −0.93* | 7.33(96) |
| F3: + γM1M1, γE1E1 free + trunc. | 282 | 11.67(52) | 2.88(46) | 14.55(27) | 8.79(94) | 241.0/260 | 0.80 | 3.00(48) | −3.25(55) | −1.07(30) | 7.22(99) |
| F3, no trunc. cov. | 282 | 10.92(27) | 3.53(26) | 14.44(7) | 7.39(53) | 303.8/260 | 0.03 | 3.20(31) | −2.84(31) | −1.69(4) | 7.01(62) |
| F3g0: fwd dσ(0) + GHP γ0 penalty + trunc. | 272 | 11.65(52) | 2.85(46) | 14.50(29) | 8.80(94) | 236.5/251 | 0.74 | 2.97(48) | −3.37(49) | −0.92(10) | 7.30(96) |
| F1 (Baldin 14.0(2)) + trunc. | 262 | 11.41(45) | 2.75(43) | 14.16(19) | 8.66(86) | 236.1/245 | 0.65 | * | * | * | * |
| F1, no trunc. cov. | 262 | 10.89(29) | 3.36(24) | 14.25(18) | 7.53(50) | 244.9/245 | 0.49 | * | * | * | * |
| F3b: F1 + γπ free + trunc. | 262 | 11.39(49) | 2.77(46) | 14.16(19) | 8.63(93) | 236.1/244 | 0.63 | 2.89(48) | −3.28(48) | −0.93* | 7.14(96) |
| F3c: F1 + γM1M1 free + trunc. | 262 | 11.34(46) | 2.80(44) | 14.14(19) | 8.54(88) | 235.6/244 | 0.64 | 2.59(53) | −3.34* | −0.58(53) | 6.89(53) |

(\* = fixed at NNLO.) Correlation coefficients: ρ(α, β) = −0.83 (main) and −0.91 (F1 + trunc.).

**What the spin-polarisability fits show (F3).** The sub-threshold data, with the A2 Σ3 and the backward
cross sections, constrain γπ to ±1.0 (and to ±0.6 without the truncation covariance), and γM1M1 to ±0.5. The
results are consistent with NNLO (γπ = 7.25, γM1M1 = 2.94) and with the MAMI dispersive γπ = 8.0(1.8). Freeing
them inflates the α, β errors only slightly (spin item of the budget). Fixing the spin polarisabilities instead
at the MRPP-2022 values (γπ = 9.9) shifts β by −0.40; at the Paudyal values, by −0.12.

**What the forward points show.**

- **dσ(0°).** Alone, the 10 forward dσ(0) points give α+β = 14.10(9), in agreement with Baldin; the NNLO forward
  energy dependence fits them with χ² = 5.8/7.
- **Σ2z(0°).** Without the truncation covariance, the nuisance model fails, with χ² = 328 for the 10 Σ2z points.
  - The NNLO g(ν)/f(ν) shape differs from GHP by 0.01–0.03 at 120–135 MeV. This is consistent with the known
    difference in the higher forward spin polarisability: NNLO γ̄0 = 1.1(5) against GHP 0.48(8).
  - With the truncation covariance, the forward χ² is 3.0 (dσ) + 4.6 (Σ2z) for 20 points.
- **Free spin polarisabilities without the covariance.** F3 without the truncation covariance pulls γ0 to −1.69
  to accommodate that shape. With the covariance, γ0 = −1.07(30).

### 4.2 Normalisations and per-set χ² (main fit)

| set | N | χ² | χ²/N | n_k | prior δ_k | (n−1)/δ |
|---|---|---|---|---|---|---|
| A2 22 dσ | 60 | 47.6 | 0.79 | 0.988(12) | 3.0% | −0.41 |
| A2 22 Σ3 | 36 | 26.0 | 0.72 | 1.008(17) | 2.9% | +0.26 |
| Baranov 74 150° | 4 | 3.5 | 0.87 | 0.894(28) | 8% | −1.33 |
| Baranov 74 90° | 3 | 2.9 | 0.97 | 1.019(27) | 8% | +0.24 |
| Bernardini 60 | 2 | 1.0 | 0.48 | 0.991(65) | 10% | −0.09 |
| Federspiel 91 | 16 | 17.9 | 1.12 | 0.998(16) | 2% | −0.10 |
| Goldansky 60 | 5 | 3.8 | 0.77 | 1.040(23) | 8% | +0.50 |
| HIγS dσ (circ.) | 8 | 8.4 | 1.05 | 0.974(19) | 2.24% | −1.14 |
| HIγS dσ∥,⊥ | 8 | 8.5 | 1.06 | 0.990(19) | 2.24% | −0.43 |
| Hallin 93 (≤150) | 13 | 8.7 | 0.67 | 1.044(26) | 3.6% | +1.22 |
| Hyman 59 | 12 | 1.6 | 0.13 | 1.001(48) | 10% | +0.01 |
| MacGibbon tagged | 8 | 5.3 | 0.66 | 1.000(30) | 4% | −0.01 |
| MacGibbon untagged | 10 | 3.5 | 0.35 | 0.999(23) | 2.9% | −0.03 |
| TAPS (Olmos de León) | 55 | 72.5 | 1.32 | 1.000(13) | 3% | −0.01 |
| Oxley 58 | 4 | 9.3 | 2.32 | 0.919(21) | 10% | −0.81 |
| Pugh 57 | 16 | 3.9 | 0.24 | 0.961(47) | 10% | −0.39 |
| Zieger 92 | 2 | 2.6 | 1.32 | 1 (fixed) | – | – |
| GHP dσ(0°) | 10 | 3.0 | 0.30 | nuisance λ_f = −0.65(43) | | |
| GHP Σ2z(0°) | 10 | 4.6 | 0.46 | nuisance λ_g = −0.05(17) | | |

χ²_data = 234.6, χ²_norm = 6.1, χ²_nuis = 0.4.

How to read the per-set χ²:

- With a correlated truncation covariance, the per-set χ² are the contributions Σ_i r_i(C⁻¹r)_i, which sum to
  χ²_data.
- The pulls in Fig. 6 are r_i/σ_pt,i.
- The old sets (Hyman, Pugh) have χ²/N ≪ 1, i.e. generous total errors.
- Oxley has χ²/N = 2.3, consistent with GMPF's remark about its tiny errors.

### 4.3 Database variants (main configuration, + trunc. cov.)

| variant | N | αE1 | βM1 | α−β | χ²/dof | p | Δα | Δβ |
|---|---|---|---|---|---|---|---|---|
| (a) pre-2015 only | 170 | 11.73(60) | 2.60(60) | 9.13(1.17) | 145.6/154 | 0.67 | +0.04 | −0.28 |
| (b) A2 2022 only | 116 | 10.95(65) | 3.52(58) | 7.43(1.20) | 78.7/110 | 0.99 | −0.75 | +0.64 |
| (c) HIγS 2022 only | 36 | 13.94(1.18) | 0.46(1.18) | 13.48(2.34) | 21.8/30 | 0.86 | +2.25 | −2.42 |
| **(d) all (main)** | 282 | 11.69(47) | 2.88(43) | 8.81(86) | 241.2/262 | 0.82 | 0 | 0 |
| (e) all − TAPS | 227 | 11.52(50) | 3.03(46) | 8.49(92) | 166.6/208 | 0.98 | −0.18 | +0.15 |
| (f) GMPF/McGovern pruning | 270 | 11.47(47) | 3.08(44) | 8.39(87) | 195.3/253 | 1.00 | −0.23 | +0.20 |
| (g) A2 per beamtime | 378 | 11.74(46) | 2.83(43) | 8.91(85) | 314.2/356 | 0.95 | +0.05 | −0.04 |
| (h) 3σ outliers removed | 281 | 11.61(47) | 2.94(43) | 8.67(86) | 219.5/261 | 0.97 | −0.08 | +0.06 |
| (i) all − HIγS 2022 | 266 | 11.31(50) | 3.20(46) | 8.11(92) | 222.3/248 | 0.88 | −0.38 | +0.32 |
| F4: E ≤ 170 MeV | 292 | 11.63(47) | 2.94(43) | 8.69(86) | 253.2/272 | 0.79 | −0.06 | +0.07 |
| F4 + GMPF pruning | 280 | 11.40(47) | 3.15(44) | 8.25(87) | 207.3/263 | 1.00 | −0.29 | +0.27 |
| F4 with Hallin ≤ 170 | 307 | 11.30(46) | 3.24(43) | 8.06(85) | 312.8/287 | 0.14 | −0.39 | +0.36 |

In the main configuration, α+β stays at 14.5–14.6 for every variant. F1 variants (Baldin), and all variants
without the truncation covariance, are in `results/tables.md`. Without the covariance, the same pattern holds,
e.g. (e) −0.25/+0.25, (f) −0.17/+0.16, (i) −0.14/+0.13 relative to (d) = 11.16/3.38. The χ² values there are
dominated by the forward Σ2z.

**Outliers (3σ rule).** In all configurations exactly one point is removed: TAPS (108.42 MeV, 133°), with pull
−4.2. Points with |pull| > 2.5 are Federspiel (44.2 MeV, 135°) at +2.7, which is the second GMPF outlier, and A2 Σ3 (129.5 MeV, 35°) at +2.5. No
second-pass outliers. (KLP's Fit 1_3σ removed four points with a PWA and fixed normalisations; floating
normalisations absorb the Oxley/Goldansky 150° and TAPS 89/155 tensions here.)

**HIγS consistency test (variant (i), requested).**

| configuration | α, β (all) | α, β (no HIγS) | α−β (no HIγS) | Δα, Δβ | χ²/dof all → no HIγS | χ²_PG (dof 2), p | χ²(HIγS predicted), n = 1 / n fitted |
|---|---|---|---|---|---|---|---|
| F2 + trunc. (main) | 11.69, 2.88 | 11.31(50), 3.20(46) | 8.11 | −0.38, +0.32 | 241.2/262 → 222.3/248 | 5.31, 0.070 | 22.7 / 19.9 |
| F2 | 11.16, 3.38 | 11.03(25), 3.51(24) | 7.52 | −0.14, +0.13 | 616.4/262 → 595.5/248 | 7.23, 0.027 | 24.5 / 21.2 |
| F1 + trunc. | 11.41, 2.75 | 11.03(48), 3.09(46) | 7.94 | −0.38, +0.34 | 236.1/245 → 216.5/231 | 5.94, 0.051 | 23.7 / 20.6 |
| F1 | 10.89, 3.36 | 10.71(30), 3.49(24) | 7.22 | −0.18, +0.13 | 244.9/245 → 223.1/231 | 8.14, 0.017 | 26.0 / 22.2 |
| F3 + trunc. | 11.67, 2.88 | 11.25(55), 3.23(49) | 8.02 | −0.42, +0.35 | 241.0/260 → 221.9/246 | 5.43, 0.066 | 23.0 / 20.1 |
| F3 | 10.92, 3.53 | 10.76(28), 3.68(27) | 7.08 | −0.16, +0.15 | 303.8/260 → 281.8/246 | 8.32, 0.016 | 26.2 / 22.3 |

How the PG test is set up:

- **Definition.** Δχ²_PG = χ²(all) − χ²(HIγS alone, α and β free, no sum constraint) − χ²(rest, with
  Baldin/forward). The Baldin or forward information is attributed to the "rest". dof = 2 + 2 − 2.
- **HIγS alone.** α = 15.18(1.74), β = 2.02(1.97), χ² = 13.6/13, against Li et al. 15.4(1.8), 2.1(2.0).
- **Fixed-sum version** (α+β = 14.0 exactly in all three fits; one shared parameter): χ²_PG = 7.6 (p = 0.006)
  without and 5.0 (p = 0.025) with the truncation covariance. α−β: all 7.32/8.57, no HIγS 7.04/7.86, HIγS 13.53.
- **Truncation-covariance caveat.** In the "+ trunc." PG numbers, the subset fits drop the theory correlations
  between HIγS and the rest (approximation).

Pulls of the 16 HIγS points predicted by the no-HIγS main fit (n_HIγS refitted with the 2.24% prior and P fixed:
CIRC 0.971, LIN 0.987):

| point | E [MeV] | θ_lab | data | prediction | pull (n = 1) | pull (n fitted) |
|---|---|---|---|---|---|---|
| CIRC 55° (×2) | 81.3 | 55 | 11.45(1.34), 11.61(1.82) | 12.08 | −0.47, −0.26 | −0.21, −0.07 |
| CIRC 90° (×3) | 81.3 | 90 | 8.87(1.11), 10.14(1.28), 9.55(0.82) | 10.70 | −1.64, −0.43, −1.39 | −1.36, −0.19, −1.02 |
| CIRC 125° (×3) | 81.3 | 125 | 11.46(0.92), 12.83(1.33), 15.08(1.48) | 13.84 | −2.58, −0.76, +0.84 | −2.15, −0.46, +1.11 |
| dσ∥ 55° | 83.4 | 55 | 5.19(1.16) | 5.69 | −0.43 | −0.37 |
| dσ∥ 90° (×2) | 83.4 | 90 | 3.10(1.01), 3.99(0.88) | 2.14 | +0.95, +2.12 | +0.98, +2.15 |
| dσ∥ 125° (×2) | 83.4 | 125 | 7.06(0.95), 7.99(0.75) | 8.56 | −1.58, −0.76 | −1.47, −0.62 |
| dσ⊥ 55°, 90°, 125° | 83.4 | 55/90/125 | 17.03(2.27), 18.24(1.18), 19.18(1.18) | 18.19, 19.26, 19.14 | −0.51, −0.86, +0.03 | −0.41, −0.66, +0.23 |

### 4.4 Energy-cut dependence (Fig. 8)

| E_max [MeV] | 90 | 100 | 110 | 120 | 130 | 140 | 150 | 170 |
|---|---|---|---|---|---|---|---|---|
| β, main config. + trunc. | 1.52(71) | 1.94(60) | 2.12(52) | 2.48(49) | 2.73(45) | 2.91(44) | 2.88(43) | 2.94(43) |
| β, main config., no trunc. | 1.57(68) | 2.01(55) | 2.14(43) | 2.51(36) | 3.14(28) | 3.38(24) | 3.38(24) | 3.47(24) |
| β, F1, no trunc. | 1.34(69) | 1.79(56) | 2.01(44) | 2.43(36) | 3.11(29) | 3.37(24) | 3.36(24) | 3.45(24) |

The rise of β between 120 and 140 MeV comes mainly from the old data:

- Hallin, MacGibbon untagged and TAPS at 130–145 MeV. Fitted alone, the pre-2015 data give β = 1.4–1.9 for
  E ≤ 110–130 MeV and 2.7 for E ≤ 140 MeV.
- The A2 data alone are stable within errors: β = 2.7(1.2), 3.3(7), 3.6(5), 4.0(4) for E ≤ 100, 110, 120, 130 MeV.

This is the main argument for including the truncation covariance, which grows rapidly towards the cusp.

### 4.5 Benchmarks, order-by-order fits and parameter variations

**Mornacchi benchmark (thesis Table 10.4, BχPT, fixed spin polarisabilities).**

- Setup: A2 per-beamtime data (192 points), four scale factors (dσ 3%, Σ3 5%/√3), Baldin 13.8(4), spin
  polarisabilities = Paudyal (−2.87, 2.70, −0.85, 2.04).
- Result: α = 10.81 ± 0.16 (stat) [± 0.45 total], β = 3.32 ± 0.21 (stat) [± 0.30 total], χ² = 153.9/187.
- Thesis: 10.66(15)(41), 3.40(20)(22), χ² = 154.64/187.
- Scale pattern: theory scales 0.988/0.979 (dσ Mar/Jul) and 0.991/1.007 (Σ3), matching his data scales
  1.008/1.016 and 1.012/0.995 (his S multiplies the data).
- Residual differences: +0.15 in α, −0.08 in β. These are compatible with code-version and constant differences.

**Lensky–McGovern-2014 type.**

- Setup: pre-2015 data, GMPF pruning, E ≤ 170 MeV with Hallin ≤ 150 MeV (148 points), Baldin 13.8(4), no
  truncation covariance.
- Result: α = 10.48(45), β = 3.40(43), χ² = 112.3/138.
- LM14: 10.6 ± 0.25 ± 0.2 ± 0.4, 3.2 ∓ 0.25 ± 0.2 ∓ 0.4.

**Order-by-order (main configuration, + trunc.).**

| amplitude | α | β | χ²/dof |
|---|---|---|---|
| NLO + δα, δβ | 10.44(45) | 5.29(42) | 457/262 |
| NLO with all spin polarisabilities matched to NNLO | 10.84(47) | 4.54(43) | 263/262 |
| NNLO | 11.69(47) | 2.88(43) | 241/262 |

The NLO → NNLO shift of the fitted values (spin-matched) is 0.85 (α) and 1.66 (β). Multiplied by δ = 0.46 this
gives the conservative next-order estimate ±0.39 (α), ±0.77 (β) quoted in the summary.

**Non-static effect of the πΔ loops.** Removing them (static polarisabilities kept fixed) shifts the fit by
−0.11 (α) and +0.13 (β).

**Parameter variations (main configuration, + trunc.).**

| change | Δα | Δβ |
|---|---|---|
| gM = 2.87 | −0.08 | +0.15 |
| gM = 3.07 | +0.08 | −0.16 |
| gA = 1.2754, fπ = 92.07 MeV | ≤ 0.01 | ≤ 0.01 |
| hA (est.) | < 0.003 | < 0.003 |

- **Why the constants hardly matter.** The constants act mostly on the static polarisabilities, which are fitted;
  without the covariance, F1 shifts by Δα = +0.17 in the prediction but ≤ 0.01 in the fit.
- **Why hA is estimated rather than varied.** hA is a COMMON initialiser in `ampl.f90`, not a namelist variable,
  so it cannot be varied at run time. The πΔ loops scale as hA², and the Γ_Δ-driven uncertainty of hA is about
  1%; hence 2% × the πΔ-loop non-static shift.

### 4.6 Error budget

**Main fit (F2 + trunc. cov.)**

| source | δαE1 | δβM1 |
|---|---|---|
| statistics (stat ⊕ p2p syst) | 0.12 | 0.12 |
| floating normalisations | 0.30 | 0.28 |
| forward points (nuisances λ_f, λ_g) | 0.17 | 0.06 |
| EFT truncation (covariance in χ²) | 0.30 | 0.30 |
| *total fit error* | *0.47* | *0.43* |
| truncation-model choice (Λ_b, δ, c̄, ℓ; rms over 7 variants) | 0.15 | 0.07 |
| database selection (half range of d, e, f, g, h, F4) | 0.14 | 0.12 |
| spin polarisabilities (γπ free, F2π) | 0.20 | 0.16 |
| γNΔ coupling gM ± 0.1 | 0.08 | 0.15 |
| constants gA, fπ → PDG | 0.01 | 0.01 |
| hA (Δ width) | 0.00 | 0.00 |
| forward-error model (½ \|nuisance − uncorrelated\|) | 0.08 | 0.02 |
| **total** | **0.56** | **0.51** |

**F1 (Baldin) variant.** α = 11.41 and β = 2.75.

| source | δα | δβ |
|---|---|---|
| stat | 0.13 | 0.16 |
| normalisations | 0.32 | 0.27 |
| Baldin | 0.14 | 0.05 |
| truncation | 0.30 | 0.30 |
| *total fit error* | *0.45* | *0.43* |
| truncation model | 0.09 | 0.11 |
| database | 0.15 | 0.14 |
| spin | 0.18 | 0.17 |
| gM | 0.14 | 0.14 |
| **total** | **0.53** | **0.51** |

Not included in the totals:

- **No-HIγS shift** (−0.38/+0.32). It is reported separately because it is a consistency test rather than a
  selection choice.
- **Conservative next-order estimate** (±0.4/±0.8).
- **Central-value shift between the fits with and without the truncation covariance** (+0.53/−0.50). The
  covariance fit is the principled one, so this shift is a diagnostic rather than an error.

**LaTeX (headline):**

```latex
\begin{tabular}{lcccc}
\hline
 & $\alpha_{E1}$ & $\beta_{M1}$ & $\alpha_{E1}+\beta_{M1}$ & $\chi^2/{\rm dof}$\\
\hline
This work (main)  & $11.69\pm0.47_{\rm fit}\pm0.31_{\rm syst}$ & $2.88\pm0.43_{\rm fit}\pm0.27_{\rm syst}$ & $14.57\pm0.27$ & $241.2/262$\\
This work (F1, Baldin) & $11.41\pm0.45\pm0.28$ & $2.75\pm0.43\pm0.29$ & $14.16\pm0.19$ & $236.1/245$\\
\hline
\end{tabular}
% fit error = stat(0.12) (+) norm(0.30/0.28) (+) forward(0.17/0.06) (+) EFT truncation(0.30/0.30)
```

### 4.7 Comparison with other extractions

| | αE1 | βM1 | comment |
|---|---|---|---|
| **this work, main** | 11.7 ± 0.6 | 2.9 ± 0.5 | BχPT NNLO + O(p⁴) LECs, world data + forward points, trunc. cov. |
| this work, F1 (Baldin) | 11.4 ± 0.5 | 2.8 ± 0.5 | |
| this work, no truncation covariance (F1) | 10.89(29) | 3.36(24) | directly comparable to LM14/MPG13 |
| NNLO prediction (LMP15) | 11.2(7) | 3.9(7) | |
| Lensky–McGovern 2014 (BχPT) | 10.6(5) | 3.2(5) | |
| McGovern–Phillips–Grießhammer 2013 (HBχPT) | 10.65(50) | 3.15(50) | |
| Mornacchi 2022 (A2), average of 3 models | 10.99(16)(47)(17)(34) | 3.14(21)(24)(20)(35) | |
| Mornacchi 2022 (A2), BχPT only | 10.65(16)(47)(4) | 3.28(21)(24)(9) | thesis 10.66(15)(41), 3.40(20)(22) |
| HIγS 2022 (HBχPT, BSR) | 13.8(1.2)(0.1)(0.3) | 0.2(1.2)(0.1)(0.3) | our BχPT, HIγS alone: 13.77(1.16), 0.24(1.16) |
| Pasquini–Pedroni–Sconfietti 2019 (DR) | 12.03(+48/−54) | 1.77(+52/−54) | |
| Mornacchi–Rodini–Pasquini–Pedroni 2022 (DR bootstrap) | 12.7(8) | 2.4(6) | α+β = 15.1(7) |
| PDG 2024 | 11.2(4) | 2.5(4) | |
| PDG 2025/26 | 11.5(4) | 2.31(29) | |

**Assessment.**

- **Main result.** It lies between the earlier χPT fits (β ≈ 3.2) and the DR fits (β ≈ 1.8–2.4). It agrees
  within 1σ with the NNLO prediction, with both PDG averages, and with MRPP22.
- **Comparison with the earlier BχPT fits.** In the "traditional" analysis without a truncation error, our
  result (10.9/3.4) reproduces LM14/MPG13. The shift towards larger α and smaller β comes from including the EFT
  truncation error, which de-weights the 120–150 MeV region.

---

## 5. Δ region (150–350 MeV)

### 5.1 χ² per data set

χ²/N with the normalisation fitted (penalty); the χ² with n = 1 is in `results/tables.md`.

| set | N | NNLO prediction | main-fit α, β (spin NNLO) | n (NNLO) |
|---|---|---|---|---|
| A2 Σ2x (Martel, 288 MeV) | 4 | 0.63 | 0.67 | – |
| A2 Σ2z (Paudyal, 275/295 MeV) | 10 | 1.42 | 1.61 | 0.90 |
| LEGS Σ3 (all; < 2π threshold) | 79; 65 | 1.86; 2.09 | 2.52; 2.89 | – |
| LEGS dσ | 77 | 6.97 | 7.92 | 1.04 |
| MAMI 97 75° c.m. (Hünger) | 12 | 1.54 | 1.59 | 0.95 |
| MAMI 97 90° c.m. (Hünger) | 9 | 2.27 | 2.60 | 1.02 |
| MAMI 99 (Wissmann) | 8 | 1.13 | 1.26 | 1.00 |
| MAMI 02 (Camen/Kossert) | 7 | 0.87 | 1.15 | 0.98 |
| SAL 93 (Hallin) > 150 MeV | 70 | 3.30 | 3.64 | 1.17 |
| TAPS 157–164 MeV | 10 | 1.36 | 1.73 | 1.03 |
| Bonn 72 (Genzel) / Illinois 67 (Gray) / Moscow 66 / Kolb / Nagashima | 14/5/10/1/1 | 2.9 / 4.9 / 12.6 / 6.7 / 8.4 | similar | |

**Sensitivity to gM.** The parameter-free χ² as a function of gM (`results/delta_region.json`, `gM_scan`) has
these minima:

- Hünger 75° ≈ 2.90;
- Hünger 90° ≈ 3.00;
- Camen ≈ 2.95;
- Wissmann ≈ 2.95–3.00;
- LEGS dσ ≈ 3.03 (but χ² ≥ 670);
- Paudyal ≈ 2.85;
- LEGS Σ3: χ² still decreasing at 3.5;
- Martel Σ2x: insensitive.

Fits of gM alone:

| data | gM | χ²/dof |
|---|---|---|
| MAMI dσ | 2.912(26) | 54/33 |
| A2 Σ2x+Σ2z | 3.23(18) | 16/14 |
| LEGS Σ3 | 3.41(5) | 71/66 |
| Hallin | 2.92(4) | 139/55 |
| LEGS dσ | 2.82(2) | 318/77 |

The cross-section sets mildly favour gM slightly below 2.97, as anticipated in LMP15 ("rather high around the Δ
peak"); LEGS Σ3 pulls the opposite way. **These numbers are limited by the missing γNΔ vertex corrections**
(LMP2015 Sect. 2.2): Σ2x is off by 0.03–0.06 at 286–334 MeV. They are exploratory only.

### 5.2 Exploratory spin-polarisability fits

Setup: α−β from the main fit, α+β = 14.0(2), γ0 = −0.929(105).

| data / constraints | γE1E1 | γM1M1 | γE1M2 | γM1E2 | γπ | χ²/dof |
|---|---|---|---|---|---|---|
| A2 Σ2x+Σ2z; γE1E1, γM1M1 free | −3.59(1.15) | 3.18(1.15) | 0.18* | 1.15* | 7.7(2.3) | 18.6/12 |
| A2 Σ2x+Σ2z + γπ = 8.0(1.8); 4 free | −3.68(99) | 2.58(88) | 0.07(1.19) | 1.95(88) | 8.1(1.8) | 16.7/11 |
| + LEGS Σ3 (< 2π), 4 free | −3.59(63) | 3.83(79) | −3.32(1.26) | 4.02(82) | 14.8(3.1) | 81.1/75 |
| **+ LEGS Σ3, γπ = 8.0(1.8), 4 free** | **−2.96(52)** | **2.58(41)** | **−1.39(76)** | **2.71(43)** | **9.6(1.6)** | 84.6/76 |
| same + vertex theory error (0.03/0.02/0.01) | −2.93(55) | 2.56(41) | −1.39(78) | 2.69(44) | 9.6(1.6) | 78.5/76 |
| same + gM floated, prior 2.97(10) → gM = 2.98(10) | −2.95(53) | 2.61(46) | −1.40(77) | 2.67(49) | 9.6(1.6) | 84.6/76 |
| same + gM floated, prior 2.97(30) → gM = 3.07(26) | −2.87(58) | 2.81(72) | −1.44(80) | 2.45(76) | 9.6(1.6) | 84.5/76 |

With A2 Σ2x+Σ2z alone and all four γ's free, the system is under-determined (γπ runs to 46 ± 10); it is shown in
the tables only for completeness.

Comparison:

| source | γE1E1 | γM1M1 | γE1M2 | γM1E2 |
|---|---|---|---|---|
| Paudyal 2020, BχPT | −2.65(43) | 2.43(42) | −1.32(72) | 2.47(42) |
| Paudyal 2020, HDPV | −3.18(52) | 2.98(43) | −0.44(67) | 1.58(43) |
| Martel 2015, DR | −3.5(1.2) | 3.16(85) | −0.7(1.2) | 1.99(29) |
| Martel 2015, BχPT | −2.6(8) | 2.7(5) | | |
| MRPP 2022 | −3.0(6)(4) | 3.7(5)(1) | −1.2(1.0)(3) | 2.0(7)(4) |
| NNLO prediction | −3.3(8) | 2.9(1.5) | 0.2(2) | 1.1(3) |

- **Agreement with Paudyal.** The γπ-constrained A2+LEGS fit agrees with Paudyal's BχPT fit within 0.3σ, so the
  machinery reproduces the A2 analysis.
- **Floating gM.** Within a reasonable prior, floating gM changes the spin polarisabilities by ≤ 0.2 (wide prior:
  γM1M1 ± 0.7).
- **Without γπ.** Removing the γπ constraint lets LEGS Σ3 (mean pull +0.9 at 80–120° c.m.) drive γπ to about 15.
  This is the Σ3 analogue of the known LEGS–MAMI γπ discrepancy.
- **Consistency with region L.** The sub-threshold F3 fit, γM1M1 = 3.00(48), γπ = 7.2(1.0), is consistent with
  the Δ-region results.

---

## 6. Figures (`work/fits/figs/*.pdf`, style `work/figstyle.mplstyle`)

| file | caption (draft) |
|---|---|
| `fig1a_cs_A2_HIGS.pdf` | dσ/dΩ_lab vs θ_lab at the A2 energies 92.2, 103.2, 113.6, 124.1, 135.3 MeV and the HIγS energy 81.3 MeV, with all data within ±2.5 MeV (distinct markers per experiment; c.m. data converted to lab). Dashed: parameter-free NNLO. Solid: main fit (F2) with the EFT truncation band (68%, scheme B). Dash-dot: F1 (Baldin). Data shown without fitted normalisations. |
| `fig1b_cs_older.pdf` | Same at the TAPS energies 58.9, 79.2, 108.4, 149.9 MeV with data within ±5 MeV (curves at the nominal energy). |
| `fig2_sigma3.pdf` | Beam asymmetry Σ3 at 92.2, 108.5, 129.5 MeV (A2 2022) and at 83.4 MeV. The 83.4 MeV values are derived from the HIγS dσ∥ and dσ⊥ (two in-plane detectors combined, normalisation cancels). Curves as in Fig. 1. |
| `fig3_higs_linear.pdf` | HIγS dσ∥ and dσ⊥ at 83.4 MeV (left) and unpolarised dσ at 81.3 MeV (right). Curves: NNLO, main fit, main fit without HIγS, and the HIγS-only fit (α = 15.2, β = 2.0; × its fitted normalisations). |
| `fig4_alpha_beta.pdf` | α–β plane: 68% joint-CL ellipses (Δχ² = 2.30) of the main fit, F2 without truncation covariance, F1 with/without, main configuration without HIγS, and A2-only, HIγS-only, pre-2015-only. Baldin band 14.0(2); NNLO prediction; other extractions (MPG13, LM14, A2 2022, HIγS 2022, PPS19, MRPP22, PDG 2024, PDG 2025/26). |
| `fig5a_delta_S2x_S2z.pdf` | Δ region: Σ2x at 288 MeV (Martel) and Σ2z at 275/295 MeV (Paudyal) vs θ_lab. NNLO with the main-fit α, β (dashed), gM = 2.85–3.10 band, spin-polarisability fits (A2+LEGS+γπ constraint, solid; A2 only with γE1E1, γM1M1 free, dash-dot). Missing γNΔ vertex corrections: Σ2x uncertain by about 0.03–0.06. |
| `fig5b_delta_LEGS_S3.pdf` | LEGS Σ3 vs θ_cm at 254, 286.5, 310 MeV with NNLO, gM band and the spin-polarisability fit. |
| `fig5c_delta_cs_vs_E.pdf` | dσ/dΩ vs E_γ up to 350 MeV at θ_cm = 75°, 90°, 115° (c.m.) and θ_lab ≈ 131–136° (lab), with all data within ±3°: Hünger, LEGS, Hallin, Wissmann, Camen/Kossert, Genzel, Gray, Baranov 66, Kolb, Nagashima, TAPS. NNLO prediction, main-fit α/β with gM = 2.97 and 2.90, gM band 2.85–3.10. |
| `fig6a_pulls_main_F2.pdf` | Pulls (y − n t)/σ_pt of the main fit per data set (labels: χ²/N), and the pull histogram of the angular data compared with N(0,1). |
| `fig6b_pulls_F1.pdf` | Same for F1. |
| `fig7_forward.pdf` | GHP forward dσ/dΩ_lab(0°) and Σ2z(0°) vs E_γ (filled: used, E ≤ 135 MeV; open: 145/150 MeV, not used), NNLO prediction, main fit with truncation band, F1. Dotted line: π⁰ threshold. |
| `fig8_Emax_scan.pdf` | Fitted α and β vs the upper energy cut E_max for F1/F2, with and without the truncation covariance. |

---

## 7. Physics conclusions

1. **The updated database is consistent with NNLO BχPT + O(p⁴) LECs.** The new A2 dσ/Σ3 and HIγS data and the
   GHP forward points all fit with χ²/dof ≈ 0.9 once the EFT truncation error is included. Normalisations stay
   within their priors.
2. **The result moves β down and α up relative to the 2014 BχPT fit.** The main result is α = 11.7(6),
   β = 2.9(5). The move is caused mostly by the theory treatment (the truncation covariance de-weights the
   120–150 MeV region, where the fitted β depends strongly on the energy window). The new data contribute less:
   A2 alone favours a larger β, while HIγS and TAPS favour a smaller one.
3. **α+β from Compton data within BχPT is 14.6(3).** This holds both with the forward points (main) and without
   any sum constraint (15.4(6)). It is 1.5–2σ above the Baldin/GHP value 14.0(2). The forward dσ(0°) points alone
   reproduce 14.1(1), so the tension is between the angular data and the sum rule within this amplitude.
   - It mirrors MRPP22's 15.1(7).
   - It is not resolved by any database variant.
4. **Spin polarisabilities.** Below threshold, the data constrain γπ = 7.2(1.0) and γM1M1 = 3.0(5), consistent with
   NNLO and with the dispersive MAMI value. The forward Σ2z points with the truncation covariance give
   γ0 = −1.07(30). In the Δ region, the A2 double-polarisation data plus LEGS Σ3 (with γπ = 8.0(1.8)) reproduce
   Paudyal's BχPT values.
5. **Tensions.**
   - HIγS vs the rest: PG test p = 0.07 (main), at most 2.2σ in the strictest version.
   - A2-only vs HIγS-only: 2.3σ in β.
   - TAPS-only gives β = 1.8(8) and contributes the only 3σ outlier.
   - LEGS dσ and Σ3 are inconsistent with MAMI and with the code in the Δ region.

---

## 8. Caveats and issues found (reported, not silently fixed)

1. **Code (work/code).**
   - **Single pion mass in the loops.** The loops use one (charged) pion mass, so the cusp sits at 149.95 MeV
     instead of 144.7 MeV (π⁰p). This is why the forward points at 145/150 MeV are unusable, and it limits the
     region near 140–150 MeV.
   - **No γNΔ vertex corrections.** This matters for Σ2x and for gM determinations from Δ-region data.
   - **hA cannot be set at run time.** It is a COMMON initialiser, so its effect was estimated rather than
     varied. Suggest exposing it in a namelist in `code_v2`.
2. **Forward points.**
   - Valid within the code only up to 135 MeV.
   - The table rounding (two decimals) sets the uncorrelated part of the nuisance model, so any NNLO shape
     deviation larger than about 0.003 in Σ2z is penalised. A theory error on the forward points is therefore
     indispensable.
   - GHP's 145/150 MeV points lie between the physical π⁰p and π⁺n thresholds.
3. **Truncation-error model.**
   - Scheme B is our construction (LEC-matched, Born + π⁰ exact). The naive scheme A double counts the fitted
     LECs.
   - The GP length scales fitted to the coefficients are short (12–26 MeV), because the matched coefficients
     grow as E³–E⁴ below mπ, i.e. they are non-stationary.
   - Results vary by up to 0.36 (α) / 0.16 (β) across the length-scale variants (ℓ = 100/100 MeV/deg being the
     extreme); this enters the budget as the rms (0.15/0.07).
4. **Database assumptions.**
   - **Old-set normalisations.** The 10% for Pugh, Hyman, Oxley and Bernardini is an assumption (the sources state
     none).
   - **Pugh and Hyman.** The true energies are reconstructed from the bin labels, two of them inferred. The KLP
     angle "50°" for Pugh/Hyman disagrees with GMPF's "45°" for Pugh; unverifiable.
   - **Bernardini.** The energies are bin labels, and KLP's c.m. → lab conversion was done at the label energy.
   - **Federspiel.** The legacy values are 3–5% below the published figure (db_old_README).
   - **Frames.** The Hallin frame (c.m.) and the Wissmann/Camen frame (lab) are inferred, and the Gray, Genzel,
     Baranov 1966 and Nagashima frames are assumed c.m. (affects only region D).
   - None of these sets dominates. The pruning variant (f), which drops Oxley and Bernardini, moves the main fit
     by −0.23/+0.20.
5. **Normalisation of asymmetries.** The A2 Σ3 5% uniform scale is approximated by a Gaussian with rms 2.89%; the
   Paudyal 10.4% scale is floated (fitted n = 0.90).
6. **No issues found in the A2/HIγS tables.**
   - The HIγS dσ∥/dσ⊥ conventions are verified by reproducing Li et al.'s HBχPT numbers with our BχPT within
     0.1–0.2.
   - The Martel/Paudyal lab-angle frame is confirmed by χ².
7. **Spin polarisabilities fixed at NNLO in the main fit.** Their uncertainty is included through the γπ-free
   variant (±0.20/±0.16). Fixing them at the MRPP22 values instead would shift β by −0.40.

---

## 9. Files

- **Code (`work/fits/`).**
  - `Makefile`, `run_all.sh`;
  - `lib/{bchpt,database,forward,theory,fitter,common}.py`;
  - `scripts/s00_validate_forward.py` … `s07_tables.py`.
- **Database (`work/fits/db/`).** `fitdb.csv` and `fitdb_pb.csv`.
- **Results (`work/fits/results/`).**
  - `tables.md` and `tables.tex`: all tables;
  - `dataset_table.{md,tex,csv}`;
  - `fits_L.json`: every fit with P, covariance, χ², per-set χ², normalisations, residuals and pids;
  - `delta_region.json`;
  - `error_budget.json`;
  - `theory_lengths.json`;
  - `forward_validation.md`;
  - `s0x_log.txt`.
- **Figures.** `work/fits/figs/*.pdf` (12 files).
