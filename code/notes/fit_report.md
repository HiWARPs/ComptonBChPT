# Proton polarisabilities from the updated world Compton database in NNLO BχPT — fit report

*Analysis for the single-author paper (V. Pascalutsa). Date: 2026-10-14.*

- **Baseline: code v2** (`work/code_v2`), which includes the γNΔ vertex corrections and physical pion and nucleon
  masses in the πN loops. The new results are in Part V.
- **Reproduction:** `cd work/fits && make` (v2 profile, about 40 s with a filled cache).
- **Tables:** all v2 tables are in `results/tables.md` and `results/tables.tex`.
- **MAIN-FIT SWITCH (2026-10-14, after the v2 baseline below was written): the main fit is now the one with the
  3σ outlier removed.** The main database ('d_all' in all keys) = all region-L data minus `MAIN_DROP =
  ('OlmosdeLeon2001_029',)` (TAPS, 108.42 MeV, 133°, pull −4.27); the all-points database is the variant 'd_allpts'.
  The new main result is α = 11.59(47)(29), β = 2.94(43)(27), χ² = 219.7/263 (Sect. V.9); the numbers quoted in
  Sects. V.1–V.6 are those of the all-points fit (now 'd_allpts', formerly 'd_all') and are superseded where they
  differ — the complete OLD → NEW list is in `work/notes/mainfit_switch.md`.

- **MAIN-DATABASE CHANGE (2026-10-17, author decision): the HIγS 2022 data (Li et al.) are no longer part of the
  main database.** 'd_all' = pre-2015 world data + A2 2022 minus MAIN_DROP (245 angular points + 22 forward
  points, N = 267); HIγS is kept as a test and as the variant 'j_withHIGS' (= the former main, 283 points). The
  new main result is α = 11.16(48)(30), β = 3.30(44)(28), α+β = 14.46(25), α−β = 7.86(89), χ² = 201.1/249
  (= the former variant "without HIγS"); every OLD → NEW number, the code changes and the regenerated figures are
  listed in `work/notes/main_without_HIGS.md`. The numbers quoted below (Sects. V.1–V.9) are those of the
  former main database with HIγS and are superseded where they differ.

- **FORWARD POINTS 2026 (2026-10-18): the forward points come from the 2026 ForwardRCS evaluation (FWD_SOURCE=2026).** Current
  main result α = 10.89(47)(39), β = 3.23(44)(30), α+β = 14.13(19), α−β = 7.66(89), χ² = 198.0/249; PWA α−β = 8.05(1.22)(1.02),
  β = 3.08(61)(52). Every OLD → NEW number is in `work/notes/forward2026_update.md`; the team audit of 2026-10-18 is in
  `work/notes/audit_2026-10-18_team.md`. Variant (k) 'main − bremsstrahlung data of 1957–1960' and the HIγS free-normalisation
  fits were added to s03 on 2026-10-18 (see `results/tables.md`).

Units: α, β in 10⁻⁴ fm³; γ in 10⁻⁴ fm⁴.

---

# Part V — Results

## V.0 Code configuration (differences from the 2015 code)

### Code (main configuration)

- **Settings.** `dvertex=t` (imaginary parts of the γNΔ vertex loops, LMP2015 Fig. 4), `pimassphys=t`,
  `nmassphys=t`, `consts='LMP15'`. The physical masses mean π⁰p loops with m_π⁰ and π⁺n loops with m_π± and M_n,
  so the cusps sit at the physical thresholds, 144.7 and 151.4 MeV.
- **Bug fixes.**
  - The πΔ Thomson-limit leak is fixed.
  - θ = 0° and 180° are now exact: the loop amplitudes are extrapolated in t at fixed s. the 2015 code extrapolated the
    observable instead, which was off by up to 5.6e-4.
  - Analytic reference polarisabilities.
  - LoopTools state leak fixed.
- **New options** used here: the order switch (`order=2|3|3.5`), h_A as a namelist variable, and PDG-2024
  constants.

### NNLO static polarisabilities of the v2 base amplitude

The physical masses add a numerically extracted shift to the πN-loop part:
(α, β, γE1E1, γM1M1, γE1M2, γM1E2) = (11.415, 3.729, −3.246, 2.939, 0.170, 1.133), with γ0 = −0.996 and
γπ = 7.15. The 2015-code values were 11.246, 3.871, −3.335, 2.943, 0.179, 1.147, with γ0 = −0.933. In all F1/F2 fits the
spin polarisabilities are fixed at these v2 values.

### Pipeline

- `FIT_PROFILE=v2` (default) selects the executable and settings (`lib/bchpt.py`).
- Database, χ² machinery, nuisance model, database variants and truncation scheme are unchanged. The fixed
  references (`PREF_NNLO`) and the base static polarisabilities are now taken from the driver's `Pref`, which in
  v2 follows the mechanisms switched on.

### Forward points

The main fit now uses the GHP points up to **145 MeV** (11 dσ(0°) + 11 Σ2z(0°), N = 284). Variants: ≤ 135 MeV
(without 145 and 150) and ≤ 150 MeV.

At 145 MeV the v2 NNLO prediction is 2.38 nb/sr against GHP 2.77(16) (2015 code: 2.17). In the main fit the 145 MeV points
have pulls of +0.05 (dσ) and −0.7 (Σ2z), relative to exp ⊕ truncation errors.

### Truncation errors from the v2 order switch

The scheme-B coefficients were recomputed with `order=2|3` (LO/NLO) and the Born + π⁰ configuration (flags at
order 3.5).

- **Consistency with the flag-based orders of the 2015 code.** `order=2` ≡ my LO flags and `order=3` ≡ my NLO flags, to all
  printed digits (8 test kinematics). `Pref` follows the order.
- **Size compared with the 2015 code.** The ratio σ_th(this code)/σ_th(2015 code) has median 0.995 over the 284 fit points. Per-set medians
  are 0.986–1.005; individual points range 0.91–1.09, at the cusp. corr(c3^2015, c3^now) = 0.9997.
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
| 2015 code, same fit | 11.69 ± 0.47 ± 0.31 | 2.88 ± 0.43 ± 0.27 | 14.57(27) | 8.81(86) | 241.2/262 (0.82) |

- **Correlation.** ρ(α, β) = −0.84.
- **Nuisance pulls.** λ_f = −0.61(42), λ_g = −0.01(17). Normalisation pulls are within ±1.4.
- **Spin polarisabilities.** γ0 = −1.00 is the fixed NNLO input (v2 physical masses).
- **F3 (main + γM1M1, γE1E1 free).** γM1M1 = 3.03(49), γE1E1 = −3.25(56), γ0 = −1.08(29) (fixed by the forward
  Σ2z(0°) points), γπ = 7.25(1.01). With γπ free at fixed γ0 (F2π): γπ = 7.32(97).
- **F1 (Baldin 14.0(2) penalty, no forward points) + trunc. cov.** α = 11.41 ± 0.45 ± 0.30, β = 2.76 ± 0.43 ± 0.27,
  χ² = 236.9/245. Without the truncation covariance (LM14-comparable): 10.87(29), 3.40(24).

### Error budget (v2 main)

| source | δα | δβ | 2015-code δα | 2015-code δβ |
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
- **Verdict.** Unchanged from the 2015-code analysis: a mild tension of about 2σ.

## V.4 E_max scan: does the β drift shrink with the physical cusps? — **No.**

| E_max [MeV] | 90 | 100 | 110 | 120 | 130 | 140 | 150 | 170 |
|---|---|---|---|---|---|---|---|---|
| β, F1 no trunc., 2015 code | 1.34(69) | 1.79(56) | 2.01(44) | 2.43(36) | 3.11(29) | 3.37(24) | 3.36(24) | 3.45(24) |
| β, F1 no trunc., **v2** | 1.32(69) | 1.78(56) | 2.00(44) | 2.43(36) | 3.13(28) | 3.40(24) | 3.40(24) | 3.49(23) |
| β, main + trunc., 2015 code | 1.52(71) | 1.94(60) | 2.12(52) | 2.48(49) | 2.73(45) | 2.91(44) | 2.88(43) | 2.94(43) |
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
  is 0.91 against GHP 1.12(9); with the 2015 code it was 0.18. That is much closer, but within the fully correlated model the
  steep threshold dependence still gives the point excessive leverage. It is kept as a variant only.
- **Without the truncation covariance**, the forward Σ2z(0°) shape problem is reduced but persists: the
  Σ2z(0°) χ² for the same ≤ 135 MeV points falls from 328 (2015 code) to 198, still unacceptable.

## V.6 Δ region with γNΔ vertex corrections

### χ²/N of the parameter-free NNLO prediction (normalisation fitted)

χ² with n = 1 in brackets.

| set | N | 2015 code | without vertex corr. | **main** |
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
(151 → 146). They worsen LEGS dσ (646 → 684). All other changes are ≤ 2%. This is consistent with the code notes
(Sect. C): the vertex corrections mostly fix Σ2x/Σ2z shapes at 300+ MeV, where there are few data.

### gM fitted to single data sets

Spin polarisabilities at NNLO; α, β constrained to the main fit; prior 2.97(50).

| data | 2015 code | without vertex corr. | **main** |
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
| 2015 code | −2.96(52) | 2.58(41) | −1.39(76) | 2.71(43) | 9.6(1.6) | 84.6/76 |
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
- **Code.** `lib/bchpt.py` (profile `v2`), `Makefile` (`FIT_PROFILE`, `FIT_RESULTS`, `FIT_FIGS`).

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
