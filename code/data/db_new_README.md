# db_new.csv: new (2014–2026) proton RCS data from A2@MAMI and HIGS

Compiled 2026-10-13. The database is `db_new.csv` (126 rows). It is built by
`scripts/build_db_new.py` and cross-checked by `scripts/verify_db_new.py`. Run both with
`python3 -I <script> <release>/code/data`.
Raw downloads are in `raw/new/<dataset>/`. Plain-text extractions of the PDFs (pdftotext -layout)
are in `extracted/`.

## Column conventions

| column | meaning |
|---|---|
| dataset | unique tag per independent data set (separate normalisation) |
| texkey | INSPIRE texkey (first/collaboration key) |
| obs | `CS` = unpolarised dσ/dΩ; `S3` = Σ3; `S2X` = Σ2x; `S2Z` = Σ2z; **extra:** `CS_PAR` / `CS_PERP` = dσ/dΩ for a linearly polarised beam whose E-vector is parallel / perpendicular to the scattering plane (HIGS only) |
| E_lab_MeV | lab photon energy: bin centre (A2 tagged bins) or intensity-weighted mean beam energy (HIGS) |
| dE_MeV | full bin width, upper − lower edge (A2). For HIGS it is the **FWHM** of the beam energy spread |
| theta_lab_deg, theta_cm_deg | photon scattering angle. All sources quote **lab** angles. θ_cm is computed with exact 2-body kinematics at E_lab, M_p = 938.272 MeV: cosθ_cm = (cosθ_L − β)/(1 − β cosθ_L), β = E/(E+M_p). The original frame is stated in `notes` |
| frame_of_cs | frame of dσ/dΩ as quoted (`lab` for all CS here); `na` for asymmetries |
| value | transcribed verbatim, with no conversion |
| stat_err | statistical error. Where the source gives a combination, `notes` says so (e.g. Paudyal "Rand." = stat ⊕ random syst; Martel = stat ⊕ pol. syst) |
| sys_err_pt | absolute systematic error per point, as tabulated. It is point-to-point uncorrelated unless `notes` says otherwise |
| norm_err_frac | correlated fractional scale error quoted in the text. **Overlap flags in notes:** `NORM_ALREADY_IN_STAT_ERR` or `NORM_ALREADY_IN_SYS_ERR_PT` mean the scale error is already contained in that column, so do not add it twice |
| units | `nb/sr` for cross sections, `1` for asymmetries |
| source | exact table and file |
| notes | bin edges, original frame, details of the error model |

## Summary

| dataset | texkey | obs | E_lab (MeV) | θ_lab (θ_cm) | N | scale error | provenance |
|---|---|---|---|---|---|---|---|
| A2_Martel2015_S2x | A2:2014iky | S2X | 273–303 (288) | 90–150° (103.6–156.2°) | 4 | P_T 2.1% (already inside error bar) | Martel PhD thesis Tab. 7.3 "Summed" = PRL Fig. 3 |
| A2_Paudyal2020_S2z | A2:2019bqm | S2Z | 265–285, 285–305 | 87.5–147.5° (100.7–154.3°) | 10 | 10.4% (target 10% ⊕ beam 2.7%), already inside "Syst." | PRC 102, 035205 Table I |
| A2_Mornacchi2022_CS | A2:2021wmm | CS | 86.3–140.4 (5 bins) | 35–145° (38.1–149.0°) | 60 | 3% | PRL 128, 132503 Suppl. Table I |
| A2_Mornacchi2022_S3 | A2:2021wmm | S3 | 86.3–140.4 (3 bins) | 35–145° (38.1–148.8°) | 36 | 5% uniform (rms 2.9%) | PRL 128, 132503 Suppl. Table II |
| HIGS_Li2022_CIRC | Li:2022vnz | CS | 81.3 | 55, 90, 125° (58.8–128.7°) | 8 | 2.24% (flux 2% ⊕ target 1%) | PRL 128, 132502 Table I |
| HIGS_Li2022_LIN | Li:2022vnz | CS_PAR (5), CS_PERP (3) | 83.4 | 55, 90, 125° (58.9–128.7°) | 8 | 2.24% | PRL 128, 132502 Table I |
| Sokhoyan 2017 Σ3 | A2:2016nio | S3 | 79–139 | 30–~155° | 21 | – | **figures only, not in the db** |

---

## 1. Martel et al. (A2), Σ2x in the Δ region

* **Reference.** P. P. Martel et al. (A2 Collaboration at MAMI), "Measurements of Double-Polarized Compton
  Scattering Asymmetries and Extraction of the Proton Spin Polarizabilities", Phys. Rev. Lett. 114, 112501 (2015).
  arXiv:1408.1576. INSPIRE texkeys A2:2014iky and Martel:2014pba.
* **Provenance.** The PRL shows Σ2x only in Fig. 3 (4 points). The numbers are tabulated in the PhD thesis the PRL
  cites as its analysis reference (P. P. Martel, UMass Amherst, 2013, "Measuring Proton Spin Polarizabilities
  with Polarized Compton Scattering"):
  https://download.uni-mainz.de/fb08-kpha2/publications/phd/pmartel_dis.pdf, Table 7.3, column "Summed"
  (φ-summing method). Local copies: `raw/new/Martel2013thesis/pmartel_dis.pdf` and `extracted/Martel2013_thesis.txt`.
  * The PRL figure file `raw/new/Martel2015/src/Sig2x_gE1E1.eps` (arXiv source) carries the internal title
    "AsySum.eps: Summed Asymmetry".
  * Its vector coordinates reproduce the thesis "Summed" values and errors to ≤0.0014, which is the EPS
    resolution. This check is `[M1]` in `verify_db_new.py`. The EPS was used only to verify the thesis numbers,
    not to digitise anything.
  * The thesis φ-fit alternative values are given in `notes`.
* **Kinematics.** E_γ = 273–303 MeV, one bin. θ_lab bins are 80–100, 100–120, 120–140 and 140–160°
  (the PRL axis reads "Compton θ_lab").
* **Definition (PRL Eq. 3, quoted).**
  Σ2x(E_γ,θ_c) = [P_T P_γ(E_γ) cos(φ_c)]^(-1) · [N^R − N^L]/[N^R + N^L].
  * φ_c is the azimuth of the scattered photon relative to the target-polarisation direction.
  * The target is polarised transversely, along x in the scattering plane on the side of the scattered photon
    when φ_c = 0.
  * R/L is the right/left photon helicity.
* **Uncertainties.**
  * The error bar is statistical plus the beam- and target-polarisation systematics. The PRL says these
    systematics are at most 9% of the total. Average P_T = 81.6 ± 1.7% and P_e = 81.9 ± 0.1%.
    `norm_err_frac` = 0.021 is flagged `NORM_ALREADY_IN_STAT_ERR`.
  * The carbon-subtraction systematic is not in the error bar (at most about 10% of the total).
  * The carbon subtraction was varied by ±20%.
* **Polarisabilities (PRL Table II), in 10⁻⁴ fm⁴.** The fits use Σ2x plus LEGS Σ3 (Blanpied 2001, below 2π
  threshold), with constraints α+β = 13.8 ± 0.4, α−β = 7.6 ± 0.9, γ0 and γπ.
  * DR fit: γE1E1 = −3.5 ± 1.2, γM1M1 = 3.16 ± 0.85, γE1M2 = −0.7 ± 1.2, γM1E2 = 1.99 ± 0.29.
  * BχPT fit: γE1E1 = −2.6 ± 0.8, γM1M1 = 2.7 ± 0.5.
  * Σ2x alone (DR): γE1E1 = −4.6 ± 1.6, γM1M1 = −7 ± 11.
  * γπ is quoted without the π0 pole.
* **Theory used.** The dispersion model is Drechsel, Pasquini and Vanderhaeghen, Phys. Rep. 378 (2003) 99.
  The BχPT calculation is Lensky and Pascalutsa, EPJC 65 (2010) 195.
* **Not in the db (thesis only, above the 2π threshold).** Thesis Table 7.5, Σ2x at 315–346 MeV, φ-summed values:
  * θ_lab = 70°: −0.160 ± 0.081
  * 90°: −0.168 ± 0.038
  * 110°: −0.202 ± 0.052
  * 130°: −0.308 ± 0.114
  * 150°: −0.074 ± 0.108

## 2. Sokhoyan et al. (A2), Σ3 below the pion threshold: NOT AVAILABLE AS NUMBERS

* **Reference.** V. Sokhoyan et al. (A2 Collaboration), "Determination of the scalar polarizabilities of the proton
  using beam asymmetry Σ3 in Compton scattering", Eur. Phys. J. A 53, 14 (2017). arXiv:1611.03769.
  INSPIRE texkeys A2:2016nio and Sokhoyan:2016yrc.
* **Correction to the task description.** This is not a Δ-region measurement. It is the first Σ3 below the pion
  threshold: E_γ = 79–98, 98–119 and 119–139 MeV, with 7 θ bins each between 30° and ~155°. That gives 21 points
  (χ²/ndf = 19.2/20).
* **Numbers are only in Fig. 4 and Fig. 5. No table exists** in the arXiv source, the EPJA article, or any public
  thesis or proceedings I could find. I checked:
  * Mornacchi's thesis
  * Li's thesis
  * Sconfietti's thesis
  * the arXiv sources of Krupina et al. 2018, Pasquini et al. 2019 and Mornacchi–Pasquini–Pedroni 2022
  * the EPJ Web Conf. 241 proceedings
  * HEPData

  Mornacchi–Pasquini–Pedroni (PRL 129, 102501) did fit these points, presumably with numbers from the collaboration.
  **The points were not digitised, so the set is not in the db.** It has been superseded by Mornacchi 2022 in the same
  energy range.
* **Definition.** Σ3 ≡ (σ∥ − σ⊥)/(σ∥ + σ⊥), where ∥/⊥ means the photon polarisation is parallel/perpendicular to the
  scattering plane. The plotted values are negative (about −0.3 to −0.7).
* **Systematics.**
  * Flux normalisation: 15% at low E, 4% at high E.
  * Background: up to 50% at extreme angles, about 3% at central angles.
  * Linear polarisation: below 2%.
* **Fit, with α+β = 14.0 (Baldin).** βM1 = 2.8 (+2.3, −2.1) in BχPT [Lensky–Pascalutsa] and 3.7 (+2.5, −2.3) in HBχPT
  [McGovern–Phillips–Grießhammer 2013], in 10⁻⁴ fm³.
* **Theory compared.** Born, DR (Drechsel 1999 / Pasquini 2007), BχPT (Lensky–Pascalutsa 2010) and HBχPT (McGovern 2013),
  all with α = 10.65 and β = 3.15.

## 3. Paudyal et al. (A2), Σ2z

* **Reference.** D. Paudyal et al. (A2 Collaboration at MAMI), "Extracting the spin polarizabilities of the proton by
  measurement of Compton double-polarization observables", Phys. Rev. C 102, 035205 (2020). arXiv:1909.02032 (v3).
  INSPIRE texkeys A2:2019bqm and Paudyal:2019mee.
* **Provenance.** Table I, parsed from the arXiv TeX source `raw/new/Paudyal2020/src/PRC_S2z_Paudyal.tex`. It matches
  the arXiv PDF text exactly (check `[P]`).
* **Kinematics.**
  * E_γ bins are 265–285 and 285–305 MeV. Two beamtimes (2014 and 2015) were combined by weighted average.
  * θ_lab centres are 87.5, 102.5, 117.5, 132.5 and 147.5°, each 15° wide. The figure axis is "θ_lab".
* **Definition (quoted).**
  Σ2z = 1/(P^γ_circ · P^t_z) · [(N^R_{+z} + N^L_{−z}) − (N^L_{+z} + N^R_{−z})] / [(N^R_{+z} + N^L_{−z}) + (N^L_{+z} + N^R_{−z})].
  * R/L are the right/left-handed beam helicities and ±z is the target polarisation along or against the beam.
  * The beam polarisation follows from the helicity-transfer (Olsen–Maximon) formula. Average P_e = 86.8 ± 0.1%.
  * At 0° the paper plots the dispersive sum-rule value Σ2z ≈ +0.78 (Gryniuk et al.). This fixes the sign convention
    used in the HDPV and BχPT codes.
* **Uncertainties.**
  * `stat_err` is the "Rand." column: statistical plus random systematic.
  * `sys_err_pt` is the "Syst." column. This is the **correlated** systematic: target polarisation 10%, beam polarisation
    2.7% and carbon subtraction 3–6%, added in quadrature. It is 11–13% of the value.
  * `norm_err_frac` = 0.104 = √(0.10² + 0.027²), flagged `NORM_ALREADY_IN_SYS_ERR_PT`.
  * Recommended treatment: put stat_err in the χ² and use norm_err_frac as a floating scale. The remaining carbon part,
    3–6%, can be treated as a correlated per-point error.
* **Polarisabilities (Table II), in 10⁻⁴ fm⁴.** The fit uses Σ2z, Σ2x (Martel) and Σ3 (LEGS), with
  α+β = 13.8 ± 0.4, α−β = 8.7 ± 0.7, γ0 = −0.929 ± 0.105 and γπ = 8 ± 1.8 (without the π0 pole of −46.7).
  * HDPV: γE1E1 = −3.18 ± 0.52, γM1M1 = 2.98 ± 0.43, γE1M2 = −0.44 ± 0.67, γM1E2 = 1.58 ± 0.43 (χ²/dof = 1.14).
  * BχPT: −2.65 ± 0.43, 2.43 ± 0.42, −1.32 ± 0.72, 2.47 ± 0.42 (χ²/dof = 1.36).
  * Weighted average: −2.87 ± 0.52, 2.70 ± 0.43, −0.85 ± 0.72, 2.04 ± 0.43.
* **Theory compared.**
  * HDPV: Holstein et al. 2000, Drechsel et al. 2003, Pasquini et al. 2007.
  * BχPT: Lensky and Pascalutsa, EPJC 65 (2010).
  * HBχPT: Grießhammer, McGovern and Phillips, EPJA 52 (2016) 139 and EPJA 54 (2018) 37.

## 4. Mornacchi et al. (A2), dσ/dΩ and Σ3 below the pion threshold

* **Reference.** E. Mornacchi et al. (A2 Collaboration at MAMI), "Measurement of Compton Scattering at MAMI for the
  Extraction of the Electric and Magnetic Polarizabilities of the Proton", Phys. Rev. Lett. 128, 132503 (2022).
  arXiv:2110.15691. INSPIRE texkeys A2:2021wmm, A2CollaborationatMAMI:2021vfy and A2:2021rcs.
* **Provenance.**
  * The Supplemental Material (APS copy is paywalled) is open-access at Glasgow Enlighten:
    https://eprints.gla.ac.uk/265979/2/265979Suppl.pdf, saved as `raw/new/Mornacchi2022/supp_glasgow/265979Suppl.pdf`
    (sha256 1248874f…7bc1).
  * Table I has 60 dσ/dΩ points and Table II has 36 Σ3 points, both combining the March and July 2018 beamtimes.
    Both tables were checked visually against the rendered PDF.
  * As an independent check, the stat-weighted mean of the per-beamtime tables in Mornacchi's PhD thesis
    (Mainz 2021, doi:10.25358/openscience-6051, App. A/B) reproduces all 60 CS values to ≤0.0011 nb/sr.
  * It reproduces 34 of the 36 Σ3 values to better than 0.01. The exceptions are 92.2 MeV at 135° and 145°: the
    published values are −0.1424 ± 0.2381 and −0.0741 ± 0.3847, against thesis means of −0.111 ± 0.216 and
    −0.071 ± 0.319. The published supplement is authoritative.
* **Kinematics.**
  * All angles and cross sections are in the **lab frame** ("The cross-section values are given in the laboratory
    frame").
  * θ_lab bin centres run from 35 to 145° in 10°-wide bins.
  * CS energy bins: 86.3–98.2, 98.1–108.4, 108.5–118.7, 118.7–130.2 and 130.3–140.4 MeV.
  * Σ3 energy bins: 86.3–98.2, 98.1–118.7 and 118.7–140.4 MeV.
  * E_lab is the bin central value quoted in the thesis: 92.2, 103.2, 113.6, 124.1, 135.3 for CS, and 92.2, 108.5,
    129.5 for Σ3. These are the bin midpoints except 124.1, where the midpoint is 124.45; midpoints are in `notes`.
* **Definition.** Σ3 = (dσ∥ − dσ⊥)/(dσ∥ + dσ⊥), following Babusci et al. 1998. ∥/⊥ refers to the photon polarisation
  vector relative to the scattering plane. **All values are negative**, which is the expected Thomson-like sign
  below threshold.
* **Uncertainties.**
  * CS: stat, plus point-to-point uncorrelated syst (background, and an extra 3% in the lowest energy bin for
    beamtime inconsistency), plus a correlated **3%** scale. The 3% comes from target density 1%, flux 2% and
    cuts/MC 2%.
  * Σ3: stat, plus point-to-point syst (background), plus a correlated **5% "uniformly distributed"** scale from the
    linear-polarisation determination. `norm_err_frac` = 0.05 is the half-width; the Gaussian-equivalent rms is
    0.05/√3 = 0.0289, which the thesis lists as "corr. 5/√3 %".
* **Polarisabilities (PRL Table I).**
  * The fits use all 96 points plus the Baldin constraint α+β = 13.8 ± 0.4. Spin polarisabilities are fixed to
    Paudyal 2020. The normalisations s_σ and s_Σ are free parameters with a penalty term.
  * HDPV: α = 11.23 ± 0.16 ± 0.46 ± 0.02, β = 2.79 ± 0.20 ± 0.23 ± 0.11.
  * BχPT: α = 10.65 ± 0.16 ± 0.47 ± 0.04, β = 3.28 ± 0.21 ± 0.24 ± 0.09.
  * HBχPT: α = 11.10 ± 0.16 ± 0.47 ± 0.17, β = 3.36 ± 0.21 ± 0.24 ± 0.20.
  * Fitted normalisations: s_σ = 1.011, 1.013, 1.043 and s_Σ = 0.994, 0.996, 1.001 (HDPV, BχPT, HBχPT). χ²/DOF ≈ 0.89 (93 DOF).
  * Average: α = 10.99 ± 0.16 ± 0.47 ± 0.17 ± 0.34 and β = 3.14 ± 0.21 ± 0.24 ± 0.20 ± 0.35
    (stat, syst, spin-pol, model), with ρ(α, β) = −0.75.
* **Theory compared.** HDPV (Drechsel 1999, Holstein 2000, Pasquini 2007), BχPT (Lensky–Pascalutsa 2010) and HBχPT (McGovern 2013).

## 5. X. Li et al. (HIGS), polarised dσ/dΩ around 81–83 MeV

* **Reference.** X. Li et al., "Proton Compton Scattering from Linearly Polarized Gamma Rays", Phys. Rev. Lett. 128,
  132502 (2022). **arXiv:2205.10533** (not 2203.xxxxx). INSPIRE texkey Li:2022vnz.
* **Provenance.** Table I, parsed from the arXiv TeX source `raw/new/Li2022/src/prl_proof_Apr2022.tex`.
  * It matches the arXiv PDF text and Tables 6.5 and 6.7 of X. Li's PhD thesis (Duke 2020, hdl 10161/22194,
    https://inspirehep.net/files/27a37b67bd44333dd21c6a176f969600) exactly (check `[L]`).
  * The thesis also gives the detector names and the breakdown of the point-to-point errors.
* **Kinematics.**
  * HINDA NaI detectors sit at θ_lab = 55, 90 and 125°. In-plane detectors are at φ = 0 and 180°, out-of-plane ones at φ = 270°.
  * The beam polarisation axis is horizontal (φ = 0/180°). So φ = 0/180° gives `CS_PAR` (E-vector in the scattering
    plane) and φ = 270° gives `CS_PERP`.
  * Circular beam: E = 81.3 MeV (6.5% FWHM). With an unpolarised target this is the unpolarised dσ/dΩ (`CS`); each
    detector is listed separately.
  * Linear beam: E = 83.4 MeV (2.7% FWHM), "nearly 100%" polarised. No polarisation correction or uncertainty is
    quoted.
  * The cross sections are lab dσ/dΩ (yield over lab effective solid angle). They are bin-centre corrected to the
    detector centre using a χEFT model.
* **Definition.** Σ3 = (σ∥ − σ⊥)/(σ∥ + σ⊥), where σ∥ is the in-plane cross section and σ⊥ the out-of-plane one.
  Σ3 values appear only in Fig. 3; they are derived from Table I and are **not** entered as separate rows.
* **Uncertainties.**
  * Point-to-point syst is 4.5–13.8%: ROI and fit window, TOF and shield cuts, and geometry for two out-of-plane detectors.
  * Normalisation is 2% (photon flux) ⊕ 1% (target thickness) = 2.24%. Absolute values are in `notes`.
  * The authors floated the circular and linear normalisations independently, so they are separate `dataset` tags.
* **Polarisabilities.** HBχPT fit (McGovern–Phillips–Grießhammer 2013) to the 16 points, in 10⁻⁴ fm³.
  * With BSR α+β = 14.0 ± 0.2: α = 13.8 ± 1.2 ± 0.1 ± 0.3 and β = 0.2 ∓ 1.2 ± 0.1 ∓ 0.3 (χ² = 14.7 for 13 dof).
  * Without BSR: α = 15.4 ± 1.8 and β = 2.1 ± 2.0 (χ² = 13.2 for 12 dof).

## 6. Search for later results (2022–Oct 2026): nothing new to add

Sources searched:
* INSPIRE: A2 and HIGS author/collaboration queries, citations of all five papers, and theses.
* arXiv API: nucl-ex, "Compton" plus "polarizabilities".
* HEPData: there are **no records** for any of the datasets above.
* Web search.

No further published proton RCS data below 400 MeV were found. Related items:
* Mornacchi, Pasquini and Pedroni, PRL 129, 102501 (2022), arXiv:2204.13491. This is a bootstrap DR **fit** (not data) to
  Σ2x, Σ2z, Σ3 (Sokhoyan, Mornacchi, HIGS, LEGS) and CS. It gives αE1 = 12.7 ± 0.8 ± 0.1, βM1 = 2.4 ± 0.6 ± 0.1,
  γE1E1 = −3.0 ± 0.6 ± 0.4 and γM1M1 = 3.7 ± 0.5 ± 0.1.
* Collicott et al. (A2) Σ3 in the Δ region (e.g. 286.9–307.1 MeV, with **positive** values in the same σ∥−σ⊥
  convention). This is preliminary, shown in figures only (Martel, EPJ Web Conf. 142 (2017) 01021; PoS CD2018), and
  cited as "in preparation" in Paudyal 2020. **It is unpublished.**
* HIGS ³He Compton (arXiv:2610.07579, 2026) is not a proton measurement. A2 2026 GDH (arXiv:2604.14385) is not Compton.

## Verification log (`verify_db_new.py`)

* [K] θ_cm from the boost agrees with θ_cm from the invariant t to within 0.0005° for all 126 rows.
* [M1] Martel: the EPS vector points match the CSV to ≤0.0014, which is the EPS resolution.
* [P] Paudyal: all 10 rows match the arXiv PDF text exactly.
* [Mo] Mornacchi CS: all 60 rows match the thesis beamtime-weighted mean to ≤0.0011. Σ3: 34 of 36 rows match; the 2
  exceptions are explained in section 4.
* [L] HIGS: all 16 rows match both the thesis and the arXiv PDF exactly.
