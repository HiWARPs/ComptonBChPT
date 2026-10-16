# db_old.csv — pre-2015 world data on real Compton scattering off the proton (γp→γp), E_lab ≈ 30–400 MeV

Compiled 2026-10-13/14 for a chiral-EFT fit of the proton polarisabilities.
File: `work/data/db_old.csv` (477 rows: 397 dσ/dΩ + 80 Σ₃). Built reproducibly by
`work/data/scripts/build_db_old.py` from the raw sources in `work/data/raw/old/<source>/`.

## 1. Columns and conventions

```
dataset,texkey,obs,E_lab_MeV,dE_MeV,theta_lab_deg,theta_cm_deg,frame_of_cs,value,stat_err,sys_err_pt,norm_err_frac,units,source,notes
```

* `obs`: `CS` = unpolarised dσ/dΩ (nb/sr); `S3` = linear beam asymmetry Σ (dimensionless, `frame_of_cs = n/a`).
* `E_lab_MeV`: lab photon energy (bin centre as tabulated). `dE_MeV`: full bin width when it is stated (or, for TAPS, estimated; see the dataset notes). Otherwise it is blank.
* Angles: the angle as given in the source is copied verbatim. The other frame's angle is computed with exact two-body kinematics for a massless outgoing particle, using M_p = 938.272 MeV:
  β = E/(E+M_p), cosθ_cm = (cosθ_lab − β)/(1 − β cosθ_lab). The `notes` column starts with `orig_angle=lab` or `orig_angle=cm`. The computed angle is given to 3 decimals. The kinematics were checked against the author's legacy file `e_059.leon`: θ_lab = 59° at 58.94 MeV gives θ_cm = 61.95°, identical.
* `frame_of_cs`: the frame of the quoted dσ/dΩ. **Values are never converted between frames.** The only rescaling is a unit change: μb/sr ×1000 or 10⁻³² cm²/sr ×10 → nb/sr, done with exact decimal arithmetic and noted per row. Conversion if needed: dσ/dΩ_cm = dσ/dΩ_lab · (ω′_cm/ω′_lab)², where ω′_lab = E/[1+E(1−cosθ_lab)/M], ω′_cm = E·M/√(M²+2ME).
* Error columns (absolute, same units as `value`):
  * `stat_err`: the statistical error. When the source gives only a combined error, that combined error goes here and the `notes` say so.
  * `sys_err_pt`: a point-to-point systematic error, filled only when the source (or the paper as quoted by the compilations) gives one per point.
  * `norm_err_frac`: the correlated scale (normalisation) uncertainty as a fraction. `0` means the source explicitly says there is no common scale. Blank means none is stated.
  * Where a source's per-point "sys" already contains the common part (Hallin), the notes say not to add both.
* Energy cut: rows with E_lab > 400 MeV are dropped. The Δ peak in the lab is ≈ 320–340 MeV, and the user asked for ~350 MeV. The 400 MeV cut matches the counting convention of Griesshammer et al. 2012, Table 4.

π⁰ threshold: E_thr = m_π0 + m_π0²/(2M_p) = 144.69 MeV.

## 2. Summary table

| dataset | texkey | obs | N | E_lab range (MeV) | θ_lab | θ_cm | frame | norm | N(E<170) | N(E≤E_thr) | N(E≤150) | N(E≤350) | provenance class |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| OlmosdeLeon2001 | deLeon:2001dnx | CS | 65 | 58.94–163.71 | 59–155 | 62.0–158.4 | lab | 0.03 | 65 | 50 | 55 | 65 | primary-equivalent (thesis table = author files) |
| Federspiel1991 | Federspiel:1991yd | CS | 16 | 33.5–70 | 60,135 | 61.7–137.7 | lab | 0.020 | 16 | 16 | 16 | 16 | **secondary (author legacy)** |
| Zieger1992 | Zieger:1992jq | CS | 2 | 98,132 | 180 | 180 | lab | 0 | 2 | 2 | 2 | 2 | **secondary (author legacy)** |
| Baranov1974_75_th90 | Baranov:1974ta; Baranov:1974ec; Baranov:1975ju | CS | 3 | 80.9–109.9 | 90 | 94.6–96.0 | lab | 0.08 | 3 | 3 | 3 | 3 | primary (HEPData + JETP Lett.) |
| Baranov1974_75_th150 | same | CS | 4 | 81.9–111.1 | 150 | 152.2–152.9 | lab | 0.08 | 4 | 4 | 4 | 4 | primary |
| Pugh1957_Frisch1967_legacy_partial | Pugh:1957zz | CS | 6 | 55–135 | 90 | 93.2–97.2 | lab | – | 6 | 6 | 6 | 6 | **secondary, partial** |
| Hyman1959_legacy_partial | Hyman:1959zz | CS | 8 | 60–128 | 90 | 93.4–96.9 | lab | – | 8 | 8 | 8 | 8 | **secondary, partial** |
| Hallin1993_TabII | Hallin:1993ft | CS | 83 | 130–291 | 24.6–136.4 | 31–141 | cm (inferred) | 0.036 | 27 | 8 | 13 | 83 | primary (HEPData) |
| Hallin1993_TabI | Hallin:1993ft | CS | 21 | 149–286 | 24.7–135.2 | 30–144 | cm (inferred) | 0.036 | 4 | 0 | 4 | 21 | primary, **not independent of TabII** |
| MacGibbon1995_tagged | MacGibbon:1995in | CS | 8 | 73.2–96.8 | 90,135 | 94.2–138.7 | lab | 0.04 | 8 | 8 | 8 | 8 | primary (arXiv TeX = HEPData) |
| MacGibbon1995_untagged | MacGibbon:1995in | CS | 10 | 105–145 | 90,135 | 95.8–140.2 | lab | 0.029 | 10 | 8 | 10 | 10 | primary |
| Blanpied2001_LEGS | Blanpied:2001ae | CS | 77 | 213.1–333.8 | 51.9–126.9 | 65–135 | cm | 0.02 | 0 | 0 | 0 | 77 | primary (LEGS-posted manuscript) |
| Blanpied2001_LEGS | Blanpied:2001ae | S3 | 79 | 213.1–333.8 | 51.9–126.9 | 65–135 | – | 0 | 0 | 0 | 0 | 79 | primary |
| Hunger1997_th75cm | Hunger:1997pt; Peise:1996zz | CS | 17 | 207.3–399.1 | 58.9–65.2 | 75 | cm | 0.044/0.06/0.08 | 0 | 0 | 0 | 12 | primary (HEPData) |
| Hunger1997_th90cm | Hunger:1997pt; Molinari:1996zw | CS | 13 | 250.4–396.3 | 72.7–77.8 | 90 | cm | 0.04/0.05 | 0 | 0 | 0 | 9 | primary (HEPData) |
| Wissmann1999 | Wissmann:1999vi | CS | 11 | 199.1–390.8 | 130.7 | 137.9–142.6 | lab (inferred) | 0.03 | 0 | 0 | 0 | 8 | primary (HEPData) |
| Camen2002_Kossert2003 | Camen:2001st; Kossert:2002ws | CS | 10 | 211.1–388.5 | 136.2 | 143.1–146.9 | lab (inferred) | 0.044 | 0 | 0 | 0 | 7 | primary (arXiv TeX = HEPData) |
| Kolb2000 | Kolb:2000ix | CS | 1 | 247 | 135 | 142.9 | lab | – (sys in col.) | 0 | 0 | 0 | 1 | primary (arXiv TeX) |
| Gray1967 | Gray:1967zzd | CS | 5 | 210–310 | 75.6–135 | 90–142 | cm (assumed) | – | 0 | 0 | 0 | 5 | primary (HEPData) |
| Genzel1972 | Genzel:1972ba; Genzel:1976wy | CS | 19 | 236.9–386.5 | 38.5–120 | 50–130 | cm (assumed) | – (sys in col.) | 0 | 0 | 0 | 14 | primary (HEPData) |
| Baranov1966 | Baranov:1966aaa | CS | 10 | 212.9–249.2 | 59.1–141.1 | 68.8–147.7 | cm (assumed) | 0.20 | 0 | 0 | 0 | 10 | primary (HEPData) |
| Baranov1966_refined | Baranov:1966oqa | CS | 6 | 235–248.3 | 46.4–140.4 | 55.4–147.5 | cm (assumed) | – | 0 | 0 | 0 | 6 | primary; **4 pts overlap Baranov1966 T4** |
| Nagashima1964 | Nagashima:1964zz | CS | 2 | 314, 351 | 74.2–75.8 | 90.0–90.3 | cm (assumed) | 0.05 | 0 | 0 | 0 | 1 | primary (HEPData) |
| Barbiellini1968 | Barbiellini:1968rlx | S3 | 1 | 317.5 (300–335) | 75.4 | 90 | – | – | 0 | 0 | 0 | 1 | primary (HEPData); **sign convention ≠ LEGS** |
| **TOTAL** | | | **477** | | | | | | **153** | **113** | **129** | **456** | |

Low-energy (E < 170 MeV) points: there are 153 in total.
* 117 are primary: OdL 65 + Hallin TabII 27 + MacGibbon 18 + Baranov 7.
* 18 are secondary: Federspiel 16 + Zieger 2.
* 14 are partial legacy points: Pugh/Frisch 6 + Hyman 8.
* 4 are Hallin Table I points.

Pasquini et al. 2019 ("FULL") have 150 points at ≤150 MeV. This file has 129 points at ≤150 MeV. The difference comes from the old sets that are missing or only partial: Oxley 4, Hyman 12, Goldansky 5, Bernardini 2, Pugh 16.

## 3. Datasets — provenance, errors, compilation usage

### Olmos de León et al. 2001 (TAPS/MAMI), EPJA 10, 207 — `deLeon:2001dnx`
* **Source:** the EPJA article is not open access, and the tables were not on arXiv or HEPData. Values are taken from V. Olmos de León, PhD thesis, Mainz 2000, Anhang A, Tabellen A.1–A.5, column "gewicht. Mittel" (weighted mean of runs 2 and 3), pp. 123–125. The thesis is at http://download.uni-mainz.de/fb08-kpha2/publications/phd/2000_olmos.ps.pdf (saved in `raw/old/OlmosdeLeon2000_thesis/`).
* **Transcription:** made by hand from 220-dpi page renders, in `transcription_weighted_mean.csv`.
* **Cross-check:** all 65 values and errors are **identical** to the author's independent legacy files `a059/a085/a133.leon01` and `Olmos_{59,85,107,135,155}.dat` (104 comparisons, 0 mismatches). This, and the matching counts (65 points; 55 at ≤150 MeV as in Pasquini 2019), indicate that the published table equals the thesis weighted means. This could not be checked against the EPJA PDF itself.
* **Coverage:** 13 energies (58.94–163.71 MeV, 20 tagger channels = 6–10 MeV per bin) × 5 lab angles (59, 85, 107, 133, 155°, ±10° bins).
* **Frame:** lab.
* **Errors:**
  * stat: the thesis statistical error. It already includes a ~3% cut-variation component.
  * `sys_err_pt`: 5% point-to-point (geometry/simulation) per the published paper, as quoted by Griesshammer 2012, McGovern 2013 ("following Wissmann") and Pasquini 2019. These works add it in quadrature to the statistical error.
  * norm: 3% (flux 2% ⊕ target 2%). The thesis itself quotes a 6% total systematic (5% from the effective solid angle).
* **Compilations:** all three main analyses use it.
  * Griesshammer 2012 and McGovern 2013 exclude the outlier at **133°, 108.42 MeV** (10.84 nb/sr).
  * The SELECTED set of Pasquini 2019 drops the same point.
  * The ≤150 MeV analyses use 55 points.

### Federspiel et al. 1991 (Illinois), PRL 67, 1511 — `Federspiel:1991yd` — SECONDARY
* **The PRL contains no table:** the cross sections appear only in Fig. 2. This was checked in the KEK scan of preprint ILL-NPL-91-012 (https://lib-extopc.kek.jp/preprints/PDF/1991/9109/9109187.pdf, saved in `raw/old/Federspiel1991_KEKscan/`).
* **Values used:** the 16 values are the author's legacy files `DATA2/a060.fed91` and `a135.fed91`, embedded in `LP_ms2/figs/dasfixA.agr` and other figures.
* **Difference from the published figure:** these values are systematically **3–5% lower** than the published figure. The digitised-looking version is `Codes/execs/CSeft/RCSdata/feder_91_th*.exp`, which has a constant error of 1.44 nb/sr; the db/exp ratio is 0.95–0.97.
* **Likely origin:** this is consistent with the later corrections "A.M. Nathan, priv. comm. 1993" used in the Baranov et al. 2001 compilation. The origin of the legacy numbers is not documented, so **verify before use.**
* **Errors:** stat. norm 2.0% (preprint); Griesshammer 2012 quotes 2.2%, and LEGS 2001 used 2.0%. Bins are 4 MeV wide.
* **Compilations:** used by all. Outlier **135°, 44.2 MeV** is excluded by Griesshammer 2012, McGovern 2013 and the Pasquini SELECTED set.

### Zieger et al. 1992 (MAMI, 180°), PLB 278, 34 — `Zieger:1992jq` — SECONDARY
* **Source:** HEPData ins338241 has only the polarisabilities. The values (98 MeV: 18.42±3.39; 132 MeV: 21.82±1.18 nb/sr) come from the author's file `RCSdata/zieg_92.exp`.
* **Consistency check:** the relative statistical errors (18.4%, 5.4%) and the systematics (6.4%, 4.3%) are exactly those quoted in Griesshammer 2012. `sys_err_pt` = 6.4% and 4.3%. norm = 0, since LEGS 2001 Table IX says the two points have no common scale. Bins are 16 MeV.
* **Frame:** taken as lab, following the author's own plot (`RCSdata/rcs.ps`, "dσ/dΩ^lab"). This was not verified in the paper. At 180°, c.m./lab ≈ 1.17–1.28.

### Baranov et al. 1974/75 (Moscow) — `Baranov:1974ta` (JETP Lett. 19, 398), `Baranov:1974ec` (PLB 52, 122), `Baranov:1975ju` (Yad.Fiz. 21, 689)
* **Source:** HEPData ins102948 and ins96876 (https://www.hepdata.net/record/ins102948). The values were checked digit by digit against the table in the JETP Lett. paper (INSPIRE fulltext, saved in `raw/old/Baranov1974_JETPLett/`): all 7 points agree.
* **Frame:** lab, θ_lab = 90°, 150°.
* **Errors:** the quoted errors include random and reduction errors plus ~1% absolutization systematic, so they are placed in stat_err. norm = 0.08 is the scale uncertainty LEGS 2001 (Table IX) used for "Moscow 81–110 MeV". The original quotes only ~1%. Baranov 2001 assigns different systematic errors to the 90° and 150° sets, but those values could not be accessed.
* **Compilations:** treated as two sets (Baranov 2001, Pasquini 2019 sets 6 and 7). **The 150° set is excluded** by Griesshammer 2012 and McGovern 2013 (χ²/dof > 3). MacGibbon 1995 excluded the whole Moscow set from its global fit.
* **Legacy check:** `a090.bar74` is identical.

### Hallin et al. 1993 (SAL), PRC 48, 1497 — `Hallin:1993ft`
* **Source:** HEPData ins369337 (https://www.hepdata.net/record/ins369337).
  * Tables 5–16 are the paper's Table II (excitation functions, 5 MeV bins), 83 points, in `Hallin1993_TabII`.
  * Tables 1–4 are the paper's Table I (angular distributions at 149/182/230/286 MeV, 15 MeV bins), 21 points, in `Hallin1993_TabI`.
* **Angles:** these are **c.m.** angles. The abstract's 25° < θ_lab < 135° maps to 31–144° c.m., and Griesshammer quote "141° (cm)".
* **Frame (inferred c.m.):** the author's legacy `a060.hal93` equals these values multiplied by the c.m.→lab Jacobian at θ_lab ≈ 60°. Magnitudes vs TAPS and LEGS also favour c.m. Verify in the paper.
* **Errors:** stat as tabulated. `sys_err_pt` = per-point "sys" as tabulated. The paper's systematics are a correlated 3.6% (flux 3% ⊕ target 2%) plus 1–2% point-to-point, so the column contains the norm part; do not add both. norm 0.036. LEGS 2001 used 3.8%.
* **Compilations:**
  * Pasquini 2019 and the Pedroni 2020 test discard Table I as possibly not independent. They use the 13 Table-II points with E ≤ 150 MeV.
  * McGovern 2013 uses Table II only at or below the cusp (≤150 MeV).
  * Griesshammer 2012 excludes all Hallin data from the final fit (27 points < 170 MeV, χ² = 41/27). They flag **141° cm, 170 MeV** (56.0±2.4) as an outlier.
  * The Griesshammer review lists "77 points" for SAL 93. HEPData Table II has 83 (incl. duplicate energies at 141°), and the difference could not be resolved.
* **Legacy check:** values identical. The author's errors are stat⊕sys.
  * One value differs: `a141.sal93` has 222 MeV = 72.0±4.0, while HEPData has 71.7±3.1(stat)±2.7(sys) (and a second point 73.0±3.4).
  * The author's angular files `e149/e182/e230/e286.sal93` exist in two variants: one with stat⊕sys, one with larger errors of unknown origin.
  * `e135.sal93` is not a published angular distribution.

### MacGibbon et al. 1995 (SAL), PRC 52, 2097 — `MacGibbon:1995in`
* **Source:** arXiv:nucl-ex/9507001 TeX source, Tables "Tagged…" and "Untagged…". These are identical to HEPData ins397101.
* **Frame:** lab, θ_lab = 90°, 135°.
* **Tagged set (73.2–96.8 MeV, 8 MeV bins):** stat errors. norm 3.0–4.0% per the paper; 0.04 is entered, which LEGS also used.
* **Untagged set (100–150 MeV, 10 MeV bins):** the error is the diagonal of a 5×5 error matrix, and adjacent bins are strongly anti-correlated. The full matrix is not published. norm 2.9%.
* **Compilations:** used as two separate sets in McGovern 2013 and Pasquini 2019 (sets 11 and 12). Only ±0.4% (tagged) and ±1.0–1.4% (untagged) of the systematics are uncorrelated (Griesshammer 2012).
* **Legacy check:** `RCSdata/sal_95_*.exp` are identical. `a090.mac95` is rounded to 0.1.

### Blanpied et al. 2001 (LEGS), PRC 64, 025203 — `Blanpied:2001ae`
* **Source:** the APS PDF is paywalled and the paper is not on arXiv or HEPData. Tables II (dσ/dΩ_cm) and III (Σ) were transcribed from the LEGS group's own posted manuscript `legs_n_delta.pdf` (74 pp., PDF dated 2001-04-18), saved via the Internet Archive: http://web.archive.org/web/20050118065130/http://www.legs.bnl.gov:80/pub/pi_comp/legs_n_delta.pdf (in `raw/old/Blanpied2001_LEGSweb/`).
* **Verification:** the abstract numbers equal the published abstract. The transcription was checked against the pdftotext token stream (24/24 rows identical).
* **Legacy check:** `a90.legs01` (12/12), `sigma286.legs01` (7/7) and `asym90.legs01` match.
  * `asym90.legs01` has one typo: 254.1 MeV error 0.04, where the table has 0.047.
  * `e_286.legs` has the same values but larger errors (e.g. 10.7 vs 7.3), probably an older version.
* **Coverage:** 77 dσ/dΩ + 79 Σ points at θ_cm = 65–135°, E = 213.1–333.8 MeV.
* **Errors:** the error is stat ⊕ point-specific systematic. The cross sections have a 2% common scale. Σ has no additional scale uncertainty.
* **Compilations:** Griesshammer 2012 and McGovern 2013 fit the Δ parameters to MAMI rather than LEGS, because of the LEGS–MAMI disagreement at θ_cm > 100° (up to ~25% at 115°, the Δ peak). LEGS itself fitted all world data with floating scales (Table IX). LEGS warns that its 1996 PRL 90° data (L8 release) must **not** be combined with these data (double counting). The 1996 release was saved (`raw/old/Blanpied1996_LEGSweb/data_release_l8_90_1.html`) but **excluded** from the db. The author's `asym90.legs96` corresponds to that superseded set.

### Hünger et al. 1997 (MAMI), NPA 620, 385 — `Hunger:1997pt` (incl. Peise 1996 `Peise:1996zz` at 75° and Molinari 1996 `Molinari:1996zw` at 90°)
* **Source:** HEPData ins458618, Tables 1 (θ_cm = 75°) and 2 (θ_cm = 90°). The c.m. frame is explicit (RF=CM).
* **Errors:** stat. The quoted systematic (75°: 4.4%/6%/8% below 350 / 350–380 / above 380 MeV; 90°: 4%/5% below/above 360 MeV) is entered as norm. It contains the common flux 2% and target 2%. LEGS used 4.4% and 5.0%.
* **Coverage:** 17 + 13 points ≤400 MeV, matching the Griesshammer counts.
* **Relation to Peise:** the Peise PhD thesis (Mainz 1996, saved) Table 6 at 75° differs by −1.6% to +3% (older analysis, 5.3% systematic). Hünger is the final version.
* **Legacy mislabel:** the author's `a90.lags97` (legend "LEGS 97" in `FIGS/*.mgr`) is **in fact these MAMI 90° c.m. data** (13/13 identical), not LEGS.

### Wissmann et al. 1999 (MAMI, CATS), NPA 660, 232 — `Wissmann:1999vi`
* **Source:** HEPData ins514218, θ_lab = 130.7°.
* **Frame:** **not stated, taken as lab.** The peak value of ~170 nb/sr agrees with the lab-frame Kossert/Camen data at 136° lab, whereas c.m. would be ×1.47. Verify.
* **Errors:** stat. norm 3% (flux, target, acceptance).
* **Compilations:** Griesshammer 2012 notes two points beyond the peak lie high.

### Camen et al. 2002 (MAMI, CATS+SENECA), PRC 65, 032202 — `Camen:2001st`; numbers from Kossert et al. 2003 — `Kossert:2002ws`
* **Source:** Camen shows the data only graphically, shifted to θ_cm = 135°. The free-proton cross sections of the same experiment are tabulated in Kossert et al., EPJA 16, 259 (2003), appendix (arXiv:nucl-ex/0210020; = HEPData ins599960 Table 1). θ_lab = 136.2°; 14 points, of which 10 are ≤400 MeV (= Griesshammer count).
* **Frame inferred lab:** Kossert Fig. 7. Camen Fig. 2 (c.m., 135°) values are 1.2–1.7× larger, rising with E as the lab→cm Jacobian does.
* **Errors:** stat. Systematic 4.4% as norm.

### Kolb et al. 2000 (SAL), PRL 85, 1388 — `Kolb:2000ix`
* **Source:** the free-proton calibration point, 94.6 ± 9.0 (stat) ± 3.8 (syst) nb/sr at 247 MeV, θ_lab = 135° (arXiv:nucl-ex/0003002).
* **Notes:** not in the Griesshammer table; it is added for completeness. Frame lab.

### Older Δ-region data (all from HEPData; frames not stated, c.m. assumed unless noted)
* **Gray & Hanson 1967**, PR 160, 1212 — `Gray:1967zzd`. HEPData ins944949 Table 3. 5 points, 210–310 MeV, with E bins. This is "Illinois 67" in Griesshammer (5 points).
* **Genzel et al. 1972**, Lett. Nuovo Cim. 4, 695 — `Genzel:1972ba`.
  * Source and coverage: HEPData ins77720, θ_cm = 50–130° (RF=CM). 19 points ≤400 MeV, the same count as Griesshammer's "Bonn 76".
  * Relation to the 1976 paper: Griesshammer cite Genzel:1976wy (Z.Phys.A 279, 399), which is not on HEPData and could not be checked.
  * Errors: the per-point % systematic (4.4–6.3%) is put in `sys_err_pt`.
* **Baranov et al. 1966**, Sov.J.Nucl.Phys. 3, 791 — `Baranov:1966aaa`. HEPData ins1392560 Tables 3 (215 MeV) and 4 (248 MeV). 20% common photon-spectrum uncertainty (norm 0.20). Griesshammer: the 215 MeV points at 91.9° and 107.2° (cm) are implausibly low.
* **Baranov et al. 1966 "refinement"**, JETP 23, 242 — `Baranov:1966oqa`. HEPData ins1387582. This is a re-analysis of Baranov 1962 (`Baranov:1961sya`, 247 MeV). The 4 backward points overlap Table 4 of `Baranov:1966aaa` (which "incorporates the revised 1962 data"), so **use only one of them.** Griesshammer's "Moscow 66" count of 12 = T3 (5) + T4 (5) + the 2 forward refined points.
* **Nagashima 1964** (Tokyo INS report 81, via the Bonn compilation) — `Nagashima:1964zz`. HEPData ins1101996. 2 points ≤400 MeV at ~90° c.m. 5% systematic.
* **Barbiellini et al. 1968** (Frascati), PR 174, 1665 — `Barbiellini:1968rlx`. HEPData ins55887. Photon asymmetry −0.35 ± 0.10 at 90° c.m., E = 300–335 MeV. **Its sign convention is opposite to LEGS:** LEGS Σ(90°, ~318 MeV) ≈ +0.38. The value is copied as given.

### Secondary partial legacy sets (only the 90° subsets survive in the author's files)
* **Hyman et al. 1959**, PRL 3, 93 — `Hyman:1959zz`. `DATA2/a090.hym59`: 8 points, 60–128 MeV.
* **Pugh et al. 1957**, PR 105, 982 — `Pugh:1957zz`, as corrected by D. Frisch 1967 (cited as "Fri67" in Baranov et al. 2001; the "MIT 67" set in Griesshammer 2012). `DATA2/a090.fri67`: 6 points.
* For both: the error type and normalisation are unknown, and the values could not be verified. The full sets have 12 and 16 points.

## 4. Requested datasets that could NOT be obtained numerically from public sources (not invented)

| dataset | texkey | status |
|---|---|---|
| Oxley 1958, PR 110, 733 | `Oxley:1958zz` | APS paywall; no HEPData. Griesshammer: 4 pts, 25–87 MeV, 70–150°. Excluded from their fits (tiny statistical errors). |
| Hyman 1959 (full set) | `Hyman:1959zz` | Only the 90° subset is in the author files (see above). Full set: 12 pts, 50° & 90°. |
| Goldansky et al. 1960, Nucl. Phys. 18, 473 | no INSPIRE record found (cited as "Goldansky:1960" / "GOLDANSKY1960473") | Elsevier paywall; no HEPData. 5 pts, 40–70 MeV, 45–150°. ~8% normalisation (MacGibbon 1995). |
| Bernardini et al. 1960, Nuovo Cim. 18, 1203 | `Bernardini:1960wya` | The HEPData record ins1184998 contains **only π⁰ photoproduction**. Griesshammer: 17 pts at 100–290 MeV (7 at 120–200 MeV per their text), not fitted. Pasquini: 2 pts ≈ 135°, 140 MeV. |
| Pugh 1957 / Frisch 1967 (full set) | `Pugh:1957zz` | Only the 90° subset is in the author files. Full set: 16 pts, 45/90/135°, 50–130 MeV. |
| DeWire et al. 1961, PR 124, 909 | `DeWire:1961zz` | APS paywall. Griesshammer: 5 pts, 275–425 MeV. |
| Galler et al. 2001, PLB 503, 245 / Wolf et al. 2001, EPJA 12, 231 (LARA) | `Galler:2001ht`, `Wolf:2001ha` | arXiv sources (nucl-ex/0102003, nucl-ex/0109013) contain **figures only** (>700 points, 294 ≤400 MeV in Griesshammer). Not on HEPData. The numbers were obtained privately by the compilers. |
| Camen 2002 (as published) | `Camen:2001st` | Only a figure, but the same data are tabulated in Kossert 2003 (included above). |
| Peise 1996 / Molinari 1996 | `Peise:1996zz`, `Molinari:1996zw` | Superseded/summarised by Hünger 1997 (included). |
| Fuhrberg 1995, NPA 591, 1 | `Fuhrberg:1995zz` | This paper is **Compton scattering on ⁴He at 87 MeV**, not the proton. No proton Compton dataset by Fuhrberg was found. |
| Blanpied 1996, PRL 76, 1023 / Blanpied 1997, PRL 79, 4337 / Tonnison 1998 | `LEGS:1996cao`, `Blanpied:1997zz` | Superseded by Blanpied 2001; LEGS says not to combine (double counting). The 1996 90° release (12 dσ/dΩ + 12 Σ) is saved in raw but excluded. |
| Baranov et al. 2001 compilation, Phys.Part.Nucl. 32, 376 | `Baranov:2001jv` | Holds the corrected old data with systematic assignments used by Griesshammer and Pasquini. Not accessible. The arXiv talk nucl-ex/0011015 has only fit results. |

With institutional access, these DOIs would complete the low-energy set: PR 105, 982; PR 110, 733; PRL 3, 93; Nucl. Phys. 18, 473; Nuovo Cim. 18, 1203; PR 124, 909. The Baranov 2001 PEPAN review and the published EPJA 10, 207 / PLB 278, 34 / PRC 48, 1497 tables would allow the inferred frames to be verified.

## 5. Compilations consulted (sources in `raw/old/arXiv_*`)
* **Griesshammer, McGovern, Phillips, Feldman, PPNP 67 (2012) 841, arXiv:1203.6834.**
  * Table 4 (`experiment.tex`) gives sets, ranges and N ≤400 MeV. The arXiv source contains **no data tables**; the "supplemental material online" contains plots only.
  * Fit set: E ≤ 170 MeV. Excludes Oxley and Bernardini, Baranov 150°, all of Hallin (2012 fit), and the outliers Federspiel 135°/44 MeV, OdL 133°/108 MeV and Hallin 141°cm/170 MeV.
  * Adds 5% point-to-point to OdL and floats the normalisation of each set.
* **McGovern, Phillips, Griesshammer, EPJA 49 (2013) 12, arXiv:1210.4104.** Same selection, but **re-includes Hallin ≤150 MeV**. MacGibbon is counted as two sets. Cut-off 170 MeV.
* **Pasquini, Pedroni, Sconfietti, PRC 98 (2018) 015204 (arXiv:1711.07401) and J.Phys.G 46 (2019) 104001 (arXiv:1903.07952), Appendix A Table.** 13 sets at ≤150 MeV.
  * FULL = 150 points; TAPS = 55.
  * SELECTED = 137. In their wording it is "FULL except Oxley, Bernardini, Baranov:1975ju, OdL 133°/108 MeV and Federspiel 135°/44 MeV". The count of 13 removed points does not match removing all 7 Baranov points, so their exact Baranov cut is unclear.
  * Uses the Baranov 2001 data selection for Oxley, Hyman, Goldansky and Pugh. Discards Hallin Table I.
  * No supplementary numeric tables.
* **Pedroni & Sconfietti, J.Phys.G 47 (2020) 054001 (arXiv:1909.03885).** The same 13 sets / 150 points. Its Table "normal" gives the fitted offsets f_k − 1 per set, which is useful for checking normalisation behaviour.
* **LEGS 2001, Table IX:** the scale uncertainties used for world data (Mainz 75° 4.4%, Mainz 90° 5.0%, SAL93 3.8%, SAL95 4.0%, MPI 4.3%/6.4%, Moscow 8.0%, Illinois 2.0%).

## 6. Cross-check against the author's legacy (xmgrace) data

* **How the dump was made:** all data sets embedded in `LP_ms2/figs/*.agr`, `Codes/execs/CSeft/{FIGS/*.mgr, NewFigs/*.agr, *.agr}` and `RCSdata/*.exp` were extracted by `scripts/dump_agr_data.py` (with `parse_agr2.py`). The result is `raw/old/author_legacy/author_agr_embedded_data.txt` (56 unique sets).
* **How the comparison was run:** with `scripts/crosscheck_agr.py`. The full output is in `raw/old/author_legacy/crosscheck_agr_output.txt`.
* **Findings:**
  * **OdL:** 65/65 identical (values and errors).
  * **OdL fixed-energy files `e_XXX.leon`:** these are the same data converted to **c.m.** (c.m. angle and c.m. cross section). Some variants carry inflated errors of unclear origin.
  * **OdL, `Codes/kmat/DATA/e_164.leon` (in `FIGS/4report.mgr`):** this variant mixes c.m. angles with lab cross sections.
  * **LEGS:**
    * `a90.legs01` and `sigma286.legs01` are identical.
    * `asym90.legs01` is identical except one error (254.1 MeV: 0.04 vs 0.047).
    * `e_286.legs` has the same values but larger errors.
    * `asym90.legs96` is the superseded 1996 data.
  * **"LEGS 97" legend in `FIGS/*.mgr` (`a90.lags97`):** in fact the MAMI Hünger/Molinari 90° c.m. data (13/13 identical ≤400 MeV).
  * **Hallin:** values are identical; the legacy errors are stat⊕sys. One value differs (141° cm, 222 MeV: legacy 72.0 vs HEPData 71.7/73.0). `e135.sal93` is not a published distribution. `a060.hal93` = HEPData θ_cm = 67° converted to lab (confirms c.m. frame of HEPData values).
  * **MacGibbon:** identical to `sal_95_*.exp`.
  * **Baranov 74:** identical.
  * **Zieger:** identical.
  * **Federspiel:** the `a0xx.fed91` values (used in the db) are 3–5% below `RCSdata/feder_91_*.exp`, which appear digitised from the PRL figure (constant ±1.44). Energies differ by ≤0.5 MeV.

## 7. Reproduction
```
cd work/data/scripts
python3 -I build_db_old.py ../raw/old ../db_old.csv
python3 -I verify_sources.py ../db_old.csv ../raw/old
python3 -I crosscheck_agr.py ../db_old.csv ../raw/old/author_legacy/author_agr_embedded_data.txt ../raw/old/author_legacy
python3 -I summarize_db.py ../db_old.csv
```
* `verify_sources.py` re-reads every HEPData-derived row from the YAML files, and the MacGibbon and Kossert rows from the alternate source (HEPData vs TeX). It finds 0 mismatches.
* The two hand transcriptions (OdL thesis, LEGS Tables II/III) are verified by `check_odl.py`, `check_legs.py` and the legacy comparison.
