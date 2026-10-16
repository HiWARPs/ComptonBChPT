# Data provenance

Every data set of the fit database, with its origin and error model.  Nothing here is new: it is condensed from the
documentation of the curated tables (`code/data/db_old_README.md`, `code/data/db_new_README.md`), from the assembly
code (`code/fits/lib/database.py`), from the per-set table of the fit pipeline (`results/fits/dataset_table.md`) and
from the figure provenance of the Δ-region plots (`code/delta_plots/results/provenance.md`).  Those files hold the
full details (bin edges, cross-checks, which compilation used which point); consult them before using a data set.

## 1. Files in `data/`

| file | content |
|---|---|
| `fitdb.csv` | the fit database: 624 points (all data sets below, A2 2022 combined) with the provenance columns described in Sec. 2 |
| `fitdb_pb.csv` | the same with the A2 2022 data per beam time (March and July 2018, from the Mornacchi thesis) instead of the combined values: 720 points |
| `db_forward_obs.csv` | the forward points: dσ/dΩ(0°) and Σ2z(0°) at 12 energies, 50–150 MeV, derived from the empirical forward amplitudes (Sec. 3, set GHP) |
| `db_forward_amp.csv` | the forward amplitudes Re f and Re g themselves (μb GeV), with the uncertainties, from which the forward points are computed |

`fitdb.csv` and `fitdb_pb.csv` are built by `code/fits/scripts/s01_build_db.py` (`code/fits/lib/database.py`) from the
curated source tables in `code/data/`: `db_old.csv` (pre-2015 data, 477 rows), `db_new.csv` (A2 and HIγS 2014–2022,
126 rows), `db_klp2018.csv` (the 138-point database of the 2018 partial-wave analysis, exported from the author's
notebook), `db_mornacchi_thesis.csv` (A2 per beam time), `db_forward_obs.csv`.  `db_old_ext450.csv` is a plot-only
extension of the old data to 450 MeV.  The scripts that produced the curated tables from the raw sources are in
`code/data/scripts/`; the raw downloads themselves (HEPData records, arXiv sources, theses; `work/data/raw/`,
155 MB) and the text extractions of the theses are not redistributed.  Their locations are given per set below.

## 2. Columns of `fitdb.csv`

| column | meaning |
|---|---|
| `pid` | unique point id (`<dataset>_<nnn>`) |
| `dataset` | independent data set (own normalization) |
| `expt` | short label for legends |
| `source` | curated table of origin (`db_old`, `db_new`, `db_klp2018`, `db_mornacchi_thesis`, `db_forward_obs`) |
| `obs` | observable as in the source (`CS`, `S3`, `S2X`, `S2Z`, `CS_PAR`, `CS_PERP`) |
| `dobs` | observable code of the amplitude driver (`CS` lab, `CSCM` c.m., `S3`, `S2X`, `S2Z`, `SPAR`, `SPERP`) |
| `frame` | frame of the angle passed to the driver (`LAB` or `CM`) |
| `E` | lab photon energy in MeV (true energy where known) |
| `theta`, `theta_lab`, `theta_cm` | angle in `frame`, and in both frames (deg; exact two-body kinematics with M_p = 938.272 MeV) |
| `value` | measured value (nb/sr for cross sections; dimensionless for asymmetries) |
| `err_stat`, `err_syspt`, `err_pt` | statistical, point-to-point systematic, and their quadrature sum (used in χ²) |
| `norm` | 1σ (Gaussian-equivalent) scale uncertainty of the set, as a fraction; 0 = not floated |
| `flags` | `GMPF_out` (excluded by the pruning of Griesshammer et al. 2012 / McGovern et al. 2013), `KLP_E` (energy reconstructed from the 2018 bin label), `FWD` (forward point) |
| `prov` | provenance and error-model remark, as written by `database.py` |

Values are never converted between frames; only the units are changed (to nb/sr).

## 3. Data sets below the pion-production threshold (fits of the polarizabilities)

N is the number of points in the database (in parentheses: at E ≤ 150 MeV, as used in the fits).  "Legacy files" =
the author's own data files used in Lensky and Pascalutsa, Eur. Phys. J. C 65, 195 (2010) and in the 2018 PWA.

| dataset | experiment and reference | N | origin of the numbers | error model, normalization |
|---|---|---|---|---|
| `OlmosdeLeon2001` | TAPS at MAMI, Olmos de León et al., Eur. Phys. J. A 10, 207 (2001) | 65 (55; 54 in the main fit, see note a) | V. Olmos de León, PhD thesis, Mainz 2000, Tables A.1–A.5 (weighted means), transcribed by hand; identical to the legacy files | stat ⊕ 5% point-to-point (as in the earlier analyses); norm 3% |
| `Federspiel1991` | Illinois, Federspiel et al., Phys. Rev. Lett. 67, 1511 (1991) | 16 (16) | legacy files (the PRL has no table; the values are 3–5% below the published figure, cf. `db_old_README.md`) | stat; norm 2% |
| `Zieger1992` | MAMI (180°), Zieger et al., Phys. Lett. B 278, 34 (1992) | 2 (2) | legacy file | stat ⊕ per-point systematic (6.4%, 4.3%); no common scale |
| `Baranov1974_75_th90`, `_th150` | Moscow (Lebedev), Baranov et al., JETP Lett. 19, 398 (1974); Phys. Lett. B 52, 122 (1974); Yad. Fiz. 21, 689 (1975) | 3 + 4 (3 + 4) | HEPData ins102948, ins96876, checked against the JETP Lett. table | as quoted; norm 8% (LEGS 2001, Table IX) |
| `Hallin1993_TabII` | SAL, Hallin et al., Phys. Rev. C 48, 1497 (1993), Table II | 28 (13; the 55 points above 170 MeV form `Hallin1993_TabII_D`, Sec. 4) | HEPData ins369337 | c.m. angle and c.m. cross section (inferred); stat ⊕ (syst ⊖ 3.6% scale); norm 3.6% |
| `MacGibbon1995_tagged`, `_untagged` | SAL, MacGibbon et al., Phys. Rev. C 52, 2097 (1995) | 8 + 10 (8 + 10) | arXiv:nucl-ex/9507001 tables (= HEPData ins397101) | stat; norm 4% / 2.9% |
| `Pugh1957` | MIT, Pugh et al., Phys. Rev. 105, 982 (1957), as corrected by Frisch (1967) | 16 (16) | database of the 2018 PWA (`db_klp2018.csv`), total errors; energies reconstructed from the 90° legacy subset and the bin pattern (flag `KLP_E`) | total error; norm 10% adopted |
| `Hyman1959` | MIT, Hyman et al., Phys. Rev. Lett. 3, 93 (1959) | 12 (12) | as Pugh | total error; norm 10% adopted |
| `Goldansky1960` | Moscow, Goldansky et al., Nucl. Phys. 18, 473 (1960) | 5 (5) | database of the 2018 PWA | total error; norm 8% (MacGibbon et al. 1995) |
| `Oxley1958` | Chicago, Oxley, Phys. Rev. 110, 733 (1958) | 4 (4) | database of the 2018 PWA | total error; norm 10% adopted |
| `Bernardini1960` | Illinois, Bernardini et al., Nuovo Cim. 18, 1203 (1960) | 2 (2) | database of the 2018 PWA (c.m. → lab at the bin-label energy) | total error; norm 10% adopted |
| `A2_2022_CS`, `A2_2022_S3` | A2 at MAMI, Mornacchi et al., Phys. Rev. Lett. 128, 132503 (2022) | 60 + 36 (60 + 36) | Supplemental Material, Tables I and II (`raw/new/Mornacchi2022`) | stat ⊕ point-to-point syst; norm 3% (dσ); 5% uniform scale → Gaussian 2.9% (Σ3) |
| `A2_2018Mar_*`, `A2_2018Jul_*` (`fitdb_pb.csv` only) | the same data per beam time | 60 + 36 each | E. Mornacchi, PhD thesis, JGU Mainz 2021, Appendices A and B | as above, per beam time |
| `HIGS_2022_CIRC`, `HIGS_2022_LIN` | HIγS, X. Li et al., Phys. Rev. Lett. 128, 132502 (2022), Table I | 8 + 8 (8 + 8) | published table (`raw/new/Li2022`) | stat ⊕ point-to-point syst; norm 2.24% per set; unpolarized dσ at 81.3 MeV (circular beam), dσ∥ (5) and dσ⊥ (3) at 83.4 MeV (linear beam) |
| `GHP_fwd_CS`, `GHP_fwd_S2Z` | forward points from the empirical forward amplitudes: Re f, Gryniuk, Hagelstein, Pascalutsa, Phys. Rev. D 92, 074031 (2015); Re g, Phys. Rev. D 94, 034043 (2016) | 12 + 12 (11 + 11: the point at 150 MeV is not used) | table of arXiv:1604.00789; dσ(0) = \|f\|² + \|g\|², Σ2z(0) = −2fg/(f² + g²); Re f = mean of fits I and II with half the difference added to the error | errors of different energies strongly correlated: two common nuisance parameters n_f, n_g (Appendix D of the paper); only the rounding of the table is uncorrelated |

Note a: the TAPS point at 108.42 MeV, 133° (10.84 ± 0.94 nb/sr) is the only point flagged by the 3σ rule and is
excluded from the main database (`MAIN_DROP` in `code/fits/lib/common.py`); it is kept in `fitdb.csv`.

## 4. Data sets in the Δ region (comparison with the parameter-free prediction, Sec. IV D and Appendix G)

| dataset | experiment and reference | N (E ≤ 350 MeV) | origin of the numbers | error model, normalization |
|---|---|---|---|---|
| `A2_Martel2015_S2x` | A2, Martel et al., Phys. Rev. Lett. 114, 112501 (2015) | 4 | P. Martel, PhD thesis (UMass 2013), Table 7.3 "Summed" (= PRL Fig. 3) | stat ⊕ beam/target-polarization syst (the 2.1% scale is inside); not floated |
| `A2_Paudyal2020_S2z` | A2, Paudyal et al., Phys. Rev. C 102, 035205 (2020) | 10 | Table I (arXiv TeX source) | "Rand." ⊕ (Syst. ⊖ 10.4% scale); scale 10.4% floated |
| `LEGS2001_CS`, `LEGS2001_S3` | LEGS, Blanpied et al., Phys. Rev. C 64, 025203 (2001) | 77 + 79 | Tables II and III, transcribed from the manuscript posted by the LEGS group (Internet Archive copy) | stat ⊕ point syst; norm 2% (dσ), none (Σ) |
| `Hunger1997_75`, `Hunger1997_90` | MAMI, Hünger et al., Nucl. Phys. A 620, 385 (1997) (incl. Peise 1996 at 75°, Molinari 1996 at 90° c.m.) | 12 + 9 | HEPData ins458618 | stat; norm = quoted syst (4.4%/6%/8%; 4%/5%) |
| `Wissmann1999` | MAMI (CATS), Wissmann et al., Nucl. Phys. A 660, 232 (1999) | 8 | HEPData ins514218 | stat; norm 3%; lab frame (inferred) |
| `Camen2002` | MAMI (CATS+SENECA), Camen et al., Phys. Rev. C 65, 032202 (2002); numbers from Kossert et al., Eur. Phys. J. A 16, 259 (2003) | 7 | arXiv:nucl-ex/0210020 (= HEPData ins599960) | stat; norm 4.4%; lab frame (inferred) |
| `Hallin1993_TabII_D` | SAL, Hallin et al. (above 170 MeV) | 55 | HEPData ins369337 | as `Hallin1993_TabII` |
| `Gray1967` | Illinois, Gray and Hanson, Phys. Rev. 160, 1212 (1967) | 5 | HEPData ins944949 | as quoted; c.m. frame assumed; norm 5% adopted |
| `Genzel1972` | Bonn, Genzel et al., Lett. Nuovo Cim. 4, 695 (1972) | 14 | HEPData ins77720 | as quoted; c.m. frame; norm 5% adopted |
| `Baranov1966` | Moscow, Baranov et al., Sov. J. Nucl. Phys. 3, 791 (1966) | 10 | HEPData ins1392560 | as quoted; c.m. frame assumed; norm 20% |
| `Kolb2000` | SAL, Kolb et al., Phys. Rev. Lett. 85, 1388 (2000) | 1 | arXiv:nucl-ex/0003002 (free-proton calibration point) | stat ⊕ syst per point; not floated |
| `Nagashima1964` | Tokyo, Nagashima (INS report 81, 1964) | 1 | HEPData ins1101996 | as quoted; c.m. frame assumed; norm 5% adopted |

## 5. Data not available in numerical form (not in the database)

Wolf et al., Eur. Phys. J. A 12, 231 (2001) and Galler et al., Phys. Lett. B 503, 245 (2001) (LARA); the pilot beam
asymmetry of Sokhoyan et al., Eur. Phys. J. A 53, 14 (2017); DeWire et al., Phys. Rev. 124, 909 (1961); the
superseded LEGS 1996 release; the Baranov et al. 2001 compilation.  Details in `code/data/db_old_README.md`, Sec. 4,
and `code/data/db_new_README.md`, Sec. 2.
