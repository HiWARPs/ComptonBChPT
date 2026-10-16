# ComptonBChPT — NNLO BχPT calculations and partial-wave analysis of proton Compton scattering

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23413172.svg)](https://doi.org/10.5281/zenodo.23413172)

Code, data and amplitude tables for the next-to-next-to-leading-order (NNLO) calculations of proton Compton
scattering in covariant baryon chiral perturbation theory (BχPT) with the Δ(1232), in the framework of

> V. Lensky and V. Pascalutsa, *Predictive powers of chiral perturbation theory in Compton scattering off protons*,
> Eur. Phys. J. C 65, 195 (2010), [doi:10.1140/epjc/s10052-009-1183-z](https://doi.org/10.1140/epjc/s10052-009-1183-z);
> V. Lensky, J. McGovern and V. Pascalutsa, *Predictions of covariant chiral perturbation theory for nucleon
> polarisabilities and polarised Compton scattering*, Eur. Phys. J. C 75, 604 (2015),
> [doi:10.1140/epjc/s10052-015-3791-0](https://doi.org/10.1140/epjc/s10052-015-3791-0),

and for the partial-wave analysis (PWA) of the proton Compton-scattering data below the pion-production threshold of

> N. Krupina, V. Lensky and V. Pascalutsa, *Partial-wave analysis of proton Compton scattering data below the
> pion-production threshold*, Phys. Lett. B 782, 34 (2018),
> [doi:10.1016/j.physletb.2018.04.066](https://doi.org/10.1016/j.physletb.2018.04.066).

The fits of the polarizabilities to the world database, the error budget, the radiative corrections and the figures
are those of

> V. Pascalutsa, *Chiral-EFT prediction for proton polarizabilities tested by the new Compton-scattering data*,
> [journal reference / arXiv identifier: to be added] (2026),

whose complete analysis can be rerun from the tabulated amplitudes in a few minutes (see below).

Author: Vladimir Pascalutsa, Institut für Kernphysik, Johannes Gutenberg-Universität Mainz.  Version 1.1
(2026-10-16).  Repository: https://github.com/vladipas/ComptonBChPT; archived on Zenodo, DOI [10.5281/zenodo.23413172](https://doi.org/10.5281/zenodo.23413172) (all versions; v1.0.1: 10.5281/zenodo.23413173).  Citation metadata: `CITATION.cff`.  License: MIT for the code (`LICENSE`), CC BY 4.0 for the data, amplitude tables and results (`LICENSE-DATA`).

The package contains (i) the world database of proton Compton-scattering data used in the paper, with the
provenance of every point, (ii) the tabulated NNLO BχPT amplitudes for every data point and for the curves of the
figures, as exact quadratic forms in the polarizabilities, (iii) the Python code of the fits, of the error budget,
of the partial-wave analysis (PWA), of the radiative corrections, of the total cross section and of the figures,
(iv) the Fortran code of the NNLO amplitude, and (v) the results.  With (ii) the complete fit analysis of the 2026 paper
can be rerun in a few minutes without compiling the Fortran code.

## Layout

```
README.md, CITATION.cff, LICENSE, LICENSE-DATA, requirements.txt
data/           fitdb.csv, fitdb_pb.csv (the fit database with provenance columns), db_forward_obs.csv,
                db_forward_amp.csv (the forward points and amplitudes), DATA_PROVENANCE.md
amplitudes/     the tabulated amplitudes (README.md there describes the format and the configurations)
code/           the analysis code; this directory is the author's former `work/` directory, so that all relative
                paths of the scripts are unchanged:
  code_v2/        Fortran code of the NNLO BχPT Compton amplitude with its fit driver (driver/README.md),
                  regression tests and validation tools; GIT_COMMIT.txt gives the version
  fits/           BχPT fits: lib/ (database assembly, driver interface, forward-point error model, truncation
                  uncertainty, χ² machinery), scripts/s00–s10 (the pipeline, see Makefile), scenarios/ (side
                  studies quoted in the paper), db/ (the database as read by the scripts), results/ (working copy)
  pwa/            partial-wave analysis: re-implementation of the 2018 PWA and its update (scenarios/
                  main_fits_legacynorm.py is the PWA fit of the paper); results/, validation/, figs/
  radcorr/        one-loop QED (radiative) corrections and their effect on the fits
  sigma_tot/      total γp → γp cross section and its contribution to the Baldin integral
  delta_plots/    the Δ-region figures and χ² tables (needs the executable, see below)
  data/           the curated source tables db_*.csv with their READMEs and build scripts
  notes/          fit_report.md and pwa_report.md: the working notes that describe the two pipelines
  figstyle*.mplstyle  shared matplotlib styles
results/        frozen reference copies of code/fits/results (fits_L.json: every fit of the paper; tables.md/.tex;
                error_budget.json; delta_region.json; dataset_table.*), code/pwa/results and code/radcorr/results.
                The copies under code/ are the working copies that a rerun overwrites; initially they are identical.
paper/figs/     empty target directory for the scripts that copy figures to the paper
```

## Requirements

Python 3.13 with numpy, scipy, pandas, matplotlib (and mpmath for `code/radcorr` and `code/code_v2/tools`);
the versions used are in `requirements.txt` (`python3 -m pip install -r requirements.txt`).  GNU make.  The
Fortran compiler and LoopTools are needed only to regenerate the amplitude tables (below).  The scripts use
`python3` from the PATH (`PY=... make` to override).

## Rerunning the fits from the amplitude tables (no Fortran)

All commands are run in `code/fits`.

```
cd code/fits
python3 scripts/s03_fits_L.py            # all fits below threshold (Table "fits", variants, scans; ~1 min);
                                         #   writes results/fits_L.json and prints the fits
python3 scenarios/nnlo_prediction_table.py   # the parameter-free NNLO prediction versus the data (Table "pred")
make                                     # the complete pipeline: database assembly (s01), driver tables (s02,
                                         #   all found in amplitudes/), forward validation (s00), fits (s03),
                                         #   truncation correlation lengths (s04), Δ region (s05), figures (s06)
                                         #   and tables incl. the error budget (s07); ~40 s
python3 scripts/s07b_error_budget_table.py   # the error-budget table of Appendix D from results/error_budget.json
```

Expected numbers (`results/tables.md`, `results/fits_L.json` entry `F2_d_all_thB`): NNLO fit α_E1 = 11.59(47),
β_M1 = 2.94(43), χ² = 219.7 for 263 degrees of freedom (283 points); Baldin fit (`F1_d_all_thB`) 11.34(45),
2.83(43); parameter-free NNLO prediction χ² = 227.2 for the 283 points with normalizations and nuisance parameters
floated and the truncation covariance (`scenarios/nnlo_prediction_table.txt`, line GLOBAL).

Other analyses that run from the shipped tables and results:

```
cd code/radcorr   && python3 run_radcorr_fit.py        # QED corrections: shifts of the fits (results/fit_shifts.md)
cd code/sigma_tot && python3 sigma_tot.py && python3 plot_sigma_tot.py   # total cross section, Baldin integral
cd code/pwa       && python3 scenarios/main_fits_legacynorm.py           # the PWA fit of the paper (~5 s)
cd code/fits      && python3 scripts/s08_fig_pred_ratio.py && python3 scripts/s09_fig_coverage.py \
                  && python3 scripts/s10_fig_pred_S2.py    # further figures of the paper
```

Scripts that need the executable: `code/delta_plots/make_all.sh` (the Δ-region figures; its own amplitude cache,
0.9 GB, is not included), `code/pwa/scenarios/multipoles_figure.py` (multipoles of the BχPT amplitude), the
`code/pwa` validation against the Fortran code, and the tools in `code/code_v2`.

## Rebuilding the amplitude tables with the Fortran driver

`code/code_v2` is the exported git tree (commit in `GIT_COMMIT.txt`; no binaries).  Requirements: a Fortran
compiler and the LoopTools library (version 2.16 was used; https://feynarts.de/looptools/), built with the same
compiler.  The primary build is `makefile.gfortran-native` (gfortran; flags `-std=legacy
-fallow-argument-mismatch -cpp -ffree-line-length-none`), the cross-check build `makefile.ifort` (Intel ifort).
Edit the compiler path (`FC`) and the LoopTools path (`LT`) at the top of the makefile, then

```
cd code/code_v2
make -f makefile.gfortran-native          # -> execs/ampl.exe
tests/run_regression.sh                   # regression suite (~20 s)
```

`code/fits/lib/bchpt.py` looks for the executable at `code/code_v2/execs/ampl.exe` (override with the environment
variable `AMPL_EXE`).  A table is recomputed only when it is missing from `amplitudes/`: delete the tables to be
regenerated (or the whole directory; `make clean-cache` in `code/fits` removes `code/fits/cache`, not
`amplitudes/`) and run the pipeline.  The tables are keyed by the driver inputs, so the executable used is
recorded in the `.out.json` sidecars rather than in the key; set `BCHPT_CODE_ID=<git commit>` when regenerating,
so that the sidecars record the code version.  Regenerating everything takes a few minutes (1.4 ms per point).

## Changes relative to the author's working tree

Applied by `patch_release.py` of the build (every change is marked `[release]` in the code):

1. `code/fits/lib/bchpt.py`: the amplitude tables are read from `<release>/amplitudes` (environment variable
   `BCHPT_CACHE` overrides); the cache key is the hash of the driver inputs only, not of the executable's path,
   size and modification time, so that the shipped tables are found without the Fortran code; a missing table
   without an executable raises an error that says so; the sidecar records the executable, the code version and
   the full driver input.
2. `code/fits/Makefile` and `code/delta_plots/make_all.sh` use `python3` instead of the author's virtual
   environment.
3. Absolute paths of the author's machine in `code/sigma_tot/*.py`, `code/pwa/{make_figures,sumrules,plb_inputs,
   build_update_db}.py`, `code/pwa/scenarios/multipoles_figure.py` and `code/fits/scenarios/pi0ff/bv2mult.py` are
   replaced by paths relative to the file (`_WORK` = `code/`, `_ROOT` = the package root); remaining mentions in
   docstrings, notes and logs read `<release>/...`.
4. An empty `paper/figs/delta` directory is provided for the scripts that copy figures to the paper.

The scripts `code/pwa/digitize_plb.py` and `code/pwa/scenarios/klp18_fig4_digitise.py` read figure files of the
2018 PWA paper that are not included; their outputs are in `code/pwa/results`.

## Not included

The raw downloads behind the curated tables (`work/data/raw`, 155 MB: HEPData records, arXiv sources, theses) and
the text extractions of the theses; the reference papers of `code/radcorr/refs`; the Δ-region amplitude cache of
`code/delta_plots` (0.9 GB); the binaries, build directories and run directories of the Fortran code (including
the copy under `code/pwa/fortran` used for the PWA validation, of which the sources are kept); the archived
version-1 results of the fits; the LaTeX sources of the paper.

## Reproducibility test of this package

The tables in `amplitudes/` were generated by replaying the complete pipeline (`make all`, the scenarios,
`sigma_tot.py`, the figure scripts) in a copy of this package with the executable built from `code/code_v2` at
the commit of `GIT_COMMIT.txt`.  The tables of the main database computed by that executable are identical, value
by value, to the tables used for the paper.  The package as shipped was then tested in a fresh copy without any
executable (2026-10-15, macOS, Python 3.13.15, package versions of `requirements.txt`):

| step | result |
|---|---|
| `python3 scripts/s03_fits_L.py` (13 s) | NNLO fit α_E1 = 11.59(47), β_M1 = 2.94(43), χ²/dof = 219.7/263; Baldin fit 11.34(45), 2.83(43), 214.8/244; all 219 fits of `fits_L.json` identical to the paper's (23 314 numbers, maximal difference 0) |
| `python3 scenarios/nnlo_prediction_table.py` | NNLO prediction: χ² = 227.2 for the 283 points (normalizations and nuisances floated, truncation covariance); table identical |
| `make all` (41 s) | database, forward validation, fits, correlation lengths, Δ region, figures, tables and error budget regenerated; `error_budget.json`, `delta_region.json`, `theory_lengths.json` identical |
| `s07b`, `s08`, `s09`, `s10`; `sigma_tot.py`; `run_radcorr_fit.py` (9 s); PWA `main_fits_legacynorm.py` (5 s) | run without error; `sigma_tot.json`, `fit_shifts.json` and `main_fits_legacynorm.json` identical to the shipped results (apart from recorded file paths) |
