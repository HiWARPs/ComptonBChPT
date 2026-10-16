# Proton polarisability fits with the NNLO BChPT code -- full pipeline.
#   make            : everything (database, driver runs, fits, Delta region, figures, tables)
#   make clean-cache: remove cached driver outputs (forces recomputation)
# Executable: chosen by FIT_PROFILE (v2: ../code_v2/execs/ampl.exe, v1: ../code/execs/ampl.exe); override with
# AMPL_EXE=...; in.dat template override: AMPL_TEMPLATE (optional).  The archived v1 results (results_v1/, figs_v1/)
# were produced by the v1 pipeline state of 2026-10-14; re-running FIT_PROFILE=v1 into a NEW directory is possible.
PY ?= python3
# FIT_PROFILE = v2 (default: work/code_v2, vertex corrections + physical masses) or v1 (work/code, archived in
# results_v1/, figs_v1/: run `make FIT_PROFILE=v1 FIT_RESULTS=results_v1 FIT_FIGS=figs_v1 all`)
export FIT_PROFILE ?= v2
export FIT_RESULTS ?= results
export FIT_FIGS ?= figs
export AMPL_EXE
export AMPL_TEMPLATE

all: db driver validate fitsL theorylen delta figs tables compare

db:
	$(PY) scripts/s01_build_db.py
driver: db
	$(PY) scripts/s02_run_driver.py
validate: driver
	$(PY) scripts/s00_validate_forward.py > /dev/null
fitsL: driver
	$(PY) scripts/s03_fits_L.py > $(FIT_RESULTS)/s03_log.txt
theorylen: driver
	$(PY) scripts/s04_theory_lengths.py > $(FIT_RESULTS)/s04_log.txt
delta: fitsL
	$(PY) scripts/s05_delta.py > $(FIT_RESULTS)/s05_log.txt
figs: fitsL delta
	$(PY) scripts/s06_figures.py > $(FIT_RESULTS)/s06_log.txt
tables: fitsL delta
	$(PY) scripts/s07_tables.py > /dev/null
compare: tables
	@if [ -f results_v1/fits_L.json ] && [ "$(FIT_PROFILE)" = "v2" ]; then $(PY) scripts/s08_compare_v1_v2.py > /dev/null; fi

clean-cache:
	rm -f cache/*.out cache/*.json

.PHONY: all db driver validate fitsL theorylen delta figs tables compare clean-cache
