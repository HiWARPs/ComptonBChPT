# Proton polarisability fits with the NNLO BChPT code -- full pipeline.
#   make            : everything (database, driver runs, fits, Delta region, figures, tables)
#   make clean-cache: remove cached driver outputs (forces recomputation)
# Executable: ../code_v2/execs/ampl.exe (profile v2); override with AMPL_EXE=...; in.dat template override:
# AMPL_TEMPLATE (optional).
PY ?= python3
# FIT_PROFILE = v2 (work/code_v2, vertex corrections + physical masses)
export FIT_PROFILE ?= v2
export FIT_RESULTS ?= results
export FIT_FIGS ?= figs
export AMPL_EXE
export AMPL_TEMPLATE

all: db driver validate fitsL theorylen delta combined figs tables

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
# main BChPT fit: alpha, beta, R_EM up to 360 MeV (needs the Delta band of work/delta_plots/results/data_with_theory.csv)
combined: fitsL
	$(PY) scenarios/combined_fit.py > $(FIT_RESULTS)/combined_log.txt
figs: fitsL delta combined
	$(PY) scripts/s06_figures.py > $(FIT_RESULTS)/s06_log.txt
tables: fitsL delta
	$(PY) scripts/s07_tables.py > /dev/null

clean-cache:
	rm -f cache/*.out cache/*.json

.PHONY: all db driver validate fitsL theorylen delta combined figs tables clean-cache
