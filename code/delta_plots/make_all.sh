#!/usr/bin/env bash
# ---------------------------------------------------------------------------------------------
# make_all.sh -- regenerate all Delta-region figures, the chi^2 table and PDF->PNG renders.
#
#   cd work/delta_plots && ./make_all.sh      # defaults: code_v2; main curve = parameter-free NNLO prediction
#                                             # (DP_MAIN=pred); the 'NNLO fit' curve reads alpha, beta from
#                                             # work/fits/results/fits_L.json ['F2_d_all_thB'] automatically
#   DP_MAIN=fit ./make_all.sh                 # legacy: main curve = NNLO fit with its band (figures until 2026-10-15)
#   DP_FITS_JSON=../fits/results/fits_L.json DP_FIT_KEY=F2_d_all_thB ./make_all.sh   # other fit file/key
#   AMPL_EXE=/path/to/ampl.exe DP_V2OPTS="pimassphys=t, nmassphys=t" ./make_all.sh
#   ./make_all.sh --check-v1      # additionally: chi^2 regression check against fit_report Sect. 5.1 (v1)
#
# Outputs: paper/figs/delta/*.pdf ; work/delta_plots/results/* ; work/delta_plots/png/pdf/*.png (renders of
# the final PDFs at 150 dpi, for inspection).  Runtime ~4 min from an empty cache (cache/ keyed by inputs).
# ---------------------------------------------------------------------------------------------
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
PY=${PY:-python3}
cd "$HERE"
"$PY" make_figs.py
if [ "${1:-}" = "--check-v1" ]; then "$PY" check_v1_chi2.py; fi
mkdir -p png/pdf
for f in "$ROOT"/paper/figs/delta/*.pdf; do
  b=$(basename "${f%.pdf}")
  pdftoppm -r 150 -png -singlefile "$f" "png/pdf/$b"
done
echo "PDFs:   $ROOT/paper/figs/delta/"
echo "tables: $HERE/results/chi2_table.md (.csv, .tex)"
echo "renders: $HERE/png/pdf/"
