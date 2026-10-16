#!/usr/bin/env bash
# ---------------------------------------------------------------------------------------------
# make_all.sh -- regenerate all Delta-region figures, the chi^2 table and PDF->PNG renders.
#
#   cd work/delta_plots && ./make_all.sh                       # defaults (code_v2, alpha=11.69, beta=2.88)
#   DP_ALPHA=11.67 DP_BETA=2.88 DP_COV_AB=0.218,-0.168,0.187 ./make_all.sh
#   DP_FITS_JSON=../fits/results/fits_L.json DP_DELTA_JSON=../fits/results/delta_region.json ./make_all.sh
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
