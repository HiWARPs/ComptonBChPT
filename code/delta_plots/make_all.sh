#!/usr/bin/env bash
# ---------------------------------------------------------------------------------------------
# make_all.sh -- regenerate all Delta-region figures, the chi^2 table and PDF->PNG renders.
#
#   cd work/delta_plots && ./make_all.sh      # defaults: code_v2; main curve = parameter-free NNLO prediction
#                                             # (DP_MAIN=pred); the 'NNLO fit' curve reads alpha, beta from
#                                             # work/fits/results/fits_L.json ['F2_d_all_thB'] automatically
#   DP_MAIN=fit ./make_all.sh                 # legacy: main curve = NNLO fit with its band (figures until 2026-10-15)
#   DP_EXTRA=novc ./make_all.sh               # legacy: dotted curve = prediction without the vertex corrections (figures until
#                                             # 2026-10-17); default DP_EXTRA=run3 = prediction with the running Delta parameters
#                                             # (code_v2 &v2 irun=3, static polarisabilities unchanged), also in the angular
#                                             # distributions (DP_EXTRA_THETA=0 switches it off there)
#   DP_FITS_JSON=../fits/results/fits_L.json DP_FIT_KEY=F2_d_all_thB ./make_all.sh   # other fit file/key
#   AMPL_EXE=/path/to/ampl.exe DP_V2OPTS="pimassphys=t, nmassphys=t" ./make_all.sh
#
# Outputs: paper/figs/delta/*.pdf ; work/delta_plots/results/* ; work/delta_plots/png/pdf/*.png (renders of
# the final PDFs at 150 dpi, for inspection).  Runtime ~5 min from an empty cache (cache/ keyed by inputs; one extra
# complete run for the running curve).
# ---------------------------------------------------------------------------------------------
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
PY=${PY:-python3}
cd "$HERE"
"$PY" make_figs.py
mkdir -p png/pdf
for f in "$ROOT"/paper/figs/delta/*.pdf; do
  b=$(basename "${f%.pdf}")
  pdftoppm -r 150 -png -singlefile "$f" "png/pdf/$b"
done
echo "PDFs:   $ROOT/paper/figs/delta/"
echo "tables: $HERE/results/chi2_table.md (.csv, .tex)"
echo "renders: $HERE/png/pdf/"
