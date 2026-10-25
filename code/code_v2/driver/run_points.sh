#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_points.sh -- fit driver for the NNLO BChPT Compton code (ampl.exe)
#
# usage:  run_points.sh POINTS_FILE [OUT_FILE] [key=value ...]
#
#   POINTS_FILE : lines  "id  E_lab[MeV]  theta[deg]  OBS  [LAB|CM]"   ('#' = comment)
#                 OBS = CS (dsigma/dOmega_lab, nb/sr) | CSCM (c.m., nb/sr) | S3 | S2X | S2Z
#                       | SPAR | SPERP (lab dsigma for 100% linear pol. parallel/perp. to plane, nb/sr)
#                 angle frame: LAB (default) or CM
#   OUT_FILE    : default = <POINTS_FILE basename>.out
#   keys (all optional):
#     alpha=  beta=   (1e-4 fm^3)    gEE= gMM= gEM= gME=  (1e-4 fm^4)
#                   polarisabilities at which the 'value' column is evaluated;
#                   default = the NNLO prediction (code value -999 = "use reference")
#     selftest=t|f  compare quadratic form with a direct full-amplitude evaluation (default f)
#     onlylex=t|f   Born + pi0 pole + contact terms only (reference polarisabilities = 0)
#     ref="refalphaE1=..,refbetaM1=..,..."  override the reference values (default: analytic NNLO)
#     exe=PATH      executable (default ../execs/ampl.exe relative to this script)
#     keep=1        keep the temporary run directory
#   [v2] options (namelist &v2 of ampl.exe; defaults reproduce the 2015 code apart from the bug fixes):
#     order=2|3|3.5|4|4.5  or  order=N2LO|N3LO|N4LO   chiral order of the base amplitude (delta counting):
#                   2 = Born (charges only), 3 = + AMM, pi0 pole, piN loops, 3.5 = N2LO = + Delta pole, piDelta
#                   loops (default; LMP2015), 4 = N3LO = + O(p^4) piN loops and the photon on the Delta charge in the
#                   piDelta loops, 4.5 = N4LO = + gamma N Delta vertex in the piN loops and the one-loop gamma N Delta
#                   vertex in the Delta pole (orders 4, 4.5: proton only; ~1 s per point at 4.5)
#     dvertex=t|f   Im parts of the gamma-N-Delta vertex loops (NLO in the Delta region, LMP2015 Fig. 4)
#     drunning=t|f|0|1|2|3  running Re gM, gE, M_Delta(s), 1+Re tau_R(s) (LMP2015 App. D; higher order, off):
#                   t = 1 = once subtracted at the Delta pole (changes the static polarisabilities),
#                   2 = doubly subtracted (pole and static point s0 = M_N^2: alpha_E1, beta_M1 unchanged),
#                   3 = also vanishing slopes at s0 (all six static pols unchanged; exploratory); numbers set &v2 irun
#     consts=LMP15|PDG24  parameter set (default LMP15 = LMP2015 Table 2)
#     v2="key=val,..."    any further &v2 entries (pimassphys, nmassphys, pi0ff, loopsafe, elex, tausafe, tzfrac, ...)
# ---------------------------------------------------------------------------
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
[ $# -ge 1 ] || { sed -n '2,31p' "$0"; exit 1; }
PTS=$1; shift
OUT=""
if [ $# -ge 1 ] && [[ "$1" != *=* ]]; then OUT=$1; shift; fi
[ -n "$OUT" ] || OUT="$(basename "${PTS%.*}").out"
ALPHA=-999; BETA=-999; GEE=-999; GMM=-999; GEM=-999; GME=-999
SELFTEST=f; ONLYLEX=f; REFS=""; EXE="$HERE/../execs/ampl.exe"; KEEP=0
ORDER=3.5; DVERTEX=""; DRUNNING=""; CONSTS="LMP15"; V2X=""
for kv in "$@"; do
  k=${kv%%=*}; v=${kv#*=}
  case $k in
    alpha) ALPHA=$v;; beta) BETA=$v;; gEE) GEE=$v;; gMM) GMM=$v;; gEM) GEM=$v;; gME) GME=$v;;
    selftest) SELFTEST=$v;; onlylex) ONLYLEX=$v;; ref) REFS=", $v";; exe) EXE=$v;; keep) KEEP=$v;;
    order) case $v in N2LO|NNLO) ORDER=3.5;; N3LO) ORDER=4;; N4LO) ORDER=4.5;; *) ORDER=$v;; esac;; dvertex) DVERTEX=", dvertex=$v";;
    drunning) case $v in [0-9]) DRUNNING=", irun=$v";; *) DRUNNING=", drunning=$v";; esac;;
    consts) CONSTS=$v;; v2) V2X=", $v";;
    *) echo "unknown key $k" >&2; exit 1;;
  esac
done
[ -x "$EXE" ] || { echo "executable $EXE not found (build with: make -f makefile.gfortran-native)" >&2; exit 1; }
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"      # absolute path (we cd into a temp dir)
PTS_ABS="$(cd "$(dirname "$PTS")" && pwd)/$(basename "$PTS")"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/cbchpt_run.XXXXXX")
cp "$HERE/baryon.dat" "$WORK/baryon.dat"
cp "$PTS_ABS" "$WORK/points.txt"
sed -e "s/@ALPHA@/$ALPHA/" -e "s/@BETA@/$BETA/" -e "s/@GEE@/$GEE/" -e "s/@GMM@/$GMM/" \
    -e "s/@GEM@/$GEM/" -e "s/@GME@/$GME/" -e "s/@SELFTEST@/$SELFTEST/" -e "s/@ONLYLEX@/$ONLYLEX/" \
    -e "s/@REFS@/$REFS/" -e "s/@ORDER@/$ORDER/" -e "s/@DVERTEX@/$DVERTEX/" -e "s/@DRUNNING@/$DRUNNING/" \
    -e "s/@CONSTS@/$CONSTS/" -e "s/@V2X@/$V2X/" "$HERE/in_driver.template" > "$WORK/in.dat"
# (executables without the &v2 namelist simply ignore the trailing group)
t0=$(date +%s.%N 2>/dev/null || date +%s)
( cd "$WORK" && "$EXE" > run.log 2>&1 ) || { echo "ampl.exe failed, see $WORK/run.log" >&2; exit 1; }
t1=$(date +%s.%N 2>/dev/null || date +%s)
cp "$WORK/points.out" "$OUT"
[ "$SELFTEST" = "t" ] && cp "$WORK/points.out.selftest" "$OUT.selftest"
grep -E "run_points|selftest" "$WORK/run.log" || true
n=$(grep -vc '^#' "$OUT" || true)
awk -v a="$t0" -v b="$t1" -v n="$n" 'BEGIN{ if(n>0) printf(" wall time %.2f s for %d points (%.2f ms/point)\n", b-a, n, 1000*(b-a)/n) }'
if [ "$KEEP" = "1" ]; then echo " run directory kept: $WORK"; else rm -rf "$WORK"; fi
echo " output: $OUT"
