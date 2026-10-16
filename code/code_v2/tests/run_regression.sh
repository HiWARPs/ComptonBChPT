#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Regression tests of the NNLO BChPT Compton code, version 2 (work/code_v2).
#
# usage: tests/run_regression.sh [--build] [a b c d e f g]
#   --build  rebuild execs/ampl.exe (gfortran) and execs/ampl_ifort.exe (ifort) first
#   a  static polarisabilities vs LMP2015 Tables 3-4 and the analytic App. C formulas
#   b  Thomson limit, Xi/mu independence (bug B1), NNLO vs Born+LEX at 0.5-20 MeV
#   c  unitarity (multipoles vs tree-level photoproduction, per isospin channel with physical
#      masses; forward Im A1/Im A3 at theta=0; gamma N Delta vertex Im parts vs PV2006/Watson)
#   d  v1 (tests/ref_bin/ampl_v1.exe) vs v2 below/above threshold; intentional changes listed
#   e  gfortran vs ifort
#   f  fit-driver exactness (selftest, driver vs main program)
#   g  pi0 gamma gamma transition form factor (&v2 pi0ff): off = default, on leaves static pols unchanged, vs LoopTools C0
# Requirements: gfortran-16 + LoopTools (makefile.gfortran-native), ifort (makefile.ifort),
# Python venv with numpy/scipy/mpmath: <release>/.venv
# Output: tests/out/regression.log ; exit code = number of FAILs.  Runtime ~1 min (+ builds).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
PY=${PYTHON:-python3}
if [ "${1:-}" = "--build" ]; then
  shift
  ( cd "$ROOT" && make -f makefile.gfortran-native >/dev/null ) || { echo "gfortran build failed"; exit 99; }
  ( cd "$ROOT" && make -f makefile.ifort >/dev/null 2>&1 ) || echo "WARNING: ifort build failed (test e will fail)"
fi
[ -x "$ROOT/execs/ampl.exe" ] || { echo "missing execs/ampl.exe (use --build)"; exit 99; }
[ -x "$HERE/ref_bin/ampl_v1.exe" ] || echo "WARNING: tests/ref_bin/ampl_v1.exe (v1 reference) missing: test d will fail"
mkdir -p "$HERE/out"
"$PY" "$HERE/regression.py" "$@" 2>&1 | tee "$HERE/out/regression.log"
rc=${PIPESTATUS[0]}
rm -rf "$HERE/out/tmp"
exit $rc
