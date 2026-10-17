#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# build_looptools_dq.sh -- LoopTools 2.16 with a quad-precision (REAL*16) core behind the usual
# double-precision Fortran interface ("LoopTools-dq"); see README.md in this directory.
#
# usage:  build_looptools_dq.sh PREFIX [LoopTools-2.16.tar.gz]
#   PREFIX  : install directory -> PREFIX/lib/libooptools.a, PREFIX/include/looptools.h
#   tarball : the official LoopTools-2.16.tar.gz (MD5 checked); downloaded from feynarts.de if omitted
#   FC=...  : Fortran compiler (default gfortran); it must support REAL*16 and be the one used for the code
#
# Steps: (1) plain --quad build, to list the public (non-lj) symbols; (2) rebuild with every public name
# prefixed by q through cpp defines appended to externals.h; (3) check that no unprefixed public symbol
# is left; (4) add the double-precision front end ltdq.f90; (5) run the stability self-check test_dq.f90.
# ---------------------------------------------------------------------------
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PREFIX=${1:?usage: $0 PREFIX [LoopTools-2.16.tar.gz]}
TAR=${2:-}
FC=${FC:-gfortran}
MD5=9503320a00e7c16b4ac1fd20749d3feb          # https://feynarts.de/looptools/ (tarball of 6 Nov 2024)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/ltdq.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
[ "$(uname -s)" = Darwin ] && export SDKROOT=${SDKROOT:-$(xcrun --sdk macosx --show-sdk-path)}

if [ -z "$TAR" ]; then
  TAR=$WORK/LoopTools-2.16.tar.gz
  curl -fsSL https://feynarts.de/looptools/LoopTools-2.16.tar.gz -o "$TAR"
fi
sum=$( (md5sum "$TAR" 2>/dev/null || md5 -r "$TAR") | cut -d' ' -f1)
[ "$sum" = "$MD5" ] || { echo "MD5 of $TAR is $sum, expected $MD5" >&2; exit 1; }

# public symbols of a library: defined or common Fortran names that are not LoopTools-internal (lj...)
pubsyms() { nm -g "$1" | awk '$2 ~ /^[TDBCS]$/ {print $3}' | sed -e 's/^_//' | grep -v '^lj' | grep '_$' | sed 's/_$//' | sort -u; }
# configure + make lib in $LT, log to $WORK/build$1.log, show its tail on failure
build() { (cd "$LT" && FC=$FC ./configure --quad --prefix="$WORK/q$1" && make lib) > "$WORK/build$1.log" 2>&1 \
          || { tail -30 "$WORK/build$1.log" >&2; exit 1; }; }

echo "== (1) plain quad build"
tar -xzf "$TAR" -C "$WORK"; LT=$WORK/LoopTools-2.16
build 1
NAMES=$(pubsyms "$(ls "$LT"/build-quad/libooptools-quad.a)")
echo "   $(echo "$NAMES" | wc -l | tr -d ' ') public names (routines and common blocks)"

echo "== (2) rebuild with the public names prefixed by q"
rm -rf "$LT"; tar -xzf "$TAR" -C "$WORK"
for n in $NAMES; do                     # Fortran is case-insensitive, cpp is not: define every spelling
  case $n in d0|e0) continue;; esac     # would also hit the exponent of literals like 1D0; remapped in defs.h
  grep -rhoiw --include='*.F' --include='*.h' --include='*.c' "$n" "$LT/src" | sort -u | sed 's/.*/#define & q&/'
done > "$LT/src/include/qrename.h"
echo '#include "qrename.h"' >> "$LT/src/include/externals.h"
for n in $NAMES; do echo "#define ${n}_ q${n}_"; done >> "$LT/src/include/cexternals.h"   # C side (common /ltcache/)
sed -i.orig -e 's/^#define XD0 D0$/#define XD0 qD0/' -e 's/^#define XE0 E0$/#define XE0 qE0/' "$LT/src/include/defs.h"
build 2
QLIB=$(ls "$LT"/build-quad/libooptools-quad.a)

echo "== (3) check: no unprefixed public symbol left (defined or referenced)"
left=$(nm "$QLIB" | awk '{print $NF}' | sed -e 's/^_//' | grep -E '_$' | sed 's/_$//' | sort -u | grep -xF "$NAMES" || true)
[ -z "$left" ] || { echo "unprefixed public symbols remain: $left" >&2; exit 1; }

echo "== (4) double-precision front end"
mkdir -p "$PREFIX/lib" "$PREFIX/include"
cp "$QLIB" "$PREFIX/lib/libooptools.a"
cp "$LT/src/include/looptools.h" "$PREFIX/include/"   # same header as the double build (types are cpp macros)
(cd "$WORK" && $FC -O2 -c "$HERE/ltdq.f90" -o ltdq.o)
ar rs "$PREFIX/lib/libooptools.a" "$WORK/ltdq.o" 2> /dev/null

echo "== (5) stability self-check"
$FC -O2 -cpp -ffree-line-length-none -I"$PREFIX/include" "$HERE/test_dq.f90" "$PREFIX/lib/libooptools.a" -o "$WORK/test_dq"
"$WORK/test_dq" | grep test_dq
echo "installed: $PREFIX/lib/libooptools.a  $PREFIX/include/looptools.h"
