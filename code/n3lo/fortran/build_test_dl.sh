#!/bin/bash
# builds test_dl.exe from fortran/p4 (Gram routines) and fortran/dl against LoopTools-dq
set -e
cd "$(dirname "$0")"
export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
FC=/opt/homebrew/bin/gfortran-16
LT=$HOME/.local/looptools-dq
mkdir -p build/dl build/p4
for f in p4/*.F90; do o=build/p4/$(basename ${f%.F90}).o; [ $o -nt $f ] || $FC -O1 -w -ffree-line-length-none -Jbuild/p4 -c $f -o $o; done
$FC -O1 -w -ffree-line-length-none -Jbuild/dl -c dl/CovDeltaLoops.F90 -o build/dl/CovDeltaLoops.o
ls dl/dlg_*.F90 | xargs -P 8 -I{} sh -c 'f={}; o=build/dl/$(basename ${f%.F90}).o; [ $o -nt $f ] || '"$FC"' -O1 -w -ffree-line-length-none -c $f -o $o'
$FC -O1 -w -Ibuild/dl -o test_dl.exe test_dl.f90 build/dl/*.o build/p4/*.o $LT/lib/libooptools.a
