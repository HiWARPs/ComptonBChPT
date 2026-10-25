#!/bin/bash
# builds the Fortran test harnesses (covp3.exe, covd.exe) against the code_v2 sources and LoopTools-dq
set -e
cd "$(dirname "$0")"
FC=/opt/homebrew/bin/gfortran-16
export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
LT=$HOME/.local/looptools-dq
FL="-O2 -std=legacy -fallow-argument-mismatch -w -ffree-form -ffree-line-length-none -cpp -I. -I../../code_v2 -I$LT/include"
mkdir -p build
(cd build && $FC $FL -c ../../../code_v2/trigd.f90 -o trigd.o && $FC $FL -c ../../../code_v2/Cov_1loop_As.F90 -o Cov_1loop_As.o && $FC $FL -c ../seccorriso.F90 -o seccorriso.o)
for p in covd covdsec; do $FC $FL -Jbuild -o $p.exe $p.f90 build/Cov_1loop_As.o build/seccorriso.o build/trigd.o $LT/lib/libooptools.a; done
