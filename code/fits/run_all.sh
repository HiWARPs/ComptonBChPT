#!/usr/bin/env bash
# One-command reproduction (same as `make`).  Usage: ./run_all.sh [path/to/ampl.exe]
set -euo pipefail
cd "$(dirname "$0")"
[ $# -ge 1 ] && export AMPL_EXE="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
make all
