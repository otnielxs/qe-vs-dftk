#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

cd "$ROOT/systems/si/qe"
mpirun -np 8 pw.x -in scf.in > scf.out

cd "$ROOT/systems/si/dftk"
julia --project="$ROOT" scf.jl

cd "$ROOT"
python scripts/compare_si.py
