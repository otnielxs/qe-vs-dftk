#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NP="${NP:-8}"

python "$ROOT/scripts/make_si_bands_inputs.py"

cd "$ROOT/systems/si_bands/qe"
mpirun -np "$NP" pw.x -in scf.in > scf.out
mpirun -np "$NP" pw.x -in bands.in > bands.out

cd "$ROOT/systems/si_bands/dftk"
julia --project="$ROOT" bands.jl

cd "$ROOT"
python scripts/compare_si_bands.py
