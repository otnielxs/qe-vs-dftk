#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NP="${NP:-8}"

python "$ROOT/scripts/make_al_inputs.py"

cd "$ROOT/systems/al/qe"
for f in *.in; do
  echo "Running QE: $f"
  mpirun -np "$NP" pw.x -in "$f" > "${f%.in}.out"
done

cd "$ROOT/systems/al/dftk"
julia --project="$ROOT" al.jl

cd "$ROOT"
python scripts/compare_al.py
