#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NP="${NP:-8}"

python "$ROOT/scripts/make_convergence_inputs.py"

cd "$ROOT/systems/convergence/qe"
for f in *.in; do
  echo "Running QE: $f"
  mpirun -np "$NP" pw.x -in "$f" > "${f%.in}.out"
done

cd "$ROOT/systems/convergence/dftk"
julia --project="$ROOT" convergence.jl

cd "$ROOT"
python scripts/compare_convergence.py
