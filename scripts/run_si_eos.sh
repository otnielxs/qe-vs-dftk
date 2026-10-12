#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

python "$ROOT/scripts/make_qe_eos_inputs.py"

cd "$ROOT/systems/si_eos/qe"
for f in a_*.in; do
  echo "Running QE: $f"
  mpirun -np 8 pw.x -in "$f" > "${f%.in}.out"
done

cd "$ROOT/systems/si_eos/dftk"
julia --project="$ROOT" eos.jl

cd "$ROOT"
python scripts/compare_si_eos.py
