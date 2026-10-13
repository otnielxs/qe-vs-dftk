#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export OMP_NUM_THREADS=1   # avoid hidden threading inside QE

python "$ROOT/scripts/make_timing_inputs.py"

read -r -a CORES <<< "$(python -c "import json; print(*json.load(open('$ROOT/systems/timing/params.json'))['cores'])")"
REPEATS="$(python -c "import json; print(json.load(open('$ROOT/systems/timing/params.json'))['repeats'])")"
WORKLOADS="si8 al"

cd "$ROOT/systems/timing/qe"
for w in $WORKLOADS; do
  for n in "${CORES[@]}"; do
    for r in $(seq 1 "$REPEATS"); do
      echo "QE: $w, $n processes, repeat $r"
      mpirun -np "$n" pw.x -nk "$n" -in "$w.in" > "${w}_np${n}_r${r}.out"
    done
  done
done

cd "$ROOT/systems/timing/dftk"
for w in $WORKLOADS; do
  for n in "${CORES[@]}"; do
    echo "DFTK (MPI): $w, $n ranks"
    "$HOME/.julia/bin/mpiexecjl" --project="$ROOT" -n "$n" julia --project="$ROOT" timing_mpi.jl "$w"
  done
done

cd "$ROOT"
python scripts/compare_timing.py
