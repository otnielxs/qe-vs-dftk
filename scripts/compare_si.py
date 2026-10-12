import json
import sys
from pathlib import Path

import numpy as np

from qevsdftk.parse_qe import parse_qe

ROOT = Path(__file__).resolve().parents[1]
RY_TO_EV = 13.605693122994
N_ATOMS = 2
VBM = 3  # Si has 8 valence electrons -> 4 occupied bands -> index 3 is the VBM at Gamma

TOL_ENERGY_MEV_ATOM = 1.0
TOL_EIG_MEV = 5.0

qe = parse_qe(ROOT / "systems/si/qe/scf.out")
dftk = json.loads((ROOT / "systems/si/dftk/results.json").read_text())

dE = (dftk["total_energy_Ry"] - qe["total_energy_Ry"]) * RY_TO_EV * 1000 / N_ATOMS

n = min(len(qe["gamma_eigenvalues_eV"]), len(dftk["gamma_eigenvalues_eV"]), 8)
e_qe = np.array(qe["gamma_eigenvalues_eV"][:n])
e_dftk = np.array(dftk["gamma_eigenvalues_eV"][:n])
d_eig = ((e_dftk - e_dftk[VBM]) - (e_qe - e_qe[VBM])) * 1000  # meV, relative to VBM

print(f"QE {qe['qe_version']}  vs  DFTK {dftk['dftk_version']}")
print(f"Total energy difference : {dE:+.4f} meV/atom")
print("Gamma eigenvalues relative to VBM (meV difference, DFTK - QE):")
for i, d in enumerate(d_eig, start=1):
    print(f"  band {i}: {d:+.3f}")
print(f"SCF iterations: QE {qe['n_scf_iterations']}, DFTK {dftk['n_scf_iterations']}")

summary = {
    "delta_E_meV_per_atom": dE,
    "delta_eig_meV": d_eig.tolist(),
    "qe": qe,
    "dftk": dftk,
}
out = ROOT / "results" / "si_comparison.json"
out.write_text(json.dumps(summary, indent=2))

ok = abs(dE) < TOL_ENERGY_MEV_ATOM and np.all(np.abs(d_eig) < TOL_EIG_MEV)
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
