import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "si_eos"
QE_DIR = SYS / "qe"

TEMPLATE = """&CONTROL
  calculation = 'scf'
  prefix = 'si_eos'
  outdir = './tmp'
  pseudo_dir = '../../../pseudos/'
  disk_io = 'none'
/
&SYSTEM
  ibrav = 2, celldm(1) = {a:.6f}, nat = 2, ntyp = 1
  ecutwfc = {ecut:.1f}
  occupations = 'fixed'
/
&ELECTRONS
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
  Si 28.0855 {pseudo}
ATOMIC_POSITIONS (crystal)
  Si 0.00 0.00 0.00
  Si 0.25 0.25 0.25
K_POINTS (automatic)
  {nk} {nk} {nk} 0 0 0
"""

params = json.loads((SYS / "params.json").read_text())
QE_DIR.mkdir(parents=True, exist_ok=True)

# Remove stale files from earlier runs with different parameters
for old in list(QE_DIR.glob("a_*.in")) + list(QE_DIR.glob("a_*.out")):
    old.unlink()

for i, a in enumerate(params["a_values_bohr"], start=1):
    text = TEMPLATE.format(
        a=a, ecut=params["ecut_ry"], nk=params["nk"], pseudo=params["pseudo"]
    )
    (QE_DIR / f"a_{i:02d}.in").write_text(text)

print(f"Wrote {len(params['a_values_bohr'])} QE inputs to {QE_DIR}")
