import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "timing"
QE_DIR = SYS / "qe"

TEMPLATE = """&CONTROL
  calculation = 'scf'
  prefix = 'timing'
  outdir = './tmp'
  pseudo_dir = '../../../pseudos/'
  disk_io = 'none'
/
&SYSTEM
  ibrav = {ibrav}, celldm(1) = {a:.6f}, nat = {nat}, ntyp = 1
  ecutwfc = {ecut:.1f}
{occ}/
&ELECTRONS
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
  {el} {mass} {pseudo}
ATOMIC_POSITIONS (crystal)
{pos}K_POINTS (automatic)
  {nk} {nk} {nk} 0 0 0
"""


def occ_lines(occ):
    if occ["kind"] == "fixed":
        return "  occupations = 'fixed'\n"
    return (f"  occupations = 'smearing'\n  smearing = '{occ['smearing']}'\n"
            f"  degauss = {occ['width_ry']:.4f}\n")


params = json.loads((SYS / "params.json").read_text())
QE_DIR.mkdir(parents=True, exist_ok=True)
for old in list(QE_DIR.glob("*.in")) + list(QE_DIR.glob("*.out")):
    old.unlink()

for name, w in params["workloads"].items():
    pos = "".join(f"  {w['element']} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n"
                  for p in w["positions_crystal"])
    text = TEMPLATE.format(
        ibrav=1 if w["cell"] == "cubic" else 2, a=w["a_bohr"],
        nat=len(w["positions_crystal"]), ecut=w["ecut_ry"],
        occ=occ_lines(w["occupations"]), el=w["element"], mass=w["mass"],
        pseudo=params["pseudos"][w["pseudo"]], pos=pos, nk=w["nk"],
    )
    (QE_DIR / f"{name}.in").write_text(text)

print(f"Wrote QE inputs for {list(params['workloads'])}")
