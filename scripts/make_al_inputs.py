import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "al"
QE_DIR = SYS / "qe"

TEMPLATE = """&CONTROL
  calculation = 'scf'
  prefix = 'al'
  outdir = './tmp'
  pseudo_dir = '../../../pseudos/'
  disk_io = 'none'
  verbosity = 'high'
/
&SYSTEM
  ibrav = 2, celldm(1) = {a:.6f}, nat = 1, ntyp = 1
  ecutwfc = {ecut:.1f}
  occupations = 'smearing'
  smearing = '{smearing}'
  degauss = {width:.4f}
/
&ELECTRONS
  conv_thr = 1.0d-10
/
ATOMIC_SPECIES
  Al 26.9815 {pseudo}
ATOMIC_POSITIONS (crystal)
  Al 0.00 0.00 0.00
K_POINTS (automatic)
  {nk} {nk} {nk} 0 0 0
"""

params = json.loads((SYS / "params.json").read_text())

cases = {}
groups = {"base": [], "smearing_scan": [], "k_scan": []}


def add(group, smearing, width, nk):
    name = f"{smearing}_w{width:.3f}_k{nk}"
    cases.setdefault(
        name, {"name": name, "smearing": smearing, "width_ry": width, "nk": nk}
    )
    groups[group].append(name)


b = params["base"]
add("base", b["smearing"], b["width_ry"], b["nk"])

s = params["smearing_scan"]
for sm in s["types"]:
    for w in s["widths_ry"]:
        add("smearing_scan", sm, w, s["nk"])

k = params["k_scan"]
for nk in k["nks"]:
    add("k_scan", k["smearing"], k["width_ry"], nk)

(SYS / "cases.json").write_text(
    json.dumps({"cases": list(cases.values()), "groups": groups}, indent=1)
)

QE_DIR.mkdir(parents=True, exist_ok=True)
for old in list(QE_DIR.glob("*.in")) + list(QE_DIR.glob("*.out")):
    old.unlink()

for c in cases.values():
    text = TEMPLATE.format(
        a=params["a_bohr"], ecut=params["ecut_ry"], pseudo=params["pseudo"],
        smearing=c["smearing"], width=c["width_ry"], nk=c["nk"],
    )
    (QE_DIR / f"{c['name']}.in").write_text(text)

print(f"Wrote {len(cases)} unique cases to {SYS}")
