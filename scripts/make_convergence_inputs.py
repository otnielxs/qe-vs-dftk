import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "convergence"
QE_DIR = SYS / "qe"

TEMPLATE = """&CONTROL
  calculation = 'scf'
  prefix = 'conv'
  outdir = './tmp'
  pseudo_dir = '../../../pseudos/'
  disk_io = 'none'
/
&SYSTEM
  ibrav = 2, celldm(1) = {a:.6f}, nat = {nat}, ntyp = 1
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
cases, groups = {}, {}


def add(system, group, ecut, nk):
    name = f"{system}_e{int(ecut):03d}_k{int(nk):02d}"
    cases.setdefault(name, {"name": name, "system": system,
                            "ecut_ry": float(ecut), "nk": int(nk)})
    groups.setdefault(system, {"ecut_scan": [], "k_scan": []})[group].append(name)


for system, s in params["systems"].items():
    for e in s["ecut_scan_ry"]:
        add(system, "ecut_scan", e, s["base"]["nk"])
    for nk in s["k_scan"]:
        add(system, "k_scan", s["base"]["ecut_ry"], nk)

(SYS / "cases.json").write_text(json.dumps({"cases": list(cases.values()), "groups": groups}, indent=1))

QE_DIR.mkdir(parents=True, exist_ok=True)
for old in list(QE_DIR.glob("*.in")) + list(QE_DIR.glob("*.out")):
    old.unlink()

for c in cases.values():
    s = params["systems"][c["system"]]
    pos = "".join(f"  {s['element']} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n"
                  for p in s["positions_crystal"])
    text = TEMPLATE.format(
        a=s["a_bohr"], nat=len(s["positions_crystal"]), ecut=c["ecut_ry"],
        occ=occ_lines(s["occupations"]), el=s["element"], mass=s["mass"],
        pseudo=params["pseudos"][c["system"]], pos=pos, nk=c["nk"],
    )
    (QE_DIR / f"{c['name']}.in").write_text(text)

print(f"Wrote {len(cases)} unique cases to {SYS}")
