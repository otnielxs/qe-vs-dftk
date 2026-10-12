import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "si_bands"
QE_DIR = SYS / "qe"

SCF_TEMPLATE = """&CONTROL
  calculation = 'scf'
  prefix = 'si_bands'
  outdir = './tmp'
  pseudo_dir = '../../../pseudos/'
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

BANDS_TEMPLATE = """&CONTROL
  calculation = 'bands'
  prefix = 'si_bands'
  outdir = './tmp'
  pseudo_dir = '../../../pseudos/'
  verbosity = 'high'
/
&SYSTEM
  ibrav = 2, celldm(1) = {a:.6f}, nat = 2, ntyp = 1
  ecutwfc = {ecut:.1f}
  occupations = 'fixed'
  nbnd = {nbnd}
/
&ELECTRONS
  conv_thr = 1.0d-10
  diago_full_acc = .true.
/
ATOMIC_SPECIES
  Si 28.0855 {pseudo}
ATOMIC_POSITIONS (crystal)
  Si 0.00 0.00 0.00
  Si 0.25 0.25 0.25
K_POINTS tpiba
{nks}
{klist}"""

params = json.loads((SYS / "params.json").read_text())
pts = {k: np.array(v, dtype=float) for k, v in params["points_tpiba"].items()}
path = params["path"]
n = params["points_per_segment"]

# Build the path: n points per segment, each vertex appears once
kpts, label_index = [], []
for i in range(len(path) - 1):
    k0, k1 = pts[path[i]], pts[path[i + 1]]
    label_index.append(len(kpts))
    for t in np.arange(n) / n:
        kpts.append(k0 + t * (k1 - k0))
label_index.append(len(kpts))
kpts.append(pts[path[-1]])
kpts = np.array(kpts)

dist = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kpts, axis=0), axis=1))])

(SYS / "kpath.json").write_text(json.dumps({
    "kpoints_tpiba": kpts.tolist(),
    "distance": dist.tolist(),
    "labels": path,
    "label_index": label_index,
}, indent=1))

QE_DIR.mkdir(parents=True, exist_ok=True)
fmt = dict(a=params["a_bohr"], ecut=params["ecut_ry"], pseudo=params["pseudo"])

(QE_DIR / "scf.in").write_text(SCF_TEMPLATE.format(nk=params["nk_scf"], **fmt))

klist = "\n".join(f"  {k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1.0" for k in kpts)
(QE_DIR / "bands.in").write_text(
    BANDS_TEMPLATE.format(nbnd=params["n_bands"], nks=len(kpts), klist=klist, **fmt)
)

print(f"Wrote kpath.json and QE inputs ({len(kpts)} k-points)")
