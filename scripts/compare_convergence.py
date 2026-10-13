import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from qevsdftk.parse_qe import parse_qe_energy

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "convergence"
RY_TO_EV = 13.605693122994
TOL_CODES_MEV = 0.1          # meV/atom, QE vs DFTK
CHECK_MIN_ECUT_RY = 40.0     # cases below this cutoff are reported but not checked
CRITERIA_MEV = [1.0, 0.1]

params = json.loads((SYS / "params.json").read_text())
cdata = json.loads((SYS / "cases.json").read_text())
cases = {c["name"]: c for c in cdata["cases"]}
groups = cdata["groups"]
dftk_data = json.loads((SYS / "dftk" / "convergence_results.json").read_text())
dftk = {r["name"]: r for r in dftk_data["results"]}
natoms = {k: len(v["positions_crystal"]) for k, v in params["systems"].items()}
systems = list(params["systems"])

rows = {}
for name, c in cases.items():
    n = natoms[c["system"]]
    eq = parse_qe_energy(SYS / "qe" / f"{name}.out")
    ed = dftk[name]["total_energy_Ry"]
    rows[name] = {**c, "E_qe_Ry": eq, "E_dftk_Ry": ed,
                  "dE_meV": (ed - eq) * RY_TO_EV * 1000 / n}


def scan(system, group, key):
    """Energies relative to the largest parameter of the scan (per code), meV/atom."""
    names = sorted(groups[system][group], key=lambda nm: rows[nm][key])
    n = natoms[system]
    ref = rows[names[-1]]
    return [{
        "name": nm, key: rows[nm][key],
        "dE_qe_meV": (rows[nm]["E_qe_Ry"] - ref["E_qe_Ry"]) * RY_TO_EV * 1000 / n,
        "dE_dftk_meV": (rows[nm]["E_dftk_Ry"] - ref["E_dftk_Ry"]) * RY_TO_EV * 1000 / n,
    } for nm in names]


def converged_at(points, key, tol):
    """Smallest value from which all larger values stay within tol (both codes)."""
    ok = [max(abs(p["dE_qe_meV"]), abs(p["dE_dftk_meV"])) < tol for p in points]
    for i in range(len(points) - 1):   # the last point is the reference itself
        if all(ok[i:]):
            return points[i][key]
    return None


summary = {}
print(f"DFTK {dftk_data['dftk_version']} vs QE")
for system in systems:
    base = params["systems"][system]["base"]
    summary[system] = {}
    for group, key, label in [("ecut_scan", "ecut_ry", "Ecut (Ry)"), ("k_scan", "nk", "k-grid")]:
        pts = scan(system, group, key)
        for tol in CRITERIA_MEV:
            val = converged_at(pts, key, tol)
            text = f"{val}" if val is not None else "not reached within the scan"
            print(f"  {system}: {label} converged to {tol} meV/atom at: {text}")
            summary[system][f"{group}_{tol}_meV"] = val
    print(f"  {system}: base settings were Ecut = {base['ecut_ry']} Ry, k-grid = {base['nk']}")

max_all = max(abs(r["dE_meV"]) for r in rows.values())
checked = [r for r in rows.values() if r["ecut_ry"] >= CHECK_MIN_ECUT_RY]
max_codes = max(abs(r["dE_meV"]) for r in checked)
worst = max(rows.values(), key=lambda r: abs(r["dE_meV"]))
print(f"Max |E_DFTK - E_QE| over all {len(rows)} cases: {max_all:.4f} meV/atom ({worst['name']})")
print(f"Max |E_DFTK - E_QE| over the {len(checked)} cases with Ecut >= {CHECK_MIN_ECUT_RY:.0f} Ry: "
      f"{max_codes:.4f} meV/atom")
# ---------------- Figure ----------------
fig, axes = plt.subplots(len(systems), 3, figsize=(15, 4.3 * len(systems)), squeeze=False)
for i, system in enumerate(systems):
    for j, (group, key, xlabel) in enumerate(
            [("ecut_scan", "ecut_ry", "Ecut (Ry)"), ("k_scan", "nk", "k-grid (n x n x n)")]):
        pts = scan(system, group, key)
        x = [p[key] for p in pts]
        ax = axes[i][j]
        ax.plot(x, [p["dE_qe_meV"] for p in pts], "o-", label="QE")
        ax.plot(x, [p["dE_dftk_meV"] for p in pts], "x--", label="DFTK")
        ax.axhspan(-1.0, 1.0, color="gray", alpha=0.15)
        ax.set_yscale("symlog", linthresh=0.1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("E - E(largest) (meV/atom)")
        ax.set_title(f"{system.upper()}: {xlabel} convergence")
        ax.legend()
    names = groups[system]["ecut_scan"] + [n for n in groups[system]["k_scan"]
                                           if n not in groups[system]["ecut_scan"]]
    ax = axes[i][2]
    ax.bar(range(len(names)), [rows[n]["dE_meV"] for n in names])
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_ylabel("E$_{DFTK}$ - E$_{QE}$ (meV/atom)")
    ax.set_title(f"{system.upper()}: difference between codes")
fig.tight_layout()
(ROOT / "results").mkdir(exist_ok=True)
fig.savefig(ROOT / "results" / "convergence.png", dpi=150)

(ROOT / "results" / "convergence_comparison.json").write_text(json.dumps({
    "dftk_version": dftk_data["dftk_version"],
    "max_abs_dE_codes_meV_per_atom": max_codes,
    "max_abs_dE_all_cases_meV_per_atom": max_all,
    "check_min_ecut_ry": CHECK_MIN_ECUT_RY,
    "converged_settings": summary,
    "cases": list(rows.values()),
}, indent=2))

ok = max_codes < TOL_CODES_MEV
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
