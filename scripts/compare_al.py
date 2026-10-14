import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qevsdftk.parse_qe import parse_qe
from qevsdftk.plotstyle import apply_style

apply_style()

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "al"
RY_TO_EV = 13.605693122994

TOL_E_MEV = 0.1    # meV/atom
TOL_BW_MEV = 1.0   # meV

cases_data = json.loads((SYS / "cases.json").read_text())
params = json.loads((SYS / "params.json").read_text())
iv = params["n_semicore_bands"]  # index of the first valence band at Gamma
cases = {c["name"]: c for c in cases_data["cases"]}
groups = cases_data["groups"]

dftk_data = json.loads((SYS / "dftk" / "al_results.json").read_text())
dftk = {r["name"]: r for r in dftk_data["results"]}

rows = []
for name, c in cases.items():
    q = parse_qe(SYS / "qe" / f"{name}.out")
    d = dftk[name]
    bw_q = q["fermi_eV"] - q["gamma_eigenvalues_eV"][iv]
    bw_d = d["fermi_eV"] - d["gamma_eigenvalues_eV"][iv]
    rows.append({
        "name": name,
        "smearing": c["smearing"], "width_ry": c["width_ry"], "nk": c["nk"],
        "E_qe_Ry": q["total_energy_Ry"], "E_dftk_Ry": d["total_energy_Ry"],
        "dE_meV": (d["total_energy_Ry"] - q["total_energy_Ry"]) * RY_TO_EV * 1000,
        "bw_qe_eV": bw_q, "bw_dftk_eV": bw_d,
        "dbw_meV": (bw_d - bw_q) * 1000,
    })
by_name = {r["name"]: r for r in rows}

print(f"DFTK {dftk_data['dftk_version']} vs QE, aluminum, {len(rows)} cases")
print(f"{'case':22s}{'E_QE (Ry)':>16s}{'dE (meV)':>11s}{'BW_QE (eV)':>12s}{'dBW (meV)':>11s}")
for r in rows:
    print(f"{r['name']:22s}{r['E_qe_Ry']:16.8f}{r['dE_meV']:+11.4f}"
          f"{r['bw_qe_eV']:12.4f}{r['dbw_meV']:+11.3f}")

max_dE = max(abs(r["dE_meV"]) for r in rows)
max_dbw = max(abs(r["dbw_meV"]) for r in rows)
print(f"Max |dE| = {max_dE:.4f} meV/atom, max |d bandwidth| = {max_dbw:.3f} meV")

# ---------------- Figure ----------------

colors = {
    "gaussian": "#2ca02c",
    "mv": "#1f77b4",
    "fd": "#ff7f0e",
}

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# Panel 1: free energy vs smearing width
mv_rows = [
    by_name[n]
    for n in groups["smearing_scan"]
    if by_name[n]["smearing"] == "mv"
]

ref = min(mv_rows, key=lambda r: r["width_ry"])["E_qe_Ry"]

types = sorted({
    by_name[n]["smearing"]
    for n in groups["smearing_scan"]
})

for sm in types:
    rs = sorted(
        (
            by_name[n]
            for n in groups["smearing_scan"]
            if by_name[n]["smearing"] == sm
        ),
        key=lambda r: r["width_ry"],
    )

    w = [r["width_ry"] for r in rs]

    ax[0].plot(
        w,
        [(r["E_qe_Ry"] - ref) * RY_TO_EV * 1000 for r in rs],
        marker="o",
        linestyle="-",
        markerfacecolor="white",
        color=colors[sm],
        label=f"{sm} QE",
    )

    ax[0].plot(
        w,
        [(r["E_dftk_Ry"] - ref) * RY_TO_EV * 1000 for r in rs],
        marker="x",
        linestyle="--",
        color=colors[sm],
        label=f"{sm} DFTK",
    )

ax[0].set_xlabel("Smearing width (Ry)")
ax[0].set_ylabel("F - F$_{ref}$ (meV/atom)")
ax[0].set_title("Free energy vs smearing")

ax[0].legend(
    frameon=False,
    loc="best",
    fontsize=7,
)

# Panel 2: k-grid convergence
ks = sorted(
    (by_name[n] for n in groups["k_scan"]),
    key=lambda r: r["nk"],
)

ref_k = ks[-1]["E_qe_Ry"]

ax[1].plot(
    [r["nk"] for r in ks],
    [(r["E_qe_Ry"] - ref_k) * RY_TO_EV * 1000 for r in ks],
    marker="o",
    linestyle="-",
    markerfacecolor="white",
    color="#1f77b4",
    label="QE",
)

ax[1].plot(
    [r["nk"] for r in ks],
    [(r["E_dftk_Ry"] - ref_k) * RY_TO_EV * 1000 for r in ks],
    marker="x",
    linestyle="--",
    color="#d62728",
    label="DFTK",
)

ax[1].set_xlabel("k-grid (n x n x n)")
ax[1].set_ylabel("F - F(24x24x24, QE) (meV/atom)")
ax[1].set_title("k-point convergence (MV, 0.02 Ry)")

ax[1].legend(
    frameon=False,
    loc="best",
)

# Panel 3: difference between codes, all cases
ax[2].bar(
    range(len(rows)),
    [r["dE_meV"] for r in rows],
    color="#1f77b4",
)

ax[2].axhline(
    TOL_E_MEV,
    color="gray",
    linestyle=":",
    alpha=0.8,
)

ax[2].axhline(
    -TOL_E_MEV,
    color="gray",
    linestyle=":",
    alpha=0.8,
)

ax[2].set_xticks(range(len(rows)))
ax[2].set_xticklabels(
    [r["name"] for r in rows],
    rotation=90,
    fontsize=6,
)

ax[2].set_ylabel("E$_{DFTK}$ - E$_{QE}$ (meV/atom)")
ax[2].set_title("Difference between codes")

fig.tight_layout()

(ROOT / "results").mkdir(exist_ok=True)

fig.savefig(
    ROOT / "results" / "al_smearing.png",
    dpi=150,
)

summary = {
    "dftk_version": dftk_data["dftk_version"],
    "max_abs_dE_meV_per_atom": max_dE,
    "max_abs_d_bandwidth_meV": max_dbw,
    "rows": rows,
}

(ROOT / "results" / "al_comparison.json").write_text(
    json.dumps(summary, indent=2)
)

ok = max_dE < TOL_E_MEV and max_dbw < TOL_BW_MEV

if not ok:
    for r in rows:
        if (
            abs(r["dE_meV"]) >= TOL_E_MEV
            or abs(r["dbw_meV"]) >= TOL_BW_MEV
        ):
            print("  outside tolerance:", r["name"])

print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
