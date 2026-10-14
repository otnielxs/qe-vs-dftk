import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qevsdftk.eos import RY_TO_EV, eos_curve, fit_eos
from qevsdftk.parse_qe import parse_qe_energy
from qevsdftk.plotstyle import apply_style

apply_style()


ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "si_eos"
N_ATOMS = 2

TOL_A0_ANG = 0.001
TOL_B0_GPA = 0.1

params = json.loads((SYS / "params.json").read_text())
a = np.array(params["a_values_bohr"])
V = a**3 / 4.0

e_qe = np.array(
    [parse_qe_energy(SYS / "qe" / f"a_{i:02d}.out") for i in range(1, len(a) + 1)]
)
dftk = json.loads((SYS / "dftk" / "eos_results.json").read_text())
assert np.allclose(dftk["a_bohr"], a), "DFTK and QE used different lattice constants"
e_dftk = np.array(dftk["energy_Ry"])

fit_qe = fit_eos(V, e_qe, N_ATOMS)
fit_dftk = fit_eos(V, e_dftk, N_ATOMS)

dE_mev_atom = (e_dftk - e_qe) * RY_TO_EV * 1000 / N_ATOMS
d_a0 = fit_dftk["a0_angstrom"] - fit_qe["a0_angstrom"]
d_b0 = fit_dftk["B0_GPa"] - fit_qe["B0_GPa"]

print(f"DFTK {dftk['dftk_version']} vs QE (E-V curve, {len(a)} points)")
print(f"{'':10s}{'a0 (A)':>12s}{'B0 (GPa)':>12s}{'B0prime':>10s}{'rms (meV/atom)':>16s}")
for name, f in [("QE", fit_qe), ("DFTK", fit_dftk)]:
    print(f"{name:10s}{f['a0_angstrom']:12.5f}{f['B0_GPa']:12.3f}"
          f"{f['B0_prime']:10.3f}{f['rms_residual_meV_per_atom']:16.5f}")
print(f"Difference (DFTK - QE): da0 = {d_a0:+.6f} A, dB0 = {d_b0:+.4f} GPa")
print(f"Max |E_DFTK - E_QE| over all volumes: {np.abs(dE_mev_atom).max():.5f} meV/atom")

# Figure
Vfine = np.linspace(V.min(), V.max(), 300)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))

colors = {
    "QE": "#1f77b4",
    "DFTK": "#d62728",
}

for name, e, fit, marker in [
    ("QE", e_qe, fit_qe, "o"),
    ("DFTK", e_dftk, fit_dftk, "x"),
]:
    ref = fit["E0_Ry"]

    ax[0].plot(
        V,
        (e - ref) * RY_TO_EV * 1000 / N_ATOMS,
        marker=marker,
        linestyle="None",
        color=colors[name],
    )

    ax[0].plot(
        Vfine,
        (eos_curve(fit, Vfine) - ref) * RY_TO_EV * 1000 / N_ATOMS,
        color=colors[name],
        lw=1.5,
        label=f"{name} BM fit",
    )

ax[0].set_xlabel("Volume per primitive cell (Bohr$^3$)")
ax[0].set_ylabel("E - E$_0$ (meV/atom)")
ax[0].set_title("Si equation of state")
ax[0].legend(
    frameon=False,
    loc="best",
)

ax[1].plot(
    a,
    dE_mev_atom,
    marker="o",
    linestyle="-",
    markerfacecolor="white",
    color="#1f77b4",
)

ax[1].set_xlabel("Lattice constant (Bohr)")
ax[1].set_ylabel("E$_{DFTK}$ - E$_{QE}$ (meV/atom)")
ax[1].set_title("Energy difference at each volume")

fig.tight_layout()

(ROOT / "results").mkdir(exist_ok=True)
fig.savefig(ROOT / "results" / "si_eos.png", dpi=150)

summary = {
    "a_bohr": a.tolist(),
    "qe": {**fit_qe, "energy_Ry": e_qe.tolist()},
    "dftk": {
        **fit_dftk,
        "energy_Ry": e_dftk.tolist(),
        "dftk_version": dftk["dftk_version"],
    },
    "delta_a0_angstrom": d_a0,
    "delta_B0_GPa": d_b0,
    "max_abs_dE_meV_per_atom": float(np.abs(dE_mev_atom).max()),
}

(ROOT / "results" / "si_eos_comparison.json").write_text(
    json.dumps(summary, indent=2)
)

ok = abs(d_a0) < TOL_A0_ANG and abs(d_b0) < TOL_B0_GPA

print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
