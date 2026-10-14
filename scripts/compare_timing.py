import json
import os
import re
import statistics as st
import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from qevsdftk.parse_qe import parse_qe_energy, parse_qe_wall_time
from qevsdftk.plotstyle import apply_style

apply_style()

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "timing"
RY_TO_EV = 13.605693122994
TOL_E_MEV = 0.1

params = json.loads((SYS / "params.json").read_text())
cores = params["cores"]
repeats = params["repeats"]


def cpu_model():
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return "unknown"


table, max_dE, qe_version, dftk_version, julia_version = {}, 0.0, None, None, None
for w, wp in params["workloads"].items():
    n_atoms = len(wp["positions_crystal"])
    for n in cores:
        qe_times, qe_energy, qe_iters = [], None, None
        for r in range(1, repeats + 1):
            f = SYS / "qe" / f"{w}_np{n}_r{r}.out"
            text = f.read_text()
            qe_times.append(parse_qe_wall_time(f))
            qe_energy = parse_qe_energy(f)
            m = re.search(r"convergence has been achieved in\s+(\d+)\s+iterations", text)
            qe_iters = int(m.group(1)) if m else None
            v = re.search(r"Program PWSCF v\.([\w.\-]+)", text)
            qe_version = v.group(1) if v else qe_version
        d = json.loads((SYS / "dftk" / f"timing_{w}_t{n}.json").read_text())
        dftk_version, julia_version = d["dftk_version"], d["julia_version"]
        dE = (d["energies_Ry"][-1] - qe_energy) * RY_TO_EV * 1000 / n_atoms
        max_dE = max(max_dE, abs(dE))
        table[(w, n)] = {
            "qe_times": qe_times, "dftk_times": d["times_s"],
            "qe_med": st.median(qe_times), "dftk_med": st.median(d["times_s"]),
            "qe_iters": qe_iters, "dftk_iters": d["n_scf_iterations"][-1],
            "dE_meV": dE,
        }

print(f"QE {qe_version} vs DFTK {dftk_version} (Julia {julia_version})")
print(f"CPU: {cpu_model()}, logical cores: {os.cpu_count()}")
for w in params["workloads"]:
    t1q, t1d = table[(w, cores[0])]["qe_med"], table[(w, cores[0])]["dftk_med"]
    print(f"\nWorkload {w}")
    print(f"{'cores':>6}{'QE (s)':>10}{'DFTK (s)':>10}{'DFTK/QE':>9}"
          f"{'QE speedup':>12}{'DFTK speedup':>14}{'iters QE/DFTK':>15}")
    for n in cores:
        r = table[(w, n)]
        print(f"{n:6d}{r['qe_med']:10.2f}{r['dftk_med']:10.2f}{r['dftk_med'] / r['qe_med']:9.2f}"
              f"{t1q / r['qe_med']:12.2f}{t1d / r['dftk_med']:14.2f}"
              f"{str(r['qe_iters']) + '/' + str(r['dftk_iters']):>15}")
print(f"\nMax |E_DFTK - E_QE| over all runs: {max_dE:.4f} meV/atom")

# ---------------- Figure ----------------
wl = list(params["workloads"])

fig, axes = plt.subplots(
    len(wl),
    2,
    figsize=(11, 4.3 * len(wl)),
    squeeze=False,
)
fig, axes = plt.subplots(
    len(wl),
    2,
    figsize=(11, 4.3 * len(wl)),
    squeeze=False,
)

for i, w in enumerate(wl):

    for label, key, marker, linestyle, color in [
        ("QE", "qe", "o", "-", "#1f77b4"),
        ("DFTK (MPI)", "dftk", "s", "--", "#d62728"),
    ]:

        med = [table[(w, n)][f"{key}_med"] for n in cores]

        lo = [
            m - min(table[(w, n)][f"{key}_times"])
            for m, n in zip(med, cores)
        ]

        hi = [
            max(table[(w, n)][f"{key}_times"]) - m
            for m, n in zip(med, cores)
        ]

        axes[i][0].errorbar(
            cores,
            med,
            yerr=[lo, hi],
            marker=marker,
            linestyle=linestyle,
            color=color,
            markerfacecolor="white" if marker == "o" else color,
            label=label,
            capsize=3,
        )

        axes[i][1].plot(
            cores,
            [med[0] / m for m in med],
            marker=marker,
            linestyle=linestyle,
            color=color,
            markerfacecolor="white" if marker == "o" else color,
            label=label,
        )

    axes[i][1].plot(
        cores,
        [n / cores[0] for n in cores],
        linestyle=":",
        color="gray",
        alpha=0.8,
        label="ideal",
    )

    axes[i][0].set_xscale(
        "log",
        base=2,
    )

    axes[i][0].set_yscale("log")

    axes[i][0].set_xticks(cores)
    axes[i][0].set_xticklabels(cores)

    axes[i][0].set_xlabel(
        "Cores (MPI processes / Julia threads)"
    )

    axes[i][0].set_ylabel("Wall time (s)")
    axes[i][0].set_title(f"{w}: time to solution")

    axes[i][0].legend(
        frameon=False,
        loc="best",
    )

    axes[i][1].set_xlabel("Cores")
    axes[i][1].set_ylabel("Speedup vs 1 core")
    axes[i][1].set_title(f"{w}: parallel speedup")

    axes[i][1].legend(
        frameon=False,
        loc="best",
    )

fig.tight_layout()

(ROOT / "results").mkdir(exist_ok=True)

fig.savefig(
    ROOT / "results" / "timing.png",
    dpi=150,
)

(ROOT / "results" / "timing_comparison.json").write_text(
    json.dumps(
        {
            "date": str(date.today()),
            "cpu": cpu_model(),
            "logical_cores": os.cpu_count(),
            "qe_version": qe_version,
            "dftk_version": dftk_version,
            "julia_version": julia_version,
            "dftk_parallel_mode": (
                "MPI (1 Julia thread and 1 BLAS thread per rank)"
            ),
            "max_abs_dE_meV_per_atom": max_dE,
            "runs": {
                f"{w}_{n}": v
                for (w, n), v in table.items()
            },
        },
        indent=2,
    )
)

ok = max_dE < TOL_E_MEV

print(
    "PASS (both codes solved the same problem)"
    if ok
    else "FAIL (energies differ)"
)

sys.exit(0 if ok else 1)
