import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qevsdftk.parse_qe import parse_qe_bands

ROOT = Path(__file__).resolve().parents[1]
SYS = ROOT / "systems" / "si_bands"

TOL_BANDS_MEV = 1.0
TOL_GAP_MEV = 1.0
VBM = 3  # 0-based index of the highest occupied band (8 electrons -> 4 bands)

params = json.loads((SYS / "params.json").read_text())
kpath = json.loads((SYS / "kpath.json").read_text())
nb = params["n_bands"]

qe = np.array(parse_qe_bands(SYS / "qe" / "bands.out", nb))
dftk_data = json.loads((SYS / "dftk" / "bands_results.json").read_text())
dftk = np.array(dftk_data["eigenvalues_eV"])
assert qe.shape == dftk.shape, f"Shape mismatch: QE {qe.shape}, DFTK {dftk.shape}"


def analyse(e):
    vbm = e[:, VBM].max()
    cbm = e[:, VBM + 1].min()
    return {
        "vbm_eV": float(vbm),
        "cbm_eV": float(cbm),
        "indirect_gap_eV": float(cbm - vbm),
        "direct_gap_eV": float((e[:, VBM + 1] - e[:, VBM]).min()),
        "k_index_vbm": int(e[:, VBM].argmax()),
        "k_index_cbm": int(e[:, VBM + 1].argmin()),
    }


a_qe, a_dftk = analyse(qe), analyse(dftk)

# Align each code to its own VBM (absolute eigenvalues differ between codes)
rel_qe = qe - a_qe["vbm_eV"]
rel_dftk = dftk - a_dftk["vbm_eV"]
diff_meV = (rel_dftk - rel_qe) * 1000.0
max_diff = float(np.abs(diff_meV).max())
d_indirect = (a_dftk["indirect_gap_eV"] - a_qe["indirect_gap_eV"]) * 1000.0
d_direct = (a_dftk["direct_gap_eV"] - a_qe["direct_gap_eV"]) * 1000.0

print(f"DFTK {dftk_data['dftk_version']} vs QE: {qe.shape[0]} k-points, {nb} bands")
print(f"{'':8s}{'indirect gap (eV)':>20s}{'direct gap (eV)':>18s}")
for name, r in [("QE", a_qe), ("DFTK", a_dftk)]:
    print(f"{name:8s}{r['indirect_gap_eV']:20.4f}{r['direct_gap_eV']:18.4f}")
print(f"Gap difference (DFTK - QE): indirect {d_indirect:+.3f} meV, direct {d_direct:+.3f} meV")
print(f"Max |band difference| over all k and bands (relative to VBM): {max_diff:.3f} meV")

# Figure
x = np.array(kpath["distance"])
ticks = [x[i] for i in kpath["label_index"]]
labels = ["Γ" if t == "G" else t for t in kpath["labels"]]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
for b in range(nb):
    ax[0].plot(x, rel_qe[:, b], "-", color="tab:blue", lw=1.8, label="QE" if b == 0 else None)
    ax[0].plot(x, rel_dftk[:, b], "--", color="tab:red", lw=1.0, label="DFTK" if b == 0 else None)
    ax[1].plot(x, diff_meV[:, b], lw=0.8)
for a_ in ax:
    for t in ticks:
        a_.axvline(t, color="gray", lw=0.5)
    a_.set_xticks(ticks)
    a_.set_xticklabels(labels)
    a_.set_xlim(x[0], x[-1])
ax[0].axhline(0, color="gray", lw=0.5)
ax[0].set_ylim(-13, 8)
ax[0].set_ylabel("E - E$_{VBM}$ (eV)")
ax[0].set_title("Si band structure (PBE)")
ax[0].legend()
ax[1].set_ylabel("DFTK - QE (meV)")
ax[1].set_title("Band differences (all 8 bands)")
fig.tight_layout()
(ROOT / "results").mkdir(exist_ok=True)
fig.savefig(ROOT / "results" / "si_bands.png", dpi=150)

summary = {
    "qe": a_qe,
    "dftk": a_dftk,
    "delta_indirect_gap_meV": d_indirect,
    "delta_direct_gap_meV": d_direct,
    "max_abs_band_difference_meV": max_diff,
    "n_kpoints": int(qe.shape[0]),
}
(ROOT / "results" / "si_bands_comparison.json").write_text(json.dumps(summary, indent=2))

ok = max_diff < TOL_BANDS_MEV and abs(d_indirect) < TOL_GAP_MEV and abs(d_direct) < TOL_GAP_MEV
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
