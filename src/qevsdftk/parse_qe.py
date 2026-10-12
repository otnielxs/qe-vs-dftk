import re
from pathlib import Path

FLOAT = r"-?\d+\.\d+"


def parse_qe(path):
    text = Path(path).read_text()

    energies = re.findall(rf"^!\s+total energy\s+=\s+({FLOAT})\s+Ry", text, re.M)
    if not energies:
        raise ValueError(f"No converged total energy found in {path}")

    # Only look at the final SCF section for eigenvalues
    final = text.split("End of self-consistent calculation")[-1]
    block = re.search(
        r"k =\s*-?0\.0000\s+-?0\.0000\s+-?0\.0000.*?bands \(ev\):(.*?)"
        r"(?:occupation numbers|\n\s*k =|the Fermi energy)",
        final,
        re.S,
    )
    if block is None:
        raise ValueError("Gamma-point eigenvalues not found (is verbosity='high'?)")
    gamma = [float(x) for x in re.findall(FLOAT, block.group(1))]

    fermi = re.search(rf"the Fermi energy is\s+({FLOAT})\s+ev", text)
    version = re.search(r"Program PWSCF v\.([\w.\-]+)", text)
    iters = re.search(r"convergence has been achieved in\s+(\d+)\s+iterations", text)

    return {
        "code": "QE",
        "qe_version": version.group(1) if version else None,
        "total_energy_Ry": float(energies[-1]),
        "fermi_eV": float(fermi.group(1)) if fermi else None,
        "gamma_eigenvalues_eV": gamma,
        "n_scf_iterations": int(iters.group(1)) if iters else None,
    }
