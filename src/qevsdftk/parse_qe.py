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
def parse_qe_energy(path):
    """Total energy (Ry) only; works for runs without verbosity='high'."""
    text = Path(path).read_text()
    energies = re.findall(rf"^!\s+total energy\s+=\s+({FLOAT})\s+Ry", text, re.M)
    if not energies:
        raise ValueError(f"No converged total energy found in {path}")
    return float(energies[-1])
def parse_qe_bands(path, n_bands):
    """Eigenvalues (eV) of a QE 'bands' run: list of nk lists with n_bands floats each.

    Requires verbosity='high' in the QE input. k-points are returned in the
    same order as in the input file.
    """
    text = Path(path).read_text()
    section = text.split("End of band structure calculation")[-1]
    parts = re.split(r"k =[^\n]*?\(\s*\d+ PWs\)\s+bands \(ev\):", section)[1:]
    if not parts:
        raise ValueError("No band energies found (is verbosity='high' set?)")

    bands = []
    for p in parts:
        values = [float(x) for x in re.findall(FLOAT, p)][:n_bands]
        if len(values) != n_bands:
            raise ValueError("Fewer eigenvalues than expected in a k-point block")
        bands.append(values)
    return bands
def _parse_clock(s):
    """Convert QE clock strings like '0.94s', '1m 5.20s' or '1h 2m 3.00s' to seconds."""
    units = {"h": 3600.0, "m": 60.0, "s": 1.0}
    return sum(float(v) * units[u] for v, u in re.findall(r"(\d+(?:\.\d+)?)\s*([hms])", s))


def parse_qe_wall_time(path):
    """Total wall time (s) from the final PWSCF timing line of a QE output."""
    text = Path(path).read_text()
    found = re.findall(r"^\s*PWSCF\s*:\s*(.*?)\s+CPU\s+(.*?)\s+WALL", text, re.M)
    if not found:
        raise ValueError(f"No PWSCF timing line found in {path}")
    return _parse_clock(found[-1][1])
