from pathlib import Path

from qevsdftk.parse_qe import parse_qe

FIXTURE = Path(__file__).parent / "fixtures" / "si_scf.out"


def test_parse_qe_basic():
    r = parse_qe(FIXTURE)
    assert r["total_energy_Ry"] < 0
    assert len(r["gamma_eigenvalues_eV"]) >= 4
    assert r["n_scf_iterations"] is not None
