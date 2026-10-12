from pathlib import Path

from qevsdftk.parse_qe import parse_qe

FIXTURE = Path(__file__).parent / "fixtures" / "al_scf.out"


def test_parse_qe_metal():
    r = parse_qe(FIXTURE)
    assert r["total_energy_Ry"] < 0
    assert r["fermi_eV"] is not None
    # Occupied bandwidth of aluminum should be of the order of 10 eV
    n_semicore = 4  # this Al pseudopotential has 2s2p semicore states (11 valence electrons)
    bandwidth = r["fermi_eV"] - r["gamma_eigenvalues_eV"][n_semicore]
    assert 5.0 < bandwidth < 15.0
