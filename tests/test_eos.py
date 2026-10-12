import numpy as np

from qevsdftk.eos import RY_BOHR3_TO_GPA, birch_murnaghan, fit_eos


def test_fit_recovers_known_parameters():
    E0, V0, B0, B0p = -15.7, 270.0, 0.006, 4.2
    V = np.linspace(0.88 * V0, 1.12 * V0, 9)
    E = birch_murnaghan(V, E0, V0, B0, B0p)

    r = fit_eos(V, E)

    assert abs(r["E0_Ry"] - E0) < 1e-8
    assert abs(r["V0_bohr3"] - V0) < 1e-4
    assert abs(r["B0_GPa"] - B0 * RY_BOHR3_TO_GPA) < 1e-3
    assert abs(r["B0_prime"] - B0p) < 1e-3
    assert abs(r["a0_bohr"] - (4.0 * V0) ** (1.0 / 3.0)) < 1e-6
