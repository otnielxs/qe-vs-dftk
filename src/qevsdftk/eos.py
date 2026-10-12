import numpy as np
from scipy.optimize import curve_fit

RY_BOHR3_TO_GPA = 14710.507848
BOHR_TO_ANG = 0.529177210903
RY_TO_EV = 13.605693122994


def birch_murnaghan(V, E0, V0, B0, B0p):
    """Third-order Birch-Murnaghan E(V). Units: V in Bohr^3, E in Ry, B0 in Ry/Bohr^3."""
    eta = (V0 / V) ** (2.0 / 3.0)
    return E0 + 9.0 * V0 * B0 / 16.0 * (
        (eta - 1.0) ** 3 * B0p + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )


def fit_eos(volumes_bohr3, energies_ry, n_atoms=2):
    V = np.asarray(volumes_bohr3, dtype=float)
    E = np.asarray(energies_ry, dtype=float)

    # Fit shifted energies: better numerical conditioning
    e_shift = E.min()
    i = int(np.argmin(E))
    p0 = [0.0, V[i], 0.006, 4.0]
    popt, _ = curve_fit(birch_murnaghan, V, E - e_shift, p0=p0, maxfev=20000)
    E0, V0, B0, B0p = popt

    resid = (E - e_shift) - birch_murnaghan(V, *popt)
    a0_bohr = (4.0 * V0) ** (1.0 / 3.0)  # fcc primitive cell: V = a^3 / 4
    return {
        "E0_Ry": float(E0 + e_shift),
        "V0_bohr3": float(V0),
        "B0_GPa": float(B0 * RY_BOHR3_TO_GPA),
        "B0_prime": float(B0p),
        "a0_bohr": float(a0_bohr),
        "a0_angstrom": float(a0_bohr * BOHR_TO_ANG),
        "rms_residual_meV_per_atom": float(
            np.sqrt(np.mean(resid**2)) * RY_TO_EV * 1000 / n_atoms
        ),
    }


def eos_curve(fit, V):
    """Evaluate the fitted curve (Ry) at volumes V, from a fit_eos result."""
    return birch_murnaghan(
        V, fit["E0_Ry"], fit["V0_bohr3"], fit["B0_GPa"] / RY_BOHR3_TO_GPA, fit["B0_prime"]
    )
