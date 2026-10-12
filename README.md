# qe-vs-dftk

Reproducible comparison of [Quantum ESPRESSO](https://www.quantum-espresso.org/) (QE)
and [DFTK.jl](https://github.com/JuliaMolSim/DFTK.jl) on identical plane-wave DFT calculations.

The goal is to verify that both codes give the same answers under the same settings,
and to help people who do not have QE installed try simple DFT simulations with DFTK.

## Quick start

Requirements: QE (`pw.x`, MPI), Julia, Python 3.9+.

```bash
git clone https://github.com/otnielxs/qe-vs-dftk.git
cd qe-vs-dftk
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python -m pip install -e ".[dev]"
./scripts/run_si.sh
```

## Stage 1: bulk silicon

| Setting | QE | DFTK |
|---|---|---|
| Functional | PBE | PBE |
| Pseudopotential | ONCV (SG15) | same file |
| Cutoff | 60 Ry | 30 Ha |
| k-grid | 6x6x6, Gamma-centered | 6x6x6 |
| Smearing | Marzari-Vanderbilt, 0.2 Ry | Marzari-Vanderbilt, 0.1 Ha |

Results (QE 7.5 vs DFTK 0.8.0):

| Quantity | Difference (DFTK - QE) |
|---|---|
| Total energy | +0.0046 meV/atom |
| Gamma eigenvalues relative to VBM (8 bands) | all within 0.04 meV |

### Notes

- Absolute Fermi levels differ (about 0.9 eV) because the zero of the electrostatic
  potential is defined differently in each code. Only energy differences are
  physically meaningful, so this repo compares eigenvalues relative to the VBM.
- SCF iteration counts are not comparable: the mixing schemes and convergence
  criteria are different.
- Both codes use the same pseudopotential, so this verifies the implementations,
  not the accuracy of the physics.

## Roadmap

- [x] Stage 1: Si total energy and Gamma eigenvalues
- [ ] Si with fixed occupations (clean band gap)
- [ ] Stage 2: E-V curve (lattice constant, bulk modulus)
- [ ] Stage 3: metal (Al) and smearing
- [ ] Stage 4: convergence and timing
- [ ] Stage 5: MoS2 monolayer

## Versions and credits

QE 7.5, DFTK 0.8.0
Pseudopotential: SG15 ONCV
