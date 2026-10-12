using DFTK
using LinearAlgebra
using JSON3

const HA_TO_EV = 27.211386245988
const HA_TO_RY = 2.0

root = normpath(joinpath(@__DIR__, "..", "..", ".."))
psp_file = joinpath(root, "pseudos", "Si_ONCV_PBE-1.2.upf")

a = 10.26  # Bohr
lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
Si = ElementPsp(:Si, DFTK.load_psp(psp_file))
atoms = [Si, Si]
positions = [zeros(3), ones(3) / 4]

model = model_DFT(lattice, atoms, positions;
                  functionals=PBE(),
                  temperature=0.1,
                  smearing=DFTK.Smearing.MarzariVanderbilt())
basis = PlaneWaveBasis(model; Ecut=30, kgrid=[6, 6, 6])

# First run includes Julia compilation time; the second run is the timed one
scfres = self_consistent_field(basis; tol=1e-10)
t = @elapsed scfres = self_consistent_field(basis; tol=1e-10, callback=identity)

# Eigenvalues at Gamma (find the k-point by its coordinates, do not assume index 1)
ik = findfirst(k -> norm(k.coordinate) < 1e-8, basis.kpoints)
eig_gamma_eV = Float64.(scfres.eigenvalues[ik]) .* HA_TO_EV

results = Dict(
    "code" => "DFTK",
    "dftk_version" => string(pkgversion(DFTK)),
    "julia_version" => string(VERSION),
    "total_energy_Ry" => scfres.energies.total * HA_TO_RY,
    "fermi_eV" => scfres.εF * HA_TO_EV,
    "gamma_eigenvalues_eV" => eig_gamma_eV,
    "n_scf_iterations" => scfres.n_iter,
    "scf_time_s" => t,
)

open(joinpath(@__DIR__, "results.json"), "w") do io
    JSON3.pretty(io, results)
end
println("Wrote results.json")
