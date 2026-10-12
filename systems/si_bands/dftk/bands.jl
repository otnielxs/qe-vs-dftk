using DFTK
using JSON3
using LinearAlgebra

const HA_TO_EV = 27.211386245988

root = normpath(joinpath(@__DIR__, "..", "..", ".."))
sysdir = normpath(joinpath(@__DIR__, ".."))
params = JSON3.read(read(joinpath(sysdir, "params.json"), String))
kpath = JSON3.read(read(joinpath(sysdir, "kpath.json"), String))

a = Float64(params.a_bohr)
n_bands = Int(params.n_bands)
nk = Int(params.nk_scf)

psp_file = joinpath(root, "pseudos", String(params.pseudo))
Si = ElementPsp(:Si, DFTK.load_psp(psp_file))
lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
atoms = [Si, Si]
positions = [zeros(3), ones(3) / 4]

model = model_DFT(lattice, atoms, positions; functionals=PBE())
basis = PlaneWaveBasis(model; Ecut=params.ecut_ry / 2, kgrid=[nk, nk, nk])
scfres = self_consistent_field(basis; tol=1e-10)

# Cartesian k (units of 2π/a) -> reduced coordinates of THIS lattice:
# k_red = A' * k_cart / (2π)  with  k_cart = 2π/a * κ   =>   k_red = A' * (κ / a)
kcoords = [lattice' * (collect(Float64, κ) ./ a) for κ in kpath.kpoints_tpiba]

# The compute_bands interface changed between DFTK versions, so try two forms
bands = try
    compute_bands(scfres.basis, ExplicitKpoints(kcoords); n_bands=n_bands, ρ=scfres.ρ)
catch err
    err isa MethodError || rethrow()
    compute_bands(scfres, ExplicitKpoints(kcoords); n_bands=n_bands)
end

eigenvalues_eV = [Float64.(e[1:n_bands]) .* HA_TO_EV for e in bands.eigenvalues]

results = Dict(
    "code" => "DFTK",
    "dftk_version" => string(pkgversion(DFTK)),
    "julia_version" => string(VERSION),
    "eigenvalues_eV" => eigenvalues_eV,
)

open(joinpath(@__DIR__, "bands_results.json"), "w") do io
    JSON3.pretty(io, results)
end
println("Wrote bands_results.json with $(length(eigenvalues_eV)) k-points")
