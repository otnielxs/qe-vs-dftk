using DFTK
using JSON3
using LinearAlgebra

const HA_TO_EV = 27.211386245988
const HA_TO_RY = 2.0

root = normpath(joinpath(@__DIR__, "..", "..", ".."))
sysdir = normpath(joinpath(@__DIR__, ".."))
params = JSON3.read(read(joinpath(sysdir, "params.json"), String))
cases = JSON3.read(read(joinpath(sysdir, "cases.json"), String))

a = Float64(params.a_bohr)
psp_file = joinpath(root, "pseudos", String(params.pseudo))
Al = ElementPsp(:Al, DFTK.load_psp(psp_file))
lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
atoms = [Al]
positions = [zeros(3)]

smearings = Dict(
    "gaussian" => DFTK.Smearing.Gaussian(),
    "fd" => DFTK.Smearing.FermiDirac(),
    "mv" => DFTK.Smearing.MarzariVanderbilt(),
)

function run_case(c)
    model = model_DFT(lattice, atoms, positions;
                      functionals=PBE(),
                      temperature=Float64(c.width_ry) / 2,   # QE degauss (Ry) -> DFTK kT (Ha)
                      smearing=smearings[String(c.smearing)])
    nk = Int(c.nk)
    basis = PlaneWaveBasis(model; Ecut=params.ecut_ry / 2, kgrid=[nk, nk, nk])
    scfres = self_consistent_field(basis; tol=1e-10, callback=identity)
    ik = findfirst(k -> norm(k.coordinate) < 1e-8, basis.kpoints)
    return Dict(
        "name" => String(c.name),
        "total_energy_Ry" => scfres.energies.total * HA_TO_RY,
        "fermi_eV" => scfres.εF * HA_TO_EV,
        "gamma_eigenvalues_eV" => Float64.(scfres.eigenvalues[ik][1:min(8, end)]) .* HA_TO_EV,
        "n_scf_iterations" => scfres.n_iter,
    )
end

results = []
for c in cases.cases
    println("Running $(c.name) ...")
    push!(results, run_case(c))
end

out = Dict(
    "code" => "DFTK",
    "dftk_version" => string(pkgversion(DFTK)),
    "julia_version" => string(VERSION),
    "results" => results,
)
open(joinpath(@__DIR__, "al_results.json"), "w") do io
    JSON3.pretty(io, out)
end
println("Wrote al_results.json with $(length(results)) cases")
