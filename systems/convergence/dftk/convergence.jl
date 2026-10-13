using DFTK
using JSON3

const HA_TO_RY = 2.0

root = normpath(joinpath(@__DIR__, "..", "..", ".."))
sysdir = normpath(joinpath(@__DIR__, ".."))
params = JSON3.read(read(joinpath(sysdir, "params.json"), String))
cases = JSON3.read(read(joinpath(sysdir, "cases.json"), String))

smearings = Dict(
    "gaussian" => DFTK.Smearing.Gaussian(),
    "fd" => DFTK.Smearing.FermiDirac(),
    "mv" => DFTK.Smearing.MarzariVanderbilt(),
)

function build_model(system_name)
    s = params.systems[Symbol(system_name)]
    psp_file = joinpath(root, "pseudos", String(params.pseudos[Symbol(system_name)]))
    el = ElementPsp(Symbol(String(s.element)), DFTK.load_psp(psp_file))
    a = Float64(s.a_bohr)
    lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
    positions = [Float64.(collect(p)) for p in s.positions_crystal]
    atoms = fill(el, length(positions))
    occ = s.occupations
    if String(occ.kind) == "smearing"
        return model_DFT(lattice, atoms, positions; functionals=PBE(),
                         temperature=Float64(occ.width_ry) / 2,
                         smearing=smearings[String(occ.smearing)])
    end
    return model_DFT(lattice, atoms, positions; functionals=PBE())
end

results = []
for c in cases.cases
    println("Running $(c.name) ...")
    model = build_model(String(c.system))
    nk = Int(c.nk)
    basis = PlaneWaveBasis(model; Ecut=Float64(c.ecut_ry) / 2, kgrid=[nk, nk, nk])
    scfres = self_consistent_field(basis; tol=1e-10, callback=identity)
    push!(results, Dict(
        "name" => String(c.name),
        "total_energy_Ry" => scfres.energies.total * HA_TO_RY,
        "n_scf_iterations" => scfres.n_iter,
    ))
end

out = Dict(
    "code" => "DFTK",
    "dftk_version" => string(pkgversion(DFTK)),
    "julia_version" => string(VERSION),
    "results" => results,
)
open(joinpath(@__DIR__, "convergence_results.json"), "w") do io
    JSON3.pretty(io, out)
end
println("Wrote convergence_results.json with $(length(results)) cases")
