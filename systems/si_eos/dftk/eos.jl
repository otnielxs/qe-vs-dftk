using DFTK
using JSON3

const HA_TO_RY = 2.0

root = normpath(joinpath(@__DIR__, "..", "..", ".."))
params = JSON3.read(read(joinpath(@__DIR__, "..", "params.json"), String))

psp_file = joinpath(root, "pseudos", String(params.pseudo))
Si = ElementPsp(:Si, DFTK.load_psp(psp_file))
atoms = [Si, Si]
positions = [zeros(3), ones(3) / 4]

Ecut_Ha = params.ecut_ry / 2     # 60 Ry = 30 Ha
nk = Int(params.nk)
a_values = Float64.(params.a_values_bohr)

function total_energy_Ry(a)
    lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
    model = model_DFT(lattice, atoms, positions; functionals=PBE())
    basis = PlaneWaveBasis(model; Ecut=Ecut_Ha, kgrid=[nk, nk, nk])
    scfres = self_consistent_field(basis; tol=1e-10, callback=identity)
    return scfres.energies.total * HA_TO_RY
end

energies = Float64[]
for a in a_values
    E = total_energy_Ry(a)
    println("a = $a Bohr   E = $E Ry")
    push!(energies, E)
end

results = Dict(
    "code" => "DFTK",
    "dftk_version" => string(pkgversion(DFTK)),
    "julia_version" => string(VERSION),
    "a_bohr" => a_values,
    "energy_Ry" => energies,
)

open(joinpath(@__DIR__, "eos_results.json"), "w") do io
    JSON3.pretty(io, results)
end
println("Wrote eos_results.json")
