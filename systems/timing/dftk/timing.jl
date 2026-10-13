using DFTK
using JSON3
using LinearAlgebra

root = normpath(joinpath(@__DIR__, "..", "..", ".."))
sysdir = normpath(joinpath(@__DIR__, ".."))
params = JSON3.read(read(joinpath(sysdir, "params.json"), String))

workload = ARGS[1]
w = params.workloads[Symbol(workload)]

smearings = Dict(
    "gaussian" => DFTK.Smearing.Gaussian(),
    "fd" => DFTK.Smearing.FermiDirac(),
    "mv" => DFTK.Smearing.MarzariVanderbilt(),
)

function build_model()
    psp_file = joinpath(root, "pseudos", String(params.pseudos[Symbol(w.pseudo)]))
    el = ElementPsp(Symbol(String(w.element)), DFTK.load_psp(psp_file))
    a = Float64(w.a_bohr)
    if String(w.cell) == "cubic"
        lattice = a * Matrix{Float64}(I, 3, 3)
    else
        lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
    end
    positions = [Float64.(collect(p)) for p in w.positions_crystal]
    atoms = fill(el, length(positions))
    occ = w.occupations
    if String(occ.kind) == "smearing"
        return model_DFT(lattice, atoms, positions; functionals=PBE(),
                         temperature=Float64(occ.width_ry) / 2,
                         smearing=smearings[String(occ.smearing)])
    end
    return model_DFT(lattice, atoms, positions; functionals=PBE())
end

function run_scf(ecut_Ha, nk)
    model = build_model()
    basis = PlaneWaveBasis(model; Ecut=ecut_Ha, kgrid=[nk, nk, nk])
    return self_consistent_field(basis; tol=1e-6, callback=identity)
end

ecut_Ha = Float64(w.ecut_ry) / 2
nk = Int(w.nk)
BLAS.set_num_threads(parse(Int, get(ENV, "DFTK_BLAS_THREADS", "1")))
println("Threads: $(Threads.nthreads()), BLAS threads: $(BLAS.get_num_threads())")
println("Warm-up (compilation, not timed) ...")
run_scf(ecut_Ha / 3, 2)

DFTK.reset_timer!(DFTK.timer)
times, energies, iters = Float64[], Float64[], Int[]
for r in 1:Int(params.repeats)
    GC.gc()
    t = @elapsed scfres = run_scf(ecut_Ha, nk)
    push!(times, t)
    push!(energies, scfres.energies.total * 2.0)
    push!(iters, scfres.n_iter)
    println("repeat $r: $(round(t, digits=2)) s, $(scfres.n_iter) iterations")
end
println(DFTK.timer)
out = Dict(
    "code" => "DFTK",
    "dftk_version" => string(pkgversion(DFTK)),
    "julia_version" => string(VERSION),
    "workload" => workload,
    "threads" => Threads.nthreads(),
    "blas_threads" => BLAS.get_num_threads(),
    "times_s" => times,
    "energies_Ry" => energies,
    "n_scf_iterations" => iters,
)
open(joinpath(@__DIR__, "timing_$(workload)_t$(Threads.nthreads()).json"), "w") do io
    JSON3.pretty(io, out)
end
