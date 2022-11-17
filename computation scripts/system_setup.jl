using Distributed
if nprocs() < 6
    addprocs(5)
end
@everywhere using Random
@everywhere include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
d = 60       # Number of time steps in the fastest chain mode
τmax = 5    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 50000   # Number of chain atoms tracked

function precompute(ωmax, τmax, lmax, d)
    δ = (1 / ωmax) / d
    n_pts = floor(τmax / δ) |> Int

    # Read in files that have the correct ωmax and d
    filenames = filter(
        x -> first(x) !== '.' && occursin("_ωmax$(ωmax)_", x) && occursin("_d$(d)_", x),
        readdir(joinpath(pwd(), "precomputed/systems/")),
    )

    # No precomputation files exist
    if isempty(filenames)
        res = mkChainSystem(ωmax, τmax, lmax, d, Matrix[])
        return save_object(
            "precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax)_τmax$(τmax).jld2",
            res,
        )
    end

    # Check if precomputation is necessary
    size_elems =
        [load_object(joinpath("precomputed/systems/", ii)).Γ |> size for ii in filenames]

    if any(x -> (lmax <= x[1] && n_pts <= x[2]) == true, size_elems)
        error("A file containing a sufficient amount of precomputed values already exists")
    else
        # Return existing matrix with most number of elements already computed
        diff_elems =
            map(x -> (lmax * n_pts) - min(lmax, x[1]) * min(n_pts, x[2]), size_elems)

        Γ_prev =
            load_object(joinpath("precomputed/systems/", filenames[argmin(diff_elems)])).Γ

        res = mkChainSystem(ωmax, τmax, lmax, d, Γ_prev)
        save_object(
            "precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax)_τmax$(τmax).jld2",
            res,
        )
    end
end

# precompute(ωmax, τmax, lmax, d)

## Prepare the ultrafine ChainSystem's by calculating the recoil term
# d = 6000     # Number of time steps in the fastest chain mode
# τmax = 20    # Simulation time
# ωmax = 10    # Maximum frequency
# lmax = 20    # Number of chain atoms tracked
#
# if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
#     res = mkChainSystem(ωmax, τmax, lmax, d)
#     save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
# end

# precompute(ωmax, τmax, lmax, d)

println("Calculating thermal trajectories")
## Precompute the thermal trajectories
τmax = 5                         # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_modes = 100000                  # Number of chain masses for simulating ρ0

qa_s = 2 * pi .* (1:n_modes) / n_modes
ωs = ω.(ωmax, qa_s ./ 2)

Random.seed!(150)
ϕs = 2 * π * rand(n_modes)

# Range of temperatures
ωTs = [0.0, 2.0, 5.0, 10.0]

for ωT in ωTs
    println("ωT is ", ωT)
    ζs = ζq.(ωs, ωT)

    if (
        !isfile(
            "precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_nmodes$(n_modes).jld2",
        )
    )
        # Populate each row of matrix
        gs = collect(1:lmax)
        full_res = @showprogress pmap(n -> real(ρH(n, δ, ζs, ϕs, ωs, gs)), 1:n_pts)
        full_res = reduce(hcat, full_res)

        res = ThermalTrajectory(ωmax, δ, full_res, ωT)
        save_object(
            "precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_nmodes$(n_modes).jld2",
            res,
        )
    end

end
