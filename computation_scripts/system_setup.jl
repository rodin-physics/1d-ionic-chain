using Random
include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
d = 60       # Number of time steps in the fastest chain mode
τmax = 20    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 100   # Number of chain atoms tracked

function precompute(ωmax, τmax, lmax, d)
    δ = (1 / ωmax) / d
    n_pts = floor(τmax / δ) |> Int

    # Read in files that have the correct ωmax and d
    filenames = filter(x -> first(x) !== '.' && occursin("_ωmax$(ωmax)_", x) && occursin("_d$(d)_", x), readdir(joinpath(pwd(), "precomputed/systems/")))

    # No precomputation files exist
    if isempty(filenames)
        res = mkChainSystem(ωmax, τmax, lmax, d, Matrix[])
        return save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax)_τmax$(τmax).jld2", res)
    end

    # Check if precomputation is necessary
    size_elems = [load_object(joinpath("precomputed/systems/", ii)).Γ |> size for ii in filenames]

    if any(x -> (lmax <= x[1] && n_pts <= x[2]) == true, size_elems)
        error("A file containing a sufficient amount of precomputed values already exists")
    else
        # Return existing matrix with most number of elements already computed
        diff_elems = map(x -> (lmax * n_pts) - min(lmax,x[1])*min(n_pts,x[2]), size_elems)

        Γ_prev = load_object(joinpath("precomputed/systems/", filenames[argmin(diff_elems)])).Γ

        res = mkChainSystem(ωmax, τmax, lmax, d, Γ_prev)
        save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax)_τmax$(τmax).jld2", res)
    end
end

precompute(ωmax, τmax, lmax, d)

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


## Precompute the thermal trajectories
τmax = 1000                     # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_modes = 50000                 # Number of chain masses for simulating ρ0

qa_s = 2 * pi .* (1:n_modes) / n_modes
ωs = ω.(ωmax, qa_s ./ 2)
lmax = 300
ε = reduce(hcat, [exp.(1im * 2 * π / n_modes .* (1:n_modes) * g) for g = 1:lmax])

# Range of temperatures
ωTs = [0, 1, 2, 5, 10, 15, 25, 50, 100, 250]

Random.seed!(150)
for ωT in ωTs
    println("ωT is ", ωT)
    ζs = ζq.(ωs, ωT)
    ϕs = 2 * π * rand(n_modes)
    if (!isfile("precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_l$(lmax)_modes$(n_modes).jld2"))
        res = zeros(lmax, n_pts)
        pr = Progress(n_pts)

        Threads.@threads for n = 1:n_pts
            res[:, n] = ρH(n * δ, ζs, ϕs, ωs, ε) |> real
            next!(pr)
        end

        traj = ThermalTrajectory(ωmax, δ, res, ωT)
        save_object(
            "precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_l$(lmax)_modes$(n_modes).jld2",
            traj,
        )
    end

end