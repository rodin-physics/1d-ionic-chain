include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
d = 60          # Number of time steps in the fastest chain mode
τmax = 300      # Simulation time
ωmax = 10       # Maximum frequency
lmax = 300      # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, 1 / (d * ωmax); batch_size = 100)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end

## Precompute the thermal trajectories
τmax = 1000                     # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_modes = 10000                 # Number of chain masses for simulating ρ0

qa_s = 2 * pi .* (1:n_modes) / n_modes
ωs = ω.(ωmax, qa_s ./ 2)
lmax = 300
ε = reduce(hcat, [exp.(1im * 2 * π / n_modes .* (1:n_modes) * g) for g = 1:lmax])

# Range of temperatures
ωTs = [0, 1, 2, 5, 10, 25, 50, 100, 250]

Random.seed!(150)
for ωT in ωTs
    println("ωT is ", ωT)
    ζs = ζq.(ωs, ωT)
    ϕs = 2 * π * rand(n_modes)
    if (!isfile("precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_l$(lmax).jld2"))
        res = zeros(lmax, n_pts)
        pr = Progress(n_pts)

        Threads.@threads for n = 1:n_pts
            res[:, n] = ρH(n * δ, ζs, ϕs, ωs, ε) |> real
            next!(pr)
        end

        traj = ThermalTrajectory(ωmax, δ, res, ωT)
        save_object(
            "precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_l$(lmax).jld2",
            traj,
        )
    end

end
