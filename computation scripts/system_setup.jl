include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
d = 60          # Number of time steps in the fastest chain mode
τmin = 0
τmax = 300      # Simulation time
ωmax = 10       # Maximum frequency
lmax = 300      # Number of chain atoms tracked



if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(τmin, τmax, 1 / (d * ωmax), 0:lmax, ωmax)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
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

# println("Calculating thermal trajectories")
# ## Precompute the thermal trajectories
# τmax = 5                         # Simulation time
# δ = (1 / ωmax) / d              # Time step
# n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
# n_modes = 100000                  # Number of chain masses for simulating ρ0

# qa_s = 2 * pi .* (1:n_modes) / n_modes
# ωs = ω.(ωmax, qa_s ./ 2)

# Random.seed!(150)
# ϕs = 2 * π * rand(n_modes)

# # Range of temperatures
# ωTs = [0.0, 2.0, 5.0, 10.0]

# for ωT in ωTs
#     println("ωT is ", ωT)
#     ζs = ζq.(ωs, ωT)

#     if (
#         !isfile(
#             "precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_nmodes$(n_modes).jld2",
#         )
#     )
#         # Populate each row of matrix
#         gs = collect(1:lmax)
#         full_res = @showprogress pmap(n -> real(ρH(n, δ, ζs, ϕs, ωs, gs)), 1:n_pts)
#         full_res = reduce(hcat, full_res)

#         res = ThermalTrajectory(ωmax, δ, full_res, ωT)
#         save_object(
#             "precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax)_nmodes$(n_modes).jld2",
#             res,
#         )
#     end

# end
