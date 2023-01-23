include("main.jl")

# τmax = 5                        # Simulation time
# δ = (1 / ωmax) / d              # Time step
# n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
# n_modes = 100000                # Number of chain masses for simulating ρ0

# qa_s = 2 * pi .* (1:n_modes) / n_modes
# ωs = ω.(ωmax, qa_s ./ 2)

# n_ζ = length(ζs)
# ε = reduce(hcat, [exp.(1im * 2 * π / n_ζ .* (1:n_ζ) * g) for g = 1:lmax])

# Random.seed!(150)
# ϕs = 2 * π * rand(n_modes)
# ζs = ζq.(ωs, ωT)
# # Range of temperatures
# ωTs = [0.0, 2.0, 5.0, 10.0]

# @time ρH(1000, ζs, ϕs, ωs, ε)
# @time ρH_TEST(1000, ζs, ϕs, ωs, ε)


# function ρH_TEST(τ, ζs, ϕs, ωs, ε)
#     n_ζ = length(ζs)
#     # Get the displacement of each mode ζ
#     f = transpose(exp.(-1im * (2 * π * τ * ωs + ϕs)) / √(n_ζ) .* ζs)
#     # Multiply vector of ζ's by the polarization matrix ε
#     res = f * ε
#     return res
# end
# a
# system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")
# d = 60
# τ = 100                             # Simulation time
# δ = system.δ                       # Time step
# α = 40                              # Distance between chain atoms
# μ = 1
# σ0 = [Int(10.5 * α)]
# n_pts = τ / δ |> floor |> Int
# nChain = 20
# ρHs = zeros(nChain, n_pts)
# tTraj = ThermalTrajectory(system.ωmax, δ, ρHs, nothing)
# mem = Inf

# bias = 0.0

# Φ0 = 2.5
# λ = 1
# σdot0 = [25]

# res =
#     motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ; box = (4.5 * α, 15.5 * α))
# lines(res.τs, res.σs[1, :])
# res.σs
# z
# function full_traj(param)
#     println(param)
#     Φ0 = param[1]
#     λ = param[2]
#     σdot0 = param[3]
#     if (
#         !isfile(
#             "data/Non_Thermal/Single_sigma0$(σ0)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(nothing)_tau$(τ).jld2",
#         )
#     )


# system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")
# system.τs

# nMasses = 10000
# ωT = 1
# ωmax = 10
# θs = pi .* (1:nMasses) / nMasses
# ωs = ω.(ωmax, θs)
# ϕs = 2 * π * rand(nMasses)
# ζs = ζq.(ωs, ωT)

# @time ρH(1, ζs, ϕs, ωs, 0:300) |> real


# 1.25 * 40 / 10 * 600

# rr = load_object("data/NumericalLoss_Φ05_λ1_ωmax10_v80_μ1_ωT1.jld2")
# mean(rr[2])
# Δ_thermal_analytic(30, 5, 1, 10, 1)
# hist(rr[2], bins = 50)
# (1 / 2 / (2 * π)^2 * (15^2))
# (1 / 2 / (2 * π)^2 * (10^2)) * 0.5
# (1 / 2 / (2 * π)^2 * (6^2))
# Δ_analytic(80, 5, 1, 10)

# zz = load_object("data/NumericalLossFast_Φ01_λ1_ωmax10_μ1_ωT5.0.jld2")
# hist(zz[2][end, :], bins = 40)
# mean(zz[2][4, :])
# Δ_thermal_analytic(24, 1, 1, 10, 1)


# function Δ_thermal_analytic(v, Φ, λ, Ω, ωT)
#     factor = π^2 * λ^2 * Φ^2 * √(π) / v / 2

#     C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)
#     # C_NN_func(τ) = 0

#     int_func(θ, τ) =
#         exp(2im * π * τ * ω(Ω, θ)) *
#         exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) *
#         (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2 * τ^2) /
#         (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)

#     res = hcubature(x -> int_func(x[1], x[2]), [0, -1], [π / 2, 1], atol = 1e-4)

#     return factor * res[1] * 2 / π |> real
# end

# rr = load_object("data/NumericalLoss_Φ05_λ1_ωmax10_v80_μ1_ωT1.jld2")
# mean(rr[2])
# Δ_thermal_analytic(80, 5, 1, 10, 10)
# Δ_analytic(80, 5, 1, 10)
# (mean(rr[2]) - Δ_thermal_analytic(80, 5, 1, 10, 1)) /
# sqrt(Δ2_thermal_analytic(80, 5, 1, 10, 1) - Δ_thermal_analytic(80, 5, 1, 10, 1)^2) *
# sqrt(length(rr[2]))


# Δ_thermal_analytic(10, 5, 1, 10, 5000)



# ts = range(0, 10, length = 1000)
# temp = 50
# r = -C_corr.(ts, 0, 10, temp) .+ C_corr(0, 0, 10, temp)
# lines(ts, r)
# C_corr(20, 0, 10, 10000)




# 2 * pi * 10
# @time Δ2_thermal_analytic(2, 5, 1, 10, 100)

# vs = range(1000, 1500, length = 20)
# vs = range(4, 5, length = 10)


# r = Δ2_thermal_analytic.(vs, 1, 1, 10, 150)

# lines(log10.(vs), log10.(real(r)))

# (log10.(real(r)[end]) - log10.(real(r)[1])) / (log10.(vs[end]) - log10.(vs[1]))



# C_NN_func(τ) = C_corr(τ, 0, 10, 10)
# v = 0.3;
# int_func(τ) =
#     v *
#     exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) *
#     (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2 * τ^2) /
#     (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)
# lines(τs, int_func.(τs) |> real)
# int_func.(τs)
