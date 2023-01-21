include("../src/main.jl")

# Potential broadening

nτs = 1000
τs = range(τmin, τmax, length = nτs)
ωTs = [1e-5, 1, 5, 10, 20]
fig =
    Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 40, figure_padding = 30)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}",
    ylabel = L"\dot{\sigma}\langle\Delta\rangle",
    xscale = log10,
    yscale = log10,
)










nωTs = 100
nτs = 1000
ωmax = 10

ωTmin = 1e-5
ωTmax = 50
τmin = 0
τmax = 20

ωTs = range(ωTmin, ωTmax, length = nωTs)
τs = range(τmin, τmax, length = nτs)

res = zeros(nτs, nωTs)

ωT = 20
res2 =
    [(C_corr(0, 0, ωmax, ωT) - C_corr(τ, 0, ωmax, ωT)) / C_corr(0, 0, ωmax, ωT) for τ in τs]
lines(τs, res2)
res = [C_corr(τ, 0, ωmax, ωT) for τ in τs, ωT in ωTs]
res2 = [C_corr(0, 0, ωmax, ωT) - C_corr(τ, 0, ωmax, ωT) for τ in τs, ωT in ωTs]

res = [C_corr(0, 0, ωmax, ωT) for ωT in ωTs]
lines(ωTs, res)
# ωT
# heatmap(τs, ωTs, res2, colormap = :oslo)
# lines(τs, res2[:,end])

# # Analytic dissipation for Gaussian potential
# function Δ_analytic(v, Φ, λ, Ω)
#     z = (2 * π * λ / v)^2
#     return (
#         4 * π^3 * Φ^2 / v^2 *
#         z *
#         exp(-z * (Ω^2 + 1) / 2) *
#         (
#             besseli(0, z * (Ω^2 - 1) / 2) +
#             (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
#         )
#     )
# end

# # Mean Delta with thermal fluctuations 
# function Δ_thermal_analytic(v, Φ, λ, Ω, ωT)
#     factor = π^2 * λ^2 * Φ^2 * √(π) / v / 2

#     C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)
#     τlim = min(30 / v, 3)
#     int_func(θ, τ) =
#         exp(2im * π * τ * ω(Ω, θ)) *
#         exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) *
#         (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2 * τ^2) /
#         (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)

#     res = hcubature(
#         x -> int_func(x[1], x[2]),
#         [0, -τlim],
#         [π / 2, τlim],
#         rtol = 1e-2,
#         initdiv = 1,
#     )

#     return factor * res[1] * 2 / π |> real
# end

# # Mean Delta with thermal fluctuations 
# function Δ2_thermal_analytic(v, Φ, λ, Ω, ωT)
#     factor = π^2 * λ^2 * Φ^2 * √(π) / v / 2

#     C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)
#     int_func(τ) =
#         exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) *
#         (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2 * τ^2) /
#         (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)
#     (sol, err) = quadgk(int_func, -Inf, Inf, atol = 1e-8)

#     return (sol * factor * 2 / 4 / pi^2 * v^2)
# end



# ## SYSTEM PARAMETERS
# α = 40
# μ = 1
# Φ0 = 1
# λ = 1
# ωmax = 10
# ωTs = [1e-5, 1, 5, 10, 25, 50, 100, 250, 1000]
# # , 2500, 10000]

# ## ANALYTIC MEAN
# Φ0 = 1

# nPts = 6000
# vMin = 0.1
# vMax = 300
# # vs = 10 .^ range(log10(vMin), log10(vMax), length = nPts)
# # scatter(log10.(vs), log10.(vs))

# vs = range(vMin, vMax, length = nPts) |> collect
# vs_idx = reduce(
#     vcat,
#     [collect(eachindex(vs))[n:Threads.nthreads():end] for n = 1:Threads.nthreads()],
# )

# for ωT in ωTs
#     println("ωT = $(ωT)")
#     if (!isfile("data/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2"))
#         result = zeros(nPts)
#         pr = Progress(nPts)
#         Threads.@threads for ii in vs_idx
#             result[ii] = Δ_thermal_analytic(vs[ii], Φ0, λ, ωmax, ωT)
#             next!(pr)
#         end
#         save_object("data/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2", result)
#     end
# end


# # nPts = 200
# # vMin = 2
# # vMax = 20
# # vs = range(vMin, vMax, length = nPts) |> collect
# # vs_idx = reduce(
# #     vcat,
# #     [collect(eachindex(vs))[n:Threads.nthreads():end] for n = 1:Threads.nthreads()],
# # )

# # for ωT in ωTs
# #     println("ωT = $(ωT)")
# #     if (!isfile("data/AnalyticLossSlow_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2"))
# #         result = zeros(nPts)
# #         pr = Progress(nPts)
# #         Threads.@threads for ii in vs_idx
# #             result[ii] = Δ_thermal_analytic(vs[ii], Φ0, λ, ωmax, ωT)
# #             next!(pr)
# #         end
# #         save_object("data/AnalyticLossSlow_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2", result)
# #     end
# # end


# # @time Δ_thermal_analytic(300, 1, 1, 10, 1)
# # @time Δ_thermal_analytic(20, 1, 1, 10, 10)
# # Δ_analytic(20, 1, 1, 10)


# function Δ_thermal_analytic_CHECK(v, Φ, λ, Ω, ωT)
#     factor = π^2 * λ^2 * Φ^2 * √(π) / v / 2

#     C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)
#     τlim = min(30 / v, 3)
#     τlim = 1
#     int_func(θ, τ) =
#         exp(2im * π * τ * ω(Ω, θ)) *
#         (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))) /
#         (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)

#     res = hcubature(
#         x -> int_func(x[1], x[2]),
#         [0, -τlim],
#         [π / 2, τlim],
#         rtol = 1e-2,
#         initdiv = 1,
#     )

#     return factor * res[1] * 2 / π |> real
# end

# Δ_thermal_analytic_CHECK(50,1, 1, 10, 2000)
# Δ_thermal_analytic(50,1, 1, 10, 2000)

# Ω = 10
# ωT = 60
# C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)
# λ = 1
# f(τ) =
#     (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))) /
#     (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)
# lines(ts, f.(ts)|>real)
# f.(ts)
