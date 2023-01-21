include("../src/main.jl")

# Potential broadening
ωmax = 10
nτs = 1000
τmin = 0
τmax = 10
τs = range(τmin, τmax, length = nτs)
ωTs = [0, 1, 5, 10]

colors = [my_vermillion, my_orange, my_green, my_sky, my_blue]

fig = Figure(
    resolution = (1200, 800),
    fonts = (; math = "CMU Serif"),
    fontsize = 40,
    figure_padding = 30,
)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\tau",
    ylabel = L"[C_0(0) - C_0(\tau)]/C_0(0)",
    xticklabelfont = :math,
    yticklabelfont = :math,
    xgridvisible = false,
    ygridvisible = false,
    xlabelpadding = -15,
    limits = (0, 10, 0, 1.5),
)

ax2 = Axis(
    fig,
    bbox = BBox(750, 750 + 360, 200, 200 + 240),
    xlabel = L"\omega_T",
    ylabel = L"C_0(0)",
    xticklabelfont = :math,
    yticklabelfont = :math,
    xgridvisible = false,
    ygridvisible = false,
    xlabelpadding = -15,
    xlabelsize = 40,
    xticklabelsize = 40,
    ylabelsize = 40,
    yticklabelsize = 40,
    limits = (-1, 11, 0, 1.1),
)

for ii in eachindex(ωTs)
    ωT = ωTs[ii]
    res = [
        (C_corr(0, 0, ωmax, ωT) - C_corr(τ, 0, ωmax, ωT)) / C_corr(0, 0, ωmax, ωT) for
        τ in τs
    ]
    lines!(
        ax1,
        τs,
        res,
        color = colors[ii],
        linewidth = 4,
        label = L"\omega_T = %$(ωTs[ii])",
    )
end
scatter!(ax2, ωTs, C_corr.(0, 0, ωmax, ωTs), markersize = 15, color = my_blue)
lines!(
    ax2,
    range(0, maximum(ωTs), length = 100),
    C_corr.(0, 0, ωmax, range(0, maximum(ωTs), length = 100)),
    linewidth = 4,
    color = my_blue,
)

axislegend(ax1, position = :rt, orientation = :horizontal)
save("Broadening_parameter.pdf", fig)
# fig

#  Analytic Delta
α = 40
μ = 1
Φ0 = 1
λ = 1
ωmax = 10
# Analytic dissipation for a non-thermal Gaussian potential
function Δ_analytic(v, Φ, λ, Ω)
    z = (2 * π * λ / v)^2
    return (
        4 * π^3 * Φ^2 / v^2 *
        z *
        exp(-z * (Ω^2 + 1) / 2) *
        (
            besseli(0, z * (Ω^2 - 1) / 2) +
            (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
        )
    )
end

ωTs = [0, 1, 5, 25, 100, 250]

colors = [my_red, my_vermillion, my_orange, my_green, my_sky, my_blue]

data = [
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT1.0e-5.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT1.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT5.0.jld2",
    # "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT10.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT25.0.jld2",
    # "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT50.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT100.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT250.0.jld2",
]

nPts = 6000
vMin = 0.1
vMax = 300
vs = range(vMin, vMax, length = nPts) |> collect

fig = Figure(
    resolution = (1200, 800),
    fonts = (; math = "CMU Serif"),
    fontsize = 40,
    figure_padding = 30,
)

ax1 = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}",
    ylabel = L"\dot{\sigma}\langle\Delta\rangle / \Phi^2_0",
    xscale = log10,
    yscale = log10,
    xticklabelfont = :math,
    yticklabelfont = :math,
    xgridvisible = false,
    ygridvisible = false,
    xlabelpadding = -15,
    limits = (1e-1, 300, 1e-2, 1),
)

for ii in eachindex(data)
    d = load_object(data[ii])
    lines!(
        ax1,
        (vs),
        (vs) .* (d),
        linewidth = 4,
        color = colors[ii],
        label = L"\omega_T = %$(ωTs[ii])",
    )
end

nPts = 1000
vMin = 2
vMax = 300
vs = range(vMin, vMax, length = nPts) |> collect

noFluctuation = Δ_analytic.(vs, Φ0, λ, ωmax)
lines!(ax1, (vs), (vs) .* (noFluctuation), linewidth = 4, color = my_black)

axislegend(position = :lb, nbanks = 2)

vHigh = range(50, 300, length = 10)
lines!(
    ax1,
    vHigh,
    vHigh .* 2 * π^3 * Φ0^2 * (2 * π * λ)^2 .* (ωmax^2 + 1) ./ vHigh .^ 4,
    linewidth = 4,
    color = my_black,
    linestyle = :dash,
)

save("Mean_loss.pdf", fig)

# Broadened potential
nτs = 1000
τmax = 10
τs = range(0, τmax, length = nτs)

fig = Figure(
    resolution = (1200, 800),
    fonts = (; math = "CMU Serif"),
    fontsize = 40,
    figure_padding = 30,
)

ax1 = Axis(
    fig[1, 1],
    xlabel = L"\tau",
    ylabel = L"\int dx \,[\tilde{\Phi}(x,\tau) / \Phi_0]^2",
    xticklabelfont = :math,
    yticklabelfont = :math,
    xgridvisible = false,
    ygridvisible = false,
    xlabelpadding = -15,
    limits = (0, τmax, 0, 1.05),
)

for ii in eachindex(ωTs)
    ωT = ωTs[ii]
    lines!(
        ax1,
        τs,
        sqrt(π) * λ^2 / 2 ./
        sqrt.(λ^2 .+ C_corr(0, 0, ωmax, ωT) .- C_corr.(τs, 0, ωmax, ωT)),
        color = colors[ii],
        linewidth = 4,
        label = L"\omega_T = %$(ωTs[ii])",
    )
end
axislegend(position = :rt, orientation = :vertical, nbanks = 3)

save("Phi_prime_integral.pdf", fig)



# function Φ_broad(x, λ, τ, ωT)
#     C = C_corr(0, 0, ωmax, ωT) - C_corr(τ, 0, ωmax, ωT)
#     res = exp(-x^2 / 2 / (C + λ^2))
#     return (λ * res / √(C + λ^2))
# end

# τs = [0, 0.5, 1, 1.5]
# ωTs = [0, 1, 5, 25, 100, 250]

# nxs = 500
# xmax = 10
# xs = range(-xmax, xmax, length = nxs)

# fig = Figure(
#     resolution = (1200, 800),
#     fonts = (; math = "CMU Serif"),
#     fontsize = 40,
#     figure_padding = 30,
# )

# ax1 = Axis(
#     fig[1, 1],
#     xlabel = L"x",
#     ylabel = L"\tilde{\Phi}(x) / \Phi_0",
#     xticklabelfont = :math,
#     yticklabelfont = :math,
#     xgridvisible = false,
#     ygridvisible = false,
#     xlabelpadding = -15,
#     limits = (-xmax, xmax, 0, 1.05),
# )

# # for ii in eachindex(ωTs)
# #     ωT = ωTs[ii]
# #     lines!(
# #         ax1,
# #         xs,
# #         Φ_broad.(xs, λ, 1.5, ωT),
# #         color = colors[ii],
# #         linewidth = 4,
# #         label = L"\omega_T = %$(ωTs[ii])",
# #     )
# # end


# for ii in eachindex(ωTs)
#     ωT = ωTs[ii]
#     lines!(
#         ax1,
#         xs,
#         Φ_broad.(xs, λ, 1.5, ωT),
#         color = colors[ii],
#         linewidth = 4,
#         label = L"\omega_T = %$(ωTs[ii])",
#     )
# end

# axislegend(position = :rt)

# fig


# @time quadgk(x -> Φ_broad(x, λ, 0, 1), -Inf, Inf)
