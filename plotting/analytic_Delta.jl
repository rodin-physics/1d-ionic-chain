include("../src/main.jl")

colors = [
    my_vermillion,
    my_orange,
    my_green,
    my_sky,
    my_blue,
    my_vermillion,
    my_orange,
    my_green,
    my_sky,
    my_blue,
    my_red,
]
ωTs = [0, 1, 5, 10, 25, 50, 100, 250, 1000, 2500, 10000]
ωmax = 10
data = [
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT1.0e-5.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT1.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT5.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT10.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT25.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT50.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT100.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT250.0.jld2",
    "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT1000.0.jld2",
    # "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT2500.0.jld2",
    # "data/AnalyticLoss_Φ01_λ1_ωmax10_ωT10000.0.jld2",
]

nPts = 6000
vMin = 0.1
vMax = 300
vs = range(vMin, vMax, length = nPts) |> collect

fig =
    Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 40, figure_padding = 30)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}",
    ylabel = L"\dot{\sigma}\langle\Delta\rangle",
    xscale = log10,
    yscale = log10,
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

α = 40
μ = 1
Φ0 = 1
λ = 1
ωmax = 10

nPts = 1000
vMin = 2
vMax = 300
vs = range(vMin, vMax, length = nPts) |> collect

noFluctuation = Δ_analytic.(vs, Φ0, λ, ωmax)
lines!(ax1, (vs), (vs) .* (noFluctuation), linewidth = 4, color = my_black)

axislegend(position = :rt, nbanks = 2)



vHigh = range(50, 300, length = 10)
lines!(
    ax1,
    vHigh,
    vHigh .* 2 * π^3 * Φ0^2 * (2 * π * λ)^2 .* (ωmax^2 + 1) ./ vHigh .^ 4,
    linewidth = 4,
    color = my_black,
    linestyle = :dash,
)

ylims!(ax1, (0.01, 1))
# z = load_object("data/NumericalLossFast_Φ01_λ1_ωmax10_μ1_ωT1.0e-5.jld2")
# # z = load_object("data/NumericalLossFast_Φ01_λ1_ωmax10_μ1_ωT1.0e-5.jld2")
# mn = [mean(z[2][ii, :]) for ii = 1 : 6]
# scatter!(ax1, log.(z[1]), log.(mn))
fig
# save("MeanLoss.pdf", fig)
# l.attributes


# d = [load_object(r)[1] for r in data]

# scatter(log.(ωTs[2:end]), log.(d[2:end]))
# # # ωTs
# # log.(ωTs)
# C_corr(0,0,10,2*pi * 10)
