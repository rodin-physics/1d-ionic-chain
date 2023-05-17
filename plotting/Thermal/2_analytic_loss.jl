include("../../src/main.jl")

## Parameters 
colors = [my_blue, my_sky, my_green, my_yellow, my_orange, my_vermillion, my_red, my_black]
ωTs = [0, 1, 5, 25, 100, 250]
ωmax = 10
data = ["data/Thermal/AnalyticLoss_Φ01_λ1_ωmax10_ωT$(ωT).jld2" for ωT in ωTs]


# Plotting 
fig = Figure(resolution = (1200, 800), fontsize = 36, figure_padding = 30)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"Particle velocity $\dot{\sigma}$",
    ylabel = L"Normalized single-pass loss $\dot{\sigma}\langle\Delta\rangle / \Phi_0^2$",
    xscale = log10,
    yscale = log10,
    xgridvisible = false,
    ygridvisible = false
)

## Analytic thermal data
nPts = 10000
vMin = 0.1
vMax = 1000
vs = range(vMin, vMax, length = nPts) |> collect


for ii in eachindex(data)
    d = load_object(data[ii])
    lines!(
        ax1,
        d[1],
        d[1] .* d[2],
        linewidth = 4,
        color = colors[ii],
        label = "$(ωTs[ii])",
    )
end


## Analytic non-thermal data
α = 40
μ = 1
Φ0 = 1
λ = 1
ωmax = 10

nPts = 3000
vMin = 2
vMax = 1000
vs = range(vMin, vMax, length = nPts) |> collect

noFluctuation = Δ_analytic.(vs, Φ0, λ, ωmax)
lines!(ax1, (vs), (vs) .* (noFluctuation), linewidth = 4, color = my_black, label = "Nonthermal")


## High speed asymptote

vHigh = range(50, 1000, length = 100)
lines!(
    ax1,
    vHigh,
    vHigh .* 2 * π^3 * Φ0^2 * (2 * π * λ)^2 .* (ωmax^2 + 1) ./ vHigh .^ 4,
    linewidth = 4,
    color = my_red,
    linestyle = :dash,
    label = "Fast limit"
)


temp_color = [LineElement(color = color, linewidth = 4) for color in colors[1:end-2]]
res_type = [LineElement(color = my_black, linewidth = 4), LineElement(color = my_red, linewidth = 4)]


Legend(fig[1,1],
    [temp_color, res_type],
    [string.(ωTs), ["Nonthermal", "Fast limit"]],
    [L"\omega_T", nothing],
    tellheight = false,
    tellwidth = false,
    halign = 0.6,
    valign = 0.2,
    margin = (10, 10, 10, 10)
    ) 


xlims!(ax1, 0.1, 1e3)
ylims!(ax1, (1e-3, 1))
fig
