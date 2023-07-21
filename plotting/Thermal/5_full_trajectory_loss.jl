include("../../src/main.jl")
using Random
using BinnedStatistics

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 4.0
ωT = 10.0

## Plotting 
fig = Figure(resolution = (2400, 800), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"Velocity midway between lattice sites $\dot{\sigma}$", ylabel = L"Loss $\Delta$", title = L"Memory $\tau_0$ = 1", yticks = -50:5:50, xgridvisible = false, ygridvisible = false)
ax2 = Axis(fig[1, 2], xticklabelsvisible = false, yticklabelsvisible = false, title = L"$\tau_0$ = 10", yticks = -50:5:50, xgridvisible = false, ygridvisible = false)
ax3 = Axis(fig[1, 3], xticklabelsvisible = false, yticklabelsvisible = false, title = L"$\tau_0$ = 100", yticks = -50:5:50, xgridvisible = false, ygridvisible = false)

## Analytics 
speeds = range(√(8*pi^2 * Φ0), √(8*pi^2 * Φ0) + 20, length = 100)
mean_data = readdlm("data/Thermal/Full_Trajectory/MeanΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")
var_data = readdlm("data/Thermal/Full_Trajectory/VarΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")
new_var = var_data[2,:]
inds = findall(x -> x>= 0, new_var)

for ax in [ax1, ax2, ax3]
    vlines!(ax,√((8*π^2*Φ0)/μ), linewidth = 4, linestyle = :dash, color = my_green, label = "Capture Speed")
    lines!(ax, speeds, ((μ/8/π^2) .* speeds.^2) .- Φ0, linestyle = :dash, color = my_black, linewidth = 3, label = "Exclusion Zone")

    band!(ax, mean_data[1,:][inds], mean_data[2,:][inds] .- sqrt.(new_var[inds]), mean_data[2,:][inds] .+ sqrt.(new_var[inds]), color = (my_blue, 0.3))
    lines!(ax, mean_data[1,:][inds], mean_data[2,:][inds], linewidth = 4, color = my_blue)

end

## Numerics
mem_ax_pairs = [(1, ax1), (10, ax2), (100, ax3)]

for pair in mem_ax_pairs
    data = readdlm("data/Thermal/Full_Trajectory/deltas_mem$(pair[1])_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")

    edges, centers, mean_res = binnedStatistic(data[1,:], data[2,:], statistic = :mean, nbins = 150)
    lines!(pair[2], centers, mean_res, linewidth = 4, color = my_red)

    edges, centers, var_res = binnedStatistic(data[1,:], data[2,:], statistic = :var, nbins = 150)

    band!(pair[2], centers, mean_res .- sqrt.(var_res), mean_res .+ sqrt.(var_res), color = (my_red, 0.4))

    xlims!(pair[2], 0, 140)
    ylims!(pair[2], -10, 10)
end

# Legend entries
axislegend(ax2, position = :rb, patchsize = (50, 50), rowgap = -10)

mean_entries = [LineElement(color = (my_blue, 0.5), linewidth = 5), LineElement(color = (my_red, 0.5), linewidth = 5)]
res_type = [PolyElement(color = (my_blue, 0.5), strokecolor = :transparent), PolyElement(color = (my_red, 0.5), strokecolor = :transparent)]

Legend(fig[1,3],
    [mean_entries, res_type],
    [["Analytics", "Numerics"], [" ", " "]],
    ["Mean", "Fluctuation"],
    nbanks = 2,
    orientation = :horizontal,
    height = 180,
    width = 440,
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :bottom,
    margin = (10, 10, 10, 10),
    patchsize = (30, 20),
    ) 

fig
