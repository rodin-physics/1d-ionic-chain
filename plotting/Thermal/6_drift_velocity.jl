include("../../src/main.jl")

## Stacked Histogram
timings = string.(0:250:1000)
timings = collect(Iterators.flatten(zip(timings, repeat([" "], 5))))
timings = vcat(L"\tau", timings)

fig = Figure(resolution = (1800, 800), fontsize = 28, figure_padding = 40)
ax1 = Axis(fig[1, 1], title = L"Thermal frequency $\omega_T$ = 0", xlabel = L"Particle speed $\dot{\sigma}", ylabel = "Probability density", yticks = 0:0.08:0.4, yticklabelsvisible = false, xgridvisible = false)
ax2 = Axis(fig[1, 2], title = L"\omega_T = 5", yticklabelsvisible = false, ygridvisible = false, yticksvisible = false, xgridvisible = false, xticklabelsvisible = false)
ax3 = Axis(fig[1, 3], title = L"\omega_T = 25", yticklabelsvisible = false, ygridvisible = false, yticksvisible = false, xgridvisible = false, xticklabelsvisible = false)
ax4 = Axis(fig[1, 4], title = L"\omega_T = 50", ygridvisible = false, yaxisposition = :right, yticklabelsvisible = false, yticksvisible = false, xgridvisible = false, xticklabelsvisible = false)
ax5 = Axis(fig[1,5], width = 30, xticklabelsvisible = false, xticksvisible = false, xgridvisible = false, ygridvisible = false, yticks = (0.0:0.04:0.4, reverse(timings)), yticksize = 20, yticksvisible = false, yticklabelpad = 20, yticklabelalign = (:center, :center))

ax_pair = [(0.0, ax1, my_blue), (5.0, ax2, my_sky), (25.0, ax3, my_green), (50.0, ax4, my_vermillion)]

for pair in ax_pair 
    data = readdlm("data/Thermal/Drift_Velocity/Speeds_1000τ_inter250_α40_Φ2.0_λ4.0_bias0.01_T$(pair[1]).dat")
    num_times = size(data, 1)

    hlines!(pair[2], 0.08:0.08:0.36, linewidth = 2, color = (my_black, 0.2))

    for t in 1:num_times
        hist!(pair[2], data[t, :], bins = 44, offset = 0.08 * (num_times - t), color = pair[3], normalization = :pdf, strokewidth = 0.3, strokecolor = pair[3])
    end

    drift = vlines!(pair[2], [40/n for n in 1:3], linestyle = :dash, color = my_orange, linewidth = 2.5, label = "Drift velocity")
    cap_speed = vlines!(pair[2], [√(8*π^2*2)], linewidth = 4, linestyle = :dash, color = my_red, label = "Capture speed")
    xlims!(pair[2], -20, 80)
    ylims!(pair[2], 0, 0.4)
end

arrows!(ax5, [0.0], [0.36], [0.0], [-0.345], arrowsize = 20, linewidth = 2)

text!(ax1, -18, 0.38, text = "(a)", fontsize = 26)
text!(ax2, -18, 0.38, text = "(b)", fontsize = 26)
text!(ax3, -18, 0.38, text = "(c)", fontsize = 26)
text!(ax4, -18, 0.38, text = "(d)", fontsize = 26)


# Legend entries
leg_entries = [LineElement(color = my_red, linewidth = 4, linestyle = :dash), LineElement(color = my_orange, linewidth = 3, linestyle = :dash)]

Legend(fig[2,3:4], 
       leg_entries,
       ["Capture speed", "Drift velocity"],
       tellwidth = false, 
       tellheight = true,
       orientation = :horizontal,
       patchsize = (36, 36),
       height = 44)

colgap!(fig.layout, 4, -25)
rowgap!(fig.layout, -45)
hidespines!(ax5)
ylims!(ax5, 0, 0.4)
fig