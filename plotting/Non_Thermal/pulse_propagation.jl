include("../../src/main.jl")

ωmax = 10
l_max = 500
times = [1, 5, 10, 15]

Γs = [Γ(τ, 0:l_max, ωmax) for τ in times]
colors = [my_vermillion, my_green, my_sky, my_blue]
speed = (ωmax - 1) / 2 * 2 * π

## Plotting 
fig = Figure(resolution = (600, 400), font = "CMU Serif", fontsize = 18)
ax1 = Axis(fig[1, 1], ylabel = L"\Gamma_n(\tau)", xlabel = L"n")
for ii = 1:length(times)
    sc = scatter!(
        ax1,
        0:l_max,
        Γs[ii],
        markersize = 5,
        color = colors[ii],
        label = L"\tau = %$(times[ii])",
    )

    ln = vlines!(
        ax1,
        times[ii] * [speed],
        linewidth = 2,
        linestyle = :dash,
        color = [colors[ii]],
    )
end
xlims!(ax1, (-0.05, 500))
ylims!(ax1, (-0.1, 0.1))
axislegend(ax1, position = :rb, labelsize = 18, nbanks = 4)

save("Pulse_Propagation.pdf", fig)
