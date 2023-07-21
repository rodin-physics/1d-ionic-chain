include("../../src/main.jl")

step_size = 20

## Plotting
fig = Figure(resolution = (2400, 1400), fontsize = 38)

ax1 = Axis(fig[1, 1], title = L"Thermal frequency $\omega_T = 0.0$", xticklabelsvisible = false, yticklabelsvisible = false, titlesize = 45, xgridvisible = false, ygridvisible = false, xticksize = 8, yticksize = 8)
ax2 = Axis(fig[2, 1], xlabel = L"Evolution time $\tau$", ylabel = L"Mobile particle position $\sigma$", xgridvisible = false, ygridvisible = false, xticksize = 8, yticksize = 8)

ax3 = Axis(fig[1, 2], title = L"\omega_T = 5.0", xticklabelsvisible = false, yticklabelsvisible = false, titlesize = 45, xgridvisible = false, ygridvisible = false, xticksize = 8, yticksize = 8)
ax4 = Axis(fig[2, 2], xticklabelsvisible = false, yticklabelsvisible = false, xgridvisible = false, ygridvisible = false, xticksize = 8, yticksize = 8)

ax5 = Axis(fig[1, 3], title = L"\omega_T = 25.0", xticklabelsvisible = false, yticklabelsvisible = false, titlesize = 45, xgridvisible = false, ygridvisible = false, xticksize = 8, yticksize = 8)
ax6 = Axis(fig[2, 3], xticklabelsvisible = false, yticklabelsvisible = false, xgridvisible = false, ygridvisible = false, xticksize = 8, yticksize = 8)

ωTs = [0.0, 5.0, 25.0]
ax_pairs = [(ax1, ax2), (ax3, ax4), (ax5, ax6)]


for ax_pair in ax_pairs
    ωT = ωTs[findfirst(x -> x == ax_pair, ax_pairs)]

    # REPULSIVE
    data = load_object("data/Thermal/Single_Particle/Single_sigma0[180]_sigmadot0[120]_MemInf_lambda4_Phi20_mu1_d60_bias0.0_omegaT$(ωT)_tau200.jld2")

    δ = data.τs[2] - data.τs[1]
    rr = reduce(hcat, [data.ρs[ii, :] .- ii * data.α for ii = 1:size(data.ρs)[1]])
    # mx = maximum(abs.(rr))
    mx = 1.5
    hm = heatmap!(
        ax_pair[1],
        data.τs[1:step_size:end],
        collect(1:size(data.ρs)[1]) .* data.α,
        rr[1:step_size:end, :],
        colormap = :RdBu,
        colorrange = (-mx, mx),
    )
    lines!(ax_pair[1], data.τs[1:3:end], vec([x[1] for x in data.σs])[1:3:end], color = my_black, linewidth = 5)
    xlims!(ax_pair[1], (0, 200))
    ylims!(ax_pair[1], (0, 10000))

    # ATTRACTIVE
    data = load_object("data/Thermal/Single_Particle/Single_sigma0[180]_sigmadot0[120]_MemInf_lambda4_Phi-20_mu1_d60_bias0.0_omegaT$(ωT)_tau200.jld2")
    δ = data.τs[2] - data.τs[1]
    rr = reduce(hcat, [data.ρs[ii, :] .- ii * data.α for ii = 1:size(data.ρs)[1]])
    # mx = maximum(abs.(rr))
    mx = 1.5
    hm = heatmap!(
        ax_pair[2],
        data.τs[1:step_size:end],
        collect(1:size(data.ρs)[1]) .* data.α,
        rr[1:step_size:end, :],
        colormap = :RdBu,
        colorrange = (-mx, mx))

lines!(ax_pair[2], data.τs[1:3:end], vec([x[1] for x in data.σs])[1:3:end], color = my_black, linewidth = 5)
xlims!(ax_pair[2], (0, 200))
ylims!(ax_pair[2], (0, 10000))
end

## Annotations
text!(ax1, 10, 8600, text = "(a)", fontsize = 45, font = :bold)
text!(ax2, 10, 8600, text = "(b)", fontsize = 45, font = :bold)
text!(ax3, 10, 8600, text = "(c)", fontsize = 45, font = :bold)
text!(ax4, 10, 8600, text = "(d)", fontsize = 45, font = :bold)
text!(ax5, 10, 8600, text = "(e)", fontsize = 45, font = :bold, color = :white)
text!(ax6, 10, 8600, text = "(f)", fontsize = 45, font = :bold, color = :white)

text!(ax1, 130, 200, text = "Repulsive", fontsize = 48)
text!(ax2, 130, 200, text = "Attractive", fontsize = 48)

Colorbar(fig[:, 4], limits = (-1.5,1.5), colormap = :RdBu; label = L"Chain mass displacement $\delta\rho$", width = 15, ticksize = 15, tickalign = 1)

rowgap!(fig.layout, 30)
colgap!(fig.layout, 45)
fig