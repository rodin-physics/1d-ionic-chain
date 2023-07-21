include("../src/main.jl")

step_size = 20

## Plotting
fig = Figure(resolution = (1200, 1600), font = "CMU Serif", fontsize = 36)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma")

## REPULSIVE

# Load data
data = load_object(
    "data/non_thermal/Single_sigma0[220]_sigmadot0[120]_MemInf_lambda4_Phi20_mu1_d60_bias0.0_omegaTnothing_tau125.jld2")

τ_max = 125                 # Maximum time plotted
idx = findall(data.τs .< τ_max)
τs = data.τs[idx]
rr = reduce(hcat, [data.ρs[ii, idx] .- ii * data.α for ii = 1:size(data.ρs)[1]])

mx = maximum(abs.(rr))      # Colorbar limits

# Heatmap
hm = heatmap!(
    ax1,
    τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)

# Particle trajectory
lines!(ax1, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)

# Sound cone
lines!(
    ax1,
    0:0.1:10,
    π * data.α * (data.ωmax - 1) * (0:0.1:10) .+ data.σs[1][1],
    color = my_black,
    linewidth = 4,
    linestyle = :dash,
)

# Colorbar
Colorbar(fig[1, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)
xlims!(ax1, (0, 125))
ylims!(ax1, (0, 10000))


## ATTRACTIVE

# Load data
data = load_object(
    "data/non_thermal/Single_sigma0[220]_sigmadot0[120]_MemInf_lambda4_Phi-20_mu1_d60_bias0.0_omegaTnothing_tau125.jld2",
)

τ_max = 125                 # Maximum time plotted
idx = findall(data.τs .< τ_max)
τs = data.τs[idx]
rr = reduce(hcat, [data.ρs[ii, idx] .- ii * data.α for ii = 1:size(data.ρs)[1]])

mx = maximum(abs.(rr))      # Colorbar limits

# Heatmap
hm = heatmap!(
    ax2,
    τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)

# Particle trajectory
lines!(ax2, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)

# Sound cone
lines!(
    ax2,
    0:0.1:10,
    π * data.α * (data.ωmax - 1) * (0:0.1:10) .+ data.σs[1][1],
    color = my_black,
    linewidth = 4,
    linestyle = :dash,
)

# Colorbar
Colorbar(fig[2, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)
xlims!(ax2, (0, 125))
ylims!(ax2, (0, 10000))

# Save figure
save("General_Example.pdf", fig)
