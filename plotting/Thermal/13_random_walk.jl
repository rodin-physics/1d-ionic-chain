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

mean_data = readdlm("data/Thermal/Full_Trajectory/MeanΔ_ωT10.0_Φ20.0_λ4.0_α40.dat")
var_data = readdlm("data/Thermal/Full_Trajectory/VarΔ_ωT10.0_Φ20.0_λ4.0_α40.dat")
speeds = range(√(8*pi^2 * Φ0), √(8*pi^2 * Φ0) + 20, length = 100)

## Truncated analytics
speeds = mean_data[1, 25:139]
cutoffs = ((μ/8/π^2) .* speeds.^2) .- Φ0
expectations = zeros(length(speeds))
variances = zeros(length(speeds))

pdf_func(x) = 1/√(2*π) * exp(-x^2 / 2)
cdf_func(x) = (1/2) * (1 + erf(x/√2))

# Calculate for each speed 
for ii in eachindex(speeds)
    mean_val = mean_data[2, (25 + ii)]
    var_val = var_data[2, (25 + ii)]
    cutoff = ((μ/8/π^2) * speeds[ii]^2) .- Φ0
    β = (cutoff - mean_val) / sqrt(var_val)

    expectations[ii] = mean_val - sqrt(var_val) * pdf_func(β) / cdf_func(β)
    variances[ii] = var_val * (1 - (β * pdf_func(β) / cdf_func(β)) - (pdf_func(β) / cdf_func(β))^2)
end


## Plotting
fig = Figure(resolution = (1200, 1400), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xgridvisible = false, ygridvisible = false, xticklabelsvisible = false, xticksize = 10, yticklabelsvisible = false)
ax2 = Axis(fig[2, 1], xlabel = L"Velocity midway between lattice sites $\dot{\sigma}$", ylabel = L"Single Pass Loss $\Delta$", xgridvisible = false, ygridvisible = false, xticksize = 10)

# Analytical prediction
for ax in [ax1]
    band!(ax, mean_data[1,:], mean_data[2,:] .- sqrt.(var_data[2,:]), mean_data[2,:] .+ sqrt.(var_data[2,:]), color = (my_blue, 0.3), label = "Analytics")
    lines!(ax, mean_data[1,:], mean_data[2,:], linewidth = 3, color = my_blue, label = "Analytics")
end


# Truncated prediction
lines!(ax2, speeds, expectations, color = my_green, linewidth = 4, label = "Truncated analytics")
band!(ax2, speeds, expectations .- sqrt.(variances), expectations .+ sqrt.(variances), label = "Truncated analytics", linestyle = :dash, linewidth = 2, color = (my_green, 0.3))


## Numerical data 
# without exclusion zone
data = readdlm("data/Thermal/random_walk.dat")
edges, centers, mean_res = binnedStatistic(data[1,:], data[2,:], statistic = :mean, nbins = 100)
lines!(ax1, centers, mean_res, linewidth = 4, color = my_red, label = "Numerics")
edges, centers, std_res = binnedStatistic(data[1,:], data[2,:], statistic = :std, nbins = 100)
band!(ax1, centers, mean_res .- std_res, mean_res .+ std_res, color = (my_red, 0.4), label = "Numerics")

# with exclusion zone 
data = readdlm("data/Thermal/random_walk_exclusion.dat")
edges, centers, mean_res = binnedStatistic(data[1,:], data[2,:], statistic = :mean, nbins = 100)
lines!(ax2, centers, mean_res, linewidth = 4, color = my_red, label = "Truncated numerics")
edges, centers, std_res = binnedStatistic(data[1,:], data[2,:], statistic = :std, nbins = 100)
band!(ax2, centers, mean_res .- std_res, mean_res .+ std_res, color = (my_red, 0.4), label = "Truncated numerics")


# Exclusion zone and capture speed
lines!(ax1, speeds, ((μ/8/π^2) .* speeds.^2) .- Φ0, linestyle = :dash, color = my_black, linewidth = 3)
vlines!(ax1, [√(8*pi^2 * Φ0)], linestyle = :dash, linewidth = 4, color = my_green, label = "Capture Speed")

lines!(ax2, speeds, ((μ/8/π^2) .* speeds.^2) .- Φ0, linestyle = :dash, color = my_black, linewidth = 3, label = "Exclusion Zone")
vlines!(ax2, [√(8*pi^2 * Φ0)], linestyle = :dash, linewidth = 4, color = my_green)


xlims!(ax1, 0, 140)
ylims!(ax1, -10, 10)

xlims!(ax2, 0, 140)
ylims!(ax2, -10, 10)

axislegend(ax1, merge = true, patchsize = (40, 40), position = :rb)
axislegend(ax2, merge = true, patchsize = (40, 40), position = :rb)
fig

