using SpecialFunctions

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 4.0
ωT = 10.0

# Load analytic data
mean_data = readdlm("data/Thermal/Full_Trajectory/MeanΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")
var_data = readdlm("data/Thermal/Full_Trajectory/VarΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")

# Calculate truncated analytics
speeds = range(20, 120, length = 100)
cutoffs = ((μ/8/π^2) .* speeds.^2) .- Φ0
expectations = zeros(length(speeds))
variances = zeros(length(speeds))

pdf_func(x) = 1/√(2*π) * exp(-x^2 / 2)
cdf_func(x) = (1/2) * (1 + erf(x/√2))

for ii in eachindex(speeds)
    mean_val = mean_data[2, (40 + ii)]
    var_val = var_data[2, (40 + ii)]
    cutoff = ((μ/8/π^2) * speeds[ii]^2) .- Φ0
    β = (cutoff - mean_val) / sqrt(var_val)

    expectations[ii] = mean_val - sqrt(var_val) * pdf_func(β) / cdf_func(β)
    variances[ii] = var_val * (1 - (β * pdf_func(β) / cdf_func(β)) - (pdf_func(β) / cdf_func(β))^2)
end

# Load numeric data
numeric_data = readdlm("data/Thermal/Full_Trajectory/deltas_mem100_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")

## Plotting
fig = Figure(resolution = (1200, 800), fontsize = 32, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"Velocity midway between lattice sites $\dot{\sigma}$", ylabel = L"Single-pass loss $\Delta$", xgridvisible = false, ygridvisible = false)


# Truncated prediction
lines!(ax1, speeds, expectations, color = my_green, linewidth = 4, label = "Truncated analytics")
band!(ax1, speeds, expectations .- sqrt.(variances), expectations .+ sqrt.(variances), label = "Truncated analytics", linestyle = :dash, linewidth = 2, color = (my_green, 0.3))


# Numerical data
edges, centers, mean_res = binnedStatistic(numeric_data[1,:], numeric_data[2,:], statistic = :mean, nbins = 150)
lines!(ax1, centers, mean_res, linewidth = 4, color = my_red, label = "Numerics")

edges, centers, var_res = binnedStatistic(numeric_data[1,:], numeric_data[2,:], statistic = :var, nbins = 150)
band!(ax1, centers, mean_res .- sqrt.(var_res), mean_res .+ sqrt.(var_res), color = (my_red, 0.4), label = "Numerics")

# Dashed lines
vlines!(ax1,√((8*π^2*Φ0)/μ), linewidth = 3, linestyle = :dash, color = my_green, label = "Capture Speed")
lines!(ax1, speeds, ((μ/8/π^2) .* speeds.^2) .- Φ0, linestyle = :dash, color = my_black, linewidth = 2.5, label = "Exclusion Zone")


xlims!(ax1, 0, 120)
ylims!(ax1, -10, 10)

axislegend(ax1, position = :rb, labelsize = 30, patchsize = (30, 30), merge = true)
fig