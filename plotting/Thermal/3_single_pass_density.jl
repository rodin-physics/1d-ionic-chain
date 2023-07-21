include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
nChain = 3
Φ0 = 2.0
λ = 4.0

σdot0 = 50.0
colors = [my_vermillion, my_green, my_blue]
ωTs = [25.0, 5.0, 0.0]


## Plotting density compared to predicted distribution
fig = Figure(resolution = (1800, 1000), fontsize = 34, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel = L"Single-pass loss $\Delta$", ylabel = "Probability Density", xgridvisible = false, ygridvisible = false, yscale = log10)

for ωT in ωTs 
    # Read in data
    Δs = readdlm("data/Thermal/ωT$(ωT)/delta_ωmax10_speed$(σdot0)_ωT$(ωT).dat") |> vec

    # Fit to histogram and plot
    hist_fit = fit(Histogram, Δs, nbins = 100)
    hist_fit = normalize(hist_fit, mode = :pdf)
    pos_idx = findall(x -> x > 1e-12, hist_fit.weights)
    barplot!(ax1, (((hist_fit.edges[1])[1:end-1] + (hist_fit.edges[1])[2:end]) ./ 2)[pos_idx], hist_fit.weights[pos_idx], fillto = 1e-10, color = (colors[findfirst(x -> x == ωT, ωTs)],  0.6), gap = 0)

    # Analytic mean and variance
    pred_mean = Δ_thermal_analytic(σdot0, Φ0, λ, ωmax, ωT)
    pred_secmom = Δ_thermal_variance(σdot0, Φ0, λ, ωmax, ωT)
    pred_std = √(pred_secmom)

    # Plot normal distribution based on analytics 
    xs = range(-3.0, 3.0, length = 1000)
    ys = 1/(pred_std * sqrt(2π)) .* exp.(-((xs .- pred_mean)).^2 / (2 .* pred_std.^2))
    pos_idx = findall(x -> x > 1e-20, ys)
    lines!(ax1, xs[pos_idx], ys[pos_idx], linewidth = 4, color = colors[findfirst(x -> x == ωT, ωTs)])
end

# Legend entries
temp_color = [PolyElement(color = color, strokecolor = :transparent) for color in reverse(colors)]
res_type = [PolyElement(color = my_black, strokecolor = :transparent), LineElement(color = my_black, linewidth = 4)]
ωTs = [0, 5, 25]

Legend(fig[1,1],
    [temp_color, res_type],
    [string.(ωTs), ["Numerics", "Analytics"]],
    [L"\omega_T", nothing],
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :top,
    margin = (10, 10, 10, 10)) 


xlims!(ax1, -3.0, 3.0)
ylims!(ax1, 1e-8, 20.0)
fig