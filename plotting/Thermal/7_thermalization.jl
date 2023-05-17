include("../../src/main.jl")

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 1.0
τ = 500
δ = ((1 / ωmax) / 60)
ωTs = [1.0, 2.0, 5.0, 10.0, 25.0, 100.0, 250.0]

colors = [my_blue, my_sky, my_green, my_yellow, my_orange, my_vermillion, my_red, my_black]

# Box Parameters
box_size = 100
left_boundary = 1.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)
nParticles = 20

## Get particle energy distribution
function get_dist(σs, ωT, ens)
    numP = size(σs)[1]
    σs_final = copy(σs)

    for particle in 1:numP
        curr_σs = σs_final[particle, :]

        # Get timesteps at which particle crosses a midpoint 
        midpoint = floor.((curr_σs .- (α / 2)) ./ α)
        idx = findall(x -> x != 0, midpoint[2:end] - midpoint[1:end-1])

        # Calculate kinetic energy at those timesteps 
        kin_en = ((curr_σs[idx] - curr_σs[idx.-1]) ./ δ) .^ 2 ./ 2 ./ ωT ./ (2 * pi)^2
        ens = vcat(ens, kin_en)
    end

    # Get normalized kinetic energy distribution of all numP particles
    hist_fit = fit(Histogram, ens |> vec, 0.0:0.02:5)
    hist_fit = normalize(hist_fit, mode = :pdf)
    return hist_fit
end

## Plotting
fig = Figure(resolution = (2400, 800), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"Scaled energy $\mathcal{E} / \omega_T$", ylabel = L"Probability distribution $\ln(P(\mathcal{E}))$", title = L"Memory $\tau_0$ = 1", xgridvisible = false, ygridvisible = false, yticks = (-8:2:0), yminorticksvisible = true)
ax2 = Axis(fig[1, 2], title = L"$\tau_0$ = 10", xgridvisible = false, ygridvisible = false, yticklabelsvisible = false, xticklabelsvisible = false)
ax3 = Axis(fig[1, 3], title = L"$\tau_0$ = 100", xgridvisible = false, ygridvisible = false, yticklabelsvisible = false, xticklabelsvisible = false)


mem_ax_pair = [(1, ax1), (10, ax2), (100, ax3)]

for pair in mem_ax_pair
    for ωT in ωTs
        println("OmegaT is $(ωT)")

        energies = Float64[]
        σs_data = readdlm("data/Thermal/Thermalization_Single/SingleParticle_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem$(pair[1]).dat")

        hist_data = get_dist(σs_data, ωT, energies)

        scatter!(
            pair[2],
            (hist_data.edges[1])[1:end-1],
            log.(hist_data.weights),
            # log.(hist_data.weights .* sqrt.((hist_data.edges[1])[1:end-1])),
            color = colors[findfirst(x -> x == ωT, ωTs)],
            label = "$(Int(ωT))",
            markersize = 12,)

    end
end

# Line with slope -1 
for ax in [ax1, ax2, ax3]
    lines!(ax, 0:5, -(0:5), color = my_black, linewidth = 4.5)
    xlims!(ax, -0.1, 5.1)
    ylims!(ax, -8, 1)
end

# Annotations
text!(ax1, 4.6, 0, text = "(a)", fontsize = 36)
text!(ax2, 4.6, 0, text = "(b)", fontsize = 36)
text!(ax3, 4.6, 0, text = "(c)", fontsize = 36)

axislegend(L"\omega_T", position = :lb)


fig
