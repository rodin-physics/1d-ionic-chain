include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 2.0
λ = 4.0
ωTs = [0.0, 1.0, 5.0, 10.0]
colors = [my_blue, my_sky, my_green, my_yellow, my_orange, my_vermillion, my_red, my_black]


## Plotting 
fig = Figure(resolution = (1200, 1600), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], ylabel = L"[C_0(0) - C_0(\tau)]/C_0(0)", xgridvisible = false, ygridvisible = false, xticklabelsvisible = false)
inset_ax = Axis(fig[1, 1], width=Relative(0.3), height=Relative(0.3), halign=0.9, valign=0.3, xlabel = L"Temperature $\omega_T$", ylabel = L"C_0(0)")
ax2 = Axis(fig[2, 1], xlabel = L"Evolution time $\tau$", ylabel = L"1000 \int dx \, [ \tilde{\Phi}'(x,\tau) / \Phi_0 ]^2", xgridvisible = false, ygridvisible = false)

τs = range(0, 10, length = 1000)
ωT_range = range(0, 10, length = 50)

## First axis
for ωT in ωTs 
    broadening_xs = map(x -> (C_corr(0, 0, ωmax, ωT) - C_corr(x, 0, ωmax, ωT)) / C_corr(0, 0, ωmax, ωT), τs)

    lines!(ax1, τs, broadening_xs, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
end

## Inset axis 
for ωT in ωTs
    scatter!(inset_ax, [ωT], [C_corr(0,0,ωmax, ωT)], color = my_blue)
end
factors = map(x -> C_corr(0,0,ωmax, x), ωT_range)
lines!(inset_ax, ωT_range, factors, color = my_blue, linewidth = 3)

## Second axis 
ωTs = [0.0, 1.0, 5.0, 25.0, 100.0, 250.0]
colors = [my_blue, my_sky, my_green, my_orange, my_vermillion, my_red, my_black]

# Temperature induced broadened potential - normalised by Φ0 and squared
function broadened_potential(x, τ, ωT)
    factor = λ^2 + (C_corr(0, 0, ωmax, ωT) - C_corr(τ, 0, ωmax, ωT))
    return (λ  / factor) * exp(- x^2 / factor / 2)
end


for ωT in ωTs 
    res = (λ^2 * √(π) / 2) .* (λ^2 .+ C_corr(0, 0, ωmax, ωT) .- C_corr.(τs, 0, ωmax, ωT)).^(-5/2)
    lines!(ax2, τs, 1000 * res, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
end

# Legend entries 
ωTs = [0, 1, 5, 10, 25, 100, 250]
colors = [my_blue, my_sky, my_green, my_yellow, my_orange, my_vermillion, my_red]

leg_entries = [LineElement(color = cs, linewidth = 4) for cs in colors]

Legend(fig[1:2,1],
      leg_entries,
      string.(ωTs),
      L"\omega_T",
      tellwidth = false,
      orientation = :horizontal,
      titleposition = :left,
      titlegap = 25)


## Annotation
text!(ax1, 1, 1.4, text = "(a)", fontsize = 36)
text!(ax2, 1, 15, text = "(b)", fontsize = 36)

xlims!(ax1, 0, 10)
ylims!(ax1, 0, nothing)

ylims!(inset_ax, 0, nothing)

xlims!(ax2, 0, 10)
ylims!(ax2, 0, 17)


translate!(inset_ax.scene, 0, 0, 10)
translate!(inset_ax.elements[:background], 0, 0, 9)
fig

