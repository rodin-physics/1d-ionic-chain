include("../../src/main.jl")
## Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 2.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
bias = 0.01

## Computation
xs = range(1,60, length = 2000)
xs2 = range(sqrt(8*π^2*Φ0/μ)+0.1,60, length = 2000)
# Δ_bias_rep = @showprogress map(x -> Δ_transport_broadened(x, Φ0, λ, ωmax, α, μ), xs2)
# Δ_bias_att = @showprogress map(x -> Δ_transport_broadened(x, -Φ0, λ, ωmax, α, μ), xs)

# writedlm("rep_Deltabar.dat", Δ_bias_rep)
# writedlm("att_Deltabar.dat", Δ_bias_att)

## Plotting Delta
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=28, figure_padding = 10)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\bar{\Delta}")

# Predicted drift velocities
vlines!(ax1, [α/n for n in 1:10], linewidth = 1.5, linestyle = :dash, color = my_black)

# Capture speed
vlines!(ax1, [sqrt(8*π^2*Φ0/μ)], linestyle = :dash, linewidth = 2, color = my_red, label = "Min. Speed")

# Bias line
hlines!(ax1, [bias], linewidth = 2, color = my_black, label = "Bias")

# Plot broadened DeltaBar
Δ_bias_rep = readdlm("rep_Deltabar.dat") |> vec
Δ_bias_att = readdlm("att_Deltabar.dat") |> vec
scatterlines!(ax1, xs, Δ_bias_att, color = my_blue, markersize = 7, label = "Attractive")
scatterlines!(ax1, xs2, Δ_bias_rep, color = my_vermillion, markersize = 7, label = "Repulsive")


xlims!(ax1, 0, 60)
ylims!(ax1, 0, nothing)
axislegend(ax1, labelsize = 22)
fig
