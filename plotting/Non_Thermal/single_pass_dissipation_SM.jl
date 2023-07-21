include("../../src/main.jl")

## Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 20.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")

# Numerically calculate energy loss after single pass
function Δ_numeric(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 2.2 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    σs = res.σs |> vec
    ρs = res.ρs

    # Find index closest to next midpoint
    chain_idx = searchsortedlast(ρs[:,1], σs[1])
    mod_val = mod(σs[1], res.α)
    mob_final = findmin(abs.(σs .- (ρs[chain_idx+1, :] .+ mod_val)))[2]

    # Calculate final kinetic energy
    v_final = (σs[mob_final] - σs[mob_final-1]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

## Plotting Delta
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=24, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta")

# Capture speed
vlines!(ax1, [sqrt(8*Φ0*π^2)], color = my_red, linewidth = 3, label = "Min. Speed")

# Analytic Δ
xs = range(20.0, 120, step = 1.0)
Δs = Δ_analytic.(xs, Φ0, λ, ωmax)
lines!(ax1, xs, Δs, color = my_green, label = "Analytic", linewidth = 4)

# Numeric losses
xs = range(ceil(√(8*π^2*Φ0)), 120, step = 0.5)
xs2 = range(1, 120, step = 0.5)
numeric_rep = map(x -> Δ_numeric(x, 5.5 * α, Φ0, λ, system), xs)
numeric_att = map(x -> Δ_numeric(x, 5.5 * α, -Φ0, λ, system), xs2)
scatter!(ax1, xs, numeric_rep, color = my_vermillion, markersize = 12, label = "Repulsive")
scatter!(ax1, xs2, numeric_att, color = my_blue, markersize = 12, label = "Attractive")

xlims!(ax1, 0, 120)
ylims!(ax1, 0, nothing)
axislegend(ax1, labelsize = 20, unique = true)
fig
# save("Single_Pass_SM.pdf", fig)
