include("../../src/main.jl")

## Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
system2 = load_object("precomputed/systems/System_ωmax10_d6000_l20.jld2")

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
fig = Figure(resolution=(1200, 1800), font="CMU Serif", fontsize=32, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta", yscale = log10)
ax2 = Axis(fig[2, 1], xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta", yscale = log10)

colors = [my_vermillion, my_orange, my_green, my_sky]
Φ0 = 0.01
λs = [0.5, 1, 2, 4]
init = [2.0, 3.0, 4.0, 8.0]

# Low Phi
for ii in eachindex(λs)
    # Analytic Δ
    xs = range(init[ii], 80, step = 0.5)
    Δs = Δ_analytic.(xs, Φ0, λs[ii], ωmax)
    lines!(ax1, xs, Δs, color = colors[ii], label = L"\lambda = %$(λs[ii])", linewidth = 4)

    # Numeric Δ
    numeric_rep = map(x -> Δ_numeric(x, 5.5 * α, Φ0, λs[ii], system), xs)
    numeric_att = map(x -> Δ_numeric(x, 5.5 * α, -Φ0, λs[ii], system), xs)
    if any(x -> x<0, numeric_att) || any(x -> x<0, numeric_rep)
        println(λs[ii])
    end
    scatter!(ax1, xs, numeric_rep, color = colors[ii], marker = :cross, markersize = 14)
    scatter!(ax1, xs, numeric_att, color = colors[ii], markersize = 14, marker = :hline)
end

# High Phi
Φ0 = 24
for ii in eachindex(λs)
    # Analytic Δ
    xs = range(100, 1250, step = 10.0)
    Δs = Δ_analytic.(xs, Φ0, λs[ii], ωmax)
    lines!(ax2, xs, Δs, color = colors[ii], label = L"\lambda = %$(λs[ii])", linewidth = 4)

    # Numeric Δ
    numeric_rep = map(x -> Δ_numeric(x, 5.5 * α, Φ0, λs[ii], system2), xs)
    numeric_att = map(x -> Δ_numeric(x, 5.5 * α, -Φ0, λs[ii], system2), xs)
    if any(x -> x<0, numeric_att)
        println(λs[ii])
    end
    scatter!(ax2, xs, numeric_rep, color = colors[ii], marker = :cross, markersize = 14)
    scatter!(ax2, xs, numeric_att, color = colors[ii], markersize = 14, marker = :hline)
end

axislegend(ax1, labelsize = 20, unique = true)
axislegend(ax2, labelsize = 20, unique = true)
save("Single_Pass.pdf",fig)
