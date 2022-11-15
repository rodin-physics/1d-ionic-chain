include("../../src/main.jl")

## Parameters 
α = 40
μ = 1
ωmax = 10
nChain = 100
Φ0 = 2.0
λ = 4.0 

n_sys = 500
n_runs = 100

speed = 40
ωT = 5.0

system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
# system = load_object("precomputed/systems/System_ωmax10_d60_l10000_τmax5.jld2")

# Numerically calculate energy loss after single pass
function Δ_numeric(σ_dot, σ0, Φ0, λ, system, tTraj)
    δ = system.δ
    τ = 1.5 * (α / σ_dot)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ)

    σs = res.σs |> vec
    ρs = res.ρs

    # Find index closest to next midpoint; take the thermal motion into account
    chain_idx = searchsortedlast(ρs[:,1], σs[1])
    mod_val = mod(σs[1], res.α)
    mob_final = findmin(abs.(σs .- (ρs[chain_idx+1, :] .+ mod_val)))[2]

    # Calculate final kinetic energy
    v_final = (σs[mob_final] - σs[mob_final-1]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

# Get distribution of energy loss values
function Δ_thermal(σ_dot, σ0, Φ0, λ, system, tTraj)

    function get_tTraj(ind, tTraj)
        return ThermalTrajectory(tTraj.ωmax, tTraj.δ, tTraj.ρHs[ind:ind+nChain, :], tTraj.ωT)
    end
    nPts = length(1:(size(tTraj.ρHs)[1] - nChain))
    res = zeros(nPts)
    p = Progress(nPts)
    Threads.@threads for ii in 1:(size(tTraj.ρHs)[1] - nChain)
        res[ii] = Δ_numeric(σ_dot, σ0, Φ0, λ, system, get_tTraj(ii, tTraj))
        next!(p)
    end 

   return res
end


## Computation loop 
# Δs = zeros(n_runs * n_sys)
# p = Progress(n_runs * n_sys)
# Threads.@threads for ii in 1:n_sys
#     # Generate thermal trajectory for n_runs + nChain masses 
#     τmax = 1.5 * (α / speed)                            # Simulation time
#     δ = system.δ                        # Time step
#     n_pts = floor(τmax / δ) |> Int      # Number of time steps given t_max and δ
#     n_modes = 10000                     # Number of chain masses for simulating ρ0

#     qa_s = 2 * pi .* (1:n_modes) / n_modes
#     ωs = ω.(ωmax, qa_s ./ 2)
#     ϕs = 2 * π * rand(n_modes)
#     ζs = ζq.(ωs, ωT)

#     gs = collect(1:(n_runs + 50))
#     full_res = [real(ρH(n, δ, ζs, ϕs, ωs, gs)) for n in 1:n_pts]
#     full_res = reduce(hcat, full_res)
#     full_tTraj = ThermalTrajectory(ωmax, δ, full_res, ωT)

#     # Get single pass losses for smaller chunks of chain mass systems 
#     for jj in 1:n_runs
#         # Get small chunk of tTraj
#         # function get_tTraj(ind, tTraj)
#         #     return ThermalTrajectory(tTraj.ωmax, tTraj.δ, tTraj.ρHs[ind:ind+nChain, :], tTraj.ωT)
#         # end
#         init_pos = (4.5 + jj) * α
#         Δs[(ii-1)*n_runs + jj] = Δ_numeric(speed, init_pos, Φ0, λ, system, full_tTraj)
#         next!(p)
#     end
# end

# writedlm("Delta_speed$(speed)_temp$(ωT)2.dat", Δs)
 
# tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ5_nmodes100000.jld2")
# Δs = Δ_thermal(speed, 4.5 * α, Φ0, λ, system, tTraj)
# writedlm("Delta_speed$(speed)_temp$(ωT)3.dat", Δs)




Δs = readdlm("Delta_speed$(speed)_temp$(ωT)2.dat")

fig = Figure(resolution=(1600, 1200), font="CMU Serif", fontsize=40, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel= "Number of Values", ylabel=L"\Delta", title = L"\omega_T = %$(ωT), \, \Phi_0 = %$(Φ0), \,\lambda = %$(λ), \, \dot{\sigma}_0 = %$(speed)")

test = map(nPts -> mean(Δs[1:nPts]), 1:100:length(Δs))
append!(test, mean(Δs))
scatter!(ax1, 1:length(test), test, markersize = 14)
hlines!(ax1, [Δ_analytic(speed, Φ0, λ, ωmax)], linewidth = 3, color = my_black)
hlines!(ax1, [Δ_thermal_analytic(speed, Φ0, λ, ωmax, ωT)], linewidth = 3, color = my_vermillion)

ylims!(ax1, 0.0, 0.03)

fig
# save("convergence_fig2.pdf", fig)