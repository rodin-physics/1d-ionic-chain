include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
nChain = 3
Φ0 = 2.0
λ = 4.0 
# Load memory kernel and thermal trajectory
system = load_object("precomputed/systems/System_ωmax10_d60_l10000_τmax5.jld2")

# Numerically calculate energy loss after single pass
function Δ_numeric(σ_dot, σ0, Φ0, λ, system, tTraj)
    δ = system.δ
    τ = 1.75 * (α / σ_dot)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ)

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

## Computation 
σdot0 = 20.0
ωT = 25.0
ωTs = [0.0]
init_pos = 2.5 * α

# for ωT in ωTs 
#     # tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ5.jld2")
#     # Δs = Δ_thermal(σdot0, init_pos, Φ0, λ, system, tTraj)
#     Δs = readdlm("data/Thermal/ωT$(ωT)/delta_ωmax$(ωmax)_speed$(σdot0)_ωT$(ωT).dat") |> vec

#     ## Plotting 
#     fig = Figure(resolution=(1200, 1200), font="CMU Serif", fontsize=40, figure_padding = 30)
#     ax1 = Axis(fig[1, 1], xlabel=L"\Delta", ylabel="Frequency", title = L"\omega_T = %$(tTraj.ωT), \, \Phi_0 = %$(Φ0), \,\lambda = %$(λ), \, \dot{\sigma}_0 = %$(σdot0)")

#     hist!(ax1, Δs, bins = 24, normalization = :pdf, bar_labels = :values, label_formatter=x-> round(x, digits=2), label_size = 20, strokewidth = 0.5, strokecolor = (:black, 0.5), color = my_green)

#     pred_val = Δ_analytic(σdot0, Φ0, λ, ωmax)
#     vlines!(ax1, [pred_val], color = my_black, label = "Analytic - $(round(pred_val, digits = 3))")
#     vlines!(ax1, [mean(Δs)], color = my_vermillion, label = "Numerical Mean - $(round(mean(Δs), digits = 3))")

#     axislegend(ax1, labelsize = 40)
#     ylims!(0, nothing)
#     # save("dist_$(tTraj.ωT).pdf", fig)
#     fig
# end


tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ5_nmodes100000.jld2")
# Δs = Δ_thermal(σdot0, init_pos, Φ0, λ, system, tTraj)
# cd("data/Thermal/ωT$(ωT)")
# writedlm("delta_ωmax$(ωmax)_speed$(σdot0)_ωT$(ωT)_Pts$(length(Δs)).dat", Δs)
# cd("../../../")
Δs = readdlm("data/Thermal/ωT25.0/delta_ωmax10_speed60.0_ωT25.0.dat") |> vec

# ## Plotting 
fig = Figure(resolution=(1200, 1200), font="CMU Serif", fontsize=40, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\Delta", ylabel="Frequency", title = L"\omega_T = %$(tTraj.ωT), \, \Phi_0 = %$(Φ0), \,\lambda = %$(λ), \, \dot{\sigma}_0 = %$(σdot0)")

hist!(ax1, Δs, bins = 24, normalization = :pdf, bar_labels = :values, label_formatter=x-> round(x, digits=2), label_size = 20, strokewidth = 0.5, strokecolor = (:black, 0.5), color = my_green)

pred_val = Δ_analytic(σdot0, Φ0, λ, ωmax)
vlines!(ax1, [pred_val], color = my_black, label = "Analytic - $(round(pred_val, digits = 3))")
vlines!(ax1, [mean(Δs)], color = my_vermillion, label = "Numerical Mean - $(round(mean(Δs), digits = 3))")

axislegend(ax1, labelsize = 40, position = :lt)
ylims!(0, nothing)
fig


# fig = Figure(resolution=(1600, 1200), font="CMU Serif", fontsize=40, figure_padding = 30)
# ax1 = Axis(fig[1, 1], xlabel= "Number of Values", ylabel=L"\Delta", title = L"\omega_T = %$(tTraj.ωT), \, \Phi_0 = %$(Φ0), \,\lambda = %$(λ), \, \dot{\sigma}_0 = %$(σdot0)")

# test = map(nPts -> mean(Δs[1:nPts]), 1:1000:length(Δs))
# append!(test, mean(Δs))
# scatter!(ax1, 1:length(test), test, markersize = 14)
# hlines!(ax1, [Δ_analytic(σdot0, Φ0, λ, ωmax)], linewidth = 3, color = my_black)

# ylims!(ax1, -0.1, 0.1)
# fig