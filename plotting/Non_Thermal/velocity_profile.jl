include("../../src/main.jl")
# Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 20.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")

# Trajectory for a single pass
function single_pass(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.0 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    return res
end


# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

# Maximum speed deviation
function speed_deviation(σdot, Φ0, μ)
    return sqrt(σdot^2 - (8*π^2*Φ0/μ))
end

## Plotting
fig = Figure(resolution=(1600, 800), font="CMU Serif", fontsize=24, figure_padding = 30)

ax1 = Axis(fig[1, 1], xlabel=L"\tau", ylabel=L"\sigma", title = "Particle Position")
ax2 = Axis(fig[2, 1], xlabel=L"\tau", ylabel=L"\dot{\sigma}", title = "Velocities")
ax3 = Axis(fig[1, 2], xlabel=L"\dot{\sigma}", ylabel=L"\dot{\sigma}_0 - \dot{\sigma}_\mathrm{ext}", title = "Velocity Fluctuations")
ax4 = Axis(fig[2, 2], xlabel=L"\tau", ylabel=L"\dot{\sigma}", title = "Velocity Over Time")

## Computation
σdot0 = 50
xs = range(-α/σdot0/2, α/σdot0/2, length = 10000)
res_rep = single_pass(σdot0, 5.5 * α, Φ0, λ, system)
res_att = single_pass(σdot0, 5.5 * α, -Φ0, λ, system)

## Particle Position
hlines!(ax1, [6 * α], linestyle = :dash, linewidth = 2, color = my_black)
vlines!(ax1, [res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))]], color = my_vermillion, linestyle = :dash, linewidth = 2)
vlines!(ax1, [res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))]], color = my_blue, linestyle = :dash, linewidth = 2)
lines!(ax1, res_rep.τs, vec(res_rep.σs), color = my_vermillion, linewidth = 4, label = "Repulsive")
lines!(ax1, res_att.τs, vec(res_att.σs), color = my_blue, linewidth = 4, label = "Attractive")


## Velocities
(xs_rep, speed_rep) = particle_speed(res_rep)
(xs_att, speed_att) = particle_speed(res_att)

# Numerical particle velocities
vlines!(ax2, [res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))]], color = my_vermillion, linestyle = :dash, linewidth = 2)
vlines!(ax2, [res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))]], color = my_blue, linestyle = :dash, linewidth = 2)
lines!(ax2, xs_rep, speed_rep, linewidth = 4, color = my_vermillion)
lines!(ax2, xs_att, speed_att, linewidth = 4, color = my_blue)

# Analytical particle velocities
lines!(ax2, xs .+ res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))], sqrt.(σdot0^2 .- 8*π^2*U_profile.(xs, Φ0, λ/σdot0)), linewidth = 4, color = my_black)
lines!(ax2, xs .+ res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))], sqrt.(σdot0^2 .- 8*π^2*U_profile.(xs, -Φ0, λ/σdot0)), linewidth = 4, color = my_black, label = "Analytic")


## Speed Fluctuations
speeds = range(sqrt(8*Φ0*π^2/μ), 150, step = 0.5)
speeds2 = range(0, 150, step = 0.5)
vlines!(ax3, [sqrt(Φ0 * 8 * π^2)], color = my_red, linewidth = 2, label = "Capture speed")
lines!(ax3, speeds, speeds .- speed_deviation.(speeds, Φ0, μ), color = my_vermillion, linewidth = 4)
lines!(ax3, speeds2, speeds2 .- speed_deviation.(speeds2, -Φ0, μ), color = my_blue, linewidth = 4)


# Full velocity profile
data_rep = load_object("data/Non_Thermal/Single_sigma0[220]_sigmadot0[120]_MemInf_lambda4_Phi20_mu1_d60_bias0.0_omegaTnothing_tau125.jld2")
data_att = load_object("data/Non_Thermal/Single_sigma0[220]_sigmadot0[120]_MemInf_lambda4_Phi-20_mu1_d60_bias0.0_omegaTnothing_tau125.jld2")
hlines!(ax4, [sqrt(Φ0 * 8 * π^2)], color = my_red, linewidth = 2, label = "Capture speed")
hlines!(ax4, [0.0], color = my_black, linewidth = 2)
(xs_rep, speed_rep) = particle_speed(data_rep)
(xs_att, speed_att) = particle_speed(data_att)
lines!(ax4, xs_rep[1:2:end], speed_rep[1:2:end], linewidth = 4, color = my_vermillion)
lines!(ax4, xs_att[1:2:end], speed_att[1:2:end], linewidth = 4, color = my_blue)


## Labels and axes limits
Label(fig[1, 1, TopLeft()], "(a)", fontsize = 30)
Label(fig[2, 1, TopLeft()], "(b)", fontsize = 30)
Label(fig[1, 2, TopLeft()], "(c)", fontsize = 30)
Label(fig[2, 2, TopLeft()], "(d)", fontsize = 30)

xlims!(ax1, 0.0, 1 * (α / σdot0))
ylims!(ax1, 220, 264)
xlims!(ax2, 0.0, 1 * (α / σdot0))
xlims!(ax3, 0.0, 150)
xlims!(ax4, 0.0, 125)
axislegend(ax1, labelsize = 22, position = :rb)
axislegend(ax2, labelsize = 22, position = :rt)
axislegend(ax3, labelsize = 22, position = :rb)
axislegend(ax4, labelsize = 22, position = :rt)

# save("Velocity_Profile.pdf", fig)
fig
