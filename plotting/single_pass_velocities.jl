include("../../src/main.jl")

## Parameters
α = 40
μ = 1
nChain = 5
ωmax = 10
Φ0 = 2.0
λ = 4.0
# Load memory kernel and thermal trajectory
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT0.0_τ10.jld2")

# Trajectory for a single pass
function single_pass(σ_dot, σ0, Φ0, λ, system, tTraj)
    τ = 1.25 * (α / σ_dot)
    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads = true)
    return res
end

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii = 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end


# Computation 
σdot0 = 50.0
init_pos = range(2.5, 990.5, step = 50.0) .* α

# Plotting 
fig = Figure(
    resolution = (1200, 1200),
    font = "CMU Serif",
    fontsize = 40,
    figure_padding = 30,
)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\tau",
    ylabel = "Speed with offset",
    title = L"\omega_T = %$(tTraj.ωT), \, \Phi_0 = %$(Φ0), \,\lambda = %$(λ), \, \dot{\sigma}_0 = %$(σdot0)",
)

for ii in eachindex(init_pos)
    data = single_pass(σdot0, init_pos[ii], Φ0, λ, system, tTraj)
    times, speeds = particle_speed(data)
    lines!(ax1, times, speeds .+ (ii * 0.1), color = my_blue, linewidth = 2)
end

fig
