include("../src/main.jl")
using Random

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 1.0
τ = 1000
ωTs = [1.0, 2.0, 5.0, 10.0, 25.0, 100.0, 250.0]


# Load memory kernel and thermal trajectory
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
tTraj_paths = ["precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ1000_lmax300_modes50000.jld2" for ωT in ωTs] 


# Particle-confining box
box_size = 250              
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)
nParticles = 25

## Computation for 25 particles simultaneously in a box
mem = 10
for path in tTraj_paths
    σ0 = (right_boundary + left_boundary - α) / 2 .* ones(nParticles) + 10 * α * rand(nParticles)
    σdot0 = zeros(nParticles)

    tTraj = load_object(path)
    res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ; threads = true, box = box)

    save_object("data/Thermal/Thermalization/Particle$(nParticles)_Φ$(Φ0)_λ$(λ)_ωT$(tTraj.ωT)_τ$(τ)_mem$(mem).jld2", res)
end


### Computation for 20 single-particle trajectories

## Obtain mobile particle positions over τ for single-particle trajectories
function get_single_particle_dist(ωT, box, numP)
    # Initialise final positions array 
    n_pts = τ / ((1 / ωmax) / 60) |> Int
    final_pos = zeros(numP, n_pts)

    # Load ThermalTrajectory
    full_tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ1000_lmax300_modes50000.jld2")

    # Solve for each particle over different parts of the chain
    @showprogress for part in 1:numP
        curr_section = (1:102) .+ (10 * (numP - 1))
        curr_tTraj = ThermalTrajectory(full_tTraj.ωmax, full_tTraj.δ, full_tTraj.ρHs[curr_section, :], full_tTraj.ωT)

        σ0 = (box[2] + box[1] - α) / 2 + 10 * α * rand()
        σdot0 = [0.0]

        res = motion_solver(system, Φ0, λ, α, [σ0], σdot0, μ, curr_tTraj, mem, τ, threads = true, box = box)

        # Save only the mobile particle position 
        final_pos[part,:] = vec(res.σs) 
    end

    return final_pos
end

# Particle-confining box
box_size = 100              
left_boundary = 1.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)
nParticles = 20

mem = 10
for ωT in ωTs
    single_traj = get_single_particle_dist(ωT, box, 20)
    writedlm("data/Thermal/Thermalization_Single/SingleParticle_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem$(mem).dat", single_traj)
end

