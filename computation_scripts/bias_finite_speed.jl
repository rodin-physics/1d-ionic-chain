include("../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 2.0
λ = 4.0
τ = 1000
bias_val = 0.01

ωT = 0.0
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
speed_range = range(15, 60, step = 1.0)

# Box parameters
box_size = 250
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)

# Number of particles per batch
nParticles = 25
total_num = nParticles * 40
per_traj = 4
particles_per_traj = nParticles * per_traj

speed_array = zeros(5, total_num)

## Unfold particle trajectories given a confining box
function traj_unfold(data, box; periodic = false)
    σs_final = copy(data.σs)
    num_P = size(data.σs, 1)

    for particle in 1:num_P
        if periodic 
            for st in 2:lastindex(σs_final[particle, :])
                if (σs_final[particle, st] - σs_final[particle, st-1]) > 0.9 * (box[2] - box[1])
                    σs_final[particle, st:end] = σs_final[particle, st:end] .- (box[2] - box[1])
                elseif (σs_final[particle, st] - σs_final[particle, st-1]) < -0.9 * (box[2] - box[1])
                    σs_final[particle, st:end] = σs_final[particle, st:end] .+ (box[2] - box[1])
                end
            end
        else 
            τ_ids = findall(x -> x < box[1] || x > box[2], σs_final[particle,:])

            for τ_id in τ_ids
                σs_final[particle, (τ_id+1):end] = 2 * σs_final[particle, τ_id] .- σs_final[particle, (τ_id+1):end]
            end
        end
        σs_final[particle,:] = σs_final[particle,:] .- σs_final[particle,1]
    end

    return σs_final
end

# Thermal trajectory parameters
d = 60
τmax = 1000                     # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_modes = 5000                 # Number of chain masses for simulating ρ0

qa_s = 2 * pi .* (1:n_modes) / n_modes
ωs = ω.(ωmax, qa_s ./ 2)
lmax = 260
ε = reduce(hcat, [exp.(1im * 2 * π / n_modes .* (1:n_modes) * g) for g = 1:lmax])


# Get distribution of speeds for total_num particles at each time interval
function generate_final_speeds(speed_data)
    # Time interval at which data should be saved 
    num_interval = size(speed_data,1)
    save_time = Int(τ/(num_interval-1))

    for traj_ind in 1:(total_num ÷ particles_per_traj)
        ## Generate thermal trajectory 
        ζs = ζq.(ωs, ωT)
        ϕs = 2 * π * rand(n_modes)
    
        res = zeros(lmax, n_pts)
        pr = Progress(n_pts)
    
        Threads.@threads for n = 1:n_pts
            res[:, n] = ρH(n * δ, ζs, ϕs, ωs, ε) |> real
            next!(pr)
        end
        
        curr_traj = ThermalTrajectory(ωmax, δ, res, ωT)
    
        for run_ind in 1:(particles_per_traj ÷ nParticles)
            # Start nParticles with random speeds and positions within the box
            init_speeds = rand(speed_range, nParticles)
            init_pos = (box[2] + box[1] - α) / 2 .* ones(nParticles) + 10 * α * randn(nParticles)
    
            data = motion_solver(system, Φ0, λ, α, init_pos, init_speeds, μ, curr_traj, 10, τ; threads = true, box = box, periodic = true, bias = bias_val)
            
            # Add initial speeds to matrix
            curr_ind = (traj_ind - 1) * per_traj + run_ind
            curr_indices = (curr_ind * nParticles) - (nParticles - 1):(curr_ind * nParticles)
            speed_data[1, curr_indices] = init_speeds


            # Get the final speeds at each time interval
            data_unfold = traj_unfold(data, box, periodic = true)
            time_indices = vcat(findall(x -> x % save_time == 0, data.τs), lastindex(data.τs))

            for ii in 2:num_interval
                final_speeds = (data_unfold[:, time_indices[ii-1]] .- data_unfold[:, time_indices[ii-1]-1]) ./ δ

                speed_data[ii, curr_indices] = final_speeds
            end

        end
    
    end

    ## Save data 
    writedlm("Speeds_$(τ)τ_inter$(save_time)_α$(α)_Φ$(Φ0)_λ$(λ)_bias$(bias_val)_T$(ωT).dat", speed_data)
end

generate_final_speeds(speed_array)