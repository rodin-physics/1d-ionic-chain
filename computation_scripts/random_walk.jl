include("../src/main.jl")
using Random

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 4.0
ωT = 10.0

function speed_after_pass(init_speed)
    pred_mean = Δ_thermal_analytic(init_speed, Φ0, λ, ωmax, ωT)
    pred_var = Δ_thermal_variance(init_speed, Φ0, λ, ωmax, ωT)

    dist = Normal(pred_mean, √(pred_var))

    return √(init_speed^2 - (8*π^2/μ)*rand(dist))
end

function delta_after_pass(init_speed)
    pred_mean = Δ_thermal_analytic(init_speed, Φ0, λ, ωmax, ωT)
    pred_var = Δ_thermal_variance(init_speed, Φ0, λ, ωmax, ωT)

    dist = Normal(pred_mean, √(pred_var))

    return rand(dist)
end


function random_walk(nPasses, nParticles, init_speed, num_ind; filter_before = false)
    speed_data = Float64[]
    Δ_data = Float64[]

    speeds = repeat([init_speed], nParticles)

    @showprogress for _ in 1:nPasses
        # Check for speeds below capture speed 
        stuck = findall(x -> x <= √(8*π^2*Φ0/μ), speeds)
        deleteat!(speeds, stuck)

        # 
        if isempty(speeds)
            break
        end

        Δs = delta_after_pass.(speeds)
        
        if filter_before == false 
            append!(speed_data, speeds)
            append!(Δ_data, Δs)
        end

        # Filter out speeds for which the particle does not pass the chain mass
        below_cutoff = findall(x -> x < Φ0, ((μ/8/π^2) .* speeds.^2 ).- Δs)
        deleteat!(speeds, below_cutoff)
        deleteat!(Δs, below_cutoff)

        if filter_before == true
            append!(speed_data, speeds)
            append!(Δ_data, Δs)
        end

        # Final speeds
        speeds = sqrt.(speeds.^2 .- (8 * π^2 / μ) .* Δs)
    end
    writedlm("data/Thermal/random_walk$(num_ind).dat", (speed_data, Δ_data))
end

Threads.@threads for batch in 1:9
    random_walk(800, 30, 120, batch, filter_before = false)
end