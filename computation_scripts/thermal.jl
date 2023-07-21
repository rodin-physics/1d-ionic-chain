include("../src/main.jl")

## FULL TRAJECTORY
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
d = 60
τ = 200                             # Simulation time
δ = system.δ                        # Time step
α = 40                              # Distance between chain atoms
μ = 1
σ0 = [Int(4.5 * α)]

n_pts = τ / δ |> floor |> Int
ωT = 0.0
tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ1000_lmax300_modes50000.jld2")
mem = Inf

speeds = [120]
bias = 0.0

param_vals = vcat(map(ii -> [(20, 4, [ii]), (-20, 4, [ii])], speeds)...)

function full_traj(param)
    println(param)
    Φ0 = param[1]
    λ = param[2]
    σdot0 = param[3]
    if (
        !isfile(
            "data/Thermal/Single/Single_sigma0$(σ0)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(ωT)_tau$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ, bias = bias, threads = true)
        save_object(
            "data/Thermal/Single/Single_sigma0$(σ0)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(ωT)_tau$(τ).jld2",
            res,
        )
    end
end

## Extract multiple trajectories by initialising particle at different positions
function multi_traj(param, numTraj)
    println(param)
    Φ0 = param[1]
    λ = param[2]
    σdot0 = param[3]

    n_pts = floor(τ / δ) |> Int 
    final_σs = zeros(numTraj, n_pts)
    init_pos_range = range(0, numTraj - 1, step = 1)

    for init_pos in init_pos_range
        res = motion_solver(system, Φ0, λ, α, σ0 .+ (init_pos * α), σdot0, μ, tTraj, mem, τ, bias = bias, threads = true)

        final_σs[init_pos + 1,:] = res.σs 
    end

    writedlm("data/Thermal/Single/Single$(numTraj)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(ωT)_tau$(τ).dat", final_σs)
end

# full_traj.(param_vals)

multi_traj.(param_vals, 5)