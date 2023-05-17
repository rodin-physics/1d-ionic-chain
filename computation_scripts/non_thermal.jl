include("../src/main.jl")
using Peaks

## FULL TRAJECTORY
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
d = 60
τ = 500                             # Simulation time
δ = system.δ                        # Time step
α = 40                              # Distance between chain atoms
μ = 1
σ0 = [Int(4.5 * α)]

n_pts = τ / δ |> floor |> Int
nChain = 800
ρHs = zeros(nChain, n_pts)
tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
mem = 10

speeds = range(14, 50, step = 2)
bias = 0.01

param_vals = vcat(map(ii -> [(2, 4, [ii])], speeds)...)

function full_traj(param)
    println(param)
    Φ0 = param[1]
    λ = param[2]
    σdot0 = param[3]
    if (
        !isfile(
            "data/Non_Thermal/Single_sigma0$(σ0)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(nothing)_tau$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ, bias = bias, threads = true)
        particle_speeds = (res.σs[:, 2:end] .- res.σs[:, 1:end-1]) ./ system.δ |> vec
        pks, vals = findmaxima(particle_speeds)

        # Save particle speeds only 
        # writedlm("data/Non_Thermal/Single_sigma0$(σ0)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(nothing)_tau$(τ).dat", (res.τs[pks], vals))

        # Save full result
        save_object(
            "data/Non_Thermal/Single_sigma0$(σ0)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(nothing)_tau$(τ).jld2",
            res,
        )
    end
end

full_traj.(param_vals)