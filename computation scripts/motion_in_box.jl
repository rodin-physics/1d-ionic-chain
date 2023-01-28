include("../src/main.jl")
system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")

τ = 100        # Simulation time
δ = system.δ    # Time step
α = 10          # Distance between chain atoms
μ = 1           # Particle mass
Φ0 = 15          # Potential magnitude                    
λ = 1           # Potential width

box_size = 100              # Particle-confining box
left_boundary = 10.5 * α
right_boundary = left_boundary + box_size * α

nParticles = 2
τ0 = 10
bias = 0.0

Random.seed!(120)
σ0 =
    (right_boundary + left_boundary - α) / 2 .* ones(nParticles) +
    15 * α * randn(nParticles)
σdot0 = zeros(nParticles)

tTraj = [
    # "precomputed/rH/rH_ωmax10_d60_ωT0_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT1_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT2_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT5_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT10_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT15_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT25_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT50_τ1000_l300.jld2",
    # "precomputed/rH/rH_ωmax10_d60_ωT100_τ1000_l300.jld2",
    "precomputed/rH/rH_ωmax10_d60_ωT250_τ1000_l300.jld2",
]

for traj in tTraj
    tr = load_object(traj)
    if (!isfile("data/Box/BoxMultiple_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_ωT$(tr.ωT)_τ$(τ).jld2"))
        # Keep only a portion of the thermal trajectory
        tr = ThermalTrajectory(tr.ωmax, tr.δ, tr.ρHs[1:box_size+20, :], tr.ωT)
        println("ωT = $(tr.ωT)")
        res = motion_solver(
            system,
            Φ0,
            λ,
            α,
            σ0,
            σdot0,
            μ,
            tr,
            τ0,
            τ;
            box = (left_boundary, right_boundary),
            threads = true,
        )
        save_object(
            "data/Box/TEST_BoxMultiple_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_ωT$(tr.ωT)_τ$(τ).jld2",
            res,
        )
    end
end
