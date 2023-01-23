include("../src/main.jl")
system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")

τ = 1000        # Simulation time
δ = system.δ    # Time step
α = 10          # Distance between chain atoms
μ = 1           # Particle mass
Φ0 = 7          # Potential magnitude                    
λ = 1           # Potential width

box_size = 100              # Particle-confining box
left_boundary = 10.5 * α
right_boundary = left_boundary + box_size * α

nParticles = 25
τ0 = 10
bias = 0.0

Random.seed!(120)
σ0 =
    (right_boundary + left_boundary - α) / 2 .* ones(nParticles) +
    15 * α * randn(nParticles)
σdot0 = zeros(nParticles)

tTraj = [
    load_object("precomputed/rH/rH_ωmax10_d60_ωT0_τ1000_l300.jld2"),
    load_object("precomputed/rH/rH_ωmax10_d60_ωT1_τ1000_l300.jld2"),
    load_object("precomputed/rH/rH_ωmax10_d60_ωT2_τ1000_l300.jld2"),
    load_object("precomputed/rH/rH_ωmax10_d60_ωT5_τ1000_l300.jld2"),
    load_object("precomputed/rH/rH_ωmax10_d60_ωT10_τ1000_l300.jld2"),
    load_object("precomputed/rH/rH_ωmax10_d60_ωT25_τ1000_l300.jld2"),
    load_object("precomputed/rH/rH_ωmax10_d60_ωT100_τ1000_l300.jld2"),
    load_object("precomputed/rH/rH_ωmax10_d60_ωT250_τ1000_l300.jld2"),
]

Threads.@threads for tr in tTraj
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
        )
        save_object("data/Box/BoxMultiple_τ0$(τ0)_λ$(λ)_Φ0$(Φ0)_ωT$(tr.ωT)_τ$(τ).jld2", res)
    end
end
