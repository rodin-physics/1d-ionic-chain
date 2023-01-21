include("../src/main.jl")
## SYSTEM PARAMETERS
α = 40
μ = 1
Φ0 = 1
λ = 1
ωmax = 10
ωTs = [1e-5, 1, 5, 10]

## ANALYTIC MEAN
Φ0 = 1

nPts = 200
vMin = 1
vMax = 30
vs = range(vMin, vMax, length = nPts) |> collect
vs_idx = reduce(
    vcat,
    [collect(eachindex(vs))[n:Threads.nthreads():end] for n = 1:Threads.nthreads()],
)

for ωT in ωTs
    println("ωT = $(ωT)")
    if (!isfile("data/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2"))
        result = zeros(nPts)
        pr = Progress(nPts)
        Threads.@threads for ii in vs_idx
            result[ii] = Δ_thermal_analytic(vs[ii], Φ0, λ, ωmax, ωT)
            next!(pr)
            save_object("data/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2", result)
        end
    end
end

## SINGLE PASS LOSS
system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")
δ = system.τs[2] - system.τs[1]

# HIGH SPEED
# Φ0 = 1
# σ_dots = 15:3:30 |> collect
# # Use the longest τ to make sure that the homogeneous trajectory
# # is sufficient even for the slowest particle
# τ = 1.5 * (α / minimum(σ_dots))
# n_pts = τ / δ |> floor |> Int
# τs = range(0, τ, length = n_pts)
# ωTs = [1e-5, 1, 5, 10]

nMasses = 10000
ωmax = 10
θs = π .* (1:nMasses) / nMasses
ωs = ω.(ωmax, θs)

numBatch = 20
batchSize = 250

# for ωT in ωTs
#     println("ωT = $(ωT)")
#     res = zeros(length(σ_dots), numBatch * batchSize)
#     if (!isfile("data/NumericalLossFast_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_μ$(μ)_ωT$(ωT).jld2"))
#         for batch = 1:numBatch
#             # Create a thermal trajectory
#             ϕs = 2 * π * rand(nMasses)
#             ζs = ζq.(ωs, ωT)
#             ρHs = [ρH(t, ζs, ϕs, ωs, 0:300) |> real for t in τs]
#             ρHs = reduce(hcat, ρHs)
#             tTraj = ThermalTrajectory(system.ωmax, δ, ρHs, ωT)
#             println("Trajectory generated for batch number $(batch)")
#             # Use the same trajectory to compute losses at different speeds
#             @showprogress for jj in eachindex(σ_dots)
#                 # println("Batch $(batch), speed index $(jj)")
#                 Threads.@threads for ii = 1:batchSize
#                     σ_dot = σ_dots[jj]
#                     σ0 = (10.5 + ii) * α
#                     # Use the actual speed for the computation length
#                     τ_calc = 1.5 * (α / σ_dot)
#                     sol = motion_solver(
#                         system,
#                         Φ0,
#                         λ,
#                         α,
#                         [σ0],
#                         [σ_dot],
#                         μ,
#                         tTraj,
#                         Inf,
#                         τ_calc,
#                     )
#                     σs = sol.σs |> vec
#                     mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
#                     v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
#                     Δ = (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
#                     res[jj, ii+(batch-1)*batchSize] = Δ
#                 end
#             end
#             println("Batch $(batch) complete")
#         end

#         save_object(
#             "data/NumericalLossFast_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_μ$(μ)_ωT$(ωT).jld2",
#             (σ_dots, res),
#         )
#     end
# end
# nMasses = 1000
# ϕs = 2 * π * rand(nMasses)
# ζs = ζq.(ωs, ωT)
# @time ρHs = [ρH(t, ζs, ϕs, ωs, 1:20) |> real for t in τs]
# ρHs = reduce(hcat, ρHs)
# zz = load_object("data/NumericalLossFast_Φ01_λ1_ωmax10_μ1.jld2")
# zz[3]
# r = zeros(num_Batch * batchSize)

σ_dot = 80
ωT = 1
Φ0 = 5

τ = 1.5 * (α / σ_dot)
n_pts = τ / δ |> floor |> Int
τs = range(0, τ, length = n_pts)
r = zeros(numBatch * batchSize)
Threads.@threads for b = 1:numBatch
    ϕs = 2 * π * rand(nMasses)
    ζs = ζq.(ωs, ωT)
    ρHs = [ρH(t, ζs, ϕs, ωs, 0:300) |> real for t in τs]
    ρHs = reduce(hcat, ρHs)
    tTraj = ThermalTrajectory(system.ωmax, δ, ρHs, ωT)

    @showprogress for ii = 1:batchSize
        σ0 = 180 + ii * α
        res = motion_solver(
            system,
            Φ0,
            λ,
            α,
            [σ0],
            [σ_dot],
            μ,
            tTraj,
            Inf,
            τ,
            threads = false,
        )
        σs = res.σs |> vec
        mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
        v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
        r[ii+(b-1)*batchSize] = (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
    end
end

save_object(
    "data/NumericalLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_v$(σ_dot)_μ$(μ)_ωT$(ωT).jld2",
    (σ_dot, r),
)


# # # (1 / 2 / (2 * π)^2 * (30^2))

# @time Δ_thermal_analytic(55, 3, λ, ωmax, 1)
# @time Δ_analytic(100, 5, λ, ωmax)

# # 1
