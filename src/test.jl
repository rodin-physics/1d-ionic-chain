include("main.jl")

r = load_object("data/Box/BoxMultiple_τ010_λ1_Φ07_ωT1.0_τ1000.jld2")

nParticles = length(σs[:, 1])
traj = zeros(size(r.σs))
for ii = 1:nParticles
    tr = r.σs[ii, :]
    for st = 2:lastindex(tr)
        if (tr[st] - tr[st-1]) > 0.9 * (right_boundary - left_boundary)
            tr[st:end] = tr[st:end] .- (right_boundary - left_boundary)
        elseif (tr[st] - tr[st-1]) < -0.9 * (right_boundary - left_boundary)
            tr[st:end] = tr[st:end] .+ (right_boundary - left_boundary)
        end
    end
    traj[ii, :] = tr .- tr[1]

end


fig = Figure(
    resolution = (1200, 800),
    fonts = (; math = "CMU Serif"),
    fontsize = 40,
    figure_padding = 30,
)

ax1 = Axis(fig[1, 1])

for ii = 1:nParticles
    lines!(ax1, r.τs, traj[ii, :])
end

fig
traj[] = traj[ii, :] .- traj[ii, 1]


std(traj[:, 200])
st = [std(traj[:, n]) for n = 1:length(r.τs)]


fig = Figure(
    resolution = (1200, 800),
    fonts = (; math = "CMU Serif"),
    fontsize = 40,
    figure_padding = 30,
)

ax1 = Axis(fig[1, 1])
lines!(ax1, (r.τs), (st))
lines!(ax1, (r.τs), 5.2 .* sqrt.(r.τs))
fig
# tr = r.σs[1, :]
# idx = findall(x -> x < left_boundary || x > right_boundary, tr)
# tr_unfold = tr



# @showprogress for n in idx
#     if tr[n] > right_boundary
#         tr_unfold[n+1:end] = 2 * tr_unfold[n] .- tr_unfold[n+1:end]
#     elseif tr[n] < left_boundary
#         tr_unfold[n+1:end] = 2 * tr_unfold[n] .- tr_unfold[n+1:end]
#     end
# end

# lines(r.τs, tr_unfold)



lines!(ax1, r.τs, r.σs[1, :])
scatter!(ax1, r.τs[idx], r.σs[1, idx])
lines!(ax1, r.τs, tr_unfold)
fig

5
# system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")

# τ = 20        # Simulation time
# δ = system.δ    # Time step
# α = 10          # Distance between chain atoms
# μ = 1           # Particle mass
# Φ0 = 7          # Potential magnitude                    
# λ = 1           # Potential width

# box_size = 10              # Particle-confining box
# left_boundary = 10.5 * α
# right_boundary = left_boundary + box_size * α

# nParticles = 10
# τ0 = 10
# bias = 0.0

# Random.seed!(120)
# σ0 =
#     (right_boundary + left_boundary - α) / 2 .* ones(nParticles) +
#     15 * α * randn(nParticles)
# σdot0 = zeros(nParticles)

# tr = load_object("precomputed/rH/rH_ωmax10_d60_ωT1_τ1000_l300.jld2")
# tr = ThermalTrajectory(tr.ωmax, tr.δ, tr.ρHs[1:box_size+20, :], tr.ωT)
# res = motion_solver(
#     system,
#     Φ0,
#     λ,
#     α,
#     σ0,
#     σdot0,
#     μ,
#     tr,
#     τ0,
#     τ;
#     box = (left_boundary, right_boundary),
#     threads = true,
# )

# res_TEST = motion_solver(
#     system,
#     Φ0,
#     λ,
#     α,
#     σ0,
#     σdot0,
#     μ,
#     tr,
#     τ0,
#     τ;
#     box = (left_boundary, right_boundary),
#     threads = false,
# )




# # println(res == res_TEST)
# res.σs == res_TEST.σs
# # motion_solver == motion_solver_TEST
# res.σs

# # struct SystemSolution
# #     ωmax::Float64               # Largest mode frequency
# #     μ::Float64                  # Mass of the mobile atoms
# #     τs::Vector{Float64}         # Time steps
# #     τ0::Float64                 # Memory
# #     α::Float64                  # Spacing between chain atoms
# #     Φ::Float64                  # Magnitude of the Gaussian potential
# #     λ::Float64                  # Standard deviation of the potential
# #     σs::Matrix{Float64}         # Positions of mobile atoms
# #     ρs::Matrix{Float64}         # Positions of the chain atoms
# #     bias::Float64               # Applied bias
# #     ωT::Union{Nothing,Float64}  # Temperature
# # end
# # @time reduce(+, [x^2 for x in 1:3])
# using BenchmarkTools

# m = 1
# r = zeros(300, m)
# # println(@benchmark [exp(-(x - y)^2) for x = 1:300, y = 1:m])


# # println()

# @btime (Threads.@threads for y = 1:m
#     r[:, y] = exp.(-((1:300) .- y) .^ 2)
# end)
# # [exp(-(x-y)^2) for x = 1:300, y = 1:m]
# # exp.(-(1:300 .- 1).^2)
# @btime ([exp(-(x - y)^2) for x = 1:300, y = 1:m])

# @btime findall(x -> x > 10, 1:20)
# @btime findall(x -> x > 20, 1:200)

# @btime (2 * ones(25) + 2 * ones(25))
# @btime (2 * (ones(25)  + ones(25) ))

res = load_object("data/Box/TEST_BoxMultiple_τ010_λ1_Φ015_ωT250.0_τ100.jld2")
tr = res.σs[1, :]

idx = findall(x -> tr[x] - tr[x-1] > 0.9 * (right_boundary - left_boundary), 2:length(tr))

for st = 2:lastindex(tr)
    if (tr[st] - tr[st-1]) > 0.9 * (right_boundary - left_boundary)
        tr[st:end] = tr[st:end] .- (right_boundary - left_boundary)
    elseif (tr[st] - tr[st-1]) < -0.9 * (right_boundary - left_boundary)
        tr[st:end] = tr[st:end] .+ (right_boundary - left_boundary)
    end
end

lines(res.τs, tr)
