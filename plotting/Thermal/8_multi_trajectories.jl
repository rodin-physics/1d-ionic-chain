include("../../src/main.jl")

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 1.0
τ = 500
δ = (1 / ωmax) / 60
n_pts = floor(τ / δ) |> Int
τs = δ .* (1:n_pts)

# Box parameters
box_size = 100
left_boundary = 1.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)


drop = 1
spacing = 100
numTraj = 20

## Unfold trajectory for single particle dataset
function traj_unfold_single(data::Matrix{Float64}, box)
    σs_final = copy(data)
    for particle in 1:size(data,1)
        τ_ids = findall(x -> x < box[1] || x > box[2], σs_final[particle,:])

        for τ_id in τ_ids
            σs_final[particle, (τ_id+1):end] = 2 * σs_final[particle, τ_id] .- σs_final[particle, (τ_id+1):end]
        end

        σs_final[particle,:] = σs_final[particle,:] .- σs_final[particle,1]
    end
    return σs_final
end


## Plotting 
fig = Figure(resolution = (1400, 1000), font = "CMU Serif", fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xticklabelsvisible = false, xgridvisible = false, ygridvisible = false)
ax2 = Axis(fig[2, 1], xlabel = L"Evolution time $\tau$", xgridvisible = false, ygridvisible = false, yticks = (-100:50:100))

ωT = 5.0
curr_color = my_green
data = readdlm("data/Thermal/Thermalization_Single/SingleParticle_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem100.dat")
σs = traj_unfold_single(data, box)

for ind in 1:numTraj
    lines!(ax1, τs[drop:spacing:end] .- τs[drop], (σs[ind, drop:spacing:end] .- σs[ind, drop]) ./ 10^2, color = curr_color, linewidth = 2)
end

ωT = 100.0
spacing = 600
curr_color = my_vermillion
data = readdlm("data/Thermal/Thermalization_Single/SingleParticle_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem100.dat")
σs = traj_unfold_single(data, box)

for ind in 1:numTraj
    lines!(ax2, τs[drop:spacing:end] .- τs[drop], (σs[ind, drop:spacing:end] .- σs[ind, drop]) ./ 10^2, color = curr_color, linewidth = 2)
end

xlims!(ax1, 0, 500)
xlims!(ax2, 0, 500)

ylims!(ax1, -1.2, 1.2)
ylims!(ax2, -120, 120)

Label(fig[1:2, 0], L"Displacement $(\sigma - \sigma_0) / 10^2$", rotation = pi/2)

text!(ax1, 10, -1.1, text = L"(a) $\omega_T = 5$", fontsize = 35)
text!(ax2, 10, 80, text = L"(b) $\omega_T = 100$", fontsize = 35)

fig