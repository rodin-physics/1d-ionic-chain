include("../../src/main.jl")

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 1.0
τ = 500
δ = ((1 / ωmax) / 60)
τs = δ .* (1:Int(τ / δ))

# Box parameters
box_size = 100
left_boundary = 1.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)

ωTs = [1.0, 2.0, 5.0, 10.0, 25.0, 100.0, 250.0]
colors = [my_blue, my_sky, my_green, my_yellow, my_orange, my_vermillion, my_red, my_black]


## Mean squared deviation of particle trajectories 
function particle_RMSD_variable_box(data::Matrix{Float64}, box)
    num_P = size(data)[1]
    σs_final = copy(data)

    for particle in 1:num_P
        τ_ids = findall(x -> x < box[1] || x > box[2], σs_final[particle,:])

        for τ_id in τ_ids
            σs_final[particle, (τ_id+1):end] = 2 * σs_final[particle, τ_id] .- σs_final[particle, (τ_id+1):end]
        end

        σs_final[particle,:] = σs_final[particle,:] .- σs_final[particle,1]
    end

    return vec(std(σs_final, dims = 1))
end


## Plotting 
fig = Figure(resolution = (1400, 1000), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"Evolution time $\tau$", ylabel = L"RMS Deviation $\sqrt{\sum_i (\sigma_i - \sigma_0)^2 / 20}$", xscale = log10, yscale = log10, xgridvisible = false, ygridvisible = false,)

start_ind = 30
end_ind = 300

start_ind2 = 550

# Zero bias data
for ωT in ωTs 
    data = readdlm("data/Thermal/Thermalization_Single/SingleParticle_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem100.dat")
    rmsd_vals = particle_RMSD_variable_box(data, box)
    
    # Plot particle rmsd over time (avoid zeros for log10 scale to work)
    lines!(ax1, τs[start_ind:100:end], rmsd_vals[start_ind:100:end], color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = "$(Int(ωT))")

    # Ballistic curves
    if ωT > 2
        lines!(ax1, τs[start_ind2:100:end], ((rmsd_vals[start_ind2] / τs[start_ind2]) .* (τs[start_ind2:100:end]).^(1)) , linewidth = 3, linestyle = :dash, color = colors[findfirst(x -> x == ωT, ωTs)])
    end
end 

# Constant force curve
lines!(ax1, τs[start_ind:100:end_ind], (80 .* (τs[start_ind:100:end_ind]).^(2)) , linewidth = 2.5, linestyle = :dash, color = my_black)

xlims!(ax1, 0.05, 500)
ylims!(ax1, 0.1, nothing)
axislegend(L"\omega_T", position = :lt)

# Legend entries 
leg_entries = [LineElement(color = my_black, linewidth = 5), LineElement(color = my_black, linewidth = 4, linestyle = :dash), LineElement(color = my_black, linestyle = :dash, linewidth = 2.5)]

Legend(fig[1,1],
    leg_entries,
    ["Numerics", "Ballistic", "Constant force"],
    tellheight = false,
    tellwidth = false,
    patchsize = (40, 20),
    halign = 0.145,
    valign = 0.985
    ) 

fig
