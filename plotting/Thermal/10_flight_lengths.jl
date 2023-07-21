include("../../src/main.jl")
using CurveFit

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 1.0
τ = 500
δ = ((1 / ωmax) / 60)
τs = δ .* (1:Int(τ / δ))

ωTs = [2.0, 5.0, 10.0, 25.0, 100.0, 250.0] |> reverse
colors = [my_sky, my_green, my_yellow, my_orange, my_vermillion, my_red] |> reverse

box_size = 100
left_boundary = 1.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)

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


## Check if a sign change has occured in an array 
function sign_change(arr)
    if length(arr) == 0
        error("Array is empty")
    end

    return !all( ==(sign(arr[1])), sign.(arr))
end


## Flight distances for every mobile particle 
function flight_distances_single(data::Matrix{Float64}; box = (-Inf, Inf))
    num_P = size(data, 1)
    ρHs = (1:105) .* α
    σs = traj_unfold_single(data, box)
    speeds = (σs[:, 2:end] .- σs[:, 1:end-1]) ./ δ

    # Initialise arrays
    curr_lens = zeros(num_P)
    dists = [Int64[] for _ in 1:num_P]
    prev_hops = ones(Int64, num_P)
    times = [Float64[] for _ in 1:num_P]

    for ind in 1:(length(τs)-1)
        
        # Determine position of mobile particles in relation to chain
        curr = searchsortedfirst.(Ref(ρHs), σs[:, ind])
        nxt = searchsortedfirst.(Ref(ρHs), σs[:, ind+1])

        # Compare initial and final positions
        particles = findall(x -> x != 0, (nxt .- curr))

        # A hopping event has occured
        if !isempty(particles)
            for part in particles 
                if curr_lens[part] == 0 || (sign_change(speeds[part, prev_hops[part]:ind]) == false)
                    
                    curr_lens[part] += 1

                # Particle did get stuck between the last hop event and now 
                elseif (τs[ind]-τs[prev_hops[part]] > 0.1)
                # else
                    # Record time it took for next hop to occur           
                    push!(times[part], τs[ind] - τs[prev_hops[part]])

                    # Record flight distance
                    push!(dists[part], curr_lens[part])

                    # Reset flight distance
                    curr_lens[part] = 1

                    # Update the latest hop event time index
                    prev_hops[part] = ind

                end
            end
        end
    end
    return dists
end


## Plotting density plots of flight lengths + times 

fig = Figure(resolution = (1200, 800), fontsize = 30, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"Flight Length $\ell$", ylabel = L"Probability Density $P(\ell)$", xgridvisible = false, ygridvisible = false, xscale = Makie.log10, yscale = Makie.log10,
xminorticksvisible = true,  xminorticks = IntervalsBetween(9), xticks = 10 .^(0:1:2))


for ωT in ωTs
    # Load data 
    data = readdlm("data/Thermal/Thermalization_Single/SingleParticle_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem100.dat")

    ## Array of flight lengths
    flight_lengths = vcat(flight_distances_single(data, box = box)...)

    if !isempty(flight_lengths)
        hist_fit = fit(Histogram, flight_lengths, 1:1:100)
        hist_fit = normalize(hist_fit, mode = :pdf)

        scatter!(
            ax1,
            ((hist_fit.edges[1])[1:end-1]),
            (hist_fit.weights),
            markersize = 18,
            color = colors[findfirst(x -> x == ωT, ωTs)]
        )

        # Numerical fit to power law 
        # Filter out small values to avoid numerical errors 
        inds = findall(x -> x >= 1e-6, hist_fit.weights)

        xs = ((hist_fit.edges[1])[1:end-1])[inds]
        ys = (hist_fit.weights)[inds]
        numerical_fit = curve_fit(PowerFit, float.(xs), float.(ys))

        # Plot for x range 
        xs = range(1, 1000, length = 1000)
        lines!(ax1, xs, numerical_fit.(xs), linewidth = 3, color =colors[findfirst(x -> x == ωT, ωTs)], linestyle = :dash)
    end
end


# Legend entries 
leg_entries = [MarkerElement(marker = :circle, color = cs) for cs in reverse(colors)]
line_entry = [LineElement(color = my_black, linestyle = :dash, linewidth = 3.5)]

Legend(fig[1,1],
       [leg_entries, line_entry],
       [string.(Int.(reverse(ωTs))), ["Power law fits"]],
       [L"\omega_T", " "],
       halign = :right,
       valign = :top,
       tellheight = false,
       tellwidth = false,
       titleposition = :left,
       orientation = :horizontal,
       margin = (10, 10, 10, 10),
       titlegap = 15,
       nbanks = 2,
       patchsize = (40, 20)
       )



xlims!(ax1, nothing, 101)
ylims!(ax1, 0.001, 1.2)

fig