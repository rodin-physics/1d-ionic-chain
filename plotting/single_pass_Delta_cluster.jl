include("../../src/main.jl")
## Parameters
α = 40
μ = 1
nChain = 3
ωmax = 10
Φ0 = 2.0
λ = 4.0

system = load_object("precomputed/systems/System_ωmax10_d60_l10000_τmax5.jld2")
ωTs = [10.0]

# Numerically calculate energy loss after single pass
function Δ_numeric(σ_dot, σ0, Φ0, λ, system, tTraj)
    δ = system.δ
    τ = 1.25 * (α / σ_dot)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ)

    σs = res.σs |> vec
    ρs = res.ρs

    # Find index closest to next midpoint
    chain_idx = searchsortedlast(ρs[:, 1], σs[1])
    mod_val = mod(σs[1], res.α)
    mob_final = findmin(abs.(σs .- (ρs[chain_idx+1, :] .+ mod_val)))[2]

    # Calculate final kinetic energy
    v_final = (σs[mob_final] - σs[mob_final-1]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

# Get distribution of energy loss values
function Δ_thermal(σ_dot, σ0, Φ0, λ, system, tTraj)

    function get_tTraj(ind, tTraj)
        return ThermalTrajectory(
            tTraj.ωmax,
            tTraj.δ,
            tTraj.ρHs[ind:ind+nChain, :],
            tTraj.ωT,
        )
    end
    nPts = length(1:(size(tTraj.ρHs)[1]-nChain))
    res = zeros(nPts)
    p = Progress(nPts)
    Threads.@threads for ii = 1:(size(tTraj.ρHs)[1]-nChain)
        res[ii] = Δ_numeric(σ_dot, σ0, Φ0, λ, system, get_tTraj(ii, tTraj))
        next!(p)
    end

    return res
end


function Δ_distribution(speed, temp)
    println(speed)
    init_pos = 2.5 * α
    tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT$(temp)_τ5_nmodes100000.jld2")
    Δs = Δ_thermal(speed, init_pos, Φ0, λ, system, tTraj)

    cd("data/Thermal/ωT$(temp)")
    writedlm("delta_ωmax$(ωmax)_speed$(speed)_ωT$(temp).dat", Δs)
    cd("../../../")
end

function Δ_per_temp(dir_name)
    filenames = filter(x -> x[1] !== '.', readdir(dir_name))
    speeds = Float32[]
    Δs = Float32[]
    sdevs = Float32[]
    for file in filenames
        # Obtain mean and error bars 
        data = readdlm(joinpath(dir_name, file)) |> vec
        append!(Δs, mean(data))
        append!(sdevs, std(data))

        # Parse the filename and get the speed 
        names = filter(x -> x[1] == 's', split(file, "_"))
        speed = parse(Float64, chopprefix(names[1], "speed"))
        append!(speeds, speed)
    end
    return (speeds, Δs, sdevs)
end


# for ωT in ωTs 
#     speeds = range(16.0, 18.0, step = 2.0)
#     map(x -> Δ_distribution(x, ωT), speeds)
# end


## Plotting 
fig = Figure(
    resolution = (1600, 1200),
    font = "CMU Serif",
    fontsize = 40,
    figure_padding = 30,
)
ax = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}_0",
    ylabel = L"\Delta",
    title = L"\Phi_0 = %$(Φ0), \,\lambda = %$(λ)",
)

ωTs = reverse([0.0, 2.0, 5.0, 10.0, 25.0])
speeds = range(10.0, 80.0, step = 1.0)
Δ_analytics = Δ_analytic.(speeds, Φ0, λ, ωmax)
lines!(speeds, Δ_analytics, color = my_black, linewidth = 5, label = "Analytic")

colors = reverse([my_blue, my_green, my_vermillion, my_red, my_yellow])
for ωT in ωTs
    dir_name = joinpath(pwd(), "data/Thermal/", "ωT$(ωT)/")
    (xs, ys, errors) = Δ_per_temp(dir_name)
    scatter!(ax, xs, ys, markersize = 20, color = colors[findfirst(x -> x == ωT, ωTs)])
    lines!(
        ax,
        xs,
        ys,
        linewidth = 4,
        color = colors[findfirst(x -> x == ωT, ωTs)],
        label = L"\omega_T = %$(ωT)",
    )
    # errorbars!(ax, xs, ys, errors, whiskerwidth = 10, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4.0)
    # vlines!([sqrt(2 * 8 * pi^2 / μ)], color = colors[findfirst(x -> x == ωT, ωTs)], linewidth =4)
end

axislegend(position = :rb)
xlims!(ax, 16, 80)
ylims!(ax, -0.02, 0.05)
fig
