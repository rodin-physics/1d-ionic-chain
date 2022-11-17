include("../../src/main.jl")

## Parameters 
α = 40
μ = 1
ωmax = 10
nChain = 3
Φ0 = 2.0
λ = 10.0

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



colors = reverse([my_blue, my_green, my_vermillion, my_red, my_yellow])
ωTs = reverse([0.0, 2.0, 5.0, 10.0, 25.0])

speeds = range(15, 30, length = 50)

## Plotting 
fig = Figure(
    resolution = (1600, 1200),
    font = "CMU Serif",
    fontsize = 40,
    figure_padding = 30,
)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}",
    ylabel = L"\Delta",
    title = L"\Phi_0 = %$(Φ0), \,\lambda = %$(λ), \, \alpha = %$(α)",
)

for ωT in ωTs
    nPts = length(speeds)
    analytic_Δs = zeros(nPts)
    p = Progress(nPts)
    Threads.@threads for ii in eachindex(speeds)
        analytic_Δs[ii] = Δ_thermal_analytic(speeds[ii], Φ0, λ, ωmax, ωT)
        next!(p)
    end
    lines!(
        speeds,
        analytic_Δs,
        linewidth = 5,
        color = colors[findfirst(x -> x == ωT, ωTs)],
        label = L"\omega_T = %$(ωT)",
    )
end

# for ωT in ωTs 
#     dir_name = joinpath(pwd(), "data/Thermal/", "ωT$(ωT)/")
#     (xs, ys, errors) = Δ_per_temp(dir_name)
#     scatter!(ax1, xs, ys, markersize = 20, color = colors[findfirst(x -> x == ωT, ωTs)])

#     # errorbars!(ax, xs, ys, errors, whiskerwidth = 10, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4.0)
# end



# lines!(ax1, speeds, Δ_analytic.(speeds, Φ0, λ, ωmax), label = "Non-Thermal", linewidth = 4 , color = my_black)

axislegend(ax1, position = :rt)
# xlims!(ax1, 16, 80)
# ylims!(ax1, -0.05, 0.05)
# save("Thermal_lambda10_low_speeds.pdf", fig)
fig
