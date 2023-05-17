include("../../src/main.jl")
using Random

## Parameters
α = 40
μ = 1
ωmax = 10
nChain = 3
Φ0 = 2.0
λ = 4.0

## Function to parse the datafiles and obtain mean and standard deviations
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

speeds = range(15, 80, step = 5.0)
ωTs = [0.0, 2.0, 5.0, 10.0, 25.0]
colors = [my_blue, my_sky, my_green, my_orange, my_vermillion]

fig = Figure(
    resolution = (1200, 1300),
    fontsize = 36,
    figure_padding = 30,
)

ax1 = Axis(
    fig[1, 1],
    ylabel = L"Mean Loss $\langle \Delta \rangle$",
    xticklabelsvisible = false,
    xgridvisible = false,
    ygridvisible = false
)

ax2 = Axis(
    fig[2, 1],
    xlabel = L"Velocity midway between lattice sites $\dot{\sigma}$",
    ylabel = L"Loss fluctuation $\sqrt{\langle \Delta^2_\mathrm{fluc} \rangle}$",
    xgridvisible = false,
    ygridvisible = false
)

vlines!(ax1, [√(8*pi^2 * Φ0)], linestyle = :dash, linewidth = 5, color = my_red)
vlines!(ax2, [√(8*pi^2 * Φ0)], linestyle = :dash, linewidth = 5, color = my_red, label = "Capture Speed")

range_len = 100
x_range_nt = range(7, 80, length = range_len)


for ωT in ωTs 
    mean_data = readdlm("data/Thermal/MeanΔ_ωT$(ωT).dat")
    var_data = readdlm("data/Thermal/VarΔ_ωT$(ωT).dat")

    lines!(ax1, mean_data[1,:], mean_data[2,:], color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = "$(Int(ωT))")
    new_var = (var_data[2,:])
    lines!(ax2, var_data[1,:], sqrt.(new_var), color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4)

    (speeds, Δs, sdevs) = Δ_per_temp("data/Thermal/ωT$(ωT)/")

    scatter!(ax1, speeds, Δs, color = colors[findfirst(x -> x == ωT, ωTs)], markersize = 16)
    scatter!(ax2, speeds, sdevs, color = colors[findfirst(x -> x == ωT, ωTs)], markersize = 16)

end

# Non-thermal results
lines!(ax1, x_range_nt, Δ_analytic.(x_range_nt, Φ0, λ, ωmax), color = my_black, linewidth = 3, label = "Nonthermal")

# Legend entries
temp_color = [LineElement(color = color, linewidth = 5) for color in colors]
res_type = [LineElement(color = my_black, linewidth = 5)]

Legend(fig[1,1],
    [temp_color, res_type],
    [string.(Int.(ωTs)), ["Nonthermal"]],
    [L"\omega_T", nothing],
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :top,
    titleposition = :left,
    margin = (0, 10, 0, 10),
    width = 260,
    titlegap = -30,
    rowgap = -2,
    ) 


axislegend(ax2, position = :rt, patchsize = (35, 35))

xlims!(ax1, 0.0, nothing)
xlims!(ax2, 0.0, nothing)

ylims!(ax1, 0.0, 0.034)
ylims!(ax2, 0.0, nothing)

fig

