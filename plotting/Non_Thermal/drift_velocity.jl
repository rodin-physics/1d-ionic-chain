include("../../src/main.jl")
# Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 2.0
λ = 4.0
bias = 0.01

# Bias drift calculation
xs = range(1,90, step = 0.001)
Δ_bias = Δ_transport.(xs, Φ0, λ, ωmax, α)

## Plotting
fig = Figure(resolution = (1800, 1000), font = "CMU Serif", fontsize = 32, figure_padding = 30)

ax1_top = fig[1,1] = Axis(fig, title = "Repulsive")
ax1 = fig[2,1] = Axis(fig, xlabel = L"\tau", ylabel = L"\dot{\sigma}")

ax2_top = fig[1,2] = Axis(fig, title = "Attractive")
ax2 = fig[2,2] = Axis(fig, xlabel = L"\tau")

ax3_top = fig[1,3] = Axis(fig)
ax3 = fig[2,3] = Axis(fig, xlabel = L"\bar{\Delta}", ylabel = L"\dot{\sigma}")


# Predicted drift velocities
hlines!(ax1, [α/n for n in 1:20], linewidth = 2, linestyle = :dash, color = my_black)
hlines!(ax1_top, [α/n for n in 1:20], linewidth = 2, linestyle = :dash, color = my_black)

hlines!(ax2, [α/n for n in 1:20], linewidth = 2, linestyle = :dash, color = my_black)
hlines!(ax2_top, [α/n for n in 1:20], linewidth = 2, linestyle = :dash, color = my_black)

hlines!(ax3, [α/n for n in 1:20], linewidth = 2, linestyle = :dash, color = my_black)

# Capture speed
hlines!(ax1, [sqrt(8*Φ0*π^2)], linewidth = 2, color = my_red, linestyle = :dash, label = "Capture speed")

# Bias value
vlines!(ax3, [bias], linewidth = 3, color = my_black, label = "Bias")

# Read in full trajectories
filenames = filter(x -> first(x) !== '.', readdir(joinpath(pwd(), "data/proc_tau800/")))

for ii in eachindex(filenames)
    data = readdlm(joinpath("data/proc_tau800/", filenames[ii]))
    τs = data[:,1]
    speeds = data[:,2]

    if occursin("Phi2", filenames[ii])
        lines!(ax1, τs, speeds, linewidth = 4, color = my_vermillion)
        lines!(ax1_top, τs, speeds, linewidth = 4, color = my_vermillion)

    elseif occursin("Phi-2", filenames[ii])
        lines!(ax2_top, τs, speeds, linewidth = 4, color = my_blue)
        lines!(ax2, τs, speeds, linewidth = 4, color = my_blue)
    end

end


hidexdecorations!(ax1_top)
xlims!(ax1_top, 0.0, 800.0)
ylims!(ax1_top, 86.0, 86.9)
ax1_top.yticks = 86:1:90


xlims!(ax1, 0.0, 800.0)
ylims!(ax1, 0.0, 56.0)
ax1.yticks = 0:10:52

rowsize!(fig.layout, 2, Relative(4.6/5))
rowgap!(fig.layout, 1, Relative(0.0))
ax1_top.bottomspinevisible=false
ax1.topspinevisible=false


hidexdecorations!(ax2_top)
hideydecorations!(ax2_top)


hideydecorations!(ax2)
ax2_top.bottomspinevisible=false
ax2.topspinevisible=false
ax2.xticks = 200:200:800

xlims!(ax2, 0.0, 800.0)
ylims!(ax2, 0.0, 56.0)
xlims!(ax2_top, 0.0, 800.0)
ylims!(ax2_top, 86.0, 86.9)


## Drift velocity
scatter!(ax3_top, Δ_bias, xs, color = my_black, linewidth = 2)
scatter!(ax3, Δ_bias, xs, color = my_black, linewidth = 2)

hidexdecorations!(ax3_top)
hideydecorations!(ax3_top, ticklabels=false)
ax3_top.yticks = 86:1:87

ax3_top.bottomspinevisible=false
ax3.topspinevisible=false
ax3.yticks = 0:10:60

xlims!(ax3, 0.0, 1.0)
ylims!(ax3, 0.0, 56.0)
xlims!(ax3_top, 0.0, 1.0)
ylims!(ax3_top, 86.0, 86.9)

axislegend(ax1, labelsize = 30, position = :rb)
axislegend(ax3, labelsize = 30, position = :rb)

# save("Drift_Velocities.pdf", fig)
fig