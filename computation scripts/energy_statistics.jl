include("../src/main.jl")

data = [
    # "data/Box/Box_τ010_λ1_Φ07_ωT0.0_τ1000.jld2",
    "data/Box/Box_τ010_λ1_Φ07_ωT1.0_τ1000.jld2",
    "data/Box/Box_τ010_λ1_Φ07_ωT2.0_τ1000.jld2",
    "data/Box/Box_τ010_λ1_Φ07_ωT5.0_τ1000.jld2",
    "data/Box/Box_τ010_λ1_Φ07_ωT10.0_τ1000.jld2",
    "data/Box/Box_τ010_λ1_Φ07_ωT25.0_τ1000.jld2",
    "data/Box/Box_τ010_λ1_Φ07_ωT100.0_τ1000.jld2",
    "data/Box/Box_τ010_λ1_Φ07_ωT250.0_τ1000.jld2",
]

colors = [my_red, my_vermillion, my_orange, my_yellow, my_green, my_sky, my_blue]

fig = Figure(
    resolution = (1200, 800),
    fonts = (; math = "CMU Serif"),
    fontsize = 40,
    figure_padding = 30,
)

ax1 = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}",
    ylabel = L"\dot{\sigma}\langle\Delta\rangle / \Phi^2_0",
    # xscale = log10,
    # yscale = log10,
    xticklabelfont = :math,
    yticklabelfont = :math,
    xgridvisible = false,
    ygridvisible = false,
    xlabelpadding = -15,
    # limits = (1e-1, 300, 1e-2, 1),
)
for ii in eachindex(data)
    d = load_object(data[ii])
    ωT = d.ωT
    Φ0 = d.Φ
    λ = d.λ
    σs = d.σs
    ρs = d.ρs
    τs = d.τs
    δ = τs[2] - τs[1]

    if (!isfile("data/Kinetic/Kinetic_τ010_λ$(λ)_Φ0$(Φ0)_ωT$(d.ωT)_τ1000.jld2"))
        prt = size(σs)[1]
        ens = Vector{Float64}([])
        for p = 1:prt
            println(p)
            kin_en = ((σs[p, 2:end] - σs[p, 1:(end-1)]) ./ δ) .^ 2 ./ 2 ./ ωT ./ (2 * pi)^2
            pot_en = [
                sum((Φ0 .* exp.(-(ρs[:, t] .- σs[p, t]) .^ 2 ./ (2 * λ^2))) ./ ωT) for
                t in eachindex(σs[p, :])
            ]
            pot_en = pot_en[2:end]
            idx = findall(x -> x < 1e-3, pot_en)
            ens = vcat(ens, kin_en[idx])
        end
        save_object("data/Kinetic/Kinetic_τ010_λ$(λ)_Φ0$(Φ0)_ωT$(d.ωT)_τ1000.jld2", ens)
    end

    ens = load_object("data/Kinetic/Kinetic_τ010_λ$(λ)_Φ0$(Φ0)_ωT$(d.ωT)_τ1000.jld2")
    hist_fit = fit(Histogram, ens |> vec, 0.0:0.05:10)
    hist_fit = normalize(hist_fit, mode = :pdf)
    scatter!(
        ax1,
        (hist_fit.edges[1])[1:end-1],
        log.(hist_fit.weights .* sqrt.((hist_fit.edges[1])[1:end-1])),
        color = colors[ii],
        label = L"\omega_T = %$(ωT)",
        markersize = 12,
    )
end

lines!(ax1, 0:5, -(0:5))
axislegend(position = :rt, orientation = :vertical, nbanks = 1)


fig
