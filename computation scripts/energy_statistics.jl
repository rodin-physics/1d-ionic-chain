include("../src/main.jl")

data = [
    # "data/Box/Box_τ010_λ1_Φ07_ωT0.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT1.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT2.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT5.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT10.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT15.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT25.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT50.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT100.0_τ1000.jld2",
    "data/Box/BoxMultiple_τ010_λ1_Φ07_ωT250.0_τ1000.jld2",
]

colors = [
    my_red,
    my_vermillion,
    my_orange,
    my_yellow,
    my_green,
    my_sky,
    my_blue,
    my_red,
    my_vermillion,
    my_orange,
    my_yellow,
    my_green,
    my_sky,
    my_blue,
]

fig = Figure(
    resolution = (1200, 800),
    fonts = (; math = "CMU Serif"),
    fontsize = 40,
    figure_padding = 30,
)

ax1 = Axis(
    fig[1, 1],
    xlabel = L"\mathcal{E} / \omega_T",
    ylabel = L"\ln[P(\mathcal{E})]",
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
    drop = 1
    σs = d.σs[:, drop:end]
    ρs = d.ρs[:, drop:end]
    τs = d.τs[drop:end]
    δ = τs[2] - τs[1]

    if (!isfile("data/Kinetic/Total_τ010_λ$(λ)_Φ0$(Φ0)_ωT$(d.ωT)_τ1000.jld2"))
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
            # idx = findall(x -> x < 1e-3, pot_en)
            ens = vcat(ens, kin_en + pot_en)
        end
        save_object("data/Kinetic/Total_τ010_λ$(λ)_Φ0$(Φ0)_ωT$(d.ωT)_τ1000.jld2", ens)
    end

    ens = load_object("data/Kinetic/Total_τ010_λ$(λ)_Φ0$(Φ0)_ωT$(d.ωT)_τ1000.jld2")
    kin_ens = Vector{Float64}([])
    for p = 1:prt
        kin_en = ((σs[p, 2:end] - σs[p, 1:(end-1)]) ./ δ) .^ 2 ./ 2 ./ ωT ./ (2 * pi)^2
        # idx = findall(x -> x < 1e-3, pot_en)
        kin_ens = vcat(kin_ens, kin_en)
    end
    # println(length(ens|>vec))
    # println(length(kin_ens|>vec))
    hist_fit = fit(Histogram, (ens |> vec), weights(sqrt.(kin_ens)), 0.0:0.025:5)
    # hist_fit = fit(Histogram, (ens |> vec), weights(sqrt.(ens)), 0.0:0.025:5)
    hist_fit = normalize(hist_fit, mode = :pdf)
    scatter!(
        ax1,
        (hist_fit.edges[1])[1:end-1],
        log.(hist_fit.weights),
        # log.(hist_fit.weights .* sqrt.((hist_fit.edges[1])[1:end-1])),
        color = colors[ii],
        label = L"\omega_T = %$(ωT)",
        markersize = 12,
    )
end

lines!(ax1, 0:5, -(0:5), linewidth = 4, color = my_black)
axislegend(position = :lb, orientation = :vertical, nbanks = 1)


fig
save("Thermalization.pdf", fig)
