data = load_object("box_test_7.jld2")

ωT = data.ωT
Φ0 = data.Φ
λ = data.λ
σs = data.σs
ρs = data.ρs
τs = data.τs
δ = τs[2] - τs[1]

drop = 50000
σs = σs[:, drop:end]
ρs = ρs[:, drop:end]

prt = size(σs)[1]

ens = Vector{Float64}([])
for p = 1:25
    println(p)
    kin_en = ((σs[p, 2:end] - σs[p, 1:(end-1)]) ./ δ) .^ 2 ./ 2 ./ ωT ./ (2 * pi)^2
    pot_en = [
        sum((Φ0 .* exp.(-(ρs[:, t] .- σs[p, t]) .^ 2 ./ (2 * λ^2))) ./ ωT) for
        t in eachindex(σs[p, :])
    ]
    pot_en = pot_en[2:end]
    idx = findall(x -> x < 1e-3, pot_en)
    ens = vcat(ens, kin_en[idx])
    # tot_en = (pot_en[2:end] + kin_en) |> vec
    # ens = vcat(ens, tot_en)
end

hist_fit = fit(Histogram, ens |> vec, 0.01:0.01:4.5)
hist_fit = normalize(hist_fit, mode = :pdf)

scatter(
    (hist_fit.edges[1])[1:end-1],
    log.(hist_fit.weights .* sqrt.((hist_fit.edges[1])[1:end-1])),
    # label = lab,
    markersize = 12,
)

length(ens)

t = 120
(Φ0 .* exp.(-(ρs[:, t] .- σs[t]) .^ 2 ./ (2 * λ^2))) ./ ωT |> sum
findmin(tot_en)
mean(tot_en)
σs[70553]
σs[70552]
σs[70551]
σs[70550]

ρs[30, 70550]
# (Φ0 .* exp.(-(5) .^ 2 ./ (2 * λ^2)))
τs[70551]
