include("main.jl")

system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")
system.τs

nMasses = 10000
ωT = 1
ωmax = 10
θs = pi .* (1:nMasses) / nMasses
ωs = ω.(ωmax, θs)
ϕs = 2 * π * rand(nMasses)
ζs = ζq.(ωs, ωT)

@time ρH(1, ζs, ϕs, ωs, 0:300) |> real


1.25 * 40 / 10 * 600

rr = load_object("data/NumericalLoss_Φ05_λ1_ωmax10_v80_μ1_ωT1.jld2")
mean(rr[2])
Δ_thermal_analytic(30, 5, 1, 10, 1)
hist(rr[2], bins = 50)
(1 / 2 / (2 * π)^2 * (15^2))
(1 / 2 / (2 * π)^2 * (10^2)) * 0.5
(1 / 2 / (2 * π)^2 * (6^2))
Δ_analytic(80, 5, 1, 10)

zz = load_object("data/NumericalLossFast_Φ01_λ1_ωmax10_μ1_ωT5.0.jld2")
hist(zz[2][end, :], bins = 40)
mean(zz[2][4, :])
Δ_thermal_analytic(24, 1, 1, 10, 1)


function Δ_thermal_analytic(v, Φ, λ, Ω, ωT)
    factor = π^2 * λ^2 * Φ^2 * √(π) / v / 2

    C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)
    # C_NN_func(τ) = 0

    int_func(θ, τ) =
        exp(2im * π * τ * ω(Ω, θ)) *
        exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) *
        (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2 * τ^2) /
        (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)

    res = hcubature(x -> int_func(x[1], x[2]), [0, -1], [π / 2, 1], atol = 1e-4)

    return factor * res[1] * 2 / π |> real
end

rr = load_object("data/NumericalLoss_Φ05_λ1_ωmax10_v80_μ1_ωT1.jld2")
mean(rr[2])
Δ_thermal_analytic(80, 5, 1, 10, 10)
Δ_analytic(80, 5, 1, 10)
(mean(rr[2]) - Δ_thermal_analytic(80, 5, 1, 10, 1)) /
sqrt(Δ2_thermal_analytic(80, 5, 1, 10, 1) - Δ_thermal_analytic(80, 5, 1, 10, 1)^2) *
sqrt(length(rr[2]))


Δ_thermal_analytic(10, 5, 1, 10, 5000)



ts = range(0, 10, length = 1000)
temp = 50
r = -C_corr.(ts, 0, 10, temp) .+ C_corr(0, 0, 10, temp)
lines(ts, r)
C_corr(20, 0, 10, 10000)




2 * pi * 10
@time Δ2_thermal_analytic(2, 5, 1, 10, 100)

vs = range(3, 100, length = 300)
r = Δ2_thermal_analytic.(vs, 1, 1, 10, 50)

lines(log10.(vs), log10.(real(r)))