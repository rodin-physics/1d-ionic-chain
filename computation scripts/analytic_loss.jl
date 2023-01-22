include("../src/main.jl")

function Δ_thermal_analytic(v, Φ, λ, ωmax, ωT)
    factor = π^2 * λ^2 * Φ^2 * √(π) / v / 2

    C_NN_func(τ) = C_corr(τ, 0, ωmax, ωT)
    τlim = min(30 / v, 3)
    int_func(θ, τ) =
        exp(2im * π * τ * ω(ωmax, θ)) *
        exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) *
        (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2 * τ^2) /
        (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5 / 2)

    res = hcubature(
        x -> int_func(x[1], x[2]),
        [0, -τlim],
        [π / 2, τlim],
        rtol = 1e-2,
        initdiv = 1,
    )

    return factor * res[1] * 2 / π |> real
end

Φ0 = 1
λ = 1
ωmax = 10

nPts = 6000
vMin = 0.1
vMax = 300

vs = range(vMin, vMax, length = nPts) |> collect
vs_idx = reduce(
    vcat,
    [collect(eachindex(vs))[n:Threads.nthreads():end] for n = 1:Threads.nthreads()],
)

ωTs = [0, 1, 5, 25, 100, 250]

for ωT in ωTs
    println("ωT = $(ωT)")
    if (!isfile("data/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2"))
        result = zeros(nPts)
        pr = Progress(nPts)
        Threads.@threads for ii in vs_idx
            result[ii] = Δ_thermal_analytic(vs[ii], Φ0, λ, ωmax, ωT)
            next!(pr)
        end
        save_object(
            "data/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2",
            (vs, result),
        )
    end
end
