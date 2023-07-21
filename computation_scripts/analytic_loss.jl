include("../src/main.jl")

Φ0 = 1
λ = 1
ωmax = 10

nPts = 10000
vMin = 0.1
vMax = 1000

vs = range(vMin, vMax, length = nPts) |> collect
vs_idx = reduce(
    vcat,
    [collect(eachindex(vs))[n:Threads.nthreads():end] for n = 1:Threads.nthreads()],
)

ωTs = [0, 1, 5, 25, 100, 250]

for ωT in ωTs
    println("ωT = $(ωT)")
    if (!isfile("data/Thermal/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2"))
        result = zeros(nPts)
        pr = Progress(nPts)
        Threads.@threads for ii in vs_idx
            result[ii] = Δ_thermal_analytic(vs[ii], Φ0, λ, ωmax, ωT)
            next!(pr)
        end
        save_object(
            "data/Thermal/AnalyticLoss_Φ0$(Φ0)_λ$(λ)_ωmax$(ωmax)_ωT$(ωT).jld2",
            (vs, result),
        )
    end
end
