using CairoMakie
using Colors
using Distributions
using JLD2
using LaTeXStrings
using LinearAlgebra
using ProgressMeter
using QuadGK
using SpecialFunctions
using Statistics
using StatsBase
using ToeplitzMatrices
using Roots
using KernelDensity
using DelimitedFiles
using HCubature
using Random

## Parameters
ϵ = 1e-12       # Minimum energy cutoff
## Friendly colors
my_red = colorant"rgba(204, 121, 167, 1.0)"
my_vermillion = colorant"rgba(213, 94, 0, 1.0)"
my_orange = colorant"rgba(230, 159, 0, 1.0)"
my_yellow = colorant"rgba(240, 228, 66, 1.0)"
my_green = colorant"rgba(0, 158, 115, 1.0)"
my_sky = colorant"rgba(86, 180, 233, 1.0)"
my_blue = colorant"rgba(0, 114, 178, 1.0)"
my_black = colorant"rgba(0, 0, 0, 1.0)"

## Types
struct ChainSystem
    ωmax::Float64               # Largest mode frequency
    δ::Float64                  # Timestep
    Γ::Matrix{Float64}          # Response array
end

struct ThermalTrajectory
    ωmax::Float64               # Largest mode frequency
    δ::Float64                  # Time step
    ρHs::Matrix{Float64}        # Thermal trajectory
    ωT::Union{Nothing,Float64}  # Temperature
end

struct SystemSolution
    ωmax::Float64               # Largest mode frequency
    μ::Float64                  # Mass of the mobile atoms
    τs::Vector{Float64}         # Time steps
    τ0::Float64                 # Memory
    α::Float64                  # Spacing between chain atoms
    Φ::Float64                  # Magnitude of the Gaussian potential
    λ::Float64                  # Standard deviation of the potential
    σs::Matrix{Float64}         # Positions of mobile atoms
    ρs::Matrix{Float64}         # Positions of the chain atoms
    bias::Float64               # Applied bias
    ωT::Union{Nothing,Float64}  # Temperature
end

## Functions
# Frequency as a function of momentum
@inline function ω(ωmax, θ)
    return sqrt(1 + (ωmax^2 - 1) * sin(θ)^2)
end

# Chain recoil function
function Γ(τs, ls, ωmax)
    int_fun(θ) = cos.(2 * θ * ls) * sin.(2 * π * τs * ω(ωmax, θ))' / ω(ωmax, θ)
    res = quadgk(int_fun, 0, π / 2, atol = 1e-8)
    return (res[1] * 2 / π)
end

# Covariance matrix between chain mass displacements
function C_corr(τs, ls, ωmax, ωT)
    int_fun(θ) =
        cos.(2 * θ * ls) * cos.(2 * π * τs * ω(ωmax, θ))' / ω(ωmax, θ) *
        coth(ω(ωmax, θ) / 2 / ωT)
    res = quadgk(int_fun, 0, π / 2, atol = 1e-8)
    return (res[1] / π)
end

# Precompute recoil term
function mkChainSystem(ωmax, τmax, lmax, δ; batch_size = 100)
    τs = range(0, τmax, step = δ)
    ls = 0:lmax
    Γ_mat = zeros(length(ls), length(τs))       # Preallocate Γ matrix

    # Reorder the time steps so that the load is equal for every thread
    # because latter times take longer to evaluate
    τs_idx = reduce(
        vcat,
        [collect(eachindex(τs))[n:Threads.nthreads():end] for n = 1:Threads.nthreads()],
    )

    τs_idx_batches = Iterators.partition(τs_idx, batch_size) |> collect
    pr = Progress(length(ls) * length(τs_idx_batches))  # Setup the progress meter

    for ll in eachindex(ls)
        Threads.@threads for ii in τs_idx_batches
            Γ_mat[ll, ii] = Γ(τs[ii], ls[ll], ωmax)
            next!(pr)
        end
    end

    return ChainSystem(ωmax, δ, Γ_mat)
end

# Mode amplitude
function ζq(ωq, ωT)
    # Subtract a small number from p. The reason is that for low ωT, p ≈ 1,
    # causing issues with the rand() generator
    η = 1e-12
    n = rand(Geometric(1 - exp(-ωq / ωT) - η))
    res = √(n + 1 / 2) * √(2 / ωq)
    return res
end

# Homogeneous displacement of the chain atoms at time τ given a set of ωs
# and the corresponding ζs and ϕs
function ρH(τ, ζs, ϕs, ωs, ε)
    n_ζ = length(ζs)
    # Get the displacement of each mode ζ
    f = transpose(exp.(-1im * (2 * π * τ * ωs + ϕs)) / √(n_ζ) .* ζs)
    # Multiply vector of ζ's by the polarization matrix ε
    res = f * ε
    return res
end

function motion_solver(
    system::ChainSystem,
    Φ0::T where {T<:Real},
    λ::T where {T<:Real},
    α::T where {T<:Real},
    σ0::Vector{T} where {T<:Real},
    σ_dot0::Vector{T} where {T<:Real},
    μ::T where {T<:Real},
    tTraj::ThermalTrajectory,
    τ0::T where {T<:Real},
    τ::T where {T<:Real};
    threads::Bool = false,
    bias::T where {T<:Real} = 0,
    box::Tuple{T,T} where {T<:Real} = (-Inf, Inf),
)
    ωmax = system.ωmax              # Maximum chain frequency
    Γ_mat = system.Γ                # Memory term
    δ = system.δ                    # Time step
    n_pts = floor(τ / δ) |> Int     # Number of time steps
    τs = δ .* (1:n_pts) |> collect  # Times
    nChain = size(tTraj.ρHs)[1]     # Number of chain particles for which the homogeneous motion is available
    F_bias = bias / α               # Force due to the bias
    τ0_pts = max(floor(τ0 / δ), 1)  # Memory time points
    # Even if τ0 == 0, we have to keep a single time point to make sure arrays work.
    # For zero memory, the recoil contribution is dropped

    # Check that the thermal trajectory is for the correct system
    if (ωmax != tTraj.ωmax || δ != tTraj.δ)
        error("Thermal trajectory describes a different system. Check your input.")
    elseif size(tTraj.ρHs)[2] < n_pts
        error("Thermal trajectory provided does not span the necessary time range.")
        # If the precomputed memory is shorter than the simulation time AND shorter
        # than the desired memory, terminate the calculation.
    elseif (size(Γ_mat)[2] < n_pts && size(Γ_mat)[2] < τ0_pts)
        error("Chosen memory and the simulation time exceed the precomputed recoil.")
        # If the desired number of chain particles is greater than what is contained in
        # the precomputed Γ, terminate the calculation. Otherwise, retain the appropriate
        # number of terms
    elseif (size(Γ_mat)[1] < nChain)
        error("The recoil term does not contain the desired number of chain masses.")
    else
        ρs = (tTraj.ρHs)[:, 1:n_pts] .+ α .* repeat(1:size(tTraj.ρHs)[1], 1, n_pts)
        σs = zeros(length(σ0), n_pts)
        τ0_pts = min(τ0_pts, n_pts) |> Int
        Γ_mat = (2 * π * δ) .* Γ_mat[1:nChain, 1:τ0_pts]
        Γ_mat = vcat(reverse(Γ_mat, dims = 1)[1:end-1, :], Γ_mat)

    end
    # Interaction terms
    @inline function dU_dρ(r)
        return (-Φ0 * exp(-r^2 / (2 * λ^2)) * r / λ^2)
    end
    ## Initial values
    σs[:, 1] = σ0
    σs[:, 2] = σ0 + δ .* σ_dot0
    σ_dot = σ_dot0
    U_pr = zeros(nChain, length(σ0))
    # for ii = 3:n_pts
    @showprogress for ii = 3:n_pts
        nxt = ii        # Next time step index
        curr = ii - 1   # Current time step index
        # Calculate the forces on all the masses
        U_pr = [dU_dρ(ρ - σ) for ρ in ρs[:, curr], σ in σs[:, curr]]
        U_pr_chain = sum(U_pr, dims = 2) |> vec
        U_pr_mob = -sum(U_pr, dims = 1) |> vec
        # Find the indices of the chain masses where the force is larger than ϵ
        idx = findall(x -> abs(x) > ϵ, U_pr_chain)

        σs[:, nxt] =
            σs[:, curr] +
            δ .* σ_dot +
            -(2 * π * δ)^2 / μ .* (U_pr_mob - F_bias .* ones(length(σ0)))

        σ_dot = (σs[:, nxt] - σs[:, curr]) / δ
        for particle in eachindex(σ_dot)
            if (σs[particle, nxt] < box[1])
                σs[particle, nxt] += (box[2] - box[1])
            elseif (σs[particle, nxt] > box[2])
                σs[particle, nxt] -= (box[2] - box[1])
            end
        end
        # for particle in eachindex(σ_dot)
        #     if (σs[particle, nxt] < box[1] || σs[particle, nxt] > box[2])
        #         σ_dot[particle] = -σ_dot[particle]
        #     end
        # end

        steps_left = n_pts - curr               # Number of time steps remaining
        step_memory = min(τ0_pts, steps_left)   # Number of steps for which ρs are affected by the impulse
        if threads == true
            steps_per_thread = step_memory ÷ Threads.nthreads() + 1
            step_alloc = [
                (x*steps_per_thread+1):min(step_memory, (x + 1) * steps_per_thread) for
                x = 0:Threads.nthreads()-1
            ]
            Threads.@threads for t = 1:Threads.nthreads()
                for n in idx
                    view(ρs, :, curr .+ step_alloc[t]) .-=
                        view(Γ_mat, nChain-n+1:2*nChain-n, step_alloc[t]) .*
                        U_pr_chain[n] .* (τ0 != 0)
                end
            end
        else
            for n in idx
                view(ρs, :, curr .+ (1:step_memory)) .-=
                    view(Γ_mat, nChain-n+1:2*nChain-n, 1:step_memory) .* U_pr_chain[n] .*
                    (τ0 != 0)
            end
        end
    end
    return SystemSolution(ωmax, μ, τs, τ0, α, Φ0, λ, σs, ρs, bias, tTraj.ωT)
end
