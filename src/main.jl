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

## Parameters
η = 1e-12                       # Small number
ϵ = 1e-20
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
    δ::Float64                  # Time step
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
@inline function ω(ωmax, x)
    return sqrt(1 + (ωmax^2 - 1) * sin(x)^2)
end

# Chain recoil function
function Γ(τ, ls, ωmax)
    int_fun(x) = cos.(2 * x * ls) * sin(2 * π * τ * ω(ωmax, x)) / ω(ωmax, x)
    res = quadgk(int_fun, 0, π / 2, atol=1e-8)
    return (res[1] * 2 / π)
end


# Covariance matrix between chain mass displacements
function C_corr(τ, ls, ωmax, ωT)
    int_fun(x) = cos.(2 * x * ls) * cos(2 * π * τ * ω(ωmax, x)) / ω(ωmax, x) * coth(ω(ωmax, x) / 2 / ωT)
    res = quadgk(int_fun, 0, π / 2, atol=1e-8)
    return (res[1] / π)
end

# Covariance matrix between chain mass displacements
function effective_T(ωmax, ωT)
    int_fun(x) = ω(ωmax, x) * coth(ω(ωmax, x) / 2 / ωT)
    res = quadgk(int_fun, 0, π / 2, atol=1e-8)
    return (res[1] / π)
end


# Precompute recoil term and take existing files into account
function mkChainSystem(ωmax, τ_max, lmax, d, Γ_prev)
    δ = (1 / ωmax) / d                          # Time step in units of t_slow
    n_pts = floor(τ_max / δ) |> Int             # Number of time steps given τ_max and δ
    Γ_dim = size(Γ_prev)                        # Existing matrix size (if any)
    Γ_mat = zeros(lmax + 1, n_pts)              # Prepare the Γ matrix

    # No precomputation files exist
    if isempty(Γ_prev)
        println("No existing file, starting precomputation")
        pr = Progress(n_pts)                        # Setup the progress meter
        pts_per_thread = n_pts ÷ Threads.nthreads() # Calculate points per thread

        # Reorder the time steps so that the load is equal for every thread because latter times take longer to evaluate
        thread_pts =
            reshape(1:(Threads.nthreads()*pts_per_thread), Threads.nthreads(), :)' |> vec |> collect

        Threads.@threads for ii in thread_pts
            Γ_mat[:, ii] = Γ(δ * ii, 0:lmax, ωmax)
            next!(pr)
        end

        # Evaluate the remaining time steps if needed
        if length(thread_pts) != n_pts
            Threads.@threads for ii in (length(thread_pts)+1:n_pts)
                Γ_mat[:, ii] = Γ(δ * ii, 0:lmax, ωmax)
                next!(pr)
            end
        end

    # A file exists with enough time steps, but not enough chain masses
    elseif lmax > Γ_dim[1] && n_pts <= Γ_dim[2]
        println("Sufficient number of time steps, calculating for more chain masses")
        Γ_mat[1:Γ_dim[1], :] = Γ_prev[:,1:n_pts]   # Fill matrix with previous data
        pr = Progress(n_pts)                        # Setup the progress meter
        pts_per_thread = n_pts ÷ Threads.nthreads() # Calculate points per thread

        # Reorder the time steps
        thread_pts =
            reshape(1:(Threads.nthreads()*pts_per_thread), Threads.nthreads(), :)' |> vec |> collect

        indices = (Γ_dim[1]+1:lmax)
        Threads.@threads for ii in thread_pts
            Γ_mat[indices, ii] = Γ(δ * ii, indices, ωmax)
            next!(pr)
        end

        # Evaluate the remaining time steps
        if length(thread_pts) != n_pts
            Threads.@threads for ii in (length(thread_pts)+1:n_pts)
                Γ_mat[indices, ii] = Γ(δ * ii, indices, ωmax)
                next!(pr)
            end
        end

    # A file exists with enough chain masses, but not enough time steps
    elseif lmax <= Γ_dim[1] && n_pts > Γ_dim[2]
        println("Sufficient number of chain atoms, calculating for more time steps")
        Γ_mat[:, 1:Γ_dim[2]] = Γ_prev[1:lmax+1,:]   # Fill matrix with previous data
        rem_pts = n_pts - Γ_dim[2]
        pr = Progress(rem_pts)                        # Setup the progress meter
        pts_per_thread = rem_pts ÷ Threads.nthreads() # Calculate points per thread

        # Reorder the time steps
        thread_pts =
            reshape(Γ_dim[2]+1:(Threads.nthreads()*pts_per_thread+Γ_dim[2]), Threads.nthreads(), :)' |>
            vec |>
            collect

        Threads.@threads for ii in thread_pts
            Γ_mat[:, ii] = Γ(δ * ii, 0:lmax, ωmax)
            next!(pr)
        end

        # Evaluate the remaining time steps
        if length(thread_pts) != rem_pts
            Threads.@threads for ii in (length(thread_pts)+1:n_pts)
                Γ_mat[:, ii] = Γ(δ * ii, 0:lmax, ωmax)
                next!(pr)
            end
        end

    # A file exists with insufficient number of chain masses and time steps
    elseif lmax > Γ_dim[1] && n_pts > Γ_dim[2]
        println("Calculating for more chain masses and time steps")
        Γ_mat[1:Γ_dim[1], 1:Γ_dim[2]] = Γ_prev   # Fill matrix with previous data
        pr = Progress(n_pts)                        # Setup the progress meter
        pts_per_thread = n_pts ÷ Threads.nthreads() # Calculate points per thread

        # Reorder the time steps
        thread_pts =
            reshape(1:(Threads.nthreads()*pts_per_thread), Threads.nthreads(), :)' |>
            vec |>
            collect

        Threads.@threads for ii in thread_pts
            indices1 = ii <= Γ_dim[2] ? (Γ_dim[1]+1:lmax+1) : (1:lmax+1)
            indices2 = ii <= Γ_dim[2] ? (Γ_dim[1]:lmax) : (0:lmax)
            Γ_mat[indices1, ii] = Γ(δ * ii, indices2, ωmax)
            next!(pr)
        end

        # Evaluate the remaining time steps
        if length(thread_pts) != n_pts
            Threads.@threads for ii in (length(thread_pts)+1:n_pts)
                indices1 = ii <= Γ_dim[2] ? (Γ_dim[1]+1:lmax+1) : (1:lmax+1)
                indices2 = ii <= Γ_dim[2] ? (Γ_dim[1]:lmax) : (0:lmax)
                Γ_mat[indices1, ii] = Γ(δ * ii, indices2, ωmax)
                next!(pr)
            end
        end
    end

    return ChainSystem(ωmax, δ, Γ_mat)
end

# Mode amplitude
function ζq(ωq, ωT)
    # Subtract a small number from p. The reason is that for low ωT, p ≈ 1,
    # causing issues with the rand() generator
    n = rand(Geometric(1 - exp(-ωq / ωT) - η))
    res = √(n + 1 / 2) * √(2 / ωq)
    return res
end

# Homogeneous displacement of the chain atom g at time step n given a set of ωs
# and the corresponding ζs and ϕs as a sum of normal coordinates.
# function ρH(n, δ, ζs, ϕs, ωs, gs)
#     n_ζ = length(ζs)
#     f = transpose(exp.(-1im * (2 * π * δ * n * ωs + ϕs)) / √(n_ζ) .* ζs)
#     r = [f * exp.(1im * 2 * π / n_ζ .* (1:n_ζ) * g) for g in gs]
#     return r
# end

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
    periodic::Bool = false,
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

        if periodic
            for particle in eachindex(σ_dot)
                if (σs[particle, nxt] < box[1])
                    σs[particle, nxt] += (box[2] - box[1])
                elseif (σs[particle, nxt] > box[2])
                    σs[particle, nxt] -= (box[2] - box[1])
                end
            end

        else 
            for particle in eachindex(σ_dot)
                if (σs[particle, nxt] < box[1] || σs[particle, nxt] > box[2])
                    σ_dot[particle] = -σ_dot[particle]
                end
            end
        end

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

# Analytic dissipation for Gaussian potential
function Δ_analytic(v, Φ, λ, Ω)
    z = (2 * π * λ / v)^2
    return (
        4 * π^3 * Φ^2 / v^2 *
        z *
        exp(-z * (Ω^2 + 1) / 2) *
        (
            besseli(0, z * (Ω^2 - 1) / 2) +
            (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
        )
    )
end

function motion_solver_test(
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
    periodic::Bool = false,
    check_stuck::Bool = false,
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
    final_ind = n_pts
    # for ii = 3:n_pts
    for ii = 3:n_pts
        nxt = ii        # Next time step index
        curr = ii - 1   # Current time step index

        # If the particle is stuck, break the loop
        if check_stuck && σ_dot[1] <= 0
            final_ind = curr
            break
        end

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

        if periodic
            for particle in eachindex(σ_dot)
                if (σs[particle, nxt] < box[1])
                    σs[particle, nxt] += (box[2] - box[1])
                elseif (σs[particle, nxt] > box[2])
                    σs[particle, nxt] -= (box[2] - box[1])
                end
            end

        else 
            for particle in eachindex(σ_dot)
                if (σs[particle, nxt] < box[1] || σs[particle, nxt] > box[2])
                    σ_dot[particle] = -σ_dot[particle]
                end
            end
        end

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
    return SystemSolution(ωmax, μ, τs[1:final_ind], τ0, α, Φ0, λ, σs[:, 1:final_ind], ρs[:, 1:final_ind], bias, tTraj.ωT)
end

# Calculate energy losses along a full trajectory
function Δ_traj(data)
    δ = data.τs[2] - data.τs[1]
    σs = data.σs |> vec
    ρs = data.ρs

    # Find the chain indices where the mobile particle starts and ends
    chain_idx = searchsortedlast(ρs[:,1], σs[1])
    max_lattice_pos = min(data.ρs[end, 1], maximum(σs)) / data.α |> floor |> Int

    # Find the indices halfway between chain masses as chain moves
    idx = [argmin(abs.(σs .- (ρs[n, :] .+ (0.5 * data.α)))) for n in chain_idx:max_lattice_pos]
    idx = filter(x -> x <= length(σs) - 1, idx)

    # Get speeds at these points and corresponding times
    σ_dots = (σs[idx.+1] - σs[idx]) ./ δ
    τs = data.τs[idx]

    # Get the kinetic energy and filter out energies less that potential height
    KE = 0.5 * data.μ * (σ_dots ./ 2 ./ π) .^ 2
    idx = findall(x -> x > data.Φ, KE)
    KE = KE[idx]
    σ_dots = σ_dots[idx]
    Δs = KE[1:end-1] - KE[2:end]
    return (σ_dots[1:end-1], Δs)
end

## Unfold particle trajectories given a confining box
function traj_unfold(data, box; periodic = false, minus_init = false)
    σs_final = copy(data.σs)
    num_P = size(data.σs, 1)

    for particle in 1:num_P
        if periodic 
            for st in 2:lastindex(σs_final[particle, :])
                if (σs_final[particle, st] - σs_final[particle, st-1]) > 0.9 * (box[2] - box[1])
                    σs_final[particle, st:end] = σs_final[particle, st:end] .- (box[2] - box[1])
                elseif (σs_final[particle, st] - σs_final[particle, st-1]) < -0.9 * (box[2] - box[1])
                    σs_final[particle, st:end] = σs_final[particle, st:end] .+ (box[2] - box[1])
                end
            end
        else 
            τ_ids = findall(x -> x < box[1] || x > box[2], σs_final[particle,:])

            for τ_id in τ_ids
                σs_final[particle, (τ_id+1):end] = 2 * σs_final[particle, τ_id] .- σs_final[particle, (τ_id+1):end]
            end
        end
        if minus_init
            σs_final[particle,:] = σs_final[particle,:] .- σs_final[particle,1]
        end
    end

    return σs_final
end


# Calculate energy losses along a full trajectory
function Δ_traj_unfold(data, box; periodic = false)
    δ = data.τs[2] - data.τs[1]
    σs = traj_unfold(data, box, periodic = periodic) |> vec

    # Find the chain indices where the mobile particle starts and ends
    chain_idx = searchsortedlast(data.ρs[:,1], σs[1])
    max_lattice_pos = maximum(σs) / data.α |> floor |> Int

    # Find the indices halfway between chain mass rest positions
    idx = [argmin(abs.(σs .- ((n + 0.5) * data.α))) for n in chain_idx:max_lattice_pos]
    idx = filter(x -> x <= length(σs) - 1, idx)

    # Get speeds at these points and corresponding times
    σ_dots = (σs[idx.+1] - σs[idx]) ./ δ

    # Get the kinetic energy a
    KE = 0.5 * data.μ * (σ_dots ./ 2 ./ π) .^ 2
    Δs = KE[1:end-2] - KE[2:end-1]

    return (σ_dots[1:end-2], Δs)
end

# Energy loss DeltaBar in the steady-state limit
function Δ_transport(v, Φ, λ, Ω, α)
    # Find zeros of function
    sols = find_zeros(x -> ω(Ω, π*x) + x*(v/α), -Ω*α/v, 0, no_pts = 55)

    ωprime(x) = π*sin(π*x)*cos(π*x)*(Ω^2 -1)/ω(Ω, π*x)

    vals = [(4*π^2*λ*Φ*ω(Ω, π*sol)/v^2)^2 * π * exp(-2*π^2*λ^2*(ω(Ω, π*sol)^2)^2/v^2) * (1 / abs(α*ωprime(sol)/v + 1)) for sol in sols]

    return isempty(vals) ? 0 : sum(vals)
end

# Gaussian interaction profile
function U_profile(r, Φ, λ)
    return Φ*exp(-(r^2)/(2*λ^2))
end

# Broadened energy loss DeltaBar in the steady-state limit
function Δ_transport_broadened(v, Φ, λ, Ω, α, μ)
    # Get range of particle velocities as it passes a chain mass
    num_speeds = 10000
    v_ext = sqrt(v^2 - (8*π^2*Φ/μ))
    v_range = range(min(v, v_ext), max(v, v_ext), length = num_speeds)

    # Get probability density of speeds as particle passes a chain mass
    xs = range(-α/v/2, α/v/2, length = num_speeds)
    kde_data = kde(sqrt.(v^2 .- 8*π^2/μ*U_profile.(xs, Φ, λ/v)), npoints = num_speeds, boundary = (min(v, v_ext), max(v, v_ext)))

    # Obtain weighted sum of Deltabar solutions
    vals = 0
    for ii in eachindex(v_range)
        sols = find_zeros(x -> ω(Ω, π*x) + x*(v_range[ii]/α), -Ω*α/v_range[ii], 0, no_pts = 60)
        ωprime(x) = π*sin(π*x)*cos(π*x)*(Ω^2 -1)/ω(Ω, π*x)

        for sol in sols
            vals += kde_data.density[ii] * (4*π^2*λ*Φ*ω(Ω, π*sol)/v_range[ii]^2)^2 * π * exp(-2*π^2*λ^2*(ω(Ω, π*sol)^2)^2/v_range[ii]^2) * (1 / abs(α*ωprime(sol)/v_range[ii] + 1))
        end
    end

    return (vals/sum(kde_data.density))
end


# Mean Delta with thermal fluctuations 
function Δ_thermal_analytic(v, Φ, λ, Ω, ωT)
    factor = π^2*λ^2*Φ^2*√(π) / v / 2

    C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)

    int_func(θ, τ) = exp(2im * π * τ * ω(Ω, θ)) * 
                    exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) * 
                    (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2*τ^2) / 
                    (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5/2)
                    
    τlim = min(30 / v, 3)
    res = hcubature(x -> int_func(x[1], x[2]), [0, -τlim], [π/2, τlim], rtol = 1e-2, initdiv = 1)

    return factor * res[1] * 2 / π |> real
end


# Delta variance with thermal fluctuations
function Δ_thermal_variance(v, Φ, λ, Ω, ωT)
    factor = v * Φ^2 * λ^2 * √(π) / 4

    C_NN_func(τ) = C_corr(τ, 0, Ω, ωT)

    int_func(τ) = exp(-v^2 * τ^2 / (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) / 4) * 
                (2(λ^2 + Complex(C_NN_func(0) - C_NN_func(τ))) - v^2*τ^2) / 
                (λ^2 + Complex(C_NN_func(0) - C_NN_func(τ)))^(5/2)


    res = quadgk(int_func, -40, 40, atol = 1e-9)
    
    return factor*res[1] |> real
end
