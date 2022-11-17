include("main.jl")

ls = collect(0:1)
ωmax = 10
τs = collect(1:(1/600):300)
# Precompute recoil term and take existing files into account
function mkChainSystem(τs, ls, ωmax; batch_size = 100)
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

    return ChainSystem(ωmax, τs, ls, Γ_mat)
end
# mkChainSystem(τs, ls, ωmax; batch_size = 50)
# show(mkChainSystem(τs, ls, ωmax))
