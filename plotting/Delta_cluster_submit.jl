using Distributed
proc_number = 6
if nprocs() < proc_number
    addprocs(5)
end
@everywhere include("single_pass_Delta_cluster.jl")

for ωT in ωTs
    speeds = range(10.0, 80.0, step = 5.0)
    map(x -> Δ_distribution(x, ωT), speeds)
end
