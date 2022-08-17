using Distributed
proc_number = 3
if nprocs() < proc_number
    addprocs(proc_number - nprocs())
end

@everywhere include("cluster_parameters.jl")

println("Starting Calculations")

pmap(multi_rep, speeds)
pmap(multi_att, speeds)
