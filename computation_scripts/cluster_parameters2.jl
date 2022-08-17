using Distributed
proc_number = 3
if nprocs() < proc_number
    addprocs(proc_number - nprocs())
end

@everywhere include("non_thermal2.jl")

pmap(full_traj, param_vals)
