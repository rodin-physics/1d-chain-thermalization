using Distributed
proc_number = 2
if nprocs() < proc_number
    addprocs(proc_number - nprocs())
end

@everywhere include("non_thermal.jl")

pmap(full_traj, param_vals)
