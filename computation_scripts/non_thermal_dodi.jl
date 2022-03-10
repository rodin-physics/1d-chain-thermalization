using Distributed

# proc_num = 12
# addprocs(proc_num - nprocs())

@everywhere include("../src/main.jl")

system = load_object("../precomputed/systems/System_K1_k20_m1_d60_l200.jld2")

d = 60

τ = 100                             # Simulation time
Ωmin = √(system.K / system.m)       # Minimum phonon frequency
tmin = 2 * π / Ωmin                 # Period of the slowest mode
δ = system.δ                        # Time step
n_pts = floor(τ * tmin / δ) |> Int  # Number of time steps
nChain = 50                         # Number of chain atoms tracked
a = 1                               # Distance between chain atoms

rHs = [a .* collect(1:nChain) for n = 1:n_pts]
tTraj = ThermalTrajectory(system.k, system.K, system.m, a, system.δ, rHs, nothing, ħ)

## Width and Depth Dependence
s = [1 / 2, 1 / 4, 1 / 8, 1 / 16]
F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]
v0 = [[1 / 2]]
M = [1 / 2]
mems = [Inf]

parameters = [(w, x, y, z, mem) for w in M, x in s, y in F, z in v0, mem in mems] |> vec

x0 = [5]

m = system.m            # Mass of the chain atoms
k = system.k            # Spring force constant
K = system.K            # Confining potential force constant
δ = system.δ            # Array of time steps
G_list = system.G       # Memory term
G_list = [x[1:nChain] for x in G_list]

system = ChainSystem(k, K, m, δ, G_list)

# Base.summarysize(system)
@showprogress pmap(parameters) do param
    M = param[1]
    s = param[2]
    F = param[3]
    v0 = param[4]
    mem = param[5]
    if (
        !isfile(
            "../data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
        save_object(
            "../data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

# mass dependence 

s = [1 / 16]
F = [1.0]    #put a decimal next time
v0 = [[1 / 2]]
M = [0.25, 0.5, 1, 1.25, 1.5, 2.0]
mems = [Inf]

parameters = [(w, x, y, z, mem) for w in M, x in s, y in F, z in v0, mem in mems] |> vec
parameters[1]
x0 = [5]

m = system.m            # Mass of the chain atoms
k = system.k            # Spring force constant
K = system.K            # Confining potential force constant
δ = system.δ            # Array of time steps
G_list = system.G       # Memory term
G_list = [x[1:nChain] for x in G_list]

system = ChainSystem(k, K, m, δ, G_list)


@showprogress pmap(parameters) do param
    M = param[1]
    s = param[2]
    F = param[3]
    v0 = param[4]
    mem = param[5]
    if (
        !isfile(
            "../data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
        save_object(
            "../data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end


#velocity dependence

s = [1 / 16]
F = [ 1.0 ]   #put a decimal next time
#v0 = [ [2/2] ]
v0 = [[1 / 4], [1 /2] , [1.0], [3 / 2], [2.0]]
M = [1.0]
mems = [Inf]
parameters = [(w, x, y, z, mem) for w in M, x in s, y in F, z in v0, mem in mems] |> vec
parameters[1]
x0 = [5]

m = system.m            # Mass of the chain atoms
k = system.k            # Spring force constant
K = system.K            # Confining potential force constant
δ = system.δ            # Array of time steps
G_list = system.G       # Memory term
G_list = [x[1:nChain] for x in G_list]

system = ChainSystem(k, K, m, δ, G_list)

@showprogress pmap(parameters) do param
    M = param[1]
    s = param[2]
    F = param[3]
    v0 = param[4]
    mem = param[5]
    if (
        !isfile(
            "../data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
        save_object(
            "../data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end