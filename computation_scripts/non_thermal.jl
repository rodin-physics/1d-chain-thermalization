using Distributed

proc_num = 9
addprocs(proc_num - nprocs())

include("../src/main.jl")
# @everywhere include("../src/main.jl")

system = load_object("precomputed/systems/System_K1_k20_m1_d60_l500.jld2")
gg = system.G

@time gg[50000][200]
d = 60

τ = 50                              # Simulation time
Ωmin = √(system.K / system.m)       # Minimum phonon frequency
tmin = 2 * π / Ωmin                 # Period of the slowest mode
δ = system.δ                        # Time step
n_pts = floor(τ * tmin / δ) |> Int  # Number of time steps
nChain = 50                         # Number of chain atoms tracked
a = 1                               # Distance between chain atoms

rHs = [a .* collect(1:nChain) for n = 1:n_pts]
tTraj = ThermalTrajectory(system.k, system.K, system.m, a, system.δ, rHs, nothing, ħ)

# Modify the system variable to reduce its size for parallelization
m = system.m            # Mass of the chain atoms
k = system.k            # Spring force constant
K = system.K            # Confining potential force constant
δ = system.δ            # Array of time steps
G_list = system.G       # Memory term
G_list = G_list[1:n_pts]
G_list = [x[1:nChain] for x in G_list]

system = ChainSystem(k, K, m, δ, G_list)

## Width and Depth Dependence
s = [1 / 2, 1 / 4, 1 / 8, 1 / 16]
F = [-1 / 2, -1 / 4, -1 / 8, -1 / 16, -1 / 32, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2]
# F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]
v0 = [[1 / 2]]
M = [1]
mems = [Inf]

parameters = [(w, x, y, z, mem) for w in M, x in s, y in F, z in v0, mem in mems] |> vec
x0 = [5.5]
@showprogress pmap(parameters) do param
    M = param[1]
    s = param[2]
    F = param[3]
    v0 = param[4]
    mem = param[5]
    if (
        !isfile(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
        save_object(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

## Speed dependence
s = [1 / 8, 1 / 4]
F = [1 / 8, 1 / 4, 1 / 2, 1, 2]
v0 = [[1 / 2], [1], [2], [4], [5]]
mems = [Inf]
M = [1]


parameters = [(w, x, y, z, mem) for w in M, x in s, y in F, z in v0, mem in mems] |> vec
x0 = [5.5]

@showprogress pmap(parameters) do param
    M = param[1]
    s = param[2]
    F = param[3]
    v0 = param[4]
    mem = param[5]
    if (
        !isfile(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
        save_object(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

s = [1 / 8]
F = [1 / 16, 1 / 8, 1 / 4, 1 / 2, 1, 2]
v0 = [[4], [6]]
mems = [Inf]
M = [Inf]


parameters = [(w, x, y, z, mem) for w in M, x in s, y in F, z in v0, mem in mems] |> vec
x0 = [5.5]

@showprogress pmap(parameters) do param
    M = param[1]
    s = param[2]
    F = param[3]
    v0 = param[4]
    mem = param[5]
    if (
        !isfile(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
        save_object(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

## Memory Dependence
s = [1 / 4]
F = [-1 / 4, 1 / 4]
v0 = [[1 / 2]]
M = [1 / 2]
mems = [0, 0.05, 0.5, 1, 10, 50, Inf]

parameters = [(w, x, y, z, mem) for w in M, x in s, y in F, z in v0, mem in mems] |> vec
x0 = [5.5]

@showprogress pmap(parameters) do param
    M = param[1]
    s = param[2]
    F = param[3]
    v0 = param[4]
    mem = param[5]
    if (
        !isfile(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )
        res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
        save_object(
            "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end


M = 1.0
s = 0.5
F = -1
v0 = [1/2]
mem = Inf
x0 = [5.0]
if (
    !isfile(
        "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
    )
)
    res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)
    save_object(
        "data/Non_Thermal/Single_x0$(x0)_v0$(v0)_Mem$(mem)_s$(s)_F$(F)_M$(M)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        res,
    )
end