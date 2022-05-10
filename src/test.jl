include("../src/main.jl")
# @everywhere include("../src/main.jl")

system = load_object("precomputed/systems/System_ωmax10_d60_l50.jld2")

d = 60
τ = 50                             # Simulation time
δ = system.δ                        # Time step
a = 1                               # Distance between chain atoms
n_pts = τ / δ |> floor |> Int
nChain = 50
ρHs = zeros(nChain, n_pts)
tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

res = motion_solver(system, 1, 1, 1, [5.5], [0.25], 1, tTraj, Inf, τ)
# res.σs
# aa = ones(500000, 500);
# bb = ones(500000, 500);
# cc = ones(500000, 500);

# @time aa .= aa .+ aa;

# @time for ii = 1:500
#     Threads.@spawn(view(bb, 1:100000, ii) .= view(bb, 1:100000, ii) .+3 .* view(bb, 1:100000, ii))
# end

# @time Threads.@threads for ii = 1:500
#     view(cc, 1:100000, ii) .= view(cc, 1:100000, ii) .+  view(cc, 1:100000, ii)
# end


# # @time aa .= aa .* 3;

# @time for ii in 1 : 500
#     Threads.@spawn(view(bb,:, ii) .= view(bb,:, ii) .* 3)
# end


# @time Threads.@threads for ii in 1 : 500
#     view(cc,:, ii) .=  view(cc,:, ii) .* 3
# end



# proc_num = 9
# addprocs(proc_num - nprocs())

