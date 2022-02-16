include("main.jl")
system = load_object("precomputed/systems/System_K1_k20_m1_d60_l500.jld2")

d = 60

τ = 30                              # Simulation time
Ωmin = √(system.K / system.m)       # Minimum phonon frequency
tmin = 2 * π / Ωmin                 # Period of the slowest mode
δ = system.δ                        # Time step
n_pts = floor(τ * tmin / δ) |> Int  # Number of time steps
nChain = 20                         # Number of chain atoms tracked
a = 1                               # Distance between chain atoms

rHs = [a .* collect(1:nChain) for n = 1:n_pts]
tTraj = ThermalTrajectory(system.k, system.K, system.m, a, system.δ, rHs, nothing, ħ)

## Width and Depth Dependence
s = 1/4
F = 2/5
v0 = [3 / 5]
M = 1/2
mem = Inf
x0 = [3]

res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ);
res_NEW = motion_solver_NEW(system, F, s, M, x0, v0, tTraj, mem, τ);
res.Rs ≈ res_NEW.Rs
res.rs[4] ≈ res_NEW.rs[4]
res.rs[3]
res_NEW.rs[3]
plot(res.ts, [x[1] for x in res.Rs])
# lines( res.ts, [x[1] for x in res.Rs])

# fig = Figure()
# ax = Axis(fig[1, 1])
# l1 = lines!(ax, res.ts, [x[1] for x in res.Rs])
# # l2 = lines!(ax, res.ts, [x[2] for x in res.Rs] , color = my_red)
# # l3 = lines!(ax, res.ts, [x[3] for x in res.Rs] , color = my_green)
# # l3 = lines!(ax, res.ts, [x[4] for x in res.Rs] , color = my_orange)
# # for ii = 1:nChain
# #     lines!(ax, res.ts, [x[ii] for x in res.rs], color = colorant"rgba(0, 0, 0, 0.35)")
# # end
# # l2 = lines!(ax, (1:length(Rs)) .* δ, [x[2] for x in res.Rs] .- 0, color = my_red)
# # l2 = lines!(ax, (1:length(Rs)) .* δ, [x[3] for x in res.Rs] .- 0, color = my_green)
# fig
