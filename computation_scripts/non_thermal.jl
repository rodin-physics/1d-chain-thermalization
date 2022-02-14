include("../src/main.jl")

system = load_object("precomputed/systems/System_K1_k20_m1_d60_l500.jld2")

d = 60

M = 1 / 2                           # Mass of the mobile particles

s = 1 / 4
F = 1

x0 = [0.5,1.5, 10.5, 15.5]
v0 = [0, 3.0, 0, 0]
# x0 = [3, 4, 5]
# v0 = [2.5, 2.5, 2.5]
mem = 1
τ = 40

Ωmin = √(system.K / system.m)                 # Minimum phonon frequency
tmin = 2 * π / Ωmin                 # Period of the slowest mode
δ = system.δ
n_pts = floor(τ * tmin / δ) |> Int  # Number of time steps
nChain = 20
a = 1
rHs = [a .* collect(0:nChain) for n = 1:n_pts]

tTraj = ThermalTrajectory(system.k, system.K, system.m, a, system.δ, rHs, nothing, ħ)
res = motion_solver(system, F, s, M, x0, v0, tTraj, mem, τ)


fig = Figure()
ax = Axis(fig[1, 1])
l1 = lines!(ax, res.ts, [x[1] for x in res.Rs] )
l2 = lines!(ax, res.ts, [x[2] for x in res.Rs] , color = my_red)
l3 = lines!(ax, res.ts, [x[3] for x in res.Rs] , color = my_green)
l3 = lines!(ax, res.ts, [x[4] for x in res.Rs] , color = my_orange)
for ii = 1:nChain
    lines!(
        ax,
        res.ts,
        [x[ii] for x in res.rs] ,
        color = colorant"rgba(0, 0, 0, 0.35)",
    )
end
# l2 = lines!(ax, (1:length(Rs)) .* δ, [x[2] for x in res.Rs] .- 0, color = my_red)
# l2 = lines!(ax, (1:length(Rs)) .* δ, [x[3] for x in res.Rs] .- 0, color = my_green)
fig

G(2,[9], K, k, m)


δ * 4.5
scatter(δ .* ( 1:700),(Rs[2:end]-Rs[1:end-1])[1:700]./δ)