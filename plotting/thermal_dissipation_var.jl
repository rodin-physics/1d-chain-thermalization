## DISTRIBUTION OF DELTA FOR ENERGY RANGE
include("../src/main.jl")
## System parameters
system_slow = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
# system_fast = load_object("precomputed/systems/System_ωmax10_d6000_l20.jld2")
filenames = ["precomputed/rH/rH_ωmax10_d60_ωT0.0_τ200.jld2",
             # "precomputed/rH/rH_ωmax10_d60_ωT1.0_τ200.jld2",
             # "precomputed/rH/rH_ωmax10_d60_ωT5.0_τ200.jld2",
             # "precomputed/rH/rH_ωmax10_d60_ωT10.0_τ200.jld2",
             # "precomputed/rH/rH_ωmax10_d60_ωT20.0_τ200.jld2"
             ]

α = 10                  # Distance between chain atoms
μ = 1                   # Dimensionless mass
σ0 = 4.5 * α            # Mobile particle initial position
nChain = 10             # Number of chain atoms
vPts = 200
λ = 0.5
ωTs = [0.0, 1.0, 5.0, 10.0, 20.0]

colors = [my_vermillion, my_orange, my_green, my_sky, my_blue]

## Numerically calculate the energy loss for given system at zero temperature
function Δ_numeric(σ_dot, σ0, Φ0, λ, system)
    # Initialize parameters
    δ = system.δ
    τ = 1.25 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    # Solve mobile particle trajectory
    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ)
    σs = res.σs |> vec

    # Find index of particle at the next midpoint
    mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
    # Calculate speed and K.E. of particle at next midpoint
    v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
    kinetic = (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
    return kinetic
end

## Numerically calculate the energy loss for given system at given T
function Δ_numeric(σ_dot, σ0, Φ0, λ, system, tTraj)
    # Initialize parameters
    δ = system.δ
    τ = 1.25 * (α / σ_dot)

    # Solve mobile particle trajectory
    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ)
    σs = res.σs |> vec

    # Find index of particle at the next midpoint
    mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]

    # Calculate speed and K.E. of particle at next midpoint
    v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
    kinetic = (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
    return kinetic
end

d = 60
Φ0 = 1.0

## Plotting range of Delta+σ0 for single energy
speeds = range(5, 100, step = 0.2)
σ_dots = range(5, 100, length = vPts)
σ0s = range(2.5, 150.5, step = 1.0) .* α
tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT0.0_τ200.jld2")

# Calculation
function mean_var_per_speed(speed, σ0s, Φ0, λ, system, tTraj)
    res = [Δ_numeric(speed, x, Φ0, λ, system, tTraj) for x in σ0s]
    return (mean(res), var(res))
end

println("Calculating attractive potential")
res_att = map(v -> mean_var_per_speed(v, σ0s, -Φ0, λ, system_slow, tTraj), speeds)
println("Calculating repulsive potential")
res_rep = map(v -> mean_var_per_speed(v, σ0s, Φ0, λ, system_slow, tTraj), speeds)

mean_att = first.(res_att)
var_att = last.(res_att)
mean_rep = first.(res_rep)
var_rep = last.(res_rep)

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\dot{\sigma}_0", ylabel = L"\Delta", title = L"\omega_T = %$(tTraj.ωT), \Phi_0 = %$(Φ0)")

# for speed in speeds
#     res_att = [Δ_numeric(speed, x, -Φ0, λ, system_slow, tTraj) for x in σ0s]
#     res_rep = [Δ_numeric(speed, x, Φ0, λ, system_slow, tTraj) for x in σ0s]
#
#     scatter!(ax1, repeat([speed], length(σ0s)), (res_att), color = my_red, marker = :hline, markersize = 15)
#     scatter!(ax1, repeat([speed], length(σ0s)), (res_rep), color = my_blue, marker = :cross, markersize = 15)
# end

# No motion case
# zero_T_att = [Δ_numeric(x, σ0, -Φ0, λ, system_slow) for x in σ_dots]
# zero_T_rep = [Δ_numeric(x, σ0, Φ0, λ, system_slow) for x in σ_dots]
#
# scatter!(ax1, σ_dots, zero_T_att, color = my_vermillion, marker = :hline, markersize = 15, label = "No motion")
# scatter!(ax1, σ_dots, zero_T_rep, color = my_vermillion, marker = :cross, markersize = 15)
band!(ax1, speeds, mean_att .- var_att, mean_att .+ var_att, color = RGBA(my_blue2, 0.3))
band!(ax1, speeds, mean_rep .- var_rep, mean_rep .+ var_rep, color = RGBA(my_red2, 0.3))

println("Calculating analytic curve")
analytic_rep = Δ_analytic.(σ_dots, Φ0, λ, system_slow.ωmax)
analytic_att = Δ_analytic.(σ_dots, -Φ0, λ, system_slow.ωmax)
lines!(ax1, σ_dots, analytic_rep, color = my_vermillion, linewidth = 4)
lines!(ax1, σ_dots, analytic_att, color = my_vermillion, linewidth = 4, label = "Analytic")

# Attractive potential
lines!(ax1, speeds, mean_att, color = my_blue, linewidth = 4, label = "Attractive")

# Repulsive potential
lines!(ax1, speeds, mean_rep, color = my_red, linewidth = 4, label = "Repulsive")


axislegend(ax1, labelsize = 20)
CairoMakie.xlims!(0.0, 100)

fig
# save("Delta_Single_Pass_$(tTraj.ωT).pdf", fig)
