## SINGLE PASS DISSIPATION FOR THERMAL MOTION
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
Φ0 = 2.0

## Plotting over energy range
σ_dots = range(20, 350, length = vPts)

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36, yscale = Makie.pseudolog10)
ax1 = Axis(fig[1,1], xlabel = L"\dot{\sigma}_0", ylabel = L"\Delta")

for ii in 1:length(filenames)
    tTraj = load_object(filenames[ii])
    ωT = tTraj.ωT
    println("ωT is ", ωT)

    res_att = [Δ_numeric(x, σ0, -Φ0, λ, system_slow, tTraj) for x in σ_dots]
    res_rep = [Δ_numeric(x, σ0, Φ0, λ, system_slow, tTraj) for x in σ_dots]

    scatter!(ax1, (σ_dots), (res_att), color = colors[ii], marker = :hline, markersize = 15, label = L"\omega_T = %$(ωT)")
    scatter!(ax1, (σ_dots), (res_rep), color = colors[ii], marker = :cross, markersize = 15)

end

zero_T_att = [Δ_numeric(x, σ0, -Φ0, λ, system_slow) for x in σ_dots]
zero_T_rep = [Δ_numeric(x, σ0, Φ0, λ, system_slow) for x in σ_dots]

scatter!(ax1, σ_dots, zero_T_att, color = my_black, marker = :hline, markersize = 15, label = "No motion")
scatter!(ax1, σ_dots, zero_T_rep, color = my_black, marker = :cross, markersize = 15)

analytic = Δ_analytic.(σ_dots, Φ0, λ, system_slow.ωmax)
lines!(ax1, σ_dots, analytic, color = my_black, linewidth = 4)
axislegend(ax1, labelsize = 20)
CairoMakie.xlims!(0.0, 350)
# save("Fast_system.pdf", fig)
fig

## Plotting range of Delta+σ0 for single energy
# speed = 100
# σ0s = range(2.5, 100.5, step = 1.0) .* α
# tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT0.0_τ200.jld2")
#
# fig2 = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
# ax2 = Axis(fig2[1,1], xlabel = L"\sigma_0", ylabel = L"\Delta", title = L"\dot{\sigma}_0 = %$speed")
#
# res_att = [Δ_numeric(speed, x, -Φ0, λ, system_slow, tTraj) for x in σ0s]
# res_rep = [Δ_numeric(speed, x, Φ0, λ, system_slow, tTraj) for x in σ0s]
#
# scatter!(ax2, (σ0s), (res_att), color = my_red, marker = :hline, markersize = 15)
# scatter!(ax2, (σ0s), (res_rep), color = my_blue, marker = :cross, markersize = 15)
#
# hlines!(ax2, [Δ_analytic(speed, Φ0, λ, system_slow.ωmax), mean(res_att), mean(res_rep)], color = my_black, linestyle = :dashdot)
# fig2
