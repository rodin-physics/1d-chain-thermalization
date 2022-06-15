include("../src/main.jl")
system_slow = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
system_fast = load_object("precomputed/systems/System_ωmax10_d6000_l10.jld2")
colors = [my_vermillion, my_orange, my_green, my_sky, my_blue]

## System parameters
α = 10                  # Distance between chain atoms
μ = 1                   # Dimensionless mass
σ0 = 4.5 * α            # Mobile particle initial position
nChain = 10             # Number of chain atoms
vPts = 100
λs = [1 / 8, 1 / 4, 1 / 2, 1]

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

## Analytic energy loss of mobile particle per encounter with chain atom
function Δ_analytic(v, Φ, λ, Ω)
    z = (2 * π * λ / v)^2
    return (
        4 * π^3 * Φ^2 / v^2 *
        z *
        exp(-z * (Ω^2 + 1) / 2) *
        (
            besseli(0, z * (Ω^2 - 1) / 2) +
            (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
        )
    )
end

## Low speed
d = 60
Φ0 = 0.025
σ_dots = range(2, 20, length = vPts)

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}_0",
    ylabel = L"\Delta",
    yscale = log10,
    # xscale = log10,
    # yminorticksvisible = true,
    # yminorgridvisible = true,
    # yminorticks = IntervalsBetween(10),
    # xminorticksvisible = true,
    # xminorgridvisible = true,
    # xminorticks = IntervalsBetween(10),
)


for ii = 1:length(λs)

    λ = λs[ii]
    res_att = [Δ_numeric(x, σ0, -Φ0, λ, system_slow) for x in σ_dots]
    res_rep = [Δ_numeric(x, σ0, Φ0, λ, system_slow) for x in σ_dots]
    analytic = Δ_analytic.(σ_dots, Φ0, λ, system_slow.ωmax)

    scatter!(ax1, (σ_dots), (res_att), color = colors[ii], marker = :hline, markersize = 20)
    scatter!(ax1, (σ_dots), (res_rep), color = colors[ii], marker = :cross, markersize = 20)
    lines!(
        ax1,
        (σ_dots),
        (analytic),
        color = colors[ii],
        linewidth = 4,
        label = L"\lambda = %$λ",
    )
end
axislegend(ax1, position = :rb)
save("Slow_dissipation.pdf", fig)


## High speed
d = 60
Φ0 = 2
σ_dots = range(20, 350, length = vPts)

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
ax1 = Axis(
    fig[1, 1],
    xlabel = L"\dot{\sigma}_0",
    ylabel = L"\Delta",
    yscale = log10,
    # xscale = log10,
    # yminorticksvisible = true,
    # yminorgridvisible = true,
    # yminorticks = IntervalsBetween(10),
    # xminorticksvisible = true,
    # xminorgridvisible = true,
    # xminorticks = IntervalsBetween(10),
)


for ii = 1:length(λs)

    λ = λs[ii]
    res_att = [Δ_numeric(x, σ0, -Φ0, λ, system_fast) for x in σ_dots]
    res_rep = [Δ_numeric(x, σ0, Φ0, λ, system_fast) for x in σ_dots]
    analytic = Δ_analytic.(σ_dots, Φ0, λ, system_fast.ωmax)

    scatter!(ax1, (σ_dots), (res_att), color = colors[ii], marker = :hline, markersize = 20)
    scatter!(ax1, (σ_dots), (res_rep), color = colors[ii], marker = :cross, markersize = 20)
    lines!(
        ax1,
        (σ_dots),
        (analytic),
        color = colors[ii],
        linewidth = 4,
        label = L"\lambda = %$λ",
    )
end
axislegend(ax1, position = :rt)
fig
save("Fast_dissipation.pdf", fig)