include("../src/main.jl")

α = 10
μ = 1
σ0 = 45
nChain = 10
vPts = 100
λs = [1 / 8, 1 / 4, 1 / 2, 1]
system_slow = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
system_fast = load_object("precomputed/systems/System_ωmax10_d6000_l10.jld2")
colors = [my_vermillion, my_orange, my_green, my_sky, my_blue]

function Δ_numeric(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.25 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ)

    σs = res.σs |> vec
    mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
    v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

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

# Low speed
d = 60
Φ0 = 0.025
σ_dots = range(2, 20, length=vPts)

fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=36)
ax1 = Axis(
    fig[1, 1],
    xlabel=L"\dot{\sigma}_0",
    ylabel=L"\Delta",
    yscale=log10,
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

    scatter!(ax1, (σ_dots), (res_att), color=colors[ii], marker=:hline, markersize=20)
    scatter!(ax1, (σ_dots), (res_rep), color=colors[ii], marker=:cross, markersize=20)
    lines!(
        ax1,
        (σ_dots),
        (analytic),
        color=colors[ii],
        linewidth=4,
        label=L"\lambda = %$λ",
    )
end
axislegend(ax1, position=:rb)
save("Slow_dissipation.pdf", fig)


# High speed
d = 60
Φ0 = 2
σ_dots = range(20, 350, length=vPts)

fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=36)
ax1 = Axis(
    fig[1, 1],
    xlabel=L"\dot{\sigma}_0",
    ylabel=L"\Delta",
    yscale=log10,
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

    scatter!(ax1, (σ_dots), (res_att), color=colors[ii], marker=:hline, markersize=20)
    scatter!(ax1, (σ_dots), (res_rep), color=colors[ii], marker=:cross, markersize=20)
    lines!(
        ax1,
        (σ_dots),
        (analytic),
        color=colors[ii],
        linewidth=4,
        label=L"\lambda = %$λ",
    )
end
axislegend(ax1, position=:rt)
fig
save("Fast_dissipation.pdf", fig)

## TRAJECTORY ENERGY LOSS

fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=36)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta")

# REPULSIVE
data = load_object(
    "data/non_thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ0.5_μ1_d60_ΩTnothing_τ120.jld2",
)
δ = data.τs[2] - data.τs[1]
σs = data.σs |> vec
σ_dots = (σs[2:end] - σs[1:end-1]) ./ δ
# Get the data before the first turning point
idx = findfirst(σ_dots .< 0)
idx = isnothing(idx) ? length(data.τs) - 1 : idx
ts = data.τs[1:idx]
σ_dots = σ_dots[1:idx]

findlocalmaxima(σ_dots)
# Get all the indices when the particle is at the minimum velocity
σ_dots_max = prepend!(maxima(σ_dots),20)
idx = maxima(σ_dots)
idx = prepend!(idx, 1)
# # Compute the kinetic energy at the maximum 
σ_dots_max = σ_dots[idx]
KE = 0.5 * (σ_dots_max ./ 2 ./ π) .^ 2
Δs = KE[1:end-1] - KE[2:end]
# lines!(ax1, σ_dots_max, KE)
# lines!(ax1, σ_dots_max[1:end-1], Δs)
scatter!(ax1, σ_dots_max[1:end-1], Δs, marker='⊕', markersize=20)
fig
# σ_dots_burst = σ_dots_max[findall(σ_dots_max .> 10)]
# Δs = [Δ_numeric(x, σs[1], data.Φ, data.λ, system_slow) for x in σ_dots_burst]
# scatter!(ax1, σ_dots_burst, Δs, marker=:cross, markersize=20)
# # ATTRACTIVE

# data = load_object(
#     "data/non_thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ-1.0_μ1_d60_ΩTnothing_τ120.jld2",
# )
# δ = data.τs[2] - data.τs[1]
# σs = data.σs |> vec
# σ_dots = (σs[2:end] - σs[1:end-1]) ./ δ
# # Get the data before the first turning point
# idx = findfirst(σ_dots .< 0)
# ts = data.τs[1:idx]
# σ_dots = σ_dots[1:idx]

# # Get all the indices when the particle is at the maximum velocity
# idx = findall(n -> σ_dots[n] < σ_dots[n-1] && σ_dots[n] < σ_dots[n+1], 2:length(σ_dots)-1)
# idx = prepend!(idx, 1)

# # Compute the kinetic energy at the maximum 
# σ_dots_max = σ_dots[idx]
# KE = 0.5 * (σ_dots_max ./ 2 ./ π) .^ 2
# Δs = KE[1:end-1] - KE[2:end]
# scatter!(ax1, σ_dots_max[1:end-1], Δs, marker='⊖', markersize=20)

# σ_dots_burst = σ_dots_max[findall(σ_dots_max .> 10)]
# Δs = [Δ_numeric(x, σs[1], data.Φ, data.λ, system_slow) for x in σ_dots_burst]
# scatter!(ax1, σ_dots_burst, Δs, marker=:hline, markersize=20)
# ## ANALYTIC
# σ_dots = range(10, 20, length=100)
# Δs = Δ_analytic.(σ_dots, data.Φ, data.λ, data.ωmax)
# lines!(ax1, σ_dots, Δs, linewidth=4)

# xlims!(ax1, (10, 20))
# ylims!(ax1, (0, 0.15))

# fig
