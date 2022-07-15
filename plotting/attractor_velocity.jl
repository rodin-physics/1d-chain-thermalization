include("../src/main.jl")

α = 40
μ = 1
nChain = 100
ωmax = 10
Φ0 = 8.0
function Δ_numeric(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.25 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    σs = res.σs |> vec
    ρs = res.ρs

    chain_idx = searchsortedlast(ρs[:,1], σs[1])
    mod_val = mod(σs[1], res.α)

    mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
    mob_final = findmin(abs.(σs .- (ρs[chain_idx+1, mid_pt_idx] + mod_val)))[2]

    # v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
    v_final = (σs[mob_final] - σs[mob_final - 1]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

function Δ_numeric2(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.25 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    return res
end

system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
## Plotting
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=34, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta")

resonance_speed(n) = (2 * α * ωmax) / ((2 * n) + 1)

xs = range(10, 150, step = 1.0)
λs = [4.0]
numeric_rep = map(x -> Δ_numeric(x, 5.5 * α, Φ0, λs[1], system), xs)
numeric_att = map(x -> Δ_numeric(x, 5.5 * α, -Φ0, λs[1], system), xs)

colors = reverse([my_blue, my_vermillion, my_green])

for ii in 1:length(λs)
    Δs = Δ_analytic.(xs, Φ0, λs[ii], ωmax)
    lines!(ax1, xs, Δs, color = colors[ii], linewidth = 5, label = L"\lambda = %$(λs[ii])")
end
hlines!(ax1, [Δ_analytic(50, Φ0, λs[1], ωmax)], linestyle = :dash, color = my_black, linewidth = 3)
text!(L"V_\mathrm{bias}", position  = (180, 0.22), textsize = 30)

# points = findall(x -> isapprox(x, Δ_analytic(50, 8.0, 4.0, 10), atol = 1e-4), Δ_analytic.(xs, 8.0, 4.0, 10))
# scatter!(ax1, [14.5, 50], Δ_analytic.([14.5, 50], 8.0, 4.0, 10), color = my_vermillion, markersize = 18)
scatter!(ax1, xs, numeric_rep, color = my_vermillion, markersize = 4, label = "Repulsive")
scatter!(ax1, xs, numeric_att, color = my_blue, markersize = 4, label = "Attractive")
# vlines!(ax1, [42, 50])
# vlines!(ax1, resonance_speed.(0:20), color = my_black, linewidth = 2, linestyle = :dashdot)

xlims!(ax1, 0, 155)
ylims!(ax1, -0.1, 0.45)
axislegend(ax1, labelsize = 20)
fig
