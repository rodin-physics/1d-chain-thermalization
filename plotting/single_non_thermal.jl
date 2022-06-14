include("../src/main.jl")

## SINGLE PASS ENERGY LOSS
slow = load_object("data/non_thermal/Single_Pass_Slow_Φ0.025_μ1.jld2")
fast = load_object("data/non_thermal/Single_Pass_Fast_Φ2.0_μ1.jld2")
colors = [my_vermillion, my_orange, my_green, my_sky]

fig = Figure(resolution=(1200, 1600), font="CMU Serif", fontsize=36)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta", yscale=log10)
ax2 = Axis(fig[2, 1], xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta", yscale=log10)
slow_keys = sort(collect(keys(slow)), by=x -> x[2])
fast_keys = sort(collect(keys(fast)), by=x -> x[2])

for ii = 1:length(slow_keys)
    k = slow_keys[ii]
    r = get(slow, k, nothing)
    scatter!(
        ax1,
        r[1],
        r[2],
        color=colors[(ii+1)÷2],
        marker=k[1] > 0 ? :cross : :hline,
        markersize=20,
    )

    if k[1] > 0
        analytic = Δ_analytic.(r[1], k[1], k[2], 10)
        lines!(
            ax1,
            r[1],
            analytic,
            color=colors[(ii+1)÷2],
            linewidth=4,
            label=L"\lambda = %$(k[2])",
        )
    end

    k = fast_keys[ii]
    r = get(fast, k, nothing)
    scatter!(
        ax2,
        r[1],
        r[2],
        color=colors[(ii+1)÷2],
        marker=k[1] > 0 ? :cross : :hline,
        markersize=20,
    )

    if k[1] > 0
        analytic = Δ_analytic.(r[1], k[1], k[2], 10)
        lines!(
            ax2,
            r[1],
            analytic,
            color=colors[(ii+1)÷2],
            linewidth=4,
            label=L"\lambda = %$(k[2])",
        )
    end


end
axislegend(ax1, position=:rb)
axislegend(ax2, position=:rt)

save("Dissipation.pdf", fig)


## TRAJECTORY ENERGY LOSS

fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=36)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta")

# REPULSIVE
data = load_object(
    "data/non_thermal/Single_σ0[55]_σdot0[30]_MemInf_λ0.5_Φ1_μ1_d60_ΩTnothing_τ250.jld2",
)

res_numeric = Δ_traj(data)

scatter!(ax1, res_numeric[1], res_numeric[2], marker=:cross, markersize=20, color=my_vermillion)

# ATTRACTIVE
data = load_object(
    "data/non_thermal/Single_σ0[55]_σdot0[30]_MemInf_λ0.5_Φ-1_μ1_d60_ΩTnothing_τ250.jld2",
)

res_numeric = Δ_traj(data)
vs = range(res_numeric[1][1], res_numeric[1][end], length=1000)
res_analytic = Δ_analytic.(vs, data.Φ, data.λ, data.ωmax)

scatter!(ax1, res_numeric[1], res_numeric[2], marker=:hline, markersize=20, color=my_blue)
lines!(ax1, vs, res_analytic,
    linewidth=4, color=my_black
)

fig
save("Trajectory_Loss.pdf", fig)


## General Example
step_size = 20

fig = Figure(resolution=(1200, 1600), font="CMU Serif", fontsize=36)
ax1 = Axis(fig[1, 1], xlabel=L"\tau", ylabel=L"\sigma")
ax2 = Axis(fig[2, 1], xlabel=L"\tau", ylabel=L"\sigma")

# REPULSIVE

data = load_object(
    "data/non_thermal/Single_σ0[55]_σdot0[30]_MemInf_λ0.5_Φ1_μ1_d60_ΩTnothing_τ250.jld2",
)
δ = data.τs[2] - data.τs[1]
τ_max = 105
idx = findall(data.τs .< τ_max)
τs = data.τs[idx]
rr = reduce(hcat, [data.ρs[ii, idx] .- ii * data.α for ii = 1:size(data.ρs)[1]])
# mx = maximum(abs.(rr)) / 1
mx = 0.03
hm = heatmap!(
    ax1,
    τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap=:RdBu,
    colorrange=(-mx, mx),
)
lines!(ax1, data.τs, [x[1] for x in data.σs] |> vec, color=my_black, linewidth=5)
xlims!(ax1, (0, 105))
ylims!(ax1, (0, 2500))

lines!(
    ax1,
    0:0.1:10,
    π * 10 * 9 * (0:0.1:10) .+ data.σs[1][1],
    color=my_black,
    linewidth=4,
    linestyle=:dash,
)
Colorbar(fig[1, 2], hm; label=L"\Delta\rho", width=15, ticksize=15, tickalign=1)
# ATTRACTIVE

data = load_object(
    "data/non_thermal/Single_σ0[55]_σdot0[30]_MemInf_λ0.5_Φ-1_μ1_d60_ΩTnothing_τ250.jld2",
)
δ = data.τs[2] - data.τs[1]
τ_max = 105
idx = findall(data.τs .< τ_max)
τs = data.τs[idx]
rr = reduce(hcat, [data.ρs[ii, idx] .- ii * data.α for ii = 1:size(data.ρs)[1]])
mx = 0.03
# mx = maximum(abs.(rr)) / 3
hm = heatmap!(
    ax2,
    τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap=:RdBu,
    colorrange=(-mx, mx),
)
lines!(ax2, data.τs, [x[1] for x in data.σs] |> vec, color=my_black, linewidth=5)
xlims!(ax2, (0, 105))
ylims!(ax2, (0, 2500))

lines!(
    ax2,
    0:0.1:10,
    π * 10 * 9 * (0:0.1:10) .+ data.σs[1][1],
    color=my_black,
    linewidth=4,
    linestyle=:dash,
)
Colorbar(fig[2, 2], hm; label=L"\Delta\rho", width=15, ticksize=15, tickalign=1)
fig
save("General_Example.pdf", fig)

