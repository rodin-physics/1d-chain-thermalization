include("../src/main.jl")

## Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 8.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")

## Load data
speed = 80
# bias = 0.21465496799193154
bias = 0.0
data_rep = [load_object("Data/Non_Thermal_Multi/$(n)_Single_σ0[220]_σdot0[$(speed)]_MemInf_λ4_Φ8_μ1_d300_bias$(bias)_ΩTnothing_τ30.jld2") for n in 1:10]

data_att = [load_object("Data/Non_Thermal_Multi/$(n)_Single_σ0[220]_σdot0[$(speed)]_MemInf_λ4_Φ-8_μ1_d300_bias$(bias)_ΩTnothing_τ30.jld2") for n in 1:10]

## Functions
resonance_speed(n) = (2 * ωmax * α) / ((2 * n) + 1)

function Δ_numeric(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 2.2 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    σs = res.σs |> vec
    ρs = res.ρs

    chain_idx = searchsortedlast(ρs[:,1], σs[1])
    mod_val = mod(σs[1], res.α)

    mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
    mob_final = findmin(abs.(σs .- (ρs[chain_idx+1, :] .+ mod_val)))[2]

    # v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
    v_final = (σs[mob_final+1] - σs[mob_final]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

function Δ_numeric2(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.8 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    return res
end


start_frac_rep = 0.5
start_frac_att = 0.5

## Plotting Delta
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=34, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta", title = "Start points - Rep: $(start_frac_rep), Att: $(start_frac_att)")

xs = range(10, 110, step = 0.05)
xs2 = range(20, 55, step = 0.001)

# Analytic Δ
Δs = Δ_analytic.(xs, Φ0, λ, ωmax)
Δnews = Δ_transport.(xs2, Φ0, λ, ωmax, α)
lines!(ax1, xs, Δs, color = my_green, linewidth = 5, label = L"\lambda = %$(λ)")
lines!(ax1, xs2, Δnews, color = my_black, label = L"\Delta_\mathrm{transport}")

# vlines!(ax1, [sqrt(8*Φ0*π^2/μ), 41,50])
# vlines!(ax1, [35, 40])
# vlines!(ax1, resonance_speed.(1:30), color = my_black, linewidth = 2, linestyle = :dash)
# hlines!(ax1, [Δ_analytic(50, Φ0, λ, ωmax)], linestyle = :dash, color = my_black, linewidth = 3)

# Numeric Δ
# numeric_rep = map(x -> Δ_numeric(x, (5.0 + start_frac_rep) * α, Φ0, λ, system), xs)
# numeric_att = map(x -> Δ_numeric(x, (5.0 + start_frac_att) * α, -Φ0, λ, system), xs2)

# scatter!(ax1, xs, numeric_rep, color = my_vermillion, markersize = 4, label = "Repulsive")
# scatter!(ax1, xs2, numeric_att, color = my_blue, markersize = 4, label = "Attractive")

# Plot full trajectory losses
for ii in 1:length(data_rep)
    (rep_xs, rep_ys) = Δ_traj(data_rep[ii])
    (att_xs, att_ys) = Δ_traj(data_att[ii])

    scatter!(ax1, rep_xs, rep_ys .+ data_rep[ii].bias, markersize=10, color = my_vermillion, label = "Repulsive")
    scatter!(ax1, att_xs, att_ys .+ data_rep[ii].bias, markersize=10, color = my_blue, label = "Attractive")
end


xlims!(ax1, 0, 110)
ylims!(ax1, -0.1, 1.0)
axislegend(ax1, labelsize = 20, unique = true)
fig


## Speed of particle over time
# function particle_speed(data)
#     σs = data.σs |> vec
#     τs = data.τs
#     δ = τs[2] - τs[1]
#     speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]
#
#     return (τs[2:end], speeds)
# end
#
# ## Find final velocity at midpoint
# function find_midpoint(σ0, data)
#     σs = data.σs |> vec
#     ρs = data.ρs
#
#     chain_idx = searchsortedlast(ρs[:,1], σs[1])
#     mod_val = mod(σs[1], data.α)
#
#     mid_pt_idx = findmin(abs.(σs .- (σ0 + data.α)))[2]
#     mob_final = findmin(abs.(σs .- (ρs[chain_idx+1, :] .+ mod_val)))[2]
#
#     return mid_pt_idx
# end

## Plot velocities
# fig = Figure(resolution=(1200, 1600), font="CMU Serif", fontsize=32, figure_padding = 30)
# ax1 = Axis(fig[1, 1], xlabel=L"\tau", ylabel=L"\dot{\sigma}", title = "Repulsive - $(start_frac_rep)a")
# ax2 = Axis(fig[2, 1], xlabel=L"\tau", ylabel=L"\dot{\sigma}", title = "Attractive - $(start_frac_att)a")
#
# xs_2 = range(5, 20, step = 0.5)
#
# for ii in xs_2
#     # res_rep = Δ_numeric2(ii, (5.0 + start_frac_rep) * α, Φ0, λ, system)
#     res_att = Δ_numeric2(ii, (5.0 + start_frac_att) * α, -Φ0, λ, system)
#
#     # (xs_rep, speed_rep) = particle_speed(res_rep)
#     (xs_att, speed_att) = particle_speed(res_att)
#
#     # mid_rep = find_midpoint((5.0 + start_frac_rep) * α, res_rep)+1
#     mid_att = find_midpoint((5.0 + start_frac_att) * α, res_att)+1
#
#     # lines!(ax1, xs_rep, speed_rep, linewidth = 2, color = my_vermillion, label = "Repulsive")
#     lines!(ax2, xs_att, speed_att, linewidth = 2, color = my_blue, label = "Attractive")
#
#     # scatter!(ax1, [xs_rep[mid_rep]], [speed_rep[mid_rep]], color = my_black, markersize = 10)
#     scatter!(ax2, [xs_att[mid_att]], [speed_att[mid_att]], color = my_black, markersize = 10)
# end
#
# # axislegend(ax1, unique = true, labelsize = 20)
# # xlims!(ax1, 0.0, 10.0)
# xlims!(ax2, 0.0, 10.0)
# fig
