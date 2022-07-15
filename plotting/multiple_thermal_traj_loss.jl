include("../src/main.jl")

# Non-thermal trajectories
no_motion = [load_object("data/Non_Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩTnothing_τ100.jld2"),
load_object("data/Non_Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩTnothing_τ100.jld2")]

# Repulsive T = 0 trajectories
T0_rep = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩT0.0_τ100.jld2")


# Attractive T = 0 trajectories
T0_att = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩT0.0_τ100.jld2")


# Repulsive T = 0.5 trajectories
T0p5_rep = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩT0.5_τ100.jld2")

# Attractive T = 0.5 trajectories
T0p5_att = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩT0.5_τ100.jld2")


# Repulsive T = 1 trajectories
T1_rep = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩT1.0_τ100.jld2")

#Attractive T = 1 trajectories
T1_att = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩT1.0_τ100.jld2")

# Repulsive T = 2 trajectories
T2_rep = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩT2.0_τ100.jld2")

#Attractive T = 2 trajectories
T2_att = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩT2.0_τ100.jld2")

# Repulsive T = 5 trajectories
T5_rep = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩT5.0_τ100.jld2")

#Attractive T = 5 trajectories
T5_att = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩT5.0_τ100.jld2")

# Repulsive T = 10.0 trajectories
T10_rep = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩT10.0_τ100.jld2")

# Attractive T = 10.0 trajectories
T10_att = load_object("data/Thermal/Single_σ0[220]_σdot0[40]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩT10.0_τ100.jld2")

data_rep = [T10_rep, T5_rep, T2_rep, T1_rep, T0p5_rep, T0_rep, no_motion[1]]
data_att = [T10_att, T5_att, T2_att, T1_att, T0p5_att, T0_att, no_motion[2]]
colors = [my_vermillion, my_orange, my_yellow, my_red, my_green, my_blue, my_black]

# ind = 1
# data_rep = [data_rep[ind]]
# data_att = [data_att[ind]]
# colors = [colors[ind]]
resonance_speed(n) = (2 * data_rep[1].ωmax * data_rep[1].α) / ((2 * n) + 1)

## Plotting
# fig = Figure(resolution = (1400, 1600), font = "CMU Serif", fontsize = 36, figure_padding = 30)
# ax1 = Axis(fig[1,1], xlabel = L"\dot{\sigma}", ylabel = L"\Delta", title = "Repulsive")
# ax2 = Axis(fig[2,1], xlabel = L"\dot{\sigma}", ylabel = L"\Delta", title = "Attractive")
#
# # Repulsive case
# for ii in 1:length(data_rep)
#     (xs, ys) = Δ_traj(data_rep[ii])
#     scatter!(ax1, xs, ys, color = colors[ii], markersize = 15, label = L"\omega_T = %$(data_rep[ii].ωT)")
# end
#
# σ_dots = range(10, 60, length=100)
# Δs = Δ_analytic.(σ_dots, T0_rep.Φ, T0_rep.λ, T0_rep.ωmax)
# lines!(ax1, σ_dots, Δs, linewidth=4, color = my_black, label = "Analytic")
# vlines!(ax1, resonance_speed.(8:20), color = my_black, linestyle = :dash, linewidth = 3)
#
# # Attractive case
# for ii in 1:length(data_att)
#     (xs, ys) = Δ_traj(data_att[ii])
#     scatter!(ax2, xs, ys, color = colors[ii], markersize = 15, label = L"\omega_T = %$(data_att[ii].ωT)")
# end
# σ_dots = range(10, 60, length=100)
# Δs = Δ_analytic.(σ_dots, T0_att.Φ, T0_att.λ, T0_att.ωmax)
# lines!(ax2, σ_dots, Δs, linewidth=4, color = my_black, label = "Analytic")
# vlines!(ax2, resonance_speed.(8:20), color = my_black, linestyle = :dash, linewidth = 3)
#
# axislegend(ax1, labelsize = 24, unique = true, position = :lt)
# axislegend(ax2, labelsize = 24, unique = true, position = :lt)
#
# xlims!(ax1, 0, 50)
# # ylims!(ax1, -1.5, 1.5)
#
# xlims!(ax2, 0, 50)
# # ylims!(ax2, -1.5, 1.5)
# fig


## Root Mean Squared Deviation
function rmsd_Delta(data)
    (xs, ys) = Δ_traj(data)
    indices = findall(x -> x >= 10, xs)
    (xs, ys) = (xs[indices], ys[indices])
    deltas = Δ_analytic.(xs, data.Φ, data.λ, data.ωmax)
    return rmsd(ys, deltas)
end

xs_rmsd = [10.0, 5.0, 2.0, 1.0, 0.5, 0.0]
rep_rmsd = map(rmsd_Delta, data_rep)
att_rmsd = map(rmsd_Delta, data_att)

fig = Figure(resolution = (1400, 1600), font = "CMU Serif", fontsize = 36, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\omega_T", ylabel = "RMSD from Analytics", title = "Repulsive")
ax2 = Axis(fig[2,1], xlabel = L"\omega_T", ylabel = "RMSD from Analytics", title = "Attractive")

scatterlines!(ax1, xs_rmsd, rep_rmsd[1:(end-1)], color = my_blue, linewidth = 3, markersize = 15)
hlines!(ax1, [rep_rmsd[end]], linestyle = :dash, color = my_black, linewidth = 4)
text!(ax1, "Non-Thermal", position = (8.5, 0.2), textsize = 30)

scatterlines!(ax2, xs_rmsd, att_rmsd[1:(end-1)], color = my_blue, linewidth = 3, markersize = 15)
hlines!(ax2, [att_rmsd[end]], linestyle = :dash, color = my_black, linewidth = 4)
text!(ax2, "Non-Thermal", position = (8.5, 0.2), textsize = 30)

xlims!(ax1, -0.5, 10.5)
ylims!(ax1, 0, nothing)

xlims!(ax2, -0.5, 10.5)
ylims!(ax2, 0, nothing)
fig
