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

## Plotting
fig = Figure(resolution = (1400, 1600), font = "CMU Serif", fontsize = 36, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\sigma - \sigma_0", title = "Repulsive")
ax2 = Axis(fig[2,1], xlabel = L"\tau", ylabel = L"\sigma - \sigma_0", title = "Attractive")

## REPULSIVE
# T = 10
for ii in 1:length(data_rep)
    y_vals = data_rep[ii].σs |> vec
    lines!(
        ax1,
        data_rep[ii].τs,
        (y_vals .- y_vals[1]),
        color = colors[ii],
        linewidth = 5,
        label = L"\omega_T = %$(data_rep[ii].ωT)",
    )
end

for ii in 1:length(data_att)
    y_vals = data_att[ii].σs |> vec
    lines!(
        ax2,
        data_att[ii].τs,
        (y_vals .- y_vals[1]),
        color = colors[ii],
        linewidth = 5,
        label = L"\omega_T = %$(data_att[ii].ωT)",
    )
end

CairoMakie.xlims!(ax1, 0, 100)
CairoMakie.ylims!(ax1, 0, 4000)

CairoMakie.xlims!(ax2, 0, 100)
CairoMakie.ylims!(ax2, 0, 4000)

axislegend(ax1, labelsize = 34, unique = true, position = :lt)
axislegend(ax2, labelsize = 34, unique = true, position = :lt)
fig


## TRAJECTORY ENERGY LOSS

# fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}_0", ylabel = L"\Delta")
#
# # REPULSIVE
# for data in T0_rep
#     res_numeric = Δ_traj(data)
#
#     scatter!(
#         ax1,
#         res_numeric[1],
#         res_numeric[2],
#         marker = :cross,
#         markersize = 20,
#         color = my_vermillion,
#         label = L"\omega_T = 0.0"
#         )
# end
#
# # ATTRACTIVE
# for data in T0_att
#     res_numeric = Δ_traj(data)
#     scatter!(
#         ax1,
#         res_numeric[1],
#         res_numeric[2],
#         marker = :hline,
#         markersize = 20,
#         color = my_blue,
#         label = L"\omega_T = 0.0"
#         )
# end
# # NON-THERMAL
# for data in [no_motion[1]]
#     res_numeric = Δ_traj(data)
#     scatter!(
#         ax1,
#         res_numeric[1],
#         res_numeric[2],
#         marker = :cross,
#         markersize = 20,
#         color = my_black,
#         )
# end
# for data in [no_motion[2]]
#     res_numeric = Δ_traj(data)
#     scatter!(
#         ax1,
#         res_numeric[1],
#         res_numeric[2],
#         marker = :hline,
#         markersize = 20,
#         color = my_black,
#         label = "Non-thermal"
#         )
# end
#
# # ANALYTIC
# # vs = range(res_numeric[1][1], res_numeric[1][end], length = 1000)
# vs = range(1, 40, length = 100)
# data = T0_att[1]
# res_analytic = Δ_analytic.(vs, data.Φ, data.λ, data.ωmax)
#
# lines!(ax1, vs, res_analytic, linewidth = 4, color = my_black, label = "Analytic")
# axislegend(ax1, labelsize = 25, position = :rb, unique = false, merge = true)
# xlims!(ax1, 0, 26)
# fig

##
