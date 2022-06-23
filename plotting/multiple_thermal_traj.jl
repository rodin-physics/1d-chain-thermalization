include("../src/main.jl")

rep = false

# Non-thermal trajectories
no_motion = [load_object("data/Non_Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩTnothing_τ100.jld2"),
load_object("data/Non_Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩTnothing_τ100.jld2")]

# Repulsive T = 0 trajectories
T0_rep = [
load_object("data/Thermal/Single_σ0[45.0]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55.0]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.0_τ100.jld2")
]

# Attractive T = 0 trajectories
T0_att = [
load_object("data/Thermal/Single_σ0[45.0]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55.0]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.0_τ100.jld2")
]

# Repulsive T = 0.5 trajectories
T0p5_rep = [
load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT0.5_τ100.jld2")
]

# Attractive T = 0.5 trajectories
T0p5_att = [
load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.5_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT0.5_τ100.jld2")
]


# Repulsive T = 1 trajectories
T1_rep = [
load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT1.0_τ100.jld2")
]

#Attractive T = 1 trajectories
T1_att = [
load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT1.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT1.0_τ100.jld2")
]

# Repulsive T = 10.0 trajectories
T10_rep = [
load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_ΩT10.0_τ100.jld2")
]
# Attractive T = 10.0 trajectories
T10_att = [
load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[65]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[75]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT10.0_τ100.jld2"),
load_object("data/Thermal/Single_σ0[85]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_ΩT10.0_τ100.jld2")
]




## Plotting
fig = Figure(resolution = (1400, 1600), font = "CMU Serif", fontsize = 36, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\sigma - \sigma_0")

if rep == true
    ## REPULSIVE
    T = 10
    for ii in T10_rep
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_vermillion, linewidth = 5, label = L"\omega_T = 10.0")
    end

    # T = 1
    for ii in T1_rep
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_red, linewidth = 5, label = L"\omega_T = 1.0")
    end

    # T = 0.5
    for ii in T0p5_rep
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_green, linewidth = 5, label = L"\omega_T = 0.5")
    end

    # T = 0
    for ii in T0_rep
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_blue, linewidth = 5, label = L"\omega_T = 0.0")
    end

    # Non-thermal motion
    for ii in [no_motion[1]]
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_black, linewidth = 7, label = "Non-thermal")
    end
else
    # ATTRACTIVE
    # T = 10
    for ii in T10_att
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_vermillion, linewidth = 5, label = L"\omega_T = 10.0")
    end

    # T = 1
    for ii in T1_att
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_red, linewidth = 5, label = L"\omega_T = 1.0")
    end

    # T = 0.5
    for ii in T0p5_att
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_green, linewidth = 5, label = L"\omega_T = 0.5")
    end

    # T = 0
    for ii in T0_att
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_blue, linewidth = 5, label = L"\omega_T = 0.0")
    end

    # Non-thermal motion
    for ii in [no_motion[2]]
        y_vals = ii.σs |> vec
        lines!(ax1, ii.τs, (y_vals .- y_vals[1]), color = my_black, linewidth = 7, label = "Non-thermal")
    end
end

CairoMakie.xlims!(0, 100)
CairoMakie.ylims!(0, 2000)
axislegend(ax1, labelsize = 40, unique = true, position = :lt)
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
