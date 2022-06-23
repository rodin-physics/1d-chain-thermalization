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

# ROOT MEAN SQUARE DEVIATION
function root_mse(all_data, non_thermal)
    σs = [vec(data.σs .- data.σs[1]) for data in all_data]
    σs_no_T = vec(non_thermal.σs .- non_thermal.σs[1])

    error = transpose(reduce(hcat, σs .- repeat([σs_no_T], length(σs))))

    return (non_thermal.τs, vec(sqrt.(mean(error.^2, dims = 1))))
end

## Plotting
fig = Figure(resolution = (1400, 1600), font = "CMU Serif", fontsize = 36, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = "Root MSE")
ax2 = Axis(fig[2,1], xlabel = L"\tau", ylabel = "Root MSE")

all_rep = [T0_rep, T0p5_rep, T1_rep, T10_rep]
all_att = [T0_att, T0p5_att, T1_att, T10_att]
colors = [my_blue, my_green, my_red, my_vermillion]

for ii in reverse(1:length(all_rep))
    data = all_rep[ii]
    res = root_mse(data, no_motion[1])
    lines!(ax1, res[1], res[2], color = colors[ii] , linewidth = 5, label = L"\omega_T = %$(data[1].ωT)")
end

for ii in reverse(1:length(all_att))
    data = all_att[ii]
    res = root_mse(data, no_motion[2])
    lines!(ax2, res[1], res[2], color = colors[ii] , linewidth = 5, label = L"\omega_T = %$(data[1].ωT)")
end

axislegend(ax1, labelsize = 25, unique = true, position = :lt)
axislegend(ax2, labelsize = 25, unique = true, position = :lt)
xlims!(ax1, 0, 100)
xlims!(ax2, 0, 100)
fig
