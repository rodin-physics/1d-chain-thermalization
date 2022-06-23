include("../src/main.jl")

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


# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end


## Plotting
fig = Figure(resolution = (1800, 1600), font = "CMU Serif", fontsize = 36, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\sigma", title = "Repulsive")
ax2 = Axis(fig[2,1], xlabel = L"\tau", ylabel = L"\sigma", title = "Attractive")

# REPULSIVE
for data in T10_rep
    res = particle_speed(data)
    lines!(ax1, res[1], res[2], linewidth = 3, color = my_vermillion, label = L"\omega_T = 10.0")
end

for data in T1_rep
    res = particle_speed(data)
    lines!(ax1, res[1], res[2], linewidth = 3, color = my_red, label = L"\omega_T = 1.0")
end

for data in T0p5_rep
    res = particle_speed(data)
    lines!(ax1, res[1], res[2], linewidth = 3, color = my_green, label = L"\omega_T = 0.5")
end

for data in T0_rep
    res = particle_speed(data)
    lines!(ax1, res[1], res[2], linewidth = 3, color = my_blue, label = L"\omega_T = 0.0")
end

for data in [no_motion[1]]
    res = particle_speed(data)
    lines!(ax1, res[1], res[2], linewidth = 3, color = my_black, label = "Non-Thermal")
end


# ATTRACTIVE
for data in T10_att
    res = particle_speed(data)
    lines!(ax2, res[1], res[2], linewidth = 3, color = my_vermillion, label = L"\omega_T = 10.0")
end

for data in T1_att
    res = particle_speed(data)
    lines!(ax2, res[1], res[2], linewidth = 3, color = my_red, label = L"\omega_T = 1.0")
end

for data in T0p5_att
    res = particle_speed(data)
    lines!(ax2, res[1], res[2], linewidth = 3, color = my_green, label = L"\omega_T = 0.5")
end

for data in T0_att
    res = particle_speed(data)
    lines!(ax2, res[1], res[2], linewidth = 3, color = my_blue, label = L"\omega_T = 0.0")
end
for data in [no_motion[2]]
    res = particle_speed(data)
    lines!(ax2, res[1], res[2], linewidth = 3, color = my_black, label = "Non-Thermal")
end

axislegend(ax1, labelsize = 24, position = :lb, unique = true)
axislegend(ax2, labelsize = 24, position = :lb, unique = true)
xlims!(ax1, 0.0, 30)
xlims!(ax2, 0.0, 45)
fig
