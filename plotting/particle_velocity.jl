include("../src/main.jl")


data_rep = load_object(
    "data/non_thermal/Single_σ0[200]_σdot0[50]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩTnothing_τ100.jld2")

data_att = load_object("data/Non_Thermal/Single_σ0[220]_σdot0[50]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩTnothing_τ100.jld2")

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

σs_rep = data_rep.σs |> vec
τs_rep = data_rep.τs
δ_rep = τs_rep[2] - τs_rep[1]

σs_att = data_att.σs |> vec
τs_att = data_att.τs
δ_att = τs_att[2] - τs_att[1]

speeds_rep = [((σs_rep[ii+1] - σs_rep[ii]) / δ_rep) for ii in 1:(length(σs_rep)-1)]
speeds_att = [((σs_att[ii+1] - σs_att[ii]) / δ_att) for ii in 1:(length(σs_att)-1)]

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 34, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}")
# ax2 = Axis(fig[2,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}", title = "Attractive")

hlines!(ax1, [speeds_rep[1]], linewidth = 3, color = my_black, linestyle = :dash, label = "Initial Speed")
lines!(ax1, τs_rep[2:end], speeds_rep, linewidth = 3, color = my_vermillion, label = "Repulsive")

hlines!(ax1, [speeds_att[1]], linewidth = 3, color = my_black, linestyle = :dash)
lines!(ax1, τs_att[2:end], speeds_att, linewidth = 3, color = my_blue, label = "Attractive")
# hlines!(ax2, resonance_speed.(3:12))

xlims!(ax1, 0, 1)
ylims!(ax1, 0, nothing)

# xlims!(ax2, 0, 100)
# ylims!(ax2, 0, 50)
axislegend(ax1, labelsize = 25, position = :lb)
fig
