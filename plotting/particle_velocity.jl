include("../src/main.jl")

data_rep = load_object("data/Thermal/Single_σ0[45]_σdot0[30]_MemInf_λ1_Φ1_μ1_d60_ΩT0.0_τ100.jld2")

data_att = load_object("data/Thermal/Single_σ0[45]_σdot0[30]_MemInf_λ1_Φ-1_μ1_d60_ΩT0.0_τ100.jld2")

σs_rep = data_rep.σs |> vec
τs_rep = data_rep.τs
δ_rep = τs_rep[2] - τs_rep[1]

σs_att = data_att.σs |> vec
τs_att = data_att.τs
δ_att = τs_att[2] - τs_att[1]

speeds_rep = [((σs_rep[ii+1] - σs_rep[ii]) / δ_rep) for ii in 1:(length(σs_rep)-1)]
speeds_att = [((σs_att[ii+1] - σs_att[ii]) / δ_att) for ii in 1:(length(σs_att)-1)]

fig = Figure(resolution = (1800, 1600), font = "CMU Serif", fontsize = 36)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\sigma", title = "Repulsive")
ax2 = Axis(fig[2,1], xlabel = L"\tau", ylabel = L"\sigma", title = "Attractive")

lines!(ax1, τs_rep[2:end], speeds_rep, linewidth = 3, color = my_blue)
lines!(ax2, τs_att[2:end], speeds_att, linewidth = 3, color = my_blue)
fig
