include("../src/main.jl")
using Peaks

data_rep = load_object(
    "data/non_thermal/Bias_Single_σ0[180]_σdot0[13]_MemInf_λ4_Φ2_μ1.0_d60_ΩTnothing_τ800.jld2")

data_att = load_object("data/Non_Thermal/Single_sigma0[220]_sigmadot0[42]_MemInf_lambda4_Phi-8_mu1_d60_bias0.2_omegaTnothing_tau200.jld2")

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end
(τs_rep, speeds_rep) = particle_speed(data_rep)
(τs_att, speeds_att) = particle_speed(data_att)


fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 34, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}")
# ax2 = Axis(fig[2,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}", title = "Attractive")

# hlines!(ax1, speeds_att[1], linewidth = 3, color = my_black, linestyle = :dash)
lines!(ax1, τs_att, speeds_att, linewidth = 3, color = my_blue, label = "Attractive")

# vlines!(ax1,  τs_att[findfirst(x -> isapprox(x,20.0, rtol = 1e-8) == false, speeds_att)])

# hlines!(ax1, speeds_rep[1], linewidth = 3, color = my_black, linestyle = :dash, label = "Initial Speed")
lines!(ax1, τs_rep, speeds_rep, linewidth = 3, color = my_vermillion, label = "Repulsive")
# pks, vals = findminima(speeds_att)
# lines!(ax1, τs_att[pks], vals, linewidth = 4, color = my_blue)
# pks, vals = findmaxima(speeds_rep)
# lines!(ax1, τs_rep[pks], vals, linewidth = 4, color = my_vermillion)
hlines!(ax1, [data_rep.α/n for n in 1:10], linewidth = 2, linestyle = :dash, color = :black)
hlines!(ax1, [37.8, 41.79], linewidth = 2, linestyle = :dash, color = :black)

# hlines!(ax1, resonance_speed.(1))

xlims!(ax1, 0, 800)
ylims!(ax1, 0, nothing)

# xlims!(ax2, 0, 30)
# ylims!(ax2, 0, 50)
axislegend(ax1, labelsize = 25, position = :lb)
fig
