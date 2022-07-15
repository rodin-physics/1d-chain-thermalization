include("../src/main.jl")
speed = 80
# bias = 0.21465496799193154
bias = 0.0

data_rep = [load_object("Data/Non_Thermal_Multi/$(n)_Single_σ0[220]_σdot0[$(speed)]_MemInf_λ4_Φ8_μ1_d60_bias$(bias)_ΩTnothing_τ30.jld2") for n in 1:10]

data_att = [load_object("Data/Non_Thermal_Multi/$(n)_Single_σ0[220]_σdot0[$(speed)]_MemInf_λ4_Φ-8_μ1_d60_bias$(bias)_ΩTnothing_τ30.jld2") for n in 1:10]
resonance_speed(n) = (2 * data_rep[1].ωmax * data_rep[1].α) / ((2 * n) + 1)

## Computation
σ_dots = range(7, 70, length=100)
Δs = Δ_analytic.(σ_dots, data_rep[1].Φ, data_rep[1].λ, data_rep[1].ωmax)

## Plotting
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=34, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta + V_{\mathrm{bias}}")
# band!(ax1, 0:20, repeat([-1.0], 21), repeat([1.0], 21), color = RGBA(0.0, 0.0,0.0, 0.1))

for ii in 1:length(data_rep)
    (rep_xs, rep_ys) = Δ_traj2(data_rep[ii])
    scatter!(ax1, rep_xs, rep_ys .+ data_rep[ii].bias, markersize=12, color = my_vermillion, label = "Repulsive")
    # if !(isempty(rep_xs))
    #     vlines!(ax1, [rep_xs[1]], linestyle = :dashdot, color = my_vermillion, linewidth = 2)
    # end
end


for ii in 1:length(data_att)
    (att_xs, att_ys) = Δ_traj2(data_att[ii])
    scatter!(ax1, att_xs, att_ys .+ data_rep[ii].bias, markersize=12, color = my_blue, label = "Attractive")
    # if !(isempty(att_xs))
    #     vlines!(ax1, [att_xs[1]], linestyle = :dashdot, color = my_blue, linewidth = 2)
    # end
end

# scatter!(ax1, rep_xs2, rep_ys2 .+ data_rep2.bias, markersize=12, color = my_vermillion, label = "Repulsive")
# scatter!(ax1, att_xs2, att_ys2.+ data_att2.bias, markersize=12, color = my_blue, label = "Attractive")

lines!(ax1, σ_dots, Δs, linewidth=4, color = my_black, label = "Analytic")
vlines!(ax1, [40])
# vlines!(ax1, resonance_speed.(1:20), color = my_black, linewidth = 2, linestyle = :dash)
# text!(ax1, L"V_{\mathrm{bias}} = 0.05", position = (1.0, -0.04), textsize = 28)
# text!(ax1, L"V_{\mathrm{bias}} = 0.2", position = (52.0, -0.04), textsize = 28)

xlims!(ax1, 0, 90)
ylims!(ax1, -0.1, 0.65)
axislegend(ax1, labelsize = 24, unique = true)
fig
