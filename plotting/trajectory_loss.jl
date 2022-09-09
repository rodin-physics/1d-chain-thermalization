include("../src/main.jl")
## Load data
data_rep = load_object("data/Non_Thermal/Single_sigma0[220]_sigmadot0[120]_MemInf_lambda4_Phi20_mu1_d60_bias0.0_omegaTnothing_tau125.jld2")

data_att = load_object("data/Non_Thermal/Single_sigma0[220]_sigmadot0[120]_MemInf_lambda4_Phi-20_mu1_d60_bias0.0_omegaTnothing_tau125.jld2")

# data_rep2 = load_object("data/Non_Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_bias0.2_ΩTnothing_τ200.jld2")
#
# data_att2 = load_object("data/Non_Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_bias0.2_ΩTnothing_τ200.jld2")

resonance_speed(n) = (2 * data_rep.ωmax * data_rep.α) / ((2 * n) + 1)

## Computation
(rep_xs, rep_ys) = Δ_traj(data_rep)
(att_xs, att_ys) = Δ_traj(data_att)

# (rep_xs2, rep_ys2) = Δ_traj(data_rep2)
# (att_xs2, att_ys2) = Δ_traj(data_att2)

σ_dots = range(7, 130, length=100)
Δs = Δ_analytic.(σ_dots, data_rep.Φ, data_rep.λ, data_rep.ωmax)

## Plotting
fig = Figure(resolution=(1600, 800), font="CMU Serif", fontsize=30, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta", xticks = 0:20:120)
inset_ax = Axis(fig[1, 1],
    width=Relative(0.45),
    height=Relative(0.45),
    halign=0.964,
    valign=0.9,
    backgroundcolor = :white,
    yticksize = 0,
    xticklabelsize = 20,
    yticklabelsize = 20,
    xlabel=L"\dot{\sigma}_0", ylabel=L"\Delta",
    xlabelsize = 22, ylabelsize = 22,
    xticks = 80:20:120)
# band!(ax1, 0:20, repeat([-1.0], 21), repeat([4.0], 21), color = RGBA(0.0, 0.0,0.0, 0.1))
vlines!(ax1, resonance_speed.(1:50), color = my_black, linewidth = 2)

vlines!(ax1, [sqrt(data_rep.Φ * π^2 * 8)], color = my_red, linewidth = 2, label = "Capture speed")
scatter!(ax1, rep_xs, rep_ys .+ data_rep.bias, markersize=12, color = my_vermillion, label = "Repulsive")
scatter!(ax1, att_xs, att_ys .+ data_rep.bias, markersize=12, color = my_blue, label = "Attractive")

# scatter!(ax1, rep_xs2, rep_ys2 .+ data_rep2.bias, markersize=12, color = my_vermillion, label = "Repulsive")
# scatter!(ax1, att_xs2, att_ys2.+ data_att2.bias, markersize=12, color = my_blue, label = "Attractive")

lines!(ax1, σ_dots, Δs, linewidth=4, color = my_black, label = "Analytic")


# text!(ax1, L"V_{\mathrm{bias}} = 0.05", position = (1.0, 1.5), textsize = 28)
# text!(ax1, L"V_{\mathrm{bias}} = 0.2", position = (52.0, -0.04), textsize = 28)

translate!(inset_ax.scene, 0, 0, 10)
# this needs separate translation as well, since it's drawn in the parent scene
translate!(inset_ax.elements[:background], 0, 0, 9)
vlines!(inset_ax, resonance_speed.(1:6), color = my_black, linewidth = 2)
lines!(inset_ax, σ_dots, Δs, linewidth=4, color = my_black)
lines!(inset_ax, att_xs, att_ys .+ data_rep.bias, color = my_blue, linewidth = 2)
lines!(inset_ax, rep_xs, rep_ys .+ data_rep.bias, color = my_vermillion, linewidth = 2)



xlims!(ax1, 0, 125)
ylims!(ax1, -0.05, 3.5)

xlims!(inset_ax, 80, 122)
ylims!(inset_ax, 0.485, 1.0)
axislegend(ax1, labelsize = 28, unique = true, position = :lt)
fig
