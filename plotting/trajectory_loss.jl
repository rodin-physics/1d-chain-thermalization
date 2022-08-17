include("../src/main.jl")
## Load data
data_rep = load_object("data/Non_Thermal/Single_σ0[220]_σdot0[50]_MemInf_λ4_Φ8_μ1_d60_bias0.0_ΩTnothing_τ100.jld2")

data_att = load_object("data/Non_Thermal/Single_σ0[220]_σdot0[50]_MemInf_λ4_Φ-8_μ1_d60_bias0.0_ΩTnothing_τ100.jld2")

# data_rep2 = load_object("data/Non_Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ1_μ1_d60_bias0.2_ΩTnothing_τ200.jld2")
#
# data_att2 = load_object("data/Non_Thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ-1_μ1_d60_bias0.2_ΩTnothing_τ200.jld2")

resonance_speed(n) = (2 * data_rep.ωmax * data_rep.α) / ((2 * n) + 1)

## Computation
(rep_xs, rep_ys) = Δ_traj(data_rep)
(att_xs, att_ys) = Δ_traj(data_att)

# (rep_xs2, rep_ys2) = Δ_traj(data_rep2)
# (att_xs2, att_ys2) = Δ_traj(data_att2)

σ_dots = range(10, 50, length=100)
Δs = Δ_analytic.(σ_dots, data_rep.Φ, data_rep.λ, data_rep.ωmax)

## Plotting
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=34, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta + V_{\mathrm{bias}}")
# band!(ax1, 0:20, repeat([-1.0], 21), repeat([1.0], 21), color = RGBA(0.0, 0.0,0.0, 0.1))
scatter!(ax1, rep_xs, rep_ys .+ data_rep.bias, markersize=12, color = my_vermillion, label = "Repulsive")
scatter!(ax1, att_xs, att_ys .+ data_rep.bias, markersize=12, color = my_blue, label = "Attractive")

# scatter!(ax1, rep_xs2, rep_ys2 .+ data_rep2.bias, markersize=12, color = my_vermillion, label = "Repulsive")
# scatter!(ax1, att_xs2, att_ys2.+ data_att2.bias, markersize=12, color = my_blue, label = "Attractive")

lines!(ax1, σ_dots, Δs, linewidth=4, color = my_black, label = "Analytic")

vlines!(ax1, resonance_speed.(1:20), color = my_black, linewidth = 2, linestyle = :dash)
# text!(ax1, L"V_{\mathrm{bias}} = 0.05", position = (1.0, -0.04), textsize = 28)
# text!(ax1, L"V_{\mathrm{bias}} = 0.2", position = (52.0, -0.04), textsize = 28)

xlims!(ax1, 10, 60)
ylims!(ax1, -0.05, 0.65)
axislegend(ax1, labelsize = 24, unique = true)
fig
