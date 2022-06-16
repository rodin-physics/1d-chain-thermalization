include("../src/main.jl")
## Load data
data_rep = load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ1_Φ1_μ1_d60_ΩT0.0_τ100.jld2")

data_att = load_object("data/Thermal/Single_σ0[45]_σdot0[20]_MemInf_λ1_Φ-1_μ1_d60_ΩT0.0_τ100.jld2")

## Computation
(rep_xs, rep_ys) = Δ_traj(data_rep)
(att_xs, att_ys) = Δ_traj(data_att)

σ_dots = range(2, 25, length=100)
Δs = Δ_analytic.(σ_dots, data.Φ, data.λ, data.ωmax)

## Plotting
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=36)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta")

scatter!(ax1, rep_xs, rep_ys, marker='⊕', markersize=20, color = my_blue)
scatter!(ax1, att_xs, att_ys, marker='⊖', markersize=20, color = my_red)

lines!(ax1, σ_dots, Δs, linewidth=4, color = my_black)

fig
