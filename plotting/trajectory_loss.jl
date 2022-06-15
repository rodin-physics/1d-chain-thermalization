## Load data
data_rep = load_object("data/non_thermal/Single_σ0[55]_σdot0[25]_MemInf_λ1_Φ1_μ1_d60_ΩTnothing_τ200.jld2")

data_att = load_object("data/non_thermal/Single_σ0[50]_σdot0[25]_MemInf_λ1_Φ-1_μ1_d60_ΩTnothing_τ200.jld2")

## Computation
(rep_xs, rep_ys) = Δ_traj(data_rep)
(att_xs, att_ys) = Δ_traj(data_att)

σ_dots = range(2, 25, length=100)
Δs = Δ_analytic.(σ_dots, data.Φ, data.λ, data.ωmax)

## Plotting
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=36)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta")

scatter!(ax1, rep_xs, rep_ys, marker='⊕', markersize=20)
scatter!(ax1, att_xs, att_ys, marker='⊖', markersize=20)
lines!(ax1, σ_dots, Δs, linewidth=4)

fig
