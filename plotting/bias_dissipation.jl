include("../src/main.jl")

## Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 8.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")
bias = 0.2
# bias = 8.1

## Resonances
resonance_speed(n) = (2 * ωmax * α) / ((2 * n) + 1)

## Computation
xs = range(1,65, step = 0.01)
xs2 = range(sqrt(8*π^2*Φ0/μ)+0.1,65, step = 0.01)
# Δ_bias = Δ_transport.(xs, Φ0, λ, ωmax, α)
# Δ_bias_rep = @showprogress map(x -> Δ_transport_smeared(x, Φ0, λ, ωmax, α), xs2)
# Δ_bias_att = @showprogress map(x -> Δ_transport_smeared(x, -Φ0, λ, ωmax, α), xs)

## Plotting Delta
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=34, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta_\mathrm{drift}", title = "Bias Dissipation - Bias = $bias")

# scatter!(ax1, xs, Δ_bias, color = my_black, markersize = 5)

scatter!(ax1, xs2, Δ_bias_rep, color = my_vermillion, markersize = 5)
scatter!(ax1, xs, Δ_bias_att, color = my_blue, markersize = 5)


vlines!(ax1, [α/n for n in 1:10], linewidth = 2, linestyle = :dash)
hlines!(ax1, [bias], color = my_blue, linewidth = 2, label = "Bias")
# vlines!(ax1, [14.168, 21.085, 40.114], linewidth = 2, linestyle = :dash)
vlines!(ax1, [sqrt(8*Φ0*π^2)], linestyle = :dash, color = my_green, linewidth = 2, label = "Min. Speed")

xlims!(ax1, 0, 65)
ylims!(ax1, -0.1, 50)
axislegend(ax1, labelsize = 26, position = :rt)
fig
