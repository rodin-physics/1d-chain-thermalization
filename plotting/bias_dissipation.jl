include("../src/main.jl")

## Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 4.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
# bias = 0.21465496799193154
bias = 8.1

## Resonances
resonance_speed(n) = (2 * ωmax * α) / ((2 * n) + 1)

## Computation
xs = range(1,80, step = 0.001)
Δ_bias = Δ_transport.(xs, Φ0, λ, ωmax, α)

## Plotting Delta
fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=34, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta_\mathrm{drift}", title = "Bias Dissipation")

lines!(ax1, xs, Δ_bias, color = my_black, linewidth = 2)
# vlines!(ax1, [α/n for n in 1:10], linewidth = 2, linestyle = :dash)
hlines!(ax1, [bias], color = my_blue, linewidth = 2)
# vlines!(ax1, [11.5, 17.2, 28.1, 44.4], linewidth = 2, linestyle = :dash)

xlims!(ax1, 0, 80)
ylims!(ax1, -0.1, 10)
fig
