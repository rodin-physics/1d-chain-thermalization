include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 2.0
λ = 4.0
ωTs = [0.0, 1.0, 5.0, 10.0]

colors = [my_blue, my_green, my_vermillion, my_yellow]


fig = Figure(resolution = (1200, 1600), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"[C_0(0) - C_0(\tau)]/C_0(0)")
inset_ax = Axis(fig[1, 1], width=Relative(0.35), height=Relative(0.35), halign=0.9, valign=0.22, xlabel = L"\omega_T", ylabel = L"C_0(0)")
ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\int dx \, [ \tilde{\Phi}'(x,\tau) / \Phi_0 ]^2")

τs = range(0, 10, length = 1000)
ωT_range = range(0, 10, length = 50)

## First axis
for ωT in ωTs 
    broadening_xs = map(x -> (C_corr(0, 0, ωmax, ωT) - C_corr(x, 0, ωmax, ωT)) / C_corr(0, 0, ωmax, ωT), τs)

    lines!(ax1, τs, broadening_xs, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
end

## Inset axis 
for ωT in ωTs
    scatter!(inset_ax, [ωT], [C_corr(0,0,ωmax, ωT)], color = my_blue)
end
factors = map(x -> C_corr(0,0,ωmax, x), ωT_range)
lines!(inset_ax, ωT_range, factors, color = my_blue, linewidth = 3)

## Second axis 
ωTs = [0.0, 5.0, 10.0, 25.0, 100.0, 250.0]
colors = [my_blue, my_green, my_sky, my_vermillion, my_red, my_yellow]

# Temperature induced broadened potential - normalised by Φ0 and squared
function broadened_potential(x, τ, ωT)
    factor = λ^2 + (C_corr(0, 0, ωmax, ωT) - C_corr(τ, 0, ωmax, ωT))
    return (λ  / factor) * exp(- x^2 / factor / 2)
end

for ωT in ωTs 
    res = [quadgk(x -> broadened_potential(x, τ, ωT), -20, 20, rtol = 1e-5)[1] for τ in τs]
    lines!(ax2, τs, res, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
end

# for ωT in ωTs 
#     res = (λ^2 * √(π) / 2) .* (λ^2 .+ C_corr(0, 0, ωmax, ωT) .- C_corr.(τs, 0, ωmax, ωT)).^(-5/2)
#     lines!(ax2, τs, res, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
# end


axislegend(ax1, position = :rt, nbanks = 2)
axislegend(ax2, position = :rt, nbanks = 3)

xlims!(ax1, 0, 10)
ylims!(ax1, 0, nothing)

ylims!(inset_ax, 0, nothing)

xlims!(ax2, 0, 10)
ylims!(ax2, 0, 0.017)

Label(fig[1,1, TopLeft()], "(a)", font = :bold)
Label(fig[2,1, TopLeft()], "(b)", font = :bold)

translate!(inset_ax.scene, 0, 0, 10)
# this needs separate translation as well, since it's drawn in the parent scene
translate!(inset_ax.elements[:background], 0, 0, 9)
fig



