include("../src/main.jl")

## Confinement-free chain
@inline function ω0(ωmax, x)
    return sqrt(ωmax^2 * sin(x)^2)
end

function Γ0(τ, ls, ωmax)
    int_fun(x) = cos.(2 * x * ls) * sin(2 * π * τ * ω0(ωmax, x)) / ω0(ωmax, x)
    res = quadgk(int_fun, 0, π / 2)
    return (res[1] * 2 / π)
end

ωmax = 10
l_max = 500
times = [1, 5, 10, 15]

Γ0s = [Γ0(τ, 0:l_max, ωmax) for τ in times]
Γs = [Γ(τ, 0:l_max, ωmax) for τ in times]
speed0 = ωmax / 2 * 2 * π
speed = (ωmax - 1) / 2 * 2 * π

fig = Figure(resolution = (600, 960), font = "CMU Serif", fontsize = 18)

ax1 = Axis(fig[1, 1], ylabel = L"\Gamma_n(1)")
ax2 = Axis(fig[2, 1], ylabel = L"\Gamma_n(5)")
ax3 = Axis(fig[3, 1], ylabel = L"\Gamma_n(10)")
ax4 = Axis(fig[4, 1], xlabel = L"n", ylabel = L"\Gamma_n(20)")
axes = [ax1, ax2, ax3, ax4]
txts = ["(a)", "(b)", "(c)", "(d)"]
fig
for ii = 1:length(times)
    ln = vlines!(
        axes[ii],
        times[ii] * [speed0, speed],
        linewidth = 2,
        linestyle = :dash,
        color = [my_red, my_blue],
    )
    sc = scatter!(axes[ii], 0:l_max, Γ0s[ii], markersize = 5, color = my_red)
    sc2 = scatter!(axes[ii], 0:l_max, Γs[ii], markersize = 5, color = my_blue)
    xlims!(axes[ii], (-0.05, 500))
    ylims!(axes[ii], (-0.05, 0.15))
    if ii != length(times)
        hidexdecorations!(axes[ii], grid = false)
    end
    text!(axes[ii], txts[ii], position = (475, 0.1), textsize = 18, font = "CMU Serif")
end

fig
save("Pulse_Propagation.pdf", fig)