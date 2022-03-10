include("../src/main.jl")

K = 1
k = 20
m = 1
l_max = 500
times = [1, 10, 25, 50, 100]

Gs0 = [G(t, 0:l_max, 0, k, m) for t in times]
Gs = [G(t, 0:l_max, K, k, m) for t in times]

speed =  √((2 * k + 0 - √(0 * (4 * k + 0))) / m) / √(2)   # speed = √(k / m) / 2
speed2 = √((2 * k + K - √(K * (4 * k + K))) / m) / √(2)

fig = Figure(resolution = (600, 1200), font = "CMU Serif", fontsize = 18)

ax1 = Axis(fig[1, 1], ylabel = L"\Delta r")
ax2 = Axis(fig[2, 1], ylabel = L"\Delta r")
ax3 = Axis(fig[3, 1], ylabel = L"\Delta r")
ax4 = Axis(fig[4, 1], ylabel = L"\Delta r")
ax5 = Axis(fig[5, 1], xlabel = L"n", ylabel = L"\Delta r")
axes = [ax1, ax2, ax3, ax4, ax5]
txts = ["(a)", "(b)", "(c)", "(d)", "(e)"]
for ii = 1:length(times)
    ln = vlines!(
        axes[ii],
        times[ii] * [speed, speed2],
        linewidth = 2,
        linestyle = :dash,
        color = [my_red, my_blue],
    )
    sc = scatter!(axes[ii], 0:l_max, Gs0[ii], markersize = 5, color = my_red)
    sc2 = scatter!(axes[ii], 0:l_max, Gs[ii], markersize = 5, color = my_blue)
    ylims!(axes[ii], (-0.05, 0.15))
    if ii != length(times)
        hidexdecorations!(axes[ii], grid = false)
    end
    text!(axes[ii], txts[ii], position = (475, 0.1), textsize = 18, font = "CMU Serif")
end

fig
save("Pulse_Propagation.pdf", fig)