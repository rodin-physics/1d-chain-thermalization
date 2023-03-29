include("../../src/main.jl")

# Plotting speeds 
# fig = Figure(resolution = (2400, 600), fontsize = 36, figure_padding = 40)
# ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}", ylabel = "Probability Density", title = L"Temperature $\omega_T = 0.0$", xticks = -40:20:90)
# ax2 = Axis(fig[1, 2], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"\omega_T = 5.0", xticks = -40:20:90)
# ax3 = Axis(fig[1, 3], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"\omega_T = 25.0", xticks = -40:20:90)
# ax4 = Axis(fig[1, 4], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"\omega_T = 50.0", xticks = -40:20:90)

# ax_pair = [(0.0,ax1, my_blue), (5.0, ax2, my_green), (25.0, ax3, my_vermillion), (50.0, ax4, my_yellow)]

# for pair in ax_pair
#     data = readdlm("Final_speeds_α40_Φ2.0_λ4.0_bias0.01_T$(pair[1]).dat") |> vec
#     vlines!(pair[2], [0.0], color = my_black, linewidth = 3)
#     hist!(pair[2], data, color = pair[3], bins = 30, normalization = :pdf)
#     vlines!(pair[2], [40/n for n in 1:3], linestyle = :dash, color = my_black, linewidth = 3)
#     vlines!(pair[2], [√(8*π^2*2)], linewidth = 4, linestyle = :dash, color = my_red)
#     xlims!(pair[2], -40, 90)
#     ylims!(pair[2], 0, 0.055)
# end

# Label(fig[1,1, TopLeft()], "(a)", font = :bold)
# Label(fig[1,2, TopLeft()], "(b)", font = :bold)
# Label(fig[1,3, TopLeft()], "(c)", font = :bold)
# Label(fig[1,4, TopLeft()], "(d)", font = :bold)


# fig

    
## Stacked Histogram
timings = [L"0\tau", " ", L"250\tau", " ", L"500\tau", " ", L"750\tau", " ", L"1000\tau", " "]

# yticks = (0:0.04:0.36, reverse(timings))

fig = Figure(resolution = (1800, 1000), fontsize = 34, figure_padding = 40)
ax1 = Axis(fig[1, 1], title = L"Temperature $\omega_T$ = 0.0", xlabel = "Speed", ylabel = "Probability density", yticks = 0:0.08:0.4, xgridvisible = false)
ax2 = Axis(fig[1, 2], title = L"\omega_T = 5.0", xlabel = "Speed", yticklabelsvisible = false, ygridvisible = false, yticksvisible = false, xgridvisible = false)
ax3 = Axis(fig[1, 3], title = L"\omega_T = 25.0", xlabel = "Speed", yticklabelsvisible = false, ygridvisible = false, yticksvisible = false, xgridvisible = false)
ax4 = Axis(fig[1, 4], title = L"\omega_T = 50.0", xlabel = "Speed", ygridvisible = false, yaxisposition = :right, yticks = (0:0.04:0.36, reverse(timings)), yticksize = 20, yticksvisible = false, xgridvisible = false)

ax_pair = [(0.0, ax1, my_blue), (5.0, ax2, my_sky), (25.0, ax3, my_green), (50.0, ax4, my_vermillion)]

for pair in ax_pair 
    data = readdlm("Speeds_1000τ_inter250_α40_Φ2.0_λ4.0_bias0.01_T$(pair[1]).dat")
    num_times = size(data, 1)

    hlines!(pair[2], 0.08:0.08:0.36, linewidth = 2, color = (my_black, 0.2))

    for t in 1:num_times
        hist!(pair[2], data[t, :], bins = 44, offset = 0.08 * (num_times - t), color = pair[3], normalization = :pdf, strokewidth = 0.3, strokecolor = pair[3])
        # density!(pair[2], data[t, :],offset = 0.08 * (num_times - t), color = pair[3])
    end

    vlines!(pair[2], [40/n for n in 1:3], linestyle = :dash, color = my_black, linewidth = 2.5)
    vlines!(pair[2], [√(8*π^2*2)], linewidth = 3, linestyle = :dash, color = my_red, label = "Capture Speed")
    xlims!(pair[2], -20, 80)
    ylims!(pair[2], 0, 0.4)
end

# axislegend(ax4, position = :lt, labelsize = 20)
colgap!(fig.layout, 60)
fig