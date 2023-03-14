include("../../src/main.jl")

# Plotting speeds 
fig = Figure(resolution = (2400, 600), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}", ylabel = "Probability Density", title = L"Temperature $\omega_T = 0.0$", xticks = -40:20:90)
ax2 = Axis(fig[1, 2], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"\omega_T = 5.0", xticks = -40:20:90)
ax3 = Axis(fig[1, 3], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"\omega_T = 25.0", xticks = -40:20:90)
ax4 = Axis(fig[1, 4], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"\omega_T = 50.0", xticks = -40:20:90)

ax_pair = [(0.0,ax1, my_blue), (5.0, ax2, my_green), (25.0, ax3, my_vermillion), (50.0, ax4, my_yellow)]

for pair in ax_pair
    data = readdlm("Final_speeds_α40_Φ2.0_λ4.0_bias0.01_T$(pair[1]).dat") |> vec
    vlines!(pair[2], [0.0], color = my_black, linewidth = 3)
    density!(pair[2], data, color = pair[3])
    vlines!(pair[2], [40/n for n in 1:3], linestyle = :dash, color = my_black, linewidth = 3)
    vlines!(pair[2], [√(8*π^2*2)], linewidth = 4, linestyle = :dash, color = my_red)
    xlims!(pair[2], -40, 90)
    ylims!(pair[2], 0, 0.035)
end

Label(fig[1,1, TopLeft()], "(a)", font = :bold)
Label(fig[1,2, TopLeft()], "(b)", font = :bold)
Label(fig[1,3, TopLeft()], "(c)", font = :bold)
Label(fig[1,4, TopLeft()], "(d)", font = :bold)


fig