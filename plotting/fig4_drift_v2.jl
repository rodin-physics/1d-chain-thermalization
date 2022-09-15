include("../src/main.jl")
using Peaks, DelimitedFiles

# Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 2.0
λ = 4.0
bias = 0.01

## Bias drift calculation
xs = range(1,90, step = 0.005)
# Δ_bias = Δ_transport.(xs, Φ0, λ, ωmax, α)


## Plotting
fig = Figure(resolution = (1800, 1200), font = "CMU Serif", fontsize = 32, figure_padding = 30)
ax1 = Axis(fig[2,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}", xticks = 0:200:800)
ax2 = Axis(fig[2,2], xlabel = L"\tau", xticks = 200:200:800, yticklabelsvisible = false)
ax3 = Axis(fig[2,3], xlabel=L"\Delta_\mathrm{drift}", ylabel=L"\dot{\sigma}")

# Top axes
ax4 = Axis(fig[1,1], xlabel = L"\dot{\sigma}", ylabel = L"Speed Diff after 200 $\tau$", xticks = 50:10:100, title = "Repulsive")
ax5 = Axis(fig[1,2], xlabel = L"\dot{\sigma}", xticks = 50:10:100, title = "Attractive", yticklabelsvisible = false)
ax6 = Axis(fig[1,3], xlabel=L"\Delta_\mathrm{drift}", ylabel=L"\dot{\sigma}")
rowsize!(fig.layout, 1, Relative(0.25))

hlines!(ax1, [α/n for n in 1:8], linewidth = 2, linestyle = :dash, color = my_black)
hlines!(ax1, [sqrt(8*Φ0*π^2)], linewidth = 2, color = my_red, linestyle = :dash, label = "Capture speed")

hlines!(ax2, [α/n for n in 1:8], linewidth = 2, linestyle = :dash, color = my_black)

hlines!(ax3, [α/n for n in 1:8], linewidth = 2, linestyle = :dash, color = my_black)
hlines!(ax4, [0.0], color = my_black, linewidth = 3)
hlines!(ax5, [0.0], color = my_black, linewidth = 3)

scatter!(ax3, Δ_bias, xs, color = my_black, linewidth = 2)
scatter!(ax6, Δ_bias, xs, color = my_black, linewidth = 2)


# Read in full trajectories
filenames = filter(x -> first(x) !== '.', readdir(joinpath(pwd(), "data/proc_tau800/")))

for ii in 1:length(filenames)
    data = readdlm(joinpath("data/proc_tau800/", filenames[ii]))
    τs = data[:,1]
    speeds = data[:,2]

    if occursin("Phi2", filenames[ii])
        lines!(ax1, τs, speeds, linewidth = 4, color = my_vermillion)

    elseif occursin("Phi-2", filenames[ii])
        lines!(ax2, τs, speeds, linewidth = 4, color = my_blue)
    end

end

# Read in full trajectories
filenames = filter(x -> first(x) !== '.', readdir(joinpath(pwd(), "data/proc_tau200/")))

for ii in 1:length(filenames)
    data = readdlm(joinpath("data/proc_tau200/", filenames[ii]))
    τs = data[:,1]
    speeds = data[:,2]
    speed_diff = speeds[end] - speeds[1]

    if occursin("Phi2", filenames[ii])
        scatter!(ax4, [speeds[1]], [speed_diff], markersize = 13, color = my_vermillion)

    elseif occursin("Phi-2", filenames[ii])
        scatter!(ax5, [speeds[1]], [speed_diff], markersize = 13, color = my_blue)
    end

end



xlims!(ax1, 0, 800)
axislegend(ax1, labelsize = 30, position = :rb)
#
xlims!(ax2, 0, 800)
# # hideydecorations!(ax2)

vlines!(ax3, [0.01], label = "Bias", color = my_black, linewidth = 2)
xlims!(ax3, 0, 0.5)

xlims!(ax4, 55, 91)
xlims!(ax5, 55, 91)
xlims!(ax6, 0, 0.1)

ylims!(ax1, 0, 94)
ylims!(ax2, 0, 94)
ylims!(ax3, 0, 50)
# ylims!(ax4, 58, 80)
# ylims!(ax5, 58, 80)

vlines!(ax6, [0.01], label = "Bias", color = my_black, linewidth = 2)
ylims!(ax6, 58, 90)
# axislegend(ax3, labelsize = 30, position = :rb)
fig
