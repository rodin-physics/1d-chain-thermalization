include("../src/main.jl")
using Peaks

# Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 8.0
λ = 4.0
bias = 0.2
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")

## Bias drift calculation
xs = range(1,65, step = 0.001)
# Δ_bias = Δ_transport.(xs, Φ0, λ, ωmax, α)


## Plotting
fig = Figure(resolution = (1800, 1000), font = "CMU Serif", fontsize = 36, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}", xticks = 0:60:750, title = "Repulsive")
ax2 = Axis(fig[1,2], xlabel = L"\tau", xticks = 0:60:750, title = "Attractive", yticklabelsvisible = false)
ax3 = Axis(fig[1,3], xlabel=L"\Delta_\mathrm{drift}", ylabel=L"\dot{\sigma}")


# Read in full trajectories
# filenames = filter(x -> first(x) !== '.', readdir(joinpath(pwd(), "data/Non_Thermal/")))

# data = [load_object(joinpath(pwd(), "data/Non_Thermal/", x)) for x in filenames]
# data_rep = filter(x -> x.Φ == Φ0 && x.α == α && x.bias == bias, data)
# data_att = filter(x -> x.Φ == -Φ0 && x.α == α && x.bias == bias, data)

# hlines!(ax1, [α/n for n in 1:8], linewidth = 2, linestyle = :dash, color = my_black)
# hlines!(ax1, [sqrt(8*Φ0*π^2)], linewidth = 2, color = my_red, linestyle = :dash, label = "Minimum speed")
# for ii in 1:length(data_rep)
#     (τs, speeds) = particle_speed(data_rep[ii])
#     # lines!(ax1, τs, speeds, linewidth = 3, color = my_vermillion)
#     if all(speeds.>=0)
#         pks, vals = findmaxima(speeds)
#         lines!(ax1, τs[pks], vals, linewidth = 4, color = my_vermillion)
#     end
#
# end

# hlines!(ax2, [α/n for n in 1:8], linewidth = 2, linestyle = :dash, color = my_black)
# for ii in 1:length(data_att)
#     (τs, speeds) = particle_speed(data_att[ii])
#     # lines!(ax2, τs, speeds, linewidth = 3, color = my_blue)
#
#     if all(speeds.>=0)
#         pks, vals = findminima(speeds)
#         lines!(ax2, τs[pks], vals, linewidth = 4, color = my_blue)
#     end
#
# end

hlines!(ax3, [α/n for n in 1:8], linewidth = 2, linestyle = :dash, color = my_black)
scatter!(ax3, Δ_bias, xs, color = my_black, linewidth = 2)


# xlims!(ax1, 0, 200)
# axislegend(ax1, labelsize = 30, position = :rb)
#
# xlims!(ax2, 0, 200)
# # hideydecorations!(ax2)

vlines!(ax3, [0.2])
xlims!(ax3, 0, 1)
ylims!(ax1, 0, 60)
ylims!(ax2, 0, 60)
ylims!(ax3, 0, 60)

fig
