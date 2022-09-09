include("../src/main.jl")
using Peaks
α= 40
Φ0=  2.0
bias = 0.01

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end


## Plotting

fig = Figure(resolution = (1800, 1200), font = "CMU Serif", fontsize = 32, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}", xticks = 0:100:800, title = "Repulsive")
ax2 = Axis(fig[1,2], xlabel = L"\tau", ylabel = L"\dot{\sigma}", xticks = 0:100:800, title = "Attractive")

# Vertical line solutions
hlines!(ax1, [α/n for n in 1:10], linewidth = 2.5, linestyle = :dash, color = my_black, label = "Vertical solutions")
hlines!(ax2, [α/n for n in 1:10], linewidth = 2.5, linestyle = :dash, color = my_black, label = "Vertical solutions")

# Minimum speeds
hlines!(ax1, [sqrt(8*Φ0*π^2)], linewidth = 2.5, color = my_red, linestyle = :dash, label = "Minimum speed")
hlines!(ax2, [1.9], linewidth = 2.5, color = my_red, linestyle = :dash, label = "Minimum speed")

# Read in full trajectories
filenames = filter(x -> first(x) !== '.', readdir(joinpath(pwd(), "data/final_tau800/")))
#
# data = [load_object(joinpath(pwd(), "data/final_tau800/", x)) for x in filenames]
# data_rep = filter(x -> x.Φ == Φ0 && x.α == α && x.bias == bias, data)
# data_att = filter(x -> x.Φ == -Φ0 && x.α == α && x.bias == bias, data)
#

for ii in 1:length(filenames)
    data = load_object(joinpath(pwd(), "data/final_tau800/", filenames[ii]))
    (τs, speeds) = particle_speed(data)

    if data.Φ > 0 && all(speeds.>=0)
        pks, vals = findmaxima(speeds)
        lines!(ax1, τs[pks], vals, linewidth = 4, color = my_vermillion)
    elseif data.Φ < 0 && all(speeds.>=0)
        pks, vals = findminima(speeds)
        lines!(ax2, τs[pks], vals, linewidth = 4, color = my_blue)
    end
end


xlims!(ax1, 0.0, 800)
xlims!(ax2, 0.0, 800)

ylims!(ax1, 0.0, nothing)
ylims!(ax2, 0, nothing)

axislegend(ax1, labelsize = 30, position = :rb, unique = true)
axislegend(ax2, labelsize = 30, position = :rb, unique = true)
save("final_Bias$(bias)_Phi$(Φ0)_alpha$(α).pdf", fig)
