include("../src/main.jl")
using Peaks
# speed = 45
bias = 0.21465496799193154
# bias = 16.0

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

## Plotting

fig = Figure(resolution = (1800, 1200), font = "CMU Serif", fontsize = 34, figure_padding = 30)
ax1 = Axis(fig[1,1], xlabel = L"\tau", ylabel = L"\dot{\sigma}", xticks = 0:60:600, title = "Repulsive")
ax2 = Axis(fig[1,2], xlabel = L"\tau", ylabel = L"\dot{\sigma}", xticks = 0:60:600, title = "Attractive")

segment = 30 .* range(0, 20, step = 1)

final_τ = 0
for speed in [29, 35, 45, 50]
    data_rep = [load_object("Data/Non_Thermal_Multi/$(n)_Single_σ0[220]_σdot0[$(speed)]_MemInf_λ4_Φ8_μ1_d60_bias$(bias)_ΩTnothing_τ30.jld2") for n in 1:15]

    data_att = [load_object("Data/Non_Thermal_Multi/$(n)_Single_σ0[220]_σdot0[$(speed)]_MemInf_λ4_Φ-8_μ1_d60_bias$(bias)_ΩTnothing_τ30.jld2") for n in 1:15]

    final_τ = 0
    for ii in 1:length(data_rep)
        (τs, speeds) = particle_speed(data_rep[ii])
        # lines!(ax1, τs .+ final_τ, speeds, linewidth = 3, color = my_vermillion)

        pks, vals = findmaxima(speeds)
        lines!(ax1, τs[vcat([1], pks)] .+ final_τ, vcat(speeds[1], vals), linewidth = 4, color = my_vermillion)
        final_τ += τs[pks[end]]
    end

    final_τ = 0
    for ii in 1:length(data_att)
        (τs, speeds) = particle_speed(data_att[ii])
        # lines!(ax2, τs .+ final_τ, speeds, linewidth = 3, color = my_blue)

        pks, vals = findminima(speeds)
        lines!(ax2, τs[vcat([1], pks)] .+ final_τ, vcat(speeds[1], vals), linewidth = 4, color = my_blue)

        final_τ += τs[pks[end]]
    end

end

# hlines!(ax1, [11.5, 17.2, 28.1, 44.4], linewidth = 3, color = my_black, linestyle = :dash)
# hlines!(ax1, [α/n for n in 1:5], linewidth = 3, linestyle = :dash, color = my_blue)
#
# hlines!(ax2, [11.5, 17.2, 28.1, 44.4], linewidth = 3, color = my_black, linestyle = :dash)
# hlines!(ax2, [α/n for n in 1:5], linewidth = 3, linestyle = :dash, color = my_blue)

xlims!(ax1, 0.0, 450)
xlims!(ax2, 0.0, 450)

ylims!(ax1, 0.0, 60)
ylims!(ax2, 0, 60)

# xlims!(ax2, 0, 100)
# ylims!(ax2, 0, 50)
# axislegend(ax1, labelsize = 25, position = :lb, unique = true)
fig
