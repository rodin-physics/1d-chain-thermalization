include("../../src/main.jl")
using Random
using BinnedStatistics

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 4.0
ωT = 10.0
# x_range = range(√((8*π^2*Φ0)/μ), 120, length = 80)
# x_range_capture = shuffle(range(1, √((8*π^2*Φ0)/μ), length = 40))

## Plotting 
fig = Figure(resolution = (2400, 800), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}", ylabel = L"\Delta", title = L"Memory $\tau_0$ = 1", yticks = -50:25:50, xgridvisible = false, ygridvisible = false)
ax2 = Axis(fig[1, 2], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"$\tau_0$ = 10", yticks = -50:25:50, xgridvisible = false, ygridvisible = false)
ax3 = Axis(fig[1, 3], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"$\tau_0$ = 100", yticks = -50:25:50, xgridvisible = false, ygridvisible = false)

# Prediction 
# mean_data = zeros(length(x_range))
# var_data = zeros(length(x_range))
# mean_data2 = zeros(length(x_range_capture))
# var_data2 = zeros(length(x_range_capture))

# p = Progress(length(x_range))
# Threads.@threads for ii in eachindex(x_range)
#     mean_data[ii] =  Δ_thermal_analytic(x_range[ii], Φ0, λ, ωmax, ωT)
#     var_data[ii] =  Δ_thermal_variance(x_range[ii], Φ0, λ, ωmax, ωT)
#     next!(p)
# end

# p = Progress(length(x_range_capture))
# Threads.@threads for ii in eachindex(x_range_capture)
#     mean_data2[ii] =  Δ_thermal_analytic(x_range_capture[ii], Φ0, λ, ωmax, ωT)
#     var_data2[ii] =  Δ_thermal_variance(x_range_capture[ii], Φ0, λ, ωmax, ωT)
#     next!(p)
# end
# indices = sortperm(x_range_capture)
# x_range_capture = x_range_capture[indices]
# mean_data2 = mean_data2[indices]
# var_data2 = var_data2[indices]

# xs = vcat(x_range_capture, x_range)
# ys = vcat(mean_data2, mean_data)
# ys2 = vcat(var_data2, var_data)

# writedlm("data/Thermal/Full_Trajectory/MeanΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat", (xs, ys))
# writedlm("data/Thermal/Full_Trajectory/VarΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat", (xs, ys2))

mean_data = readdlm("data/Thermal/Full_Trajectory/MeanΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")
var_data = readdlm("data/Thermal/Full_Trajectory/VarΔ_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")

for ax in [ax1, ax2, ax3]
    vlines!(ax,√((8*π^2*Φ0)/μ), linewidth = 4, linestyle = :dash, color = my_black, label = "Capture Speed")
    band!(ax, mean_data[1,:], mean_data[2,:] .- var_data[2,:], mean_data[2,:] .+ var_data[2,:], color = (my_blue, 0.3), label = "Analytic")
    lines!(ax, mean_data[1,:], mean_data[2,:], linewidth = 4, color = my_blue, label = "Analytic")
 
    # band!(ax, x_range_capture, mean_data2 .- var_data2, mean_data2 .+ var_data2, color = (my_blue, 0.3))
    # lines!(ax, x_range_capture, mean_data2, linewidth = 3, color = my_blue)
end


# # Numerical data 
mem_ax_pairs = [(1, ax1), (10, ax2), (100, ax3)]


for pair in mem_ax_pairs
    data = readdlm("data/Thermal/Full_Trajectory/deltas_mem$(pair[1])_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")

    edges, centers, mean_res = binnedStatistic(data[1,:], data[2,:], statistic = :mean, nbins = 150)
    lines!(pair[2], centers, mean_res, linewidth = 4, color = my_red, label = "Numerical")

    edges, centers, var_res = binnedStatistic(data[1,:], data[2,:], statistic = :var, nbins = 150)

    band!(pair[2], centers, mean_res .- var_res, mean_res .+ var_res, color = (my_red, 0.4), label = "Numerical")

    xlims!(pair[2], 0, 140)
    ylims!(pair[2], -52, 52)
end

Label(fig[1,1, TopLeft()], "(a)", font = :bold, padding = (0, 0, 10, 0))
Label(fig[1,2, TopLeft()], "(b)", font = :bold, padding = (0, 0, 10, 0))
Label(fig[1,3, TopLeft()], "(c)", font = :bold, padding = (0, 0, 10, 0))


axislegend(ax3, merge = true, labelsize = 32, patchsize = (50, 50), height = 185)
# save("full_trajectory.pdf", fig)
fig

