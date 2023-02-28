include("../../src/main.jl")
using Random

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 5.0
λ = 1.0
ωT = 10.0
# x_range = range(√((8*π^2*Φ0)/μ), 120, length = 80)
# x_range_capture = shuffle(range(1, √((8*π^2*Φ0)/μ), length = 40))

## Plotting 
fig = Figure(resolution = (2400, 800), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}", ylabel = L"\Delta", title = L"Memory $\tau_0$ = 1")
ax2 = Axis(fig[1, 2], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"$\tau_0$ = 10")
ax3 = Axis(fig[1, 3], xlabel = L"\dot{\sigma}", yticklabelsvisible = false, title = L"$\tau_0$ = 100")

# Prediction 
mean_data = zeros(length(x_range))
var_data = zeros(length(x_range))
mean_data2 = zeros(length(x_range_capture))
var_data2 = zeros(length(x_range_capture))

p = Progress(length(x_range))
Threads.@threads for ii in eachindex(x_range)
    mean_data[ii] =  Δ_thermal_analytic(x_range[ii], Φ0, λ, ωmax, ωT)
    var_data[ii] =  Δ_thermal_variance(x_range[ii], Φ0, λ, ωmax, ωT)
    next!(p)
end

p = Progress(length(x_range_capture))
Threads.@threads for ii in eachindex(x_range_capture)
    mean_data2[ii] =  Δ_thermal_analytic(x_range_capture[ii], Φ0, λ, ωmax, ωT)
    var_data2[ii] =  Δ_thermal_variance(x_range_capture[ii], Φ0, λ, ωmax, ωT)
    next!(p)
end
indices = sortperm(x_range_capture)
x_range_capture = x_range_capture[indices]
mean_data2 = mean_data2[indices]
var_data2 = var_data2[indices]

for ax in [ax1, ax2, ax3]
    band!(ax, x_range, mean_data .- var_data, mean_data .+ var_data, color = (my_blue, 0.3))
    lines!(ax, x_range, mean_data, linewidth = 3, color = my_blue)

    band!(ax, x_range_capture, mean_data2 .- var_data2, mean_data2 .+ var_data2, color = (my_black, 0.3))
    lines!(ax, x_range_capture, mean_data2, linewidth = 3, color = my_black)
end


# Numerical data 
mem_ax_pairs = [(1, ax1), (10, ax2)]

for pair in mem_ax_pairs
    data = readdlm("data/Thermal/Full_Trajectory/deltas_mem$(pair[1])_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat")

    pick_range = range(1, length(data[1,:]), step = 1)
    sampled = sample(pick_range, 20000)
    scatter!(pair[2], data[1,:][sampled], data[2,:][sampled], markersize = 8, color = (my_black, 0.3))
end

Label(fig[1,1, TopLeft()], "(a)", font = :bold)
Label(fig[1,2, TopLeft()], "(b)", font = :bold)
Label(fig[1,3, TopLeft()], "(c)", font = :bold)

xlims!(ax1, 0, 120)
ylims!(ax1, -16, 16)

xlims!(ax2, 0, 120)
ylims!(ax2, -16, 16)

xlims!(ax3, 0, 120)
ylims!(ax3, -16, 16)


fig

