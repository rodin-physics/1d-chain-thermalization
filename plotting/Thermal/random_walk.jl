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
x_range = range(√((8*π^2*Φ0)/μ), 120, length = 80)
x_range_capture = shuffle(range(1, √((8*π^2*Φ0)/μ), length = 40))


function speed_after_pass(init_speed)
    pred_mean = Δ_thermal_analytic(init_speed, Φ0, λ, ωmax, ωT)
    pred_var = Δ_thermal_variance(init_speed, Φ0, λ, ωmax, ωT)

    dist = Normal(pred_mean, √(pred_var))

    return √(init_speed^2 - (8*π^2/μ)*rand(dist))
end

function delta_after_pass(init_speed)
    pred_mean = Δ_thermal_analytic(init_speed, Φ0, λ, ωmax, ωT)
    pred_var = Δ_thermal_variance(init_speed, Φ0, λ, ωmax, ωT)

    dist = Normal(pred_mean, √(pred_var))

    return rand(dist)
end


## Trajectory loss based on random walk model
fig = Figure(resolution = (1600, 1200), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}", ylabel = L"\Delta", xgridvisible = false, ygridvisible = false)

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

mean_data = readdlm("data/Thermal/Full_Trajectory/MeanΔ_ωT10.0_Φ20.0_λ4.0_α40.dat")
var_data = readdlm("data/Thermal/Full_Trajectory/VarΔ_ωT10.0_Φ20.0_λ4.0_α40.dat")


band!(ax1, mean_data[1,:], mean_data[2,:] .- var_data[2,:], mean_data[2,:] .+ var_data[2,:], color = (my_blue, 0.3), label = "Analytic")
lines!(ax1, mean_data[1,:], mean_data[2,:], linewidth = 3, color = my_blue, label = "Analytic")


function random_walk(nPasses, nParticles, init_speed, num_ind)
    speed_data = Float64[]
    Δ_data = Float64[]

    speeds = repeat([init_speed], nParticles)

    @showprogress for _ in 1:nPasses
        stuck = findall(x -> x <= √(8*π^2*Φ0/μ), speeds)
        deleteat!(speeds, stuck)

        if isempty(speeds)
            break
        end

        Δs = delta_after_pass.(speeds)
        append!(speed_data, speeds)
        append!(Δ_data, Δs)
        scatter!(ax1, speeds, Δs, markersize=8, color=(my_black, 0.4))

        square_speeds = speeds.^2 .- (8 * π^2 / μ) .* Δs
        below_zero = findall(x -> x < 0, square_speeds)
        deleteat!(square_speeds, below_zero)

        speeds = sqrt.(square_speeds)
    end
    writedlm("data/Thermal/random_walk$(num_ind).dat", (speed_data, Δ_data))
end

# Threads.@threads for batch in 1:9
#     random_walk(800, 20, 120, batch)
# end

data = readdlm("data/Thermal/random_walk.dat")
edges, centers, mean_res = binnedStatistic(data[1,:], data[2,:], statistic = :mean, nbins = 100)
lines!(ax1, centers, mean_res, linewidth = 4, color = my_red, label = "Numerical")
edges, centers, var_res = binnedStatistic(data[1,:], data[2,:], statistic = :var, nbins = 100)
band!(ax1, centers, mean_res .- var_res, mean_res .+ var_res, color = (my_red, 0.4), label = "Numerical")

vlines!(ax1, [√(8*pi^2 * Φ0)], linestyle = :dash, linewidth = 4, color = my_black, label = "Capture Speed")
xlims!(ax1, 0, 140)
ylims!(ax1, -50, 50)
# save("random_walk.pdf", fig)

axislegend(ax1, merge = true, patchsize = (40, 40))
fig