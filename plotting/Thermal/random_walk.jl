include("../../src/main.jl")
using Random
## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 5.0
λ = 1.0
ωT = 10.0

x_range = range(20, 120, length = 70)
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
ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}", ylabel = L"\Delta")

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

band!(ax1, x_range, mean_data .- var_data, mean_data .+ var_data, color = (my_blue, 0.3))
lines!(ax1, x_range, mean_data, linewidth = 3, color = my_blue)

band!(ax1, x_range_capture, mean_data2 .- var_data2, mean_data2 .+ var_data2, color = (my_black, 0.3))
lines!(ax1, x_range_capture, mean_data2, linewidth = 3, color = my_black)


function random_walk(nPasses, nParticles, init_speed)
    speeds = repeat([init_speed], nParticles)

    @showprogress for _ in 1:nPasses
        stuck = findall(x -> x <= √(8*π^2*Φ0/μ), speeds)
        deleteat!(speeds, stuck)

        if isempty(speeds)
            break
        end

        Δs = delta_after_pass.(speeds)

        scatter!(ax1, speeds, Δs, markersize=8, color=(my_black, 0.4))

        square_speeds = speeds.^2 .- (8 * π^2 / μ) .* Δs
        below_zero = findall(x -> x < 0, square_speeds)
        deleteat!(square_speeds, below_zero)

        speeds = sqrt.(square_speeds)
    end

end

Threads.@threads for batch in 1:9
    random_walk(500, 20, 80)
end


xlims!(ax1, 0, 120)
ylims!(ax1, -16, 16)
save("random_walk.pdf", fig)
# fig