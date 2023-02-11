include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
nChain = 3
Φ0 = 2.0
λ = 4.0
# Load memory kernel and thermal trajectory
system = load_object("precomputed/systems/System_ωmax10_d60_l10000_τmax5.jld2")

# Numerically calculate energy loss after single pass
function Δ_numeric(σ_dot, σ0, Φ0, λ, system, tTraj)
    δ = system.δ
    τ = 1.75 * (α / σ_dot)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ)

    σs = res.σs |> vec
    ρs = res.ρs

    # Find index closest to next midpoint
    chain_idx = searchsortedlast(ρs[:, 1], σs[1])
    mod_val = mod(σs[1], res.α)
    mob_final = findmin(abs.(σs .- (ρs[chain_idx+1, :] .+ mod_val)))[2]

    # Calculate final kinetic energy
    v_final = (σs[mob_final] - σs[mob_final-1]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

# Get distribution of energy loss values
function Δ_thermal(σ_dot, σ0, Φ0, λ, system, tTraj)

    function get_tTraj(ind, tTraj)
        return ThermalTrajectory(
            tTraj.ωmax,
            tTraj.δ,
            tTraj.ρHs[ind:ind+nChain, :],
            tTraj.ωT,
        )
    end
    nPts = length(1:(size(tTraj.ρHs)[1]-nChain))
    res = zeros(nPts)
    p = Progress(nPts)
    Threads.@threads for ii = 1:(size(tTraj.ρHs)[1]-nChain)
        res[ii] = Δ_numeric(σ_dot, σ0, Φ0, λ, system, get_tTraj(ii, tTraj))
        next!(p)
    end

    return res
end

## Computation 
σdot0 = 50.0
ωT = 0.0
colors = reverse([my_blue, my_green, my_vermillion, my_red, my_yellow])
ωTs = reverse([0.0, 2.0, 5.0, 10.0, 25.0])
init_pos = 2.5 * α



# tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ1000_lmax300_modes50000.jld2")
# Δs = Δ_thermal(σdot0, init_pos, Φ0, λ, system, tTraj)
# cd("data/Thermal/ωT$(ωT)")
# writedlm("delta_ωmax$(ωmax)_speed$(σdot0)_ωT$(ωT)_Pts$(length(Δs)).dat", Δs)
# cd("../../../")

Δs = readdlm("data/Thermal/ωT$(ωT)/delta_ωmax10_speed$(σdot0)_ωT$(ωT).dat") |> vec

## Plotting density compared to predicted distribution
fig = Figure(resolution = (1200, 1200), font = "CMU Serif", fontsize = 40, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel = L"\Delta", ylabel = "Density",
    title = L"\omega_T = %$(ωT), \, \Phi_0 = %$(Φ0), \,\lambda = %$(λ), \, \dot{\sigma}_0 = %$(σdot0)")

density!(ax1, Δs, npoints = 200, label = "Numerical")

pred_mean = Δ_thermal_analytic(σdot0, Φ0, λ, ωmax, ωT)
pred_var = Δ_thermal_variance(σdot0, Φ0, λ, ωmax, ωT)

density!(ax1, rand(Normal(pred_mean, √(pred_var)), 1000000), color = :transparent, strokearound = true, strokewidth = 4, strokecolor = my_black, label = "Predicted")

axislegend(ax1, labelsize = 35, position = :rt)
ylims!(0, nothing)
fig

## Plotting mean and variance 

# function Δ_per_temp(dir_name)
#     filenames = filter(x -> x[1] !== '.', readdir(dir_name))
#     speeds = Float32[]
#     Δs = Float32[] 
#     sdevs = Float32[]
#     for file in filenames  
#         # Obtain mean and error bars 
#         data = readdlm(joinpath(dir_name, file)) |> vec
#         append!(Δs, mean(data))
#         append!(sdevs, std(data))

#         # Parse the filename and get the speed 
#         names = filter(x -> x[1] == 's', split(file, "_"))
#         speed = parse(Float64, chopprefix(names[1], "speed"))
#         append!(speeds, speed)
#     end
#     return (speeds, Δs, sdevs)
# end

# speeds = range(15, 80, step = 5.0)
# ωTs = [0.0, 2.0, 5.0, 10.0, 25.0]

# fig = Figure(
#     resolution = (1200, 1600),
#     font = "CMU Serif",
#     fontsize = 40,
#     figure_padding = 30,
# )
# ax1 = Axis(
#     fig[1, 1],
#     # xlabel = "Speed",
#     ylabel = "Mean Loss",
#     title = L"\Phi_0 = %$(Φ0), \,\lambda = %$(λ)",
# )

# ax2 = Axis(
#     fig[2, 1],
#     xlabel = "Speed",
#     ylabel = "Std Dev of Loss",
# )

# vlines!(ax1, [√(8*pi^2 * 2)], linestyle = :dash, linewidth = 4, color = my_black)
# vlines!(ax2, [√(8*pi^2 * 2)], linestyle = :dash, linewidth = 4, color = my_black, label = "Capture Speed")

# x_range = range(5, 80, length = 50)
# for ωT in ωTs 
#     mean_data = zeros(length(x_range))
#     var_data = zeros(length(x_range))
#     p = Progress(length(x_range))

#     Threads.@threads for ii in eachindex(x_range)
#         mean_data[ii] =  Δ_thermal_analytic(x_range[ii], Φ0, λ, ωmax, ωT)
#         var_data[ii] =  Δ_thermal_variance(x_range[ii], Φ0, λ, ωmax, ωT)
#         next!(p)
#     end

#     lines!(ax1, x_range, mean_data, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
#     lines!(ax2, x_range, sqrt.(var_data), color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4)

#     (speeds, Δs, sdevs) = Δ_per_temp("data/Thermal/ωT$(ωT)/")

#     scatter!(ax1, speeds, Δs, color = colors[findfirst(x -> x == ωT, ωTs)], markersize = 16)
#     scatter!(ax2, speeds, sdevs, color = colors[findfirst(x -> x == ωT, ωTs)], markersize = 16)

# end

# axislegend(ax1, position = :rt)
# axislegend(ax2, position = :rt)

# xlims!(ax1, 0.0, nothing)
# xlims!(ax2, 0.0, nothing)

# ylims!(ax1, 0.0, nothing)
# ylims!(ax2, 0.0, nothing)

# save("Loss_Mean_Std.pdf", fig)
# fig

