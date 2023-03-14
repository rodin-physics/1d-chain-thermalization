include("../../src/main.jl")
using Random

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 10.0
λ = 1.0
τ = 500
# ωTs = reverse([0.0, 5.0, 25.0, 100.0, 250.0])
ωTs = [1.0, 2.0, 5.0, 10.0, 25.0] |> reverse
# ωTs = [1.0, 5.0, 10.0, 25.0, 100.0] |> reverse

# Load memory kernel and thermal trajectory
colors = [my_blue, my_green, my_sky, my_orange, my_vermillion] |> reverse
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
tTraj_paths = ["precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ1000_lmax300_modes50000.jld2" for ωT in ωTs] 


## Check if a sign change has occured in an array 
function sign_change(arr)
    if length(arr) == 0
        error("Array is empty")
    end

    return !all( ==(sign(arr[1])), sign.(arr))
end


# Speed of particles over time
function particle_speed(data; unfold = false, box = (-Inf, Inf))

    if unfold 
        σs = traj_unfold(data, box)
    else
        σs = data.σs
    end

    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = (σs[:, 2:end] .- σs[:, 1:end-1]) ./ δ

    return speeds
end


## Flight distances for every mobile particle 
function flight_distances(sol::SystemSolution; unfold = false, box = (-Inf, Inf), ret_times = false)
    num_P = size(sol.σs, 1)
    τs = sol.τs
    ρs = sol.ρs
    σs = sol.σs
    speeds = particle_speed(sol, unfold = unfold, box = box)

    # Initialise arrays
    curr_lens = zeros(num_P)
    dists = [Int64[] for _ in 1:num_P]
    prev_hops = ones(Int64, num_P)
    times = [Float64[] for _ in 1:num_P]

    for ind in 1:(length(τs)-1)
        
        # Determine position of mobile particles in relation to chain
        curr = searchsortedfirst.(Ref(ρs[:,ind]), σs[:, ind])
        nxt = searchsortedfirst.(Ref(ρs[:, (ind + 1)]), σs[:, ind+1])

        # Compare initial and final positions
        particles = findall(x -> x != 0, (nxt .- curr))

        # A hopping event has occured
        if !isempty(particles)
            for part in particles 
                if curr_lens[part] == 0 || (sign_change(speeds[part, prev_hops[part]:ind]) == false)
                    
                    curr_lens[part] += 1

                # Particle did get stuck between the last hop event and now 
                elseif (τs[ind]-τs[prev_hops[part]] > 0.1)
                # else
                    # Record time it took for next hop to occur           
                    push!(times[part], τs[ind] - τs[prev_hops[part]])

                    # Record flight distance
                    push!(dists[part], curr_lens[part])

                    # Reset flight distance
                    curr_lens[part] = 1

                    # Update the latest hop event time index
                    prev_hops[part] = ind

                end
            end
        end
    end

    if ret_times 
        return times 
    else
        return dists
    end
end

## Cluster sorted vector based on tolerance 
function cluster_sizes(arr, tol)
    lens = Int64[]
    times = [arr[1]]
    curr_len = 0

    for ind in eachindex(arr)[2:end]
        curr_len += 1
        if (arr[ind]-arr[ind-1]) > tol
            push!(times, arr[ind])
            push!(lens, curr_len)
            curr_len = 0 
        end
    end

    if (arr[end] - arr[end - 1]) > tol 
        push!(lens, 1)
    else 
        push!(lens, curr_len + 1)
    end

    # return (times, lens)
    return lens
end


## Mean squared deviation of particle trajectories 
function particle_RMSD(data, box)
    σs_unfolded = traj_unfold(data, box, minus_init = true)
    return (vec(std(σs_unfolded, dims = 1)), data.τs)
end


## Mean displacement of particles
function particle_mean_displacement(data, box; periodic = false)
    σs_unfolded = traj_unfold(data, box, periodic = periodic, minus_init = true)
    return (vec(mean(σs_unfolded, dims = 1)), data.τs)
end

## Plotting particle trajectories
box_size = 250
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)
data = load_object("data/Thermal/Thermalization/Particle25_Φ7.0_λ1.0_ωT1.0_τ1000_mem10.jld2")

# fig = Figure(resolution = (1600, 1200), font = "CMU Serif", fontsize = 36, figure_padding = 40)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma", title = L"\Phi = %$(data.Φ), \, \lambda = %$(data.λ), \, \alpha = %$(data.α)")

# for ωT in ωTs 
#     data = load_object("data/Thermal/Thermalization/Particle25_ωT$(ωT)_τ500_mem10.jld2")
#     τs = data.τs
#     num_P = size(data.σs, 1)
#     σs_test = traj_unfold(data, box)

#     for particle in 1:num_P
#         lines!(ax1, τs, σs_test[particle, :], linewidth = 2, color = colors[findfirst(x -> x == ωT, ωTs)], label = L"\omega_T = %$(ωT)")
#     end
# end

# axislegend(ax1, unique = true, merge = true)

# xlims!(ax1, 0, 500)
# ylims!(ax1, -25000, 25000)
# fig


## Plotting particle trajectory variance / MSD 
fig = Figure(resolution = (1600, 1200), fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = "Std Dev",)

D0 = 30.0

# for ωT in ωTs 
#     data = load_object("data/Thermal/Thermalization/Particle25_Φ7.0_λ1.0_ωT$(ωT)_τ1000_mem10.jld2")
#     (msd_vals, times) = particle_RMSD(data, box)
#     lines!(ax1, times, msd_vals, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")

#     lines!(ax1, times[1:100:end], D0 * sqrt(data.ωT) * exp(-data.Φ / data.ωT) .* sqrt.(times[1:100:end]), linewidth = 3, linestyle = :dash, color = colors[findfirst(x -> x == ωT, ωTs)])
# end

for ωT in ωTs 
    data = load_object("data/Thermal/Thermalization/Particle25_Φ7.0_λ1.0_ωT$(ωT)_τ1000_mem10.jld2")
    (msd_vals, times) = particle_mean_displacement(data, box)
    lines!(ax1, times, msd_vals, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")

    # lines!(ax1, times[1:100:end], D0 * sqrt(data.ωT) * exp(-data.Φ / data.ωT) .* sqrt.(times[1:100:end]), linewidth = 3, linestyle = :dash, color = colors[findfirst(x -> x == ωT, ωTs)])
end

axislegend(ax1, position = :lt)

# xlims!(ax1, 0, 1000)
# ylims!(ax1, 0, nothing)
fig


## Plotting single particle trajectory
# test_data = load_object("data/Thermal/Thermalization/Particle25_ωT0.0_τ500_mem10.jld2")

# fig = Figure(resolution = (1600, 1200), font = "CMU Serif", fontsize = 36, figure_padding = 40)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma", title = L"\Phi = %$(test_data.Φ), \, \lambda = %$(test_data.λ), \, \alpha = %$(test_data.α), \, \omega_T = %$(test_data.ωT)")

# part_ind = 1
# hops = hopping_freqs(test_data)
# flights = cluster_sizes(hops[1], 0.5)

# vlines!(ax1, hops[part_ind])
# lines!(ax1, test_data.τs, test_data.σs[part_ind,:], linewidth = 3, color = my_black)

# xlims!(ax1, 80, 85)
# fig


## Plotting density plots of flight lengths + times 

# fig = Figure(resolution = (1600, 1200), fontsize = 36, figure_padding = 40)
# ax1 = Axis(fig[1, 1], xlabel = "Log(Flight Length)", ylabel = "Log(Density)")
# # ax1 = Axis(fig[1, 1], xlabel = "Time until Activation (τ)", ylabel = "Density")

# for ωT in ωTs
#     data = load_object("data/Thermal/Thermalization/Particle25_Φ7.0_λ1.0_ωT$(ωT)_τ1000_mem10.jld2")

#     ## FLIGHT LENGTHS
#     flight_lengths = vcat(flight_distances(data, unfold = true, box = box)...)

#     # density!(ax1, flight_lengths, label = L"\omega_T = %$(ωT)", color = (colors[findfirst(x -> x == ωT, ωTs)], 0.4))
#     hist_fit = fit(Histogram, flight_lengths, 1:2:100)
#     hist_fit = normalize(hist_fit, mode = :pdf)
#     scatter!(
#             ax1,
#             # ((hist_fit.edges[1])[1:end-1]).^(-2/3) .* exp.(-(0.5 / ωT) * (((hist_fit.edges[1])[1:end-1]).^(1/3))),
#             log.((hist_fit.edges[1])[1:end-1]),
#             log.(hist_fit.weights),
#             markersize = 18,
#             label = L"\omega_T = %$(ωT)",
#             color = colors[findfirst(x -> x == ωT, ωTs)]
#         )

    # hist!(ax1, flight_lengths, strokewidth = 2, normalization = :pdf, 
    # bins = 1:1:50, label = L"\omega_T = %$(ωT)", color = colors[findfirst(x -> x == ωT, ωTs)], scale_to=0.6, offset=findfirst(x -> x == ωT, ωTs))

    ## FLIGHT TIMES
    # flight_times = vcat(flight_distances(data, unfold = true, box = box, ret_times = true)...)
    # density!(ax1, flight_times, label = L"\omega_T = %$(ωT)", color = (colors[findfirst(x -> x == ωT, ωTs)], 0.4), strokearound = true, strokewidth = 4, strokecolor = colors[findfirst(x -> x == ωT, ωTs)], npoints = 300)
    
    # hist_fit = fit(Histogram, flight_times, 0.01:0.5:20)
    # hist_fit = normalize(hist_fit, mode = :pdf)
    # scatter!(
    #         ax1,
    #         # ((hist_fit.edges[1])[1:end-1]).^(-2/3) .* exp.(-(0.5 / ωT) * (((hist_fit.edges[1])[1:end-1]).^(1/3))),
    #         log.((hist_fit.edges[1])[1:end-1]),
    #         log.(hist_fit.weights),
    #         markersize = 18,
    #         label = L"\omega_T = %$(ωT)",
    #         color = colors[findfirst(x -> x == ωT, ωTs)]
    #     )

# end
# lines!(ax1, range(0, 4.3, length = 100), -1.5 .* range(0, 4.3, length = 100) .- 1.1, linewidth = 4, color = my_black, linestyle = :dash)


# axislegend(ax1, position = :rt)
# # xlims!(ax1, 0, 20)
# # ylims!(ax1, 0, nothing)
# # save("Loglog_plot.pdf", fig)
# fig