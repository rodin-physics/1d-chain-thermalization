include("../../src/main.jl")

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 2.0
λ = 1.0
τ = 1000

box_size = 250
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)

ωTs = [1.0, 5.0, 10.0, 25.0] |> reverse
colors = [my_blue, my_sky, my_orange, my_vermillion] |> reverse


## Mean squared deviation of particle trajectories 
function particle_RMSD(data, box; periodic = false)
    σs_unfolded = traj_unfold(data, box, periodic = periodic, minus_init = true)
    return (vec(std(σs_unfolded, dims = 1)), data.τs)
end


## Mean displacement of particles
function particle_mean_displacement(data, box; periodic = false)
    σs_unfolded = traj_unfold(data, box, periodic = periodic, minus_init = true)
    return (vec(mean(σs_unfolded, dims = 1)), data.τs)
end


## Plotting 
# fig = Figure(resolution = (2000, 1200), fontsize = 36, figure_padding = 40)
# ax1 = Axis(fig[1, 1], ylabel = "Mean Displacement", xticklabelsvisible = false, title = "No bias")
# ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = "Standard Deviation of Displacement")
# ax3 = Axis(fig[1, 2], xticklabelsvisible = false, title = "Bias = 0.0025")
# ax4 = Axis(fig[2, 2], xlabel = L"\tau",)

# D0 = 30.0

# # Zero bias data
# for ωT in ωTs 
#     data = load_object("data/Thermal/Thermalization/Particle25_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem10.jld2")
#     (rmsd_vals, times) = particle_RMSD(data, box)
#     (means, times) = particle_mean_displacement(data, box)

#     lines!(ax1, times, means, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
#     lines!(ax2, times, rmsd_vals, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")

#     lines!(ax2, times[1:100:end], D0 * sqrt(data.ωT) * exp(-data.Φ / data.ωT) .* sqrt.(times[1:100:end]), linewidth = 3, linestyle = :dash, color = colors[findfirst(x -> x == ωT, ωTs)])
# end 

# # Bias data
# bias_val = 0.0025 
# for ωT in ωTs
#     data = load_object("data/Thermal/Bias/Bias0.0025_speed0.0_ωT$(ωT)_mem10.jld2")
#     (rmsd_vals, times) = particle_RMSD(data, box, periodic = true)
#     (means, times) = particle_mean_displacement(data, box, periodic = true)

#     lines!(ax3, times, means, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
#     lines!(ax4, times, rmsd_vals, color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")

#     lines!(ax4, times[1:100:end], D0 * sqrt(data.ωT) * exp(-data.Φ / data.ωT) .* sqrt.(times[1:100:end]), linewidth = 3, linestyle = :dash, color = colors[findfirst(x -> x == ωT, ωTs)])
# end



# fig[:, 3] = Legend(fig, ax1, "Temperature")

# xlims!(ax1, 0, τ)
# xlims!(ax2, 0, τ)
# xlims!(ax3, 0, τ)
# xlims!(ax4, 0, τ)

# Label(fig[1,1, TopLeft()], "(a)", font = :bold)
# Label(fig[2,1, TopLeft()], "(b)", font = :bold)
# Label(fig[1,2, TopLeft()], "(c)", font = :bold)
# Label(fig[2,2, TopLeft()], "(d)", font = :bold)

# fig

test = load_object("data/Thermal/Thermalization/Particle25_Φ2.0_λ1.0_ωT25.0_τ1000_mem10.jld2")
fig = Figure()
ax = Axis(fig[1,1])
test_traj = traj_unfold(test, box, minus_init = true)
test_means = vec(mean(test_traj, dims = 1))

for ii in 1:25
    lines!(ax, test.τs, test_means)
    
    # lines!(ax, test.τs, test_traj[ii,:] |> vec)
end

fig