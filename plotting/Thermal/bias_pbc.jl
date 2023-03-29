include("../../src/main.jl")

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 2.0
λ = 1.0
τ = 1000
# ωTs = reverse([0.0, 5.0, 25.0, 100.0, 250.0])
ωTs = reverse([1.0, 5.0, 10.0, 25.0])

# Load memory kernel and thermal trajectory
colors = reverse([my_red, my_vermillion, my_green, my_blue, my_black])
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
tTraj_paths = ["precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ1000_lmax300_modes50000.jld2" for ωT in ωTs] 


# Speed of particles over time
function particle_speed(data; unfold = false, box = (-Inf, Inf), periodic = false)

    if unfold 
        σs = traj_unfold(data, box, periodic = periodic)
    else
        σs = data.σs
    end

    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = (σs[:, 2:end] .- σs[:, 1:end-1]) ./ δ

    return speeds
end

## Mean squared deviation of particle trajectories 
function particle_mean_displacement(data, box; periodic = false)
    σs_unfolded = traj_unfold(data, box, periodic = periodic, minus_init = true)
    return (vec(std(σs_unfolded, dims = 1)))
end



## Plotting particle trajectories
box_size = 250
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)

σ0 = [box[1] + α]
bias_val = 0.01

# Finite speeds with bias
# for tTraj_path in tTraj_paths
# tTraj = load_object(tTraj_path)
#     for speed in range(10, 60, step = 10)
#         σdot0 = [speed]
#         data = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, 10, τ; threads = true, box = box, periodic = true, bias = bias_val)
#         save_object("data/Thermal/Bias/Bias$(bias_val)_speed$(speed)_ωT$(tTraj.ωT)_mem10.jld2", data)
#     end
# end


# Finite speeds with bias with newly generated tTrajs

## Precompute the thermal trajectories
d = 60
τmax = 150                     # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_modes = 5000                 # Number of chain masses for simulating ρ0

qa_s = 2 * pi .* (1:n_modes) / n_modes
ωs = ω.(ωmax, qa_s ./ 2)
lmax = 300
ε = reduce(hcat, [exp.(1im * 2 * π / n_modes .* (1:n_modes) * g) for g = 1:lmax])

ωT = 5.0
speed_range = range(15, 60, step = 5.0)

# for tTraj_path in tTraj_paths
# tTraj = load_object(tTraj_path)
#     for speed in range(10, 60, step = 10)
#         σdot0 = [speed]
#         data = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, 10, τ; threads = true, box = box, periodic = true, bias = bias_val)
#         save_object("data/Thermal/Bias/Bias$(bias_val)_speed$(speed)_ωT$(tTraj.ωT)_mem10.jld2", data)
#     end
# end

nParticles = 25
## Zero speed with bias 
for tTraj_path in tTraj_paths
    tTraj = load_object(tTraj_path)
    σ0 = (right_boundary + left_boundary - α) / 2 .* ones(nParticles) + 10 * α * randn(nParticles)
    σdot0 = zeros(nParticles)

    data = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, 10, τ; threads = true, box = box, periodic = true, bias = bias_val)
    save_object("data/Thermal/Bias/Bias$(bias_val)_speed0.0_ωT$(tTraj.ωT)_mem10_τ$(τ).jld2", data)
end
    

# fig = Figure(resolution = (1600, 1200), font = "CMU Serif", fontsize = 36, figure_padding = 40)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\dot{\sigma}", title = L"\Phi = %$(Φ0), \, \lambda = %$(λ), \, \alpha = %$(α), \, \mathrm{Bias} = %$(bias_val)")

# # for speed in range(10, 60, step = 10)
# #     data = load_object("data/Thermal/Bias/Bias$(bias_val)_speed$(speed)_ωT25.0_memInf.jld2")
# #     lines!(ax1, data.τs[2:end], vec(particle_speed(data, unfold = true, box = box, periodic = true)))
# # end

# for ωT in ωTs 
#     data = load_object("data/Thermal/Bias/Bias$(bias_val)_speed0.0_ωT$(ωT)_mem10.jld2")
#     lines!(ax1, data.τs, particle_mean_displacement(data, box, periodic = true))
# end

# fig

# # Plotting particle trajectory variance / MSD 
# fig = Figure(resolution = (1600, 1200), font = "CMU Serif", fontsize = 36, figure_padding = 40)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = "Std Dev", title = L"\Phi = %$(Φ0), \, \lambda = %$(λ), \, \alpha = %$(α)")


# for ωT in ωTs 
#     data = load_object("data/Thermal/Bias/Bias$(bias_val)_speed0.0_ωT$(ωT)_mem10.jld2")
#     # data = load_object("data/Thermal/Thermalization/Particle25_Φ7.0_λ1.0_ωT$(ωT)_τ1000_mem1.jld2")
#     means = particle_mean_displacement(data, box, periodic = true)
#     # speeds = particle_speed(data; unfold = true, box = box, periodic = true)

#     # for particle in 1:25
#     #     lines!(ax1,  data.τs[2:end], speeds[particle, :], linewidth = 2, color = colors[findfirst(x -> x == ωT, ωTs)], label = L"\omega_T = %$(ωT)")
#     # end
#     lines!(ax1, data.τs, means, linewidth = 4, color = colors[findfirst(x -> x == ωT, ωTs)], label = L"\omega_T = %$(ωT)")
# end

# axislegend(ax1, position = :lt)

# xlims!(ax1, 0, nothing)
# fig