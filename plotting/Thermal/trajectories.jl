include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 20.0
λ = 4.0
τ = 400
# ωTs = reverse([0.0, 5.0, 25.0, 100.0, 250.0])
ωTs = reverse([1.0, 5.0, 10.0, 25.0, 100.0])

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


# Define box
box_size = 250
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)

# num_trajs = 996 * 0.5 |> Int
num_trajs = 249
particles_per_traj = box_size - 1
batches = ones(num_trajs ÷ particles_per_traj) * particles_per_traj
if num_trajs % particles_per_traj != 0 
    batches = vcat(batches, num_trajs % particles_per_traj)
end

## Precompute the thermal trajectories
d = 60
τmax = 800                     # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_modes = 5000                 # Number of chain masses for simulating ρ0

qa_s = 2 * pi .* (1:n_modes) / n_modes
ωs = ω.(ωmax, qa_s ./ 2)
lmax = 300
ε = reduce(hcat, [exp.(1im * 2 * π / n_modes .* (1:n_modes) * g) for g = 1:lmax])

ωT = 10.0
init_speed = 120
x_range = range(20, 100, length = 70)
speed_range = range(30, 80, step = 1.0)

function full_trajectory(mem)
    # ## Plotting 
    fig = Figure(resolution = (1600, 1200), font = "CMU Serif", fontsize = 36, figure_padding = 40)
    ax1 = Axis(fig[1, 1], xlabel = L"\dot{\sigma}", ylabel = L"\Delta", title = L"\Phi = %$(Φ0), \, \lambda = %$(λ), \, \alpha = %$(α)")

    mean_data = zeros(length(x_range))
    var_data = zeros(length(x_range))
    p = Progress(length(x_range))
    Threads.@threads for ii in eachindex(x_range)
        mean_data[ii] =  Δ_thermal_analytic(x_range[ii], Φ0, λ, ωmax, ωT)
        var_data[ii] =  Δ_thermal_variance(x_range[ii], Φ0, λ, ωmax, ωT)
        next!(p)
    end

    lines!(ax1, x_range, mean_data, linewidth = 3, color = my_blue)
    band!(ax1, x_range, mean_data .- var_data, mean_data .+ var_data, color = (my_blue, 0.3))

    data_xs = Float64[]
    data_ys = Float64[] 
    for traj_ind in eachindex(batches)
        # Generate thermal trajectory
        ζs = ζq.(ωs, ωT)
        ϕs = 2 * π * rand(n_modes)

        res = zeros(lmax, n_pts)
        pr = Progress(n_pts)

        Threads.@threads for n = 1:n_pts
            res[:, n] = ρH(n * δ, ζs, ϕs, ωs, ε) |> real
            next!(pr)
        end

        curr_traj = ThermalTrajectory(ωmax, δ, res, ωT)

        @showprogress for start in (left_boundary .+ (1:batches[traj_ind]) .* α)
            res = motion_solver_test(system, Φ0, λ, α, [start], [init_speed], μ, curr_traj, mem, τmax, periodic = true, check_stuck = true, box = box, threads = true)

            (xs, ys) = Δ_traj_unfold(res, box; periodic = true)
            append!(data_xs, xs)
            append!(data_ys, ys)
            scatter!(ax1, xs, ys, markersize = 8, color = (my_black, 0.3))
        end
    end

    writedlm("data/Thermal/Full_Trajectory/deltas_mem$(mem)_ωT$(ωT)_Φ$(Φ0)_λ$(λ)_α$(α).dat", (data_xs, data_ys))

    save("traj_fig_mem$(mem)_test.pdf", fig)

end

full_trajectory(100)
