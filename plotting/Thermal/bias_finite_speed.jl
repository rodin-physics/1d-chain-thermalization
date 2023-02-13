include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
Φ0 = 2.0
λ = 4.0
τ = 150
bias_val = 0.01

ωT = 0.0
ωTs = reverse([1.0, 5.0, 10.0, 25.0, 100.0])
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
speed_range = range(15, 60, step = 1.0)

box_size = 250
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)

nParticles = 25
total_num = nParticles * 40
particles_per_traj = nParticles * 4

v_initials = Float64[]
v_finals = Float64[]

## Unfold particle trajectories given a confining box
function traj_unfold(data, box; periodic = false)
    σs_final = copy(data.σs)
    num_P = size(data.σs, 1)

    for particle in 1:num_P
        if periodic 
            for st in 2:lastindex(σs_final[particle, :])
                if (σs_final[particle, st] - σs_final[particle, st-1]) > 0.9 * (box[2] - box[1])
                    σs_final[particle, st:end] = σs_final[particle, st:end] .- (box[2] - box[1])
                elseif (σs_final[particle, st] - σs_final[particle, st-1]) < -0.9 * (box[2] - box[1])
                    σs_final[particle, st:end] = σs_final[particle, st:end] .+ (box[2] - box[1])
                end
            end
        else 
            τ_ids = findall(x -> x < box[1] || x > box[2], σs_final[particle,:])

            for τ_id in τ_ids
                σs_final[particle, (τ_id+1):end] = 2 * σs_final[particle, τ_id] .- σs_final[particle, (τ_id+1):end]
            end
        end
        σs_final[particle,:] = σs_final[particle,:] .- σs_final[particle,1]
    end

    return σs_final
end

# Thermal trajectory parameters
d = 60
τmax = 150                     # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_modes = 5000                 # Number of chain masses for simulating ρ0

qa_s = 2 * pi .* (1:n_modes) / n_modes
ωs = ω.(ωmax, qa_s ./ 2)
lmax = 260
ε = reduce(hcat, [exp.(1im * 2 * π / n_modes .* (1:n_modes) * g) for g = 1:lmax])


function generate_final_speeds(v_init, v_final)
    for traj_ind in 1:(total_num ÷ particles_per_traj)
        ## Generate thermal trajectory 
        ζs = ζq.(ωs, ωT)
        ϕs = 2 * π * rand(n_modes)
    
        res = zeros(lmax, n_pts)
        pr = Progress(n_pts)
    
        Threads.@threads for n = 1:n_pts
            res[:, n] = ρH(n * δ, ζs, ϕs, ωs, ε) |> real
            next!(pr)
        end
        
        curr_traj = ThermalTrajectory(ωmax, δ, res, ωT)
    
        for run_ind in 1:(particles_per_traj ÷ nParticles)
            # Start nParticles with random speeds and positions within the box
            init_speeds = rand(speed_range, nParticles)
            init_pos = (box[2] + box[1] - α) / 2 .* ones(nParticles) + 10 * α * randn(nParticles)
    
            data = motion_solver(system, Φ0, λ, α, init_pos, init_speeds, μ, curr_traj, 10, τ; threads = true, box = box, periodic = true, bias = bias_val)
    
            # Get the final speeds after 150τ
            data_unfold = traj_unfold(data, box, periodic = true)
            final_speeds = (data_unfold[:, end] .- data_unfold[:, end-1]) ./ δ
    
            v_init = vcat(v_init, init_speeds)
            v_final = vcat(v_final, final_speeds)
        
        end
    
    end

    ## Save data 
    writedlm("Initial_speeds_α$(α)_Φ$(Φ0)_λ$(λ)_bias$(bias_val)_T$(ωT).dat", v_init)
    writedlm("Final_speeds_α$(α)_Φ$(Φ0)_λ$(λ)_bias$(bias_val)_T$(ωT).dat", v_final)   

end

generate_final_speeds(v_initials, v_finals)



