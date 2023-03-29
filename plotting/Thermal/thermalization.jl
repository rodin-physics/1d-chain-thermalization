include("../../src/main.jl")
using Random

## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 2.0
λ = 1.0
τ = 1000
ωTs = [1.0, 2.0, 5.0, 10.0, 25.0, 100.0, 250.0]

# Load memory kernel and thermal trajectory
colors = [my_blue, my_green, my_sky, my_orange, my_vermillion, my_red, my_yellow, my_black]
system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
tTraj_paths = ["precomputed/rH/rH_ωmax10_d60_ωT$(ωT)_τ1000_lmax300_modes50000.jld2" for ωT in ωTs] 

# Get every particle's energy distribution 
function get_particle_energy_dist(data, box)
    δ = data.τs[2] - data.τs[1]
    num_P = size(data.σs, 1)
    ens = Float64[]
    σs_unfold = traj_unfold(data, box)

    for particle in 1:num_P
        # Get times at which particle is at midpoints for unfolded trajectory 
        curr_σs = σs_unfold[particle, :]
        # first_idx = searchsortedlast(data.ρs[:,1], curr_σs[1])
        # last_idx = maximum(curr_σs) / data.α |> floor |> Int

        # idx = [argmin(abs.(curr_σs .- ((n + 0.5) * data.α))) for n in first_idx:last_idx]
        # idx = filter(x -> x <= length(curr_σs) - 1, idx)

        # σ_dots = (curr_σs[idx.+1] - curr_σs[idx]) ./ δ
        # kin_en = σ_dots .^ 2 ./ 2 ./ data.ωT ./ (2 * pi)^2

        midpoint = floor.(curr_σs .- data.α / 2)
        idx = findall(x -> x != 0, midpoint[2:end] - midpoint[1:end-1])
        kin_en = ((curr_σs[idx] - curr_σs[idx.-1]) ./ δ) .^ 2 ./ 2 ./ effective_T(data.ωmax, data.ωT) ./ (2 * pi)^2

        
        # pot_en = [sum((data.Φ .* exp.(-(data.ρs[:, t] .- data.σs[particle, t]) .^ 2 ./ (2 * data.λ^2))) ./ data.ωT) for t in eachindex(data.σs[particle, :])]
        # pot_en = pot_en
        # idx = findall(x -> x < 1e-3, pot_en)
        # ens = kin_en[idx]
        # ens = vcat(ens, kin_en + pot_en)
        ens = vcat(ens, kin_en)
    end

    hist_fit = fit(Histogram, ens |> vec, 0.0:0.02:5)
    hist_fit = normalize(hist_fit, mode = :pdf)

    return hist_fit
end



## Computation

box_size = 250              # Particle-confining box
left_boundary = 4.5 * α
right_boundary = left_boundary + box_size * α
box = (left_boundary, right_boundary)
nParticles = 25

for path in tTraj_paths
    σ0 = (right_boundary + left_boundary - α) / 2 .* ones(nParticles) + 10 * α * randn(nParticles)
    σdot0 = zeros(nParticles)

    tTraj = load_object(path)
    res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, 100, τ; threads = true, box = box)

    save_object("data/Thermal/Thermalization/Particle$(nParticles)_Φ$(Φ0)_λ$(λ)_ωT$(tTraj.ωT)_τ$(τ)_mem100.jld2", res)
end



## Plotting
# fig = Figure(resolution = (1600, 1200), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(fig[1, 1], xlabel = L"\mathcal{E} / \omega_T", ylabel = L"\ln(P(\mathcal{E}))", title = L"\Phi = %$(Φ0), \, \lambda = %$(λ), \, \alpha = %$(α)")

# # Line with slope -1 
# lines!(ax1, 0:5, -(0:5), color = my_black, linewidth = 4)

# for ωT in ωTs 
#     hist_data = get_particle_energy_dist(load_object("data/Thermal/Thermalization/Particle$(nParticles)_Φ$(Φ0)_λ$(λ)_ωT$(ωT)_τ$(τ)_mem10.jld2"), box)

#     scatter!(
#         ax1,
#         (hist_data.edges[1])[1:end-1],
#         log.(hist_data.weights),
#         # log.(hist_data.weights .* sqrt.((hist_data.edges[1])[1:end-1])),
#         color = colors[findfirst(x -> x == ωT, ωTs)],
#         label = L"\omega_T = %$(ωT)",
#         markersize = 12,)

# end


# axislegend(ax1)

# # save("thermalization_mem100_Φ$(Φ0)_λ$(λ).pdf", fig)
# fig