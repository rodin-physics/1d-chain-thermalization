include("../../src/main.jl")

ωTs = [0.0, 5.0, 25.0]
colors = [my_blue, my_vermillion, my_yellow] |> reverse

system = load_object("precomputed/systems/System_ωmax10_d60_l1000.jld2")
d = 60
τ = 200                             # Simulation time
δ = system.δ                        # Time step
α = 40                              # Distance between chain atoms
μ = 1
mem = Inf
σdot0 = [120]
bias = 0.0
Φ0 = 20
λ = 4

n_pts = floor(τ / δ) |> Int
τs = δ .* (1:n_pts)

numTraj = 5


## Plotting 
fig = Figure(resolution = (1200, 1600), font = "CMU Serif", fontsize = 36, figure_padding = 40)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma - \sigma_0", yticks = 0:5000:10000)
ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma - \sigma_0", yticks = 0:5000:10000)

for ωT in ωTs
    data_rep = readdlm("data/Thermal/Single/Single$(numTraj)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(ωT)_tau$(τ).dat")

    data_att = readdlm("data/Thermal/Single/Single$(numTraj)_sigmadot0$(σdot0)_Mem$(mem)_lambda$(λ)_Phi-$(Φ0)_mu$(μ)_d$(d)_bias$(bias)_omegaT$(ωT)_tau$(τ).dat")

    for traj_ind in 1:numTraj
        lines!(ax1, τs, data_rep[traj_ind, :] .- data_rep[traj_ind, 1], color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
        lines!(ax2, τs, data_att[traj_ind, :] .- data_att[traj_ind, 1], color = colors[findfirst(x -> x == ωT, ωTs)], linewidth = 4, label = L"\omega_T = %$(ωT)")
    end

end

axislegend(ax1, position = :rb, merge = true, unique = true)
axislegend(ax2, position = :rb, merge = true, unique = true)


xlims!(ax1, 0, τ)
xlims!(ax2, 0, τ)

ylims!(ax1, 0, 12000)
ylims!(ax2, 0, 12000)

fig