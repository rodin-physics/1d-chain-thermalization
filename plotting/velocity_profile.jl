include("../src/main.jl")

# Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 8.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")

function Δ_numeric2(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.25 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    return res
end

## Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

function speed_deviation(σdot, Φ0, μ)
    return sqrt(σdot^2 - (8*π^2*Φ0/μ))
end


function U_profile(r, Φ)
    return Φ*exp(-(r^2)/2*(λ^2))
end

## Plotting
fig = Figure(resolution=(1600, 800), font="CMU Serif", fontsize=24, figure_padding = 30)

# ax1 = Axis(fig[1, 1], xlabel=L"\rho", ylabel=L"U", title = "Potential Profile")
ax2 = Axis(fig[1, 1], xlabel=L"\tau", ylabel=L"\sigma", title = "Particle Position")
ax3 = Axis(fig[2, 1], xlabel=L"\tau", ylabel=L"\dot{\sigma}", title = "Velocities")
right_ax = Axis(fig[:, 2], xlabel=L"\dot{\sigma}", ylabel=L"\dot{\sigma} - \dot{\sigma}_\mathrm{ext}", title = "Velocity Fluctuations")

# Potential Profile
# ρs = range(-20, 20, step = 0.01)
# lines!(ax1, ρs, U_profile.(ρs, Φ0), color = my_vermillion, linewidth = 3, label = "Repulsive")
# lines!(ax1, ρs, U_profile.(ρs, -Φ0), color = my_blue, linewidth = 3, label = "Attractive")

# Computation
σdot0 = 40
res_rep = Δ_numeric2(σdot0, 5.5 * α, Φ0, λ, system)
res_att = Δ_numeric2(σdot0, 5.5 * α, -Φ0, λ, system)

# Particle Position
hlines!(ax2, [6 * α], linestyle = :dash, linewidth = 2, color = my_black)
vlines!(ax2, [res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))]], color = my_vermillion, linestyle = :dash, linewidth = 2)
vlines!(ax2, [res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))]], color = my_blue, linestyle = :dash, linewidth = 2)
lines!(ax2, res_rep.τs, vec(res_rep.σs), color = my_vermillion, linewidth = 3)
lines!(ax2, res_att.τs, vec(res_att.σs), color = my_blue, linewidth = 3)
Label(fig[1, 1, TopLeft()], "(a)", textsize = 30)

# Velocities
(xs_rep, speed_rep) = particle_speed(res_rep)
(xs_att, speed_att) = particle_speed(res_att)

vlines!(ax3, [res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))]], color = my_vermillion, linestyle = :dash, linewidth = 2)
vlines!(ax3, [res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))]], color = my_blue, linestyle = :dash, linewidth = 2)
hlines!(ax3, [sqrt(σdot0^2 - 8*Φ0*π^2), sqrt(σdot0^2 + 8*Φ0*π^2)])
lines!(ax3, xs_rep, speed_rep, linewidth = 3, color = my_vermillion, label = "Repulsive")
lines!(ax3, xs_att, speed_att, linewidth = 3, color = my_blue, label = "Attractive")
Label(fig[2, 1, TopLeft()], "(b)", textsize = 30)

# Speed Fluctuations
speeds = range(sqrt(8*Φ0*π^2/μ), 150, step = 0.5)
speeds2 = range(0, 150, step = 0.5)
lines!(right_ax, speeds, speeds .- speed_deviation.(speeds, Φ0, μ), color = my_vermillion, label = "Repulsive", linewidth = 3)
lines!(right_ax, speeds2, speeds2 .- speed_deviation.(speeds2, -Φ0, μ), color = my_blue, label = "Attractive", linewidth = 3)
Label(fig[:, 2, TopLeft()], "(c)", textsize = 30)

# xlims!(ax1, -20, 20)
xlims!(ax2, 0.0, 1 * (α / σdot0))
xlims!(ax3, 0.0, 1 * (α / σdot0))
xlims!(right_ax, 0.0, 150)
axislegend(right_ax, labelsize = 30, position = :rt)
fig
