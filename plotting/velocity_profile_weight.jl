include("../src/main.jl")
using LambertW

# Parameters
α = 40
μ = 1
nChain = 10
ωmax = 10
Φ0 = 8.0
λ = 4.0
system = load_object("precomputed/systems/System_ωmax10_d1000_l20.jld2")

function Δ_numeric2(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.5 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)

    return res
end

## Speed of particle over time
function particle_speed_filtered(data, σ_next)
    σs = data.σs |> vec
    mid_pos = findmin(abs.(σs .- σ_next))[2]

    σs = σs[1:mid_pos]
    τs = data.τs[1:mid_pos]
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

function speed_deviation(σdot, Φ0, μ)
    return sqrt(σdot^2 - (8*π^2*Φ0/μ))
end


function U_profile(r, Φ, λ)
    return Φ*exp(-r^2/(2*λ^2))
end

function productlog_approx(τ, σdot, Φ)
    return sqrt(σdot^2 + (2*λ^2/τ^2)*lambertw(-8*Φ*π^2/μ *τ^2/(2*λ^2) * exp(-σdot^2*τ^2/(2*λ^2))))
end


## Plotting
σdot0 = 26
fig = Figure(resolution=(1600, 800), font="CMU Serif", fontsize=24, figure_padding = 30, supertitle = L"\Phi_0 = $(Φ0), \lambda = $(λ), \dot{\sigma}_0 = $(σdot0)")

ax2 = Axis(fig[1, 2], xlabel=L"\dot{\sigma}", ylabel="Density", title = "Attractive Density")
ax3 = Axis(fig[2, 2], xlabel=L"\dot{\sigma}", ylabel="Density", title = "Repulsive Density")
ax1 = Axis(fig[:, 1], xlabel=L"\tau", ylabel=L"\dot{\sigma}", title = "Particle Velocity")
Label(fig[1, 1:2, Top()], L"\Phi_0 = %$(Φ0), \, \lambda = %$(λ), \, \dot{\sigma}_0 = %$(σdot0)", textsize = 40, padding = (0, 0, 35, 0))
# Computation
xs = range(-α/σdot0/2, α/σdot0/2, length = 10000)
xs2 = range(-α/2, α/2, length = 500)
res_rep = Δ_numeric2(σdot0, 5.5 * α, Φ0, λ, system)
res_att = Δ_numeric2(σdot0, 5.5 * α, -Φ0, λ, system)

# Velocities
(xs_rep, speed_rep) = particle_speed_filtered(res_rep, 6.5*α)
(xs_att, speed_att) = particle_speed_filtered(res_att, 6.5*α)

# Vertical lines
vlines!(ax1, [res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))]], color = my_vermillion, linestyle = :dash, linewidth = 2)
vlines!(ax1, [res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))]], color = my_blue, linestyle = :dash, linewidth = 2)

# Maximum and minimum deviation
hlines!(ax1, [sqrt(σdot0^2 + 8*Φ0*π^2)], color = my_blue, linestyle = :dash, linewidth = 2)
hlines!(ax1, [sqrt(σdot0^2 - 8*Φ0*π^2)], color = my_vermillion, linestyle = :dash, linewidth = 2)

# Repulsive
lines!(ax1, xs_rep, speed_rep, linewidth = 4, color = my_vermillion, label = "Repulsive")
# lines!(ax1, xs .+ res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))], productlog_approx.(xs,σdot0, Φ0), color = my_green, linewidth = 2.5, label = "ProductLog Sol")
# lines!(ax1, xs .+ res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))], U_profile.(xs, sqrt(σdot0^2 - 8*Φ0*π^2)-σdot0, λ/(σdot0-Φ0)) .+ σdot0 , linewidth = 2.5, color = my_black)
test_rep = map(x -> sqrt(σdot0^2 - 8*π^2*U_profile(x, Φ0, λ/σdot0)), xs)
lines!(ax1, xs .+ res_rep.τs[argmin(abs.(vec(res_rep.σs) .- 6*α))], test_rep, linewidth = 2.5, color = my_black)

# Attractive
# lines!(ax1, xs .+ res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))], productlog_approx.(xs,σdot0, -Φ0), color = my_green, linewidth = 2.5)
test_att = map(x -> sqrt(σdot0^2 - 8*π^2*U_profile(x, -Φ0, λ/σdot0)), xs)
lines!(ax1, xs .+ res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))], test_att, linewidth = 2.5, color = my_black)

# lines!(ax1, xs .+ res_att.τs[argmin(abs.(vec(res_att.σs) .- 6*α))], U_profile.(xs, sqrt(σdot0^2 + 8*Φ0*π^2)-σdot0, λ/(σdot0+Φ0)) .+ σdot0 , linewidth = 2.5, color = my_black, label = "Approximation")
lines!(ax1, xs_att, speed_att, linewidth = 4, color = my_blue, label = "Attractive")

# Density kernels
density!(ax2, U_profile.(xs, sqrt(σdot0^2 + 8*Φ0*π^2)-σdot0, λ/(σdot0+Φ0)) .+ σdot0 , color = my_black, npoints = 200, boundary = (σdot0, sqrt(σdot0^2 + 8*Φ0*π^2)), label = "Approximation")
# density!(ax2, productlog_approx.(xs,σdot0, -Φ0), color = my_green, boundary = (σdot0, sqrt(σdot0^2 + 8*Φ0*π^2)), label = "ProductLog Sol")
density!(ax2, speed_att, color = my_blue, npoints = 200, boundary = (minimum(speed_att), maximum(speed_att)), label = "Numerics")


density!(ax3, U_profile.(xs, sqrt(σdot0^2 - 8*Φ0*π^2)-σdot0, λ/(σdot0-Φ0)) .+ σdot0 , color = my_black, npoints = 200, boundary = (sqrt(σdot0^2 - 8*Φ0*π^2), σdot0), label = "Approximation")
# density!(ax3, productlog_approx.(xs,σdot0, Φ0), color = my_green, boundary = (sqrt(σdot0^2 - 8*Φ0*π^2), σdot0), label = "ProductLog Sol")
density!(ax3, speed_rep, color = my_vermillion, npoints = 200, boundary = (minimum(speed_rep), maximum(speed_rep)), label = "Numerics")
# vlines!(ax2, [σdot0], color = my_black, linestyle = :dash, linewidth = 4)

# ylims!(ax2, -0.05, 1.0)
axislegend(ax1, labelsize = 30, position = :rt)
axislegend(ax2, labelsize = 24, position = :cc)
axislegend(ax3, labelsize = 24, position = :cc)
fig
