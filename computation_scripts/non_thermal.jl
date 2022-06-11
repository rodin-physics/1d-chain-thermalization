include("../src/main.jl")

## SINGLE PASS LOSS
α = 10
μ = 1
σ0 = 45
nChain = 10
vPts = 100
λs = [1 / 8, 1 / 4, 1 / 2, 1]
system_slow = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
system_fast = load_object("precomputed/systems/System_ωmax10_d6000_l10.jld2")
slow_res = Dict{Tuple{Float64,Float64},Tuple{Vector{Float64},Vector{Float64}}}()
fast_res = Dict{Tuple{Float64,Float64},Tuple{Vector{Float64},Vector{Float64}}}()

function Δ_numeric(σ_dot, σ0, Φ0, λ, system)
    δ = system.δ
    τ = 1.25 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads = true)

    σs = res.σs |> vec
    mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
    v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
    return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
end

# Low speed
Φ0s = 0.025
par = [(λ, Φ0) for λ in λs, Φ0 in [Φ0s, -Φ0s]] |> vec
σ_dots = range(2, 20, length = vPts)


if (!isfile("data/non_thermal/Single_Pass_Slow_Φ$(Φ0s)_μ$(μ).jld2"))

    for p in par
        λ = p[1]
        Φ0 = p[2]
        res = [Δ_numeric(x, σ0, Φ0, λ, system_slow) for x in σ_dots]
        merge!(slow_res, Dict((Φ0, λ) => (σ_dots, res)))
    end

    save_object("data/non_thermal/Single_Pass_Slow_Φ$(Φ0s)_μ$(μ).jld2", slow_res)
end

# High speed
Φ0s = 2.0
par = [(λ, Φ0) for λ in λs, Φ0 in [Φ0s, -Φ0s]] |> vec
σ_dots = range(20, 350, length = vPts)

if (!isfile("data/non_thermal/Single_Pass_Fast_Φ$(Φ0s)_μ$(μ).jld2"))

    for p in par
        λ = p[1]
        Φ0 = p[2]
        res = [Δ_numeric(x, σ0, Φ0, λ, system_fast) for x in σ_dots]
        merge!(fast_res, Dict((Φ0, λ) => (σ_dots, res)))
    end

    save_object("data/non_thermal/Single_Pass_Fast_Φ$(Φ0s)_μ$(μ).jld2", fast_res)
end




# d = 60
# σ_dots = range(20, 350, length = vPts)

# fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(
#     fig[1, 1],
#     xlabel = L"\dot{\sigma}_0",
#     ylabel = L"\Delta",
#     yscale = log10,
#     # xscale = log10,
#     # yminorticksvisible = true,
#     # yminorgridvisible = true,
#     # yminorticks = IntervalsBetween(10),
#     # xminorticksvisible = true,
#     # xminorgridvisible = true,
#     # xminorticks = IntervalsBetween(10),
# )


# for ii = 1:length(λs)

#     λ = λs[ii]
#     res_att = [Δ_numeric(x, σ0, -Φ0, λ, system_fast) for x in σ_dots]
#     res_rep = [Δ_numeric(x, σ0, Φ0, λ, system_fast) for x in σ_dots]
#     analytic = Δ_analytic.(σ_dots, Φ0, λ, system_fast.ωmax)

#     scatter!(ax1, (σ_dots), (res_att), color = colors[ii], marker = :hline, markersize = 20)
#     scatter!(ax1, (σ_dots), (res_rep), color = colors[ii], marker = :cross, markersize = 20)
#     lines!(
#         ax1,
#         (σ_dots),
#         (analytic),
#         color = colors[ii],
#         linewidth = 4,
#         label = L"\lambda = %$λ",
#     )
# end
# axislegend(ax1, position = :rt)
# fig

# colors = [my_vermillion, my_orange, my_green, my_sky, my_blue]


# function Δ_analytic(v, Φ, λ, Ω)
#     z = (2 * π * λ / v)^2
#     return (
#         4 * π^3 * Φ^2 / v^2 *
#         z *
#         exp(-z * (Ω^2 + 1) / 2) *
#         (
#             besseli(0, z * (Ω^2 - 1) / 2) +
#             (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
#         )
#     )
# end


# fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(
#     fig[1, 1],
#     xlabel = L"\dot{\sigma}_0",
#     ylabel = L"\Delta",
#     yscale = log10,
#     # xscale = log10,
#     # yminorticksvisible = true,
#     # yminorgridvisible = true,
#     # yminorticks = IntervalsBetween(10),
#     # xminorticksvisible = true,
#     # xminorgridvisible = true,
#     # xminorticks = IntervalsBetween(10),
# )


# for ii = 1:length(λs)

#     analytic = Δ_analytic.(σ_dots, Φ0, λ, system_slow.ωmax)

#     scatter!(ax1, (σ_dots), (res_att), color = colors[ii], marker = :hline, markersize = 20)
#     scatter!(ax1, (σ_dots), (res_rep), color = colors[ii], marker = :cross, markersize = 20)
#     lines!(
#         ax1,
#         (σ_dots),
#         (analytic),
#         color = colors[ii],
#         linewidth = 4,
#         label = L"\lambda = %$λ",
#     )
# end
# axislegend(ax1, position = :rb)
# save("Slow_dissipation.pdf", fig)


# # # GENERAL EXAMPLE
# system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
# d = 60
# τ = 200                             # Simulation time
# δ = system.δ                        # Time step
# α = 10                              # Distance between chain atoms
# μ = 1

# n_pts = τ / δ |> floor |> Int
# nChain = 200
# ρHs = zeros(nChain, n_pts)
# tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
# σdot0 = [25]
# mem = Inf

# σ0 = [55]
# Φ0s = [1 / 2, -1 / 2, -1, 1]
# # Φ0 = [1 / 2, 1, 2, -1 / 2, -1, -2]
# λs = [1 / 2, 1]
# params = [(f, l) for f in Φ0s, l in λs] |> vec
# println("Starting Calculations")
# for param in params
#     println(param)
#     Φ0 = param[1]
#     λ = param[2]
#     if (
#         !isfile(
#             "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
#         )
#     )

#         res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ)
#         save_object(
#             "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
#             res,
#         )
#     end
# end
# @showprogress pmap(params) do param
#     Φ0 = param[1]
#     λ = param[2]
#     if (
#         !isfile(
#             "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
#         )
#     )

#         res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ)
#         save_object(
#             "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
#             res,
#         )
#     end
# end


# ## Burst Energy Loss
# system = load_object("precomputed/systems/System_ωmax10_d3000_l10.jld2")

# d = 3000
# δ = system.δ
# α = 10
# Φ0 = 1/2
# σ0 = 45
# nChain = 10
# vPts = 100
# σ_dots = range(20, 350, length = vPts)
# Δs = zeros(vPts)


# @showprogress for ii = 1:vPts
#     σ_dot = σ_dots[ii]
#     τ = 1.1 * (α / σ_dot)
#     n_pts = τ / δ |> floor |> Int
#     ρHs = zeros(nChain, n_pts)
#     tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

#     res = motion_solver(system, Φ0, 1/4, α, [σ0], [σ_dot], 1, tTraj, Inf, τ)
#     σs = res.σs |> vec
#     mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
#     v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
#     Δs[ii] = μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2)
# end

# λ = 1/4
# Φ = 1/2
# function en_loss(v)
#     z = (2 * π * λ / v)^2
#     return (
#         4 * π^3 * Φ^2 / v^2 *
#         z *
#         exp(-z * (Ω^2 + 1) / 2) *
#         (
#             besseli(0, z * (Ω^2 - 1) / 2) +
#             (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
#         )
#     )
# end
# analytic = en_loss.(σ_dots)

# fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")

# # lines!(ax1,log.(σ_dots), log.(Δs))
# # # lines!(ax1,log.(σ_dots), log.(Δs_2))
# # lines!(ax1, log.(σ_dots), log.(analytic))
# # Δs- Δs_2


# lines!(ax1,(σ_dots), (Δs))
# # lines!(ax1,(σ_dots), (Δs_2))
# lines!(ax1, (σ_dots), (analytic))
# fig


# μs = [1 / 5, 1, 5]
# λs = [1 / 8, 1 / 4, 1 / 2, 1, 2]

# σ_dot = 25


# res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, Inf, τ)
# σs = res.σs |> vec
# v = (σs[2:end] - σs[1:end-1]) / (res.τs[2] - res.τs[1])
# scatter((res.τs)[2:end], v)


# mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
# v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
# Δ_en = μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2)
# μ
# σs[mid_pt_idx]


# # system = load_object("precomputed/systems/System_ωmax10_d3000_l10.jld2")

# # d = 3000
# # τ = .4                            # Simulation time
# # δ = system.δ                        # Time step
# # a = 10                               # Distance between chain atoms
# # n_pts = τ / δ |> floor |> Int
# # nChain = 10
# # ρHs = zeros(nChain, n_pts)
# # tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
# # λ = 1;
# # Φ = 1;
# # μ = 25/1;
# # @time res = motion_solver(system, Φ, λ, a, [25], [25], μ, tTraj, Inf, τ)
# # s = (res.σs |> vec)
# # v = (s[2:end] - s[1:end-1]) / (res.τs[2] - res.τs[1])
# # scatter((res.τs)[2:end], v)
# # kin_en = μ * v.^ 2 / 2 / (2 * π)^2
# # kin_en[end] - kin_en[1]
