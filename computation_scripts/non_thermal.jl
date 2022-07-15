include("../src/main.jl")

## SINGLE PASS LOSS
# α = 10
# μ = 1
# σ0 = 45
# nChain = 10
# vPts = 100
# λs = [1 / 8, 1 / 4, 1 / 2, 1]
# system_slow = load_object("precomputed/systems/System_ωmax10_d60_l300.jld2")
# system_fast = load_object("precomputed/systems/System_ωmax10_d6000_l10.jld2")
# slow_res = Dict{Tuple{Float64,Float64},Tuple{Vector{Float64},Vector{Float64}}}()
# fast_res = Dict{Tuple{Float64,Float64},Tuple{Vector{Float64},Vector{Float64}}}()
#
# function Δ_numeric(σ_dot, σ0, Φ0, λ, system)
#     δ = system.δ
#     τ = 1.25 * (α / σ_dot)
#     n_pts = τ / δ |> floor |> Int
#     ρHs = zeros(nChain, n_pts)
#     tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
#
#     res = motion_solver(system, Φ0, λ, α, [σ0], [σ_dot], μ, tTraj, Inf, τ, threads=true)
#
#     σs = res.σs |> vec
#     mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
#     v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
#     return (μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2))
# end
#
# # Low speed
# Φ0s = 0.025
# par = [(λ, Φ0) for λ in λs, Φ0 in [Φ0s, -Φ0s]] |> vec
# σ_dots = range(2, 20, length=vPts)
#
#
# if (!isfile("data/non_thermal/Single_Pass_Slow_Φ$(Φ0s)_μ$(μ).jld2"))
#
#     for p in par
#         λ = p[1]
#         Φ0 = p[2]
#         res = [Δ_numeric(x, σ0, Φ0, λ, system_slow) for x in σ_dots]
#         merge!(slow_res, Dict((Φ0, λ) => (σ_dots, res)))
#     end
#
#     save_object("data/non_thermal/Single_Pass_Slow_Φ$(Φ0s)_μ$(μ).jld2", slow_res)
# end
#
# # High speed
# Φ0s = 2.0
# par = [(λ, Φ0) for λ in λs, Φ0 in [Φ0s, -Φ0s]] |> vec
# σ_dots = range(20, 350, length=vPts)
#
# if (!isfile("data/non_thermal/Single_Pass_Fast_Φ$(Φ0s)_μ$(μ).jld2"))
#
#     for p in par
#         λ = p[1]
#         Φ0 = p[2]
#         res = [Δ_numeric(x, σ0, Φ0, λ, system_fast) for x in σ_dots]
#         merge!(fast_res, Dict((Φ0, λ) => (σ_dots, res)))
#     end
#
#     save_object("data/non_thermal/Single_Pass_Fast_Φ$(Φ0s)_μ$(μ).jld2", fast_res)
# end

## FULL TRAJECTORY
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
d = 60
τ = 100                             # Simulation time
δ = system.δ                        # Time step
α = 40                              # Distance between chain atoms
μ = 1
σ0 = [Int(5.5 * α)]

n_pts = τ / δ |> floor |> Int
nChain = 200
ρHs = zeros(nChain, n_pts)
tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
mem = Inf

params = [(8, 4, [50]), (-8, 4, [50])]
# bias = Δ_analytic(20, 8, 4, system.ωmax)
bias = 0.0

println("Starting Calculations")
for param in params
    println(param)
    Φ0 = param[1]
    λ = param[2]
    σdot0 = param[3]
    if (
        !isfile(
            "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_bias$(bias)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ, threads=true, bias=bias)
        save_object(
            "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_bias$(bias)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

## SUPERSONIC
# system = load_object("precomputed/systems/System_ωmax10_d1000_l300.jld2")
# d = 1000
# τ = 10                              # Simulation time
# δ = system.δ                        # Time step
# α = 10                              # Distance between chain atoms
# μ = 1
# σ0 = [55]
#
# n_pts = τ / δ |> floor |> Int
# nChain = 300
# ρHs = zeros(nChain, n_pts)
# tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
# mem = Inf
#
# params = [(1, 1, [300]), (-1, 1, [300]),
#     (5, 1, [300]), (-5, 1, [300]),
#     (5, 1, [600]), (-5, 1, [600]),
#     (10, 1, [600]), (-10, 1, [600]),
#     (15, 1, [600]), (-15, 1, [600]),
#     (20, 1, [600]), (-20, 1, [600])]
#
# println("Starting Supersonic Calculations")
# for param in params
#     println(param)
#     Φ0 = param[1]
#     λ = param[2]
#     σdot0 = param[3]
#     if (
#         !isfile(
#             "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
#         )
#     )
#
#         res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ, threads=true)
#         save_object(
#             "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
#             res,
#         )
#     end
# end
