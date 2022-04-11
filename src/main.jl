using CairoMakie
using Colors
using Distributions
using JLD2
using LaTeXStrings
using LinearAlgebra
using ProgressMeter
using QuadGK
using Statistics
using StatsBase
using ToeplitzMatrices

## Parameters
η = 1e-12                       # Small number
## Friendly colors
my_red = colorant"rgba(204, 121, 167, 0.75)"
my_vermillion = colorant"rgba(213, 94, 0, 0.75)"
my_orange = colorant"rgba(230, 159, 0, 0.75)"
my_yellow = colorant"rgba(240, 228, 66, 0.75)"
my_green = colorant"rgba(0, 158, 115, 0.75)"
my_sky = colorant"rgba(86, 180, 233, 0.75)"
my_blue = colorant"rgba(0, 114, 178, 0.75)"
my_black = colorant"rgba(0, 0, 0, 0.35)"

## Types
struct ChainSystem
    ωmax::Float64               # Largest mode frequency
    δ::Float64                  # Time step
    Γ::Matrix{Float64}          # Response array
end

struct ThermalTrajectory
    ωmax::Float64               # Largest mode frequency
    δ::Float64                  # Time step
    ρHs::Matrix{Float64}        # Thermal trajectory
    ωT::Union{Nothing,Float64}  # Temperature
end

struct SystemSolution
    ωmax::Float64               # Largest mode frequency
    μ::Float64                  # Mass of the mobile atoms
    τs::Vector{Float64}         # Time steps
    τ0::Float64                 # Memory
    Φ::Float64                  # Magnitude of the Gaussian potential
    λ::Float64                  # Standard deviation of the potential
    σs::Matrix{Float64}         # Positions of mobile atoms
    ρs::Matrix{Float64}         # Positions of the chain atoms
    ωT::Union{Nothing,Float64}  # Temperature
end

## Functions
# Frequency as a function of momentum
@inline function Ω(ωmax, x)
    return sqrt(1 + (ωmax^2 - 1) * sin(x)^2)
end

# Chain recoil function
function Γ(τ, ls, ωmax)
    int_fun(x) = cos.(2 * x * ls) * sin(2 * π * τ * Ω(ωmax, x)) / Ω(ωmax, x)
    res = quadgk(int_fun, 0, π / 2)
    return (res[1] * 2 / π)
end

function mkChainSystem(ωmax, τ_max, ls, d)
    δ = (1 / ωmax) / d                  # Time step in units of t_M
    n_pts = floor(τ_max / δ) |> Int     # Number of time steps given τ_max and δ
    Γ_mat = zeros(ls, n_pts)
    @showprogress Threads.@threads for ii = 1:n_pts
        Γ_mat[:, ii] = Γ(δ * ii, ls, ωmax)
    end
    return ChainSystem(ωmax, δ, Γ_mat)
end

# # Mode amplitude
# function ζq(ωq, ωT, μ)
#     # Subtract a small number from p. The reason is that for low ωT, p ≈ 1,
#     # causing issues with the rand() generator
#     n = rand(Geometric(1 - exp(-ωq / ωT) - η))
#     res = √(n + 1 / 2) * √(2 / ωq)
#     return res
# end

# # Homogeneous displacement of the active chain atom at time step n given a set of ωs
# # and the corresponding ζs and ϕs as a sum of normal coordinates. 
# function ζH(n, δ, ζs, ϕs, ωs)
#     n_ζ = length(ζs)
#     res = [ζs[x] * cos(2 * π * δ * n * ωs[x] + ϕs[x]) / √(n_ζ) for x = 1:n_ζ] |> sum
#     return res
# end

# Function that calculates the trajectories of the mobile atoms and the chain particle.
# mem determines the memory length and τ is the simulation period, both in units of
# the trap period.
function motion_solver(
    system::ChainSystem,
    Φ0::T where {T<:Real},
    λ::T where {T<:Real},
    a::T where {T<:Real},
    σ0::Vector{T} where {T<:Real},
    σ_dot0::Vector{T} where {T<:Real},
    μ::T where {T<:Real},
    tTraj::ThermalTrajectory,
    τ0::T where {T<:Real},
    τ::T where {T<:Real},
)
    δ = system.δ        # Array of time steps
    Γ_mat = system.Γ    # Memory term
    # Check that the thermal trajectory is for the correct system
    if (ωmax != tTraj.ωmax || δ != tTraj.δ)
        error("Thermal trajectory describes a different system. Check your input.")
    else
        ρHs = tTraj.ρHs
    end

    n_pts = floor(τ / δ) |> Int         # Number of time steps
    τs = δ .* (1:n_pts) |> collect      # Times
    if size(ρHs)[2] < n_pts
        error("Thermal trajectory provided does not span the necessary time range.")
    end

    τ0_pts = max(floor(τ0 / δ), 1)    # Memory time points.
    # Even if τ0 == 0, we have to keep a single time point to make sure arrays work.
    # For zero memory, the recoil contribution is dropped (see line 156)

    # If the precomputed memory is shorter than the simulation time AND shorter
    # than the desired memory, terminate the calculation.
    if (size(Γ_mat)[2] < n_pts && size(Γ_mat)[2] < τ0_pts)
        error("Chosen memory and the simulation time exceed the precomputed range.")
    else
        # The number of memory pts can be limited by the total simulation time.
        τ0_pts = min(τ0_pts, n_pts) |> Int
        Γ_mat = (2 * π * δ) .* Γ_mat[:, 1:τ0_pts]
    end

    # Interaction terms
    function dU_dρ(r)
        return (-Φ0 * exp(-r^2 / (2 * λ^2)) * r / λ^2)
    end

    ## Arrays

    σs = zeros(n_pts, length(σ0))           # Mobile mass position
    ρs = ρHs[1:n_pts]                       # Chain mass position

    # Initial values
    σs[1, :] = σ0
    σs[2, :] = σ0

    curr = 1
    nxt = curr + 2
    U_pr = dU_dρ.(ρs[curr] .- σs[curr, :])
    U_pr_chain = sum(U_pr)
    U_pr_mob = -U_pr
    ρs[nxt:nxt+min(τ0_pts - 2, n_pts - nxt)] -=
        Γ_list[2:min(τ0_pts, n_pts - curr)] .* (U_pr_chain * (τ0 != 0))
    @showprogress for ii = 3:n_pts
        nxt = ii        # Next time step index
        curr = ii - 1   # Current time step index
        U_pr = dU_dρ.(ρs[curr] .- σs[curr, :])
        U_pr_chain = sum(U_pr)
        U_pr_mob = -U_pr
        ρs[nxt:nxt+min(τ0_pts - 1, n_pts - nxt)] -=
            Γ_list[1:min(τ0_pts, n_pts - curr)] .* (U_pr_chain * (τ0 != 0))
        σs[nxt, :] =
            -(2 * π * δ)^2 .* (U_pr_mob + σs[curr, :]) + 2 .* σs[curr, :] - σs[curr-1, :]
    end
    return SystemSolution(ωmin, ωmax, μ, τs, τ0, Φ0, λ, σs, ρs, tTraj.ωT)
end















#     nChain = length(rHs[1])

#     Ωmin = √(K / m)                     # Minimum phonon frequency
#     tmin = 2 * π / Ωmin                 # Period of the slowest mode
#     n_pts = floor(τ * tmin / δ) |> Int  # Number of time steps
#     ts = δ .* (1:n_pts) |> collect      # Times
#     if length(rHs) < n_pts
#         error("Thermal trajectory provided does not span the necessary time range.")
#     end

#     mem_pts = max(floor(mem * tmin / δ), 1)  # Memory time points.
#     # Even if mem == 0, we have to keep a single time point to make sure arrays work.
#     # For zero memory, the recoil contribution is dropped (see line 175 and 184)

#     # If the precomputed memory is shorter than the simulation time AND shorter
#     # than the desired memory, terminate the calculation.
#     if (length(G_list) < n_pts && length(G_list) < mem_pts)
#         error("Chosen memory and the simulation time exceed the precomputed range.")
#     else
#         # The number of memory pts can be limited by the total simulation time.
#         mem_pts = min(mem_pts, n_pts) |> Int
#         G_list = δ .* G_list[1:mem_pts]
#     end

#     # If the desired number of chain particles is greater than what is contained in
#     # the precomputed G, terminate the calculation. Otherwise, retain the appropriate
#     # number of terms
#     if (length(G_list[1]) < nChain)
#         error("The recoil term does not contain the desired number of chain masses.")
#     else
#         G_list = [SymmetricToeplitz(x[1:nChain]) for x in G_list]
#     end

#     function dU_dr(r)
#         return (-F * exp(-r^2 / (2 * s^2)) * r / s^2)
#     end
#     ## Arrays
#     Rs = [zeros(length(x0)) for _ = 1:n_pts]# Mobile mass position
#     rs = rHs[1:n_pts]                       # Chain mass position

#     Rs[1] = x0
#     Rs[2] = x0 + v0 * δ

#     curr = 1
#     nxt = curr + 2
#     U_pr = [dU_dr(r - R) for r in rs[curr], R in Rs[curr]]
#     U_pr_chain = sum(U_pr, dims = 2) |> vec
#     U_pr_mob = -sum(U_pr, dims = 1) |> vec
#     rs[nxt:nxt+min(mem_pts - 2, n_pts - nxt)] -=
#         G_list[2:min(mem_pts, n_pts - curr)] .* Ref(U_pr_chain .* (mem != 0))

#     @showprogress for ii = 3:n_pts
#         nxt = ii        # Next time step index
#         curr = ii - 1   # Current time step index
#         U_pr = [dU_dr(r - R) for r in rs[curr], R in Rs[curr]]
#         U_pr_chain = sum(U_pr, dims = 2) |> vec
#         U_pr_mob = -sum(U_pr, dims = 1) |> vec
#         rs[nxt:nxt+min(mem_pts - 1, n_pts - nxt)] -=
#             G_list[1:min(mem_pts, n_pts - curr)] .* Ref(U_pr_chain .* (mem != 0))
#         Rs[nxt] = -δ^2 / M .* U_pr_mob + 2 .* Rs[curr] - Rs[curr-1]
#     end
#     return SystemSolution(k, K, m, ts, mem, tTraj.a, M, F, s, Rs, rs, tTraj.ΩT, tTraj.ħ)
# end

function motion_solver_NEW(
    system::ChainSystem,
    F::T where {T<:Real},
    s::T where {T<:Real},
    M::T where {T<:Real},
    x0::Vector{T} where {T<:Real},
    v0::Vector{T} where {T<:Real},
    tTraj::ThermalTrajectory,
    mem::T where {T<:Real},
    τ::T where {T<:Real},
)
    m = system.m        # Mass of the chain atoms
    k = system.k        # Spring force constant
    K = system.K        # Confining potential force constant
    δ = system.δ        # Array of time steps
    G_list = system.G   # Memory term
    # Check that the thermal trajectory is for the correct system
    if (m != tTraj.m || k != tTraj.k || K != tTraj.K || δ != tTraj.δ)
        error("Thermal trajectory describes a different system. Check your input.")
    else
        rHs = tTraj.rHs
    end

    nChain = length(rHs[1])

    Ωmin = √(K / m)                     # Minimum phonon frequency
    tmin = 2 * π / Ωmin                 # Period of the slowest mode
    n_pts = floor(τ * tmin / δ) |> Int  # Number of time steps
    ts = δ .* (1:n_pts) |> collect      # Times
    if length(rHs) < n_pts
        error("Thermal trajectory provided does not span the necessary time range.")
    end

    mem_pts = max(floor(mem * tmin / δ), 1)  # Memory time points.
    # Even if mem == 0, we have to keep a single time point to make sure arrays work.
    # For zero memory, the recoil contribution is dropped (see line 175 and 184)

    # If the precomputed memory is shorter than the simulation time AND shorter
    # than the desired memory, terminate the calculation.
    if (length(G_list) < n_pts && length(G_list) < mem_pts)
        error("Chosen memory and the simulation time exceed the precomputed range.")
    else
        # The number of memory pts can be limited by the total simulation time.
        mem_pts = min(mem_pts, n_pts) |> Int
        G_list = δ .* G_list[1:mem_pts]
    end

    # If the desired number of chain particles is greater than what is contained in
    # the precomputed G, terminate the calculation. Otherwise, retain the appropriate
    # number of terms
    if (length(G_list[1]) < nChain)
        error("The recoil term does not contain the desired number of chain masses.")
    else
        G_list = [SymmetricToeplitz(x[1:nChain]) for x in G_list]
    end

    function dU_dr(r)
        return (-F * exp(-r^2 / (2 * s^2)) * r / s^2)
    end
    ## Arrays
    Rs = [zeros(length(x0)) for _ = 1:n_pts]# Mobile mass position
    rs = rHs[1:n_pts]                       # Chain mass position

    Rs[1] = x0
    Rs[2] = x0 + v0 * δ

    curr = 1
    nxt = curr + 2
    U_pr = [dU_dr(r - R) for r in rs[curr], R in Rs[curr]]
    U_pr_chain = sum(U_pr, dims = 2) |> vec
    U_pr_mob = -sum(U_pr, dims = 1) |> vec
    rs[nxt:nxt+min(mem_pts - 2, n_pts - nxt)] -=
        G_list[2:min(mem_pts, n_pts - curr)] .* Ref(U_pr_chain .* (mem != 0))

    @showprogress for ii = 3:n_pts
        nxt = ii        # Next time step index
        curr = ii - 1   # Current time step index
        U_pr = [dU_dr(r - R) for r in rs[curr], R in Rs[curr]]
        U_pr_chain = sum(U_pr, dims = 2) |> vec
        U_pr_mob = -sum(U_pr, dims = 1) |> vec
        rs[nxt:nxt+min(mem_pts - 1, n_pts - nxt)] -=
            G_list[1:min(mem_pts, n_pts - curr)] .* Ref(U_pr_chain .* (mem != 0))
        Rs[nxt] = -δ^2 / M .* U_pr_mob + 2 .* Rs[curr] - Rs[curr-1]
    end
    return SystemSolution(k, K, m, ts, mem, tTraj.a, M, F, s, Rs, rs, tTraj.ΩT, tTraj.ħ)
end