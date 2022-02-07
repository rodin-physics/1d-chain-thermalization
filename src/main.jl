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
using Symbolics
## Parameters
η = 1e-12                       # Small number for integration
ħ = 1 / 200                     # Planck's constant
my_red = colorant"rgba(215, 67, 84, 0.75)"
my_green = colorant"rgba(106, 178, 71, 0.75)"
my_blue = colorant"rgba(100, 101, 218, 0.75)"
my_violet = colorant"rgba(169, 89, 201, 0.75)"
my_orange = colorant"rgba(209, 135, 46, 0.75)"

## Types
struct ChainSystem
    k::Float64                  # Spring force constant
    K::Float64                  # Confining force constant
    m::Float64                  # Mass of the chain atoms
    δ::Float64                  # Time step
    G::Vector{Vector{Float64}}  # Response array. External vector runs over time,
    # internal over mass separation
end

struct ThermalTrajectory
    k::Float64                  # Spring force constant
    K::Float64                  # Confining force constant
    m::Float64                  # Mass of the chain atoms
    δ::Float64                  # Time step
    rHs::Vector{Vector{Float64}}# Thermal Trajectory. External vector runs over time,
    # internal over masses
    ΩT::Union{Nothing,Float64}  # Temperature
    ħ::Float64                  # Planck's constant
end

# struct SystemSolution
#     k::Float64                  # Spring force constant
#     K::Float64                  # Confining force constant
#     m::Float64                  # Mass of the chain atoms
#     ts::Vector{Float64}         # Time steps
#     mem::Float64                # Response array
#     K_M::Float64                # Confining potential
#     M::Float64                  # Mass of the mobile atoms
#     F::Float64                  # Magnitude of the Gaussian potential
#     s::Float64                  # Standard deviation of the potential
#     Rs::Matrix{Float64}         # Positions of the mobile atoms
#     rs::Vector{Float64}         # Positions of the chain atom
#     rHs::Vector{Float64}        # Homogeneous position of the chain atom
#     U_pr_mob::Matrix{Float64}   # Derivatives of the potential for mobile atoms
#     U_pr_chain::Vector{Float64} # Derivative of the potential for the chain atom
#     ΩT::Union{Nothing,Float64}  # Temperature
#     ħ::Float64                  # Planck's constant
# end

## Functions
# Frequency as a function of momentum
@inline function Ω(K, k, m, x)
    return sqrt(4 * k / m * sin(x)^2 + K / m)
end

# Chain recoil function
function G(t, l_max, K, k, m)
    ls = 0:l_max
    int_fun(x) = cos.(2 * x * ls) * sin(t * Ω(K, k, m, x)) / Ω(K, k, m, x)
    res = quadgk(int_fun, 0, π / 2, atol = 1e-5)
    return (res[1] * 2 / π / m)
end

# Function for assembiling a ChainSystem
function mkChainSystem(K, k, m, t_max, l_max, d)
    Ωmax = sqrt(4 * k / m + K / m)      # Largest chain frequency
    δ = (2 * π / Ωmax) / d              # Time step
    n_pts = floor(t_max / δ) |> Int     # Number of time steps given t_max and δ
    # Precomputing the memory term
    G_list = @showprogress pmap(jj -> G(δ * jj, l_max, K, k, m), 1:n_pts) |> reverse
    # Note the reversal of G_list. This is done to facilitate the
    # convolution of G with dU/dr.
    return ChainSystem(k, K, m, δ, G_list)
end

# # Mode amplitude
# function ζq(Ωq, ΩT, ħ)
#     # Subtract a small number from p. The reason is that for low ΩT, p ≈ 1,
#     # causing issues with the rand() generator
#     n = rand(Geometric(1 - exp(-Ωq / ΩT) - η))
#     res = √(n + 1 / 2) * √(2 * ħ / Ωq)
#     return res
# end

# # Homogeneous displacement of the active chain atom at time step n given a set of Ωs
# # and the corresponding ζs and ϕs as a sum of normal coordinates. 
# function ζH(n, δ, ζs, ϕs, Ωs)
#     n_ζ = length(ζs)
#     res = [ζs[x] * cos(δ * n * Ωs[x] + ϕs[x]) / √(n_ζ) for x = 1:n_ζ] |> sum
#     return res
# end

# # Function that calculates the trajectories of the mobile atoms and the chain particle.
# # mem determines the memory length and τ is the simulation period, both in units of
# # the trap period.
# function motion_solver(
#     system::ChainSystem,
#     F::T where {T<:Real},
#     s::T where {T<:Real},
#     M::T where {T<:Real},
#     K_M::T where {T<:Real},
#     x0::Vector{T} where {T<:Real},
#     tTraj::ThermalTrajectory,
#     mem::T where {T<:Real},
#     τ::T where {T<:Real},
# )
#     m = system.m            # Mass of the chain atoms
#     k = system.k            # Spring force constant
#     K = system.K            # Confining potential force constant
#     δ = system.δ            # Array of time steps
#     G_list = system.G       # Memory term
#     # Check that the thermal trajectory is for the correct system
#     if (m != tTraj.m || k != tTraj.k || K != tTraj.K || δ != tTraj.δ)
#         error("Thermal trajectory describes a different system. Check your input.")
#     else
#         rHs = tTraj.rHs
#     end

#     ΩM = √(K_M / M)                     # Trap frequency
#     t_M = 2 * π / ΩM                    # Period of the trapped mass
#     n_pts = floor(τ * t_M / δ) |> Int   # Number of time steps
#     ts = δ .* (1:n_pts) |> collect      # Times
#     if length(rHs) < n_pts
#         error("Thermal trajectory provided does not span the necessary time range.")
#     end

#     mem_pts = max(floor(mem * t_M / δ), 1)  # Memory time points.
#     # Even if mem == 0, we have to keep a single time point to make sure arrays work.
#     # For zero memory, the recoil contribution is dropped (see line 208)

#     # If the precomputed memory is shorter than the simulation time AND shorter
#     # than the desired memory, terminate the calculation.
#     if (length(G_list) < n_pts && length(G_list) < mem_pts)
#         error("Chosen memory and the simulation time exceed the precomputed range.")
#     else
#         # The number of memory pts can be limited by the total simulation time.
#         mem_pts = min(mem_pts, n_pts) |> Int
#         G_list_ = G_list[(end-mem_pts+1):end]
#     end

#     # Interaction terms
#     @variables r R

#     Dr = Differential(r)
#     DR = Differential(R)

#     function U(r, R)
#         return (F * exp(-(r - R)^2 / (2 * s^2)))
#     end

#     der_r = expand_derivatives(Dr(U(r, R)))
#     der_R = expand_derivatives(DR(U(r, R)))

#     function dU_dr(r_, R_)
#         return (substitute.(der_r, (Dict(r => r_, R => R_),))[1])
#     end
#     function dU_dR(r_, R_)
#         return (substitute.(der_R, (Dict(r => r_, R => R_),))[1])
#     end

#     ## Arrays
#     Rs = zeros(n_pts, length(x0))         # Mobile mass position
#     rs = zeros(n_pts)                     # Chain mass position
#     U_pr_mob = zeros(n_pts, length(x0))   # Force on the mobile particles
#     U_pr_chain = zeros(n_pts)             # Force on the chain atom

#     ## Initial values
#     Rs[1, :] = x0              # Initialize the starting R
#     Rs[2, :] = x0              # Starting at rest
#     # Initialize the interaction terms
#     U_pr_mob[1, :] = Symbolics.value.(dU_dR.(rHs[1], x0))
#     U_pr_mob[2, :] = Symbolics.value.(dU_dR.(rHs[2], x0))
#     U_pr_chain[1] = sum(Symbolics.value.(dU_dr.(rHs[1], x0)))
#     U_pr_chain[2] = sum(Symbolics.value.(dU_dr.(rHs[2], x0)))
#     # Initialize the positions of the mobile particle
#     rs[1] = rHs[1]
#     rs[2] = rHs[2]

#     @showprogress for ii = 3:n_pts
#         nxt = ii        # Next time step index
#         curr = ii - 1   # Current time step index
#         # Get the memory term elements. If the number of elapsed steps is smaller
#         # than the number of elements in the memory term, take the corresponding
#         # number of the most recent elements
#         Gs = G_list_[(mem_pts-min(mem_pts, curr)+1):end]
#         # Get the interaction force elements. If the number of elapsed steps is
#         # smaller than the number of elements in the memory term, take the available
#         # values.
#         Us = U_pr_chain[(curr-min(mem_pts, curr)+1):curr]
#         # If mem == 0, drop the recoil contribution
#         rs[nxt] = rHs[nxt] - δ * dot(Gs, Us) * (mem != 0)
#         Rs[nxt, :] =
#             δ^2 / M .* (-U_pr_mob[curr, :] - K_M .* Rs[curr, :]) + 2 .* Rs[curr, :] -
#             Rs[curr-1, :]
#         U_pr_chain[nxt] = sum(Symbolics.value.(dU_dr.(rs[nxt], Rs[nxt, :])))
#         U_pr_mob[nxt, :] = Symbolics.value.(dU_dR.(rs[nxt], Rs[nxt, :]))
#     end

#     return SystemSolution(
#         k,
#         K,
#         m,
#         ts,
#         mem,
#         K_M,
#         M,
#         F,
#         s,
#         Rs,
#         rs,
#         rHs,
#         U_pr_mob,
#         U_pr_chain,
#         tTraj.ΩT,
#         tTraj.ħ,
#     )
# end

# function auto_corr(l, lag)
#     res = cor(l[1:(end-lag), :], l[(1+lag):end, :])
#     return res
# end



1000 * 2 * pi / Ω(1, 20, 1, 0) / (2 * pi / Ω(1, 20, 1, pi / 2) / 60)



1000 * 2 * pi / Ω(1, 20, 1, 0 * pi / 2)


Ω(1, 20, 1, pi / 2) / Ω(1, 20, 1, 0)
1000 * 2 * pi / sqrt(1 / 5) / (2 * pi / Ω(1, 20, 1 / 1, pi / 2) / 60)
2 * pi / Ω(1, 20, 1 / 1, pi / 2) / 60


Ω(1 / 20, 1, 1 / 1, pi / 2)



r = G(3000, 6000, 1 / 2000, 1, 2)

lines(0:6000, r)
2 * pi / Ω(1 / 20, 1, 1 / 2, 1 * pi / 2) / 60
500 * 2 * pi / Ω(1 / 20, 1, 1 / 2, 1 * pi / 2)
@time map(x -> G(15000, x, 1, 20, 1), collect(0:500))
1
