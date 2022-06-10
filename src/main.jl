using CairoMakie
using Colors
using Distributions
using JLD2
using LaTeXStrings
using LinearAlgebra
using ProgressMeter
using QuadGK
using SparseArrays
using SpecialFunctions
using Statistics
using StatsBase
using ToeplitzMatrices

## Parameters
η = 1e-12                       # Small number
ϵ = 1e-20
## Friendly colors
my_red = colorant"rgba(204, 121, 167, 0.75)"
my_vermillion = colorant"rgba(213, 94, 0, 0.75)"
my_orange = colorant"rgba(230, 159, 0, 0.75)"
my_yellow = colorant"rgba(240, 228, 66, 0.75)"
my_green = colorant"rgba(0, 158, 115, 0.75)"
my_sky = colorant"rgba(86, 180, 233, 0.75)"
my_blue = colorant"rgba(0, 114, 178, 0.75)"
my_black = colorant"rgba(0, 0, 0, 0.75)"

## Types
struct ChainSystem
    ωmax::Float64               # Largest mode frequency
    δ::Float64                  # Time step
    Γ::Matrix{Float64}          # Response array
end

struct Impulse
    n::Int                      # Index of the chain mass
    τ::Int                      # Time step index
    p::Float64                  # Impulse magnitude
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
    α::Float64                 # Spacing between chain atoms
    Φ::Float64                  # Magnitude of the Gaussian potential
    λ::Float64                  # Standard deviation of the potential
    σs::Matrix{Float64}         # Positions of mobile atoms
    ρs::Matrix{Float64}         # Positions of the chain atoms
    ωT::Union{Nothing,Float64}  # Temperature
end

## Functions
# Frequency as a function of momentum
@inline function ω(ωmax, x)
    return sqrt(1 + (ωmax^2 - 1) * sin(x)^2)
end

# Chain recoil function
function Γ(τ, ls, ωmax)
    int_fun(x) = cos.(2 * x * ls) * sin(2 * π * τ * ω(ωmax, x)) / ω(ωmax, x)
    res = quadgk(int_fun, 0, π / 2)
    return (res[1] * 2 / π)
end

function mkChainSystem(ωmax, τ_max, lmax, d)
    δ = (1 / ωmax) / d                          # Time step in units of t_slow
    n_pts = floor(τ_max / δ) |> Int             # Number of time steps given τ_max and δ
    Γ_mat = zeros(lmax + 1, n_pts)              # Prepare the Γ matrix
    pr = Progress(n_pts)                        # Setup the progress meter
    pts_per_thread = n_pts ÷ Threads.nthreads() # Calculate points per thread

    # Reorder the time steps so that the load is equal for every thread because
    # latter times take longer to evaluate
    thread_pts =
        reshape(1:(Threads.nthreads()*pts_per_thread), Threads.nthreads(), :)' |>
        vec |>
        collect

    Threads.@threads for ii in thread_pts
        Γ_mat[:, ii] = Γ(δ * ii, 0:lmax, ωmax)
        next!(pr)
    end

    # Evaluate the remaining time steps
    if thread_pts != n_pts
        Threads.@threads for ii in (length(thread_pts)+1:n_pts)
            Γ_mat[:, ii] = Γ(δ * ii, 0:lmax, ωmax)
            next!(pr)
        end
    end

    return ChainSystem(ωmax, δ, Γ_mat)
end

# # # # Mode amplitude
# # # function ζq(ωq, ωT, μ)
# # #     # Subtract a small number from p. The reason is that for low ωT, p ≈ 1,
# # #     # causing issues with the rand() generator
# # #     n = rand(Geometric(1 - exp(-ωq / ωT) - η))
# # #     res = √(n + 1 / 2) * √(2 / ωq)
# # #     return res
# # # end

# # # # Homogeneous displacement of the active chain atom at time step n given a set of ωs
# # # # and the corresponding ζs and ϕs as a sum of normal coordinates. 
# # # function ζH(n, δ, ζs, ϕs, ωs)
# # #     n_ζ = length(ζs)
# # #     res = [ζs[x] * cos(2 * π * δ * n * ωs[x] + ϕs[x]) / √(n_ζ) for x = 1:n_ζ] |> sum
# # #     return res
# # # end

# function motion_solver(
#     system::ChainSystem,
#     Φ0::T where {T<:Real},
#     λ::T where {T<:Real},
#     α::T where {T<:Real},
#     σ0::Vector{T} where {T<:Real},
#     σ_dot0::Vector{T} where {T<:Real},
#     μ::T where {T<:Real},
#     tTraj::ThermalTrajectory,
#     τ0::T where {T<:Real},
#     τ::T where {T<:Real},
# )
#     ωmax = system.ωmax              # Maximum chain frequency
#     δ = system.δ                    # Time step
#     Γ_mat = system.Γ                # Memory term
#     n_pts = floor(τ / δ) |> Int     # Number of time steps
#     τs = δ .* (1:n_pts) |> collect  # Times
#     nChain = size(tTraj.ρHs)[1]     # Number of chain particles for which the homogeneous motion is available
#     τ0_pts = max(floor(τ0 / δ), 1)  # Memory time points
#     # Even if τ0 == 0, we have to keep a single time point to make sure arrays work.
#     # For zero memory, the recoil contribution is dropped (see line 156)

#     # Check that the thermal trajectory is for the correct system
#     if (ωmax != tTraj.ωmax || δ != tTraj.δ)
#         error("Thermal trajectory describes a different system. Check your input.")
#     elseif size(tTraj.ρHs)[2] < n_pts
#         error("Thermal trajectory provided does not span the necessary time range.")
#         # If the precomputed memory is shorter than the simulation time AND shorter
#         # than the desired memory, terminate the calculation.
#     elseif (size(Γ_mat)[2] < n_pts && size(Γ_mat)[2] < τ0_pts)
#         error("Chosen memory and the simulation time exceed the precomputed recoil.")
#         # If the desired number of chain particles is greater than what is contained in
#         # the precomputed Γ, terminate the calculation. Otherwise, retain the appropriate
#         # number of terms
#     elseif (size(Γ_mat)[1] < nChain)
#         error("The recoil term does not contain the desired number of chain masses.")
#     else
#         ρs = (tTraj.ρHs)[:, 1:n_pts] .+ α .* repeat(1:size(tTraj.ρHs)[1], 1, n_pts)
#         σs = zeros(length(σ0), n_pts)
#         τ0_pts = min(τ0_pts, n_pts) |> Int
#         Γ_mat = (2 * π * δ) .* Γ_mat[1:nChain, 1:τ0_pts]
#         Γ_mat = vcat(reverse(Γ_mat, dims=1)[1:end-1, :], Γ_mat)
#     end
#     # Interaction terms
#     @inline function dU_dρ(r)
#         return (-Φ0 * exp(-r^2 / (2 * λ^2)) * r / λ^2)
#     end
#     ## Initial values
#     σs[:, 1] = σ0
#     σs[:, 2] = σ0 + δ .* σ_dot0

#     @showprogress for ii = 3:n_pts
#         nxt = ii        # Next time step index
#         curr = ii - 1   # Current time step index
#         # Calculate the forces on all the masses
#         U_pr = [dU_dρ(ρ - σ) for ρ in ρs[:, curr], σ in σs[:, curr]]
#         U_pr_chain = sum(U_pr, dims=2) |> vec
#         U_pr_mob = -sum(U_pr, dims=1) |> vec
#         # Find the indices of the chain masses where the force is larger than ϵ
#         idx = findall(x -> abs(x) > ϵ, U_pr_chain)
#         # For each of the impulses, get the appropriate slice of Γ_mat, multiply
#         # it by the impulse and use the result to modify the chain positions
#         for n in idx
#             Γ_curr = view(Γ_mat, nChain-n+1:2*nChain-n, 1:min(τ0_pts, n_pts - ii + 1))
#             ρs_upd = view(ρs, :, nxt:nxt+size(Γ_curr)[2]-1)
#             ρs_upd .-= Γ_curr .* U_pr_chain[n]
#         end
#         σs[:, nxt] = -(2 * π * δ)^2 / μ .* U_pr_mob + 2 .* σs[:, curr] - σs[:, curr-1]
#     end
#     return SystemSolution(ωmax, μ, τs, τ0, α, Φ0, λ, σs, ρs, tTraj.ωT)
# end
function motion_solver(
    system::ChainSystem,
    Φ0::T where {T<:Real},
    λ::T where {T<:Real},
    α::T where {T<:Real},
    σ0::Vector{T} where {T<:Real},
    σ_dot0::Vector{T} where {T<:Real},
    μ::T where {T<:Real},
    tTraj::ThermalTrajectory,
    τ0::T where {T<:Real},
    τ::T where {T<:Real},
)
    ωmax = system.ωmax              # Maximum chain frequency
    δ = system.δ                    # Time step
    Γ_mat = system.Γ                # Memory term
    n_pts = floor(τ / δ) |> Int     # Number of time steps
    τs = δ .* (1:n_pts) |> collect  # Times
    nChain = size(tTraj.ρHs)[1]     # Number of chain particles for which the homogeneous motion is available
    τ0_pts = max(floor(τ0 / δ), 1)  # Memory time points
    # Even if τ0 == 0, we have to keep a single time point to make sure arrays work.
    # For zero memory, the recoil contribution is dropped (see line 156)

    # Check that the thermal trajectory is for the correct system
    if (ωmax != tTraj.ωmax || δ != tTraj.δ)
        error("Thermal trajectory describes a different system. Check your input.")
    elseif size(tTraj.ρHs)[2] < n_pts
        error("Thermal trajectory provided does not span the necessary time range.")
        # If the precomputed memory is shorter than the simulation time AND shorter
        # than the desired memory, terminate the calculation.
    elseif (size(Γ_mat)[2] < n_pts && size(Γ_mat)[2] < τ0_pts)
        error("Chosen memory and the simulation time exceed the precomputed recoil.")
        # If the desired number of chain particles is greater than what is contained in
        # the precomputed Γ, terminate the calculation. Otherwise, retain the appropriate
        # number of terms
    elseif (size(Γ_mat)[1] < nChain)
        error("The recoil term does not contain the desired number of chain masses.")
    else
        ρs = (tTraj.ρHs)[:, 1:n_pts] .+ α .* repeat(1:size(tTraj.ρHs)[1], 1, n_pts)
        σs = zeros(length(σ0), n_pts)
        τ0_pts = min(τ0_pts, n_pts) |> Int
        Γ_mat = (2 * π * δ) .* Γ_mat[1:nChain, 1:τ0_pts]
        Γ_mat = vcat(reverse(Γ_mat, dims=1)[1:end-1, :], Γ_mat)

    end
    # Interaction terms
    @inline function dU_dρ(r)
        return (-Φ0 * exp(-r^2 / (2 * λ^2)) * r / λ^2)
    end
    ## Initial values
    σs[:, 1] = σ0
    σs[:, 2] = σ0 + δ .* σ_dot0

    @showprogress for ii = 3:n_pts
        nxt = ii        # Next time step index
        curr = ii - 1   # Current time step index
        # Calculate the forces on all the masses
        U_pr = [dU_dρ(ρ - σ) for ρ in ρs[:, curr], σ in σs[:, curr]]
        U_pr_chain = sum(U_pr, dims=2) |> vec
        U_pr_mob = -sum(U_pr, dims=1) |> vec
        # Find the indices of the chain masses where the force is larger than ϵ
        idx = findall(x -> abs(x) > ϵ, U_pr_chain)
        # For each of the impulses, get the appropriate slice of Γ_mat, multiply
        # it by the impulse and use the result to modify the chain positions
        steps_left = n_pts - curr
        steps_per_thread = steps_left ÷ Threads.nthreads()
        step_alloc = [(x*(steps_per_thread+1)+1):min(steps_left, (x + 1) * (steps_per_thread + 1))
                      for x in 0:Threads.nthreads()-1]
        Threads.@threads for t in 1:Threads.nthreads()
            for n in idx
                view(ρs, :, curr .+ step_alloc[t]) .-= view(Γ_mat, nChain-n+1:2*nChain-n, step_alloc[t]) .* U_pr_chain[n]
            end
        end
        σs[:, nxt] = -(2 * π * δ)^2 / μ .* U_pr_mob + 2 .* σs[:, curr] - σs[:, curr-1]
    end



    
    return SystemSolution(ωmax, μ, τs, τ0, α, Φ0, λ, σs, ρs, tTraj.ωT)
end