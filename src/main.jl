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
const ħ = 1 / 200                     # Planck's constant
const ϵ = 1e-150                      # Threshold to determine dUdr truncation
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
    a::Float64                  # Chain mass spacing
    δ::Float64                  # Time step
    rHs::Vector{Vector{Float64}}# Thermal Trajectory. External vector runs over time,
    # internal over masses
    ΩT::Union{Nothing,Float64}  # Temperature
    ħ::Float64                  # Planck's constant
end

struct SystemSolution
    k::Float64                  # Spring force constant
    K::Float64                  # Confining force constant
    m::Float64                  # Mass of the chain atoms
    ts::Vector{Float64}         # Time steps
    mem::Float64                # Memory in the units of the slowest chain mode
    a::Float64                  # Chain mass spacing
    M::Float64                  # Mass of the mobile atoms
    F::Float64                  # Magnitude of the Gaussian potential
    s::Float64                  # Standard deviation of the potential
    Rs::Vector{Vector{Float64}} # Positions of mobile atoms.
    # External vector runs over time, internal over masses
    rs::Vector{Vector{Float64}} # Positions of the chain atom
    # External vector runs over time, internal over masses
    ΩT::Union{Nothing,Float64}  # Temperature
    ħ::Float64                  # Planck's constant
end

struct SystemSolutionTest
    k::Float64                  # Spring force constant
    K::Float64                  # Confining force constant
    m::Float64                  # Mass of the chain atoms
    δ::Float64                  # Spacing between time steps
    nPts::Int64                 # Number of time steps
    mem::Float64                # Memory in the units of the slowest chain mode
    a::Float64                  # Chain mass spacing
    M::Float64                  # Mass of the mobile atoms
    F::Float64                  # Magnitude of the Gaussian potential
    s::Float64                  # Standard deviation of the potential
    Rs::Vector{Vector{Float64}} # Positions of mobile atoms.
    # External vector runs over time, internal over masses
    rs::Vector{Vector{Float64}} # Positions of the chain atom
    # External vector runs over time, internal over masses
    ΩT::Union{Nothing,Float64}  # Temperature
    ħ::Float64                  # Planck's constant
end

## Functions
# Frequency as a function of momentum
@inline function Ω(K, k, m, x)
    return sqrt(4 * k / m * sin(x)^2 + K / m)
end

# Chain recoil function
function G(t, ls, K, k, m)
    int_fun(x) = cos.(2 * x * ls) * sin(t * Ω(K, k, m, x)) / Ω(K, k, m, x)
    res = quadgk(int_fun, 0, π / 2)
    return (res[1] * 2 / π / m)
end

# Function for assembling a ChainSystem
function mkChainSystem(K, k, m, t_max, ls, d)
    Ωmax = sqrt(4 * k / m + K / m)      # Largest chain frequency
    δ = (2 * π / Ωmax) / d              # Time step
    n_pts = floor(t_max / δ) |> Int     # Number of time steps given t_max and δ
    # Precomputing the memory term
    G_list = @showprogress pmap(jj -> G(δ * jj, ls, K, k, m), 1:n_pts)
    return ChainSystem(k, K, m, δ, G_list)
end

# Mode amplitude
function ζq(Ωq, ΩT, ħ)
    # Subtract a small number from p. The reason is that for low ΩT, p ≈ 1,
    # causing issues with the rand() generator
    n = rand(Geometric(1 - exp(-Ωq / ΩT) - η))
    res = √(n + 1 / 2) * √(2 * ħ / Ωq)
    return res
end

# Homogeneous displacement of the active chain atom at time step n given a set of Ωs
# and the corresponding ζs and ϕs as a sum of normal coordinates.
function ζH(n, δ, ζs, ϕs, Ωs)
    n_ζ = length(ζs)
    res = [ζs[x] * cos(δ * n * Ωs[x] + ϕs[x]) / √(n_ζ) for x = 1:n_ζ] |> sum
    return res
end

# Function that calculates the trajectories of the mobile atoms and the chain particle.
# mem determines the memory length and τ is the simulation period, both in units of
# the trap period.
function motion_solver(
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
    # Rs = [zeros(length(x0)) for _ = 1:n_pts]# Mobile mass position
    Rs = fill(zeros(length(x0)), n_pts)
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
        curr = ii - 1   # Current time step index
        U_pr = [dU_dr(r - R) for r in rs[curr], R in Rs[curr]]
        U_pr_chain = sum(U_pr, dims = 2) |> vec
        U_pr_mob = -sum(U_pr, dims = 1) |> vec
        rs[ii:ii+min(mem_pts - 1, n_pts - ii)] -=
            G_list[1:min(mem_pts, n_pts - curr)] .* Ref(U_pr_chain .* (mem != 0))
        Rs[ii] = -δ^2 / M .* U_pr_mob + 2 .* Rs[curr] - Rs[curr-1]
    end
    return SystemSolution(k, K, m, ts, mem, tTraj.a, M, F, s, Rs, rs, tTraj.ΩT, tTraj.ħ)
end

## Helper functions for truncated motion_solver function

# General dU_dr function for chain masses
function dU_dr_func(r, F, s)
    return (-F * exp(-r^2 / (2 * s^2)) * r / s^2)
end

# Function to sparsify a vector in place for values below a threshold
function sparsify!(x, thr)
    return x[abs.(x) .< thr] .= 0
end

# Function to return the index of closest value to target for a sorted list
function searchsortednearest(a, x; by=identity, lt=isless, rev=false, distance=(a,b)->abs(a-b))
    i = searchsortedfirst(a, x; by, lt, rev)
    if i == 1
    elseif i > length(a)
        i = length(a)
    elseif a[i] == x
    else
        i = lt(distance(by(a[i]), by(x)), distance(by(a[i - 1]), by(x))) ? i : i - 1
    end
    return i
end

# Function that calculates number of chain masses to track given a dU_dr func
function num_tracked(a::Float64, dU_func::Function, thr::Float64)
    dUs = map(dU_func, a .* range(1, 50, step = 1))
    return searchsortedfirst(dUs, sign(dUs[1]) * thr; rev = (true * sign(dUs[1]) == 1))
end

# Function that returns the indices of chain masses to be tracked
function tracked_range(mob::Vector{Float64}, track_num::Int64, chain_max::Int64, a::Float64)
    chain_pos = a .* 1:chain_max
    mob_pos = map(x -> searchsortednearest(chain_pos, x), mob)
    all_ranges = map(x -> range(max(1, x-track_num), min(x+track_num, chain_max), step = 1), mob_pos)

    # return the set of all ranges and include the first index if not already included
    return all_ranges[1] == 1 ? unique!(sort!(vcat(all_ranges...))) : unique!(sort!(vcat([1], vcat(all_ranges...))))
end

# Calculates trajectories of chain masses and mobile particles
function motion_solver_test(
    system::ChainSystem,
    dU_func::Function,
    F::Real,
    s::Real,
    M::Real,
    x0::Vector{T} where {T<:Real},
    v0::Vector{T} where {T<:Real},
    tTraj::ThermalTrajectory,
    mem::Real,
    τ::Real)

    m = system.m        # Mass of the chain atoms
    k = system.k        # Spring force constant
    K = system.K        # Confining potential force constant
    δ = system.δ        # Array of time steps
    G_list = system.G   # Memory term

    # Check that the thermal trajectory is for the correct system
    if m != tTraj.m || k != tTraj.k || K != tTraj.K || δ != tTraj.δ
        error("Thermal trajectory describes a different system. Check your input.")
    else
        rs = tTraj.rHs
    end

    nChain = length(rs[1])
    Ωmin = √(K / m)                         # Minimum phonon frequency
    tmin = 2 * π / Ωmin                     # Period of the slowest mode
    n_pts = floor(τ * tmin / δ) |> Int      # Number of time steps

    if length(rs) < n_pts
        error("Thermal trajectory provided does not span the necessary time range.")
    end

    mem_pts = max(floor(mem * tmin / δ), 1)  # Memory time points.

    # If the precomputed memory is shorter than the simulation time AND shorter
    # than the desired memory, terminate the calculation.
    if length(G_list) < n_pts && length(G_list) < mem_pts
        error("Chosen memory and the simulation time exceed the precomputed range.")
    else
        # The number of memory pts can be limited by the total simulation time.
        mem_pts = min(mem_pts, n_pts) |> Int
        if mem_pts == n_pts
            @. G_list = δ * G_list
        else
            G_list = δ .* G_list[1:mem_pts]
        end
    end

    # Check number of atoms in G_list
    if length(G_list[1]) < nChain
        error("The recoil term does not contain the desired number of chain masses.")
    elseif length(G_list[1]) == nChain
        G_list = map(SymmetricToeplitz, G_list)
    else
        G_list = map(x -> SymmetricToeplitz(x[1:nChain]), G_list)
    end

    dU_dr(r) = dU_func(r, F, s)

    # Determine number of chain atoms to track from dUdr function
    track_num = num_tracked(tTraj.a, dU_dr, ϵ)

    ## Initialize mobile particle and chain mass arrays
    Rs = fill(zeros(length(x0)), n_pts)
    rs = rs[1:n_pts]

    ## Function to carry out the vector substraction for chain mass vectors
    function add_to_chain_vectors!(
        chain::Vector{Vector{Float64}},
        Gs::Vector{SymmetricToeplitz{Float64}},
        U_list::Vector{Float64},
        ind::Integer;
        tracked = nothing)

        if tracked == nothing
            chain_ind = ind:(ind+min(mem_pts - 2, n_pts - ind))
            g_ind = 2:min(mem_pts, n_pts - 1)
            G_curr = Gs
        else
            chain_ind = ind:(ind+min(mem_pts - 1, n_pts - ind))
            g_ind = 1:min(mem_pts, n_pts - (ind - 1))
            G_curr = view.(Gs, :, Ref(tracked))
        end

        if length(chain_ind) != length(g_ind) && length(G_list) != length(g_ind)
            error("Dimensions do not match")
        end

        # In place replacement of values in chain array
        for (ii, jj) in zip(chain_ind, g_ind)
            chain[ii] = chain[ii] - (G_curr[jj] * U_list)
        end

        return nothing
    end

    ## Function to initialize positions of mobile particles and chain masses
    function initialize_all_positions!()

        # Calculate positions without any truncation
        U_list = [dU_dr(r - R) for r in rs[1], R in Rs[1]]
        U_chain = sum(U_list, dims = 2) |> vec

        add_to_chain_vectors!(rs, G_list, U_chain, 3)

        # Initialize positions of mobile particles
        Rs[1] = x0
        Rs[2] = x0 + v0 * δ

        return nothing
    end

    ## Function to update positions of mobile particles and chain masses
    function update_all_positions!(ind::Integer)
        # Indices of tracked chain atoms
        tracked_indices = tracked_range(Rs[ind-1], track_num, nChain, tTraj.a)

        # dU terms for tracked chain atoms
        U_list = [dU_dr(r - R) for r in getindex(rs[ind-1], tracked_indices), R in Rs[ind-1]]

        # U_list = [dU_dr(r - R) for r in rs[ind-1], R in Rs[ind-1]]

        # Modified dU terms
        U_chain = sum(U_list, dims = 2) |> vec
        U_mob = -sum(U_list, dims = 1) |> vec


        add_to_chain_vectors!(rs, G_list, U_chain, ind; tracked = tracked_indices)
        Rs[ind] = -δ^2 / M .* U_mob + 2 .* Rs[ind - 1] - Rs[ind - 2]

        return nothing
    end


    ## Function containing the loop
    function solve_positions()
        initialize_all_positions!()

        @showprogress for ii in 3:n_pts
            if mem != 0
                update_all_positions!(ii)
            end
        end
    end

    # Solve for positions at each time step
    solve_positions()

    return SystemSolutionTest(k, K, m, δ, n_pts, mem, tTraj.a, M, F, s, Rs, rs, tTraj.ΩT, tTraj.ħ)
end
