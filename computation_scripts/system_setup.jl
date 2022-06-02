using Random
include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
d = 60       # Number of time steps in the fastest chain mode
τmax = 200    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 200    # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, d)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end

## Prepare the ultrafine ChainSystem's by calculating the recoil term
d = 6000     # Number of time steps in the fastest chain mode
τmax = 4     # Simulation time
ωmax = 10    # Maximum frequency
lmax = 10    # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, d)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end


# ## Precompute the thermal trajectories
# m = 1                           # Mass of the chain atoms
# τ = 1000;
# t_max = τ * t_M                 # Simulation time

# Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency
# δ = (2 * π / Ωmax) / d          # Time step
# n_pts = floor(t_max / δ) |> Int # Number of time steps given t_max and δ

# ## Thermal Trajectory
# n_masses = 1000000              # Number of chain masses for simulating r0
# qs = range(0, π / 2, length = round(n_masses / 2) |> Integer)
# Ωs = Ω.(K, k, m, qs)
# ΩTs = [1e-5, 1, 2, 5, 10, 20, 50, 100, 150, 250, 500, 1000]
# for ii = 1:length(ΩTs)
#     println(ii)
#     ΩT = ΩTs[ii]

#     if (!isfile("precomputed/rH/rH_K$(K)_k$(k)_m$(m)_d$(d)_ΩT$(ΩT)_τ$(τ)_hbar$(ħ).jld2"))
#         # Seeding the RNG
#         Random.seed!(150)
#         ϕs = 2 * π * rand(length(qs))
#         ζs = ζq.(Ωs, ΩT, ħ)
#         rHs = @showprogress pmap(n -> ζH(n, δ, ζs, ϕs, Ωs) / √(m), 1:n_pts)
#         res = ThermalTrajectory(k, K, m, δ, rHs, ΩT, ħ)
#         save_object(
#             "precomputed/rH/rH_K$(K)_k$(k)_m$(m)_d$(d)_ΩT$(ΩT)_τ$(τ)_hbar$(ħ).jld2",
#             res,
#         )
#     end

# end
