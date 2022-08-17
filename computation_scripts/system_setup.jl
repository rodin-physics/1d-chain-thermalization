@everywhere using Random
@everywhere include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
d = 60       # Number of time steps in the fastest chain mode
τmax = 500    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 300    # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, d)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end

## Prepare the ultrafine ChainSystem's by calculating the recoil term
# d = 6000     # Number of time steps in the fastest chain mode
# τmax = 20    # Simulation time
# ωmax = 10    # Maximum frequency
# lmax = 20    # Number of chain atoms tracked
#
# if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
#     res = mkChainSystem(ωmax, τmax, lmax, d)
#     save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
# end

# println("Calculating thermal trajectories")
# ## Precompute the thermal trajectories
# τmax = 200                         # Simulation time
# δ = (1 / ωmax) / d              # Time step
# n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
# n_modes = 10000                  # Number of chain masses for simulating ρ0
#
# qa_s = 2 * pi .* (1:n_modes) / n_modes
# ωs = ω.(ωmax, qa_s ./ 2)
#
# Random.seed!(150)
# ϕs = 2 * π * rand(n_modes)
#
# # Range of temperatures
# ωTs = [5.0]
#
# for ωT in ωTs
#     println("ωT is ", ωT)
#     ζs = ζq.(ωs, ωT)
#
#     if (!isfile("precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax).jld2"))
#         # Populate each row of matrix
#         gs = collect(1:lmax)
#         full_res = @showprogress pmap(n -> real(ρH(n, δ, ζs, ϕs, ωs, gs)), 1:n_pts)
#         full_res = reduce(hcat, full_res)
#
#         res = ThermalTrajectory(ωmax, δ, full_res, ωT)
#         save_object("precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τmax).jld2", res)
#     end
#
# end
