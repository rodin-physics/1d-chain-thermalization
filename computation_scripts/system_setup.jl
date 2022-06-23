using Random
include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
d = 60       # Number of time steps in the fastest chain mode
τmax = 300    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 300    # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, d)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end

## Prepare the ultrafine ChainSystem's by calculating the recoil term
d = 6000     # Number of time steps in the fastest chain mode
τmax = 20    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 20    # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, d)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end

## Prepare a fine ChainSystem's by calculating the recoil term
d = 1000     # Number of time steps in the fastest chain mode
τmax = 20    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 300   # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, d)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end


## Precompute the thermal trajectories
# τmax = 50                         # Simulation time
# δ = (1 / ωmax) / d              # Time step
# n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
# n_masses = 1000                  # Number of chain masses for simulating ρ0

# qa_s = 2 * pi .* (1:n_masses) / n_masses
# ωs = ω.(ωmax, qa_s ./ 2)
# lines(qa_s, ωs)
# # Generate random phases
# # Random.seed!(150)
# ωT = 1e-5

# ζs = ζq.(ωs, ωT)
# ϕs = 2 * π * rand(n_masses)


# @time ζH(1, δ, ζs, ϕs, ωs, 1:100)



# r1 = [ζH(n, δ, ζs, ϕs, ωs, 10) for n in 1:n_pts]
# r2 = [ζH(n, δ, ζs, ϕs, ωs, 11) for n in 1:n_pts]
# mean(abs.(r1) .* abs.(r1))
# @time pos_corr(0, 0, 0, ωmax)
# @time pos_corr_2(0, 1, ωT, ωmax)
# (std(real(r1))|>sqrt)/sqrt(2)

# lines(1:1000, real(r1)[1:1000])
# std(real(r1)[1:1000])*sqrt(2)
# # Range of temperatures
# ωTs = range(0.0, 1.0, length = 5)

# for ωT in ωTs
#     println("ωT is ", ωT)
#     # Prepare ρHs matrix
#     ρHs = zeros(n_masses, n_pts)
#     ζs = ζq.(ωs, ωT)

#     if (!isfile("precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τ).jld2"))
#         # Populate each row of matrix
#         @showprogress for ii = 1:n_masses
#             ρHs[ii, :] = map(n -> ζH(n, δ, ζs, ϕs, ωs, qs, ii), 1:n_pts)
#         end
#         res = ThermalTrajectory(ωmax, δ, ρHs, ωT)
#         save_object("precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τ).jld2", res)
#     end

# end
