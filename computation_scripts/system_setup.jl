using Random
include("../src/main.jl")

## Prepare the ChainSystem's by calculating the recoil term
# d = 60       # Number of time steps in the fastest chain mode
# τmax = 200    # Simulation time
# ωmax = 10    # Maximum frequency
# lmax = 200    # Number of chain atoms tracked
#
# if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
#     res = mkChainSystem(ωmax, τmax, lmax, d)
#     save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
# end
# println("Next computation")
## Prepare the ultrafine ChainSystem's by calculating the recoil term
d = 6000     # Number of time steps in the fastest chain mode
τmax = 20    # Simulation time
ωmax = 10    # Maximum frequency
lmax = 20    # Number of chain atoms tracked

if (!isfile("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2"))
    res = mkChainSystem(ωmax, τmax, lmax, d)
    save_object("precomputed/systems/System_ωmax$(ωmax)_d$(d)_l$(lmax).jld2", res)
end


## Precompute the thermal trajectories
τ = 10                         # Simulation time
δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
n_masses = 20                  # Number of chain masses for simulating ρ0

qs = range(0, π / 2, length = round(n_masses/2) |> Integer)
ωs = ω.(ωmax, qs)

# Generate random phases
Random.seed!(150)
ϕs = 2 * π * rand(length(qs))
# Range of temperatures
ωTs = range(0.0, 1.0, length = 5)

for ωT in ωTs
    println("ωT is ", ωT)
    # Prepare ρHs matrix
    ρHs = zeros(n_masses, n_pts)
    ζs = ζq.(ωs, ωT)

    if (!isfile("precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τ).jld2"))
        # Populate each row of matrix
        @showprogress for ii in 1:n_masses
            ρHs[ii, :] = map(n -> ζH(n, δ, ζs, ϕs, ωs, qs, ii), 1:n_pts)
        end
        res = ThermalTrajectory(ωmax, δ, ρHs, ωT)
        save_object(
            "precomputed/rH/rH_ωmax$(ωmax)_d$(d)_ωT$(ωT)_τ$(τ).jld2",
            res,
        )
    end

end
