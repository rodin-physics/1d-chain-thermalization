include("../src/main.jl")
using Peaks

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

## STITCHING TRAJECTORIES
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
d = 60
τmax = 450                          # Simulation time
τ = 30
total = ceil(τmax / τ) |> Int
δ = system.δ                        # Time step
α = 40                              # Distance between chain atoms
μ = 1
σ0 = [220]

n_pts = τ / δ |> floor |> Int
nChain = 100
ρHs = zeros(nChain, n_pts)
tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
mem = Inf

params = vcat([[(8, 4, [x]), (-8, 4, [x])] for x in [29, 35, 45, 50]]...)
bias = Δ_analytic(50, 8.0, 4.0, system.ωmax)

println("Starting Calculations")
for param in params
    println(param)
    Φ0 = param[1]
    λ = param[2]
    σdot0 = param[3]

    for ind in 1:total
        println(ind)
        if ind == 1 && !isfile(
            "Data/Non_Thermal_Multi/1_Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_bias$(bias)_ΩT$(nothing)_τ$(τ).jld2")

            res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ, threads=true, bias=bias)
            save_object(
            "data/non_thermal_multi/$(ind)_Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_bias$(bias)_ΩT$(nothing)_τ$(τ).jld2",
            res)
        elseif !isfile(
            "Data/Non_Thermal_Multi/$(ind)_Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_bias$(bias)_ΩT$(nothing)_τ$(τ).jld2")

            prev_data = load_object("data/Non_Thermal_Multi/$(ind-1)_Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_bias$(bias)_ΩT$(nothing)_τ$(τ).jld2")

            # Take final max/min speed
            (τs, speeds) = particle_speed(prev_data)

            if sign(Φ0) == 1
                pks, vals = findmaxima(speeds)
            elseif sign(Φ0) == -1
                pks, vals = findminima(speeds)
            end

            final_speed = vals[end]

            # new_σ0 = [mod(last(prev_data.σs), prev_data.α) + (5 * prev_data.α)]

            res = motion_solver(system, Φ0, λ, α, [5.5 * α], [final_speed], μ, tTraj, mem, τ, threads=true, bias=bias)
            save_object(
            "data/non_thermal_multi/$(ind)_Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_bias$(bias)_ΩT$(nothing)_τ$(τ).jld2",
            res)
        end
    end
end
