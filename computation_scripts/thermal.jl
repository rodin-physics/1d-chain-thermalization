include("../src/main.jl")

# # GENERAL EXAMPLE
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
d = 60
τ = 100                             # Simulation time
δ = system.δ                        # Time step
α = 10                              # Distance between chain atoms
μ = 1

n_pts = τ / δ |> floor |> Int
nChain = 50

tTraj = load_object("precomputed/rH/rH_ωmax10_d60_ωT500.0_τ200.jld2")
ωT = tTraj.ωT

σdot0 = [0]
mem = Inf

σ0 = [225]
Φ0 = [-1, 1]
λ = [1]
params = [(f, l) for f in Φ0, l in λ] |> vec
println("Starting Calculations")
Threads.@threads for param in params
    println(param)
    Φ0 = param[1]
    λ = param[2]
    if (
        !isfile(
            "data/thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(ωT)_τ$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ)
        save_object(
            "data/thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(ωT)_τ$(τ).jld2",
            res,
        )
    end
end
