using Distributed
## TODO
# General Example
# Past independence: different starts with infinite memory
# Memory independence
# Scaling: show that for "Instantaneous interactions", narrow potentials are worse at braking.
# Scaling: show that increasing the potential width past a certain point slows down the dissipaiton
# Scaling with Φ for instantaneous interactions
# Scaling with μ?
# Non-instantaneous interaction scaling? General formula?
# "Same potential landscape" from different λ and Φ
# 
proc_num = 6
addprocs(proc_num - nprocs())

# include("../src/main.jl")
@everywhere include("../src/main.jl")

# GENERAL EXAMPLE
system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")
d = 60
τ = 100                             # Simulation time
δ = system.δ                        # Time step
α = 10                              # Distance between chain atoms
μ = 1

n_pts = τ / δ |> floor |> Int
nChain = 100
ρHs = zeros(nChain, n_pts)
tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
σdot0 = [20]
mem = Inf

σ0 = [55]
Φ0 = [1 / 2, 1, 2]
λ = [1 / 2, 1]
params = [(f, l) for f in Φ0, l in λ] |> vec

@showprogress pmap(params) do param
    Φ0 = param[1]
    λ = param[2]
    if (
        !isfile(
            "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ)
        save_object(
            "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

σ0 = [60]
Φ0 = -[1 / 2, 1, 2]
λ = [1 / 2, 1]
params = [(f, l) for f in Φ0, l in λ] |> vec

@showprogress pmap(params) do param
    Φ0 = param[1]
    λ = param[2]
    if (
        !isfile(
            "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
        )
    )

        res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, mem, τ)
        save_object(
            "data/non_thermal/Single_σ0$(σ0)_σdot0$(σdot0)_Mem$(mem)_λ$(λ)_Φ$(Φ0)_μ$(μ)_d$(d)_ΩT$(nothing)_τ$(τ).jld2",
            res,
        )
    end
end

## Burst Energy Loss
system = load_object("precomputed/systems/System_ωmax10_d3000_l10.jld2")

d = 3000
δ = system.δ
α = 10
Φ0 = 1/2
σ0 = 45
nChain = 10
vPts = 100
σ_dots = range(20, 350, length = vPts)
Δs = zeros(vPts)


@showprogress for ii = 1:vPts
    σ_dot = σ_dots[ii]
    τ = 1.1 * (α / σ_dot)
    n_pts = τ / δ |> floor |> Int
    ρHs = zeros(nChain, n_pts)
    tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

    res = motion_solver(system, Φ0, 1/4, α, [σ0], [σ_dot], 1, tTraj, Inf, τ)
    σs = res.σs |> vec
    mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
    v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
    Δs[ii] = μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2)
end

λ = 1/4
Φ = 1/2
function en_loss(v)
    z = (2 * π * λ / v)^2
    return (
        4 * π^3 * Φ^2 / v^2 *
        z *
        exp(-z * (Ω^2 + 1) / 2) *
        (
            besseli(0, z * (Ω^2 - 1) / 2) +
            (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
        )
    )
end
analytic = en_loss.(σ_dots)

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")

# lines!(ax1,log.(σ_dots), log.(Δs))
# # lines!(ax1,log.(σ_dots), log.(Δs_2))
# lines!(ax1, log.(σ_dots), log.(analytic))
# Δs- Δs_2


lines!(ax1,(σ_dots), (Δs))
# lines!(ax1,(σ_dots), (Δs_2))
lines!(ax1, (σ_dots), (analytic))
fig


μs = [1 / 5, 1, 5]
λs = [1 / 8, 1 / 4, 1 / 2, 1, 2]

σ_dot = 25


res = motion_solver(system, Φ0, λ, α, σ0, σdot0, μ, tTraj, Inf, τ)
σs = res.σs |> vec
v = (σs[2:end] - σs[1:end-1]) / (res.τs[2] - res.τs[1])
scatter((res.τs)[2:end], v)


mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
Δ_en = μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2)
μ
σs[mid_pt_idx]


# system = load_object("precomputed/systems/System_ωmax10_d3000_l10.jld2")

# d = 3000
# τ = .4                            # Simulation time
# δ = system.δ                        # Time step
# a = 10                               # Distance between chain atoms
# n_pts = τ / δ |> floor |> Int
# nChain = 10
# ρHs = zeros(nChain, n_pts)
# tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
# λ = 1;
# Φ = 1;
# μ = 25/1;
# @time res = motion_solver(system, Φ, λ, a, [25], [25], μ, tTraj, Inf, τ)
# s = (res.σs |> vec)
# v = (s[2:end] - s[1:end-1]) / (res.τs[2] - res.τs[1])
# scatter((res.τs)[2:end], v)
# kin_en = μ * v.^ 2 / 2 / (2 * π)^2
# kin_en[end] - kin_en[1]
