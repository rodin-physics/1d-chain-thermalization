
## TRAJECTORY ENERGY LOSS

fig = Figure(resolution=(1200, 800), font="CMU Serif", fontsize=36)
ax1 = Axis(fig[1, 1], xlabel=L"\dot{\sigma}", ylabel=L"\Delta")

# REPULSIVE
data = load_object(
    "data/non_thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ0.5_μ1_d60_ΩTnothing_τ120.jld2",
)
δ = data.τs[2] - data.τs[1]
σs = data.σs |> vec
σ_dots = (σs[2:end] - σs[1:end-1]) ./ δ

# Get the data before the first turning point
idx = findfirst(σ_dots .< 0)
idx = isnothing(idx) ? length(data.τs) - 1 : idx
ts = data.τs[1:idx]
σ_dots = σ_dots[1:idx]

# lines(ts, σ_dots)

σs = data.σs[1:idx] |> vec
max_pos =  Int(maximum(σs) ÷ 10)
idx = [argmin(abs.(σs .- n * 10 .- 5)) for n in 5:max_pos]
σ_dots_max = σ_dots[idx]
KE = 0.5 * (σ_dots_max ./ 2 ./ π) .^ 2
Δs = KE[1:end-1] - KE[2:end]



# findlocalmaxima(σ_dots)
# # Get all the indices when the particle is at the minimum velocity
# σ_dots_max = prepend!(maxima(σ_dots),20)
# idx = maxima(σ_dots)
# idx = prepend!(idx, 1)
# # Compute the kinetic energy at the maximum 
# σ_dots_max = σ_dots[idx]

# lines!(ax1, σ_dots_max, KE)
# lines!(ax1, σ_dots_max[1:end-1], Δs)
scatter!(ax1, σ_dots_max[1:end-1], Δs, marker='⊕', markersize=20)
fig
# σ_dots_burst = σ_dots_max[findall(σ_dots_max .> 10)]
# Δs = [Δ_numeric(x, σs[1], data.Φ, data.λ, system_slow) for x in σ_dots_burst]
# scatter!(ax1, σ_dots_burst, Δs, marker=:cross, markersize=20)
# # ATTRACTIVE

# data = load_object(
#     "data/non_thermal/Single_σ0[55]_σdot0[20]_MemInf_λ0.5_Φ-1.0_μ1_d60_ΩTnothing_τ120.jld2",
# )
# δ = data.τs[2] - data.τs[1]
# σs = data.σs |> vec
# σ_dots = (σs[2:end] - σs[1:end-1]) ./ δ
# # Get the data before the first turning point
# idx = findfirst(σ_dots .< 0)
# ts = data.τs[1:idx]
# σ_dots = σ_dots[1:idx]

# # Get all the indices when the particle is at the maximum velocity
# idx = findall(n -> σ_dots[n] < σ_dots[n-1] && σ_dots[n] < σ_dots[n+1], 2:length(σ_dots)-1)
# idx = prepend!(idx, 1)

# # Compute the kinetic energy at the maximum 
# σ_dots_max = σ_dots[idx]
# KE = 0.5 * (σ_dots_max ./ 2 ./ π) .^ 2
# Δs = KE[1:end-1] - KE[2:end]
# scatter!(ax1, σ_dots_max[1:end-1], Δs, marker='⊖', markersize=20)

# σ_dots_burst = σ_dots_max[findall(σ_dots_max .> 10)]
# Δs = [Δ_numeric(x, σs[1], data.Φ, data.λ, system_slow) for x in σ_dots_burst]
# scatter!(ax1, σ_dots_burst, Δs, marker=:hline, markersize=20)
# ## ANALYTIC
σ_dots = range(6, 25, length=100)
Δs = Δ_analytic.(σ_dots, data.Φ, data.λ, data.ωmax)
lines!(ax1, σ_dots, Δs, linewidth=4)

# xlims!(ax1, (10, 20))
# ylims!(ax1, (0, 0.15))

fig



data = load_object(
    "data/non_thermal/Single_σ0[55]_σdot0[25]_MemInf_λ1.0_Φ1.0_μ1_d60_ΩTnothing_τ200.jld2",
)

ρs = data.ρs
r = ρs[7,10000:13000]
length(r)
lines(data.τs[10000:13000],r)

ρs = data.ρs
r = ρs[7,1:3000]
length(r)
lines(data.τs[1:3000],r)
